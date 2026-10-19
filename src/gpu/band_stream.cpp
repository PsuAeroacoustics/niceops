#include "gpu/band_stream.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

namespace niceops::gpu {

BandStream::BandStream(Device& device, std::size_t batch, Done done)
    : batch_(std::max<std::size_t>(batch, 1)), done_(std::move(done)), energies_ {BandEnergies(device), BandEnergies(device)}
{
    thread_ = std::thread([this] { run(); });
}

BandStream::~BandStream()
{
    stop();
}

void BandStream::push(std::size_t tag, PointRaySet&& set)
{
    auto packed = std::make_unique<Packed>();
    packed->tag = tag;
    packed->set = std::move(set);
    std::unique_lock<std::mutex> guard(lock_);
    space_.wait(guard, [&] { return queue_.size() < 2 * batch_ || failure_; });
    if (failure_)
        throw std::runtime_error("the GPU stopped");
    queue_.push_back(std::move(packed));
    ready_.notify_one();
}

void BandStream::finish()
{
    stop();
    if (failure_)
        std::rethrow_exception(failure_);
}

void BandStream::stop()
{
    {
        std::lock_guard<std::mutex> guard(lock_);
        producing_ = false;
    }
    ready_.notify_all();
    if (thread_.joinable())
        thread_.join();
}

void BandStream::run()
{
    try
    {
        std::vector<std::unique_ptr<Packed>> held[2];
        auto take = [&](std::vector<std::unique_ptr<Packed>>& into) {
            std::unique_lock<std::mutex> guard(lock_);
            ready_.wait(guard, [&] { return queue_.size() >= batch_ || !producing_; });
            while (!queue_.empty() && into.size() < batch_)
            {
                into.push_back(std::move(queue_.front()));
                queue_.pop_front();
            }
            space_.notify_all();
            return !into.empty();
        };
        auto finish = [&](std::vector<std::unique_ptr<Packed>>& from, BandEnergies& batch) {
            const std::vector<double>& energy = batch.wait();
            std::size_t q = 0;
            for (const auto& packed : from)
            {
                done_(packed->tag, energy.data() + q);
                q += packed->set.problems();
            }
            from.clear();
            batch.clear();
        };
        std::size_t current = 0;
        bool in_flight = false;
        while (take(held[current]))
        {
            for (const auto& packed : held[current])
                energies_[current].add(packed->set);
            energies_[current].submit();
            if (in_flight)
                finish(held[1 - current], energies_[1 - current]);
            in_flight = true;
            current = 1 - current;
        }
        if (in_flight)
            finish(held[1 - current], energies_[1 - current]);
    }
    catch (...)
    {
        std::lock_guard<std::mutex> guard(lock_);
        failure_ = std::current_exception();
        space_.notify_all();
    }
}

} // namespace niceops::gpu
