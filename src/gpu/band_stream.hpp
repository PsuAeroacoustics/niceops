#pragma once

// Ray sets streamed from worker threads to the GPU (BandEnergies): workers
// push each set as they pack it and carry on; one thread takes them a batch
// at a time, keeps two batches in flight, and hands each set's energies back
// as its batch finishes.  The queue is bounded, so a GPU that falls behind
// holds the workers back rather than the memory growing; a GPU failure wakes
// them, and push() throws.

#include <condition_variable>
#include <cstddef>
#include <deque>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>

#include "gpu/band_energy.hpp"
#include "gpu/point_rays.hpp"

namespace niceops::gpu {

class BandStream
{
public:
    // `done(tag, energy)` is called on the stream's thread for each set, its
    // energies by problem in the set's order; `batch` sets per batch.
    using Done = std::function<void(std::size_t tag, const double* energy)>;
    BandStream(Device& device, std::size_t batch, Done done);
    ~BandStream();
    BandStream(const BandStream&) = delete;
    BandStream& operator=(const BandStream&) = delete;

    // Hands a set over; blocks while the queue is full.  Throws
    // std::runtime_error where the GPU has failed.
    void push(std::size_t tag, PointRaySet&& set);
    // Waits for every set pushed to be done; rethrows a GPU failure.
    void finish();

private:
    struct Packed
    {
        std::size_t tag = 0;
        PointRaySet set;
    };
    void run();
    void stop();

    std::size_t batch_;
    Done done_;
    BandEnergies energies_[2];
    std::mutex lock_;
    std::condition_variable ready_, space_;
    std::deque<std::unique_ptr<Packed>> queue_;
    bool producing_ = true;
    std::exception_ptr failure_;
    std::thread thread_;
};

} // namespace niceops::gpu
