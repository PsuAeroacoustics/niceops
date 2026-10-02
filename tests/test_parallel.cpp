#include <atomic>
#include <catch2/catch_test_macros.hpp>
#include <stdexcept>
#include <string>
#include <vector>

#include "util/parallel.hpp"

using namespace niceops; // NOLINT

TEST_CASE("parallel_chunks visits every index exactly once", "[parallel]")
{
    for (unsigned threads : {1u, 2u, 3u, 8u, 64u})
    {
        for (std::size_t count : {std::size_t {0}, std::size_t {1}, std::size_t {7},
                     std::size_t {100}, std::size_t {1001}})
        {
            std::vector<std::atomic<int>> visits(count);
            parallel_chunks(count, threads, [&](std::size_t begin, std::size_t end) {
                REQUIRE(begin < end);
                for (std::size_t i = begin; i < end; i++)
                    visits[i]++;
            });
            for (std::size_t i = 0; i < count; i++)
                REQUIRE(visits[i] == 1);
        }
    }
}

TEST_CASE("parallel_chunks carries a worker's exception to the caller", "[parallel]")
{
    // An exception cannot leave a std::thread's function, so before it was
    // caught in the worker this one terminated the process at any thread count
    // above one.  The test binary surviving to the REQUIRE_THROWS is the point.
    for (unsigned threads : {1u, 2u, 4u, 14u})
    {
        std::atomic<std::size_t> started {0};
        REQUIRE_THROWS_AS(parallel_chunks(10000, threads,
                                  [&](std::size_t, std::size_t end) {
                                      started++;
                                      if (end > 2000)
                                          throw std::runtime_error("worker failed");
                                  }),
                std::runtime_error);
        // With more than one thread the remaining chunks are abandoned rather
        // than run to the end after the failure; one thread stops where it threw.
        REQUIRE(started.load() >= 1);
    }
}

TEST_CASE("parallel_chunks rethrows the exception whichever chunk throws it", "[parallel]")
{
    // Every chunk throws, on every thread at once: still exactly one exception
    // reaches the caller and the pool is joined cleanly.
    std::string message;
    try
    {
        parallel_chunks(1000, 8, [](std::size_t, std::size_t) { throw std::runtime_error("all fail"); });
    }
    catch (const std::runtime_error& e)
    {
        message = e.what();
    }
    REQUIRE(message == "all fail");
}
