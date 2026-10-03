#pragma once

// Minimal chunked parallel-for.  The observer loop is the only thing in
// NICE-OPS that needs one: every observer's time history is computed from
// read-only data (the database, and the per-state SourceState array built
// before the loop) and written to its own slot, so the work splits with no
// synchronization at all.
//
// std::thread rather than <execution>: libc++ ships the parallel algorithms
// only against TBB, which would be a whole new dependency for one loop.
//
// The body is handed a [begin, end) range rather than a single index so the
// caller can hoist per-worker scratch -- the level time history, in practice --
// out of the inner loop instead of allocating it per observer.
//
// An exception thrown by the body is carried back to the caller: the first one
// is rethrown from parallel_chunks once every worker has stopped, exactly as it
// would have propagated had the body run inline.  (This used to terminate the
// process whenever there was more than one thread, because an exception cannot
// leave a std::thread's function -- which defeated the caller's per-trajectory
// handling of input_error, the one the observer loop reaches when a history
// file cannot be opened.)  Remaining chunks are abandoned; their output slots
// are whatever the caller initialised them to, and the caller is expected to
// discard the lot.

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace niceops {

// How many workers to use when the caller does not say.  hardware_concurrency
// is allowed to return 0 when it cannot tell.
inline unsigned default_thread_count()
{
    const unsigned hardware = std::thread::hardware_concurrency();
    return hardware > 0 ? hardware : 1u;
}

// How many chunks to cut per worker.  Enough that a slow core holds up only its
// own share rather than a whole equal slice of the work, few enough that the
// per-chunk setup a body does (allocating its scratch history) stays noise.
inline constexpr std::size_t chunks_per_worker = 8;

// Split [0, count) into contiguous chunks and invoke body(begin, end) on each,
// handing them out to `threads` workers as they finish.  Chunks are contiguous
// rather than interleaved so every worker walks memory forwards; the calling
// thread pulls chunks itself rather than sitting idle.
//
// They are dealt out on demand rather than divided up front because the cores
// are not interchangeable.  An equal split gives every worker the same slice,
// so the wall time is the slowest core's -- on a big.LITTLE machine (every
// Apple silicon part, and most recent x86 laptops) the efficiency cores run
// this loop several times slower than the performance ones, and a run on all
// 14 threads of an M4 Pro came out no faster than one on its 10 performance
// cores alone.  Pulling chunks lets the fast cores simply take more of them.
//
// The result does not depend on the split: every index writes only its own
// output, so runs at any thread count stay byte-for-byte identical.
//
// With one thread the body runs inline on the caller, so an exception simply
// propagates; with more, the first exception any worker throws is rethrown on
// the calling thread after the pool has been joined.
template <typename Body>
void parallel_chunks(std::size_t count, unsigned threads, Body body)
{
    if (count == 0)
        return;
    if (threads <= 1)
    {
        body(std::size_t {0}, count);
        return;
    }

    const std::size_t workers = std::min<std::size_t>(threads, count);
    const std::size_t slices = workers * chunks_per_worker;
    // Round up, so the chunks cover `count` and the last one is the short one.
    // count >= 1 here, so the ceiling is already at least 1.
    const std::size_t chunk = (count + slices - 1) / slices;

    std::atomic<std::size_t> next {0};
    std::atomic<bool> stop {false};
    std::exception_ptr failure;
    std::mutex failure_mutex;
    auto pull = [&body, &next, &stop, &failure, &failure_mutex, count, chunk] {
        for (;;)
        {
            if (stop.load(std::memory_order_relaxed))
                return;
            const std::size_t begin = next.fetch_add(chunk, std::memory_order_relaxed);
            if (begin >= count)
                return;
            try
            {
                body(begin, std::min(begin + chunk, count));
            }
            catch (...)
            {
                {
                    std::lock_guard<std::mutex> lock(failure_mutex);
                    if (!failure)
                        failure = std::current_exception();
                }
                stop.store(true, std::memory_order_relaxed);
                return;
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(workers - 1);
    for (std::size_t w = 0; w + 1 < workers; w++)
        pool.emplace_back(pull);
    pull();
    for (auto& worker : pool)
        worker.join();
    if (failure)
        std::rethrow_exception(failure);
}

} // namespace niceops
