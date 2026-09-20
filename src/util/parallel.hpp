#pragma once

// Minimal chunked parallel-for.  The observer loop is the only thing in
// NICE-OPS that needs one: every observer's time history is computed from
// read-only data (the database, and the per-state SourceState array built
// before the loop) and written to its own slot, so the work splits with no
// synchronisation at all.
//
// std::thread rather than <execution>: libc++ ships the parallel algorithms
// only against TBB, which would be a whole new dependency for one loop.
//
// The body is handed a [begin, end) range rather than a single index so the
// caller can hoist per-worker scratch -- the level time history, in practice --
// out of the inner loop instead of allocating it per observer.
//
// Note that an exception escaping the body terminates the process.  That
// matches what already happens: fatal() aborts, and nothing here throws.

#include <algorithm>
#include <atomic>
#include <cstddef>
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
    // Round up, so the chunks cover `count` and the last one is the short one.
    const std::size_t chunk
            = std::max<std::size_t>(1, (count + workers * chunks_per_worker - 1)
                            / (workers * chunks_per_worker));

    std::atomic<std::size_t> next {0};
    auto pull = [&body, &next, count, chunk] {
        for (;;)
        {
            const std::size_t begin = next.fetch_add(chunk, std::memory_order_relaxed);
            if (begin >= count)
                return;
            body(begin, std::min(begin + chunk, count));
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(workers - 1);
    for (std::size_t w = 0; w + 1 < workers; w++)
        pool.emplace_back(pull);
    pull();
    for (auto& worker : pool)
        worker.join();
}

} // namespace niceops
