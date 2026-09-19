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

// Split [0, count) into contiguous chunks, one per worker, and invoke
// body(begin, end) on each.  Chunks are contiguous rather than interleaved so
// every worker walks memory forwards; the calling thread takes the last chunk
// itself rather than sitting idle.
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
    const std::size_t per_worker = count / workers;
    const std::size_t remainder = count % workers;

    std::vector<std::thread> pool;
    pool.reserve(workers - 1);
    std::size_t begin = 0;
    for (std::size_t w = 0; w < workers; w++)
    {
        // Spread the remainder over the first few workers rather than piling it
        // all onto the last one.
        const std::size_t end = begin + per_worker + (w < remainder ? 1 : 0);
        if (w + 1 < workers)
            pool.emplace_back([&body, begin, end] { body(begin, end); });
        else
            body(begin, end);
        begin = end;
    }
    for (auto& worker : pool)
        worker.join();
}

} // namespace niceops
