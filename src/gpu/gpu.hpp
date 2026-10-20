#pragma once

// A GPU, through whichever backend this build has (NICEOPS_GPU; Metal on
// macOS today).  The kernels are written once, in the dialect of
// src/gpu/kernels -- each backend's prelude maps it onto that backend's
// language -- and compiled at run time from source embedded in the binary, so
// a build needs no shader toolchain.  Buffers are shared with the CPU.
//
// Nothing here is a backend's type, so the code that drives a kernel (packing
// its inputs, dispatching it, reading its results) is written once too.

#include <cstddef>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace niceops::gpu {

class Device;
class Batch;

// Memory the CPU and the GPU both see.  data() is the CPU's view: write a
// buffer's inputs before the batch that reads them is committed, and read its
// outputs after that batch's wait().
class Buffer
{
public:
    Buffer() = default;
    std::size_t size() const;
    void* data() const;
    template <typename T>
    T* as() const
    {
        return static_cast<T*>(data());
    }
    explicit operator bool() const
    {
        return impl_ != nullptr;
    }
    struct Impl;

private:
    friend class Device;
    friend class Batch;
    std::shared_ptr<Impl> impl_;
};

// One argument of a dispatch, in the kernel's slot order: a buffer, or a
// small value (at most 4 KB) copied when it is set.
struct Arg
{
    const Buffer* buffer = nullptr;
    std::size_t offset = 0;
    std::vector<unsigned char> bytes;

    Arg(const Buffer& b, std::size_t byte_offset = 0) : buffer(&b), offset(byte_offset) {}
    template <typename T>
    static Arg value(const T& v)
    {
        Arg a;
        a.bytes.resize(sizeof(T));
        std::memcpy(a.bytes.data(), &v, sizeof(T));
        return a;
    }

private:
    Arg() = default;
};

// A compiled kernel entry point.
class Kernel
{
public:
    Kernel() = default;
    // The most threads a group of this kernel may have on this device.
    std::size_t max_threads_per_group() const;
    struct Impl;

private:
    friend class Device;
    friend class Batch;
    std::shared_ptr<Impl> impl_;
};

// Dispatches run on the GPU in the order given, after commit(); wait() blocks
// until they are done and throws if the GPU reported an error.
class Batch
{
public:
    void dispatch(const Kernel& kernel, const std::vector<Arg>& args, std::size_t groups, std::size_t threads_per_group);
    void commit();
    void wait();
    // GPU time between the batch's first dispatch starting and its last
    // finishing, s; valid after wait().
    double gpu_seconds() const;
    struct Impl;

private:
    friend class Device;
    std::shared_ptr<Impl> impl_;
};

class Device
{
public:
    // The default device, or null where this build has no GPU backend or the
    // machine no device.  Created on first use and kept for the process.
    static Device* get();
    // The backend's name ("metal") or "none".
    static const char* backend();
    // get()'s device, or a std::runtime_error saying why there is none, for
    // the option that asked for it ("--gpu").
    static Device& require(const std::string& option)
    {
        if (Device* device = get())
            return *device;
        throw std::runtime_error(option + ": "
                + (std::string(backend()) == "none" ? "this build has no GPU backend (NICEOPS_GPU)" : "no GPU device"));
    }

    std::string name() const;
    // The largest buffer the device takes, and the memory it can keep in use
    // at once, bytes: what a caller checks before choosing the GPU for a big
    // problem.
    std::size_t max_buffer_bytes() const;
    std::size_t working_set_bytes() const;
    // Zeroed shared memory of `bytes`.
    Buffer buffer(std::size_t bytes);
    // Entry point `entry` of the program made of the named kernel sources in
    // order, after the backend's prelude and common.h.  Compiled once and
    // cached; throws std::runtime_error with the compiler's message if it does
    // not compile.
    Kernel kernel(const std::vector<std::string>& sources, const std::string& entry);
    Batch batch();
    struct Impl;

private:
    std::shared_ptr<Impl> impl_;
};

// The kernel sources embedded at build time (cmake/embed_kernels.cmake), by
// file name; null for a name not embedded.
const char* kernel_source(const std::string& name);

// x as a float-float pair: the nearest float, and what it leaves, which a
// kernel adds back where a float alone would lose it (common.h's nk_phase).
inline void split(double x, float& hi, float& lo)
{
    hi = static_cast<float>(x);
    lo = static_cast<float>(x - static_cast<double>(hi));
}

} // namespace niceops::gpu
