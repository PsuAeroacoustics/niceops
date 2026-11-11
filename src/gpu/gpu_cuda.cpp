// The CUDA backend of gpu.hpp.  The kernels are compiled from source at run
// time by NVRTC to PTX, which the driver loads; both libraries are opened
// with dlopen when the device is first asked for, so a build with this backend
// links neither and runs, on the CPU, where there is no driver.
//
// Buffers are managed memory, which the CPU and the GPU both address, as
// Metal's shared buffers are: the callers write a batch's inputs while the
// previous batch runs, so a device must allow concurrent managed access (any
// NVIDIA GPU from Pascal on Linux); one that does not is no device.  A batch
// records its dispatches and launches them on the device's one stream at
// commit(), between two events, so batches run in commit order as a Metal
// queue's command buffers do.

#include "gpu/gpu.hpp"

#include <cuda.h>
#include <dlfcn.h>
#include <nvrtc.h>

#include <algorithm>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>

namespace niceops::gpu {

namespace {

// The entry points used, looked up by the name cuda.h maps each to (its _v2
// where there is one).
#define NK_STR(x) #x
#define NK_XSTR(x) NK_STR(x)
#define NK_CUDA_FUNCTIONS(X)                                                                                         \
    X(cuInit) X(cuGetErrorString) X(cuDeviceGet) X(cuDeviceGetCount) X(cuDeviceGetName) X(cuDeviceGetAttribute)      \
    X(cuDeviceTotalMem) X(cuDevicePrimaryCtxRetain) X(cuCtxSetCurrent) X(cuMemGetInfo) X(cuMemAllocManaged)            \
    X(cuMemFree) X(cuMemsetD8Async) X(cuStreamCreate) X(cuStreamSynchronize) X(cuEventCreate) X(cuEventRecord)         \
    X(cuEventSynchronize) X(cuEventElapsedTime) X(cuEventDestroy) X(cuModuleLoadDataEx) X(cuModuleGetFunction)         \
    X(cuFuncGetAttribute) X(cuLaunchKernel)
#define NK_NVRTC_FUNCTIONS(X)                                                                                        \
    X(nvrtcCreateProgram) X(nvrtcCompileProgram) X(nvrtcGetProgramLogSize) X(nvrtcGetProgramLog) X(nvrtcGetPTXSize)    \
    X(nvrtcGetPTX) X(nvrtcDestroyProgram) X(nvrtcGetErrorString) X(nvrtcGetNumSupportedArchs) X(nvrtcGetSupportedArchs)

// Each function's pointer type, named outside Api: a member may not take its
// type from the function it is named after.
#define NK_TYPE(name) using name##_type = decltype(&name);
NK_CUDA_FUNCTIONS(NK_TYPE)
NK_NVRTC_FUNCTIONS(NK_TYPE)
#undef NK_TYPE

struct Api
{
#define NK_MEMBER(name) name##_type name = nullptr;
    NK_CUDA_FUNCTIONS(NK_MEMBER)
    NK_NVRTC_FUNCTIONS(NK_MEMBER)
#undef NK_MEMBER
    // Optional: placing a buffer on the device ahead of a batch.  Absent from
    // some drivers' exports, and a hint only.
    CUresult (*prefetch)(CUdeviceptr, size_t, CUdevice, CUstream) = nullptr;
    // Why the API is not there, or "" where it is.
    std::string missing;
};

void* open_first(std::initializer_list<std::string> names)
{
    for (const std::string& name : names)
        if (void* h = dlopen(name.c_str(), RTLD_NOW | RTLD_LOCAL))
            return h;
    return nullptr;
}

const Api& api()
{
    static const Api loaded = [] {
        Api a;
        void* cuda = open_first({"libcuda.so.1", "libcuda.so"});
        if (!cuda)
        {
            a.missing = "the CUDA driver library (libcuda.so.1) is not installed";
            return a;
        }
        // The NVRTC of the headers this was built with first: its ABI is stable
        // within a major version.
        void* nvrtc = open_first({"libnvrtc.so." + std::to_string(CUDA_VERSION / 1000), "libnvrtc.so"});
        if (!nvrtc)
        {
            a.missing = "the CUDA run-time compiler (libnvrtc.so) is not on the library path";
            return a;
        }
#define NK_LOAD_FROM(library)                                                                                        \
    [&](auto& slot, const char* symbol) {                                                                            \
        slot = reinterpret_cast<std::remove_reference_t<decltype(slot)>>(dlsym(library, symbol));                     \
        if (!slot && a.missing.empty())                                                                              \
            a.missing = std::string("no ") + symbol + " in the CUDA libraries";                                       \
    }
        auto from_cuda = NK_LOAD_FROM(cuda);
        auto from_nvrtc = NK_LOAD_FROM(nvrtc);
#undef NK_LOAD_FROM
#define NK_LOAD_CUDA(name) from_cuda(a.name, NK_XSTR(name));
#define NK_LOAD_NVRTC(name) from_nvrtc(a.name, NK_XSTR(name));
        NK_CUDA_FUNCTIONS(NK_LOAD_CUDA)
        NK_NVRTC_FUNCTIONS(NK_LOAD_NVRTC)
#undef NK_LOAD_CUDA
#undef NK_LOAD_NVRTC
        a.prefetch = reinterpret_cast<decltype(a.prefetch)>(dlsym(cuda, "cuMemPrefetchAsync"));
        return a;
    }();
    return loaded;
}

std::string cuda_error(CUresult r)
{
    const char* text = nullptr;
    if (api().cuGetErrorString && api().cuGetErrorString(r, &text) == CUDA_SUCCESS && text)
        return text;
    return "CUDA error " + std::to_string(static_cast<int>(r));
}

void check(CUresult r, const char* what)
{
    if (r != CUDA_SUCCESS)
        throw std::runtime_error(std::string("GPU: ") + what + ": " + cuda_error(r));
}

// The program's PTX for compute capability `arch` (major * 10 + minor), or a
// std::runtime_error carrying NVRTC's log.  The architecture is the highest
// this NVRTC knows at or below the device's, so an older toolkit still serves
// a newer card through the driver's PTX JIT.
std::string compile(const std::string& program, const std::string& key, int arch)
{
    const Api& a = api();
    int count = 0;
    std::vector<int> archs;
    if (a.nvrtcGetNumSupportedArchs(&count) == NVRTC_SUCCESS && count > 0)
    {
        archs.resize(static_cast<std::size_t>(count));
        a.nvrtcGetSupportedArchs(archs.data());
    }
    int target = arch;
    if (!archs.empty())
    {
        int best = archs.front();
        for (int x : archs)
            if (x <= arch)
                best = std::max(best, x);
        target = best;
    }
    nvrtcProgram prog = nullptr;
    if (a.nvrtcCreateProgram(&prog, program.c_str(), key.c_str(), 0, nullptr, nullptr) != NVRTC_SUCCESS)
        throw std::runtime_error("GPU: cannot create an NVRTC program");
    // No fast math, and no contraction into fused multiply-adds: the
    // compensated sums and the float-float phases mean what they say, as
    // Metal's safe math keeps them (docs/internals/gpu.md).
    const std::string arch_flag = "--gpu-architecture=compute_" + std::to_string(target);
    const char* options[] = {arch_flag.c_str(), "--std=c++14", "--fmad=false", "--ftz=false", "--prec-div=true",
            "--prec-sqrt=true"};
    const nvrtcResult result = a.nvrtcCompileProgram(prog, static_cast<int>(std::size(options)), options);
    std::string log;
    std::size_t size = 0;
    if (a.nvrtcGetProgramLogSize(prog, &size) == NVRTC_SUCCESS && size > 1)
    {
        log.resize(size);
        a.nvrtcGetProgramLog(prog, log.data());
        log.resize(size - 1);
    }
    if (result != NVRTC_SUCCESS)
    {
        a.nvrtcDestroyProgram(&prog);
        throw std::runtime_error("GPU: kernel sources " + key + " do not compile:\n" + log);
    }
    std::string ptx;
    a.nvrtcGetPTXSize(prog, &size);
    ptx.resize(size);
    a.nvrtcGetPTX(prog, ptx.data());
    a.nvrtcDestroyProgram(&prog);
    return ptx;
}

// The prelude, common.h and the named sources, in order.
std::string assemble(const std::vector<std::string>& sources, std::string& key)
{
    key.clear();
    for (const std::string& s : sources)
        key += s + ";";
    std::string program;
    for (const std::string& name : std::vector<std::string> {"prelude_cuda.h", "common.h"})
        program += std::string(kernel_source(name)) + "\n";
    for (const std::string& name : sources)
    {
        const char* text = kernel_source(name);
        if (!text)
            throw std::runtime_error("GPU: no kernel source " + name);
        program += std::string(text) + "\n";
    }
    return program;
}

} // namespace

struct Buffer::Impl
{
    CUcontext context = nullptr;
    CUdeviceptr pointer = 0;
    std::size_t size = 0;
    std::size_t allocated = 0;
    ~Impl()
    {
        if (pointer)
        {
            api().cuCtxSetCurrent(context);
            api().cuMemFree(pointer);
        }
    }
};

struct Kernel::Impl
{
    CUfunction function = nullptr;
    std::size_t max_threads = 0;
};

struct Device::Impl
{
    CUdevice device = 0;
    CUcontext context = nullptr;
    CUstream stream = nullptr;
    // For zeroing new buffers without waiting on the batches in flight.
    CUstream setup = nullptr;
    int arch = 0;
    std::string name;
    std::size_t total = 0;
    std::mutex lock;
    std::map<std::string, CUmodule> modules;
    std::map<std::string, std::shared_ptr<Kernel::Impl>> kernels;

    void bind() const
    {
        check(api().cuCtxSetCurrent(context), "cuCtxSetCurrent");
    }
};

struct Batch::Impl
{
    struct Launch
    {
        std::shared_ptr<Kernel::Impl> kernel;
        std::vector<std::vector<unsigned char>> args;
        unsigned groups = 0, threads = 0;
    };
    Device::Impl* device = nullptr;
    std::vector<Launch> launches;
    // Every buffer a dispatch names, kept alive until the batch is done.
    std::set<std::shared_ptr<Buffer::Impl>> buffers;
    CUevent start = nullptr, end = nullptr;
    bool committed = false, waited = false;

    // A batch dropped after its commit and before its wait() still runs on the
    // buffers it holds; it is waited for here, its error ignored.
    ~Impl()
    {
        if (!device)
            return;
        api().cuCtxSetCurrent(device->context);
        if (committed && !waited && end)
            api().cuEventSynchronize(end);
        if (start)
            api().cuEventDestroy(start);
        if (end)
            api().cuEventDestroy(end);
    }
};

std::size_t Buffer::size() const
{
    return impl_ ? impl_->size : 0;
}

void* Buffer::data() const
{
    return impl_ ? reinterpret_cast<void*>(impl_->pointer) : nullptr;
}

std::size_t Kernel::max_threads_per_group() const
{
    return impl_ ? impl_->max_threads : 0;
}

void Batch::dispatch(const Kernel& kernel, const std::vector<Arg>& args, std::size_t groups,
        std::size_t threads_per_group)
{
    if (impl_->committed)
        throw std::logic_error("gpu::Batch::dispatch after commit");
    if (groups == 0)
        return;
    if (threads_per_group == 0 || threads_per_group > kernel.max_threads_per_group())
        throw std::runtime_error("gpu::Batch::dispatch: " + std::to_string(threads_per_group)
                + " threads a group, where this kernel takes at most " + std::to_string(kernel.max_threads_per_group())
                + " on this device");
    if (groups > 0x7fffffffu)
        throw std::runtime_error("gpu::Batch::dispatch: " + std::to_string(groups) + " groups, past CUDA's grid");
    Impl::Launch launch;
    launch.kernel = kernel.impl_;
    launch.groups = static_cast<unsigned>(groups);
    launch.threads = static_cast<unsigned>(threads_per_group);
    launch.args.reserve(args.size());
    for (const Arg& a : args)
    {
        if (a.buffer)
        {
            if (!a.buffer->impl_)
                throw std::logic_error("gpu::Batch::dispatch: a buffer that was never allocated");
            const CUdeviceptr at = a.buffer->impl_->pointer + a.offset;
            std::vector<unsigned char> bytes(sizeof(at));
            std::memcpy(bytes.data(), &at, sizeof(at));
            launch.args.push_back(std::move(bytes));
            impl_->buffers.insert(a.buffer->impl_);
        }
        else
            launch.args.push_back(a.bytes);
    }
    impl_->launches.push_back(std::move(launch));
}

void Batch::commit()
{
    if (impl_->committed)
        return;
    Device::Impl& d = *impl_->device;
    d.bind();
    const Api& a = api();
    check(a.cuEventCreate(&impl_->start, CU_EVENT_DEFAULT), "cuEventCreate");
    check(a.cuEventCreate(&impl_->end, CU_EVENT_DEFAULT), "cuEventCreate");
    // The inputs the CPU wrote go over ahead of the kernels that read them,
    // rather than page by page as they fault.
    if (a.prefetch)
        for (const auto& b : impl_->buffers)
            a.prefetch(b->pointer, b->allocated, d.device, d.stream);
    check(a.cuEventRecord(impl_->start, d.stream), "cuEventRecord");
    std::vector<void*> params;
    for (Impl::Launch& launch : impl_->launches)
    {
        params.clear();
        for (auto& bytes : launch.args)
            params.push_back(bytes.data());
        check(a.cuLaunchKernel(launch.kernel->function, launch.groups, 1, 1, launch.threads, 1, 1, 0, d.stream,
                      params.data(), nullptr),
                "cuLaunchKernel");
    }
    check(a.cuEventRecord(impl_->end, d.stream), "cuEventRecord");
    impl_->committed = true;
}

void Batch::wait()
{
    commit();
    impl_->device->bind();
    const CUresult r = api().cuEventSynchronize(impl_->end);
    impl_->waited = true;
    if (r != CUDA_SUCCESS)
        throw std::runtime_error("GPU batch failed: " + cuda_error(r));
}

double Batch::gpu_seconds() const
{
    float ms = 0.0f;
    if (!impl_->waited || api().cuEventElapsedTime(&ms, impl_->start, impl_->end) != CUDA_SUCCESS)
        return 0.0;
    return 1e-3 * ms;
}

namespace {

// Why there is no device, or "" where there is one.
std::string& reason()
{
    static std::string r;
    return r;
}

} // namespace

Device* Device::get()
{
    static std::unique_ptr<Device> device = []() -> std::unique_ptr<Device> {
        const Api& a = api();
        if (!a.missing.empty())
        {
            reason() = a.missing;
            return nullptr;
        }
        auto none = [](const std::string& why) -> std::unique_ptr<Device> {
            reason() = why;
            return nullptr;
        };
        CUresult r = a.cuInit(0);
        if (r != CUDA_SUCCESS)
            return none("the CUDA driver did not start (" + cuda_error(r) + ")");
        int count = 0;
        if (a.cuDeviceGetCount(&count) != CUDA_SUCCESS || count == 0)
            return none("the CUDA driver sees no device");
        auto impl = std::make_shared<Impl>();
        if ((r = a.cuDeviceGet(&impl->device, 0)) != CUDA_SUCCESS)
            return none("cuDeviceGet: " + cuda_error(r));
        int managed = 0, concurrent = 0, major = 0, minor = 0;
        a.cuDeviceGetAttribute(&managed, CU_DEVICE_ATTRIBUTE_MANAGED_MEMORY, impl->device);
        a.cuDeviceGetAttribute(&concurrent, CU_DEVICE_ATTRIBUTE_CONCURRENT_MANAGED_ACCESS, impl->device);
        a.cuDeviceGetAttribute(&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, impl->device);
        a.cuDeviceGetAttribute(&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, impl->device);
        char name[256] = {};
        a.cuDeviceGetName(name, sizeof(name) - 1, impl->device);
        impl->name = name;
        impl->arch = 10 * major + minor;
        if (!managed || !concurrent)
            return none(impl->name + " has no concurrent managed memory, which the CUDA backend needs");
        std::size_t total = 0;
        a.cuDeviceTotalMem(&total, impl->device);
        impl->total = total;
        if ((r = a.cuDevicePrimaryCtxRetain(&impl->context, impl->device)) != CUDA_SUCCESS)
            return none("cuDevicePrimaryCtxRetain: " + cuda_error(r));
        if ((r = a.cuCtxSetCurrent(impl->context)) != CUDA_SUCCESS
                || (r = a.cuStreamCreate(&impl->stream, CU_STREAM_NON_BLOCKING)) != CUDA_SUCCESS
                || (r = a.cuStreamCreate(&impl->setup, CU_STREAM_NON_BLOCKING)) != CUDA_SUCCESS)
            return none("cannot set up the CUDA context: " + cuda_error(r));
        // The primary context and the streams are kept for the process.
        std::unique_ptr<Device> d(new Device());
        d->impl_ = impl;
        return d;
    }();
    return device.get();
}

const char* Device::backend()
{
    return "cuda";
}

std::string Device::why_none()
{
    get();
    return reason();
}

std::string Device::name() const
{
    return impl_->name;
}

std::size_t Device::max_buffer_bytes() const
{
    return impl_->total;
}

std::size_t Device::working_set_bytes() const
{
    // What is free now, less a margin for the driver and the kernels' own
    // memory; managed memory could oversubscribe, but a march that pages
    // would be slower than the CPU's.
    impl_->bind();
    std::size_t free = 0, total = 0;
    check(api().cuMemGetInfo(&free, &total), "cuMemGetInfo");
    const std::size_t margin = std::min<std::size_t>(free / 10, std::size_t {512} << 20);
    return free - margin;
}

Buffer Device::buffer(std::size_t bytes)
{
    impl_->bind();
    Buffer b;
    b.impl_ = std::make_shared<Buffer::Impl>();
    b.impl_->context = impl_->context;
    b.impl_->size = bytes;
    // CUDA will not make an empty allocation; a word stands in.
    b.impl_->allocated = std::max<std::size_t>(bytes, 4);
    const CUresult r = api().cuMemAllocManaged(&b.impl_->pointer, b.impl_->allocated, CU_MEM_ATTACH_GLOBAL);
    if (r != CUDA_SUCCESS)
    {
        b.impl_->pointer = 0;
        throw std::runtime_error("GPU: cannot allocate a buffer of " + std::to_string(bytes) + " bytes ("
                + cuda_error(r) + ")");
    }
    // Zeroed on the device, so that a field the GPU alone uses is never
    // paged through the CPU; the stream is not the batches', which this does
    // not wait for.
    check(api().cuMemsetD8Async(b.impl_->pointer, 0, b.impl_->allocated, impl_->setup), "cuMemsetD8Async");
    check(api().cuStreamSynchronize(impl_->setup), "cuStreamSynchronize");
    return b;
}

Kernel Device::kernel(const std::vector<std::string>& sources, const std::string& entry)
{
    std::string key;
    const std::string program = assemble(sources, key);
    std::lock_guard<std::mutex> guard(impl_->lock);
    Kernel k;
    if (auto it = impl_->kernels.find(key + entry); it != impl_->kernels.end())
    {
        k.impl_ = it->second;
        return k;
    }
    impl_->bind();
    const Api& a = api();
    CUmodule& module = impl_->modules[key];
    if (!module)
    {
        const std::string ptx = compile(program, key, impl_->arch);
        const CUresult r = a.cuModuleLoadDataEx(&module, ptx.c_str(), 0, nullptr, nullptr);
        if (r != CUDA_SUCCESS)
        {
            module = nullptr;
            throw std::runtime_error("GPU: the driver cannot load kernel sources " + key + ": " + cuda_error(r));
        }
    }
    auto impl = std::make_shared<Kernel::Impl>();
    if (a.cuModuleGetFunction(&impl->function, module, entry.c_str()) != CUDA_SUCCESS)
        throw std::runtime_error("GPU: no kernel " + entry + " in " + key);
    int threads = 0;
    check(a.cuFuncGetAttribute(&threads, CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK, impl->function),
            "cuFuncGetAttribute");
    impl->max_threads = static_cast<std::size_t>(threads);
    impl_->kernels[key + entry] = impl;
    k.impl_ = impl;
    return k;
}

Batch Device::batch()
{
    Batch b;
    b.impl_ = std::make_shared<Batch::Impl>();
    b.impl_->device = impl_.get();
    return b;
}

bool compile_only(const std::vector<std::string>& sources)
{
    const Api& a = api();
    if (!a.nvrtcCompileProgram)
        throw std::runtime_error("GPU: " + (a.missing.empty() ? std::string("no NVRTC") : a.missing));
    std::string key;
    const std::string program = assemble(sources, key);
    // The oldest architecture with concurrent managed memory: what any device
    // this backend takes can run.
    compile(program, key, 60);
    return true;
}

} // namespace niceops::gpu
