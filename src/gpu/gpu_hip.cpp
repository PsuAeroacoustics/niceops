// The HIP backend of gpu.hpp, for AMD GPUs.  The kernels are compiled from
// source at run time by hiprtc, through prelude_cuda.h unchanged, to a code
// object for the device's own target; libamdhip64 and hiprtc are opened with
// dlopen when the device is first asked for, so a build with this backend links
// neither and runs, on the CPU, where ROCm is not installed.
//
// gpu.hpp's buffers are shared with the CPU.  Here each is device memory and a
// pinned copy in host memory, which data() is; managed memory, which CUDA's
// backend uses, moved a batch's inputs and results between the two a page at
// a time, and made a --buildings run slower on an MI210 than on the CPU.  The
// copies follow the CPU's own use, which goes through data(): what it touched
// goes over at the next commit, and what a batch writes comes back as the
// batch ends if the CPU uses it, or when data() asks for it.  So the CPU must
// not write a buffer while a batch that names it is in flight, nor keep
// data()'s pointer from one batch to the next, which no caller does
// (docs/internals/gpu.md).  A batch records its dispatches and launches them on
// the device's one stream at commit(), so batches run in commit order as a
// Metal queue's command buffers do.

#include "gpu/gpu.hpp"

#include <dlfcn.h>
#include <hip/hip_runtime_api.h>
#include <sys/mman.h>
#include <hip/hiprtc.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <thread>
#include <unistd.h>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>

namespace niceops::gpu {

namespace {

// The entry points used, looked up by the name hip_runtime_api.h maps each to
// (hipGetDevicePropertiesR0600, say, whose struct is the headers').
#define NK_STR(x) #x
#define NK_XSTR(x) NK_STR(x)
#define NK_HIP_FUNCTIONS(X)                                                                                          \
    X(hipInit) X(hipGetErrorString) X(hipGetDeviceCount) X(hipGetDeviceProperties) X(hipSetDevice) X(hipMemGetInfo)  \
    X(hipFree) X(hipHostFree) X(hipHostRegister) X(hipHostUnregister) X(hipMemcpyAsync) X(hipMemsetAsync) X(hipStreamSynchronize) X(hipStreamCreateWithFlags) X(hipEventCreate) X(hipEventRecord)      \
    X(hipEventSynchronize) X(hipEventElapsedTime) X(hipEventDestroy) X(hipModuleLoadData) X(hipModuleGetFunction)    \
    X(hipFuncGetAttribute) X(hipModuleLaunchKernel)
#define NK_HIPRTC_FUNCTIONS(X)                                                                                       \
    X(hiprtcCreateProgram) X(hiprtcCompileProgram) X(hiprtcGetProgramLogSize) X(hiprtcGetProgramLog)                 \
    X(hiprtcGetCodeSize) X(hiprtcGetCode) X(hiprtcDestroyProgram) X(hiprtcVersion)

// Each function's pointer type, named outside Api: a member may not take its
// type from the function it is named after.
#define NK_TYPE(name) using name##_type = decltype(&name);
NK_HIP_FUNCTIONS(NK_TYPE)
NK_HIPRTC_FUNCTIONS(NK_TYPE)
#undef NK_TYPE

struct Api
{
#define NK_MEMBER(name) name##_type name = nullptr;
    NK_HIP_FUNCTIONS(NK_MEMBER)
    NK_HIPRTC_FUNCTIONS(NK_MEMBER)
#undef NK_MEMBER
    // Named by hand: the header overloads them for typed pointers.
    hipError_t (*hipMalloc)(void**, size_t) = nullptr;
    hipError_t (*hipHostMalloc)(void**, size_t, unsigned) = nullptr;
    // Why the API is not there, or "" where it is.
    std::string missing;
};

// Each name as the loader finds it, then in the ROCm the build was configured
// against (NICEOPS_ROCM_LIB_DIR), for a ROCm that is not on the library path.
void* open_first(std::initializer_list<std::string> names)
{
    for (const std::string& name : names)
        if (void* h = dlopen(name.c_str(), RTLD_NOW | RTLD_LOCAL))
            return h;
#ifdef NICEOPS_ROCM_LIB_DIR
    for (const std::string& name : names)
        if (void* h = dlopen((std::string(NICEOPS_ROCM_LIB_DIR) + "/" + name).c_str(), RTLD_NOW | RTLD_LOCAL))
            return h;
#endif
    return nullptr;
}

const Api& api()
{
    static const Api loaded = [] {
        Api a;
        // The libraries of the headers this was built with first: the HIP
        // ABI is stable within a major version.
        const std::string major = std::to_string(HIP_VERSION_MAJOR);
        void* hip = open_first({"libamdhip64.so." + major, "libamdhip64.so"});
        if (!hip)
        {
            a.missing = "the HIP runtime (libamdhip64.so." + major + ") is not installed";
            return a;
        }
        void* hiprtc = open_first({"libhiprtc.so." + major, "libhiprtc.so"});
        if (!hiprtc)
        {
            a.missing = "the HIP run-time compiler (libhiprtc.so." + major + ") is not on the library path";
            return a;
        }
#define NK_LOAD_FROM(library)                                                                                        \
    [&](auto& slot, const char* symbol) {                                                                            \
        slot = reinterpret_cast<std::remove_reference_t<decltype(slot)>>(dlsym(library, symbol));                     \
        if (!slot && a.missing.empty())                                                                              \
            a.missing = std::string("no ") + symbol + " in the HIP libraries";                                        \
    }
        auto from_hip = NK_LOAD_FROM(hip);
        auto from_hiprtc = NK_LOAD_FROM(hiprtc);
#undef NK_LOAD_FROM
#define NK_LOAD_HIP(name) from_hip(a.name, NK_XSTR(name));
#define NK_LOAD_HIPRTC(name) from_hiprtc(a.name, NK_XSTR(name));
        NK_HIP_FUNCTIONS(NK_LOAD_HIP)
        NK_HIPRTC_FUNCTIONS(NK_LOAD_HIPRTC)
#undef NK_LOAD_HIP
#undef NK_LOAD_HIPRTC
        from_hip(a.hipMalloc, "hipMalloc");
        from_hip(a.hipHostMalloc, "hipHostMalloc");
        return a;
    }();
    return loaded;
}

std::string hip_error(hipError_t r)
{
    if (api().hipGetErrorString)
        if (const char* text = api().hipGetErrorString(r))
            return text;
    return "HIP error " + std::to_string(static_cast<int>(r));
}

void check(hipError_t r, const char* what)
{
    if (r != hipSuccess)
        throw std::runtime_error(std::string("GPU: ") + what + ": " + hip_error(r));
}

// hiprtc's options for target `arch` (a target ID, as the device names it:
// gfx90a:sramecc+:xnack-).  No fast math, no contraction into fused
// multiply-adds (hiprtc's default contracts), denormals kept and division and
// square roots correctly rounded: the compensated sums and the float-float
// phases mean what they say, as Metal's safe math keeps them
// (docs/internals/gpu.md).
std::vector<std::string> compile_options(const std::string& arch)
{
    return {"--offload-arch=" + arch, "-std=c++14", "-O3", "-ffp-contract=off", "-fno-fast-math",
            "-fno-gpu-flush-denormals-to-zero", "-fhip-fp32-correctly-rounded-divide-sqrt"};
}

// The program's code object for target `arch`, or a std::runtime_error
// carrying hiprtc's log.
std::vector<char> compile(const std::string& program, const std::string& key, const std::string& arch)
{
    const Api& a = api();
    hiprtcProgram prog = nullptr;
    if (a.hiprtcCreateProgram(&prog, program.c_str(), key.c_str(), 0, nullptr, nullptr) != HIPRTC_SUCCESS)
        throw std::runtime_error("GPU: cannot create a hiprtc program");
    const std::vector<std::string> strings = compile_options(arch);
    std::vector<const char*> options;
    for (const std::string& o : strings)
        options.push_back(o.c_str());
    const hiprtcResult result = a.hiprtcCompileProgram(prog, static_cast<int>(options.size()), options.data());
    std::string log;
    std::size_t size = 0;
    if (a.hiprtcGetProgramLogSize(prog, &size) == HIPRTC_SUCCESS && size > 1)
    {
        log.resize(size);
        (void)a.hiprtcGetProgramLog(prog, log.data());
        log.resize(size - 1);
    }
    if (result != HIPRTC_SUCCESS)
    {
        (void)a.hiprtcDestroyProgram(&prog);
        throw std::runtime_error("GPU: kernel sources " + key + " do not compile:\n" + log);
    }
    std::vector<char> code;
    (void)a.hiprtcGetCodeSize(prog, &size);
    code.resize(size);
    (void)a.hiprtcGetCode(prog, code.data());
    (void)a.hiprtcDestroyProgram(&prog);
    return code;
}

// The directory of compiled code objects: NICEOPS_KERNEL_CACHE, else
// $XDG_CACHE_HOME/niceops/kernels, else ~/.cache/niceops/kernels; empty for
// none (NICEOPS_KERNEL_CACHE=off, or no home).
std::filesystem::path cache_directory()
{
    if (const char* d = std::getenv("NICEOPS_KERNEL_CACHE"))
        return std::string(d) == "off" ? std::filesystem::path() : std::filesystem::path(d);
    if (const char* x = std::getenv("XDG_CACHE_HOME"); x && *x)
        return std::filesystem::path(x) / "niceops" / "kernels";
    if (const char* h = std::getenv("HOME"); h && *h)
        return std::filesystem::path(h) / ".cache" / "niceops" / "kernels";
    return {};
}

// The program's code object, from the cache where an entry was made by this
// hiprtc for the same program, options and target, else compiled and
// entered.  hiprtc took 1.2 s over point_rays.kern on every run, and ROCm's own
// cache (AMD_COMGR_CACHE) saved none of it.  An entry holds its whole key, so
// a load is checked against it byte for byte, not by its hash; a cache that
// cannot be read or written only costs the compile.
std::vector<char> cached_compile(const std::string& program, const std::string& key, const std::string& arch)
{
    int major = 0, minor = 0;
    (void)api().hiprtcVersion(&major, &minor);
    std::string full = "hiprtc " + std::to_string(major) + "." + std::to_string(minor) + "\n";
    for (const std::string& o : compile_options(arch))
        full += o + "\n";
    full += program;
    // FNV-1a, for the file's name only.
    std::uint64_t hash = 1469598103934665603ull;
    for (unsigned char c : full)
        hash = (hash ^ c) * 1099511628211ull;
    const std::filesystem::path directory = cache_directory();
    std::filesystem::path file;
    if (!directory.empty())
    {
        std::ostringstream name;
        name << std::hex << hash << ".hsaco";
        file = directory / name.str();
        std::ifstream in(file, std::ios::binary);
        std::uint64_t length = 0;
        if (in.read(reinterpret_cast<char*>(&length), sizeof(length)) && length == full.size())
        {
            std::string stored(length, '\0');
            if (in.read(stored.data(), static_cast<std::streamsize>(length)) && stored == full)
            {
                std::vector<char> code((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
                if (!code.empty())
                    return code;
            }
        }
    }
    std::vector<char> code = compile(program, key, arch);
    if (!file.empty())
    {
        std::error_code error;
        std::filesystem::create_directories(directory, error);
        std::ostringstream temporary;
        temporary << file.string() << ".tmp." << getpid() << "." << std::hash<std::thread::id>()(std::this_thread::get_id());
        {
            std::ofstream out(temporary.str(), std::ios::binary);
            const std::uint64_t length = full.size();
            out.write(reinterpret_cast<const char*>(&length), sizeof(length));
            out.write(full.data(), static_cast<std::streamsize>(full.size()));
            out.write(code.data(), static_cast<std::streamsize>(code.size()));
            if (!out)
                error = std::make_error_code(std::errc::io_error);
        }
        if (!error)
            std::filesystem::rename(temporary.str(), file, error);
        if (error)
            std::filesystem::remove(temporary.str(), error);
    }
    return code;
}

// By slot, whether kernel `entry` of `program` writes the buffer there: its
// signature's NK_OUT arguments, the dialect's only writable pointers (an
// NK_IN is const).  Empty where the signature is not found or not in the
// dialect's form, which commit() takes as every buffer written.
std::vector<bool> written_slots(const std::string& program, const std::string& entry)
{
    const std::string head = "NK_KERNEL(" + entry + ")(";
    const std::size_t at = program.find(head);
    if (at == std::string::npos)
        return {};
    std::vector<std::string> args;
    std::string arg;
    int depth = 0;
    for (std::size_t i = at + head.size(); i < program.size(); i++)
    {
        const char ch = program[i];
        if (ch == '(')
            depth++;
        else if (ch == ')' && depth-- == 0)
        {
            args.push_back(arg);
            break;
        }
        if (ch == ',' && depth == 0)
        {
            args.push_back(arg);
            arg.clear();
        }
        else
            arg += ch;
    }
    std::vector<bool> writes;
    for (std::string a : args)
    {
        a.erase(0, a.find_first_not_of(" \t\n"));
        if (a.rfind("NK_GROUP_ARGS", 0) == 0)
            break;
        const bool in = a.rfind("NK_IN(", 0) == 0, out = a.rfind("NK_OUT(", 0) == 0;
        if (!in && !out && a.rfind("NK_VALUE(", 0) != 0)
            return {};
        // The slot is the argument's last, and must be its position: binding
        // is by position (prelude_cuda.h).
        const std::size_t comma = a.rfind(','), close = a.rfind(')');
        if (comma == std::string::npos || close == std::string::npos || close < comma)
            return {};
        try
        {
            if (std::stoul(a.substr(comma + 1, close - comma - 1)) != writes.size())
                return {};
        }
        catch (const std::exception&)
        {
            return {};
        }
        writes.push_back(out);
        // The last argument before NK_GROUP_ARGS carries no comma after it.
        if (a.find("NK_GROUP_ARGS") != std::string::npos)
            break;
    }
    return writes;
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

// A buffer's memory: on the device, and its pinned copy on the host.
struct Block
{
    void* on_device = nullptr;
    void* host = nullptr;
    std::size_t capacity = 0;
    // The host's mapping, where it is this backend's own (host_memory), else
    // null and `host` is hipHostMalloc's.
    void* mapping = nullptr;
    std::size_t mapped = 0;
};

// Host memory for a block of `capacity` bytes, pinned for the copies, zeroed.
// Mapped here on 2 MB boundaries and advised into transparent huge pages,
// then registered: hipHostMalloc pins 4 KB pages, and took 213 ms for 256 MB
// with its zeroing, and 62 ms to free, where this took 35 and 3.  Where the
// mapping or its registration fails (a memlock limit, say), hipHostMalloc's.
hipError_t host_memory(Block& b)
{
    const std::size_t huge = std::size_t {2} << 20;
    const std::size_t length = (b.capacity + huge - 1) / huge * huge;
    void* map = mmap(nullptr, length + huge, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (map != MAP_FAILED)
    {
        // Trim the mapping to a 2 MB boundary at each end.
        auto* start = static_cast<unsigned char*>(map);
        auto* aligned = reinterpret_cast<unsigned char*>(
                (reinterpret_cast<std::uintptr_t>(start) + huge - 1) / huge * huge);
        if (aligned > start)
            munmap(start, static_cast<std::size_t>(aligned - start));
        const std::size_t tail = static_cast<std::size_t>(start + length + huge - (aligned + length));
        if (tail > 0)
            munmap(aligned + length, tail);
        (void)madvise(aligned, length, MADV_HUGEPAGE);
        // Touched before it is pinned, so that the pages are the huge ones.
        std::memset(aligned, 0, length);
        if (api().hipHostRegister(aligned, length, hipHostRegisterDefault) == hipSuccess)
        {
            b.host = b.mapping = aligned;
            b.mapped = length;
            return hipSuccess;
        }
        munmap(aligned, length);
    }
    const hipError_t r = api().hipHostMalloc(&b.host, b.capacity, hipHostMallocDefault);
    if (r == hipSuccess)
        std::memset(b.host, 0, b.capacity);
    return r;
}

void release_block(int device, const Block& b)
{
    (void)api().hipSetDevice(device);
    if (b.on_device)
        (void)api().hipFree(b.on_device);
    if (b.mapping)
    {
        (void)api().hipHostUnregister(b.mapping);
        munmap(b.mapping, b.mapped);
    }
    else if (b.host)
        (void)api().hipHostFree(b.host);
}

// Blocks of freed buffers, by capacity, for the next buffer of that size.
// Pinning host memory costs milliseconds an allocation, and the callers make
// buffers by the batch (a --buildings run's every block of states): unpooled,
// allocating took three times the kernels' time.  A block's capacity is a
// power of two up to 64 MB and an eighth of one past it (at most an eighth
// wasted, so that the FDTD's fields still fit the device as it judged them),
// and a buffer takes the smallest idle block that holds it and is not a
// quarter larger.  At most `idle` bytes of blocks are kept.
//
// With the pool, what every buffer shares: the lock that orders the copies
// with the batches, and the stream data() reads a buffer back on.
struct Memory
{
    static constexpr std::size_t idle = std::size_t {1} << 30;
    int device = 0;
    std::mutex flight;
    hipStream_t readback = nullptr;
    std::mutex lock;
    std::multimap<std::size_t, Block> free;
    std::size_t held = 0;

    static std::size_t capacity_for(std::size_t bytes)
    {
        std::size_t c = 4096;
        while (c < bytes)
            c *= 2;
        if (c <= (std::size_t {64} << 20))
            return c;
        const std::size_t step = c / 16;
        return (bytes + step - 1) / step * step;
    }
    bool take(std::size_t bytes, Block& b)
    {
        std::lock_guard<std::mutex> guard(lock);
        auto it = free.lower_bound(bytes);
        if (it == free.end() || it->first > bytes + bytes / 4)
            return false;
        b = it->second;
        held -= it->first;
        free.erase(it);
        return true;
    }
    void give(const Block& b)
    {
        {
            std::lock_guard<std::mutex> guard(lock);
            if (held + b.capacity <= idle)
            {
                free.emplace(b.capacity, b);
                held += b.capacity;
                return;
            }
        }
        release_block(device, b);
    }
};

} // namespace

// Under Memory::flight, where a buffer's contents are.  data() is the CPU's
// only way to them, so it is what moves them: it brings the host's copy up
// to date where a finished batch wrote the device's, and marks it touched
// (`host_dirty`), which the next commit copies over, and the CPU's
// (`cpu_visible`), which every batch that writes it then copies back as it
// ends.  A buffer the CPU never asks for -- a kernel's scratch, the FDTD's
// fields -- never crosses: copied whole either way, the point model's
// amplitudes and the FDTD's fields made most of a batch's traffic.
struct Buffer::Impl
{
    std::shared_ptr<Memory> memory;
    Block block;
    std::size_t size = 0;
    // Batches committed and not yet waited for that name this buffer: while
    // any is, the device's copy is ahead of the host's.
    int in_flight = 0;
    bool host_dirty = false, cpu_visible = false, device_newer = false;
    ~Impl()
    {
        if (block.on_device || block.host)
            memory->give(block);
    }
};

struct Kernel::Impl
{
    hipFunction_t function = nullptr;
    std::size_t max_threads = 0;
    // By slot, whether the kernel may write the buffer there (NK_OUT), from
    // its signature; empty where that could not be read, and then every
    // buffer is taken as written.
    std::vector<bool> writes;
};

struct Device::Impl
{
    int device = 0;
    hipStream_t stream = nullptr;
    std::string arch;
    std::string name;
    std::size_t total = 0;
    std::shared_ptr<Memory> memory = std::make_shared<Memory>();
    std::mutex lock;
    std::map<std::string, hipModule_t> modules;
    std::map<std::string, std::shared_ptr<Kernel::Impl>> kernels;

    // HIP's current device is a thread's: the GPU thread of a BandStream is
    // not the one that made the device.
    void bind() const
    {
        check(api().hipSetDevice(device), "hipSetDevice");
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
    // Every buffer a dispatch names, kept alive until the batch is done;
    // those a kernel may write; and of those, the ones copied back as the
    // batch ends (the CPU's, at commit).
    std::set<std::shared_ptr<Buffer::Impl>> buffers;
    std::set<Buffer::Impl*> written;
    std::set<Buffer::Impl*> back;
    // Around the kernels, and after the copies back.
    hipEvent_t start = nullptr, end = nullptr, copied = nullptr;
    bool committed = false, waited = false;

    // Once the batch is done: what it wrote and did not copy back is newer on
    // the device.
    void release()
    {
        std::lock_guard<std::mutex> guard(device->memory->flight);
        for (const auto& b : buffers)
            b->in_flight--;
        for (Buffer::Impl* b : written)
            b->device_newer = back.count(b) == 0;
    }

    // A batch dropped after its commit and before its wait() still runs on the
    // buffers it holds; it is waited for here, its error ignored.
    ~Impl()
    {
        if (!device)
            return;
        (void)api().hipSetDevice(device->device);
        if (committed && !waited)
        {
            if (copied)
                (void)api().hipEventSynchronize(copied);
            release();
        }
        for (hipEvent_t e : {start, end, copied})
            if (e)
                (void)api().hipEventDestroy(e);
    }
};

std::size_t Buffer::size() const
{
    return impl_ ? impl_->size : 0;
}

void* Buffer::data() const
{
    if (!impl_)
        return nullptr;
    Memory& m = *impl_->memory;
    std::lock_guard<std::mutex> guard(m.flight);
    if (impl_->device_newer && impl_->size > 0)
    {
        check(api().hipSetDevice(m.device), "hipSetDevice");
        check(api().hipMemcpyAsync(impl_->block.host, impl_->block.on_device, impl_->size, hipMemcpyDeviceToHost,
                      m.readback),
                "hipMemcpyAsync");
        check(api().hipStreamSynchronize(m.readback), "hipStreamSynchronize");
    }
    impl_->device_newer = false;
    impl_->host_dirty = impl_->cpu_visible = true;
    return impl_->block.host;
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
    // An AMD dispatch counts its grid in threads, in 32 bits.
    if (groups * threads_per_group > 0xffffffffu)
        throw std::runtime_error("gpu::Batch::dispatch: " + std::to_string(groups) + " groups of "
                + std::to_string(threads_per_group) + ", past HIP's grid of 2^32 threads");
    Impl::Launch launch;
    launch.kernel = kernel.impl_;
    launch.groups = static_cast<unsigned>(groups);
    launch.threads = static_cast<unsigned>(threads_per_group);
    launch.args.reserve(args.size());
    const std::vector<bool>& writes = kernel.impl_->writes;
    for (std::size_t slot = 0; slot < args.size(); slot++)
    {
        const Arg& a = args[slot];
        if (a.buffer)
        {
            if (!a.buffer->impl_)
                throw std::logic_error("gpu::Batch::dispatch: a buffer that was never allocated");
            void* at = static_cast<unsigned char*>(a.buffer->impl_->block.on_device) + a.offset;
            std::vector<unsigned char> bytes(sizeof(at));
            std::memcpy(bytes.data(), &at, sizeof(at));
            launch.args.push_back(std::move(bytes));
            impl_->buffers.insert(a.buffer->impl_);
            if (writes.empty() || slot >= writes.size() || writes[slot])
                impl_->written.insert(a.buffer->impl_.get());
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
    for (hipEvent_t* e : {&impl_->start, &impl_->end, &impl_->copied})
        check(a.hipEventCreate(e), "hipEventCreate");
    std::lock_guard<std::mutex> guard(d.memory->flight);
    // What the CPU touched goes over, unless a batch in flight names it or a
    // finished one wrote it since: then the device's copy is the newer, and
    // the CPU, which may not write a buffer in flight, only read it.  The
    // counts are taken before anything can throw, so that the destructor's
    // release always matches them.
    std::vector<Buffer::Impl*> over;
    for (const auto& b : impl_->buffers)
    {
        if (b->in_flight++ == 0 && b->host_dirty && !b->device_newer && b->size > 0)
            over.push_back(b.get());
        if (b->in_flight == 1 && b->device_newer)
            b->host_dirty = false;
    }
    for (Buffer::Impl* b : over)
        b->host_dirty = false;
    for (Buffer::Impl* b : impl_->written)
        if (b->cpu_visible)
            impl_->back.insert(b);
    impl_->committed = true;
    for (Buffer::Impl* b : over)
        check(a.hipMemcpyAsync(b->block.on_device, b->block.host, b->size, hipMemcpyHostToDevice, d.stream),
                "hipMemcpyAsync");
    check(a.hipEventRecord(impl_->start, d.stream), "hipEventRecord");
    std::vector<void*> params;
    for (Impl::Launch& launch : impl_->launches)
    {
        params.clear();
        for (auto& bytes : launch.args)
            params.push_back(bytes.data());
        check(a.hipModuleLaunchKernel(launch.kernel->function, launch.groups, 1, 1, launch.threads, 1, 1, 0, d.stream,
                      params.data(), nullptr),
                "hipModuleLaunchKernel");
    }
    check(a.hipEventRecord(impl_->end, d.stream), "hipEventRecord");
    // What the kernels may have written to the CPU's buffers comes back before
    // the next batch's kernels touch it, so that wait() leaves the CPU this
    // batch's results; a later batch in flight copies back what it changes in
    // its turn.
    for (Buffer::Impl* b : impl_->back)
        if (b->size > 0)
            check(a.hipMemcpyAsync(b->block.host, b->block.on_device, b->size, hipMemcpyDeviceToHost, d.stream),
                    "hipMemcpyAsync");
    check(a.hipEventRecord(impl_->copied, d.stream), "hipEventRecord");
}

void Batch::wait()
{
    commit();
    if (impl_->waited)
        return;
    impl_->device->bind();
    const hipError_t r = api().hipEventSynchronize(impl_->copied);
    impl_->waited = true;
    impl_->release();
    if (r != hipSuccess)
        throw std::runtime_error("GPU batch failed: " + hip_error(r));
}

double Batch::gpu_seconds() const
{
    float ms = 0.0f;
    if (!impl_->waited || api().hipEventElapsedTime(&ms, impl_->start, impl_->end) != hipSuccess)
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
        hipError_t r = a.hipInit(0);
        if (r != hipSuccess)
            return none("the HIP runtime did not start (" + hip_error(r) + ")");
        int count = 0;
        if (a.hipGetDeviceCount(&count) != hipSuccess || count == 0)
            return none("the HIP runtime sees no device");
        auto impl = std::make_shared<Impl>();
        hipDeviceProp_t properties {};
        if ((r = a.hipGetDeviceProperties(&properties, impl->device)) != hipSuccess)
            return none("hipGetDeviceProperties: " + hip_error(r));
        impl->name = properties.name;
        impl->arch = properties.gcnArchName;
        impl->total = properties.totalGlobalMem;
        impl->memory->device = impl->device;
        if ((r = a.hipSetDevice(impl->device)) != hipSuccess
                || (r = a.hipStreamCreateWithFlags(&impl->stream, hipStreamNonBlocking)) != hipSuccess
                || (r = a.hipStreamCreateWithFlags(&impl->memory->readback, hipStreamNonBlocking)) != hipSuccess)
            return none("cannot set up the HIP device: " + hip_error(r));
        // The streams are kept for the process.
        std::unique_ptr<Device> d(new Device());
        d->impl_ = impl;
        return d;
    }();
    return device.get();
}

const char* Device::backend()
{
    return "hip";
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
    // What is free now, less a margin for the runtime and the kernels' own
    // memory.
    impl_->bind();
    std::size_t free = 0, total = 0;
    check(api().hipMemGetInfo(&free, &total), "hipMemGetInfo");
    const std::size_t margin = std::min<std::size_t>(free / 10, std::size_t {512} << 20);
    return free - margin;
}

Buffer Device::buffer(std::size_t bytes)
{
    impl_->bind();
    Buffer b;
    b.impl_ = std::make_shared<Buffer::Impl>();
    b.impl_->memory = impl_->memory;
    b.impl_->size = bytes;
    // HIP will not make an empty allocation; a pooled block stands in.
    Block& block = b.impl_->block;
    const std::size_t needed = std::max<std::size_t>(bytes, 4);
    Memory& m = *impl_->memory;
    if (!m.take(needed, block))
    {
        const std::size_t capacity = Memory::capacity_for(needed);
        block.capacity = capacity;
        hipError_t r = api().hipMalloc(&block.on_device, capacity);
        if (r == hipSuccess)
            r = host_memory(block);
        if (r != hipSuccess)
        {
            release_block(impl_->device, block);
            block = Block {};
            throw std::runtime_error("GPU: cannot allocate a buffer of " + std::to_string(bytes) + " bytes ("
                    + hip_error(r) + ")");
        }
    }
    // Zeroed on both sides: a buffer only the GPU uses is never copied over.
    std::memset(block.host, 0, bytes);
    check(api().hipMemsetAsync(block.on_device, 0, needed, m.readback), "hipMemsetAsync");
    check(api().hipStreamSynchronize(m.readback), "hipStreamSynchronize");
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
    hipModule_t& module = impl_->modules[key];
    if (!module)
    {
        std::vector<char> code = cached_compile(program, key, impl_->arch);
        hipError_t r = a.hipModuleLoadData(&module, code.data());
        // An entry the runtime will not load is compiled afresh.
        if (r != hipSuccess)
        {
            code = compile(program, key, impl_->arch);
            r = a.hipModuleLoadData(&module, code.data());
        }
        if (r != hipSuccess)
        {
            module = nullptr;
            throw std::runtime_error("GPU: the runtime cannot load kernel sources " + key + ": " + hip_error(r));
        }
    }
    auto impl = std::make_shared<Kernel::Impl>();
    if (a.hipModuleGetFunction(&impl->function, module, entry.c_str()) != hipSuccess)
        throw std::runtime_error("GPU: no kernel " + entry + " in " + key);
    int threads = 0;
    check(a.hipFuncGetAttribute(&threads, HIP_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK, impl->function),
            "hipFuncGetAttribute");
    impl->max_threads = static_cast<std::size_t>(threads);
    impl->writes = written_slots(program, entry);
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
    if (!a.hiprtcCompileProgram)
        throw std::runtime_error("GPU: " + (a.missing.empty() ? std::string("no hiprtc") : a.missing));
    std::string key;
    const std::string program = assemble(sources, key);
    // The MI200's target: hiprtc needs no device for it.
    compile(program, key, "gfx90a");
    return true;
}

} // namespace niceops::gpu
