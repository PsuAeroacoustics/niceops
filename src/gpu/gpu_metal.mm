// The Metal backend of gpu.hpp.  Kernels are compiled from source at run time
// (newLibraryWithSource), which needs only the Metal framework, not Xcode.

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "gpu/gpu.hpp"

#include <map>
#include <mutex>
#include <stdexcept>

namespace niceops::gpu {

struct Buffer::Impl
{
    id<MTLBuffer> buffer = nil;
    std::size_t size = 0;
};

struct Kernel::Impl
{
    id<MTLComputePipelineState> pipeline = nil;
};

struct Batch::Impl
{
    id<MTLCommandBuffer> commands = nil;
    id<MTLComputeCommandEncoder> encoder = nil;
    bool committed = false;
};

struct Device::Impl
{
    id<MTLDevice> device = nil;
    id<MTLCommandQueue> queue = nil;
    std::mutex lock;
    std::map<std::string, id<MTLLibrary>> libraries;
    std::map<std::string, std::shared_ptr<Kernel::Impl>> kernels;
};

std::size_t Buffer::size() const
{
    return impl_ ? impl_->size : 0;
}

void* Buffer::data() const
{
    return impl_ ? impl_->buffer.contents : nullptr;
}

std::size_t Kernel::max_threads_per_group() const
{
    return impl_ ? impl_->pipeline.maxTotalThreadsPerThreadgroup : 0;
}

void Batch::dispatch(const Kernel& kernel, const std::vector<Arg>& args, std::size_t groups,
        std::size_t threads_per_group)
{
    if (impl_->committed)
        throw std::logic_error("gpu::Batch::dispatch after commit");
    if (groups == 0)
        return;
    // The kernels' group sums and tiles assume the size they are given, so a
    // pipeline that cannot run it is an error here, not garbage later.
    if (threads_per_group == 0 || threads_per_group > kernel.max_threads_per_group())
        throw std::runtime_error("gpu::Batch::dispatch: " + std::to_string(threads_per_group)
                + " threads a group, where this kernel takes at most " + std::to_string(kernel.max_threads_per_group())
                + " on this device");
    for (const Arg& a : args)
        if (a.buffer && !a.buffer->impl_)
            throw std::logic_error("gpu::Batch::dispatch: a buffer that was never allocated");
    @autoreleasepool
    {
        if (!impl_->encoder)
            impl_->encoder = [impl_->commands computeCommandEncoder];
        id<MTLComputeCommandEncoder> encoder = impl_->encoder;
        [encoder setComputePipelineState:kernel.impl_->pipeline];
        for (std::size_t i = 0; i < args.size(); i++)
        {
            const Arg& a = args[i];
            if (a.buffer)
                [encoder setBuffer:a.buffer->impl_->buffer offset:a.offset atIndex:i];
            else
                [encoder setBytes:a.bytes.data() length:a.bytes.size() atIndex:i];
        }
        [encoder dispatchThreadgroups:MTLSizeMake(groups, 1, 1) threadsPerThreadgroup:MTLSizeMake(threads_per_group, 1, 1)];
    }
}

void Batch::commit()
{
    if (impl_->committed)
        return;
    @autoreleasepool
    {
        if (impl_->encoder)
            [impl_->encoder endEncoding];
        impl_->encoder = nil;
        [impl_->commands commit];
        impl_->committed = true;
    }
}

void Batch::wait()
{
    commit();
    @autoreleasepool
    {
        [impl_->commands waitUntilCompleted];
        if (impl_->commands.status == MTLCommandBufferStatusError)
        {
            NSError* error = impl_->commands.error;
            throw std::runtime_error(std::string("GPU batch failed: ")
                    + (error ? error.localizedDescription.UTF8String : "unknown error"));
        }
    }
}

double Batch::gpu_seconds() const
{
    return impl_->commands.GPUEndTime - impl_->commands.GPUStartTime;
}

Device* Device::get()
{
    static std::unique_ptr<Device> device = [] {
        std::unique_ptr<Device> d;
        @autoreleasepool
        {
            id<MTLDevice> metal = MTLCreateSystemDefaultDevice();
            if (!metal)
                return d;
            d.reset(new Device());
            d->impl_ = std::make_shared<Impl>();
            d->impl_->device = metal;
            d->impl_->queue = [metal newCommandQueue];
        }
        return d;
    }();
    return device.get();
}

const char* Device::backend()
{
    return "metal";
}

std::string Device::name() const
{
    return impl_->device.name.UTF8String;
}

Buffer Device::buffer(std::size_t bytes)
{
    Buffer b;
    b.impl_ = std::make_shared<Buffer::Impl>();
    b.impl_->size = bytes;
    @autoreleasepool
    {
        // Metal will not make an empty buffer; a one-byte one stands in.
        b.impl_->buffer = [impl_->device newBufferWithLength:std::max<std::size_t>(bytes, 1)
                                                     options:MTLResourceStorageModeShared];
    }
    if (!b.impl_->buffer)
        throw std::runtime_error("GPU: cannot allocate a buffer of " + std::to_string(bytes) + " bytes");
    std::memset(b.data(), 0, bytes);
    return b;
}

Kernel Device::kernel(const std::vector<std::string>& sources, const std::string& entry)
{
    std::string key;
    for (const std::string& s : sources)
        key += s + ";";
    std::lock_guard<std::mutex> guard(impl_->lock);
    Kernel k;
    if (auto it = impl_->kernels.find(key + entry); it != impl_->kernels.end())
    {
        k.impl_ = it->second;
        return k;
    }
    @autoreleasepool
    {
        id<MTLLibrary> library = impl_->libraries[key];
        if (!library)
        {
            std::string program;
            for (const std::string& name : std::vector<std::string> {"prelude_metal.h", "common.h"})
                program += std::string(kernel_source(name)) + "\n";
            for (const std::string& name : sources)
            {
                const char* text = kernel_source(name);
                if (!text)
                    throw std::runtime_error("GPU: no kernel source " + name);
                program += std::string(text) + "\n";
            }
            MTLCompileOptions* options = [MTLCompileOptions new];
            // Safe math: the kernels' compensated sums and their explicit
            // fast/precise choices mean what they say.
            if (@available(macOS 15.0, *))
                options.mathMode = MTLMathModeSafe;
            else
            {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
                options.fastMathEnabled = NO;
#pragma clang diagnostic pop
            }
            NSError* error = nil;
            library = [impl_->device newLibraryWithSource:@(program.c_str()) options:options error:&error];
            if (!library)
                throw std::runtime_error("GPU: kernel sources " + key + " do not compile:\n"
                        + (error ? error.localizedDescription.UTF8String : "unknown error"));
            impl_->libraries[key] = library;
        }
        id<MTLFunction> function = [library newFunctionWithName:@(entry.c_str())];
        if (!function)
            throw std::runtime_error("GPU: no kernel " + entry + " in " + key);
        NSError* error = nil;
        auto impl = std::make_shared<Kernel::Impl>();
        impl->pipeline = [impl_->device newComputePipelineStateWithFunction:function error:&error];
        if (!impl->pipeline)
            throw std::runtime_error("GPU: cannot build kernel " + entry + ": "
                    + (error ? error.localizedDescription.UTF8String : "unknown error"));
        impl_->kernels[key + entry] = impl;
        k.impl_ = impl;
    }
    return k;
}

Batch Device::batch()
{
    Batch b;
    b.impl_ = std::make_shared<Batch::Impl>();
    @autoreleasepool
    {
        b.impl_->commands = [impl_->queue commandBuffer];
    }
    return b;
}

} // namespace niceops::gpu
