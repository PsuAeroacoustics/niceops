// gpu.hpp for a build without a GPU backend (NICEOPS_GPU=OFF): no device.
// Nothing below get() is reachable, since there is never a Device to call it
// on; each still says so rather than doing something plausible.

#include "gpu/gpu.hpp"

#include <stdexcept>

namespace niceops::gpu {

struct Buffer::Impl
{
};
struct Kernel::Impl
{
};
struct Batch::Impl
{
};
struct Device::Impl
{
};

namespace {

[[noreturn]] void none()
{
    throw std::logic_error("this build has no GPU backend (NICEOPS_GPU=OFF)");
}

} // namespace

std::size_t Buffer::size() const
{
    return 0;
}

void* Buffer::data() const
{
    return nullptr;
}

std::size_t Kernel::max_threads_per_group() const
{
    return 0;
}

void Batch::dispatch(const Kernel&, const std::vector<Arg>&, std::size_t, std::size_t)
{
    none();
}

void Batch::commit()
{
    none();
}

void Batch::wait()
{
    none();
}

double Batch::gpu_seconds() const
{
    return 0.0;
}

Device* Device::get()
{
    return nullptr;
}

const char* Device::backend()
{
    return "none";
}

std::string Device::name() const
{
    none();
}

std::size_t Device::max_buffer_bytes() const
{
    none();
}

std::size_t Device::working_set_bytes() const
{
    none();
}

Buffer Device::buffer(std::size_t)
{
    none();
}

Kernel Device::kernel(const std::vector<std::string>&, const std::string&)
{
    none();
}

Batch Device::batch()
{
    none();
}

} // namespace niceops::gpu
