// Run at configure time (CMakeLists.txt) with the CUDA backend: exits 0 where
// the driver starts and its first device takes the backend (concurrent managed
// memory), as gpu_cuda.cpp's Device::get() asks, and prints why not otherwise.
// It decides whether the tests that need a device are registered.
#include <cuda.h>
#include <dlfcn.h>

#include <cstdio>

int main()
{
    void* h = dlopen("libcuda.so.1", RTLD_NOW | RTLD_LOCAL);
    if (!h)
    {
        std::printf("the CUDA driver library (libcuda.so.1) is not installed");
        return 1;
    }
    auto init = reinterpret_cast<CUresult (*)(unsigned)>(dlsym(h, "cuInit"));
    auto count = reinterpret_cast<CUresult (*)(int*)>(dlsym(h, "cuDeviceGetCount"));
    auto get = reinterpret_cast<CUresult (*)(CUdevice*, int)>(dlsym(h, "cuDeviceGet"));
    auto attribute = reinterpret_cast<CUresult (*)(int*, CUdevice_attribute, CUdevice)>(dlsym(h, "cuDeviceGetAttribute"));
    if (!init || !count || !get || !attribute)
    {
        std::printf("libcuda.so.1 lacks the driver API");
        return 1;
    }
    const CUresult r = init(0);
    if (r != CUDA_SUCCESS)
    {
        std::printf("the CUDA driver did not start (error %d)", static_cast<int>(r));
        return 1;
    }
    int n = 0;
    CUdevice device = 0;
    if (count(&n) != CUDA_SUCCESS || n == 0 || get(&device, 0) != CUDA_SUCCESS)
    {
        std::printf("the CUDA driver sees no device");
        return 1;
    }
    int concurrent = 0;
    attribute(&concurrent, CU_DEVICE_ATTRIBUTE_CONCURRENT_MANAGED_ACCESS, device);
    if (!concurrent)
    {
        std::printf("the device has no concurrent managed memory");
        return 1;
    }
    std::printf("found");
    return 0;
}
