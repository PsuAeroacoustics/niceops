// Run at configure time (CMakeLists.txt) with the HIP backend: exits 0 where
// the HIP runtime starts and sees a device, as gpu_hip.cpp's Device::get()
// asks, and prints why not otherwise.  It decides whether the tests that need
// a device are registered.
#include <dlfcn.h>
#include <hip/hip_runtime_api.h>

#include <cstdio>
#include <string>

#define NK_STR(x) #x
#define NK_XSTR(x) NK_STR(x)

int main()
{
    const std::string major = std::to_string(HIP_VERSION_MAJOR);
    void* h = dlopen(("libamdhip64.so." + major).c_str(), RTLD_NOW | RTLD_LOCAL);
#ifdef NICEOPS_ROCM_LIB_DIR
    if (!h)
        h = dlopen((std::string(NICEOPS_ROCM_LIB_DIR) + "/libamdhip64.so." + major).c_str(), RTLD_NOW | RTLD_LOCAL);
#endif
    if (!h)
    {
        std::printf("the HIP runtime (libamdhip64.so.%s) is not installed", major.c_str());
        return 1;
    }
    auto init = reinterpret_cast<decltype(&hipInit)>(dlsym(h, "hipInit"));
    auto count = reinterpret_cast<decltype(&hipGetDeviceCount)>(dlsym(h, "hipGetDeviceCount"));
    auto properties = reinterpret_cast<decltype(&hipGetDeviceProperties)>(dlsym(h, NK_XSTR(hipGetDeviceProperties)));
    if (!init || !count || !properties)
    {
        std::printf("libamdhip64 lacks the runtime API");
        return 1;
    }
    if (init(0) != hipSuccess)
    {
        std::printf("the HIP runtime did not start");
        return 1;
    }
    int n = 0;
    if (count(&n) != hipSuccess || n == 0)
    {
        std::printf("the HIP runtime sees no device");
        return 1;
    }
    hipDeviceProp_t p {};
    if (properties(&p, 0) != hipSuccess)
    {
        std::printf("the HIP runtime cannot describe its device");
        return 1;
    }
    std::printf("found %s (%s)", p.name, p.gcnArchName);
    return 0;
}
