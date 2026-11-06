// The kernel dialect on CUDA and HIP, compiled by NVRTC or hiprtc, each of
// which takes this same source (gpu_cuda.cpp, gpu_hip.cpp).  See prelude_metal.h for what the dialect
// allows; this maps the same names.  There is no standard library here: NVRTC
// and hiprtc supply the float overloads of sqrt, exp, log, cos, atan2, fma, rint,
// floor, fabs, fmax, isfinite and the integer min and max as built-ins.

#define NK_DIALECT 1

typedef unsigned int uint;
typedef unsigned int uint32_t;
typedef unsigned char nk_u8;

// Arguments bind by position, which the slots already are; the slot is
// checked by nothing here, so a kernel's slots must run 0, 1, 2, ...
#define NK_KERNEL(name) extern "C" __global__ void name
#define NK_IN(type, name, slot) const type* name
#define NK_OUT(type, name, slot) type* name
#define NK_VALUE(type, name, slot) const type name
#define NK_GROUP_ARGS
#define NK_GROUP blockIdx.x
#define NK_THREAD threadIdx.x
#define NK_THREADS blockDim.x

#define NK_CONSTANT __constant__

// One generic address space: the qualifiers say nothing but const.
#define NK_FN static __device__ inline
#define NK_GLOBAL_IN const
#define NK_GLOBAL_OUT
#define NK_LOCAL
#define NK_SHARED
#define NK_SHARED_ARRAY(type, name, count) __shared__ type name[count]
#define NK_BARRIER() __syncthreads()

// sin and cos of x together, to full float accuracy (the backend compiles
// without fast math).
static __device__ inline float nk_sincos(float x, float& c)
{
    float s;
    sincosf(x, &s, &c);
    return s;
}
