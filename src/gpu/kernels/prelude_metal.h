// The kernel dialect on Metal.  Kernels in src/gpu/kernels are written against
// these names only, so that another backend's prelude -- CUDA's __global__,
// blockIdx and __shared__, say -- takes the same source.  What the dialect
// allows: C++14 without the standard library, exceptions, recursion or double;
// float, 32-bit integers and nk_u8; plain structs shared with the host through
// the *_types.h headers (fixed-width types and float only, so the layouts
// agree).
#include <metal_stdlib>
using namespace metal;

// Defined by every prelude: the *_types.h headers include <cstdint> only
// outside a kernel program.
#define NK_DIALECT 1

// An unsigned byte, for per-cell flags where a word would cost bandwidth.
typedef uchar nk_u8;

// A kernel: NK_KERNEL(name)(NK_IN(float, x, 0), NK_OUT(float, y, 1) NK_GROUP_ARGS)
// -- the group arguments carry their own leading comma.
#define NK_KERNEL(name) kernel void name
#define NK_IN(type, name, slot) device const type* name [[buffer(slot)]]
#define NK_OUT(type, name, slot) device type* name [[buffer(slot)]]
#define NK_VALUE(type, name, slot) constant type& name [[buffer(slot)]]
#define NK_GROUP_ARGS , uint nk_group [[threadgroup_position_in_grid]], \
    uint nk_thread [[thread_position_in_threadgroup]], uint nk_threads [[threads_per_threadgroup]]
#define NK_GROUP nk_group
#define NK_THREAD nk_thread
#define NK_THREADS nk_threads

// A table at program scope: NK_CONSTANT float table[4] = {...};
#define NK_CONSTANT constant

// Helper functions, and the address spaces their pointers and references name:
// an input buffer, an output buffer, a local variable, threadgroup memory.
#define NK_FN static inline
#define NK_GLOBAL_IN device const
#define NK_GLOBAL_OUT device
#define NK_LOCAL thread
#define NK_SHARED threadgroup
// Threadgroup memory: NK_SHARED_ARRAY(float, scratch, 256); and the barrier
// every thread of the group must reach.
#define NK_SHARED_ARRAY(type, name, count) threadgroup type name[count]
#define NK_BARRIER() threadgroup_barrier(mem_flags::mem_threadgroup)

// sin and cos of x together, to the backend's full float accuracy.
static inline float nk_sincos(float x, thread float& c)
{
    return precise::sincos(x, c);
}
