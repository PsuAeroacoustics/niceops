// The FDTD's grid as its kernels take it (fdtd.kern), shared with the host
// (fdtd.cpp): fixed-width integers and floats.
#ifndef NICEOPS_FDTD_TYPES_H
#define NICEOPS_FDTD_TYPES_H

#ifndef NK_DIALECT
#include <cstdint>
#endif

// nx by ny by nz cells, x fastest; the PML's thickness on each side, cells
// (0 for a rigid side); the steps' constants, dt / (rho h) for the velocity
// and rho c^2 dt / h for the pressure; and, where `walls` is not 0, the
// solids' faces' impedance update, v' = keep v -+ push p (fdtd.cpp).
struct FdtdGridGpu
{
    uint32_t nx;
    uint32_t ny;
    uint32_t nz;
    uint32_t x0;
    uint32_t x1;
    uint32_t y0;
    uint32_t y1;
    uint32_t z0;
    uint32_t z1;
    float push_v;
    float push_p;
    float keep;
    float push_wall;
    uint32_t walls;
};

// The cells a step marches: ni by nj by nk of them from (i, j, k).
struct FdtdBox
{
    uint32_t i;
    uint32_t j;
    uint32_t k;
    uint32_t ni;
    uint32_t nj;
    uint32_t nk;
};

// A cell's flags: whether it is air, whether the face on its low side along
// x, y and z is open (air on both sides and not the domain's edge), whether
// that face is a solid's, air on one side only, and whether the cell is clean.
#define FDTD_AIR 1u
#define FDTD_OPEN_X 2u
#define FDTD_OPEN_Y 4u
#define FDTD_OPEN_Z 8u
#define FDTD_WALL_X 16u
#define FDTD_WALL_Y 32u
#define FDTD_WALL_Z 64u
#define FDTD_FACES 126u
// The cell and the two either side of it along every axis are air with all
// their low faces open, and it lies in no PML slab: every difference there is
// the plain one, with no mirror to look for.
#define FDTD_CLEAN 128u

#endif
