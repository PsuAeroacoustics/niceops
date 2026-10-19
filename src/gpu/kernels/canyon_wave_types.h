// The canyon wave solver's grid as its kernels take it (canyon_wave.kern),
// shared with the host (canyon_wave.cpp): fixed-width integers and floats.
#ifndef NICEOPS_CANYON_WAVE_TYPES_H
#define NICEOPS_CANYON_WAVE_TYPES_H

#ifndef NK_DIALECT
#include <cstdint>
#endif

// canyon_wave.cpp's Grid and its run's constants: nx columns by nz rows,
// pressure at the cell centres, vy on the faces between columns, vz on those
// between rows.
struct CanyonGridGpu
{
    uint32_t nx;
    uint32_t nz;
    uint32_t pml;
    uint32_t wall_rows;
    uint32_t record_row;
    uint32_t record_lo;
    uint32_t record_hi;
    uint32_t source_i;
    uint32_t source_j;
    uint32_t energy_top;
    // Whether wall A and wall B are impedance walls.
    uint32_t wall_a;
    uint32_t wall_b;
    float inv_rho_h;
    float rc2_h;
    float inject;
    float keep_a;
    float push_a;
    float keep_b;
    float push_b;
    float unused;
};

#endif
