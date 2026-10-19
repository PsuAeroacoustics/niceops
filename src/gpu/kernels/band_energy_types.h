// The band-energy kernels' problem record, shared by the host (band_energy.cpp)
// and the kernels: fixed-width integers and floats only, so the two layouts agree.
#ifndef NICEOPS_BAND_ENERGY_TYPES_H
#define NICEOPS_BAND_ENERGY_TYPES_H

#ifndef NK_DIALECT
#include <cstdint>
#endif

// One band on one set of rays.  The set's rays are [rays, rays + n) of the
// per-ray arrays (length, rel), sorted by length; the problem's own entries
// are [first, first + n) of the per-problem arrays (amplitude, phasor, trig).
struct BandProblem
{
    uint32_t rays;
    uint32_t first;
    uint32_t n;
    // Gauss-Legendre nodes, and where their (t, w) pairs start in the rule
    // array; 0 nodes for the pair sum.
    uint32_t nodes;
    uint32_t rule;
    // The band's middle wavenumber and its half width, 1/m, each as a
    // float-float pair (high, low), for the phases; h again as the sums use it.
    float k_hi, k_lo, h_hi, h_lo;
    float h;
    // The pair sum's cutoff in path-length difference, 30 / h, m.
    float reach;
    // The quadrature's centre in the rel lengths, m.
    float centre;
    // Where the GPU forms the amplitudes itself (point_rays.kern): the band
    // the problem is in, for the rays' reflections, and the air's attenuation
    // of an amplitude per metre of path, Np/m (0 for none).
    uint32_t band;
    float absorption;
    // How its rays add (band_pair_sum): NK_SUM_COHERENT, all with their
    // phases; NK_SUM_CONE, with their phases within `cone` degrees of each
    // other's launch (the per-ray launch array), a cos^2 taper over a quarter
    // of it either side, in energy beyond; NK_SUM_ENERGY, all in energy.  Only
    // coherent problems may take the quadrature.
    uint32_t sum;
    float cone;
};

#define NK_SUM_COHERENT 0u
#define NK_SUM_CONE 1u
#define NK_SUM_ENERGY 2u

#endif
