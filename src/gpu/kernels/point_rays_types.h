// The point-source model's rays as the GPU forms their amplitudes
// (point_rays.kern), shared by the host (gpu/point_rays.cpp) and the kernels:
// fixed-width integers and floats only, so the two layouts agree.  Angles'
// parts are worked out in double on the host.
#ifndef NICEOPS_POINT_RAYS_TYPES_H
#define NICEOPS_POINT_RAYS_TYPES_H

#ifndef NK_DIALECT
#include <cstdint>
#endif

// One wedge on a ray (WedgeTerm): its angles' parts (WedgeAngles), n,
// sin(beta), its distance parameter, and (grazing ? 1/2 : 1) / sqrt(sin beta).
// A wedge term (WedgeTerm): its angles' part, n, sin(beta), distance and
// scale, and its edge's ends (EdgeEnds): each end's path excess, m, -1 for
// none, where the stationary point lies (-1, 0, +1), and past an end the
// factor's scale.
struct GpuWedge
{
    float cot[4];
    float c[4];
    float n;
    float sin_beta;
    float distance;
    float scale;
    float end_lo;
    float end_hi;
    float end_past;
    float end_scale;
};

// One ray, in length order: its kind (0 specular, 1 to NK_KIND_FACE - 1 its
// wedge terms, NK_KIND_FACE over a shared face), its first wedge term or its face, where its
// amplitude after reflections starts (reflected_amplitude, real per band, or
// times the floor's factor complex per problem: BandProblem::slot), and
// its length, m, for the air's absorption.
#define NK_KIND_FACE 15u
struct GpuRay
{
    uint32_t kind;
    uint32_t term;
    uint32_t reflected;
    float length;
};

// One end's aperture of a pair of parallel edges (ApertureGeometry): its
// panels [panel, panel + panels), the end's weight (0 where it has none), its
// extent scale, and its end correction's angle below the line between the
// edges (sin, cos), Keller coefficient's angles (n1 the first wedge's) and
// share of that coefficient (keller_scale).
struct GpuAperture
{
    uint32_t panel;
    uint32_t panels;
    float weight;
    float extent_scale;
    float sin_g;
    float cos_g;
    float keller_cot[4];
    float keller_c[4];
    float n1;
    float keller_scale;
};

// A ray over a shared face: its two ends' apertures, sin(beta), and the
// reference path less the ray's length in the plane, as a float-float pair
// (the phase the aperture's lengths leave out).  Sampled
// (PointPathOptions::sample_shared_faces), its apertures hold only the probes'
// wavenumbers, and `cascade` is the first of its two wedge terms, which
// ray_amplitudes takes at every problem (the first's scale carrying the
// cascade's spreading over the face's); taken at every wavenumber, `cascade`
// is no_cascade.  `ends` is the first of the same two terms where either
// edge has an end (EdgeEnds), no_cascade where neither does: the integral
// takes each one's coefficient with its ends over without
// (shared_face_amplitude).
struct GpuFace
{
    GpuAperture end[2];
    float sin_beta;
    float offset_hi;
    float offset_lo;
    uint32_t cascade;
    uint32_t ends;
};

#define NK_NO_CASCADE 0xffffffffu

// A sampled shared face's ray, for face_ratios: the problems of its set
// [first, first + problems), the ray within the set, its probes
// [probe, probe + probes) of the probe arrays (each probe's problem within
// the set, ascending in k, and the integral face_amplitudes put there), and
// the problem within the set from which the cascade stands alone (problems
// for none).
struct GpuFaceSample
{
    uint32_t first;
    uint32_t problems;
    uint32_t ray;
    uint32_t probe;
    uint32_t probes;
    uint32_t free;
    uint32_t unused[2];
};

// One Filon panel of an aperture: its three nodes' places up the aperture;
// which waves count (bit j for wave j); per wave the start node's path length
// less the reference (a float-float pair) and the other nodes' distances from
// it, formed in double; per wave and node the amplitude over the rate, less
// launch, taper and the second wedge's coefficient; and per node that
// coefficient's angles and distance parameter, and the second wedge's n.
struct GpuPanel
{
    float t[3];
    uint32_t waves;
    float u0[4];
    float u0_lo[4];
    float d1[4];
    float d2[4];
    float g[4][3];
    float cot[3][4];
    float c[3][4];
    float distance[3];
    float n2;
};

// A shared face's ray at one problem: the problem (its index in the batch),
// the ray (its index within its own set, which face_amplitudes offsets by the
// problem's first ray) and its face (its index in the batch); `probe` is 0,
// or one more than the probe slot the integral goes to instead of the
// amplitude (a sampled face).
struct GpuFaceItem
{
    uint32_t problem;
    uint32_t ray;
    uint32_t face;
    uint32_t probe;
};

#endif
