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
struct GpuWedge
{
    float cot[4];
    float c[4];
    float n;
    float sin_beta;
    float distance;
    float scale;
};

// One ray, in length order: its kind (0 specular, 1 or 2 its wedge terms, 3
// over a shared face), its first wedge term or its face, where its amplitude
// after reflections (reflected_amplitude) starts, one per band, and its
// length, m, for the air's absorption.
struct GpuRay
{
    uint32_t kind;
    uint32_t term;
    uint32_t reflected;
    float length;
};

// One end's aperture of a shared face (ApertureGeometry): its panels
// [panel, panel + panels), the end's weight (0 where it has none), its extent
// scale, and its end correction's angle below the face (sin, cos) and Keller
// coefficient's angles (n1 the first wedge's).
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
    float unused;
};

// A ray over a shared face: its two ends' apertures, sin(beta), and the
// reference path less the ray's length in the plane, as a float-float pair
// (the phase the aperture's lengths leave out).
struct GpuFace
{
    GpuAperture end[2];
    float sin_beta;
    float offset_hi;
    float offset_lo;
    float unused;
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
// problem's first ray) and its face (its index in the batch).
struct GpuFaceItem
{
    uint32_t problem;
    uint32_t ray;
    uint32_t face;
    uint32_t unused;
};

#endif
