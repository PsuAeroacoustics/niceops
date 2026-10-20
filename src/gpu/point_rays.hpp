#pragma once

// A point-source field's rays packed for the GPU to form their amplitudes
// itself (gpu/kernels/point_rays.kern), at every sub-band of every band, as
// PointSourceField::sub_band_problems hands them to the CPU.  What the
// amplitudes owe to the wavenumber is left to the GPU; everything else is
// worked out here, in double, once per field: each ray's amplitude after its
// reflections per sub-band (complex, a floor of finite impedance's Q
// included), its wedge terms' angles, and over a shared face the
// aperture's nodes for all the run's wavenumbers at once
// (shared_face_geometry).

#include <cstdint>
#include <vector>

#include "gpu/kernels/band_energy_types.h"
#include "gpu/kernels/point_rays_types.h"

namespace niceops {
class PointSourceField;
struct FloorModel;
}

namespace niceops::gpu {

struct PointRaySet
{
    // The rays' lengths, ascending, and the rays in that order.  A ray's term
    // and reflected offsets, and a face's panel offsets, index this set's own
    // arrays.
    std::vector<double> length;
    std::vector<GpuRay> rays;
    std::vector<GpuWedge> wedges;
    std::vector<GpuFace> faces;
    std::vector<GpuPanel> panels;
    // Each ray's amplitude after its reflections, the floor's included, and
    // the caller's weights, per problem: [ray][problem][re, im].
    std::vector<float> reflected;
    // Which rays are over a shared face whose integral is taken at every
    // problem.
    std::vector<uint32_t> face_rays;
    // The sampled shared faces (PointPathOptions::sample_shared_faces): one
    // record each, its problems, ray and probes set-local; each probe's
    // problem within the set; and the integrals to take, at the probes, each
    // into its probe's slot (one more than its index).
    std::vector<GpuFaceSample> samples;
    std::vector<uint32_t> probe_problems;
    std::vector<GpuFaceItem> probe_items;
    // The problems: per band, its three sub-bands in order, and the air's
    // attenuation in each (Np/m, 0 for none).
    std::vector<double> k1, k2;
    std::vector<float> absorption;
    // As RaySet::nodes (band_energy.hpp).
    std::vector<uint32_t> nodes;
    // How the rays add (BandProblem::sum, cone), and each ray's launch, a unit
    // vector, three floats a ray, for the cone.
    uint32_t sum = NK_SUM_COHERENT;
    float cone = 1.0f;
    std::vector<float> launch;

    std::size_t problems() const
    {
        return k1.size();
    }
};

// What a caller may add to a set beyond band_levels' arguments: `weights`,
// each ray's amplitude in each band, [ray][band] in rays()'s order (a
// directional source along its launch, as band_levels takes it), and how the
// rays add (PointRaySet::sum, cone).
struct PointRayPacking
{
    const double* weights = nullptr;
    uint32_t sum = NK_SUM_COHERENT;
    double cone = 1.0;
};

// Packs `field` for the bands (Hz), the facades' reflection coefficients
// [facade][band] and the air's absorption (sub_band_absorption, or empty) as
// band_levels takes them, and plans the problems.  Throws
// std::invalid_argument for no bands.
void pack_point_rays(const PointSourceField& field, const std::vector<double>& reflection, const FloorModel& floor,
        const std::vector<double>& bands, double sound_speed, PointRaySet& set,
        const std::vector<double>& absorption = {}, const PointRayPacking& packing = {});

} // namespace niceops::gpu
