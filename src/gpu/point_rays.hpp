#pragma once

// A point-source field's rays packed for the GPU to form their amplitudes
// itself (gpu/kernels/point_rays.kern), at every sub-band of every band, as
// PointSourceField::sub_band_problems hands them to the CPU.  What the
// amplitudes owe to the wavenumber is left to the GPU; everything else is
// worked out here, in double, once per field: each ray's amplitude after its
// reflections per band, its wedge terms' angles, and over a shared face the
// aperture's nodes for all the run's wavenumbers at once
// (shared_face_geometry).

#include <cstdint>
#include <vector>

#include "gpu/kernels/point_rays_types.h"

namespace niceops {
class PointSourceField;
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
    std::vector<float> reflected; // [ray][band]
    // Which rays are over a shared face.
    std::vector<uint32_t> face_rays;
    // The problems: per band, its three sub-bands in order, and the air's
    // attenuation in each (Np/m, 0 for none).
    std::vector<uint32_t> band;
    std::vector<double> k1, k2;
    std::vector<float> absorption;
    // As RaySet::nodes (band_energy.hpp).
    std::vector<uint32_t> nodes;

    std::size_t problems() const
    {
        return k1.size();
    }
};

// Packs `field` for the bands (Hz), the facades' reflection coefficients
// [facade][band] and the air's absorption (sub_band_absorption, or empty) as
// band_levels takes them, and plans the problems.  Throws
// std::invalid_argument for no bands.
void pack_point_rays(const PointSourceField& field, const std::vector<double>& reflection, double floor_reflection,
        const std::vector<double>& bands, double sound_speed, PointRaySet& set,
        const std::vector<double>& absorption = {});

} // namespace niceops::gpu
