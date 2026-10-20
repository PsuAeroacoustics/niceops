#include "gpu/point_rays.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>

#include "double_wedge.hpp"
#include "gpu/band_energy.hpp"
#include "gpu/gpu.hpp"
#include "point_paths.hpp"

namespace niceops::gpu {

namespace {

void put_angles(const WedgeAngles& a, float* cot, float* c)
{
    for (int i = 0; i < 4; i++)
    {
        cot[i] = static_cast<float>(a.cot[i]);
        c[i] = static_cast<float>(a.c[i]);
    }
}

// One end's aperture: its panels, three nodes each, appended to `panels`.
GpuAperture pack_aperture(const ApertureGeometry& g, double weight, std::vector<GpuPanel>& panels)
{
    GpuAperture a {};
    a.panel = static_cast<uint32_t>(panels.size());
    a.weight = static_cast<float>(weight);
    a.extent_scale = static_cast<float>(g.extent_scale);
    a.sin_g = static_cast<float>(g.sin_g);
    a.cos_g = static_cast<float>(g.cos_g);
    put_angles(g.keller, a.keller_cot, a.keller_c);
    a.n1 = static_cast<float>(g.keller.n);
    for (std::size_t s = 0; s + 1 < g.segments.size(); s++)
        for (std::size_t i = g.segments[s]; i + 2 < g.segments[s + 1]; i += 2)
        {
            const ApertureNode* n = &g.nodes[i];
            GpuPanel p {};
            p.waves = 0;
            for (int m = 0; m < 3; m++)
            {
                p.t[m] = static_cast<float>(n[m].t);
                put_angles(n[m].angles, p.cot[m], p.c[m]);
                p.distance[m] = static_cast<float>(n[m].distance);
            }
            p.n2 = static_cast<float>(n[0].angles.n);
            for (int j = 0; j < 4; j++)
            {
                if (n[0].weight[j] == 0.0 && n[2].weight[j] == 0.0)
                    continue;
                p.waves |= 1u << j;
                split(n[0].length[j], p.u0[j], p.u0_lo[j]);
                p.d1[j] = static_cast<float>(n[1].length[j] - n[0].length[j]);
                p.d2[j] = static_cast<float>(n[2].length[j] - n[0].length[j]);
                for (int m = 0; m < 3; m++)
                    p.g[j][m] = static_cast<float>(n[m].weight[j] / n[m].rate[j]);
            }
            panels.push_back(p);
        }
    a.panels = static_cast<uint32_t>(panels.size() - a.panel);
    return a;
}

} // namespace

void pack_point_rays(const PointSourceField& field, const std::vector<double>& reflection, const FloorModel& floor,
        const std::vector<double>& bands, double sound_speed, PointRaySet& set, const std::vector<double>& absorption,
        const PointRayPacking& packing)
{
    if (bands.empty())
        throw std::invalid_argument("pack_point_rays: no bands");
    if (!absorption.empty() && absorption.size() != 3 * bands.size())
        throw std::invalid_argument("pack_point_rays: the absorption is not three per band");
    set = PointRaySet();
    set.sum = packing.sum;
    set.cone = static_cast<float>(packing.cone);
    const std::vector<PointRay>& rays = field.rays();
    const std::size_t nb = bands.size(), facades = nb ? reflection.size() / nb : 0;
    std::vector<std::size_t> order(rays.size());
    for (std::size_t i = 0; i < order.size(); i++)
        order[i] = i;
    // The order sub_band_problems sums them in.
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return rays[a].length < rays[b].length; });

    // The problems, and every wavenumber an aperture is wanted at.
    std::vector<double> centres;
    for (std::size_t band = 0; band < nb; band++)
    {
        const std::array<SubBand, 3> parts = sub_bands(bands[band], sound_speed);
        for (int part = 0; part < 3; part++)
        {
            const SubBand& s = parts[part];
            set.k1.push_back(s.k1);
            set.k2.push_back(s.k2);
            set.absorption.push_back(absorption.empty() ? 0.0f : static_cast<float>(absorption[3 * band + part]));
            centres.push_back(0.5 * (s.k1 + s.k2));
        }
    }
    std::vector<std::vector<double>> band_reflection(nb, std::vector<double>(facades));
    for (std::size_t band = 0; band < nb; band++)
        for (std::size_t f = 0; f < facades; f++)
            band_reflection[band][f] = reflection[f * nb + band];

    if (floor.impedance() && floor.admittance.size() != centres.size())
        throw std::invalid_argument("pack_point_rays: the floor's admittance is not three per band");
    const double floor_reflection = floor.impedance() ? 1.0 : floor.reflection;
    std::vector<std::complex<double>> floor_factor(centres.size(), 1.0), plate(floor.plate ? 2 * centres.size() : 0);
    std::vector<double> in_plane;
    for (std::size_t index : order)
    {
        const PointRay& ray = rays[index];
        GpuRay r {};
        r.reflected = static_cast<uint32_t>(set.reflected.size());
        r.length = static_cast<float>(ray.length);
        if (floor.plate)
        {
            double elevation, azimuth;
            point_ray_arrival(ray, elevation, azimuth);
            floor.plate(elevation, azimuth, plate.data(), plate.data() + centres.size());
        }
        floor_factors(ray, floor, centres, floor.plate ? plate.data() : nullptr, floor_factor.data());
        for (std::size_t band = 0; band < nb; band++)
        {
            const double amplitude = reflected_amplitude(ray, band_reflection[band], floor_reflection)
                    * (packing.weights != nullptr ? packing.weights[index * nb + band] : 1.0);
            for (std::size_t part = 0; part < 3; part++)
            {
                const std::complex<double> a = amplitude * floor_factor[3 * band + part];
                set.reflected.push_back(static_cast<float>(a.real()));
                set.reflected.push_back(static_cast<float>(a.imag()));
            }
        }
        const Vec3 launch = point_ray_launch(ray);
        for (double v : {launch.x, launch.y, launch.z})
            set.launch.push_back(static_cast<float>(v));
        auto put_wedge = [&](const WedgeTerm& t, double factor) {
            GpuWedge g {};
            put_angles(t.angles, g.cot, g.c);
            g.n = static_cast<float>(t.angles.n);
            g.sin_beta = static_cast<float>(t.sin_beta);
            g.distance = static_cast<float>(t.distance);
            g.scale = static_cast<float>(factor * (t.grazing ? 0.5 : 1.0) / std::sqrt(t.sin_beta));
            set.wedges.push_back(g);
        };
        if (ray.shared)
        {
            const SharedFaceTerm& t = ray.shared_face;
            // Sampled, the integral is wanted at the probes only, and the
            // cascade, its spreading's ratio to the face's on its first
            // term, at every problem.
            SharedFaceSamples sampled;
            const bool sample = field.samples_shared_faces();
            if (sample)
            {
                sampled = shared_face_samples(t, centres);
                in_plane.clear();
                for (std::size_t p : sampled.probes)
                    in_plane.push_back(centres[p] * t.sin_beta);
            }
            else
            {
                in_plane.assign(centres.begin(), centres.end());
                for (double& k : in_plane)
                    k *= t.sin_beta;
            }
            GpuFace f {};
            f.sin_beta = static_cast<float>(t.sin_beta);
            f.cascade = NK_NO_CASCADE;
            if (!in_plane.empty())
            {
                const SharedFaceGeometry g = shared_face_geometry(t.face, t.sx, t.sz, t.rx, t.rz, in_plane);
                f.end[0] = g.forward_weight > 0.0 ? pack_aperture(g.forward, g.forward_weight, set.panels) : GpuAperture {};
                f.end[1] = g.reverse_weight > 0.0 ? pack_aperture(g.reverse, g.reverse_weight, set.panels) : GpuAperture {};
                const double reference = g.forward_weight > 0.0 ? g.forward.reference : g.reverse.reference;
                split(reference - ray.length * t.sin_beta, f.offset_hi, f.offset_lo);
            }
            const uint32_t face = static_cast<uint32_t>(set.faces.size());
            const uint32_t self = static_cast<uint32_t>(set.rays.size());
            if (sample)
            {
                f.cascade = static_cast<uint32_t>(set.wedges.size());
                put_wedge(ray.wedge[0], ray.cascade_spreading / ray.spreading);
                put_wedge(ray.wedge[1], 1.0);
                GpuFaceSample record {};
                record.problems = static_cast<uint32_t>(centres.size());
                record.ray = self;
                record.probe = static_cast<uint32_t>(set.probe_problems.size());
                record.probes = static_cast<uint32_t>(sampled.probes.size());
                record.free = static_cast<uint32_t>(sampled.free);
                for (std::size_t p : sampled.probes)
                {
                    set.probe_problems.push_back(static_cast<uint32_t>(p));
                    set.probe_items.push_back({static_cast<uint32_t>(p), self, face,
                            static_cast<uint32_t>(set.probe_problems.size())});
                }
                set.samples.push_back(record);
            }
            else
                set.face_rays.push_back(self);
            r.kind = NK_KIND_FACE;
            r.term = face;
            set.faces.push_back(f);
        }
        else
        {
            r.kind = static_cast<uint32_t>(ray.wedges);
            r.term = static_cast<uint32_t>(set.wedges.size());
            for (int w = 0; w < ray.wedges; w++)
                put_wedge(ray.wedge[w], 1.0);
        }
        set.rays.push_back(r);
        set.length.push_back(ray.length);
    }
    plan_nodes(set.length, set.k1, set.k2, set.nodes);
    // The quadrature integrates every pair with its phases.
    if (set.sum != NK_SUM_COHERENT)
        std::fill(set.nodes.begin(), set.nodes.end(), 0u);
}

} // namespace niceops::gpu
