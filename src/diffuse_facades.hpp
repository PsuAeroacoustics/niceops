#pragma once

// Diffuse reflection by the facades, for the urban tables' general builder
// (building_paths.hpp; THEORY 6.5).  A facade reflects the fraction rho of the
// energy that reaches it, and scatters the fraction s of that diffusely, by
// Lambert's law; the rest, rho (1 - s), it reflects specularly, which the
// rays carry.  The diffuse energy is radiosity: the facades within `radius`
// of the microphone are cut into patches, each patch re-emits what reaches it
// from the plane wave and from the other patches, and the microphone gathers
// what reaches it from all of them.
//
// Two simplifications set the model:
//
//   - A patch is lit by the plane wave and by every wave the floor and the
//     facades reflect specularly on the way to it, each carrying the product
//     of rho (1 - s) over its reflections (building_paths.hpp's specular
//     traces, from the patch): so every specular reflection, not just the
//     first, scatters its share.  Diffracted waves light nothing.
//   - Diffuse energy stays diffuse.  Energy a patch receives from other
//     patches it re-emits diffusely with its whole rho: a field arriving from
//     every direction reflects, specularly or not, into every direction.
//
// With B the patches' exitance, E the plane wave's irradiance, F the form
// factors (the irradiance at one patch per unit exitance of another, through
// the air and off the floor), S = diag(rho s) and D = diag(rho),
//
//   B = S E + D F B,   so the microphone gathers   g.B = w.E,
//   w = S (I - F^T D)^-1 g,
//
// with g the microphone's gathering from each patch.  Everything but E
// depends only on the geometry and the band, so the reverse solve for w is
// made once per microphone and band, and a direction costs only the test of
// which patches its wave lights.  The direct and the floor-reflected paths
// from a patch to the microphone interfere as the open ground's two rays do,
// band-averaged; between patches they add in energy.

#include <cstddef>
#include <vector>

#include "buildings.hpp"
#include "vec3.hpp"

namespace niceops {

struct DiffuseOptions
{
    // The facades within this horizontal distance of the microphone (m).
    double radius = 100.0;
    // The patches' size, at most (m).
    double patch_size = 4.0;
    // The reverse solve stops when no weight moves by more than this part of
    // the largest; it gives up after max_iterations.
    double tolerance = 1e-10;
    std::size_t max_iterations = 2000;
    // The specular waves that light a patch, as BuildingPathOptions traces
    // them for the microphone.
    std::size_t max_order = 64;
    std::size_t max_turns = 4;
    std::size_t max_departures = 4096;
};

class DiffuseFacades
{
public:
    struct Patch
    {
        Vec3 center;          // a millimetre out from the face
        Vec3 normal;          // unit, horizontal, outward
        double area = 0.0;
        double width = 0.0;   // along the facade
        double height = 0.0;
        Vec3 along;           // unit, the facade's base direction
        std::size_t facade = 0;
    };

    // `reflection` and `scattering` are [facade][band] for buildings.facades();
    // `floor_reflection` is the floor's (specular) energy reflection.  The
    // microphone is in the local frame, z above the floor.  The form factors
    // are shared out over `threads`.
    DiffuseFacades(const Buildings& buildings, const Vec3& microphone, const std::vector<double>& reflection,
            const std::vector<double>& scattering, double floor_reflection, const std::vector<double>& bands,
            double sound_speed, const DiffuseOptions& options, unsigned threads = 1);

    // The diffuse energy at the microphone per band, relative to the plane
    // wave's intensity, for the wave from `u` (unit, toward the source).
    void energy(const Vec3& u, double* out) const;

    // The same by the forward solve, B = S E + D F B per band: the reference the
    // reverse solve is tested against, and far slower.
    void energy_forward(const Vec3& u, double* out) const;

    // E: the irradiance on every patch per band, [patch][band], relative to
    // the plane wave's intensity.
    void irradiance(const Vec3& u, std::vector<double>& e) const
    {
        irradiance(u, all_, e);
    }

    const std::vector<Patch>& patches() const
    {
        return patches_;
    }

    // The patch's form factors: (other patch, F) pairs, F being the irradiance
    // at this patch per unit exitance of the other.
    const std::vector<std::pair<std::size_t, double>>& form_factors(std::size_t patch) const
    {
        return form_factors_[patch];
    }

    // Whether every band's reverse solve met its tolerance.
    bool converged() const
    {
        return converged_;
    }

private:
    // The irradiance on `which` patches per band, [patch][band]: each specular
    // arrival's cos(incidence) times its reflections' rho (1 - s) and the
    // floor's reflection.
    void irradiance(const Vec3& u, const std::vector<std::size_t>& which, std::vector<double>& e) const;

    const Buildings* buildings_;
    DiffuseOptions options_;
    std::size_t bands_ = 0;
    double floor_reflection_ = 1.0;
    std::vector<double> specular_;    // [facade][band], rho (1 - s)
    std::vector<std::size_t> all_;    // every patch
    std::vector<std::size_t> planes_; // the orientations within reach
    std::vector<Patch> patches_;
    std::vector<std::vector<std::pair<std::size_t, double>>> form_factors_;
    std::vector<double> gather_;      // [patch][band]
    std::vector<double> rho_;         // [patch][band]
    std::vector<double> rho_s_;       // [patch][band]
    std::vector<double> weight_;      // [patch][band], w
    std::vector<std::size_t> active_; // patches whose weight is not negligible
    bool converged_ = true;
};

} // namespace niceops
