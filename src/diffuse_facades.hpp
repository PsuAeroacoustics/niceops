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
#include <functional>
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
    // The images of a patch a point source's specular paths may reach it
    // through (point_energy), at most.
    std::size_t max_images = 20000;
    // The air's attenuation per band, Np/m of amplitude, on the paths from the
    // patches to the microphone; empty for none.  Between patches it is left
    // out.
    std::vector<double> air;
    // Whether the roofs scatter too, cut into patches facing up, each with
    // the rho and s of its building's first facade, as a point source's roof
    // reflection takes them.  A building then scatters as a whole, walls and
    // roof, where any of its walls is within the radius.  For a model whose
    // rays reflect in the roofs.
    bool roofs = false;
};

class DiffuseFacades
{
public:
    struct Patch
    {
        Vec3 center;          // a millimetre out from the face
        Vec3 normal;          // unit, outward: horizontal on a wall, up on a roof
        double area = 0.0;
        double width = 0.0;   // along `along`
        double height = 0.0;  // along `across`
        Vec3 along;           // unit, the facade's base direction
        Vec3 across {0.0, 0.0, 1.0}; // unit, in the face across `along`
        // The facade it lies on, or for a roof (`roof`) the facade whose rho
        // and s the roof takes (Buildings::roof_facade).
        std::size_t facade = 0;
        bool roof = false;
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

    // A source that is not the same in every direction: per band, the power
    // it sends along a lighting path, relative to the unit source's, given
    // the path's launch (unit, from the source) and its length.
    using SourcePower = std::function<void(const Vec3& launch, double length, double* power)>;

    // The same for a point source (point_paths.hpp): a unit source, whose
    // free field is |p|^2 = 1/r^2, lights each patch by every specular path
    // from it to `max_order` reflections (the patch's images, beam traced as
    // the microphone's are), each cos(incidence) / L^2 times its reflections'
    // rho (1 - s) and the floor's reflection, and exp(-2 a L) of the air
    // (`air`, Np/m of amplitude per band; empty for none).  Diffracted waves
    // light nothing, as for the plane wave.  The energy is in |p|^2, to add to
    // the rays' band energy.  The patches are shared out over `threads`.
    // With `power`, each path also carries the source's power along its
    // launch, which is the straight line to the patch's image.
    void point_energy(const Vec3& source, std::size_t max_order, const std::vector<double>& air, double* out,
            unsigned threads = 1, const SourcePower* power = nullptr) const;
    // The irradiance point_energy gathers, [patch][band].
    void point_irradiance(const Vec3& source, std::size_t max_order, const std::vector<double>& air,
            std::vector<double>& e, unsigned threads = 1, const SourcePower* power = nullptr) const;
    // Trace each active patch's images to `max_order` once, for every source
    // point_energy is then given at that order: what a moving source would
    // otherwise trace again at every position.  Not thread-safe against
    // point_energy; call it before.
    void prepare_point_lighting(std::size_t max_order, unsigned threads = 1);

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

    // Whether a facade is within the radius, cut into patches: one whose
    // scattered share this field returns.  A facade outside it reflects its
    // whole rho specularly.  With DiffuseOptions::roofs, every facade of a
    // building any of whose facades is in reach, and its roof with them.
    bool scatters(std::size_t facade) const
    {
        return facade < in_reach_.size() && in_reach_[facade] != 0;
    }

private:
    // The roofs' patches (DiffuseOptions::roofs), once the walls' reach is
    // known.
    void add_roof_patches();
    // What the microphone gathers from an irradiance [patch][band].
    void gathered(const std::vector<double>& e, double* out) const;
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
    std::vector<char> in_reach_;      // [facade], within the radius
    bool converged_ = true;
    // prepare_point_lighting's images of each active patch (in active_'s
    // order), and their order.
    struct PatchImage
    {
        Vec3 point;
        std::vector<std::size_t> facades;
    };
    std::vector<std::vector<PatchImage>> patch_images_;
    std::size_t patch_image_order_ = 0;
    // The images of active patch `a` to `max_order`: the prepared ones, or
    // traced into `scratch`.
    const std::vector<PatchImage>& patch_images(std::size_t a, std::size_t max_order,
            std::vector<PatchImage>& scratch) const;
};

} // namespace niceops
