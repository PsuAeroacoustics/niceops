#pragma once

// --buildings: the point-source urban model in a run (THEORY 6.5).  Every
// microphone hears the source's own rays among LoD1 buildings
// (point_paths.hpp) -- reflected in the facades and the floor, diffracted at
// the roof edges and corners once or at two in turn, a facade between the two
// or not, or at three roof edges -- each reading the
// sphere along the direction it leaves the source (point_ray_launch), not
// along the blocked straight line, and summed over three sub-bands as
// niceops_point_field sums them, each pair's cross term weighted by how
// coherent the source is taken to be between their launch directions
// (source_coherence.hpp): by default, coherent within a cone and in energy
// beyond.
// A ray keeps the sphere's
// conventions (eval_spectrum): it spreads from the sphere's radius, so one
// shorter than that gains nothing, and absorbs over its length beyond it.
// With no buildings and a floor that absorbs everything, the one ray left is
// the straight path, and the level is the run's without a ground model.
//
// Walls and roofs that scatter reflect specularly only rho (1 - s) of what
// reaches them, and the rest of rho diffusely, as niceops_point_field takes
// them: each microphone's radiosity field (diffuse_facades.hpp), lit by the
// source's specular paths to each patch, each path carrying the sphere along
// its own launch, adds its energy to the rays' in each band.
//
// The district is opaque here: the run's translation unit cannot see
// buildings.hpp, whose Ring is not contour.hpp's.

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "source_coherence.hpp"
#include "sphere_database.hpp"
#include "vec3.hpp"

namespace niceops {

namespace gpu {
class Device;
}

class PointSourceField;

// The q-quantile of the gamma distribution with mean `mean` and variance
// `variance`: what the interval's ends take (PhaseInterval).
double gamma_quantile(double mean, double variance, double q);

// The interval over the source's unknown phase relations: `percent` of it
// (90 for the 5th to the 95th percentile), and the levels at its ends,
// [state][mic] A-weighted dB.
struct PhaseInterval
{
    double percent = 90.0;
    std::vector<double> lower, upper;
};

// The band energies of a field's rays partially coherent (THEORY 6.5): each
// pair's cross term, band-averaged over its sub-band as band_averaged_energy
// forms it, times the mean coherency of their launch directions
// (point_ray_launch): the cone's, or the measured model's with the main
// rotor's delay added to its path difference in the bands it holds -- about
// `frame`'s normal, from its nose.  The
// arguments up to `weights` are PointSourceField::band_levels's, `bands`
// starting at the database's band `first_band`.  `energy` gets three per
// band, as band_level() takes them.  `variance`, when given, gets each band's
// variance of its band-averaged energy over the source's unknown phase
// relations: rays whose mean coherency at the band's centre is 0.99 or more,
// by single linkage, deviate together, so each pair of such groups adds half
// the square of its summed cross terms, each weighted by the root of its
// pair's spread.  Under RaySum::coherent every pair's coherency is 1 at no
// delay; under RaySum::energy no pair adds.
void partially_coherent_energies(const PointSourceField& field, const std::vector<double>& reflection,
        double floor_reflection, const std::vector<double>& bands, std::size_t first_band, double sound_speed,
        const std::vector<double>& air, const double* weights, const SourceCoherence& coherence,
        const SphereFrame& frame, double* energy, double* variance = nullptr);

class UrbanDistrict
{
public:
    // The run's frame and units, and the coefficients where a building gives
    // none: a CSV of footprint vertices or a GeoJSON projected about `origin`
    // (latitude, longitude; empty for none) with +x along `frame_bearing`.
    // Throws std::runtime_error, its message naming the problem.
    struct Input
    {
        std::string units = "ft";
        double to_meters = 0.3048;
        std::vector<double> origin;
        double frame_bearing = 90.0;
        double wall_absorption = 0.0;
        double floor_absorption = 0.0;
        // The share of what the walls and roofs reflect that they scatter
        // diffusely, where a building gives none, and a CSV of band_hz and s
        // that multiplies it band by band (read_scattering_spectrum; "" for
        // none).
        double scattering = 0.0;
        std::string scattering_file;
        // The diffuse field's patches' size at most, and the horizontal
        // distance from a microphone within which a building scatters, in
        // metres; 0 for the rays' own reach (PointPathOptions::radius).
        double patch_size = 4.0;
        double diffuse_radius = 0.0;
        // A shared face's aperture at every wavenumber, not sampled
        // (PointPathOptions::sample_shared_faces).
        bool exact_shared_faces = false;
    };
    static std::shared_ptr<const UrbanDistrict> read(const std::string& file, const Input& input);
    ~UrbanDistrict();

    std::size_t building_count() const;
    std::size_t facade_count() const;
    // Whether any building's coefficient is above 0, so that its walls and
    // roof scatter wherever the spectrum does.
    bool scatters() const;
    // The scattering spectrum's factor at each of `bands` (1 without one).
    // Throws std::runtime_error where the spectrum does not reach a band.
    std::vector<double> band_scattering(const std::vector<double>& bands) const;
    // The diffuse field's patch size and reach, metres.
    double patch_size() const;
    double diffuse_radius() const;
    // "" where a point may stand, or why not: "below the floor", "inside
    // building <id>".
    std::string misplaced(const Vec3& point) const;

    // Each state's level at every microphone, from the states' sources (their
    // positions, sphere frames and blends) and the run's absorption (dB/m per
    // database band) and gains (a_weighting_gains): `levels` gets
    // [state][mic] A-weighted dB, -infinity where no ray arrives, and
    // `band_levels`, when given, [state][mic][band] each band's unweighted
    // level, dB re 20 uPa, over every database band, those outside
    // [first_band, last_band) -infinity.  `packed`, when given, is the run's
    // packed spectra per state.  `interval`, when given, gets each level's
    // interval over the source's phases (the gamma distribution with the
    // A-weighted energy's mean and variance); -infinity where no ray arrives.
    //
    // With a `state_tolerance` (dB), not every state is evaluated: per
    // microphone, the buildings' transfer -- each band's level less the direct
    // ray's in free field from the sphere at that state -- is taken at states
    // 32 apart and at the middle of any stretch where it taken linearly
    // between the stretch's ends misses the level there by more than the
    // tolerance, or a band within 30 dB of the loudest A-weighted one by 15
    // times it; the other states are their direct ray's free-field levels plus
    // the transfer interpolated (THEORY 6.5).  0 evaluates every state.
    // Returns the pairs evaluated.
    //
    // With a `device`, each pair's rays are formed and summed on it
    // (gpu/point_rays.hpp): for rays adding coherently, in a cone or in
    // energy, without an interval -- gpu_supports() says.
    std::size_t levels(const Database& database, const std::vector<Database::SourceState>& states,
            const std::vector<Vec3>& microphones, const double* absorption, const std::vector<double>& gains,
            std::size_t first_band, std::size_t last_band, double sound_speed, unsigned threads,
            const SourceCoherence& coherence, std::vector<double>& levels,
            std::vector<double>* band_levels = nullptr,
            const std::vector<const Database::PackedSpectra*>* packed = nullptr,
            PhaseInterval* interval = nullptr, double state_tolerance = 0.0, gpu::Device* device = nullptr) const;
    // Whether levels() can take a GPU for these rays' sum, and an interval or
    // not: the measured coherence and the interval are the CPU's.
    static bool gpu_supports(const SourceCoherence& coherence, bool interval);

private:
    UrbanDistrict();
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace niceops
