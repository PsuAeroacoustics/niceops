#pragma once

// --buildings: the point-source urban model in a run (THEORY 6.5).  Every
// microphone hears the source's own rays among LoD1 buildings
// (point_paths.hpp) -- reflected in the facades and the floor, diffracted at
// the roof edges and corners once or at two in turn -- each reading the
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

class PointSourceField;

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
    };
    static std::shared_ptr<const UrbanDistrict> read(const std::string& file, const Input& input);
    ~UrbanDistrict();

    std::size_t building_count() const;
    std::size_t facade_count() const;
    // Whether any building gave a scattering coefficient, which a run leaves
    // out.
    bool scatters() const;
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
    void levels(const Database& database, const std::vector<Database::SourceState>& states,
            const std::vector<Vec3>& microphones, const double* absorption, const std::vector<double>& gains,
            std::size_t first_band, std::size_t last_band, double sound_speed, unsigned threads,
            const SourceCoherence& coherence, std::vector<double>& levels,
            std::vector<double>* band_levels = nullptr,
            const std::vector<const Database::PackedSpectra*>* packed = nullptr,
            PhaseInterval* interval = nullptr) const;

private:
    UrbanDistrict();
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace niceops
