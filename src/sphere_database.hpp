#pragma once

// Port of sphere_database.d.  Loads a netCDF noise database (.nod): one group
// per operating condition (mu advance ratio, gamma flight path angle, CT
// thrust coefficient), each holding scattered (phi, theta) channels of dBA
// level and EAA excess-attenuation factor.  eval() computes ground level:
// geometry -> sphere angles (lookup_angles), barycentric interpolation between
// operating conditions, periodic bilinear interpolation on the sphere grid,
// then spherical spreading and excess attenuation, EAA * (range - r_sphere) /
// 1000 -- EAA being the total A-weighted attenuation over the 1000 m reference
// distance the database was built at.  No ground effect: see --ground_effect,
// which is on the spectral path.  Only same_grid == true databases are
// supported (interpolation path is unimplemented).

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <complex>

#include "ground_effect.hpp"
#include "util/constants.hpp"
#include "vec3.hpp"

#include "interpolator.hpp"

namespace niceops {


// Check if points "wrap around" sphere, and swap inputs to bilinear
// interpolation if necessary
double periodic_bilinear_interp(double q11, double q12, double q21, double q22,
        double x1, double x2, double y1, double y2, double x, double y);

// Find the pair of grid indices bracketing val (wrapping around the ends)
std::array<std::size_t, 2> bracket(const std::vector<double>& a, double val);

// One axis of the sphere grid: the sorted unique angles, plus the spacing when
// the axis turns out to be uniform (every shipped database's is).  Knowing the
// spacing turns bracket() into a multiply and a compare rather than a binary
// search -- eval brackets twice, so a default run does that ~190 million times,
// and the searches were around a third of the evaluation loop.
class GridAxis
{
public:
    void reset(std::vector<double> values);

    const std::vector<double>& values() const
    {
        return values_;
    }

    std::size_t size() const
    {
        return values_.size();
    }

    double operator[](std::size_t i) const
    {
        return values_[i];
    }

    // Identical in every case to bracket(values(), val).
    std::array<std::size_t, 2> bracket(double val) const;

    // Where `val` sits on the axis: the same bracketing pair bracket() returns,
    // plus how far along the cell it lies.  `fraction` is the weight the
    // bilinear interpolation puts on `hi`, so the weight on `lo` is 1 -
    // fraction; that is exactly what periodic_bilinear_interp() derives from
    // the two knot values, wrap-around included, and
    // tests/test_sphere_database.cpp sweeps the two against each other.
    //
    // It exists because eval() wants the weights, not the knots: on a uniform
    // axis both indices and the fraction fall out of one multiply, where going
    // through bracket() meant four more loads from values_, four subtractions
    // and a divide by the cell area.
    struct Cell
    {
        std::size_t lo;
        std::size_t hi;
        double fraction;
    };

    Cell locate(double val) const;

private:
    std::vector<double> values_;
    double first_ = 0.0;
    // Reciprocal of the uniform spacing, or 0 when the axis is not uniform and
    // the binary search has to be used after all.
    double inv_step_ = 0.0;
    // Set only when every knot sits within a few ULP of where a perfectly
    // uniform axis would put it -- tighter than inv_step_'s tolerance, because
    // locate() derives the interpolation weight from the spacing rather than
    // from the knots and so has no chance to correct for a knot that is merely
    // close.  Every shipped database's axes are exact multiples of 10 and 5
    // degrees and take this path.
    bool exactly_uniform_ = false;
};

double clamp01(double x);

// acos with its argument clamped into [-1, 1], so round-off just past the ends
// does not produce NaN.
double acosc(double x);

// The sphere's orientation at one instant: the unit TPP normal, the "nose"
// (velocity projected into the TPP) and starboard.  None of it depends on
// where the observer is, so it is built once per trajectory state and reused
// across the whole observer grid.
struct SphereFrame
{
    Vec3 normal;    // n, unit TPP normal
    Vec3 nose;      // vn, velocity projected into the TPP, unit
    Vec3 starboard; // pn = vn x n
};

SphereFrame sphere_frame(const Vec3& velocity, const Vec3& normal);

// Sphere angles (phi, theta in degrees) of an already-normalized propagation
// direction.
std::array<double, 2> lookup_angles(const SphereFrame& frame, const Vec3& direction);

// Sphere angles (phi, theta in degrees) of the propagation direction
std::array<double, 2> lookup_angles(const Vec3& source,
        const Vec3& velocity, const Vec3& normal,
        const Vec3& observer);

// Repair the azimuth labels of a sphere written by a pre-versioned builder
// (database_version < 1), whose upper surface was completed with `-phi` instead
// of `180 - phi`.  Rewrites the mirrored half of `phis_scattered` in place and
// returns true if the sphere needed repair.  Leaves anything that does not match
// the broken shape untouched.
bool repair_legacy_hover_directivity(std::vector<double>& phis_scattered);

// Faithful port of D's std.math.approxEqual(value, reference, maxRelDiff,
// maxAbsDiff): relative term is against `reference`, defaults 1e-2 / 1e-5.
bool approx_equal(double value, double reference, double max_rel_diff = 1e-2,
        double max_abs_diff = 1e-5);

// An upper bound on the ground level a trajectory state can produce at an
// observer, good enough to prove that the state cannot hold the peak of that
// observer's time history -- without evaluating it.
//
// This is what lets the observer loop find the peak and then expand out to the
// 10-dB-down points, instead of marching through every emission time.  The SEL
// window is contiguous by construction (peak_range walks out from the peak and
// stops at the first sample below the threshold), so nothing outside it needs
// a bound at all; the bound is needed only to be sure the walk started on the
// right peak.
//
// The bound is `max sphere level in this direction - spherical spreading`.
// Dropping the excess attenuation is what makes it an upper bound, and that is
// only true when no excess attenuation factor in the database is negative --
// `usable` is false when one is, and the caller must then evaluate everything.
//
// Binning by direction is what makes it tight enough to be worth having.  The
// whole-sphere maximum is 11.4 dB above the level actually seen, at the median,
// which is *wider than the 10 dB window*: the peak search would then reach
// outside the window and there would be nothing to save.  Splitting the sphere
// into bins of (theta, azimuth) brings the median slack to 8.6 dB, under the
// window, and the peak search collapses inside it.
class DirectionBounds
{
public:
    // Bins along the sphere's theta axis, and azimuth sectors across it.
    static constexpr int theta_bins = 16;
    static constexpr int azimuth_bins = 8;
    static constexpr int bin_count = theta_bins * azimuth_bins;

    // Which bin a propagation direction falls in, from the *unnormalised* dot
    // products with the sphere frame and the squared range.
    //
    // Deliberately free of both a square root and any inverse trig: a state
    // this rejects is never evaluated, so everything spent here is overhead on
    // the rejection.  The azimuth test is scale-free already, and the theta bin
    // uses along^2/r^2, which needs one divide where the angle would need a
    // root and an acos.
    static int bin(double along, double up, double across, double squared_range);

    // The narrowest theta bin and the azimuth sector width, both in degrees --
    // what Database::direction_bounds pads the stencil by.  They live here with
    // bin() so that changing the binning scheme is one edit rather than three.
    static double narrowest_theta_bin_degrees();

    static constexpr double azimuth_bin_degrees()
    {
        return 180.0 / azimuth_bins;
    }

    bool usable() const
    {
        return usable_;
    }

    // 10^(bound/10) * r_sphere^2 for the state's blend in this direction, so
    // the test "bound >= level" is "numerator >= 10^(level/10) * range^2" and
    // the hot path needs no logarithm.
    double numerator(std::size_t state, int bin_index) const
    {
        return numerators_[blend_of_state_[state] * bin_count + bin_index];
    }

private:
    friend class Database;
    std::vector<double> numerators_;
    std::vector<std::size_t> blend_of_state_;
    bool usable_ = false;
};

namespace detail {

// tan of each azimuth sector edge, so the sector test in bin() is a multiply
// and a compare rather than an atan2.  An inline variable rather than a
// function-local static: bin() runs ~94 million times a run and a function
// static would test its guard on every one of them.
inline const std::array<double, DirectionBounds::azimuth_bins / 2 - 1> azimuth_edge_tangents = [] {
    std::array<double, DirectionBounds::azimuth_bins / 2 - 1> edges {};
    for (std::size_t k = 0; k < edges.size(); k++)
        edges[k] = std::tan(pi * static_cast<double>(k + 1) / DirectionBounds::azimuth_bins);
    return edges;
}();

} // namespace detail

// Defined here rather than in the .cpp because the observer loop calls it ~94
// million times from another translation unit: out of line it inlines only when
// NICEOPS_LTO is on and the toolchain supports it, and it is fifteen
// instructions that should never depend on a build option.
inline int DirectionBounds::bin(double along, double up, double across, double squared_range)
{
    // Theta: uniform in along^2, which is the square of the cosine, so no root.
    // The sign of `along` picks the hemisphere.
    constexpr int half_theta = theta_bins / 2;
    int theta_index = static_cast<int>((along * along / squared_range) * half_theta);
    if (theta_index > half_theta - 1)
        theta_index = half_theta - 1;
    theta_index = along >= 0.0 ? half_theta + theta_index : half_theta - 1 - theta_index;

    // Azimuth: |across| against |up| scaled by tan of each sector edge, which
    // needs neither the angle nor a normalised direction.
    constexpr int half_azimuth = azimuth_bins / 2;
    const double magnitude = std::fabs(across);
    const double reference = std::fabs(up);
    int sector = 0;
    while (sector < half_azimuth - 1
            && magnitude > reference * detail::azimuth_edge_tangents[sector])
        sector++;
    return theta_index * azimuth_bins
            + (across >= 0.0 ? half_azimuth + sector : half_azimuth - 1 - sector);
}

class Database
{
public:
    // The four database spheres bracketing one operating condition, and the
    // barycentric weights that blend them.  Computing this is by far the most
    // expensive part of an evaluation -- it is a CGAL Delaunay locate() walk --
    // and it depends only on (mu, gamma, CT), so it is hoisted out of any loop
    // that holds the operating condition fixed.
    struct Blend
    {
        std::array<double, condition_dims + 1> weights;
        std::array<std::size_t, condition_dims + 1> indices;
    };

    // Everything about one trajectory state that every observer shares: where
    // the aircraft is, how its sphere is oriented, and which database spheres
    // to blend.
    struct SourceState
    {
        Vec3 position;
        SphereFrame frame;
        Blend blend;
    };

    explicit Database(const std::string& datafile_path);

    Blend blend_at(double advance_ratio, double flight_path_angle,
            double thrust_coefficient) const;

    SourceState source_state(const Vec3& source, const Vec3& velocity,
            const Vec3& normal, double advance_ratio, double flight_path_angle,
            double thrust_coefficient) const;

    double eval(double range, double phi, double theta, const Blend& blend) const;

    double eval(double range, double phi, double theta, double advance_ratio,
            double flight_path_angle, double thrust_coefficient) const;

    double level_at_observer(const Vec3& observer, const SourceState& state) const;

    // Ground level from geometry that has already been formed: the squared
    // range, and the three *unnormalised* dot products of the displacement with
    // the sphere frame.  The observer loop's bound test computes exactly these,
    // so sharing them is what makes skipping states pay -- a state that
    // survives the test costs no more than it would have without one.
    double level_from_geometry(double squared_range, double along, double up, double across,
            const Blend& blend) const;

    double level_at_observer(const Vec3& observer,
            const Vec3& source, const Vec3& velocity,
            const Vec3& normal, double advance_ratio,
            double flight_path_angle, double thrust_coefficient) const;

    // ---- Spectral propagation -------------------------------------------
    //
    // The path above carries one broadband dBA from the sphere to the ground
    // and attenuates it with the database's own EAA factor, which was fitted
    // for whatever spectrum the source had.  These carry the spectrum itself:
    // each band is spread and absorbed at its own rate, and only then summed
    // into a metric.  The two answer different questions and will not agree --
    // see SpectralSetup in weighting.hpp for what the difference is made of.
    //
    // Both require has_spectra(); a database without spectra cannot use this
    // path at all.

    // Per-band ground levels in dB at one direction, written into `levels`,
    // which must have room for band_count() entries.  `absorption` is the
    // dB/m coefficient for each band, from absorption_coefficients().
    void eval_spectrum(double range, double phi, double theta, const Blend& blend,
            const double* absorption, double* levels, std::size_t first_band = 0,
            std::size_t last_band = static_cast<std::size_t>(-1)) const;

    // A-weighted, band-limited ground level at an observer, propagated per
    // band.  `gains` comes from a_weighting_gains() and must match the band
    // count.  Returns -infinity when no band is inside the limits.
    double spectral_level_at_observer(const Vec3& observer, const SourceState& state,
            const double* absorption, const std::vector<double>& gains,
            const GroundModel* ground = nullptr,
            const std::vector<std::complex<double>>* impedance = nullptr,
            double sound_speed = 343.2, std::size_t first_band = 0,
            std::size_t last_band = static_cast<std::size_t>(-1)) const;

    // Band centres carried by the database, empty when it has no spectra.
    const std::vector<double>& band_centres() const
    {
        return frequencies_;
    }

    double level_along_bearing(double range, double bearing, double altitude,
            const Vec3& velocity, const Vec3& normal,
            double advance_ratio, double flight_path_angle,
            double thrust_coefficient) const;

    double zero_range_level(double altitude, const Vec3& velocity,
            const Vec3& normal, double advance_ratio,
            double flight_path_angle, double thrust_coefficient) const;

    double find_range(double contour_level, double bearing, double altitude,
            const Vec3& velocity, const Vec3& normal,
            double advance_ratio, double flight_path_angle,
            double thrust_coefficient, double threshold, bool zero_range_check) const;

    void self_test(const std::string& datafile_path);

    // The dBA level and EAA excess-attenuation factor of one (phi, theta) grid
    // cell at one operating condition.  Kept adjacent because eval always reads
    // the pair, so one cache line serves both.
    struct Sample
    {
        double dBA;
        double EAA;
    };

    // A whole noise sphere at one operating condition: the grid axes and the
    // blended level and excess attenuation at every cell.  This is the same
    // blend eval() applies at a single direction, carried out over the entire
    // grid at once.
    struct Sphere
    {
        double advance_ratio = 0.0;
        double flight_path_angle = 0.0;
        double thrust_coefficient = 0.0;
        double radius = 0.0; // dimensional, metres
        // Mechanical shaft power for this condition, when the database carries
        // it.  NaN means the database has none.
        double shaft_power = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> phis;
        std::vector<double> thetas;
        // Row-major over (phi, theta): cells[i_phi * thetas.size() + i_theta].
        std::vector<Sample> cells;
        // Band centres, and the per-band level at every cell laid out as
        // amplitudes[cell * frequencies.size() + band].  Both empty when the
        // database carries no spectra.
        std::vector<double> frequencies;
        std::vector<double> amplitudes;
    };

    // Main rotor tip speed in m/s as recorded at the root of the database, or
    // NaN when the file does not carry it.  It turns a non-dimensional advance
    // ratio back into an airspeed without needing a matching vehicle file.
    double tip_speed() const
    {
        return tip_speed_;
    }

    // Whether the database carries per-band spectra alongside the overall level.
    bool has_spectra() const
    {
        return !frequencies_.empty();
    }

    // Whether the database carries a per-condition shaft power.
    bool has_shaft_power() const
    {
        return !shaft_power_.empty();
    }

    // Shaft power at a blended operating condition, or NaN when the database
    // has none.  Blended with the same barycentric weights as the levels.
    double shaft_power(const Blend& blend) const;

    // Blend the database down to one operating condition over the whole grid.
    Sphere interpolate_sphere(double advance_ratio, double flight_path_angle,
            double thrust_coefficient) const;

    // Range beyond which a source at this operating condition provably falls
    // below `cutoff` dBA, or infinity when no such bound can be established.
    double cull_range(const Blend& blend, double cutoff) const;

    // Directional level bounds for a trajectory's states.  One table per
    // *distinct* blend rather than per state -- a shipped track resolves to a
    // few dozen tetrahedra over thousands of states -- so the whole thing is a
    // few tens of kilobytes and stays in cache under the observer loop.
    DirectionBounds direction_bounds(const std::vector<SourceState>& states) const;

    double sphere_radius() const
    {
        return sphere_radius_;
    }


private:
    // Blend the four bracketing operating conditions at one grid cell.
    Sample blend_cell(std::size_t i_phi, std::size_t i_theta, const Blend& blend) const;

    GridAxis phis;
    GridAxis thetas;
    // Flat (phi, theta, condition) grid, indexed
    // (i_phi * thetas.size() + i_theta) * n_conditions + condition.  This was a
    // vector<vector<vector<double>>>, which put two dependent pointer loads in
    // front of every one of the 32 reads an eval makes -- the largest single
    // cost of a run once the CGAL locate had been hoisted out of the loop.
    std::vector<Sample> grid;
    std::size_t n_conditions = 0;
    double sphere_radius_ = 0.0;
    // 1 / sphere_radius_, so eval's spherical spreading term is a multiply.
    double inv_sphere_radius_ = 0.0;

    // Whether dropping the excess attenuation term bounds a level from above.
    // A level is sphere_level - 20*log10(r/r0) - EAA*(r - r0)/1000, so dropping
    // that term is conservative only if no excess factor in the database is
    // negative.  The D original took abs() of the factor precisely because it
    // could be.  Both cull_range() and direction_bounds() rest on this, so it
    // is asked in one place.
    bool excess_attenuation_is_a_loss() const
    {
        return min_excess_factor_ >= 0.0;
    }

    double advance_ratio_min = 0.0, advance_ratio_max = 0.0;
    double flight_path_angle_min = 0.0, flight_path_angle_max = 0.0;
    double thrust_coefficient_min = 0.0, thrust_coefficient_max = 0.0;
    // The triangulation over this database's own operating conditions.  Owning
    // it makes several Databases independent; it used to be a file-scope global,
    // so constructing a second one re-armed the interpolator under the first.
    ConditionInterpolator interpolator_;
    double max_level_ = 0.0;
    // Highest dBA on each operating condition's sphere, and the lowest excess
    // attenuation factor anywhere in the database.  Together they bound how
    // loud a blended state can be at a given range.
    std::vector<double> sphere_max_;
    // Band centres shared by every sphere, and the per-band levels laid out as
    // spectra_[(cell * frequencies_.size() + band) * n_conditions + condition].
    //
    // Condition innermost, matching grid above and for the same reason: the
    // observer loop blends four conditions at one (cell, band), so they want to
    // be on one cache line.  The netCDF file stores them the other way round and
    // the loader transposes once.  Both empty unless the database carries
    // spectra.
    double tip_speed_ = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> frequencies_;
    std::vector<double> spectra_;

    // Per-condition shaft power, empty when the database does not carry it.
    // Shaft power is a property of the flight condition, not of a direction on
    // the sphere, so it is one value per sphere rather than per channel.
    std::vector<double> shaft_power_;
    double min_excess_factor_ = 0.0;
};

// Write a single interpolated sphere as a one-group .nod database, the same
// format Database reads, so an exported sphere can be loaded straight back.
void write_sphere(const Database::Sphere& sphere, const std::string& path);

// Write an interpolated sphere in the AAM/RNM hemisphere format -- the format
// the source spheres themselves arrive in.  Requires spectra: the format's
// payload is AMPLITUDE(PHI, THETA, FREQUENCY), and there is nothing meaningful
// to put there when all the database holds is an overall A-weighted level.
void write_aam_sphere(const Database::Sphere& sphere, const std::string& path, double tip_speed);

} // namespace niceops
