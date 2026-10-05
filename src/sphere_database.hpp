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
#include <cstdint>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include <complex>

#include "ground_effect.hpp"
#include "plate_receiver.hpp"
#include "propagation.hpp"
#include "util/constants.hpp"
#include "vec3.hpp"

#include "condition_selection.hpp"
#include "interpolator.hpp"

namespace niceops {

// Which velocity a database's flight conditions -- mu, the path angle and the
// sphere's azimuth reference -- were taken from: the aircraft's over the
// ground, or its through the air.  Recorded at the root of the database as the
// text attribute speed_reference; see Database::speed_reference().
enum class SpeedReference
{
    ground,
    air,
};

// "ground" or "air"; nullopt for anything else.
std::optional<SpeedReference> parse_speed_reference(const std::string& name);
const char* to_string(SpeedReference reference);


// Bilinear interpolation on the sphere grid.  A cell whose knots come in
// decreasing order wraps the seam of a periodic axis: with `x_period` set, its
// upper knot is one period past x2 and x is read into the cell.  Without a
// period the knots are simply swapped, as before.
double periodic_bilinear_interp(double q11, double q12, double q21, double q22,
        double x1, double x2, double y1, double y2, double x, double y,
        double x_period = 0.0);

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

    // Make the axis periodic (phi: 360 degrees).  The cell across the seam then
    // runs from the last knot up to the first knot one period on, and a value
    // past either end is read into the axis before it is located.  Without a
    // period the seam is never crossed.
    void set_period(double period);

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
    // Period of a periodic axis, or 0 when the axis does not wrap.
    double period_ = 0.0;
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
// direction.  phi runs over the whole circle: a direction above the tip-path
// plane reads the database's upper half.  `fold_upper` reads the mirror image
// instead, for a database that stores no upper half (see
// Database::upper_hemisphere_measured()).
std::array<double, 2> lookup_angles(const SphereFrame& frame, const Vec3& direction,
        bool fold_upper = false);

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
// This is what lets the observer loop find the peak and then the 10-dB-down
// interval without evaluating every emission time.  The interval runs from the
// first to the last sample within 10 dB of the peak (ten_db_down_interval), so
// the bound serves twice: to be sure the peak found is the real one, and to
// skip states it proves are below the threshold while scanning in from either
// end for the interval's first and last samples.
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

    // Which bin a propagation direction falls in, from the *unnormalized* dot
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
        return numerators_[blend_of_state_[state] * bin_count + bin_index] * state_scale_[state];
    }

private:
    friend class Database;
    std::vector<double> numerators_;
    std::vector<std::size_t> blend_of_state_;
    // 10^(level_correction_db/10) per state.  The fixed_load_factor correction
    // shifts every direction and band alike, so it factors out of the table;
    // keeping it in the grouping key made nearly every state its own group,
    // since it varies continuously with CT.
    std::vector<double> state_scale_;
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
    // needs neither the angle nor a normalized direction.
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
    // expensive part of an evaluation -- it is a Delaunay locate() walk --
    // and it depends only on (mu, gamma, CT), so it is hoisted out of any loop
    // that holds the operating condition fixed.
    //
    // A database built with panam's load_factors=None stores only the LF=1
    // reference condition at each (mu, gamma): CT carries no information (it
    // is the same constant on every stored sphere), so the triangulation this
    // blends against is 2D, over (mu, gamma) alone; weights[3]/indices[3] is
    // then padding (weight 0, matching the existing convention for a
    // degenerate blend), and level_correction_db carries what CT alone would
    // have contributed -- see the fixed_load_factor handling in the .cpp.
    struct Blend
    {
        std::array<double, condition_dims + 1> weights;
        std::array<std::size_t, condition_dims + 1> indices;
        // Additive dB correction for a query load factor away from what the
        // database stores, from panam's own scaling law (dBA, and every band
        // of a spectral database, are exactly linear in 20*log10(load_factor)
        // -- see add_sphere_group in panam's flight_acoustics.py).  Zero for
        // an ordinary database.  Under the default one-sided correction it is
        // also zero for a query at or below the reference, and only ever adds
        // level, matching panam's own convention of leaving load factors at
        // or below 1 g unscaled from the stored reference; under
        // LoadFactorCorrection::two_sided it is negative there, down to the
        // floor BlendOptions::min_load_factor sets.  Either way it is one
        // uniform shift of every level, which is what keeps the bounds exact.
        double level_correction_db = 0.0;
        // Queries outside the stored axis ranges are clamped before blending.
        // These flags describe axis coverage, not coverage of the convex hull
        // or the measured angular directions.
        bool advance_ratio_clamped = false;
        bool flight_path_angle_clamped = false;
        bool thrust_coefficient_clamped = false;
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

    // How blend_at() picks and weights its spheres, and how it applies a
    // fixed_load_factor database's law (condition_selection.hpp).  Declared
    // there, at namespace scope: a nested struct with member initializers
    // cannot be a defaulted argument of its enclosing class's members.
    using BlendOptions = niceops::BlendOptions;

    explicit Database(const std::string& datafile_path);

    Blend blend_at(double advance_ratio, double flight_path_angle,
            double thrust_coefficient, const BlendOptions& options = {}) const;

    SourceState source_state(const Vec3& source, const Vec3& velocity,
            const Vec3& normal, double advance_ratio, double flight_path_angle,
            double thrust_coefficient, const BlendOptions& options = {}) const;

    double eval(double range, double phi, double theta, const Blend& blend) const;

    double eval(double range, double phi, double theta, double advance_ratio,
            double flight_path_angle, double thrust_coefficient) const;

    double level_at_observer(const Vec3& observer, const SourceState& state) const;

    // Ground level from geometry that has already been formed: the squared
    // range, and the three *unnormalized* dot products of the displacement with
    // the sphere frame.  The observer loop's bound test computes exactly these,
    // so sharing them is what makes skipping states pay -- a state that
    // survives the test costs no more than it would have without one.
    double level_from_geometry(double squared_range, double along, double up, double across,
            const Blend& blend) const;

    // Whether every condition group records which channels the builder
    // measured: a `coverage` variable, 1 for a measured channel and 0 for one
    // the gate zeroed or that was filled from elsewhere.  The share below is
    // only defined when it does.
    bool has_coverage() const
    {
        return has_coverage_;
    }

    // The share of a direction's blended level that comes from what the
    // builder did not measure: the bilinear weight on each corner, times the
    // blend's weight on each condition, times the unmeasured share of that
    // cell (with spectra, of its A-weighted band energy at the sphere; with
    // channel coverage, all or nothing).  `angles` are (phi, theta) in
    // degrees, as lookup_angles() gives them.
    double unmeasured_share(const std::array<double, 2>& angles, const Blend& blend) const;

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
    // dB/m coefficient for each band, from absorption_coefficients(), or the
    // mean over the heights the path spans (LayeredAbsorption).
    //
    // `packed`, when given, is pack_spectra(blend.indices, ...): the same
    // values gathered per blend.  The result is identical to the bit; it only
    // changes where they are read, and lets a stencil with nothing empty take
    // a branch-free loop across the bands.
    struct PackedSpectra;
    void eval_spectrum(double range, double phi, double theta, const Blend& blend,
            const double* absorption, double* levels, std::size_t first_band = 0,
            std::size_t last_band = static_cast<std::size_t>(-1),
            const PackedSpectra* packed = nullptr) const;
    // The same, with spreading taken over `range` and absorption over
    // `path_length - r_sphere`: a curved ray's tube and its arc (RayPath).  The
    // form above is this one with the two equal.
    void eval_spectrum(double range, double path_length, double phi, double theta,
            const Blend& blend, const double* absorption, double* levels,
            std::size_t first_band = 0, std::size_t last_band = static_cast<std::size_t>(-1),
            const PackedSpectra* packed = nullptr) const;

    // One blend's spectra, laid out [cell][k][band] for the k-th of its four
    // indices, with a bit mask per (cell, k) of the bands in [first_band,
    // last_band) that are empty there (databases with more than 64 bands
    // skip the fast path).  spectra_ interleaves every condition, so a
    // blend's four values for one (cell, band) sit on four cache lines of a
    // table tens of MB long; here each (cell, k) is one contiguous run of
    // bands, which is what the vectorized loop in eval_spectrum reads.
    // Size: cells * 4 * bands doubles.
    struct PackedSpectra
    {
        std::vector<double> values;
        std::vector<std::uint64_t> empty_bands;
        std::size_t first_band = 0;
        std::size_t last_band = 0;
    };
    PackedSpectra pack_spectra(const std::array<std::size_t, condition_dims + 1>& indices,
            std::size_t first_band = 0,
            std::size_t last_band = static_cast<std::size_t>(-1)) const;

    // A-weighted, band-limited ground level at an observer, propagated per
    // band.  With `plate` (which needs `ground` and `impedance`), the observer
    // is a flush microphone on a ground plate whose offset points along
    // `plate_mic_azimuth` (radians in the local frame, from +x toward +y), and
    // the plate's level replaces the point-receiver ground effect.  `gains` comes from a_weighting_gains() and must match the band
    // count.  Returns -infinity when no band is inside the limits.
    double spectral_level_at_observer(const Vec3& observer, const SourceState& state,
            const double* absorption, const std::vector<double>& gains,
            const GroundModel* ground = nullptr,
            const std::vector<std::complex<double>>* impedance = nullptr,
            double sound_speed = 343.2, std::size_t first_band = 0,
            std::size_t last_band = static_cast<std::size_t>(-1),
            const PackedSpectra* packed = nullptr, const PlateReceiver* plate = nullptr,
            double plate_mic_azimuth = 0.0, const GroundTable* ground_table = nullptr,
            double* band_levels = nullptr) const;

    // The same along a path the caller has already formed (propagation.hpp):
    // spectral_level_at_observer is this with straight_ray().  `extra_loss`,
    // when given, is dB per band taken off after propagation -- the shadow
    // zone's, under --refraction.  `ground_weight`, when given, scales the
    // ground effect per band in dB (0 drops it): ShadowTransition hands the
    // ground over to its own field near a shadow boundary.
    //
    // `band_levels`, when given, receives each band's unweighted level at the
    // observer, dB re 20 uPa, with every term the metric sums: spreading,
    // absorption, `extra_loss` and the ground or plate.  Only bands in
    // [first_band, last_band) are written (--band_history).
    double spectral_level_along(const RayPath& path, const SourceState& state,
            const double* absorption, const std::vector<double>& gains,
            const GroundModel* ground = nullptr,
            const std::vector<std::complex<double>>* impedance = nullptr,
            double sound_speed = 343.2, std::size_t first_band = 0,
            std::size_t last_band = static_cast<std::size_t>(-1),
            const PackedSpectra* packed = nullptr, const PlateReceiver* plate = nullptr,
            double plate_mic_azimuth = 0.0, const GroundTable* ground_table = nullptr,
            const double* extra_loss = nullptr, const double* ground_weight = nullptr,
            double* band_levels = nullptr) const;

    // Band centers carried by the database, empty when it has no spectra.
    const std::vector<double>& band_centers() const
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
        double radius = 0.0; // dimensional, meters
        // Mechanical shaft power for this condition, when the database carries
        // it.  NaN means the database has none.
        double shaft_power = std::numeric_limits<double>::quiet_NaN();
        // What the database records at its root and the loader interprets, carried
        // so an exported sphere reloads as what it was cut from: its flight
        // conditions' reference (air or ground), and the rotor it was built for.
        // NaN where the source had none.
        SpeedReference speed_reference = SpeedReference::ground;
        double tip_speed = std::numeric_limits<double>::quiet_NaN();
        double rotor_radius = std::numeric_limits<double>::quiet_NaN();
        double vehicle_weight_newtons = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> phis;
        std::vector<double> thetas;
        // Row-major over (phi, theta): cells[i_phi * thetas.size() + i_theta].
        std::vector<Sample> cells;
        // Band centers, and the per-band level at every cell laid out as
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

    // The rotor radius (m) and the weight (N) the database was built for, from
    // the root variables main_rotor_radius_meters and vehicle_weight_newtons, or
    // NaN when the file does not carry them.  With the tip speed these are what
    // the vehicle file has to agree with for mu and CT to index the right
    // conditions; only the tip speed used to be compared.
    double rotor_radius() const
    {
        return rotor_radius_;
    }
    double vehicle_weight() const
    {
        return vehicle_weight_newtons_;
    }

    // Which velocity the flight conditions were taken from, as the root
    // attribute speed_reference records it.  Absent means the ground, which is
    // what every database built before the attribute existed used.
    // True when every condition group stores channels above the tip-path
    // plane (|phi| > 90 degrees).  Then a direction above the plane reads that
    // stored half; otherwise it reads the mirror image, as before.
    bool upper_hemisphere_measured() const
    {
        return upper_hemisphere_measured_;
    }

    // The atmosphere the broadband EAA was built in, and the rotor speed the
    // spheres were built at, when the builder recorded them (panam does; older
    // databases carry neither).  Absent means the database says nothing.
    struct BuildAtmosphere
    {
        double temperature_kelvin;
        double pressure_kpa;
        double relative_humidity_percent;
    };
    const std::optional<BuildAtmosphere>& build_atmosphere() const
    {
        return build_atmosphere_;
    }
    const std::optional<double>& rotor_rpm() const
    {
        return rotor_rpm_;
    }

    SpeedReference speed_reference() const
    {
        return speed_reference_;
    }

    // Whether the database carries per-band spectra alongside the overall level.
    bool has_spectra() const
    {
        return !frequencies_.empty();
    }

    // Whether this is a fixed_load_factor build, and the thrust coefficient of
    // its reference condition (load factor 1).  Together they fix the analytic
    // load-factor law, 20 lg(C_T / C_T,ref).
    bool fixed_load_factor() const
    {
        return fixed_load_factor_;
    }
    double reference_thrust_coefficient() const
    {
        return reference_thrust_coefficient_;
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
            double thrust_coefficient, const BlendOptions& options = {}) const;

    // How many thrust-coefficient layers BlendPolicy::aam searches: one on a
    // fixed_load_factor database (the law carries C_T), otherwise one per
    // distinct stored C_T, and a state takes the nearest layer.
    std::size_t aam_layer_count() const
    {
        return aam_layers_.size();
    }

    // Range beyond which a source at this operating condition provably falls
    // below `cutoff` dBA, or infinity when no such bound can be established.
    double cull_range(const Blend& blend, double cutoff) const;

    // Directional level bounds for a trajectory's states.  One table per
    // *distinct* blend rather than per state -- a shipped track resolves to a
    // few dozen tetrahedra over thousands of states -- so the whole thing is a
    // few tens of kilobytes and stays in cache under the observer loop.
    DirectionBounds direction_bounds(const std::vector<SourceState>& states) const;

    // The same tables for the spectral path's A-weighted, band-limited level
    // (spectral_level_at_observer with these `gains` and no ground effect),
    // under the same contract: numerator >= 10^(L/10) r^2 outside the sphere.
    // Dropping the per-band absorption is what makes it a bound, and ISO 9613-1
    // absorption is never negative, so unlike direction_bounds() it does not
    // depend on the EAA factors.  The ground effect can add up to ~6 dB and is
    // not covered: a caller propagating with it must not use these.
    DirectionBounds spectral_direction_bounds(const std::vector<SourceState>& states,
            const std::vector<double>& gains) const;

    double sphere_radius() const
    {
        return sphere_radius_;
    }

    // dBA of -inf legitimately means "no energy in this direction" (ambient
    // gating emptied it), but a zero barycentric weight against -inf is a
    // NaN under IEEE 754 that poisons an entire blended sum regardless of
    // weight -- see the constructor's grid-fill comment. Every non-finite
    // dBA (and spectral band) read from a file is floored to this value once,
    // at load time, so self_test() must apply the same floor to the raw value
    // it re-reads before comparing it against eval().  The blends treat an
    // entry at the floor as "not measured" and leave it out (blend_present in
    // sphere_database.cpp), so it only comes back where every entry is empty.
    // Public so those file-local blends can name it.
    static constexpr double no_energy_floor_db = -1000.0;

private:

    // What spectral_level_at_observer and spectral_level_along share; see the
    // definition for why it is forced inline.
    double spectral_level_impl(const RayPath& path, const SourceState& state,
            const double* absorption, const std::vector<double>& gains,
            const GroundModel* ground, const std::vector<std::complex<double>>* impedance,
            double sound_speed, std::size_t first_band, std::size_t last_band,
            const PackedSpectra* packed, const PlateReceiver* plate, double plate_mic_azimuth,
            const GroundTable* ground_table, const double* extra_loss,
            const double* ground_weight = nullptr, double* band_levels = nullptr) const;

    // Group states by the conditions they blend, filling the bounds' per-state
    // group and load-factor scale.  Both direction_bounds() and
    // spectral_direction_bounds() build one table per group.
    void group_blends(const std::vector<SourceState>& states, DirectionBounds& bounds,
            std::vector<std::array<std::size_t, condition_dims + 1>>& blends) const;

    // The DirectionBounds bin a grid cell's center direction falls in, and how
    // far a bilinear stencil reading that cell can reach in bins either way.
    int home_bin(std::size_t i_phi, std::size_t i_theta) const;

    // The bins a cell is filed under before padding: its home bin, as {first,
    // count} of consecutive bins -- one, except at the poles.  A pole cell's center
    // direction has no azimuth (up = across = 0), so home_bin gives it an
    // arbitrary one, while a query a degree off the axis can arrive from any
    // azimuth and weights the pole cells at up to 0.9.  Filing the cell under
    // every azimuth bin of its theta bin is what keeps the bound a bound there.
    std::array<int, 2> home_bins(std::size_t i_phi, std::size_t i_theta) const;
    std::array<int, 2> bin_padding() const;

    // Blend the four bracketing operating conditions at one grid cell.
    Sample blend_cell(std::size_t i_phi, std::size_t i_theta, const Blend& blend) const;
    // BlendPolicy::aam's selection into `blend` (weights and indices only).
    void blend_aam(double advance_ratio, double flight_path_angle, double thrust_coefficient,
            const BlendOptions& options, Blend& blend) const;

    GridAxis phis;
    GridAxis thetas;
    // Flat (phi, theta, condition) grid, indexed
    // (i_phi * thetas.size() + i_theta) * n_conditions + condition.  This was a
    // vector<vector<vector<double>>>, which put two dependent pointer loads in
    // front of every one of the 32 reads an eval makes -- the largest single
    // cost of a run once the Delaunay locate had been hoisted out of the loop.
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
    // From the root `fixed_load_factor` flag (absent or false in every
    // database predating panam's load_factors=None): every stored sphere is
    // the same LF=1 reference condition, so CT is constant across them
    // rather than a real axis, and thrust_coefficient_min/max above are not
    // used for it -- normalizing against them would divide by zero.  The
    // triangulation is then built over (mu, gamma) alone (dims 2, not
    // condition_dims), and blend_at() derives Blend::level_correction_db
    // from the query CT against reference_thrust_coefficient_ analytically
    // instead.
    bool fixed_load_factor_ = false;
    double reference_thrust_coefficient_ = 0.0;
    std::size_t active_condition_dims_ = condition_dims;
    // The triangulation over this database's own operating conditions.  Owning
    // it makes several Databases independent; it used to be a file-scope global,
    // so constructing a second one re-armed the interpolator under the first.
    ConditionInterpolator interpolator_;
    // The conditions as the triangulation sees them -- normalized, row-major,
    // active_condition_dims_ to a row -- for BlendPolicy::nearest.
    std::vector<double> normalized_conditions_;
    // BlendPolicy::aam's search: one AamSelector per stored thrust coefficient
    // (one in all on a fixed_load_factor database), over raw advance ratio and
    // flight path angle.
    struct AamLayer
    {
        double thrust_coefficient;
        AamSelector selector;
    };
    std::vector<AamLayer> aam_layers_;
    double max_level_ = 0.0;
    // Highest dBA on each operating condition's sphere, and the lowest excess
    // attenuation factor anywhere in the database.  Together they bound how
    // loud a blended state can be at a given range.
    std::vector<double> sphere_max_;
    // Band centers shared by every sphere, and the per-band levels laid out as
    // spectra_[(cell * frequencies_.size() + band) * n_conditions + condition].
    //
    // Condition innermost, matching grid above and for the same reason: the
    // observer loop blends four conditions at one (cell, band), so they want to
    // be on one cache line.  The netCDF file stores them the other way round and
    // the loader transposes once.  Both empty unless the database carries
    // spectra.
    double tip_speed_ = std::numeric_limits<double>::quiet_NaN();
    double rotor_radius_ = std::numeric_limits<double>::quiet_NaN();
    double vehicle_weight_newtons_ = std::numeric_limits<double>::quiet_NaN();
    SpeedReference speed_reference_ = SpeedReference::ground;
    bool upper_hemisphere_measured_ = false;
    std::optional<BuildAtmosphere> build_atmosphere_;
    std::optional<double> rotor_rpm_;
    // Per grid cell and condition, parallel to grid, the share of the cell
    // that was measured: with spectra, of its A-weighted band energy at the
    // sphere; with channel coverage, 1 or 0.  Empty without coverage.
    std::vector<float> measured_;
    bool has_coverage_ = false;
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
