#pragma once

// A horizontally stratified atmosphere -- temperature, humidity and pressure
// against height -- and the layered absorption it implies.
//
// ISO 9613-1's coefficient depends on all three, and between a helicopter at a
// few hundred meters and the ground they are not one value.  Integrating the
// coefficient along the path rather than taking it at the ground moved LAE by
// -0.6 to +1.3 dB(A) at 500 m in sonAIR's comparison (Zellmann & Wunderli,
// Inter-Noise 2014).  It is the first step of the propagation roadmap because
// it needs no ray model: along a straight path, what a band loses is the
// coefficient's mean over the heights the path spans, times the path length.
//
// Heights are above the observer plane -- the same z the trajectory uses -- so
// the ground is z = 0 whatever --ground_elevation says.
//
// Wind, when the file gives it, only matters to refraction (--refraction):
// absorption does not depend on it.

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "atmosphere.hpp"
#include "util/constants.hpp"

namespace niceops {

struct AtmosphereProfile
{
    // The file it came from, for messages.
    std::string source;
    // Meters above the observer plane, strictly ascending.
    std::vector<double> height;
    std::vector<double> temperature; // K
    std::vector<double> relative_humidity; // percent
    // kPa.  Empty when the file gave no pressure; with_hydrostatic_pressure()
    // fills it.
    std::vector<double> pressure;
    // m/s toward the east and toward the north (the meteorological u and v).
    // Empty when the file gives no wind: still air.
    std::vector<double> wind_east;
    std::vector<double> wind_north;

    bool has_pressure() const
    {
        return !pressure.empty();
    }
    bool has_wind() const
    {
        return !wind_east.empty();
    }
    double bottom() const
    {
        return height.front();
    }
    double top() const
    {
        return height.back();
    }

    // The atmosphere at height z, in meters.  Temperature and humidity are
    // linear between levels.  Pressure is linear in its logarithm, since it
    // falls off exponentially.  Every quantity is held at its end value outside
    // the levels.  Needs has_pressure().
    Atmosphere at(double z) const;

    // The wind at height z, m/s east and north, interpolated like temperature.
    // Zero without wind columns.
    std::array<double, 2> wind_at(double z) const;

    // A copy with the pressure filled in: the hydrostatic equation integrated
    // through the temperature profile from `surface_pressure` kPa at z = 0,
    // treating the air as dry.  Humidity would raise the virtual temperature,
    // which moves the pressure at 1 km by well under 0.1%, and absorption is
    // the least sensitive of the three inputs to pressure.
    AtmosphereProfile with_hydrostatic_pressure(double surface_pressure) const;
};

// Turns a wind's east and north components into the local frame, whose +x
// axis points along compass bearing `frame_bearing` degrees and +y ninety
// degrees to its left, as LocalFrame::track_aligned has it.  The ray models and
// the air-referenced flight conditions all take the wind through this one
// rotation.
struct CompassToLocal
{
    explicit CompassToLocal(double frame_bearing)
        : sin_b(std::sin(frame_bearing * pi / 180.0))
        , cos_b(std::cos(frame_bearing * pi / 180.0))
    {
    }

    std::array<double, 2> operator()(const std::array<double, 2>& east_north) const
    {
        return {east_north[0] * sin_b + east_north[1] * cos_b,
                -east_north[0] * cos_b + east_north[1] * sin_b};
    }

    // The inverse: a vector in the local frame as (east, north).
    std::array<double, 2> east_north(const std::array<double, 2>& local) const
    {
        return {local[0] * sin_b - local[1] * cos_b, local[0] * cos_b + local[1] * sin_b};
    }

    double sin_b;
    double cos_b;
};

// Read a profile.  Columns are found by name (read_labeled_csv), in any order,
// one row per level:
//   height above the observers   z, agl                   in `unit`; z_ft / z_m
//                                                          are checked against it
//   or altitude above sea level  alt, altitude, height    in `unit`, less
//                                                          `ground_elevation` (m)
//   temperature                  T_K, T_C                 or T / temperature,
//                                                          unit told from values
//   relative humidity, percent   RH, humidity
//   pressure (optional)          p_kPa, p_hPa, p_Pa       or p / pressure,
//                                                          unit told from values
//   wind (optional), either      wind_speed_mps / _kt     and wind_direction, the
//                                                          compass bearing it
//                                                          blows FROM, degrees
//                                or u_mps / u_kt, v_...   toward east and north
// Wind columns must name their unit: nothing about the values tells m/s from
// knots.  Other columns are ignored and named in the log.  A single row is a uniform
// atmosphere.  Rows may come in either height order.
AtmosphereProfile read_atmosphere_profile(const std::string& path, LengthUnit unit,
        double ground_elevation);

struct LabeledCsv;

// Where a profile file keeps each quantity: the column-finding half of
// read_atmosphere_profile, with the same names and unit rules, shared with the
// profile builder's reader (atmosphere_builder.hpp) so the two accept the same
// headers.  Every index is npos when the file has no such column.
struct ProfileColumns
{
    static constexpr std::size_t npos = static_cast<std::size_t>(-1);
    // thermo: a file without temperature or humidity is refused, as
    // --atmosphere refuses it.  none: any quantity may be absent.
    enum class Require
    {
        thermo,
        none,
    };

    std::size_t height = npos;
    // The height is altitude above sea level, to have the ground elevation
    // taken off.
    bool above_sea_level = false;
    std::size_t temperature = npos;
    // False when the header did not say kelvin or Celsius and the values must.
    bool temperature_unit_named = true;
    bool celsius = false;
    std::size_t humidity = npos;
    std::size_t pressure = npos;
    bool pressure_unit_named = true;
    double pressure_to_kpa = 1.0;
    // Wind as speed and the bearing it blows from, or as u and v.
    std::size_t speed = npos;
    std::size_t direction = npos;
    double speed_to_mps = 1.0;
    std::size_t u = npos;
    std::size_t v = npos;
    double component_to_mps = 1.0;
    // Per header column: claimed by one of the above.
    std::vector<bool> used;

    // The fields a row needs to reach every column found.
    std::size_t needed() const;
};

// Throws input_error, naming `what` and `path`, for a file with no height
// column, a wind column without its unit, or half a wind.
ProfileColumns find_profile_columns(const LabeledCsv& csv, LengthUnit unit,
        const std::string& path, ProfileColumns::Require require,
        const std::string& what = "--atmosphere");

// The levels synthetic_profile places from `anchors` (meters, above the
// ground): the anchors sorted, those within 1e-9 of their upper neighbor
// merged, and between each consecutive pair ceil(n log10(ratio)) geometric
// steps, n = levels_per_decade.  Returns the levels from the lowest anchor up,
// the anchors exact.  The profile builder shares it, so a change here moves
// both, and build_atmosphere_surface_is_wind says so.
std::vector<double> geometric_levels(std::vector<double> anchors, int levels_per_decade);

// Write a profile in the form read_atmosphere_profile reads: z in `unit`, T_K,
// RH, p_kPa when it has pressure, and u_mps, v_mps when it has wind, every value
// to 17 significant digits so that reading it back gives the same doubles (the
// heights to an ULP of the unit conversion).  --write_atmosphere uses it, so a
// run made with --wind can be repeated, and inspected, through --atmosphere.
void write_atmosphere_profile(const AtmosphereProfile& profile, const std::string& path,
        LengthUnit unit);

// How the --wind profile's speed grows with height.  `neutral` is the
// logarithmic law of the neutral surface layer, the one measured against the
// 2017 data.  The others are the Klug-Manier dispersion classes of TA Luft
// (1986), each a power law u ~ z^m with a temperature gradient: unvalidated
// here, and a run that takes one says so.
enum class Stability
{
    neutral,
    very_stable, // I    (Pasquill F)
    stable, // II        (Pasquill E)
    neutral_stable, // III/1 (Pasquill D)
    neutral_unstable, // III/2 (Pasquill C)
    unstable, // IV      (Pasquill B)
    very_unstable, // V  (Pasquill A)
};

// The names --stability takes: neutral, the Klug-Manier classes I, II, III/1,
// III/2, IV, V, or the Pasquill letters A to F they correspond to.
const std::vector<std::string>& stability_names();
// Throws input_error for a name not in stability_names().
Stability parse_stability(const std::string& name);
// "neutral (log law)", "II (stable; Pasquill E)", ...
std::string describe(Stability stability);

// One surface wind observation, and the law that carries it up through the
// surface layer: the inputs of synthetic_profile().
struct SurfaceWind
{
    // m/s, at `height`.
    double speed = 0.0;
    // Degrees, the compass bearing the wind blows FROM (the meteorological
    // convention, as an --atmosphere file's wind_direction column).
    double direction = 0.0;
    // Meters above the observers: the anemometer.  10 m is the WMO standard
    // exposure for surface wind.
    double height = 10.0;
    // Meters.  0.03 m is open flat terrain -- short grass, an airfield --
    // in the Davenport-Wieringa classification (WMO-No. 8, Annex 5.B).
    double roughness_length = 0.03;
    // Meters: the law holds up to here, and the wind is held at its value
    // there above it.  The neutral surface layer is about a tenth of a
    // boundary layer one to two kilometers deep.
    double surface_layer_depth = 200.0;
    Stability stability = Stability::neutral;
    // Levels per factor of ten in height.  The profile is linear between its
    // levels; at 16 the plate run over the Coyle array moves by at most 0.016 dB
    // when they are refined to 64, at 8 by 0.053.  --profile_levels sets it,
    // which is how a run checks that nothing depends on it.
    int levels_per_decade = 16;

    // The law: m/s at height z, meters.  Zero at and below the ground (and,
    // for the log law, at and below the roughness length, where no level is
    // ever placed), held above surface_layer_depth.
    double speed_at(double z) const;
    // dT/dz in the surface layer, K/m: dry adiabatic (-9.8 K/km) for neutral,
    // the class's own otherwise.  Above the surface layer every class lapses
    // dry-adiabatically.
    double temperature_gradient() const;
    // The lowest level above the ground: ten roughness lengths, or the
    // anemometer if that is lower.  Between it and z = 0 the wind is linear.
    double first_level() const;
};

// The dry-adiabatic lapse rate, K/m: g / c_p rounded as meteorology quotes it.
inline constexpr double dry_adiabatic_lapse = -0.0098;

// Why a SurfaceWind cannot be built into a profile, or empty when it can.
// Lengths in the message are in `unit`.
std::string surface_wind_problem(const SurfaceWind& wind, LengthUnit unit);

// Parse --wind's SPEED,DIRECTION[,HEIGHT]: the speed with its unit named
// (12kt, 6.2mps, 14mph; nothing about the number tells them apart, which is
// why an --atmosphere file's wind columns name theirs too), the direction in
// degrees it blows from, and the anemometer height in `unit` (10 m when
// omitted).  Throws input_error.
SurfaceWind parse_surface_wind(const std::string& text, LengthUnit unit);

// A profile from one surface observation (--wind), in place of a file.
//
// Levels: z = 0 with no wind (the law's limit at the ground), then geometric
// from first_level() to the top at levels_per_decade, with the anemometer
// height, the top of the surface layer and the top itself as exact levels.
// The top is max(surface_layer_depth, 1.05 * highest), `highest` being the
// highest point any path reaches.  Above that the effective sound speed falls
// with height in every direction (the wind is held and the air lapses), so no
// ray turns back from air the profile leaves out.
//
// No level above the ground is ever calm: --refraction stratified takes a
// calm level just above z = 0 as a thin layer that put levels 22-28 dB low.
// A calm observation (speed 0) gives a profile without wind instead.
//
// Temperature: surface.temperature at z = 0, then temperature_gradient()
// through the surface layer and dry adiabatic above.  Pressure: hydrostatic
// from surface.pressure (with_hydrostatic_pressure).  Humidity: the molar
// concentration of water vapor -- what ISO 9613-1's absorption takes -- is
// held at its surface value, as in a well-mixed layer, so the relative
// humidity rises as the air cools; capped at 100 percent.
AtmosphereProfile synthetic_profile(const SurfaceWind& wind, const Atmosphere& surface,
        double highest);

// Per-band absorption integrated through a profile, tabulated in height so the
// observer loop never evaluates ISO 9613-1.
//
// The table holds alpha_b(z) at nodes a meter apart and its running integral
// F_b(z).  Between nodes alpha is linear and F its exact integral, a
// quadratic.  That is why (F_b(z_a) - F_b(z_b)) / (z_a - z_b) tends smoothly to
// alpha_b as the two heights meet, rather than jumping at a node.
//
// The observer loop needs F at a state's height and at the observer's.  Both
// are fixed while the other varies, so they are looked up once per state and
// once per observer.  Per pair, a band then costs a subtraction and a multiply.
// On a grid, where every observer shares a height, the caller forms each
// state's mean once and pays nothing per pair.
//
// Straight paths only.  For a curved ray the same mean is a first-order
// approximation to the integral along the arc.  The mean is taken over the
// whole height span, including the stretch inside the sphere radius that the
// path length excludes.  For the shipped databases that stretch is a few
// meters against hundreds.
class LayeredAbsorption
{
public:
    LayeredAbsorption() = default;
    // Tabulate over heights [z_low, z_high] m, for bands [first_band, last_band)
    // of `frequencies`.  The profile must carry pressure.
    LayeredAbsorption(const AtmosphereProfile& profile, const std::vector<double>& frequencies,
            std::size_t first_band, std::size_t last_band, double z_low, double z_high);

    // F_b(z), dB: the absorption integral from z_low up to z.  Written to
    // out[band] for each band in range; entries outside it are left alone.
    void integral(double z, double* out) const;

    // alpha_b(z), dB/m, written like integral().
    void coefficients(double z, double* out) const;

    // The mean coefficient over the heights between z_a and z_b, given
    // integral() at each: (F_a - F_b) / (z_a - z_b).  Within a millimeter of
    // each other it takes the coefficient at their midpoint instead, where the
    // difference quotient would lose its digits.  Written like integral().
    void mean_coefficients(double z_a, const double* integral_a, double z_b,
            const double* integral_b, double* out) const;

    bool empty() const
    {
        return alpha_.empty();
    }
    std::size_t first_band() const
    {
        return first_band_;
    }
    std::size_t last_band() const
    {
        return last_band_;
    }
    double z_low() const
    {
        return z_low_;
    }
    double z_high() const
    {
        return z_low_ + step_ * static_cast<double>(nodes_ - 1);
    }

private:
    // The node below z, and how far toward the next one z lies (0 to 1 inside
    // the table; outside it the end cell extrapolates, which covers round-off
    // at the edges).
    void locate(double z, std::size_t& node, double& fraction) const;

    double z_low_ = 0.0;
    double step_ = 1.0;
    double inv_step_ = 1.0;
    std::size_t nodes_ = 0;
    std::size_t first_band_ = 0;
    std::size_t last_band_ = 0;
    std::size_t bands_ = 0;
    // [node][band - first_band_]
    std::vector<double> alpha_;
    std::vector<double> integral_;
};

} // namespace niceops
