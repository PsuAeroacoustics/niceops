#pragma once

// An --atmosphere profile built from measured columns: a met tower, a sodar or
// LIDAR, a temperature string or tethered balloon, a radiosonde, or a weather
// model's column, with or without a surface station.  niceops_build_atmosphere
// (tools/build_atmosphere.cpp) is its front end.
//
// The run reads whatever profile it is given and only documents two things a
// profile should have: a level at the ground, and no calm level above it
// (internals/atmosphere.md).  The builder guarantees both.  Below the lowest
// wind measurement the wind comes from the law --wind uses, through a
// SurfaceWind built on that measurement, so the roughness length is the only
// thing the fill adds, and it has to be given.
//
// Everything here is SI and meters above the observer plane (z = 0), as in
// AtmosphereProfile.

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "atmosphere_profile.hpp"

namespace niceops {

// One variable's measurements, in height order, with the file they came from.
struct ColumnSeries
{
    std::string source;
    std::vector<double> z; // m above the observers; may be below 0
    std::vector<double> value;
};

// Wind measurements: m/s toward the east and the north.
struct WindSeries
{
    std::string source;
    std::vector<double> z;
    std::vector<double> east;
    std::vector<double> north;
};

// What one column file carries.  Any variable may be absent.
struct ColumnSource
{
    std::string path;
    std::optional<ColumnSeries> temperature; // K
    std::optional<ColumnSeries> humidity; // percent, as read (above 100 is left to the build)
    std::optional<ColumnSeries> pressure; // kPa
    std::optional<WindSeries> wind;
    // Cells equal to a --missing value, dropped per variable and level.
    std::size_t sentinels_dropped = 0;
    // Cells outside the plausible range, dropped under drop_implausible.
    std::size_t implausible_dropped = 0;
};

// Read a column file.  Its columns are named as an --atmosphere file's are
// (find_profile_columns), and besides relative humidity it may give the
// dew point (Td_C, Td_K) or the specific humidity (q_kgkg, with pressure in
// the same row), which are turned into relative humidity through ISO 9613-1's
// saturation pressure.  A cell equal to one of `missing` is dropped for its
// variable at its level.  Any other value outside the plausible ranges is
// refused, with the file and row named, unless `drop_implausible`, which drops
// it and counts it.  Throws input_error.
ColumnSource read_column_source(const std::string& path, LengthUnit unit,
        double ground_elevation, const std::vector<double>& missing, bool drop_implausible);

// Every variable's measurements, one source each (the caller settles which
// file gives what).
struct MeasuredColumns
{
    std::optional<ColumnSeries> temperature;
    std::optional<ColumnSeries> humidity;
    std::optional<ColumnSeries> pressure;
    std::optional<WindSeries> wind;
};

// A surface station: the air at z = 0, and a wind at its anemometer's height.
struct SurfaceAir
{
    std::optional<double> temperature; // K
    std::optional<double> humidity; // percent
    std::optional<double> pressure; // kPa
    // Speed, direction and height; the build's rules set the rest.
    std::optional<SurfaceWind> wind;
};

struct BuildRules
{
    // Below the lowest wind measurement: the --wind law through it (log), its
    // value held down to the first level (hold), or linear to calm at the
    // ground (none).
    enum class Fill
    {
        log,
        hold,
        none,
    };
    // Above a variable's highest measurement: held, or (temperature only)
    // lapsing.
    enum class Above
    {
        hold,
        lapse,
    };
    // Humidity above its highest measurement: held, or the water vapor's
    // molar concentration held, as --wind does.
    enum class HumidityAloft
    {
        interpolate,
        held_concentration,
    };

    Fill fill = Fill::log;
    // m.  Required: there is no default.
    double roughness_length = 0.0;
    // m above the observers.  Empty: geometric levels, levels_per_decade a
    // decade, through every measured level.  Otherwise these levels, every
    // measured level, the ground and the top.
    std::vector<double> levels;
    int levels_per_decade = 16;
    // m above the observers.  Required.
    double top = 0.0;
    Above above = Above::hold;
    // With a temperature from the surface alone and Above::lapse: the
    // gradient through the surface layer, and the power law in place of the
    // log law below a single wind level.
    Stability stability = Stability::neutral;
    // With one wind level: the law continues up to here, and is held above,
    // exactly as --wind does.  Refused with several wind levels.
    std::optional<double> surface_layer_depth;
    HumidityAloft humidity_aloft = HumidityAloft::interpolate;
    // Write a pressure column (hydrostatic from the surface pressure) even
    // when nothing else needs one.
    bool write_pressure = false;
    // A measured calm level above the ground: dropped (true) or refused.
    bool drop_calm = false;
    // A relative humidity above 100 percent: clipped to 100 (true) or refused.
    bool clip_humidity = false;
};

struct BuildReport
{
    std::vector<std::string> warnings;
    // Per variable ("temperature", "humidity", "pressure", "wind"), the lowest
    // and highest measurement used, m above the observers.
    std::map<std::string, std::pair<double, double>> measured_span;
    // Stretches of height whose wind no instrument saw, m: below the lowest
    // wind level, and above the highest.
    std::vector<std::pair<double, double>> filled_wind;
    std::size_t below_ground_dropped = 0;
    std::size_t calm_dropped = 0;
    std::size_t humidity_clipped = 0;
};

// Why the rules cannot build a profile, or empty when they can.  Lengths in
// the message are in `unit`.
std::string build_rules_problem(const BuildRules& rules, LengthUnit unit);

// Build the profile.
//
// Ground (z = 0): the surface station's value where it gives one (and that
// variable's measured levels at or below the ground are then dropped);
// otherwise the measurements interpolated to z = 0; otherwise the lowest held
// down, with a warning.  Levels below the ground only set the value there.
// The wind at the ground is calm.
//
// Between measured levels temperature and humidity are linear in height, the
// wind linear in its components, and pressure linear in its logarithm,
// matching AtmosphereProfile.  Pressure outside its measurements, or from the
// surface pressure alone, is hydrostatic (with_hydrostatic_pressure).
//
// No level above the ground is calm: a measured calm level is refused or
// dropped, a level the fill would leave calm is refused, and a column calm at
// every level gives a profile with no wind.  The profile is checked as
// read_atmosphere_profile would check it.  Throws input_error.
AtmosphereProfile build_profile(const MeasuredColumns& measured, const SurfaceAir& surface,
        const BuildRules& rules, LengthUnit unit, BuildReport& report);

} // namespace niceops
