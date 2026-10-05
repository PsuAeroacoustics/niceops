#include "atmosphere_profile.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "util/csv.hpp"
#include "util/logger.hpp"
#include "util/paths.hpp"

namespace niceops {

namespace {

// Standard gravity and the dry-air gas constant.  The latter is the value
// trajectory.cpp takes pressure from density with, so the two agree.
constexpr double standard_gravity = 9.80665; // m/s^2
constexpr double dry_air_gas_constant = 287.05; // J/(kg K)

// Plausible ranges, in the stored units, for telling a unit from the values
// and for catching a column labeled with the wrong one.  They do not overlap:
// no air is 150 K cold or 70 C hot.
constexpr double kelvin_low = 150.0;
constexpr double kelvin_high = 350.0;
constexpr double celsius_low = -100.0;
constexpr double celsius_high = 70.0;

enum class Unit
{
    named, // the header said which
    inferred, // the header did not, and the values did
};

// Linear interpolation of `values` against ascending `heights`, held at the
// end values outside them.
double interpolate(const std::vector<double>& heights, const std::vector<double>& values, double z)
{
    if (heights.size() == 1 || z <= heights.front())
        return values.front();
    if (z >= heights.back())
        return values.back();
    const auto upper = std::upper_bound(heights.begin(), heights.end(), z);
    const std::size_t i = static_cast<std::size_t>(upper - heights.begin()) - 1;
    const double t = (z - heights[i]) / (heights[i + 1] - heights[i]);
    return values[i] + t * (values[i + 1] - values[i]);
}

} // namespace

Atmosphere AtmosphereProfile::at(double z) const
{
    Atmosphere air;
    air.temperature = interpolate(height, temperature, z);
    air.relative_humidity = interpolate(height, relative_humidity, z);
    if (height.size() == 1 || z <= height.front())
    {
        air.pressure = pressure.front();
    }
    else if (z >= height.back())
    {
        air.pressure = pressure.back();
    }
    else
    {
        const auto upper = std::upper_bound(height.begin(), height.end(), z);
        const std::size_t i = static_cast<std::size_t>(upper - height.begin()) - 1;
        const double t = (z - height[i]) / (height[i + 1] - height[i]);
        air.pressure = pressure[i] * std::pow(pressure[i + 1] / pressure[i], t);
    }
    return air;
}

std::array<double, 2> AtmosphereProfile::wind_at(double z) const
{
    if (!has_wind())
        return {0.0, 0.0};
    return {interpolate(height, wind_east, z), interpolate(height, wind_north, z)};
}

AtmosphereProfile AtmosphereProfile::with_hydrostatic_pressure(double surface_pressure) const
{
    // dp/p = -(g / R T) dz, so p(z) = p(0) exp(-(g/R) G(z)) with G the
    // integral of 1/T from 0 to z.  T is linear between levels and held
    // constant outside them, so G is exact segment by segment: over a stretch
    // where T runs linearly from T1 to T2 it is dz ln(T2/T1) / (T2 - T1).
    auto segment = [](double dz, double t1, double t2) {
        if (std::fabs(t2 - t1) < 1e-9 * t1)
            return dz / t1;
        return dz * std::log(t2 / t1) / (t2 - t1);
    };
    // Integral of 1/T from the lowest level up to z (negative below it).
    auto from_bottom = [&](double z) {
        if (z <= height.front())
            return (z - height.front()) / temperature.front();
        double sum = 0.0;
        for (std::size_t i = 0; i + 1 < height.size(); i++)
        {
            if (z <= height[i + 1])
            {
                const double t_z = temperature[i]
                        + (z - height[i]) / (height[i + 1] - height[i])
                                * (temperature[i + 1] - temperature[i]);
                return sum + segment(z - height[i], temperature[i], t_z);
            }
            sum += segment(height[i + 1] - height[i], temperature[i], temperature[i + 1]);
        }
        return sum + (z - height.back()) / temperature.back();
    };

    AtmosphereProfile filled = *this;
    filled.pressure.resize(height.size());
    const double at_ground = from_bottom(0.0);
    const double scale = standard_gravity / dry_air_gas_constant;
    for (std::size_t i = 0; i < height.size(); i++)
        filled.pressure[i] = surface_pressure * std::exp(-scale * (from_bottom(height[i]) - at_ground));
    return filled;
}

AtmosphereProfile read_atmosphere_profile(const std::string& path, LengthUnit unit,
        double ground_elevation)
{
    LabeledCsv csv;
    CsvReport report;
    try
    {
        csv = read_labeled_csv(path, &report);
    }
    catch (const std::exception& e)
    {
        input_error("--atmosphere: ", e.what());
    }
    if (report.dropped > 0)
    {
        // A level with one unparseable cell used to vanish from the profile,
        // quietly changing every gradient fitted through it.
        input_error("--atmosphere ", path, ": line ", report.first_dropped_line,
                " is not a row of numbers", report.dropped > 1
                        ? " (and " + std::to_string(report.dropped - 1) + " more are not)"
                        : std::string(),
                ". A text line is only allowed above the first row of data.");
    }

    const bool metric = unit.to_meters == 1.0;
    std::vector<bool> used(csv.header.size(), false);
    auto take = [&used](std::size_t column) {
        if (column != LabeledCsv::npos)
            used[column] = true;
        return column;
    };

    // Height.  As with a geodetic track's altitude, only the unit spellings
    // that agree with --units are matched: a column that says "z_m" while the
    // run speaks feet fails rather than being read three times too large.
    bool above_sea_level = false;
    std::size_t height_column = csv.column({"z", "agl", "z agl", "height agl"});
    if (height_column == LabeledCsv::npos)
    {
        height_column = metric ? csv.column({"z m", "agl m", "z agl m", "height agl m"})
                               : csv.column({"z ft", "agl ft", "z agl ft", "height agl ft"});
    }
    if (height_column == LabeledCsv::npos)
    {
        above_sea_level = true;
        height_column = csv.column({"alt", "altitude", "alt msl", "altitude msl", "height",
                                          "height msl"});
        if (height_column == LabeledCsv::npos)
        {
            height_column = metric
                    ? csv.column({"alt m", "altitude m", "alt msl m", "height m", "height msl m"})
                    : csv.column({"alt ft", "altitude ft", "alt msl ft", "height ft",
                                        "height msl ft"});
        }
    }
    if (height_column == LabeledCsv::npos)
    {
        input_error("--atmosphere ", path, " names no height column. Give height above the "
                "observers as \"z\" (or \"agl\"), or altitude above sea level as \"alt\", in ",
                unit.name, "; a column whose name carries a unit (z_ft, z_m) has to agree "
                "with --units.");
    }
    take(height_column);

    Unit temperature_unit = Unit::named;
    bool celsius = false;
    std::size_t temperature_column = csv.column({"t k", "temp k", "temperature k"});
    if (temperature_column == LabeledCsv::npos)
    {
        temperature_column = csv.column({"t c", "temp c", "temperature c", "t degc",
                                               "temperature degc"});
        celsius = temperature_column != LabeledCsv::npos;
    }
    if (temperature_column == LabeledCsv::npos)
    {
        temperature_column = csv.column({"t", "temp", "temperature"});
        temperature_unit = Unit::inferred;
    }
    if (temperature_column == LabeledCsv::npos)
    {
        input_error("--atmosphere ", path, " names no temperature column. Give \"T_K\" or "
                "\"T_C\".");
    }
    take(temperature_column);

    const std::size_t humidity_column
            = take(csv.column({"rh", "relative humidity", "humidity", "rh pct", "rh percent"}));
    if (humidity_column == LabeledCsv::npos)
    {
        input_error("--atmosphere ", path, " names no relative humidity column. Give \"RH\", "
                "in percent.");
    }

    // Pressure is optional: without it, it is integrated hydrostatically from
    // the surface pressure (with_hydrostatic_pressure).
    Unit pressure_unit = Unit::named;
    double pressure_to_kpa = 1.0;
    std::size_t pressure_column = csv.column({"p kpa", "pressure kpa", "pres kpa"});
    if (pressure_column == LabeledCsv::npos)
    {
        pressure_column = csv.column({"p hpa", "pressure hpa", "pres hpa", "p mb", "pressure mb",
                                            "p mbar", "pressure mbar"});
        if (pressure_column != LabeledCsv::npos)
            pressure_to_kpa = 0.1;
    }
    if (pressure_column == LabeledCsv::npos)
    {
        pressure_column = csv.column({"p pa", "pressure pa", "pres pa"});
        if (pressure_column != LabeledCsv::npos)
            pressure_to_kpa = 1e-3;
    }
    if (pressure_column == LabeledCsv::npos)
    {
        pressure_column = csv.column({"p", "pressure", "pres"});
        pressure_unit = Unit::inferred;
    }
    take(pressure_column);
    const bool have_pressure = pressure_column != LabeledCsv::npos;

    // Wind, as speed and the direction it blows from, or as east and north
    // components.  Its unit has to be named: m/s and knots cover the same
    // values, so unlike temperature nothing in the numbers can tell them apart.
    constexpr double knot = 1852.0 / 3600.0; // m/s
    constexpr std::size_t npos = LabeledCsv::npos;
    double speed_to_mps = 1.0;
    std::size_t speed_column = csv.column({"wind speed mps", "wind speed ms", "wspd mps",
            "wspd ms", "ws mps", "ws ms"});
    if (speed_column == npos)
    {
        speed_column = csv.column({"wind speed kt", "wind speed kts", "wind speed knots",
                "wspd kt", "wspd kts", "ws kt", "ws kts"});
        speed_to_mps = knot;
    }
    const std::size_t direction_column = csv.column({"wind direction", "wind dir", "wdir", "wd",
            "wind direction deg", "wdir deg", "wd deg"});
    double component_to_mps = 1.0;
    std::size_t u_column = csv.column({"u mps", "u ms"});
    std::size_t v_column = csv.column({"v mps", "v ms"});
    if (u_column == npos && v_column == npos)
    {
        u_column = csv.column({"u kt", "u kts"});
        v_column = csv.column({"v kt", "v kts"});
        component_to_mps = knot;
    }
    if (speed_column == npos && u_column == npos && v_column == npos
            && csv.column({"wind speed", "wspd", "ws", "wind", "u", "v"}) != npos)
    {
        input_error("--atmosphere ", path, ": name the wind's unit in its column, as "
                "wind_speed_mps or wind_speed_kt (with wind_direction), or u_mps and v_mps.");
    }
    if ((speed_column == npos) != (direction_column == npos))
    {
        input_error("--atmosphere ", path, ": wind speed and wind direction come together; "
                "this file gives one without the other.");
    }
    if ((u_column == npos) != (v_column == npos))
    {
        input_error("--atmosphere ", path, ": the wind components u and v come together; this "
                "file gives one without the other.");
    }
    if (speed_column != npos && u_column != npos)
    {
        input_error("--atmosphere ", path, " gives the wind both as speed and direction and as "
                "components; give one.");
    }
    const bool polar_wind = speed_column != npos;
    const bool component_wind = u_column != npos;
    take(speed_column);
    take(direction_column);
    take(u_column);
    take(v_column);

    std::size_t needed = 1
            + std::max({height_column, temperature_column, humidity_column,
                    have_pressure ? pressure_column : 0});
    for (const std::size_t column : {speed_column, direction_column, u_column, v_column})
    {
        if (column != npos)
            needed = std::max(needed, column + 1);
    }
    if (csv.rows.empty())
        input_error("--atmosphere ", path, " holds no numeric rows under its header.");
    // nan and inf parse as numbers; the range checks below compare against them
    // and a NaN passes every one of those, then runs through the ray fits.
    for (std::size_t row = 0; row < csv.rows.size(); row++)
    {
        for (const double value : csv.rows[row])
        {
            if (!std::isfinite(value))
            {
                input_error("--atmosphere ", path, ": data row ", row + 1,
                        " holds a value that is not a number.");
            }
        }
    }

    // Gather, then order by height: a sounding may well run top down.
    std::vector<std::size_t> order(csv.rows.size());
    std::iota(order.begin(), order.end(), 0);
    for (std::size_t row = 0; row < csv.rows.size(); row++)
    {
        if (csv.rows[row].size() < needed)
        {
            input_error("Data row ", row + 1, " of ", path, " has ", csv.rows[row].size(),
                    " fields, too few for the columns its header names.");
        }
    }
    auto height_of = [&](std::size_t row) {
        const double raw = csv.rows[row][height_column] * unit.to_meters;
        return above_sea_level ? raw - ground_elevation : raw;
    };
    std::sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b) { return height_of(a) < height_of(b); });

    AtmosphereProfile profile;
    profile.source = path;
    for (const std::size_t row : order)
    {
        const auto& data = csv.rows[row];
        profile.height.push_back(height_of(row));
        profile.temperature.push_back(data[temperature_column]);
        profile.relative_humidity.push_back(data[humidity_column]);
        if (have_pressure)
            profile.pressure.push_back(data[pressure_column]);
        if (polar_wind)
        {
            // The meteorological convention: the bearing the wind comes from.
            const double speed = data[speed_column] * speed_to_mps;
            const double from = data[direction_column] * pi / 180.0;
            if (!(speed >= 0.0) || speed > 150.0)
            {
                input_error("--atmosphere ", path, ": wind speed ", speed,
                        " m/s is not plausible.");
            }
            profile.wind_east.push_back(-speed * std::sin(from));
            profile.wind_north.push_back(-speed * std::cos(from));
        }
        else if (component_wind)
        {
            const double east = data[u_column] * component_to_mps;
            const double north = data[v_column] * component_to_mps;
            if (!(std::hypot(east, north) <= 150.0))
            {
                input_error("--atmosphere ", path, ": wind of ", std::hypot(east, north),
                        " m/s is not plausible.");
            }
            profile.wind_east.push_back(east);
            profile.wind_north.push_back(north);
        }
    }
    for (std::size_t i = 1; i < profile.height.size(); i++)
    {
        if (!(profile.height[i] > profile.height[i - 1]))
        {
            input_error("--atmosphere ", path, " gives two levels at ",
                    profile.height[i] / unit.to_meters, " ", unit.name, ".");
        }
    }

    const auto [t_min, t_max]
            = std::minmax_element(profile.temperature.begin(), profile.temperature.end());
    if (temperature_unit == Unit::inferred)
    {
        if (*t_min >= kelvin_low && *t_max <= kelvin_high)
            celsius = false;
        else if (*t_min >= celsius_low && *t_max <= celsius_high)
            celsius = true;
        else
        {
            input_error("--atmosphere ", path, ": temperatures from ", *t_min, " to ", *t_max,
                    " are neither kelvin nor Celsius. Name the column T_K or T_C.");
        }
    }
    const double t_low = celsius ? celsius_low : kelvin_low;
    const double t_high = celsius ? celsius_high : kelvin_high;
    if (*t_min < t_low || *t_max > t_high)
    {
        input_error("--atmosphere ", path, ": temperatures from ", *t_min, " to ", *t_max,
                " are not plausible in ", celsius ? "Celsius" : "kelvin", ".");
    }
    if (celsius)
    {
        for (double& t : profile.temperature)
            t += 273.15;
    }

    const auto [rh_min, rh_max] = std::minmax_element(profile.relative_humidity.begin(),
            profile.relative_humidity.end());
    if (*rh_min < 0.0 || *rh_max > 100.0)
    {
        input_error("--atmosphere ", path, ": relative humidity from ", *rh_min, " to ", *rh_max,
                " is outside 0 to 100 percent.");
    }
    if (*rh_max <= 1.0)
    {
        // Nearly dry air all the way up is not impossible, but a fraction read
        // as a percentage is far more likely, and it would propagate through
        // absorption with no other sign.
        input_error("--atmosphere ", path, ": relative humidity never exceeds ", *rh_max,
                " percent. It is read in percent; a fraction (0 to 1) has to be scaled by 100.");
    }

    if (have_pressure)
    {
        const double p_max = *std::max_element(profile.pressure.begin(), profile.pressure.end());
        if (pressure_unit == Unit::inferred)
        {
            if (p_max <= 120.0)
                pressure_to_kpa = 1.0;
            else if (p_max <= 1200.0)
                pressure_to_kpa = 0.1;
            else if (p_max <= 120000.0)
                pressure_to_kpa = 1e-3;
            else
                input_error("--atmosphere ", path, ": pressure up to ", p_max,
                        " is not kPa, hPa or Pa. Name the column p_kPa, p_hPa or p_Pa.");
        }
        for (double& p : profile.pressure)
            p *= pressure_to_kpa;
        const auto [p_lo, p_hi] = std::minmax_element(profile.pressure.begin(), profile.pressure.end());
        if (!(*p_lo > 0.0) || *p_hi > 120.0)
        {
            input_error("--atmosphere ", path, ": pressure from ", *p_lo, " to ", *p_hi,
                    " kPa is not plausible.");
        }
        for (std::size_t i = 1; i < profile.pressure.size(); i++)
        {
            if (profile.pressure[i] > profile.pressure[i - 1])
            {
                input_error("--atmosphere ", path, ": pressure rises with height at ",
                        profile.height[i] / unit.to_meters, " ", unit.name, ".");
            }
        }
    }

    std::ostringstream ignored;
    for (std::size_t i = 0; i < csv.header.size(); i++)
    {
        if (used[i])
            continue;
        const auto& name = csv.header[i];
        const auto first = name.find_first_not_of(" \t");
        const auto last = name.find_last_not_of(" \t");
        ignored << (ignored.tellp() > 0 ? ", " : "")
                << (first == std::string::npos ? name : name.substr(first, last - first + 1));
    }
    std::ostringstream levels;
    if (profile.height.size() == 1)
        levels << "1 level, uniform at every height";
    else
        levels << profile.height.size() << " levels from " << profile.bottom() / unit.to_meters
               << " to " << profile.top() / unit.to_meters << " " << unit.name
               << " above the observers";
    log("Atmosphere profile ", path, ": ", levels.str(),
            above_sea_level ? " (read as altitude, less --ground_elevation)" : "",
            "; temperature in ", celsius ? "Celsius" : "kelvin",
            temperature_unit == Unit::inferred ? " (told from the values)" : "",
            have_pressure ? (pressure_unit == Unit::inferred ? "; pressure unit told from the values"
                                                             : "")
                          : "; no pressure column, so it is integrated hydrostatically from the "
                            "surface pressure",
            polar_wind ? (speed_to_mps == 1.0 ? "; wind as speed (m/s) and direction"
                                              : "; wind as speed (knots) and direction")
                    : component_wind ? (component_to_mps == 1.0 ? "; wind as u, v (m/s)"
                                                                : "; wind as u, v (knots)")
                                     : "; no wind",
            ".");
    if (ignored.tellp() > 0)
    {
        log("  ignoring columns: ", ignored.str(), ".");
    }
    return profile;
}

void write_atmosphere_profile(const AtmosphereProfile& profile, const std::string& path,
        LengthUnit unit)
{
    ensure_parent_directory(path);
    std::ofstream out(path);
    if (!out)
        input_error("--write_atmosphere: cannot open ", path, " for writing.");
    const bool metric = unit.to_meters == 1.0;
    out << (metric ? "z_m" : "z_ft") << ",T_K,RH";
    if (profile.has_pressure())
        out << ",p_kPa";
    if (profile.has_wind())
        out << ",u_mps,v_mps";
    out << '\n' << std::setprecision(17);
    for (std::size_t i = 0; i < profile.height.size(); i++)
    {
        out << profile.height[i] / unit.to_meters << ',' << profile.temperature[i] << ','
            << profile.relative_humidity[i];
        if (profile.has_pressure())
            out << ',' << profile.pressure[i];
        if (profile.has_wind())
            out << ',' << profile.wind_east[i] << ',' << profile.wind_north[i];
        out << '\n';
    }
    out.close();
    if (!out)
        input_error("--write_atmosphere: writing ", path, " failed.");
}

namespace {

struct StabilityClass
{
    Stability stability;
    const char* klug_manier; // TA Luft's name
    const char* pasquill;
    // The TA Luft (1986) wind profile exponent: u ~ z^m.
    double exponent;
    // dT/dz, K/m: the middle of the class's range in US NRC Regulatory Guide
    // 1.23's delta-T classification (Pasquill A, below -1.9 K/100 m, at -2.0).
    double temperature_gradient;
    const char* character;
};

constexpr StabilityClass stability_classes[] = {
        {Stability::very_stable, "I", "F", 0.42, 0.0275, "very stable"},
        {Stability::stable, "II", "E", 0.35, 0.005, "stable"},
        {Stability::neutral_stable, "III/1", "D", 0.22, -0.010, "neutral to stable"},
        {Stability::neutral_unstable, "III/2", "C", 0.20, -0.016, "neutral to unstable"},
        {Stability::unstable, "IV", "B", 0.15, -0.018, "unstable"},
        {Stability::very_unstable, "V", "A", 0.09, -0.020, "very unstable"},
};

const StabilityClass& class_of(Stability stability)
{
    for (const auto& entry : stability_classes)
    {
        if (entry.stability == stability)
            return entry;
    }
    throw std::logic_error("no stability class for neutral: it is the log law");
}

std::string lower(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
}

std::string trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t");
    if (first == std::string::npos)
        return "";
    const auto last = text.find_last_not_of(" \t");
    return text.substr(first, last - first + 1);
}

// A number at the front of `text`, and whatever follows it, trimmed.  False
// when there is no number there.
bool leading_number(const std::string& text, double& value, std::string& rest)
{
    const char* begin = text.c_str();
    char* end = nullptr;
    value = std::strtod(begin, &end);
    if (end == begin)
        return false;
    rest = trim(std::string(end));
    return true;
}

std::string length_text(double meters, LengthUnit unit)
{
    std::ostringstream text;
    text << meters / unit.to_meters << " " << unit.name;
    return text.str();
}

} // namespace

const std::vector<std::string>& stability_names()
{
    static const std::vector<std::string> names = [] {
        std::vector<std::string> list {"neutral"};
        for (const auto& entry : stability_classes)
            list.emplace_back(entry.klug_manier);
        for (const auto& entry : stability_classes)
            list.emplace_back(entry.pasquill);
        return list;
    }();
    return names;
}

Stability parse_stability(const std::string& name)
{
    const std::string wanted = lower(trim(name));
    if (wanted == "neutral")
        return Stability::neutral;
    for (const auto& entry : stability_classes)
    {
        if (wanted == lower(entry.klug_manier) || wanted == lower(entry.pasquill))
            return entry.stability;
    }
    input_error("--stability takes neutral, a TA Luft class (I, II, III/1, III/2, IV, V) or a "
                "Pasquill class (A to F), not '", name, "'.");
}

std::string describe(Stability stability)
{
    if (stability == Stability::neutral)
        return "neutral (log law)";
    const auto& entry = class_of(stability);
    std::ostringstream text;
    text << "TA Luft class " << entry.klug_manier << " (" << entry.character << "; Pasquill "
         << entry.pasquill << "): power law, exponent " << entry.exponent << ", "
         << entry.temperature_gradient * 1000.0 << " K/km in the surface layer";
    return text.str();
}

double SurfaceWind::speed_at(double z) const
{
    if (!(z > 0.0) || speed == 0.0)
        return 0.0;
    const double held = std::min(z, surface_layer_depth);
    if (stability == Stability::neutral)
    {
        if (held <= roughness_length)
            return 0.0;
        return speed * std::log(held / roughness_length) / std::log(height / roughness_length);
    }
    return speed * std::pow(held / height, class_of(stability).exponent);
}

double SurfaceWind::temperature_gradient() const
{
    return stability == Stability::neutral ? dry_adiabatic_lapse
                                           : class_of(stability).temperature_gradient;
}

double SurfaceWind::first_level() const
{
    return std::min(10.0 * roughness_length, height);
}

std::string surface_wind_problem(const SurfaceWind& wind, LengthUnit unit)
{
    std::ostringstream problem;
    if (!std::isfinite(wind.speed) || wind.speed < 0.0 || wind.speed > 50.0)
    {
        problem << "a surface wind of " << wind.speed << " m/s is not plausible (0 to 50 m/s, "
                << "97 kt, are accepted).";
    }
    else if (!std::isfinite(wind.direction) || wind.direction < 0.0 || wind.direction > 360.0)
    {
        problem << "a wind direction of " << wind.direction << " degrees is outside 0 to 360 "
                << "(the compass bearing the wind blows from).";
    }
    else if (!std::isfinite(wind.roughness_length) || !(wind.roughness_length > 0.0)
            || wind.roughness_length > 3.0)
    {
        problem << "a roughness length of " << length_text(wind.roughness_length, unit)
                << " is not plausible: it runs from 0.0002 m over calm water to about 2 m over "
                << "city centers (more than 0 and up to 3 m are accepted).";
    }
    else if (!std::isfinite(wind.height) || !(wind.height > 0.0))
    {
        problem << "the anemometer height has to be above the ground, not "
                << length_text(wind.height, unit) << ".";
    }
    else if (wind.stability == Stability::neutral && !(wind.height > wind.roughness_length))
    {
        problem << "the log law is zero at the roughness length, so the anemometer ("
                << length_text(wind.height, unit) << ") has to stand above it ("
                << length_text(wind.roughness_length, unit) << ").";
    }
    else if (!std::isfinite(wind.surface_layer_depth)
            || !(wind.surface_layer_depth > wind.height))
    {
        problem << "the surface layer (" << length_text(wind.surface_layer_depth, unit)
                << " deep) has to reach above the anemometer ("
                << length_text(wind.height, unit) << ").";
    }
    else if (wind.levels_per_decade < 1 || wind.levels_per_decade > 100)
    {
        problem << "levels per factor of ten in height run from 1 to 100, not "
                << wind.levels_per_decade << ".";
    }
    else if (wind.speed_at(wind.surface_layer_depth) > 150.0)
    {
        problem << "the law carries " << wind.speed << " m/s at "
                << length_text(wind.height, unit) << " to "
                << wind.speed_at(wind.surface_layer_depth) << " m/s at "
                << length_text(wind.surface_layer_depth, unit)
                << ", which is not plausible. Is the anemometer that close to the roughness "
                << "length?";
    }
    return problem.str();
}

SurfaceWind parse_surface_wind(const std::string& text, LengthUnit unit)
{
    std::vector<std::string> fields;
    std::stringstream stream(text);
    std::string field;
    while (std::getline(stream, field, ','))
        fields.push_back(trim(field));
    if (!text.empty() && text.back() == ',')
        fields.emplace_back();
    const std::string usage = "--wind takes SPEED,DIRECTION[,HEIGHT]: the speed with its unit "
                              "(12kt, 6.2mps or 14mph), the compass bearing in degrees the wind "
                              "blows from, and the anemometer's height in ";
    if (fields.size() < 2 || fields.size() > 3)
        input_error(usage, unit.name, " (default 10 m). Got '", text, "'.");

    SurfaceWind wind;
    std::string rest;
    if (!leading_number(fields[0], wind.speed, rest))
        input_error(usage, unit.name, " (default 10 m). '", fields[0], "' is not a speed.");
    const std::string speed_unit = lower(rest);
    if (speed_unit == "mps" || speed_unit == "m/s" || speed_unit == "ms")
    {
    }
    else if (speed_unit == "kt" || speed_unit == "kts" || speed_unit == "knot"
            || speed_unit == "knots")
        wind.speed *= 1852.0 / 3600.0;
    else if (speed_unit == "mph")
        wind.speed *= 1609.344 / 3600.0;
    else if (speed_unit.empty())
    {
        input_error("--wind: name the speed's unit, as ", fields[0], "kt, ", fields[0], "mps or ",
                fields[0], "mph. Nothing about the number tells them apart, and it is not "
                "taken from --units.");
    }
    else
    {
        input_error("--wind: '", rest, "' is not a speed unit; give kt, mps or mph.");
    }

    if (!leading_number(fields[1], wind.direction, rest) || !rest.empty())
    {
        input_error("--wind: the direction is the compass bearing, in degrees, that the wind "
                    "blows from; '", fields[1], "' is not one.");
    }
    if (fields.size() == 3)
    {
        double height = 0.0;
        if (!leading_number(fields[2], height, rest) || !rest.empty())
        {
            input_error("--wind: the anemometer height is a length in ", unit.name, " (--units), "
                    "without a unit of its own; '", fields[2], "' is not one.");
        }
        wind.height = height * unit.to_meters;
    }
    return wind;
}

AtmosphereProfile synthetic_profile(const SurfaceWind& wind, const Atmosphere& surface,
        double highest)
{
    const double depth = wind.surface_layer_depth;
    const double top = std::max(depth, 1.05 * highest);
    // Exact levels: the profile passes through the observation itself, and
    // the kink where the surface layer ends is a level rather than smeared
    // over whichever layer straddles it.
    std::vector<double> anchors {wind.first_level(), wind.height, depth, top};
    std::sort(anchors.begin(), anchors.end());
    anchors.erase(std::unique(anchors.begin(), anchors.end(),
                          [](double a, double b) { return b - a <= 1e-9 * b; }),
            anchors.end());

    AtmosphereProfile profile;
    profile.source = "the --wind profile";
    profile.height.push_back(0.0);
    profile.height.push_back(anchors.front());
    for (std::size_t i = 0; i + 1 < anchors.size(); i++)
    {
        const double ratio = anchors[i + 1] / anchors[i];
        const int steps = std::max(1,
                static_cast<int>(std::ceil(wind.levels_per_decade * std::log10(ratio) - 1e-9)));
        for (int j = 1; j < steps; j++)
            profile.height.push_back(anchors[i] * std::pow(ratio, static_cast<double>(j) / steps));
        profile.height.push_back(anchors[i + 1]);
    }

    const double gradient = wind.temperature_gradient();
    const double east = -std::sin(wind.direction * pi / 180.0);
    const double north = -std::cos(wind.direction * pi / 180.0);
    const bool calm = wind.speed == 0.0;
    for (const double z : profile.height)
    {
        profile.temperature.push_back(surface.temperature + gradient * std::min(z, depth)
                + dry_adiabatic_lapse * std::max(z - depth, 0.0));
        profile.relative_humidity.push_back(surface.relative_humidity);
        if (calm)
            continue;
        const double speed = wind.speed_at(z);
        // The thin-layer guard: a calm level above the ground is what made
        // --refraction stratified read 22-28 dB low.  It cannot happen with
        // first_level() above the roughness length, and this says so.
        if (z > 0.0 && !(speed > 0.0))
            throw std::logic_error("synthetic_profile placed a calm level above the ground");
        // + 0.0 turns the ground level's -0 into 0 in a written profile; it
        // leaves every other value as it is.
        profile.wind_east.push_back(speed * east + 0.0);
        profile.wind_north.push_back(speed * north + 0.0);
    }

    profile = profile.with_hydrostatic_pressure(surface.pressure);
    // A well-mixed layer keeps its water vapor: the molar concentration, which
    // is what the absorption takes, is held, and the relative humidity follows
    // the cooling air and falling pressure.  The ground level keeps the surface
    // value exactly rather than its round trip.
    const double concentration = surface.molar_concentration_water_vapor();
    for (std::size_t i = 1; i < profile.height.size(); i++)
    {
        Atmosphere air;
        air.temperature = profile.temperature[i];
        air.pressure = profile.pressure[i];
        profile.relative_humidity[i]
                = std::min(100.0, concentration * air.pressure / air.saturation_pressure());
    }
    return profile;
}

LayeredAbsorption::LayeredAbsorption(const AtmosphereProfile& profile,
        const std::vector<double>& frequencies, std::size_t first_band, std::size_t last_band,
        double z_low, double z_high)
    : z_low_(z_low), first_band_(first_band), last_band_(std::min(last_band, frequencies.size()))
{
    if (!profile.has_pressure())
        throw std::logic_error("LayeredAbsorption needs a profile with pressure");
    if (!(first_band_ < last_band_))
        throw std::logic_error("LayeredAbsorption needs at least one band");
    if (!(z_high >= z_low))
        throw std::logic_error("LayeredAbsorption needs z_high >= z_low");
    bands_ = last_band_ - first_band_;

    // A meter apart: alpha changes by parts in a thousand over that, so the
    // linear interpolant between nodes is exact to far below anything the
    // level can show.  Past 20 km of span the step widens rather than the
    // table growing without limit, which no helicopter trajectory reaches.
    constexpr std::size_t max_nodes = 20001;
    const double span = z_high - z_low;
    nodes_ = std::max<std::size_t>(2, static_cast<std::size_t>(std::ceil(span)) + 1);
    step_ = 1.0;
    if (nodes_ > max_nodes)
    {
        nodes_ = max_nodes;
        step_ = span / static_cast<double>(max_nodes - 1);
    }
    inv_step_ = 1.0 / step_;

    alpha_.resize(nodes_ * bands_);
    integral_.resize(nodes_ * bands_);
    for (std::size_t node = 0; node < nodes_; node++)
    {
        const double z = z_low_ + step_ * static_cast<double>(node);
        const auto alpha = absorption_coefficients(profile.at(z), frequencies);
        std::copy(alpha.begin() + static_cast<std::ptrdiff_t>(first_band_),
                alpha.begin() + static_cast<std::ptrdiff_t>(last_band_),
                alpha_.begin() + static_cast<std::ptrdiff_t>(node * bands_));
    }
    // Trapezoids: the exact integral of the piecewise-linear alpha.
    std::fill(integral_.begin(), integral_.begin() + static_cast<std::ptrdiff_t>(bands_), 0.0);
    for (std::size_t node = 1; node < nodes_; node++)
    {
        for (std::size_t k = 0; k < bands_; k++)
        {
            integral_[node * bands_ + k] = integral_[(node - 1) * bands_ + k]
                    + 0.5 * step_ * (alpha_[(node - 1) * bands_ + k] + alpha_[node * bands_ + k]);
        }
    }
}

void LayeredAbsorption::locate(double z, std::size_t& node, double& fraction) const
{
    const double u = (z - z_low_) * inv_step_;
    node = u <= 0.0 ? 0 : std::min(static_cast<std::size_t>(u), nodes_ - 2);
    fraction = u - static_cast<double>(node);
}

void LayeredAbsorption::integral(double z, double* out) const
{
    std::size_t node;
    double t;
    locate(z, node, t);
    const double* a0 = alpha_.data() + node * bands_;
    const double* a1 = a0 + bands_;
    const double* f0 = integral_.data() + node * bands_;
    for (std::size_t k = 0; k < bands_; k++)
        out[first_band_ + k] = f0[k] + step_ * t * (a0[k] + 0.5 * t * (a1[k] - a0[k]));
}

void LayeredAbsorption::coefficients(double z, double* out) const
{
    std::size_t node;
    double t;
    locate(z, node, t);
    const double* a0 = alpha_.data() + node * bands_;
    const double* a1 = a0 + bands_;
    for (std::size_t k = 0; k < bands_; k++)
        out[first_band_ + k] = a0[k] + t * (a1[k] - a0[k]);
}

void LayeredAbsorption::mean_coefficients(double z_a, const double* integral_a, double z_b,
        const double* integral_b, double* out) const
{
    const double dz = z_a - z_b;
    if (std::fabs(dz) < 1e-3)
    {
        coefficients(0.5 * (z_a + z_b), out);
        return;
    }
    const double inverse = 1.0 / dz;
    for (std::size_t band = first_band_; band < last_band_; band++)
        out[band] = (integral_a[band] - integral_b[band]) * inverse;
}

} // namespace niceops
