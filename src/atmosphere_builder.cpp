#include "atmosphere_builder.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "util/csv.hpp"
#include "util/logger.hpp"

namespace niceops {

namespace {

// The plausible ranges read_atmosphere_profile applies, in the stored units.
constexpr double kelvin_low = 150.0;
constexpr double kelvin_high = 350.0;
constexpr double celsius_low = -100.0;
constexpr double celsius_high = 70.0;
// A station reading a few percent over saturation is common in humid air; the
// build clips or refuses it, by --clip_humidity.  Past this it is not humidity.
constexpr double humidity_high = 110.0;
// The dry-air gas constant with_hydrostatic_pressure integrates with.
constexpr double dry_air_gas_constant = 287.05; // J/(kg K)

std::string length_text(double meters, LengthUnit unit)
{
    std::ostringstream text;
    text << meters / unit.to_meters << " " << unit.name;
    return text.str();
}

// Linear interpolation, holding the end values outside the knots.
double interpolate(const std::vector<double>& z, const std::vector<double>& value, double at)
{
    if (z.size() == 1 || at <= z.front())
        return value.front();
    if (at >= z.back())
        return value.back();
    const auto upper = std::upper_bound(z.begin(), z.end(), at);
    const std::size_t i = static_cast<std::size_t>(upper - z.begin()) - 1;
    const double t = (at - z[i]) / (z[i + 1] - z[i]);
    return value[i] + t * (value[i + 1] - value[i]);
}

// Linear in the logarithm, as AtmosphereProfile::at takes pressure.
double interpolate_log(const std::vector<double>& z, const std::vector<double>& value, double at)
{
    if (z.size() == 1 || at <= z.front())
        return value.front();
    if (at >= z.back())
        return value.back();
    const auto upper = std::upper_bound(z.begin(), z.end(), at);
    const std::size_t i = static_cast<std::size_t>(upper - z.begin()) - 1;
    const double t = (at - z[i]) / (z[i + 1] - z[i]);
    return value[i] * std::pow(value[i + 1] / value[i], t);
}

double saturation_pressure(double kelvin)
{
    Atmosphere air;
    air.temperature = kelvin;
    return air.saturation_pressure();
}

void note(BuildReport& report, const std::string& message)
{
    warn(message);
    report.warnings.push_back(message);
}

// Rows in height order, with the file's row number for messages.
struct Level
{
    double z;
    double value;
    std::size_t row;
};

template <typename Series>
void refuse_duplicates(const Series& series, const char* what, LengthUnit unit)
{
    for (std::size_t i = 1; i < series.z.size(); i++)
    {
        if (!(series.z[i] > series.z[i - 1]))
        {
            input_error("--column ", series.source, " gives two ", what, " levels at ",
                    length_text(series.z[i], unit), ".");
        }
    }
}

ColumnSeries series_from(const std::string& path, std::vector<Level> levels)
{
    std::stable_sort(levels.begin(), levels.end(),
            [](const Level& a, const Level& b) { return a.z < b.z; });
    ColumnSeries series;
    series.source = path;
    for (const auto& level : levels)
    {
        series.z.push_back(level.z);
        series.value.push_back(level.value);
    }
    return series;
}

} // namespace

ColumnSource read_column_source(const std::string& path, LengthUnit unit,
        double ground_elevation, const std::vector<double>& missing, bool drop_implausible)
{
    constexpr std::size_t npos = LabeledCsv::npos;
    LabeledCsv csv;
    CsvReport csv_report;
    try
    {
        csv = read_labeled_csv(path, &csv_report);
    }
    catch (const std::exception& e)
    {
        input_error("--column: ", e.what());
    }
    if (csv_report.dropped > 0)
    {
        input_error("--column ", path, ": line ", csv_report.first_dropped_line,
                " is not a row of numbers. A text line is only allowed above the first row of "
                "data.");
    }
    ProfileColumns columns
            = find_profile_columns(csv, unit, path, ProfileColumns::Require::none, "--column");

    // The humidity's other forms.
    bool dewpoint_celsius = true;
    std::size_t dewpoint_column = csv.column({"td c", "dewpoint c", "dew point c", "td degc"});
    if (dewpoint_column == npos)
    {
        dewpoint_column = csv.column({"td k", "dewpoint k", "dew point k"});
        dewpoint_celsius = false;
    }
    double specific_to_kgkg = 1.0;
    std::size_t specific_column = csv.column({"q kgkg", "specific humidity kgkg"});
    if (specific_column == npos)
    {
        specific_column = csv.column({"q gkg", "specific humidity gkg"});
        specific_to_kgkg = 1e-3;
    }
    const int humidity_forms = (columns.humidity != npos) + (dewpoint_column != npos)
            + (specific_column != npos);
    if (humidity_forms > 1)
    {
        input_error("--column ", path, " gives the humidity more than one way (relative "
                "humidity, dew point, specific humidity); give one.");
    }
    if (dewpoint_column != npos && columns.temperature == npos)
        input_error("--column ", path, ": a dew point needs the temperature in the same rows.");
    if (specific_column != npos && columns.pressure == npos)
    {
        input_error("--column ", path, ": a specific humidity needs the pressure in the same "
                "rows.");
    }
    std::size_t needed = columns.needed();
    for (const std::size_t column : {dewpoint_column, specific_column})
    {
        if (column != npos)
        {
            columns.used[column] = true;
            needed = std::max(needed, column + 1);
        }
    }
    const bool has_wind = columns.speed != npos || columns.u != npos;
    if (columns.temperature == npos && humidity_forms == 0 && columns.pressure == npos
            && !has_wind)
    {
        input_error("--column ", path, " names no temperature, humidity, pressure or wind "
                "column.");
    }
    if (csv.rows.empty())
        input_error("--column ", path, " holds no numeric rows under its header.");
    for (std::size_t row = 0; row < csv.rows.size(); row++)
    {
        if (csv.rows[row].size() < needed)
        {
            input_error("--column ", path, ": data row ", row + 1, " has ",
                    csv.rows[row].size(), " fields, too few for the columns its header names.");
        }
        for (const double value : csv.rows[row])
        {
            if (!std::isfinite(value))
            {
                input_error("--column ", path, ": data row ", row + 1,
                        " holds a value that is not a number.");
            }
        }
    }

    ColumnSource source;
    source.path = path;
    auto is_missing = [&](double value) {
        return std::find(missing.begin(), missing.end(), value) != missing.end();
    };
    auto implausible = [&](std::size_t row, const std::string& what) {
        if (!drop_implausible)
        {
            input_error("--column ", path, ": data row ", row + 1, " gives ", what,
                    ", which is not plausible. List a fill value with --missing, or drop such "
                    "cells with --drop_implausible.");
        }
        source.implausible_dropped++;
    };
    auto height_of = [&](std::size_t row) {
        const double raw = csv.rows[row][columns.height] * unit.to_meters;
        return columns.above_sea_level ? raw - ground_elevation : raw;
    };

    // Temperature, with its unit told from the values when the header does
    // not say, over the cells that are not fill values.
    std::vector<Level> temperatures;
    std::vector<double> kelvin(csv.rows.size(), std::numeric_limits<double>::quiet_NaN());
    if (columns.temperature != npos)
    {
        bool celsius = columns.celsius;
        if (!columns.temperature_unit_named)
        {
            double low = std::numeric_limits<double>::infinity();
            double high = -low;
            for (const auto& row : csv.rows)
            {
                if (is_missing(row[columns.temperature]) || is_missing(row[columns.height]))
                    continue;
                low = std::min(low, row[columns.temperature]);
                high = std::max(high, row[columns.temperature]);
            }
            if (low >= kelvin_low && high <= kelvin_high)
                celsius = false;
            else if (low >= celsius_low && high <= celsius_high)
                celsius = true;
            else
            {
                input_error("--column ", path, ": temperatures from ", low, " to ", high,
                        " are neither kelvin nor Celsius. Name the column T_K or T_C.");
            }
        }
        for (std::size_t row = 0; row < csv.rows.size(); row++)
        {
            const double value = csv.rows[row][columns.temperature];
            if (is_missing(value) || is_missing(csv.rows[row][columns.height]))
            {
                source.sentinels_dropped++;
                continue;
            }
            const bool ok = celsius ? value >= celsius_low && value <= celsius_high
                                    : value >= kelvin_low && value <= kelvin_high;
            if (!ok)
            {
                std::ostringstream what;
                what << "a temperature of " << value << (celsius ? " C" : " K");
                implausible(row, what.str());
                continue;
            }
            kelvin[row] = celsius ? value + 273.15 : value;
            temperatures.push_back({height_of(row), kelvin[row], row});
        }
        source.temperature = series_from(path, temperatures);
    }

    std::vector<Level> pressures;
    std::vector<double> kpa(csv.rows.size(), std::numeric_limits<double>::quiet_NaN());
    if (columns.pressure != npos)
    {
        double to_kpa = columns.pressure_to_kpa;
        if (!columns.pressure_unit_named)
        {
            double high = 0.0;
            for (const auto& row : csv.rows)
            {
                if (!is_missing(row[columns.pressure]))
                    high = std::max(high, row[columns.pressure]);
            }
            if (high <= 120.0)
                to_kpa = 1.0;
            else if (high <= 1200.0)
                to_kpa = 0.1;
            else if (high <= 120000.0)
                to_kpa = 1e-3;
            else
            {
                input_error("--column ", path, ": pressure up to ", high,
                        " is not kPa, hPa or Pa. Name the column p_kPa, p_hPa or p_Pa.");
            }
        }
        for (std::size_t row = 0; row < csv.rows.size(); row++)
        {
            const double raw = csv.rows[row][columns.pressure];
            if (is_missing(raw) || is_missing(csv.rows[row][columns.height]))
            {
                source.sentinels_dropped++;
                continue;
            }
            const double value = raw * to_kpa;
            if (!(value > 0.0) || value > 120.0)
            {
                std::ostringstream what;
                what << "a pressure of " << value << " kPa";
                implausible(row, what.str());
                continue;
            }
            kpa[row] = value;
            pressures.push_back({height_of(row), value, row});
        }
        source.pressure = series_from(path, pressures);
    }

    if (humidity_forms == 1)
    {
        std::vector<Level> humidities;
        for (std::size_t row = 0; row < csv.rows.size(); row++)
        {
            const auto& data = csv.rows[row];
            if (is_missing(data[columns.height]))
            {
                source.sentinels_dropped++;
                continue;
            }
            double rh = 0.0;
            if (columns.humidity != npos)
            {
                if (is_missing(data[columns.humidity]))
                {
                    source.sentinels_dropped++;
                    continue;
                }
                rh = data[columns.humidity];
            }
            else if (dewpoint_column != npos)
            {
                // RH = 100 e_sat(Td) / e_sat(T).
                const double raw = data[dewpoint_column];
                if (is_missing(raw))
                {
                    source.sentinels_dropped++;
                    continue;
                }
                const double td = dewpoint_celsius ? raw + 273.15 : raw;
                if (!(td >= kelvin_low && td <= kelvin_high))
                {
                    std::ostringstream what;
                    what << "a dew point of " << raw << (dewpoint_celsius ? " C" : " K");
                    implausible(row, what.str());
                    continue;
                }
                if (!std::isfinite(kelvin[row]))
                    continue; // its temperature is missing or dropped, and counted there
                rh = 100.0 * saturation_pressure(td) / saturation_pressure(kelvin[row]);
            }
            else
            {
                // The vapor pressure e = q p / (0.622 + 0.378 q), then
                // RH = 100 e / e_sat(T).
                const double raw = data[specific_column];
                if (is_missing(raw))
                {
                    source.sentinels_dropped++;
                    continue;
                }
                const double q = raw * specific_to_kgkg;
                if (!(q >= 0.0 && q < 0.05))
                {
                    std::ostringstream what;
                    what << "a specific humidity of " << q << " kg/kg";
                    implausible(row, what.str());
                    continue;
                }
                if (!std::isfinite(kelvin[row]) || !std::isfinite(kpa[row]))
                {
                    if (columns.temperature == npos)
                    {
                        input_error("--column ", path, ": a specific humidity needs the "
                                "temperature in the same rows.");
                    }
                    continue;
                }
                const double e = q * kpa[row] / (0.622 + 0.378 * q);
                rh = 100.0 * e / saturation_pressure(kelvin[row]);
            }
            if (!(rh >= 0.0 && rh <= humidity_high))
            {
                std::ostringstream what;
                what << "a relative humidity of " << rh << " percent";
                implausible(row, what.str());
                continue;
            }
            humidities.push_back({height_of(row), rh, row});
        }
        source.humidity = series_from(path, humidities);
    }

    if (has_wind)
    {
        std::vector<Level> east, north;
        for (std::size_t row = 0; row < csv.rows.size(); row++)
        {
            const auto& data = csv.rows[row];
            double u = 0.0, v = 0.0;
            if (columns.speed != npos)
            {
                if (is_missing(data[columns.speed]) || is_missing(data[columns.direction])
                        || is_missing(data[columns.height]))
                {
                    source.sentinels_dropped++;
                    continue;
                }
                const double speed = data[columns.speed] * columns.speed_to_mps;
                const double direction = data[columns.direction];
                if (!(speed >= 0.0) || speed > 150.0)
                {
                    std::ostringstream what;
                    what << "a wind speed of " << speed << " m/s";
                    implausible(row, what.str());
                    continue;
                }
                if (!(direction >= 0.0 && direction <= 360.0))
                {
                    std::ostringstream what;
                    what << "a wind direction of " << direction << " degrees";
                    implausible(row, what.str());
                    continue;
                }
                // The meteorological convention, as read_atmosphere_profile
                // has it: the bearing the wind comes from.
                const double from = direction * pi / 180.0;
                u = -speed * std::sin(from);
                v = -speed * std::cos(from);
            }
            else
            {
                if (is_missing(data[columns.u]) || is_missing(data[columns.v])
                        || is_missing(data[columns.height]))
                {
                    source.sentinels_dropped++;
                    continue;
                }
                u = data[columns.u] * columns.component_to_mps;
                v = data[columns.v] * columns.component_to_mps;
                if (!(std::hypot(u, v) <= 150.0))
                {
                    std::ostringstream what;
                    what << "a wind of " << std::hypot(u, v) << " m/s";
                    implausible(row, what.str());
                    continue;
                }
            }
            east.push_back({height_of(row), u, row});
            north.push_back({height_of(row), v, row});
        }
        const auto e = series_from(path, east);
        const auto n = series_from(path, north);
        source.wind = WindSeries {path, e.z, e.value, n.value};
    }

    if (source.temperature)
        refuse_duplicates(*source.temperature, "temperature", unit);
    if (source.humidity)
        refuse_duplicates(*source.humidity, "humidity", unit);
    if (source.pressure)
        refuse_duplicates(*source.pressure, "pressure", unit);
    if (source.wind)
        refuse_duplicates(*source.wind, "wind", unit);

    std::ostringstream ignored;
    for (std::size_t i = 0; i < csv.header.size(); i++)
    {
        if (!columns.used[i])
            ignored << (ignored.tellp() > 0 ? ", " : "") << csv.header[i];
    }
    log("Column ", path, ":", source.temperature ? " temperature" : "",
            source.humidity ? " humidity" : "", source.pressure ? " pressure" : "",
            source.wind ? " wind" : "", "; ", source.sentinels_dropped, " fill values and ",
            source.implausible_dropped, " implausible values dropped",
            ignored.tellp() > 0 ? "; ignoring columns: " + ignored.str() : std::string(), ".");
    return source;
}

std::string build_rules_problem(const BuildRules& rules, LengthUnit unit)
{
    std::ostringstream problem;
    if (!std::isfinite(rules.roughness_length) || !(rules.roughness_length > 0.0)
            || rules.roughness_length > 3.0)
    {
        problem << "a roughness length of " << length_text(rules.roughness_length, unit)
                << " is not plausible: it runs from 0.0002 m over calm water to about 2 m over "
                << "city centers (more than 0 and up to 3 m are accepted).";
    }
    else if (!std::isfinite(rules.top) || !(rules.top > 0.0))
        problem << "the top (" << length_text(rules.top, unit) << ") has to be above the ground.";
    else if (rules.levels_per_decade < 1 || rules.levels_per_decade > 100)
    {
        problem << "levels per factor of ten in height run from 1 to 100, not "
                << rules.levels_per_decade << ".";
    }
    else if (rules.surface_layer_depth
            && (!std::isfinite(*rules.surface_layer_depth) || !(*rules.surface_layer_depth > 0.0)))
    {
        problem << "the surface layer depth (" << length_text(*rules.surface_layer_depth, unit)
                << ") has to be above the ground.";
    }
    else
    {
        for (const double z : rules.levels)
        {
            if (!std::isfinite(z) || z < 0.0)
            {
                problem << "a listed level (" << length_text(z, unit)
                        << ") is below the ground.";
                break;
            }
        }
    }
    return problem.str();
}

AtmosphereProfile build_profile(const MeasuredColumns& measured, const SurfaceAir& surface,
        const BuildRules& rules, LengthUnit unit, BuildReport& report)
{
    if (const auto problem = build_rules_problem(rules, unit); !problem.empty())
        input_error("niceops_build_atmosphere: ", problem);
    const double top = rules.top;

    // --- Temperature, humidity: their knots, starting at the ground. -------
    struct Knots
    {
        std::vector<double> z;
        std::vector<double> value;
        std::string source;
    };
    auto ground_knots = [&](const std::optional<ColumnSeries>& series,
                                const std::optional<double>& at_ground, const char* surface_name,
                                const char* name, bool logarithmic) {
        Knots knots;
        if (at_ground)
        {
            knots.z.push_back(0.0);
            knots.value.push_back(*at_ground);
            knots.source = surface_name;
            if (series)
            {
                for (std::size_t i = 0; i < series->z.size(); i++)
                {
                    if (series->z[i] > 0.0)
                    {
                        knots.z.push_back(series->z[i]);
                        knots.value.push_back(series->value[i]);
                    }
                    else
                        report.below_ground_dropped++;
                }
                knots.source += std::string(" and ") + series->source;
            }
            return knots;
        }
        if (!series || series->z.empty())
            return knots;
        knots.source = series->source;
        const auto& z = series->z;
        double ground;
        if (z.front() > 0.0)
        {
            ground = series->value.front();
            note(report, std::string("No ") + name + " at or below the ground in "
                            + series->source + ": the lowest, at " + length_text(z.front(), unit)
                            + ", is held down to it.");
        }
        else if (z.back() < 0.0)
        {
            ground = series->value.back();
            note(report, std::string("Every ") + name + " in " + series->source
                            + " is below the ground: the highest, at "
                            + length_text(z.back(), unit) + ", is held up from it.");
        }
        else
            ground = logarithmic ? interpolate_log(z, series->value, 0.0)
                                 : interpolate(z, series->value, 0.0);
        knots.z.push_back(0.0);
        knots.value.push_back(ground);
        for (std::size_t i = 0; i < z.size(); i++)
        {
            if (z[i] > 0.0)
            {
                knots.z.push_back(z[i]);
                knots.value.push_back(series->value[i]);
            }
            else if (z[i] < 0.0)
                report.below_ground_dropped++;
        }
        return knots;
    };

    const Knots temperature = ground_knots(measured.temperature, surface.temperature,
            "--temperature", "temperature", false);
    if (temperature.z.empty())
    {
        input_error("niceops_build_atmosphere: no temperature. Give --temperature, or a column "
                "with T_K or T_C.");
    }
    Knots humidity = ground_knots(measured.humidity, surface.humidity, "--humidity",
            "humidity", false);
    if (humidity.z.empty())
    {
        input_error("niceops_build_atmosphere: no humidity. Give --humidity, or a column with "
                "RH, a dew point or a specific humidity.");
    }
    for (std::size_t i = 0; i < humidity.z.size(); i++)
    {
        if (humidity.value[i] > 100.0)
        {
            if (!rules.clip_humidity)
            {
                input_error("niceops_build_atmosphere: ", humidity.source,
                        " gives a relative humidity of ", humidity.value[i], " percent at ",
                        length_text(humidity.z[i], unit), ", above saturation. Clip such values "
                        "to 100 with --clip_humidity.");
            }
            humidity.value[i] = 100.0;
            report.humidity_clipped++;
        }
    }
    if (report.humidity_clipped > 0)
    {
        note(report, "Clipped " + std::to_string(report.humidity_clipped)
                        + " relative humidities above 100 percent to 100.");
    }

    const bool surface_only_temperature = temperature.z.size() == 1;
    if (rules.stability != Stability::neutral
            && !(surface_only_temperature && rules.above == BuildRules::Above::lapse))
    {
        input_error("niceops_build_atmosphere: --stability sets the temperature gradient above a "
                "temperature measured at the ground alone, with --above lapse; this build "
                "has ", surface_only_temperature ? "--above hold" : "a measured temperature profile",
                ".");
    }

    // --- Wind. -------------------------------------------------------------
    std::vector<double> wind_z, wind_east, wind_north;
    // The bearing of the lowest level, when it is the surface station's.
    bool lowest_is_surface = false;
    std::string wind_source;
    {
        struct WindLevel
        {
            double z, east, north;
            bool surface;
        };
        std::vector<WindLevel> levels;
        if (measured.wind)
        {
            wind_source = measured.wind->source;
            for (std::size_t i = 0; i < measured.wind->z.size(); i++)
            {
                if (measured.wind->z[i] > 0.0)
                {
                    levels.push_back({measured.wind->z[i], measured.wind->east[i],
                            measured.wind->north[i], false});
                }
                else
                    report.below_ground_dropped++;
            }
        }
        if (surface.wind)
        {
            const SurfaceWind& station = *surface.wind;
            if (!(station.height > 0.0))
            {
                input_error("niceops_build_atmosphere: the anemometer height has to be above the "
                        "ground, not ", length_text(station.height, unit), ".");
            }
            const double from = station.direction * pi / 180.0;
            levels.push_back({station.height, -station.speed * std::sin(from),
                    -station.speed * std::cos(from), true});
            wind_source = wind_source.empty() ? "--surface_wind" : wind_source + " and --surface_wind";
        }
        std::stable_sort(levels.begin(), levels.end(),
                [](const WindLevel& a, const WindLevel& b) { return a.z < b.z; });
        for (std::size_t i = 1; i < levels.size(); i++)
        {
            if (!(levels[i].z > levels[i - 1].z))
            {
                input_error("niceops_build_atmosphere: ", wind_source, " give two winds at ",
                        length_text(levels[i].z, unit), ".");
            }
        }
        const bool all_calm = std::all_of(levels.begin(), levels.end(),
                [](const WindLevel& level) { return level.east == 0.0 && level.north == 0.0; });
        if (all_calm)
        {
            if (!levels.empty())
            {
                note(report, "The wind from " + wind_source + " is calm at every level: the "
                                "profile has no wind.");
            }
            levels.clear();
        }
        for (const auto& level : levels)
        {
            if (level.east == 0.0 && level.north == 0.0)
            {
                if (!rules.drop_calm)
                {
                    input_error("niceops_build_atmosphere: ", wind_source, " gives a calm wind at ",
                            length_text(level.z, unit), ", above the ground. The stratified ray "
                            "model takes a calm level above the ground as a thin layer that "
                            "puts levels 20 dB or more low; drop such levels with --calm_levels "
                            "drop.");
                }
                report.calm_dropped++;
                continue;
            }
            wind_z.push_back(level.z);
            wind_east.push_back(level.east);
            wind_north.push_back(level.north);
            if (wind_z.size() == 1)
                lowest_is_surface = level.surface;
        }
        if (report.calm_dropped > 0)
        {
            note(report, "Dropped " + std::to_string(report.calm_dropped)
                            + " calm wind levels above the ground.");
        }
    }
    const bool has_wind = !wind_z.empty();
    if (rules.surface_layer_depth && wind_z.size() > 1)
    {
        input_error("niceops_build_atmosphere: --surface_layer_depth carries the law above a "
                "single wind level; with ", wind_z.size(), " wind levels it would kink at the "
                "highest. Leave it out, and the highest is held.");
    }

    // The law below the lowest wind level: --wind's, through that level.
    SurfaceWind law;
    std::array<double, 2> law_direction {0.0, 0.0};
    if (has_wind)
    {
        law.height = wind_z.front();
        law.roughness_length = rules.roughness_length;
        law.stability = rules.stability;
        law.levels_per_decade = rules.levels_per_decade;
        law.surface_layer_depth = rules.surface_layer_depth.value_or(law.height);
        if (lowest_is_surface)
        {
            // The station's own speed and bearing, so a build from it alone is
            // --wind's profile to the bit.
            law.speed = surface.wind->speed;
            law.direction = surface.wind->direction;
            law_direction = {-std::sin(law.direction * pi / 180.0),
                    -std::cos(law.direction * pi / 180.0)};
        }
        else
        {
            law.speed = std::hypot(wind_east.front(), wind_north.front());
            law_direction = {wind_east.front() / law.speed, wind_north.front() / law.speed};
        }
        if (rules.fill != BuildRules::Fill::none && !(law.height > law.roughness_length))
        {
            input_error("niceops_build_atmosphere: the lowest wind level (",
                    length_text(law.height, unit), ") has to stand above the roughness length (",
                    length_text(law.roughness_length, unit), "), where the log law is zero.");
        }
        if (rules.surface_layer_depth && !(*rules.surface_layer_depth > law.height))
        {
            input_error("niceops_build_atmosphere: the surface layer (",
                    length_text(*rules.surface_layer_depth, unit), " deep) has to reach above the "
                    "wind level (", length_text(law.height, unit), ").");
        }
        if (law.speed_at(top) > 150.0)
        {
            input_error("niceops_build_atmosphere: the law carries the wind to ", law.speed_at(top),
                    " m/s, which is not plausible.");
        }
    }
    const bool law_continues = has_wind && rules.surface_layer_depth.has_value();

    // --- Levels. -----------------------------------------------------------
    std::vector<double> measured_levels;
    auto add_measured = [&](const std::vector<double>& z) {
        for (const double level : z)
        {
            if (level > 0.0 && level < top)
                measured_levels.push_back(level);
        }
    };
    add_measured(temperature.z);
    add_measured(humidity.z);
    if (measured.pressure)
        add_measured(measured.pressure->z);
    add_measured(wind_z);

    AtmosphereProfile profile;
    profile.source = "niceops_build_atmosphere";
    profile.height.push_back(0.0);
    if (rules.levels.empty())
    {
        std::vector<double> anchors = measured_levels;
        if (has_wind && rules.fill != BuildRules::Fill::none)
            anchors.push_back(law.first_level());
        if (law_continues && *rules.surface_layer_depth < top)
            anchors.push_back(*rules.surface_layer_depth);
        anchors.push_back(top);
        const auto levels = geometric_levels(anchors, rules.levels_per_decade);
        profile.height.insert(profile.height.end(), levels.begin(), levels.end());
    }
    else
    {
        std::vector<double> levels = measured_levels;
        for (const double z : rules.levels)
        {
            if (z > 0.0 && z < top)
                levels.push_back(z);
        }
        levels.push_back(top);
        std::sort(levels.begin(), levels.end());
        levels.erase(std::unique(levels.begin(), levels.end(),
                             [](double a, double b) { return b - a <= 1e-9 * b; }),
                levels.end());
        profile.height.insert(profile.height.end(), levels.begin(), levels.end());
    }

    // --- Values at the levels. --------------------------------------------
    // --wind's air, for a temperature from the surface alone: the class's
    // gradient through the surface layer (200 m unless given), dry adiabatic
    // above it.
    SurfaceWind air_law;
    air_law.stability = rules.stability;
    const double layer = rules.surface_layer_depth.value_or(air_law.surface_layer_depth);
    const double gradient = air_law.temperature_gradient();
    for (const double z : profile.height)
    {
        double t;
        if (rules.above == BuildRules::Above::lapse && z > temperature.z.back())
        {
            if (surface_only_temperature)
            {
                // Spelled as synthetic_profile spells it.
                t = temperature.value.front() + gradient * std::min(z, layer)
                        + dry_adiabatic_lapse * std::max(z - layer, 0.0);
            }
            else
                t = temperature.value.back() + dry_adiabatic_lapse * (z - temperature.z.back());
        }
        else
            t = interpolate(temperature.z, temperature.value, z);
        profile.temperature.push_back(t);
        profile.relative_humidity.push_back(interpolate(humidity.z, humidity.value, z));

        if (!has_wind)
            continue;
        double east = 0.0, north = 0.0;
        if (z > 0.0 && (law_continues || z < law.height))
        {
            if (rules.fill == BuildRules::Fill::log || law_continues)
            {
                const double speed = law.speed_at(z);
                east = speed * law_direction[0];
                north = speed * law_direction[1];
            }
            else if (rules.fill == BuildRules::Fill::hold)
            {
                east = wind_east.front();
                north = wind_north.front();
            }
            else
            {
                east = wind_east.front() * (z / law.height);
                north = wind_north.front() * (z / law.height);
            }
        }
        else if (z > 0.0)
        {
            east = interpolate(wind_z, wind_east, z);
            north = interpolate(wind_z, wind_north, z);
        }
        if (z > 0.0 && east == 0.0 && north == 0.0)
        {
            input_error("niceops_build_atmosphere: the level at ", length_text(z, unit),
                    " would be calm: it is at or below the roughness length (",
                    length_text(law.roughness_length, unit), "), where the log law is zero. "
                    "List no level there.");
        }
        // + 0.0 turns the ground level's -0 into 0 in the written profile.
        profile.wind_east.push_back(east + 0.0);
        profile.wind_north.push_back(north + 0.0);
    }
    if (has_wind)
    {
        if (profile.height.size() > 1 && profile.height[1] < law.height)
            report.filled_wind.emplace_back(0.0, law.height);
        if (!law_continues && wind_z.back() < top)
            report.filled_wind.emplace_back(wind_z.back(), top);
    }

    // --- Pressure. ---------------------------------------------------------
    std::optional<double> ground_pressure = surface.pressure;
    const bool held_concentration
            = rules.humidity_aloft == BuildRules::HumidityAloft::held_concentration;
    if (measured.pressure && !measured.pressure->z.empty())
    {
        // Knots: the ground's (from the station, or the measurements across
        // it) and every measurement above it.  Outside them, hydrostatic from
        // the end knot, which is a level.
        const auto& series = *measured.pressure;
        std::vector<double> z, p;
        if (surface.pressure)
        {
            z.push_back(0.0);
            p.push_back(*surface.pressure);
        }
        else if (series.z.front() <= 0.0)
        {
            z.push_back(0.0);
            if (series.z.back() >= 0.0)
                p.push_back(interpolate_log(series.z, series.value, 0.0));
            else
            {
                // Every measurement below the ground: hydrostatic up from the
                // highest, the temperature held at the ground's below it, as
                // with_hydrostatic_pressure holds it.
                p.push_back(series.value.back()
                        * std::exp(-standard_gravity / dry_air_gas_constant * (0.0 - series.z.back())
                                / profile.temperature.front()));
            }
        }
        for (std::size_t i = 0; i < series.z.size(); i++)
        {
            if (series.z[i] > 0.0)
            {
                z.push_back(series.z[i]);
                p.push_back(series.value[i]);
            }
            else if (surface.pressure)
                report.below_ground_dropped++;
        }
        const auto relative = profile.with_hydrostatic_pressure(1.0).pressure;
        auto level_of = [&](double height) {
            const auto it = std::lower_bound(profile.height.begin(), profile.height.end(), height);
            std::size_t i = static_cast<std::size_t>(it - profile.height.begin());
            if (i == profile.height.size()
                    || (i > 0 && height - profile.height[i - 1] < profile.height[i] - height))
                i--;
            return i;
        };
        const std::size_t front = level_of(z.front());
        const std::size_t back = level_of(std::min(z.back(), top));
        profile.pressure.resize(profile.height.size());
        for (std::size_t i = 0; i < profile.height.size(); i++)
        {
            const double h = profile.height[i];
            if (h < z.front())
                profile.pressure[i] = p.front() * relative[i] / relative[front];
            else if (h > z.back())
                profile.pressure[i] = p.back() * relative[i] / relative[back];
            else
                profile.pressure[i] = interpolate_log(z, p, h);
        }
        for (std::size_t i = 1; i < profile.pressure.size(); i++)
        {
            if (profile.pressure[i] > profile.pressure[i - 1])
            {
                input_error("niceops_build_atmosphere: the pressure from ", series.source,
                        " rises with height at ", length_text(profile.height[i], unit), ".");
            }
        }
    }
    else if (ground_pressure && (rules.write_pressure || held_concentration))
        profile = profile.with_hydrostatic_pressure(*ground_pressure);
    else if (rules.write_pressure || held_concentration)
    {
        input_error("niceops_build_atmosphere: ", held_concentration
                        ? "--humidity_aloft held_concentration" : "--write_pressure",
                " needs a pressure: give --pressure, or a column with one.");
    }

    // --- Humidity above its highest measurement. ---------------------------
    if (held_concentration)
    {
        // The molar concentration of water vapor at the highest humidity
        // measurement is held above it, as a well-mixed layer holds it, and
        // the relative humidity follows the cooling air; capped at 100.
        const double highest = humidity.z.back();
        std::size_t from = 0;
        while (from + 1 < profile.height.size() && profile.height[from + 1] <= highest)
            from++;
        Atmosphere air;
        air.temperature = profile.temperature[from];
        air.relative_humidity = profile.relative_humidity[from];
        air.pressure = profile.pressure[from];
        const double concentration = air.molar_concentration_water_vapor();
        for (std::size_t i = from + 1; i < profile.height.size(); i++)
        {
            Atmosphere level;
            level.temperature = profile.temperature[i];
            level.pressure = profile.pressure[i];
            profile.relative_humidity[i]
                    = std::min(100.0, concentration * level.pressure / level.saturation_pressure());
        }
    }

    // --- What was measured. ------------------------------------------------
    auto span = [&](const char* name, const std::vector<double>& z) {
        if (!z.empty())
            report.measured_span[name] = {z.front(), z.back()};
    };
    if (measured.temperature && !measured.temperature->z.empty())
        span("temperature", measured.temperature->z);
    if (measured.humidity && !measured.humidity->z.empty())
        span("humidity", measured.humidity->z);
    if (measured.pressure && !measured.pressure->z.empty())
        span("pressure", measured.pressure->z);
    span("wind", wind_z);
    if (report.below_ground_dropped > 0)
    {
        log("Dropped ", report.below_ground_dropped, " measurements at or below the ground once "
                "they had set its value.");
    }

    // The run's own checks, so anything it would refuse is refused here.
    for (std::size_t i = 0; i < profile.height.size(); i++)
    {
        const double t = profile.temperature[i];
        const double rh = profile.relative_humidity[i];
        if (!(t >= kelvin_low && t <= kelvin_high))
        {
            input_error("niceops_build_atmosphere: the temperature at ",
                    length_text(profile.height[i], unit), " comes to ", t, " K, which is not "
                    "plausible.");
        }
        if (!(rh >= 0.0 && rh <= 100.0))
        {
            input_error("niceops_build_atmosphere: the relative humidity at ",
                    length_text(profile.height[i], unit), " comes to ", rh, " percent.");
        }
    }
    if (*std::max_element(profile.relative_humidity.begin(), profile.relative_humidity.end())
            <= 1.0)
    {
        input_error("niceops_build_atmosphere: the relative humidity never exceeds 1 percent; it "
                "is read in percent, and --atmosphere would refuse the profile.");
    }
    return profile;
}

namespace {

// One run's gates: ln z and the speed there.
struct RunProfile
{
    std::vector<double> log_z;
    std::vector<double> speed;
};

// The slope A_j and the residual sum for one run at ln z0 = ell.
struct RunTerms
{
    double slope;
    double squares;
    double gradient; // A_j * sum_i (U_ij - A_j x_i): half of dS_j/d ell
};

RunTerms run_terms(const RunProfile& run, double ell)
{
    double xu = 0.0, xx = 0.0;
    for (std::size_t i = 0; i < run.log_z.size(); i++)
    {
        const double x = run.log_z[i] - ell;
        xu += x * run.speed[i];
        xx += x * x;
    }
    const double slope = xu / xx;
    double squares = 0.0, residual_sum = 0.0;
    for (std::size_t i = 0; i < run.log_z.size(); i++)
    {
        const double r = run.speed[i] - slope * (run.log_z[i] - ell);
        squares += r * r;
        residual_sum += r;
    }
    return {slope, squares, slope * residual_sum};
}

constexpr double fit_low = -9.210340371976182; // ln 1e-4 m
const double fit_high = std::log(3.0); // ln 3 m

// ln z0 minimizing the pooled squares over `runs` (indices into `profiles`),
// and whether it came to an end of the search.  A scan brackets the minimum;
// bisection on dS/d ell = 0 (by the envelope theorem, 2 sum_j A_j sum_i r_ij)
// then takes it to the last bit, so exact profiles give z0 back to rounding.
std::pair<double, bool> fit_log_z0(const std::vector<RunProfile>& profiles,
        const std::vector<std::size_t>& runs)
{
    auto squares = [&](double ell) {
        double sum = 0.0;
        for (const std::size_t j : runs)
            sum += run_terms(profiles[j], ell).squares;
        return sum;
    };
    auto gradient = [&](double ell) {
        double sum = 0.0;
        for (const std::size_t j : runs)
            sum += run_terms(profiles[j], ell).gradient;
        return sum;
    };
    constexpr int scan = 240;
    int best = 0;
    double best_value = std::numeric_limits<double>::infinity();
    for (int k = 0; k <= scan; k++)
    {
        const double ell = fit_low + (fit_high - fit_low) * k / scan;
        const double value = squares(ell);
        if (value < best_value)
        {
            best_value = value;
            best = k;
        }
    }
    double lo = fit_low + (fit_high - fit_low) * std::max(best - 1, 0) / scan;
    double hi = fit_low + (fit_high - fit_low) * std::min(best + 1, scan) / scan;
    // The squares fall while the gradient is negative (dS/d ell = 2 g).  At
    // an end of the search there is no sign change to bisect.
    if (!(gradient(lo) < 0.0) || !(gradient(hi) > 0.0))
        return {fit_low + (fit_high - fit_low) * best / scan, best == 0 || best == scan};
    for (int iteration = 0; iteration < 200; iteration++)
    {
        const double mid = 0.5 * (lo + hi);
        if (mid == lo || mid == hi)
            break;
        if (gradient(mid) < 0.0)
            lo = mid;
        else
            hi = mid;
    }
    return {0.5 * (lo + hi), false};
}

double percentile(std::vector<double> values, double q)
{
    std::sort(values.begin(), values.end());
    const double position = q * static_cast<double>(values.size() - 1);
    const auto below = static_cast<std::size_t>(std::floor(position));
    const std::size_t above = std::min(below + 1, values.size() - 1);
    const double t = position - static_cast<double>(below);
    return values[below] + t * (values[above] - values[below]);
}

} // namespace

RoughnessFit fit_roughness_length(const std::vector<ColumnSource>& runs,
        const std::vector<std::string>& labels, const FitRules& rules)
{
    if (runs.size() != labels.size())
        throw std::logic_error("fit_roughness_length needs one label per run");
    if (!(rules.gate_low > 0.0) || !(rules.gate_high > rules.gate_low))
        input_error("fit_z0: the gates run from a height above the ground to a higher one.");
    if (rules.draws < 100)
        input_error("fit_z0: a bootstrap needs at least 100 draws, not ", rules.draws, ".");

    RoughnessFit fit;
    std::vector<RunProfile> profiles;
    for (std::size_t j = 0; j < runs.size(); j++)
    {
        const ColumnSource& run = runs[j];
        auto exclude = [&](const std::string& why) { fit.excluded.emplace_back(labels[j], why); };
        if (!run.wind)
        {
            exclude("no wind");
            continue;
        }
        // Gates in the band, in height order.
        std::vector<std::size_t> gates;
        for (std::size_t i = 0; i < run.wind->z.size(); i++)
        {
            const double z = run.wind->z[i];
            if (z >= rules.gate_low * (1.0 - 1e-9) && z <= rules.gate_high * (1.0 + 1e-9))
                gates.push_back(i);
        }
        std::sort(gates.begin(), gates.end(),
                [&](std::size_t a, std::size_t b) { return run.wind->z[a] < run.wind->z[b]; });
        if (gates.size() < 2)
        {
            exclude("fewer than two wind gates between the gate heights");
            continue;
        }
        const std::size_t low = gates.front();
        const std::size_t high = gates.back();
        const double z_low = run.wind->z[low];
        const double z_high = run.wind->z[high];
        const double u_low = std::hypot(run.wind->east[low], run.wind->north[low]);
        const double u_high = std::hypot(run.wind->east[high], run.wind->north[high]);
        if (u_high < rules.min_speed)
        {
            std::ostringstream why;
            why << "light wind: " << u_high << " m/s at the highest gate, under " << rules.min_speed;
            exclude(why.str());
            continue;
        }
        double veer = 0.0;
        for (const std::size_t i : gates)
        {
            const double cross = run.wind->east[low] * run.wind->north[i]
                    - run.wind->north[low] * run.wind->east[i];
            const double dot = run.wind->east[low] * run.wind->east[i]
                    + run.wind->north[low] * run.wind->north[i];
            veer = std::max(veer, std::fabs(std::atan2(cross, dot)) * 180.0 / pi);
        }
        if (veer > rules.max_veer)
        {
            std::ostringstream why;
            why << "the wind turns " << veer << " degrees across the gates, more than "
                << rules.max_veer;
            exclude(why.str());
            continue;
        }
        // The bulk Richardson number between the lowest and highest gate:
        // (g / theta) (d theta / dz) / (dU / dz)^2, theta = T + 0.0098 z.
        const auto& temperature = run.temperature;
        const bool screened = temperature && !temperature->z.empty()
                && temperature->z.front() <= z_low && temperature->z.back() >= z_high;
        if (screened)
        {
            const double t_low = interpolate(temperature->z, temperature->value, z_low);
            const double t_high = interpolate(temperature->z, temperature->value, z_high);
            const double theta_low = t_low - dry_adiabatic_lapse * z_low;
            const double theta_high = t_high - dry_adiabatic_lapse * z_high;
            const double dz = z_high - z_low;
            const double shear = (u_high - u_low) / dz;
            if (!(shear > 0.0))
            {
                exclude("the speed does not grow with height across the gates");
                continue;
            }
            const double richardson = standard_gravity / (0.5 * (theta_low + theta_high))
                    * ((theta_high - theta_low) / dz) / (shear * shear);
            if (std::fabs(richardson) > rules.max_richardson)
            {
                std::ostringstream why;
                why << "not neutral: bulk Richardson number " << richardson << ", beyond "
                    << rules.max_richardson;
                exclude(why.str());
                continue;
            }
        }
        else if (!rules.assume_neutral)
        {
            exclude("no temperature at both gate heights to screen it (--assume_neutral keeps "
                    "it)");
            continue;
        }
        else
        {
            warn("fit_z0: ", labels[j], " has no temperature at both gate heights; kept as "
                    "neutral (--assume_neutral).");
        }

        RunProfile profile;
        for (const std::size_t i : gates)
        {
            profile.log_z.push_back(std::log(run.wind->z[i]));
            profile.speed.push_back(std::hypot(run.wind->east[i], run.wind->north[i]));
        }
        profiles.push_back(profile);
        fit.used.push_back(labels[j]);
        // The run alone: U = a + b ln z, so ln z0 = -a / b.
        double mx = 0.0, mu = 0.0;
        for (std::size_t i = 0; i < profile.log_z.size(); i++)
        {
            mx += profile.log_z[i];
            mu += profile.speed[i];
        }
        mx /= static_cast<double>(profile.log_z.size());
        mu /= static_cast<double>(profile.log_z.size());
        double sxy = 0.0, sxx = 0.0;
        for (std::size_t i = 0; i < profile.log_z.size(); i++)
        {
            sxy += (profile.log_z[i] - mx) * (profile.speed[i] - mu);
            sxx += (profile.log_z[i] - mx) * (profile.log_z[i] - mx);
        }
        const double b = sxy / sxx;
        fit.per_run_z0.push_back(b > 0.0 ? std::exp(mx - mu / b)
                                         : std::numeric_limits<double>::quiet_NaN());
    }
    if (fit.used.size() < rules.min_runs)
    {
        input_error("fit_z0: ", fit.used.size(), " runs pass the screen, fewer than the ",
                rules.min_runs, " a fit needs (--min_runs); ", fit.excluded.size(),
                " were left out.");
    }

    std::vector<std::size_t> all(profiles.size());
    std::iota(all.begin(), all.end(), 0);
    const auto [ell, at_bound] = fit_log_z0(profiles, all);
    fit.z0 = std::exp(ell);
    fit.at_bound = at_bound;
    double squares = 0.0;
    std::size_t points = 0;
    for (const auto& profile : profiles)
    {
        squares += run_terms(profile, ell).squares;
        points += profile.log_z.size();
    }
    fit.residual_rms = std::sqrt(squares / static_cast<double>(points));
    std::vector<double> finite;
    for (const double z0 : fit.per_run_z0)
    {
        if (std::isfinite(z0))
            finite.push_back(z0);
    }
    fit.median_run_z0 = finite.empty() ? std::numeric_limits<double>::quiet_NaN()
                                       : percentile(finite, 0.5);

    // The interval: runs drawn with replacement, the fit repeated.
    std::mt19937_64 engine(rules.seed);
    std::vector<double> draws;
    draws.reserve(static_cast<std::size_t>(rules.draws));
    std::vector<std::size_t> sample(profiles.size());
    for (int draw = 0; draw < rules.draws; draw++)
    {
        for (auto& index : sample)
            index = portable_index(engine, profiles.size());
        draws.push_back(std::exp(fit_log_z0(profiles, sample).first));
    }
    fit.z0_low = percentile(draws, 0.025);
    fit.z0_high = percentile(draws, 0.975);
    return fit;
}

} // namespace niceops
