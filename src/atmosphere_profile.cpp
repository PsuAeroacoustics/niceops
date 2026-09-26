#include "atmosphere_profile.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "util/csv.hpp"
#include "util/logger.hpp"

namespace niceops {

namespace {

// Standard gravity and the dry-air gas constant.  The latter is the value
// trajectory.cpp takes pressure from density with, so the two agree.
constexpr double standard_gravity = 9.80665; // m/s^2
constexpr double dry_air_gas_constant = 287.05; // J/(kg K)

// Plausible ranges, in the stored units, for telling a unit from the values
// and for catching a column labelled with the wrong one.  They do not overlap:
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
    LabelledCsv csv;
    try
    {
        csv = read_labelled_csv(path);
    }
    catch (const std::exception& e)
    {
        input_error("--atmosphere: ", e.what());
    }

    const bool metric = unit.to_metres == 1.0;
    std::vector<bool> used(csv.header.size(), false);
    auto take = [&used](std::size_t column) {
        if (column != LabelledCsv::npos)
            used[column] = true;
        return column;
    };

    // Height.  As with a geodetic track's altitude, only the unit spellings
    // that agree with --units are matched: a column that says "z_m" while the
    // run speaks feet fails rather than being read three times too large.
    bool above_sea_level = false;
    std::size_t height_column = csv.column({"z", "agl", "z agl", "height agl"});
    if (height_column == LabelledCsv::npos)
    {
        height_column = metric ? csv.column({"z m", "agl m", "z agl m", "height agl m"})
                               : csv.column({"z ft", "agl ft", "z agl ft", "height agl ft"});
    }
    if (height_column == LabelledCsv::npos)
    {
        above_sea_level = true;
        height_column = csv.column({"alt", "altitude", "alt msl", "altitude msl", "height",
                                          "height msl"});
        if (height_column == LabelledCsv::npos)
        {
            height_column = metric
                    ? csv.column({"alt m", "altitude m", "alt msl m", "height m", "height msl m"})
                    : csv.column({"alt ft", "altitude ft", "alt msl ft", "height ft",
                                        "height msl ft"});
        }
    }
    if (height_column == LabelledCsv::npos)
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
    if (temperature_column == LabelledCsv::npos)
    {
        temperature_column = csv.column({"t c", "temp c", "temperature c", "t degc",
                                               "temperature degc"});
        celsius = temperature_column != LabelledCsv::npos;
    }
    if (temperature_column == LabelledCsv::npos)
    {
        temperature_column = csv.column({"t", "temp", "temperature"});
        temperature_unit = Unit::inferred;
    }
    if (temperature_column == LabelledCsv::npos)
    {
        input_error("--atmosphere ", path, " names no temperature column. Give \"T_K\" or "
                "\"T_C\".");
    }
    take(temperature_column);

    const std::size_t humidity_column
            = take(csv.column({"rh", "relative humidity", "humidity", "rh pct", "rh percent"}));
    if (humidity_column == LabelledCsv::npos)
    {
        input_error("--atmosphere ", path, " names no relative humidity column. Give \"RH\", "
                "in percent.");
    }

    // Pressure is optional: without it, it is integrated hydrostatically from
    // the surface pressure (with_hydrostatic_pressure).
    Unit pressure_unit = Unit::named;
    double pressure_to_kpa = 1.0;
    std::size_t pressure_column = csv.column({"p kpa", "pressure kpa", "pres kpa"});
    if (pressure_column == LabelledCsv::npos)
    {
        pressure_column = csv.column({"p hpa", "pressure hpa", "pres hpa", "p mb", "pressure mb",
                                            "p mbar", "pressure mbar"});
        if (pressure_column != LabelledCsv::npos)
            pressure_to_kpa = 0.1;
    }
    if (pressure_column == LabelledCsv::npos)
    {
        pressure_column = csv.column({"p pa", "pressure pa", "pres pa"});
        if (pressure_column != LabelledCsv::npos)
            pressure_to_kpa = 1e-3;
    }
    if (pressure_column == LabelledCsv::npos)
    {
        pressure_column = csv.column({"p", "pressure", "pres"});
        pressure_unit = Unit::inferred;
    }
    take(pressure_column);
    const bool have_pressure = pressure_column != LabelledCsv::npos;

    const std::size_t needed = 1
            + std::max({height_column, temperature_column, humidity_column,
                    have_pressure ? pressure_column : 0});
    if (csv.rows.empty())
        input_error("--atmosphere ", path, " holds no numeric rows under its header.");

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
        const double raw = csv.rows[row][height_column] * unit.to_metres;
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
    }
    for (std::size_t i = 1; i < profile.height.size(); i++)
    {
        if (!(profile.height[i] > profile.height[i - 1]))
        {
            input_error("--atmosphere ", path, " gives two levels at ",
                    profile.height[i] / unit.to_metres, " ", unit.name, ".");
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
                        profile.height[i] / unit.to_metres, " ", unit.name, ".");
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
        levels << profile.height.size() << " levels from " << profile.bottom() / unit.to_metres
               << " to " << profile.top() / unit.to_metres << " " << unit.name
               << " above the observers";
    log("Atmosphere profile ", path, ": ", levels.str(),
            above_sea_level ? " (read as altitude, less --ground_elevation)" : "",
            "; temperature in ", celsius ? "Celsius" : "kelvin",
            temperature_unit == Unit::inferred ? " (told from the values)" : "",
            have_pressure ? (pressure_unit == Unit::inferred ? "; pressure unit told from the values"
                                                             : "")
                          : "; no pressure column, so it is integrated hydrostatically from the "
                            "surface pressure",
            ".");
    if (ignored.tellp() > 0)
    {
        log("  ignoring columns: ", ignored.str(),
                ". Wind, if among them, is not modelled yet.");
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

    // A metre apart: alpha changes by parts in a thousand over that, so the
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
