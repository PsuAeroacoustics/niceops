#include "trajectory.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <map>
#include <sstream>

#include <nlohmann/json.hpp>

#include "smoothing_spline.hpp"
#include "util/constants.hpp"
#include "util/csv.hpp"
#include "util/logger.hpp"
#include "util/paths.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

// The vehicle JSONs store every value as a quoted string, and the D original
// read them with `.str` followed by `to!double`.  nlohmann's get<double>()
// throws on a string, so convert explicitly (tolerating a real JSON number).
double json_to_double(const nlohmann::json& value)
{
    if (value.is_string())
        return std::stod(value.get<std::string>());
    return value.get<double>();
}

// Where the labeled WGS84 format's columns sit in a particular file.  That
// format is read by column name rather than by position, so whatever a planner
// already exports can be used without being reordered first, and a file carries
// only the columns it actually has.
// The names a WGS84 track's heading column answers to.
const std::vector<const char*> heading_aliases = {"heading", "hdg", "true heading", "heading true", "heading deg"};

struct GeodeticColumns
{
    std::size_t time = LabeledCsv::npos;
    std::size_t utc = LabeledCsv::npos;
    std::size_t latitude = LabeledCsv::npos;
    std::size_t longitude = LabeledCsv::npos;
    std::size_t altitude = LabeledCsv::npos;
    std::size_t east_speed = LabeledCsv::npos;
    std::size_t north_speed = LabeledCsv::npos;
    std::size_t up_speed = LabeledCsv::npos;
    // Optional: the nose's compass heading, degrees from true north.
    std::size_t heading = LabeledCsv::npos;

    // A file is a WGS84 track when it says where the aircraft is.  Time can
    // come from either of two columns and velocity is optional, so neither of
    // those decides it.
    bool positioned() const
    {
        return latitude != LabeledCsv::npos && longitude != LabeledCsv::npos
                && altitude != LabeledCsv::npos;
    }

    // Velocity is taken from the file only when all three components are there;
    // two of them plus a differentiated third would mix two different estimates
    // of the same vector.
    bool velocity_given() const
    {
        return east_speed != LabeledCsv::npos && north_speed != LabeledCsv::npos
                && up_speed != LabeledCsv::npos;
    }

    // The widest index this file will be read at, so a short row is caught
    // where it can be named rather than as an out-of-range read.
    std::size_t width() const
    {
        std::size_t widest = 0;
        for (std::size_t index : {time, utc, latitude, longitude, altitude, east_speed,
                     north_speed, up_speed, heading})
        {
            if (index != LabeledCsv::npos)
                widest = std::max(widest, index + 1);
        }
        return widest;
    }
};

// A header that declares a unit -- "Alt (ft)", "Climb Rate [fpm]", "Time (ms)" --
// is checked against what this run reads that column in.  LabeledCsv::normalize
// drops the parenthesized part so the name matches whatever unit it carries,
// which left the declaration unchecked: a foot value read as meters put the
// aircraft 3.28 times too high, and a climb rate in feet per minute was read as
// feet per second, 60 times too large, with exit status 0.  Only tokens that are
// recognized and conflict are refused; an unfamiliar one ("deg", "msl") passes.
void check_declared_units(const std::string& filename, const LabeledCsv& csv,
        const GeodeticColumns& columns, LengthUnit unit)
{
    const bool metric = unit.to_meters == 1.0;
    const auto in = [](const std::string& token, std::initializer_list<const char*> list) {
        for (const char* candidate : list)
            if (token == candidate)
                return true;
        return false;
    };
    const auto check = [&](std::size_t index, const char* what, bool conflicts,
                               const std::string& expected) {
        if (index == LabeledCsv::npos || !conflicts)
            return;
        input_error("Column \"", csv.header[index], "\" of ", filename, " declares its ", what,
                " in ", csv.unit_token(index), ", but this run reads it in ", expected,
                ". Nothing is converted: give the column in ", expected,
                ", or change --units for the lengths.");
    };
    const auto length_conflicts = [&](const std::string& t) {
        return metric ? in(t, {"ft", "feet", "foot"}) : in(t, {"m", "meter", "meters", "meter", "meters"});
    };
    const auto time_conflicts = [&](const std::string& t) {
        return in(t, {"ms", "msec", "millisecond", "milliseconds", "us", "min", "mins", "minute",
                "minutes", "h", "hr", "hrs", "hour", "hours"});
    };
    const auto speed_conflicts = [&](const std::string& t) {
        return in(t, {"fpm", "ft/min", "m/min", "kt", "kts", "knot", "knots", "mph", "km/h", "kph",
                       "kmh", "kmph"})
                || (metric ? in(t, {"ft/s", "fps", "feet/s", "ft/sec"})
                           : in(t, {"m/s", "mps", "meters/s", "m/sec"}));
    };
    const std::string seconds = "seconds";
    const std::string speed_unit = metric ? "m/s" : "ft/s";
    if (columns.altitude != LabeledCsv::npos)
        check(columns.altitude, "altitude", length_conflicts(csv.unit_token(columns.altitude)), unit.name);
    if (columns.time != LabeledCsv::npos)
        check(columns.time, "time", time_conflicts(csv.unit_token(columns.time)), seconds);
    if (columns.utc != LabeledCsv::npos)
        check(columns.utc, "time", time_conflicts(csv.unit_token(columns.utc)), seconds);
    for (std::size_t index : {columns.east_speed, columns.north_speed, columns.up_speed})
        if (index != LabeledCsv::npos)
            check(index, "velocity", speed_conflicts(csv.unit_token(index)), speed_unit);
    // The heading is held to more than a deny-list: a magnetic heading is
    // degrees too, and read as true it turns the nose by the declination (about
    // 12 degrees at the 2017 site, worth up to 2 dB) with nothing said.  So any
    // unit but degrees is refused, a header naming a magnetic reference is
    // refused, and so is a file offering two heading columns.
    //
    // A heading whose reference or unit is spelled into the name itself
    // ("magnetic heading", "heading_rad") matches no alias, and would be
    // passed over for the velocity's nose without a word; so every column
    // naming a heading is looked at, not only the one that matched.
    for (std::size_t i = 0; i < csv.header.size(); i++)
    {
        const std::string name = LabeledCsv::normalize(csv.header[i]);
        if (name.find("heading") == std::string::npos && name.find("hdg") == std::string::npos)
            continue;
        if (name.find("mag") != std::string::npos)
            input_error("Column \"", csv.header[i], "\" of ", filename, " is a magnetic heading, "
                    "and the heading is read in degrees from true north. Add the site's declination "
                    "and name the column heading.");
        // The unit as a word of the header: "heading_rad", "hdg turns".
        std::string word;
        for (const char c : csv.header[i] + " ")
        {
            if (std::isalnum(static_cast<unsigned char>(c)))
            {
                word.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
                continue;
            }
            if (in(word, {"rad", "rads", "radian", "radians", "turn", "turns", "grad", "grads", "gon", "mil", "mils"}))
                input_error("Column \"", csv.header[i], "\" of ", filename, " declares its heading in ", word,
                        ", but the heading is read in degrees from true north. Nothing is converted.");
            word.clear();
        }
    }
    if (columns.heading != LabeledCsv::npos)
    {
        std::string raw = csv.header[columns.heading];
        std::transform(raw.begin(), raw.end(), raw.begin(),
                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (raw.find("mag") != std::string::npos)
            input_error("Column \"", csv.header[columns.heading], "\" of ", filename, " is a magnetic heading, "
                    "and the heading is read in degrees from true north. Add the site's declination "
                    "and name the column heading.");
        const std::string token = csv.unit_token(columns.heading);
        if (!token.empty() && !in(token, {"deg", "degs", "degree", "degrees", "degt", "true", "t"}))
            input_error("Column \"", csv.header[columns.heading], "\" of ", filename, " declares its heading in ",
                    token, ", but the heading is read in degrees from true north. Nothing is converted.");
        std::vector<std::string> candidates;
        for (std::size_t i = 0; i < csv.header.size(); i++)
            if (csv.column_matches(i, heading_aliases))
                candidates.push_back(csv.header[i]);
        if (candidates.size() > 1)
            input_error(filename, " has ", candidates.size(), " heading columns (\"", candidates[0], "\" and \"",
                    candidates[1], "\"): keep the one in degrees from true north.");
    }
}

GeodeticColumns find_geodetic_columns(const LabeledCsv& csv, LengthUnit unit)
{
    // Aliases are matched after LabeledCsv::normalize, which lower-cases,
    // drops a parenthesized unit and removes punctuation -- so "Alt (ft)",
    // "alt_ft" and "ALT" are all the one name, and a spelling here needs no
    // punctuation variants of its own.
    GeodeticColumns columns;
    columns.time = csv.column({"t", "time", "elapsed", "elapsed time", "seconds", "sec"});
    columns.utc = csv.column({"utcsec", "utc", "utc seconds", "utc time"});
    columns.latitude = csv.column({"lat", "latitude", "lat deg", "latitude deg"});
    columns.longitude = csv.column({"lon", "long", "longitude", "lon deg", "longitude deg"});
    // Only the spellings that agree with the unit actually being read.  A file
    // that says "alt_ft" while --units says meters is not matched, and failing
    // to parse is the right answer: matching it would read one unit as the
    // other and quietly fly the aircraft three times too high or too low.  The
    // bare spellings carry no claim and are always accepted.
    columns.altitude = csv.column({"alt", "altitude", "alt msl", "altitude msl", "height"});
    if (columns.altitude == LabeledCsv::npos)
    {
        columns.altitude = unit.to_meters == 1.0
                ? csv.column({"alt m", "altitude m", "alt msl m"})
                : csv.column({"alt ft", "altitude ft", "alt msl ft"});
    }
    columns.east_speed = csv.column({"ve", "v east", "veast", "east velocity"});
    columns.north_speed = csv.column({"vn", "v north", "vnorth", "north velocity"});
    columns.up_speed = csv.column({"vu", "vup", "v up", "up velocity", "climb rate", "climb"});
    columns.heading = csv.column(heading_aliases);
    return columns;
}

// Fit one coordinate against time.  within_tolerance bounds the *summed*
// squared residual, so an RMS bound in meters is scaled by the sample count to
// mean the same thing whatever the sample rate; zero interpolates, which is
// what a list of waypoints wants.
//
// A dense track (under half a second between samples, the ones that carry an
// aircraft's own motion) gets its ends freed from the natural condition: see
// SmoothingSpline::within_tolerance_extended.  Natural ends pin the acceleration
// at both ends of the track to zero and bias it for several seconds inward, which
// on an approach that ends in a stop is the part nearest the observers who hear it
// loudest.  Sparse waypoints are a statement of intent and are left as they were,
// and so is an exact interpolation (a zero budget has no bias to remove).
constexpr double end_window_seconds = 3.0;
constexpr double end_extension_seconds = 10.0;
constexpr double dense_track_spacing = 0.5; // seconds

SmoothingSpline fit_path(const std::vector<double>& t, const std::vector<double>& f, double rms_meters)
{
    const double spacing = (t.back() - t.front()) / static_cast<double>(t.size() - 1);
    if (rms_meters > 0.0 && t.size() >= 3 && spacing <= dense_track_spacing)
    {
        return SmoothingSpline::within_tolerance_extended(
                t, f, rms_meters, end_window_seconds, end_extension_seconds);
    }
    return SmoothingSpline::within_tolerance(
            t, f, rms_meters * rms_meters * static_cast<double>(t.size()));
}

// Read a labeled WGS84 trajectory into local Cartesian states.  The file gives
// position against time; the velocity and the load factor that sets CT come
// from differentiating the fitted path, unless the file carries velocity
// columns of its own.
void read_geodetic_states(const std::string& filename, const LabeledCsv& csv,
        const GeodeticColumns& columns, const dimensional_data& vehicle,
        const GeodeticInput& input, LengthUnit unit, std::optional<LocalFrame>& frame,
        std::vector<State>& states)
{
    const std::size_t rows = csv.rows.size();
    if (rows < 3)
    {
        input_error("Geodetic trajectory ", filename, " holds ", rows,
                " rows; three is the fewest a path and its curvature can be fitted to.");
    }
    if (!(input.resample_dt > 0.0))
        input_error("Resampling step must be positive, got ", input.resample_dt, " s.");

    const bool have_utc = columns.utc != LabeledCsv::npos;
    const bool from_utc = columns.time == LabeledCsv::npos;
    if (from_utc && !have_utc)
    {
        input_error("Geodetic trajectory ", filename, " names no time column. Add \"time\" "
                "(seconds from the start of the track) or \"utcsec\" (seconds since UTC "
                "midnight).");
    }
    check_declared_units(filename, csv, columns, unit);
    const std::size_t time_column = from_utc ? columns.utc : columns.time;
    const std::size_t needed = columns.width();
    const bool velocity_given = columns.velocity_given();

    std::vector<double> time;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;
    std::vector<double> heading;
    time.reserve(rows);
    x.reserve(rows);
    y.reserve(rows);
    z.reserve(rows);

    double first_time = 0.0;
    for (std::size_t row = 0; row < rows; row++)
    {
        const auto& data = csv.rows[row];
        if (data.size() < needed)
        {
            input_error("Data row ", row + 1, " of ", filename, " has ", data.size(),
                    " fields, too few for the columns its header names.");
        }
        const double latitude = data[columns.latitude];
        const double longitude = data[columns.longitude];
        if (!(std::fabs(latitude) <= 90.0) || !(std::fabs(longitude) <= 180.0))
        {
            input_error("Data row ", row + 1, " of ", filename,
                    " is not a WGS84 position: latitude ", latitude, ", longitude ", longitude,
                    ". Latitude and longitude are read in degrees.");
        }
        if (!frame)
        {
            // Nothing settled the frame for us, so this lone trajectory centers
            // one on its own first position, x east and y north.  A run with
            // several trajectories settles it up front instead -- they share
            // the observers, so they have to share the frame.
            frame = LocalFrame(Geodetic {latitude, longitude, 0.0});
            info("Projecting ", filename, " about its first position, ", latitude, ", ",
                    longitude, ", with x east and y north.");
        }

        // lat and lon were range-checked above, which a NaN fails; the rest of
        // what is read is not, and one NaN then ran through the spline fit into
        // every footprint cell.
        for (std::size_t index : {time_column, columns.altitude, columns.east_speed,
                     columns.north_speed, columns.up_speed, columns.heading})
        {
            if (index != LabeledCsv::npos && !std::isfinite(data[index]))
            {
                input_error("Data row ", row + 1, " of ", filename, " holds a value that is not a "
                        "number (", csv.header[index], ").");
            }
        }

        const double raw_time = data[time_column];
        if (row == 0)
            first_time = raw_time;
        time.push_back(from_utc ? raw_time - first_time : raw_time);

        // Altitude is above mean sea level while z is height above the plane
        // the observers stand on, so the ground's elevation comes off first.
        const double height
                = data[columns.altitude] * unit.to_meters - input.ground_elevation;
        const Vec3 local = frame->to_local(Geodetic {latitude, longitude, height});
        x.push_back(local(0));
        y.push_back(local(1));
        z.push_back(local(2));

        if (columns.heading != LabeledCsv::npos)
        {
            // Unwrapped as it is read, so that interpolating between two rows
            // either side of north turns through north rather than the long
            // way round through south.
            double h = data[columns.heading];
            if (!heading.empty())
                h = heading.back() + std::remainder(h - heading.back(), 360.0);
            heading.push_back(h);
        }
        if (velocity_given)
        {
            // East/north/up in ft/s, turned onto the frame's axes.  A velocity
            // is a direction, so only the rotation applies to it.
            const Vec3 enu(data[columns.east_speed] * unit.to_meters,
                    data[columns.north_speed] * unit.to_meters,
                    data[columns.up_speed] * unit.to_meters);
            const Vec3 velocity = frame->from_east_north_up(enu);
            vx.push_back(velocity(0));
            vy.push_back(velocity(1));
            vz.push_back(velocity(2));
        }
    }

    for (std::size_t row = 1; row < rows; row++)
    {
        if (!(time[row] > time[row - 1]))
        {
            input_error("Time does not increase between data rows ", row, " and ", row + 1,
                    " of ", filename, ": ", time[row - 1], " then ", time[row], ".");
        }
    }

    const double span = time.back() - time.front();
    const double mean_step = span / static_cast<double>(rows - 1);
    double worst = 0.0;
    for (std::size_t row = 1; row < rows; row++)
        worst = std::max(worst, std::fabs((time[row] - time[row - 1]) - mean_step));
    const bool uniform = worst <= 1e-6 * mean_step;

    // Never coarsen.  The sound exposure integral is a Riemann sum over
    // whatever step the states end up with, so throwing away samples a file
    // already has would cost accuracy for nothing.
    const double step = std::min(input.resample_dt, mean_step);
    std::vector<double> sample_times;
    if (uniform && mean_step <= input.resample_dt * (1.0 + 1e-9))
    {
        sample_times = time;
    }
    else
    {
        // Counted by index rather than accumulated, so the spacing cannot
        // drift over a long track.
        const auto count = static_cast<std::size_t>(std::floor(span / step)) + 1;
        sample_times.resize(count);
        for (std::size_t i = 0; i < count; i++)
            sample_times[i] = time.front() + step * static_cast<double>(i);
        info("Resampling ", filename, " from ", rows, " ",
                uniform ? "evenly" : "unevenly", " spaced samples onto a ", step,
                " s grid of ", count, " states.");
    }

    if (!velocity_given && input.path_smoothing <= 0.0 && mean_step < 0.5)
    {
        // Differentiating an interpolant twice amplifies whatever resolution the
        // coordinates were written at by 1/step^2, and the load factor is that
        // second derivative.  The shipped tracks quantize latitude at 1e-8
        // degrees, about a millimeter, which at their 0.02 s spacing moves the
        // footprint by 11.6 dB at the median against the same flight's measured
        // accelerations.  The 10 ft default takes that to 0.085 dB; asking for
        // zero is asking for the 11.6.
        log("Geodetic trajectory ", filename, " is sampled every ", mean_step,
                " s and --path_smoothing 0 was given, so it is flown through exactly. The "
                "load factor is then the curvature of a path through every given point, "
                "coordinate rounding included, which on a track this dense is mostly noise: "
                "expect the footprint to be badly wrong. Drop the flag to use the 10 ft "
                "default.");
    }
    if (!velocity_given && mean_step > 1.0)
    {
        log("Geodetic trajectory ", filename, " is sampled every ", mean_step,
                " s on average. Velocity and load factor are the first and second derivatives "
                "of the fitted path at that spacing, so CT in particular reflects the "
                "interpolant rather than anything the file measured.");
    }
    const double integration_step = sample_times.size() > 1
            ? (sample_times[1] - sample_times[0])
            : 0.0;
    if (integration_step > 0.25)
    {
        log("States come out ", integration_step, " s apart, which is a coarse step for the "
                "sound exposure integral. Pass a smaller -s to resample more finely.");
    }

    const auto fit_x = fit_path(time, x, input.path_smoothing);
    const auto fit_y = fit_path(time, y, input.path_smoothing);
    const auto fit_z = fit_path(time, z, input.path_smoothing);
    // Empty, or the three velocity components; a vector rather than three
    // optionals because either all of them are fitted or none is.
    std::vector<SmoothingSpline> fit_velocity;
    if (velocity_given)
    {
        // Velocity gets its own tolerance, not path_smoothing.  That parameter
        // is a distance, and passing it here silently reused the same number as
        // a *speed* -- one value meaning both meters and meters per second.
        // Harmless while the default was zero; at the current default it would
        // hand every velocity column a ~3 m/s (6 kt) smoothing budget nobody
        // asked for.  Zero interpolates the given velocities, which is what
        // this path did by default before and what a file that bothers to carry
        // velocity columns means by them.
        const double velocity_tolerance = 0.0;
        fit_velocity.push_back(fit_path(time, vx, velocity_tolerance));
        fit_velocity.push_back(fit_path(time, vy, velocity_tolerance));
        fit_velocity.push_back(fit_path(time, vz, velocity_tolerance));
    }

    // A file with its own utcsec column dates the flight itself, unless an
    // operations file's local_time replaces it; otherwise the run says when
    // t = 0 was, and without that the states have no clock (NaN).
    const double utc_at_start = have_utc && !input.clock_replaces_track ? csv.rows.front()[columns.utc]
                                                                         : input.start_utc;

    states.clear();
    states.reserve(sample_times.size());
    for (const double t : sample_times)
    {
        State now {};
        now.t = t;
        now.utcsec = utc_at_start + (t - time.front());
        now.x = fit_x(t);
        now.y = fit_y(t);
        now.z = fit_z(t);
        if (velocity_given)
        {
            now.vx = fit_velocity[0](t);
            now.vy = fit_velocity[1](t);
            now.vz = fit_velocity[2](t);
            now.nx = fit_velocity[0].derivative(t, 1) * mps2g;
            now.ny = fit_velocity[1].derivative(t, 1) * mps2g;
            now.nz = fit_velocity[2].derivative(t, 1) * mps2g + 1.0;
        }
        else
        {
            now.vx = fit_x.derivative(t, 1);
            now.vy = fit_y.derivative(t, 1);
            now.vz = fit_z.derivative(t, 1);
            // The load factor is the acceleration in g, with gravity added back
            // on the vertical: the same convention adapt() uses when it derives
            // states from a smoothed path.
            now.nx = fit_x.derivative(t, 2) * mps2g;
            now.ny = fit_y.derivative(t, 2) * mps2g;
            now.nz = fit_z.derivative(t, 2) * mps2g + 1.0;
        }
        if (!heading.empty())
        {
            // Linear between rows, never a spline: a heading is a measured
            // angle, and a spline through a sharp turn overshoots it.
            const auto above = std::upper_bound(time.begin(), time.end(), t);
            const std::size_t hi = std::min<std::size_t>(
                    static_cast<std::size_t>(above - time.begin()), time.size() - 1);
            const std::size_t lo = hi == 0 ? 0 : hi - 1;
            const double fraction = hi == lo ? 0.0 : (t - time[lo]) / (time[hi] - time[lo]);
            double h = std::fmod(heading[lo] + fraction * (heading[hi] - heading[lo]), 360.0);
            if (h < 0.0)
                h += 360.0;
            now.heading = h;
        }
        derive_flight_parameters(now, vehicle);
        states.push_back(now);
        trace("Parsed geodetic state: ", states.back());
    }
}

} // namespace

namespace niceops {

std::optional<Geodetic> first_geodetic_position(const std::string& filename, LengthUnit unit)
{
    // Parses the whole file to read one row of it, which is the price of having
    // a single parser rather than a second, subtly different one.  It happens
    // once per trajectory at start-up, against a database load of about a
    // second, so it does not show.
    LabeledCsv csv;
    try
    {
        csv = read_labeled_csv(filename);
    }
    catch (const std::exception&)
    {
        // Unreadable, or not a CSV at all.  Say nothing here and let the real
        // load report it, so the failure is described once and in one place.
        return std::nullopt;
    }
    if (csv.rows.empty() || csv.rows.front().size() == 27)
        return std::nullopt; // Empty, or the positional AC state format.
    const auto columns = find_geodetic_columns(csv, unit);
    if (!columns.positioned() || csv.rows.front().size() < columns.width())
        return std::nullopt;
    return Geodetic {csv.rows.front()[columns.latitude], csv.rows.front()[columns.longitude], 0.0};
}

dimensional_data configure_vehicle(const std::string& filename)
{
    info("Reading vehicle configuration: ", filename);
    std::ifstream in(filename);
    if (!in)
        input_error("Unable to open vehicle configuration file: ", filename);
    auto config = nlohmann::json::parse(in);

    dimensional_data vehicle;
    vehicle.radius = json_to_double(config["Main Rotor"]["radius"]);
    vehicle.tip_speed = json_to_double(config["Main Rotor"]["tip speed"]);
    vehicle.weight = json_to_double(config["Vehicle"]["weight"]) * standard_gravity; // Actually a mass (kg), so convert to Newtons
    vehicle.rho = json_to_double(config["Atmosphere"]["density"]);

    // Atmospheric state, for the optional spectral path.  Only density is
    // required -- it is what CT needs -- so each of these falls back rather
    // than failing, and a file written before any of this existed still loads.
    const auto& air = config["Atmosphere"];
    if (air.contains("temperature"))
    {
        vehicle.atmosphere.temperature = json_to_double(air["temperature"]);
        vehicle.temperature_given = true;
    }
    if (air.contains("pressure"))
    {
        vehicle.atmosphere.pressure = json_to_double(air["pressure"]);
        vehicle.pressure_given = true;
    }
    else
    {
        // No pressure given: take the one implied by the density and
        // temperature already in the file, p = rho R T, rather than assuming a
        // standard day the file may not be describing.  S-76D.json pairs
        // 1.225 kg/m^3 with 300.71 K, which is not an ISA point and gives
        // 105.7 kPa, so assuming 101.325 would contradict its own numbers.
        constexpr double dry_air_gas_constant = 287.05; // J/(kg K)
        vehicle.atmosphere.pressure
                = vehicle.rho * dry_air_gas_constant * vehicle.atmosphere.temperature / 1000.0;
    }
    if (air.contains("relative humidity"))
    {
        vehicle.atmosphere.relative_humidity = json_to_double(air["relative humidity"]);
        vehicle.humidity_given = true;
    }

    // What the file gave has to be air.  A temperature of 0 or a pressure of 0
    // made every level NaN, and Celsius, Pa or a humidity fraction shifted the
    // footprint several dB, all with exit status 0.  (The defaults and the
    // pressure derived from density are checked too: a density in the wrong unit
    // shows up there.)
    const auto require_air = [&](const std::string& problem) {
        if (!problem.empty())
            input_error("Vehicle file ", filename, ", Atmosphere section: ", problem);
    };
    require_air(temperature_problem(vehicle.atmosphere.temperature));
    require_air(pressure_problem(vehicle.atmosphere.pressure));
    require_air(humidity_problem(vehicle.atmosphere.relative_humidity));

    trace("Vehicle data structure: ", vehicle);
    return vehicle;
}

namespace {

// The two flight path angles, from the velocity the rotor sees (air_velocity:
// the ground velocity until relative_to_air() sets a wind).  Needs nx/ny/nz
// already normalized.
void derive_path_angles(State& now)
{
    const auto [ux, uy, uz] = air_velocity(now);
    const double v = std::hypot(ux, uy, uz);
    if (v > 0.0)
    {
        auto vn_hat = (ux * now.nx + uy * now.ny + uz * now.nz) / v;
        now.fpa_e = 180.0 / pi * std::asin(vn_hat);
        now.fpa = 180.0 / pi * std::asin(uz / v);
    }
    else
    {
        // True hover: both flight path angles are undefined rather than
        // infinite, and dividing by v would put a NaN into mu's sphere lookup
        // and from there into every level.  Zero is the D original's answer
        // (niceops commit 0b05d3a, "Fixed issue giving NaN flight path angles
        // in hover"), and it is the right one -- a stationary aircraft is
        // neither climbing nor descending.
        now.fpa_e = 0.0;
        now.fpa = 0.0;
    }
}

} // namespace

// Derive the non-dimensional flight parameters from a state's raw kinematics.
// The trajectory constructor and adapt() both need exactly this, from different
// sources, and previously each carried its own copy.  Expects nx/ny/nz to hold
// the un-normalized load-factor components on entry.
void derive_flight_parameters(State& now, const dimensional_data& vehicle)
{
    // Compute load factor
    now.LF = std::hypot(now.nx, now.ny, now.nz);
    // Normalize
    now.nx /= now.LF;
    now.ny /= now.LF;
    now.nz /= now.LF;
    // Non-dimensionalize
    now.CT = now.LF * vehicle.weight / (vehicle.rho * vehicle.tip_speed * vehicle.tip_speed * pi
            * vehicle.radius * vehicle.radius);
    const auto [ux, uy, uz] = air_velocity(now);
    now.mu = std::hypot(ux, uy, uz) / vehicle.tip_speed;
    derive_path_angles(now);
}

// The AC state format does not agree with itself about the sign of vz.  The
// 2019 tracks (tracks/B205/190337AC.csv) store it positive up, as z is; the
// 2017 Noise Abatement tracks store it positive DOWN while z stays positive up
// (vz against dz/dt correlates at -0.996 on B407 run 286270) -- except the 49
// files of EC130B4 test day 298, which are positive up.  Read as up, every
// 2017 descent became a climb of the same angle, selecting climb spheres and
// losing the descent's blade-vortex interaction.  az is positive up in both,
// so only vz is affected.
//
// Decide from the track itself: correlate vz with the derivative of z.  A track
// with too little vertical motion to tell is left alone -- there the sign
// barely matters, and guessing would be worse than not.  Returns true if the
// sign was corrected.
bool correct_vertical_speed_sign(std::vector<State>& states, const std::string& filename)
{
    if (states.size() < 3)
        return false;
    // Uncentered -- the cosine between vz and dz/dt, not their correlation -- so
    // a steady descent, whose dz/dt is constant and has no variance to
    // correlate, still decides: vz . dz/dt is +|dz/dt|^2 when the two agree.
    double s_vz_dz = 0.0, s_vz2 = 0.0, s_dz2 = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 1; i + 1 < states.size(); i++)
    {
        const double dt = states[i + 1].t - states[i - 1].t;
        if (!(dt > 0.0))
            continue;
        const double dz = (states[i + 1].z - states[i - 1].z) / dt;
        const double vz = states[i].vz;
        s_vz_dz += vz * dz;
        s_vz2 += vz * vz;
        s_dz2 += dz * dz;
        n++;
    }
    if (n < 2)
        return false;
    // Below ~0.15 m/s rms of vertical speed there is nothing to decide from.
    if (s_vz2 / n < 0.0225 || s_dz2 / n < 0.0225)
        return false;
    const double agreement = s_vz_dz / std::sqrt(s_vz2 * s_dz2);
    if (agreement > -0.5)
        return false;
    warn("Trajectory ", filename, ": vz runs opposite to dz/dt (cosine ", agreement,
            "), so this file stores vertical speed positive down; negating it.");
    for (auto& now : states)
    {
        now.vz = -now.vz;
        // Only the path angles depend on vz's sign.  The speed, and so mu, do
        // not, and nx/ny/nz were already normalized when the state was derived.
        derive_path_angles(now);
    }
    return true;
}

namespace {

// A static operation: a file whose header names a duration, one row per hold.
// The aircraft stays at one place for each row's duration, with its nose on the
// row's compass heading, and the rows follow one another.  Every hold has to sit
// at the same point: a move is a trajectory, and a jump between holds would be
// neither.  Sampled at 50 Hz, like the AC tracks.  The rotor is at 1 g with zero
// velocity, so mu, both flight path angles and CT come from the same derivation
// as a flown state, and the sphere is read at hover (its nose is the heading).
// A hold has no clock of its own: its UTC seconds are start_utc + t, NaN when
// --start_utc gave none, so the day-night metrics refuse it rather than date it.
void read_static_states(const std::string& filename, const LabeledCsv& csv,
        const dimensional_data& vehicle, LengthUnit unit, double start_utc,
        std::vector<State>& states)
{
    const auto duration_column = csv.column({"duration"});
    const auto x_column = csv.column({"x"});
    const auto y_column = csv.column({"y"});
    const auto z_column = csv.column({"z"});
    const auto heading_column = csv.column({"heading"});
    if (duration_column == LabeledCsv::npos || x_column == LabeledCsv::npos
            || y_column == LabeledCsv::npos || z_column == LabeledCsv::npos
            || heading_column == LabeledCsv::npos)
    {
        input_error(filename, " is a static operation, and needs columns named duration, x, y, z "
                "and heading (x, y and z in --units; heading is compass degrees).");
    }
    constexpr double dt = 0.02;
    const std::size_t needed = std::max({duration_column, x_column, y_column, z_column,
            heading_column}) + 1;
    double t = 0.0;
    double x0 = 0.0, y0 = 0.0, z0 = 0.0;
    for (std::size_t row = 0; row < csv.rows.size(); row++)
    {
        const auto& data = csv.rows[row];
        if (data.size() < needed)
            input_error("Row ", row + 1, " of ", filename, " has ", data.size(),
                    " values, and the static columns need ", needed, ".");
        const double duration = data[duration_column];
        if (!std::isfinite(duration) || !(duration > 0.0))
            input_error("Hold ", row + 1, " of ", filename, " has duration ", duration,
                    " s; it has to be positive.");
        const double x = data[x_column] * unit.to_meters;
        const double y = data[y_column] * unit.to_meters;
        const double z = data[z_column] * unit.to_meters;
        if (row == 0)
        {
            x0 = x;
            y0 = y;
            z0 = z;
        }
        else if (std::fabs(x - x0) > 1e-6 || std::fabs(y - y0) > 1e-6 || std::fabs(z - z0) > 1e-6)
        {
            input_error("Hold ", row + 1, " of ", filename, " is not at the first hold's point. "
                    "Static holds share one position; a move between points is a trajectory.");
        }
        const std::size_t count = std::max<std::size_t>(1,
                static_cast<std::size_t>(std::llround(duration / dt)));
        for (std::size_t k = 0; k < count; k++)
        {
            State now{};
            now.t = t + static_cast<double>(k) * dt;
            now.utcsec = start_utc + now.t;
            now.x = x;
            now.y = y;
            now.z = z;
            now.nz = 1.0; // un-normalized, as derive_flight_parameters expects
            now.heading = data[heading_column];
            derive_flight_parameters(now, vehicle);
            states.push_back(now);
        }
        t += static_cast<double>(count) * dt;
    }
}

} // namespace

Trajectory::Trajectory(const std::string& filename, const std::string& vehicle_filename,
        const GeodeticInput& geodetic, LengthUnit unit)
{
    // Get vehicle data
    vehicle = configure_vehicle(vehicle_filename);

    // D: filename.splitter('.').retro.take(1).array()[0] -- the last '.' segment
    auto dotpos = filename.rfind('.');
    std::string extension = (dotpos == std::string::npos) ? filename : filename.substr(dotpos + 1);
    // Through unsigned char: tolower of a negative char (a non-ASCII byte in the
    // extension) is undefined.
    std::for_each(extension.begin(), extension.end(), [](char& c) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    });
    info("Opening file ", filename, " with extension ", extension);
    if (extension == "csv")
    {
        trace("Loading ", filename, " as CSV.");
        // D reads this with csvReader!(double)(contents, null), which discards
        // the first record as a header regardless of content.  The demo tracks
        // all carry a text header so both paths agreed by luck; a headerless
        // trajectory would have gained a spurious first state here.  The
        // labeled reader keeps that first record instead of discarding it,
        // because the simpler formats below are read by column name -- it
        // drops the same rows, so the AC path is unchanged.
        CsvReport report;
        const auto csv = read_labeled_csv(filename, &report);
        if (report.dropped > 0)
        {
            // A row that did not parse was dropped silently: a time gap in a
            // track, a missing point in the fit.  Text before the first number
            // is a header and was always fine.
            input_error("Line ", report.first_dropped_line, " of ", filename,
                    " is not a row of numbers", report.dropped > 1
                            ? " (and " + std::to_string(report.dropped - 1) + " more are not)"
                            : std::string(),
                    ". A text line is only allowed above the first row of data.");
        }

        // An empty file falls through to the "no usable states" check below,
        // which names the format it expected; front() on no rows would not.
        if (csv.rows.empty())
        {
            input_error("No usable rows in ", filename,
                    ": expected the 27-column NASA/FAA AC state format, or a WGS84 track whose "
                    "header names latitude, longitude, altitude and time columns.");
        }

        // 27 columns is the positional NASA/FAA AC state format; anything else
        // is read by the names in its header.  The AC format's own header
        // happens to name lat/lon/alt columns, so the width has to be tested
        // first or a state file would be read as a WGS84 track -- and its x/y
        // origin, which it does not record, is not its first position.
        if (csv.rows.front().size() == 27)
        {
            std::size_t data_row = 0;
            for (const auto& data : csv.rows)
            {
                data_row++;
                if (data.size() == 27)
                {
                    // One NaN or infinity in a column the states are built from
                    // used to run through the whole acoustics and come out as
                    // NaN footprint cells, with exit status 0.
                    static const struct
                    {
                        std::size_t column;
                        const char* name;
                    } used[] = {{1, "time"}, {5, "x"}, {6, "y"}, {7, "z"}, {11, "vx"}, {13, "vy"},
                            {15, "vz"}, {21, "ax"}, {22, "ay"}, {23, "az"}};
                    for (const auto& field : used)
                    {
                        if (!std::isfinite(data[field.column]))
                        {
                            input_error("Data row ", data_row, " of ", filename, " holds a value "
                                    "that is not a number (column ", field.column + 1, ", ",
                                    field.name, ").");
                        }
                    }
                    // Assume NASA/FAA AC state data format
                    State now {};
                    // Convert to metric
                    now.utcsec = data[0];
                    now.t = data[1];
                    now.x = data[5] * unit.to_meters;
                    now.y = data[6] * unit.to_meters;
                    now.z = data[7] * unit.to_meters;
                    now.vx = data[11] * unit.to_meters;
                    now.vy = data[13] * unit.to_meters;
                    now.vz = data[15] * unit.to_meters;
                    now.heading = data[8];
                    // Convert to "g"
                    // Accelerations carry the file's length unit too.  Chosen
                    // as a single constant rather than unit.to_meters * mps2g,
                    // so the feet path multiplies by exactly the same fps2g it
                    // always did -- the product of two constants is not the
                    // same double as two multiplies, and the state golden is
                    // pinned at 1e-9.
                    const double accel_to_g = unit.to_meters == 1.0 ? mps2g : fps2g;
                    now.nx = data[21] * accel_to_g;
                    now.ny = data[22] * accel_to_g;
                    now.nz = (data[23] * accel_to_g + 1.0);
                    if (std::isfinite(data[2]) && std::isfinite(data[3]) && std::isfinite(data[4]))
                    {
                        recorded.push_back({data[2], data[3], data[4] * unit.to_meters, now.x, now.y,
                                now.z});
                    }
                    derive_flight_parameters(now, vehicle);
                    states.push_back(now);
                    trace("Parsed state: ", states.back());
                }
                else
                {
                    input_error("Wrong number of columns in ", filename, ".");
                }
            }
            // The geodetic reader refuses a clock that does not advance; this one
            // did not, and a constant or descending time column gave a zero or
            // negative integration step and an all -inf or all NaN footprint.
            for (std::size_t i = 1; i < states.size(); i++)
            {
                if (!(states[i].t > states[i - 1].t))
                {
                    input_error("Time does not increase between data rows ", i, " and ", i + 1,
                            " of ", filename, ": ", states[i - 1].t, " then ", states[i].t, ".");
                }
            }
            correct_vertical_speed_sign(states, filename);
            // An operations file's local_time dates the flight in place of
            // the file's own UTC column.
            if (geodetic.clock_replaces_track)
            {
                const double t0 = states.front().t;
                for (auto& state : states)
                    state.utcsec = geodetic.start_utc + (state.t - t0);
            }
        }
        else if (csv.column({"duration"}) != LabeledCsv::npos)
        {
            read_static_states(filename, csv, vehicle, unit, geodetic.start_utc, states);
        }
        else
        {
            const auto columns = find_geodetic_columns(csv, unit);
            if (!columns.positioned())
            {
                input_error("Unrecognized trajectory format in ", filename,
                        ": it is neither the 27-column NASA/FAA AC state format nor a WGS84 "
                        "track. A WGS84 track needs a header naming latitude and longitude "
                        "columns in degrees, an altitude column in ", unit.name, ", and a time "
                        "or utcsec "
                        "column in seconds.");
            }
            frame_ = geodetic.frame;
            read_geodetic_states(filename, csv, columns, vehicle, geodetic, unit, frame_, states);
        }
    }
    else if (extension == "nc")
    {
        input_error("netCDF trajectories not yet supported: ", filename);
    }
    else
    {
        input_error("File extension ", extension, " not recognized: ", filename);
    }

    // Fail where the cause is legible.  An empty trajectory used to surface as
    // a segfault deep in the acoustics, when the level history turned out to
    // have no peak to integrate around.
    if (states.empty())
    {
        input_error("No usable states parsed from ", filename,
                "; expected the 27-column NASA/FAA AC state format with a header row.");
    }
    info("Parsed ", states.size(), " trajectory states.");
}

void Trajectory::write_state(const std::string& state_file) const
{
    info("Writing state data to ", state_file);
    // Create the directory first, as every other writer here does.  state/ is
    // not tracked, so on a fresh checkout `-S state/ci_state.csv` had nowhere
    // to write and fatal()'d -- which aborted with the message only in the log
    // file, before warnings and errors were echoed to stderr as well.
    ensure_parent_directory(state_file);
    std::ofstream file(state_file);
    if (!file)
        input_error("Unable to open state file for writing: ", state_file);
    // Write header from field names (mirrors __traits(allMembers, State))
    file << "t,x,y,z,vx,vy,vz,nx,ny,nz,LF,CT,mu,fpa_e,fpa\n";
    file << std::setprecision(17);
    for (const auto& state : states)
    {
        file << state.t << "," << state.x << "," << state.y << "," << state.z << "," << state.vx << ","
            << state.vy << "," << state.vz << "," << state.nx << "," << state.ny << "," << state.nz << ","
            << state.LF << "," << state.CT << "," << state.mu << "," << state.fpa_e << "," << state.fpa
            << "\n";
    }
}

void Trajectory::adapt(double dt, double smoothing, double tolerance)
{
    if (tolerance >= 0.0)
        info("Adapting trajectory to dt ", dt, ", smoothing to within a residual of ", tolerance);
    else
        info("Adapting trajectory to dt ", dt, " with smoothing parameter ", smoothing);
    assert(dt > 0.0);

    // Unpack states
    std::vector<double> t;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    for (const auto& now : states)
    {
        t.push_back(now.t);
        x.push_back(now.x);
        y.push_back(now.y);
        z.push_back(now.z);
    }
    // The resampled instants, counted by index rather than accumulated so the
    // spacing cannot drift over a long track.
    std::vector<double> ti;
    const double span = t.back() - t.front();
    const auto steps = static_cast<std::size_t>(std::floor(span / dt));
    ti.reserve(steps + 1);
    for (std::size_t i = 0; i <= steps; i++)
        ti.push_back(t.front() + dt * static_cast<double>(i));

    trace("Smoothing x, y and z.");
    const auto fit = [&t, smoothing, tolerance](const std::vector<double>& f) {
        // A negative tolerance means "use the smoothing parameter"; otherwise
        // the penalty is solved for so the fit stays within `tolerance`.
        return tolerance >= 0.0 ? SmoothingSpline::within_tolerance(t, f, tolerance)
                                : SmoothingSpline::with_smoothing(t, f, smoothing);
    };
    const auto fit_x = fit(x);
    const auto fit_y = fit(y);
    const auto fit_z = fit(z);

    // UTC seconds ride along with the resampled time.  They used to be dropped:
    // the repacked State was value-initialized and nothing refilled utcsec, so
    // adapting a trajectory silently moved the whole flight to UTC midnight and
    // with it the day/night classification behind --dnlfile.
    const double utc_at_start = states.front().utcsec;
    const double t_at_start = t.front();

    // Repack data into states
    states.clear();
    trace("States reset to ", states);
    states.reserve(ti.size());
    for (const double now_t : ti)
    {
        State now {};
        now.t = now_t;
        now.utcsec = utc_at_start + (now_t - t_at_start);
        now.x = fit_x(now_t);
        now.y = fit_y(now_t);
        now.z = fit_z(now_t);
        // Velocity and acceleration come from the spline itself rather than
        // from finite differences over the resampled values.  Fitting a curve
        // and then differencing its samples throws away the exact derivatives
        // the fit already carries, and charges the resampling step's truncation
        // error for the privilege -- the load factor, which is the second
        // derivative, paid for that twice.  read_geodetic_states derives its
        // states the same way.
        now.vx = fit_x.derivative(now_t, 1);
        now.vy = fit_y.derivative(now_t, 1);
        now.vz = fit_z.derivative(now_t, 1);
        // Convert to "g", with gravity added back on the vertical
        now.nx = fit_x.derivative(now_t, 2) * mps2g;
        now.ny = fit_y.derivative(now_t, 2) * mps2g;
        now.nz = fit_z.derivative(now_t, 2) * mps2g + 1.0;
        derive_flight_parameters(now, vehicle);
        states.push_back(now);
        trace("Adapted state: ", states.back());
    }
}

std::vector<double> Trajectory::column(const std::string& varname) const
{
    static const std::map<std::string, double State::*> fields = {{"t", &State::t}, {"x", &State::x},
        {"y", &State::y}, {"z", &State::z}, {"vx", &State::vx}, {"vy", &State::vy}, {"vz", &State::vz},
        {"nx", &State::nx}, {"ny", &State::ny}, {"nz", &State::nz}, {"LF", &State::LF}, {"CT", &State::CT},
        {"mu", &State::mu}, {"fpa_e", &State::fpa_e}, {"fpa", &State::fpa}};
    auto it = fields.find(varname);
    if (it == fields.end())
        input_error("Unknown state variable: ", varname);
    std::vector<double> x;
    x.reserve(states.size());
    for (const auto& state : states)
    {
        x.push_back(state.*(it->second));
    }
    return x;
}

double Trajectory::relative_to_air(const std::function<std::array<double, 2>(double)>& wind)
{
    double strongest = 0.0;
    for (auto& now : states)
    {
        const auto w = wind(now.z);
        now.wx = w[0];
        now.wy = w[1];
        strongest = std::max(strongest, std::hypot(w[0], w[1]));
        // nx/ny/nz, LF and CT stay as the constructor derived them: normalized,
        // and inertial.
        const auto [ux, uy, uz] = air_velocity(now);
        now.mu = std::hypot(ux, uy, uz) / vehicle.tip_speed;
        derive_path_angles(now);
    }
    return strongest;
}

void Trajectory::watts_trim(double xstart_meters)
{
    auto mu = column("mu");
    auto z = column("z");
    auto imin_z = std::min_element(z.begin(), z.end()) - z.begin();
    auto imin_v = std::min_element(mu.begin(), mu.end()) - mu.begin();
    auto ixmax = imin_z < imin_v ? imin_z : imin_v;
    auto x = column("x");
    std::size_t ixmin = 0;
    for (std::size_t i = 0; i < x.size(); i++)
    {
        if (x[i] > xstart_meters)
            break;
        ixmin = i;
    }
    // The trim runs from the last state at or before `xstart` to the track's
    // lowest or slowest point.  A start past that point built a vector from a
    // reversed iterator range, which threw "vector" with no hint of the cause on
    // libc++, and a start on it left a one-state trajectory that has no sample
    // interval to integrate over.
    const auto last = static_cast<std::size_t>(ixmax);
    if (ixmin >= last)
    {
        input_error("--trim ", xstart_meters, " m leaves ",
                ixmin > last ? std::string("nothing") : std::string("one state"),
                " of the trajectory: it is cropped from the last state at or before the trim "
                "start to the track's lowest or slowest point, which is at x = ", x[ixmax],
                " m (the track runs from x = ", x.front(), " m to ", x.back(),
                " m). Start the trim earlier along the track.");
    }
    info("Trimming trajectory from x = ", x[ixmin], " to ", x[ixmax]);
    std::vector<State> sub(states.begin() + ixmin, states.begin() + ixmax + 1);
    states.swap(sub);
}

std::optional<NoseSource> parse_nose_source(const std::string& name)
{
    if (name == "velocity")
        return NoseSource::velocity;
    if (name == "ground_velocity")
        return NoseSource::ground_velocity;
    if (name == "heading")
        return NoseSource::heading;
    return std::nullopt;
}

const char* to_string(NoseSource source)
{
    switch (source)
    {
    case NoseSource::ground_velocity:
        return "ground_velocity";
    case NoseSource::heading:
        return "heading";
    case NoseSource::velocity:
        break;
    }
    return "velocity";
}

double nose_velocity_weight(double speed)
{
    if (speed <= heading_nose_speed)
        return 0.0;
    if (speed >= velocity_nose_speed)
        return 1.0;
    const double s = (speed - heading_nose_speed) / (velocity_nose_speed - heading_nose_speed);
    return s * s * (3.0 - 2.0 * s);
}

std::optional<double> heading_offset(const std::vector<State>& states, double frame_bearing,
        bool ground)
{
    const bool has_heading = std::any_of(states.begin(), states.end(),
            [](const State& s) { return !std::isnan(s.heading); });
    if (!has_heading)
        return std::nullopt;
    // A known frame wins.  The fit below absorbs the track's mean crab, which on
    // a crabbed approach ending in hover turned the hover nose 2-7 degrees off
    // the heading the file gives.
    if (!std::isnan(frame_bearing))
        return frame_bearing * pi / 180.0;
    constexpr double fast = 20.0 * knot; // m/s
    constexpr std::size_t enough = 50;
    double sum_cos = 0.0, sum_sin = 0.0;
    std::size_t count = 0;
    for (const auto& s : states)
    {
        if (std::isnan(s.heading))
            continue;
        const auto v = ground ? ground_velocity(s) : air_velocity(s);
        if (std::hypot(v[0], v[1]) < fast)
            continue;
        const double offset = std::atan2(v[1], v[0]) + s.heading * pi / 180.0;
        sum_cos += std::cos(offset);
        sum_sin += std::sin(offset);
        count++;
    }
    if (count >= enough)
        return std::atan2(sum_sin, sum_cos);
    return std::nullopt;
}

std::array<double, 3> nose_direction(const State& state, std::optional<double> offset,
        NoseSource source)
{
    const auto v = source == NoseSource::ground_velocity ? ground_velocity(state)
                                                         : air_velocity(state);
    if (!offset || std::isnan(state.heading))
        return v;
    const double speed = std::hypot(v[0], v[1]);
    const double heading_angle = *offset - state.heading * pi / 180.0;
    if (source == NoseSource::heading)
    {
        // The heading's direction with the velocity's climb: parallel to v
        // itself when the heading is along the track.  The weight takes the
        // climb out at hover, where vz over a drift speed means nothing.
        const double w = nose_velocity_weight(speed);
        return {std::cos(heading_angle), std::sin(heading_angle), w > 0.0 ? w * v[2] / speed : 0.0};
    }
    if (speed >= velocity_nose_speed)
        return v;
    if (speed <= heading_nose_speed)
        return {std::cos(heading_angle), std::sin(heading_angle), 0.0};
    // The horizontal angle turns from the heading to the velocity along the
    // short way round, by a smoothstep weight on the speed, so the nose moves
    // continuously through the band.  The vertical part of the velocity is
    // carried by the same weight, scaled to the horizontal unit length.
    const double w = nose_velocity_weight(speed);
    const double turn = std::remainder(std::atan2(v[1], v[0]) - heading_angle, 2.0 * pi);
    const double angle = heading_angle + w * turn;
    return {std::cos(angle), std::sin(angle), w * v[2] / speed};
}

} // namespace niceops
