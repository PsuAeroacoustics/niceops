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
#include "util/csv.hpp"
#include "util/logger.hpp"
#include "util/paths.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

constexpr double PI = 3.14159265358979323846;

void diff(double dx, const std::vector<double>& f, std::vector<double>& df, std::vector<double>& ddf);

void smoother(const std::vector<double>& x, const std::vector<double>& f,
        const std::vector<double>& xi, std::vector<double>& fi, std::vector<double>& dfi,
        std::vector<double>& ddfi, double smoothing, double tolerance)
{
    // A negative tolerance means "use the smoothing parameter"; otherwise the
    // penalty is solved for so the fit stays within `tolerance` of the data.
    const auto spline = tolerance >= 0.0 ? SmoothingSpline::within_tolerance(x, f, tolerance)
                                         : SmoothingSpline::with_smoothing(x, f, smoothing);
    // Use spline for interpolation
    fi = spline(xi);
    // Use finite differences to get velocity and acceleration
    const double dx = xi[1] - xi[0];
    diff(dx, fi, dfi, ddfi);
}

void diff(double dx, const std::vector<double>& f, std::vector<double>& df, std::vector<double>& ddf)
{
    // Create storage for first derivative
    df = f;
    // Forward diff first point
    df[0] = (f[1] - f[0]) / dx;
    // Central diff interior points
    for (std::size_t i = 1; i + 1 < f.size(); i++)
        df[i] = 0.5 * (f[i + 1] - f[i - 1]) / dx;
    // Backward diff last point
    df[f.size() - 1] = (f[f.size() - 1] - f[f.size() - 2]) / dx;
    // Create storage for second derivative
    ddf = df;
    // Forward diff first point
    ddf[0] = (f[2] - 2 * f[1] + f[0]) / (dx * dx);
    // Central diff interior points
    for (std::size_t i = 1; i + 1 < f.size(); i++)
        ddf[i] = (f[i + 1] - 2 * f[i] + f[i - 1]) / (dx * dx);
    // Backward diff last point
    ddf[f.size() - 1] = (f[f.size() - 1] - 2 * f[f.size() - 2] + f[f.size() - 3]) / (dx * dx);
}

// The vehicle JSONs store every value as a quoted string, and the D original
// read them with `.str` followed by `to!double`.  nlohmann's get<double>()
// throws on a string, so convert explicitly (tolerating a real JSON number).
double json_to_double(const nlohmann::json& value)
{
    if (value.is_string())
        return std::stod(value.get<std::string>());
    return value.get<double>();
}

// Where the labelled WGS84 format's columns sit in a particular file.  That
// format is read by column name rather than by position, so whatever a planner
// already exports can be used without being reordered first, and a file carries
// only the columns it actually has.
struct GeodeticColumns
{
    std::size_t time = LabelledCsv::npos;
    std::size_t utc = LabelledCsv::npos;
    std::size_t latitude = LabelledCsv::npos;
    std::size_t longitude = LabelledCsv::npos;
    std::size_t altitude = LabelledCsv::npos;
    std::size_t east_speed = LabelledCsv::npos;
    std::size_t north_speed = LabelledCsv::npos;
    std::size_t up_speed = LabelledCsv::npos;

    // A file is a WGS84 track when it says where the aircraft is.  Time can
    // come from either of two columns and velocity is optional, so neither of
    // those decides it.
    bool positioned() const
    {
        return latitude != LabelledCsv::npos && longitude != LabelledCsv::npos
                && altitude != LabelledCsv::npos;
    }

    // Velocity is taken from the file only when all three components are there;
    // two of them plus a differentiated third would mix two different estimates
    // of the same vector.
    bool velocity_given() const
    {
        return east_speed != LabelledCsv::npos && north_speed != LabelledCsv::npos
                && up_speed != LabelledCsv::npos;
    }

    // The widest index this file will be read at, so a short row is caught
    // where it can be named rather than as an out-of-range read.
    std::size_t width() const
    {
        std::size_t widest = 0;
        for (std::size_t index : {time, utc, latitude, longitude, altitude, east_speed,
                     north_speed, up_speed})
        {
            if (index != LabelledCsv::npos)
                widest = std::max(widest, index + 1);
        }
        return widest;
    }
};

GeodeticColumns find_geodetic_columns(const LabelledCsv& csv)
{
    // Aliases are matched after LabelledCsv::normalise, which lower-cases,
    // drops a parenthesised unit and removes punctuation -- so "Alt (ft)",
    // "alt_ft" and "ALT" are all the one name, and a spelling here needs no
    // punctuation variants of its own.
    GeodeticColumns columns;
    columns.time = csv.column({"t", "time", "elapsed", "elapsed time", "seconds", "sec"});
    columns.utc = csv.column({"utcsec", "utc", "utc seconds", "utc time"});
    columns.latitude = csv.column({"lat", "latitude", "lat deg", "latitude deg"});
    columns.longitude = csv.column({"lon", "long", "longitude", "lon deg", "longitude deg"});
    // Only the spellings that agree with the unit actually read.  "alt_ft" is
    // matched and "alt_m" deliberately is not: an unmatched altitude column
    // makes the file fail to parse at all, which is the right answer, where
    // matching it would read metres as feet and quietly fly the aircraft three
    // times too high.
    columns.altitude = csv.column({"alt", "altitude", "alt msl", "altitude msl", "height",
            "alt ft", "altitude ft", "alt msl ft"});
    columns.east_speed = csv.column({"ve", "v east", "veast", "east velocity"});
    columns.north_speed = csv.column({"vn", "v north", "vnorth", "north velocity"});
    columns.up_speed = csv.column({"vu", "vup", "v up", "up velocity", "climb rate", "climb"});
    return columns;
}

// Fit one coordinate against time.  within_tolerance bounds the *summed*
// squared residual, so an RMS bound in metres is scaled by the sample count to
// mean the same thing whatever the sample rate; zero interpolates, which is
// what a list of waypoints wants.
SmoothingSpline fit_path(const std::vector<double>& t, const std::vector<double>& f, double rms_metres)
{
    return SmoothingSpline::within_tolerance(
            t, f, rms_metres * rms_metres * static_cast<double>(t.size()));
}

// Read a labelled WGS84 trajectory into local Cartesian states.  The file gives
// position against time; the velocity and the load factor that sets CT come
// from differentiating the fitted path, unless the file carries velocity
// columns of its own.
void read_geodetic_states(const std::string& filename, const LabelledCsv& csv,
        const GeodeticColumns& columns, const dimensional_data& vehicle,
        const GeodeticInput& input, std::optional<LocalFrame>& frame,
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

    const bool have_utc = columns.utc != LabelledCsv::npos;
    const bool from_utc = columns.time == LabelledCsv::npos;
    if (from_utc && !have_utc)
    {
        input_error("Geodetic trajectory ", filename, " names no time column. Add \"time\" "
                "(seconds from the start of the track) or \"utcsec\" (seconds since UTC "
                "midnight).");
    }
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
            // Nothing settled the frame for us, so this lone trajectory centres
            // one on its own first position, x east and y north.  A run with
            // several trajectories settles it up front instead -- they share
            // the observers, so they have to share the frame.
            frame = LocalFrame(Geodetic {latitude, longitude, 0.0});
            info("Projecting ", filename, " about its first position, ", latitude, ", ",
                    longitude, ", with x east and y north.");
        }

        const double raw_time = data[time_column];
        if (row == 0)
            first_time = raw_time;
        time.push_back(from_utc ? raw_time - first_time : raw_time);

        // Altitude is above mean sea level while z is height above the plane
        // the observers stand on, so the ground's elevation comes off first.
        const double height = (data[columns.altitude] - input.ground_elevation_ft) * ft2m;
        const Vec3 local = frame->to_local(Geodetic {latitude, longitude, height});
        x.push_back(local(0));
        y.push_back(local(1));
        z.push_back(local(2));

        if (velocity_given)
        {
            // East/north/up in ft/s, turned onto the frame's axes.  A velocity
            // is a direction, so only the rotation applies to it.
            const Vec3 enu(data[columns.east_speed] * ft2m, data[columns.north_speed] * ft2m,
                    data[columns.up_speed] * ft2m);
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
        // second derivative.  The shipped tracks quantise latitude at 1e-8
        // degrees, about a millimetre, which at their 0.02 s spacing is worth
        // a few tenths of a g -- enough to move CT by 17% and the footprint by
        // 3 dB in the worst cell.  A path smoothing of 0.02 m takes that to
        // 0.03 dB median and 2.4 dB worst against the same flight's measured
        // accelerations.
        log("Geodetic trajectory ", filename, " is sampled every ", mean_step,
                " s and is being flown through exactly (--path_smoothing 0). The load factor "
                "is then the curvature of a path through every given point, coordinate "
                "rounding included. For measured positions rather than planned waypoints, "
                "pass --path_smoothing (0.02 m suits the shipped tracks' resolution).");
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
        fit_velocity.push_back(fit_path(time, vx, input.path_smoothing));
        fit_velocity.push_back(fit_path(time, vy, input.path_smoothing));
        fit_velocity.push_back(fit_path(time, vz, input.path_smoothing));
    }

    // A file with its own utcsec column dates the flight itself; otherwise the
    // run says when t = 0 was, and the default of zero puts it at UTC midnight.
    const double utc_at_start = have_utc ? csv.rows.front()[columns.utc] : input.start_utc;

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
        derive_flight_parameters(now, vehicle);
        states.push_back(now);
        trace("Parsed geodetic state: ", states.back());
    }
}

} // namespace

namespace niceops {

std::optional<Geodetic> first_geodetic_position(const std::string& filename)
{
    // Parses the whole file to read one row of it, which is the price of having
    // a single parser rather than a second, subtly different one.  It happens
    // once per trajectory at start-up, against a database load of about a
    // second, so it does not show.
    LabelledCsv csv;
    try
    {
        csv = read_labelled_csv(filename);
    }
    catch (const std::exception&)
    {
        // Unreadable, or not a CSV at all.  Say nothing here and let the real
        // load report it, so the failure is described once and in one place.
        return std::nullopt;
    }
    if (csv.rows.empty() || csv.rows.front().size() == 27)
        return std::nullopt; // Empty, or the positional AC state format.
    const auto columns = find_geodetic_columns(csv);
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
    vehicle.weight = json_to_double(config["Vehicle"]["weight"]) * g; // Actually a mass (kg), so convert to Newtons
    vehicle.rho = json_to_double(config["Atmosphere"]["density"]);

    // Atmospheric state, for the optional spectral path.  Only density is
    // required -- it is what CT needs -- so each of these falls back rather
    // than failing, and a file written before any of this existed still loads.
    const auto& air = config["Atmosphere"];
    if (air.contains("temperature"))
        vehicle.atmosphere.temperature = json_to_double(air["temperature"]);
    if (air.contains("pressure"))
    {
        vehicle.atmosphere.pressure = json_to_double(air["pressure"]);
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
        vehicle.atmosphere.relative_humidity = json_to_double(air["relative humidity"]);

    trace("Vehicle data structure: ", vehicle);
    return vehicle;
}

// Derive the non-dimensional flight parameters from a state's raw kinematics.
// The trajectory constructor and adapt() both need exactly this, from different
// sources, and previously each carried its own copy.  Expects nx/ny/nz to hold
// the un-normalised load-factor components on entry.
void derive_flight_parameters(State& now, const dimensional_data& vehicle)
{
    // Compute load factor
    now.LF = std::hypot(now.nx, now.ny, now.nz);
    // Normalize
    now.nx /= now.LF;
    now.ny /= now.LF;
    now.nz /= now.LF;
    // Non-dimensionalize
    now.CT = now.LF * vehicle.weight / (vehicle.rho * vehicle.tip_speed * vehicle.tip_speed * PI
            * vehicle.radius * vehicle.radius);
    auto v = std::hypot(now.vx, now.vy, now.vz);
    now.mu = v / vehicle.tip_speed;
    if (v > 0.0)
    {
        auto vn_hat = (now.vx * now.nx + now.vy * now.ny + now.vz * now.nz) / v;
        now.fpa_e = 180.0 / PI * std::asin(vn_hat);
        now.fpa = 180.0 / PI * std::asin(now.vz / v);
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

} // namespace niceops

Trajectory::Trajectory(const std::string& filename, const std::string& vehicle_filename,
        const GeodeticInput& geodetic)
{
    // Get vehicle data
    vehicle = configure_vehicle(vehicle_filename);

    // D: filename.splitter('.').retro.take(1).array()[0] -- the last '.' segment
    auto dotpos = filename.rfind('.');
    std::string extension = (dotpos == std::string::npos) ? filename : filename.substr(dotpos + 1);
    std::for_each(extension.begin(), extension.end(), [](char& c) { c = std::tolower(c); });
    info("Opening file ", filename, " with extension ", extension);
    if (extension == "csv")
    {
        trace("Loading ", filename, " as CSV.");
        // D reads this with csvReader!(double)(contents, null), which discards
        // the first record as a header regardless of content.  The demo tracks
        // all carry a text header so both paths agreed by luck; a headerless
        // trajectory would have gained a spurious first state here.  The
        // labelled reader keeps that first record instead of discarding it,
        // because the simpler formats below are read by column name -- it
        // drops the same rows, so the AC path is unchanged.
        const auto csv = read_labelled_csv(filename);

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
            for (const auto& data : csv.rows)
            {
                if (data.size() == 27)
                {
                    // Assume NASA/FAA AC state data format
                    State now {};
                    // Convert to metric
                    now.utcsec = data[0];
                    now.t = data[1];
                    now.x = data[5] * ft2m;
                    now.y = data[6] * ft2m;
                    now.z = data[7] * ft2m;
                    now.vx = data[11] * ft2m;
                    now.vy = data[13] * ft2m;
                    now.vz = data[15] * ft2m;
                    // Convert to "g"
                    now.nx = data[21] * fps2g;
                    now.ny = data[22] * fps2g;
                    now.nz = (data[23] * fps2g + 1.0);
                    derive_flight_parameters(now, vehicle);
                    states.push_back(now);
                    trace("Parsed state: ", states.back());
                }
                else
                {
                    input_error("Wrong number of columns in ", filename, ".");
                }
            }
        }
        else
        {
            const auto columns = find_geodetic_columns(csv);
            if (!columns.positioned())
            {
                input_error("Unrecognised trajectory format in ", filename,
                        ": it is neither the 27-column NASA/FAA AC state format nor a WGS84 "
                        "track. A WGS84 track needs a header naming latitude and longitude "
                        "columns in degrees, an altitude column in feet, and a time or utcsec "
                        "column in seconds.");
            }
            frame_ = geodetic.frame;
            read_geodetic_states(filename, csv, columns, vehicle, geodetic, frame_, states);
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
    // to write and fatal()'d -- which aborts silently, since the message goes
    // to the log file and not to stderr.
    ensure_parent_directory(state_file);
    std::ofstream file(state_file);
    if (!file)
        fatal("Unable to open state file for writing: ", state_file);
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
    // Build interpolants
    std::vector<double> ti;
    ti.push_back(t[0]);
    double tnext = ti[0] + dt;
    while (tnext <= t[t.size() - 1])
    {
        ti.push_back(tnext);
        tnext += dt;
    }
    trace("Smoothing x.");
    std::vector<double> xi;
    std::vector<double> vxi;
    std::vector<double> axi;
    smoother(t, x, ti, xi, vxi, axi, smoothing, tolerance);

    trace("Smoothing y.");
    std::vector<double> yi;
    std::vector<double> vyi;
    std::vector<double> ayi;
    smoother(t, y, ti, yi, vyi, ayi, smoothing, tolerance);

    trace("Smoothing z.");
    std::vector<double> zi;
    std::vector<double> vzi;
    std::vector<double> azi;
    smoother(t, z, ti, zi, vzi, azi, smoothing, tolerance);

    // Repack data into states
    states.clear();
    trace("States reset to ", states);
    for (std::size_t i = 0; i < ti.size(); i++)
    {
        State now {};
        now.t = ti[i];
        now.x = xi[i];
        now.y = yi[i];
        now.z = zi[i];
        now.vx = vxi[i];
        now.vy = vyi[i];
        now.vz = vzi[i];
        // Convert to "g"
        now.nx = axi[i] * mps2g;
        now.ny = ayi[i] * mps2g;
        now.nz = azi[i] * mps2g + 1.0;
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
        fatal("Unknown state variable: ", varname);
    std::vector<double> x;
    x.reserve(states.size());
    for (const auto& state : states)
    {
        x.push_back(state.*(it->second));
    }
    return x;
}

void Trajectory::watts_trim(double xstart_ft)
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
        if (x[i] > xstart_ft * ft2m)
            break;
        ixmin = i;
    }
    info("Trimming trajectory from x = ", x[ixmin], " to ", x[ixmax]);
    std::vector<State> sub(states.begin() + ixmin, states.begin() + ixmax + 1);
    states.swap(sub);
}
