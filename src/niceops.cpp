// Port of niceops.d.  Pipeline: load Database + Trajectory -> for each
// observer and each trajectory state, compute instantaneous dBA via
// Database::level_at_observer -> integrate the time history into a sound
// exposure level (SEL) using the 10-dB-down-from-peak window -> write CSV
// output and optional matplotlib figures.

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>
#include "vec3.hpp"

#include "plotting.hpp"
#include <complex>

#include "atmosphere.hpp"
#include "ground_effect.hpp"
#include "weighting.hpp"
#include "sphere_database.hpp"
#include "trajectory.hpp"
#include "util/csv.hpp"
#include "util/logger.hpp"
#include "util/paths.hpp"
#include "util/parallel.hpp"
#include "util/sel.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

// ft2m comes from trajectory.hpp (niceops::ft2m); defining a second one
// here made every use ambiguous under the using-directive above.

// D's iota(begin, end, step): evenly spaced values from `begin`, EXCLUSIVE of
// `end`.  Counted by index rather than accumulated, so the spacing cannot drift.
std::vector<double> iota_exclusive(double begin, double end, double step)
{
    const auto count = static_cast<std::size_t>(std::ceil((end - begin) / step));
    std::vector<double> values(count);
    for (std::size_t i = 0; i < count; i++)
        values[i] = begin + step * static_cast<double>(i);
    return values;
}

// Open a file for writing, creating its parent directory first.  All three
// writers below need exactly this, and Logger does the same for its log file.
std::ofstream open_output(const std::string& path, std::ios::openmode mode = std::ios::out)
{
    ensure_parent_directory(path);
    std::ofstream file(path, mode);
    if (!file)
        input_error("Unable to open file for writing: ", path);
    file << std::setprecision(17);
    return file;
}

void write_levels(const std::vector<double>& x_obs, const std::vector<double>& y_obs,
        const std::vector<double>& levels, const std::string& out_file, bool scattered)
{
    info("Writing output to ", out_file);
    auto file = open_output(out_file);

    if (scattered)
    {
        info("Scattered output.");
        assert(levels.size() == x_obs.size() && levels.size() == y_obs.size());
        file << "x(ft), y(ft), Level(dB)\n";
        for (std::size_t i = 0; i < levels.size(); i++)
        {
            file << x_obs[i] / ft2m << "," << y_obs[i] / ft2m << "," << levels[i] << "\n";
        }
    }
    else
    {
        info("Gridded output.");
        const std::size_t nx = x_obs.size();
        const std::size_t ny = y_obs.size();
        assert(nx * ny == levels.size());

        // Write header
        file << "Line 1: x(ft), Line 2: y(ft),  Lines 3+ levels (x cols y rows)\n";
        // Comma-separated, newline at the end of each row.
        auto write_row = [&file](const std::vector<double>& values, std::size_t first,
                                 std::size_t count, double scale) {
            for (std::size_t i = 0; i < count; i++)
                file << (i ? "," : "") << values[first + i] / scale;
            file << "\n";
        };
        write_row(x_obs, 0, nx, ft2m);
        write_row(y_obs, 0, ny, ft2m);
        for (std::size_t row = 0; row * nx < levels.size(); row++)
            write_row(levels, row * nx, nx, 1.0);
        // The D original writes one extra newline at the end; keep it.
        file << "\n";
    }
}

// Observer positions back as WGS84, so a footprint can be put on a map.  One
// row per observer whatever the observer set is: a grid that is regular in the
// local frame is not regular in latitude and longitude, so there is no gridded
// form of this file.
void write_geodetic_levels(const std::vector<double>& x_obs, const std::vector<double>& y_obs,
        const std::vector<double>& levels, const LocalFrame& frame, const std::string& out_file,
        bool scattered)
{
    info("Writing geodetic output to ", out_file);
    auto file = open_output(out_file);
    file << "latitude(deg), longitude(deg), Level(dB)\n";
    // Observers sit on the ground plane, so the altitude that comes back is the
    // ground's own and the level is what an observer standing there hears.
    auto write_observer = [&file, &frame](double x, double y, double level) {
        const auto position = frame.to_geodetic(Vec3(x, y, 0.0));
        file << position.latitude << "," << position.longitude << "," << level << "\n";
    };
    if (scattered)
    {
        assert(levels.size() == x_obs.size() && levels.size() == y_obs.size());
        for (std::size_t i = 0; i < levels.size(); i++)
            write_observer(x_obs[i], y_obs[i], levels[i]);
    }
    else
    {
        const std::size_t nx = x_obs.size();
        const std::size_t ny = y_obs.size();
        assert(nx * ny == levels.size());
        for (std::size_t iy = 0; iy < ny; iy++)
            for (std::size_t ix = 0; ix < nx; ix++)
                write_observer(x_obs[ix], y_obs[iy], levels[ix + iy * nx]);
    }
}

// Where a level file's WGS84 companion goes: alongside it, under the same name
// with "_latlon" on the end.  Derived rather than asked for, so --latlon needs
// no second output name per trajectory and cannot fall out of step with -o.
std::string latlon_companion(const std::string& level_file)
{
    const std::filesystem::path path(level_file);
    const std::string extension = path.extension().string();
    return (path.parent_path()
            / (path.stem().string() + "_latlon" + (extension.empty() ? ".csv" : extension)))
            .string();
}

void write_time_history_file(const std::vector<double>& time_history, double dt, double x, double y,
        const std::string& filename)
{
    info("Writing time history to ", filename);
    auto file = open_output(filename);
    file << "Level time history. Next row [dt (s), x (ft), y(ft)], remaining rows, dB.\n";
    file << dt << "," << x / ft2m << "," << y / ft2m << "\n";
    for (auto level : time_history)
    {
        file << level << "\n";
    }
}

void append_metrics_file(const std::string& filename, const std::string& casename,
        const std::vector<double>& levels, double threshold = 65.0)
{
    auto file = open_output(filename, std::ios::app);
    double average = 0;
    double maximum = -std::numeric_limits<double>::infinity();
    std::size_t count = 0;
    // Iterate over range and compute average and maximum value
    for (auto level : levels)
    {
        average += level;
        maximum = std::max(maximum, level);
        if (level >= threshold)
        {
            count++;
        }
    }
    average /= levels.size();
    auto area = double(count) / levels.size();
    file << casename << "," << average << "," << maximum << "," << area << "\n";
}

// Sound exposure level at one observer without marching every emission time.
//
// The integral only ever uses the window from the peak down to 10 dB below it,
// which on a shipped case is a little over half the track -- the rest is
// computed and thrown away.  peak_range's window is *contiguous by
// construction*: it walks out from the peak and stops at the first sample below
// the threshold.  So nothing outside the window needs a bound at all; the walk
// stops there by definition.  A bound is needed for one thing only -- to be
// sure the walk started on the right peak -- and DirectionBounds is tight
// enough (median slack 8.6 dB, inside the 10 dB window) that proving that costs
// almost nothing beyond the window itself.
//
// The result is identical to evaluating every state, not an approximation: a
// skipped state is provably neither the peak nor inside the window.
double windowed_sound_exposure_level(const Database& database,
        const std::vector<Database::SourceState>& states, const DirectionBounds& bounds,
        const Vec3& observer, double step, std::vector<double>& levels,
        std::vector<unsigned>& visited, unsigned generation, std::size_t& evaluated)
{
    const std::size_t nt = states.size();
    const double sphere_radius_squared = database.sphere_radius() * database.sphere_radius();
    // A coarse sweep seeds the peak search; the stride is not sensitive -- 16
    // and 1024 land within a couple of percent of each other.
    const std::size_t seed_stride = std::max<std::size_t>(1, nt / 128);

    // The displacement's components in the sphere frame.  The bound test needs
    // exactly what an evaluation needs, which is why they are formed once and
    // handed on rather than recomputed -- doing it twice was measured at 1.48x
    // down to 1.08x.
    struct Geometry
    {
        double squared_range, along, up, across;
    };
    auto geometry_of = [&](std::size_t it) {
        const Vec3 delta = observer - states[it].position;
        return Geometry {delta.squaredNorm(), states[it].frame.nose.dot(delta),
                states[it].frame.normal.dot(delta), states[it].frame.starboard.dot(delta)};
    };
    auto evaluate_with = [&](std::size_t it, const Geometry& g) {
        visited[it] = generation;
        evaluated++;
        return levels[it] = database.level_from_geometry(g.squared_range, g.along, g.up, g.across,
                       states[it].blend);
    };
    auto evaluate = [&](std::size_t it) {
        if (visited[it] == generation)
            return levels[it];
        return evaluate_with(it, geometry_of(it));
    };

    // Seed with a coarse sweep so the bound has something to beat from the
    // first state onwards.  The level history has one dominant peak, so any
    // coarse sample near it is enough.
    double best = -std::numeric_limits<double>::infinity();
    std::size_t peak = 0;
    for (std::size_t it = 0; it < nt; it += seed_stride)
    {
        const double level = evaluate(it);
        if (level > best)
        {
            best = level;
            peak = it;
        }
    }
    double best_energy = std::pow(10.0, 0.1 * best);

    // One pass: form the geometry, test the bound against the best level so far,
    // and evaluate whatever it cannot rule out -- handing that same geometry
    // straight to the evaluation, so a state that survives pays nothing for the
    // test.  `best` only rises, so a state skipped early stays correctly
    // skipped, and no second pass is needed.
    for (std::size_t it = 0; it < nt; it++)
    {
        if (visited[it] == generation)
            continue;
        const Geometry g = geometry_of(it);
        // The bound subtracts spherical spreading, which is only a loss outside
        // the sphere; inside it, evaluate rather than trust it.
        if (g.squared_range > sphere_radius_squared
                && bounds.numerator(it,
                           DirectionBounds::bin(g.along, g.up, g.across, g.squared_range))
                        < best_energy * g.squared_range)
            continue;
        const double level = evaluate_with(it, g);
        if (level > best)
        {
            best = level;
            peak = it;
            best_energy = std::pow(10.0, 0.1 * best);
        }
    }

    // `best` is now provably the peak of the whole history, so the window and
    // the integral are the ordinary ones -- expand_window() and
    // sel_over_window() from sel.hpp, reading through `evaluate` instead of a
    // filled array.  Sharing them is what keeps this path's window identical to
    // peak_range's, tie rule included, rather than merely intended to be.
    const auto index = expand_window(peak, nt, best - sel_window_db, evaluate);
    return sel_over_window(index[0], index[1], step, evaluate);
}

void read_observer_file(const std::string& observer_file, std::vector<Vec3>& observers)
{
    trace("Reading observers from ", observer_file, ".");
    // D reads this with csvReader!(double)(contents, null), which consumes the
    // first record as a header and discards it whether or not it is numeric.
    auto records = read_numeric_csv(observer_file, /*skip_first_record=*/true);
    for (const auto& data : records)
    {
        if (data.size() == 3)
        {
            observers.push_back(Vec3(data[0] * ft2m, data[1] * ft2m, data[2] * ft2m));
        }
        else if (data.size() == 4)
        {
            // Skip channel ID in first column
            observers.push_back(Vec3(data[1] * ft2m, data[2] * ft2m, data[3] * ft2m));
        }
        else
        {
            input_error("Observer file format unreadable.");
        }
    }
    if (observers.empty())
    {
        // Every row failing to parse leaves this empty, and the run then
        // "succeeds" having evaluated nothing and written a lone header line.
        // The first record is consumed as a header whether or not it holds
        // numbers, so a one-line file lands here too.
        input_error("Observer file ", observer_file, " yielded no observers. The first record "
                "is taken as a header and discarded; every row after it must be x,y,z or "
                "id,x,y,z in feet.");
    }
}

} // namespace

namespace {

// Everything about a run that does not vary between trajectories.
struct RunOptions
{
    bool qsam = true;
    double dt = 0.5;
    double ambient_level = 0.0;
    unsigned thread_count = 0;
    double db_cutoff = 0.0;
    // Spectral propagation: carry the database's per-band spectra to the ground
    // and integrate a metric from them, rather than propagating one broadband
    // dBA with the database's fitted EAA factor.
    bool spectral = false;
    BandLimits band_limits;
    // Excess ground attenuation, off unless asked for: it is the single most
    // expensive thing in the evaluation.
    bool ground_effect = false;
    GroundModel ground;
    // Overrides for the atmosphere read from the vehicle JSON; NaN means "use
    // whatever the file said".
    double humidity_override = std::numeric_limits<double>::quiet_NaN();
    double temperature_override = std::numeric_limits<double>::quiet_NaN();
    double pressure_override = std::numeric_limits<double>::quiet_NaN();
    double trim_start_ft = 0.0;
    double utc_offset = 0.0;
    bool scattered = false;
    bool want_dnl = false;
    bool write_time_history = false;
    bool write_state = false;
    std::string time_history_prefix;
    std::string metric_file;
    std::string power_file;
    // How WGS84 trajectory input is projected, and the frame it is projected
    // onto.  The frame is settled once for the whole run: the observers are
    // common to every trajectory in it, so their coordinates have to mean the
    // same thing in all of them.
    GeodeticInput geodetic;
    // Write a WGS84 companion beside every level file.
    bool latlon = false;
};

// The files one trajectory is computed from and written to, and where it sits
// in the run.
struct PathJob
{
    std::string trajectory;
    std::string vehicle;
    std::string output;
    std::string state;
    std::size_t count = 1;
};

// Compute one trajectory's footprint and write its outputs, returning the
// trajectory itself for the single-path plotting block.  Anything that throws
// in here condemns this trajectory alone; the caller catches and moves on.
Trajectory process_trajectory(const Database& database, const PathJob& job,
        const RunOptions& options, const std::vector<Vec3>& observers,
        const std::vector<double>& x_obs, const std::vector<double>& y_obs,
        std::vector<double>& integrated_levels, std::vector<double>& dnl_energy)
{
    // Named so the body below reads as it did when it was inline in main.
    const auto& path_file = job.trajectory;
    const auto& vehicle_file = job.vehicle;
    const auto& out_file = job.output;
    const std::size_t path_count = job.count;
    const bool qsam = options.qsam;
    const double dt = options.dt;
    const double ambient_level = options.ambient_level;
    const unsigned thread_count = options.thread_count;
    const double db_cutoff = options.db_cutoff;
    const double trim_start_ft = options.trim_start_ft;
    const double utc_offset = options.utc_offset;
    const bool scattered = options.scattered;
    const bool want_dnl = options.want_dnl;
    const bool write_time_history = options.write_time_history;
    const bool write_state = options.write_state;
    const auto& time_history_prefix = options.time_history_prefix;
    const auto& metric_file = options.metric_file;
    const auto& power_file = options.power_file;


    // Create trajectory from CSV data
    Trajectory trajectory(path_file, vehicle_file, options.geodetic);

    // mu and CT are non-dimensionalised by the vehicle JSON's tip speed, and
    // they index the database's operating conditions -- which were generated at
    // whatever tip speed the database records.  Disagree, and every lookup lands
    // at the wrong condition with nothing downstream able to notice.  Worth
    // checking because the two are easy to state in different units:
    // S-76D_M3.nod was built at 220 m/s, and 600 is that same figure read as
    // ft/s, which puts the tip at Mach 1.7.
    if (!std::isnan(database.tip_speed())
            && !approx_equal(trajectory.vehicle.tip_speed, database.tip_speed()))
    {
        log("Main rotor tip speed disagrees: ", vehicle_file, " gives ",
                trajectory.vehicle.tip_speed, " m/s but the database was built at ",
                database.tip_speed(), " m/s. mu and CT will index the wrong conditions.");
    }

    // Trim trajectory (per Mike): crop from the given x to whichever comes
    // first, the lowest or the slowest point.  Off unless asked for -- it
    // changes which part of the flight is integrated, and so the levels.
    if (!std::isnan(trim_start_ft))
    {
        trajectory.watts_trim(trim_start_ft);
    }

    // Trajectory adaptation parameters, for the adapt() call below.  Both the
    // call and its arguments are disabled, matching the D original.
    //const int degree = 5;
    //const double smooth = 1e6;
    //trajectory.adapt(dt, degree, smooth);

    // And time series data
    const std::size_t nt = trajectory.states.size();

    // The sound exposure level integral weights each sample by the trajectory's
    // own time step, so take that from the states rather than from -s.  Those
    // two coincide only when adapt() resamples onto the -s grid, and adapt() is
    // disabled below; every shipped track is sampled at 0.02 s, so the 0.5 s
    // default would overstate every level by 10*log10(0.5/0.02) = 14 dB.
    double integration_step = dt;
    if (nt >= 2)
    {
        std::vector<double> sample_times;
        sample_times.reserve(nt);
        for (const auto& state : trajectory.states)
            sample_times.push_back(state.t);
        double relative_spread = 0.0;
        integration_step = uniform_sample_interval(sample_times, &relative_spread);
        info("Integrating with the trajectory's sample interval of ", integration_step, " s.");
        if (relative_spread > 1e-6)
        {
            log("Trajectory sampling is not uniform (worst interval deviates by ",
                    100.0 * relative_spread, "% from the mean); the sound exposure level "
                    "integral assumes the mean step.");
        }
    }

    // The source geometry and the blend of database operating conditions depend
    // only on the trajectory state, so build them once here rather than once
    // per (observer, state) pair.  The blend is a CGAL Delaunay locate() walk
    // and was by far the largest single cost in the program; this turns
    // observers*states of them into states of them.
    // Night classification for the DNL penalty, from the AC format's UTC
    // seconds column shifted into local time.
    std::vector<char> night(want_dnl ? nt : 0);
    for (std::size_t it = 0; it < night.size(); it++)
        night[it] = is_night(trajectory.states[it].utcsec, utc_offset) ? 1 : 0;

    std::vector<Database::SourceState> source_states;
    source_states.reserve(nt);
    for (const auto& state : trajectory.states)
    {
        source_states.push_back(database.source_state(Vec3(state.x, state.y, state.z),
                Vec3(state.vx, state.vy, state.vz),
                Vec3(state.nx, state.ny, state.nz), state.mu,
                qsam ? state.fpa_e : state.fpa, state.CT));
    }

    if (!power_file.empty())
    {
        // Shaft power is a property of the flight condition, so it falls out of
        // the per-state blends with no extra geometry.  It is written per
        // trajectory rather than mixed into the state file, whose columns are
        // pinned by the state golden.
        if (!database.has_shaft_power())
        {
            log("--powerfile given but this database carries no ShaftPower; skipping ",
                    path_file, ".");
        }
        else
        {
            const auto name = path_count == 1
                    ? power_file
                    : (std::filesystem::path(power_file).parent_path()
                              / (std::filesystem::path(path_file).stem().string() + "_power.csv"))
                              .string();
            auto file = open_output(name);
            info("Writing shaft power to ", name);
            file << "t(s),advance_ratio,flight_path_angle(deg),CT,shaft_power\n";
            for (std::size_t it = 0; it < nt; it++)
            {
                const auto& state = trajectory.states[it];
                file << state.t << "," << state.mu << "," << (qsam ? state.fpa_e : state.fpa) << ","
                     << state.CT << "," << database.shaft_power(source_states[it].blend) << "\n";
            }
        }
    }
    // ambient_level stays NaN unless -a was given, which is how D's default-NaN
    // double suppressed the ambient mix-in.  The reference pressure it converts
    // to is loop-invariant; it used to be a pow() per evaluation.
    const bool add_ambient = !std::isnan(ambient_level);
    const double ambient_pressure = add_ambient ? std::pow(10.0, ambient_level / 20.) : 0.0;

    // Observers are independent of one another: each reads the database and the
    // source_states array, which are finished and const by now, and writes only
    // its own level.  Nothing in here reaches the CGAL triangulation -- the
    // blends were all resolved above -- so the loop splits with no locking.
    const unsigned workers = thread_count > 0 ? thread_count : default_thread_count();
    info("Evaluating ", observers.size(), " observers on ", workers,
            workers == 1 ? " thread." : " threads.");

    if (write_time_history)
    {
        // Create the output directory once, here.  Left to the writers it would
        // be a create_directories() race between workers, and a stat per file
        // besides.
        const auto directory = std::filesystem::path(time_history_prefix).parent_path();
        if (!directory.empty())
            std::filesystem::create_directories(directory);
    }

    // With --db_cutoff, skip any (observer, state) pair whose level is provably
    // below the cutoff.  The bound is per state, from the loudest sphere the
    // state blends and spherical spreading alone, so it never discards a sample
    // that could have been louder -- unlike the D original, which sorted
    // observers by range and broke out of the loop at the first quiet one, a
    // shortcut that ignores directivity and can drop a distant observer the
    // rotor is pointed straight at.
    //
    // The threshold sits 10 dB below the cutoff, as the D original's did, and
    // that margin is what makes the result exact where it matters: the sound
    // exposure level integrates the window from the peak down to peak - 10, so
    // for any observer whose peak reaches the cutoff, every sample in its window
    // is at or above cutoff - 10 and none of them is ever culled.  Observers
    // that never reach the cutoff are the only ones whose level changes, and
    // they are precisely the ones the cutoff says are not of interest.
    // Spectral propagation setup.  Both of these depend only on the atmosphere
    // and the database's band centres, so they are computed once here rather
    // than per observer: the coefficients cost a pow() and two exp() per band.
    std::vector<double> absorption;
    std::vector<double> band_gains;
    const GroundModel* ground_model = nullptr;
    std::vector<std::complex<double>> ground_impedance_bands;
    double sound_speed = 343.2;
    // The contiguous run of bands inside the metric's limits. Band centres are
    // sorted, so the limits always select a range rather than a scattered set.
    std::size_t first_band = 0;
    std::size_t last_band = 0;
    if (options.spectral)
    {
        if (!database.has_spectra())
        {
            input_error("--spectral needs a database carrying per-band spectra; this one "
                    "has none. Rebuild it with spectra, or drop --spectral to propagate the "
                    "broadband level.");
        }

        Atmosphere air = trajectory.vehicle.atmosphere;
        if (!std::isnan(options.temperature_override))
            air.temperature = options.temperature_override;
        if (!std::isnan(options.pressure_override))
            air.pressure = options.pressure_override;
        if (!std::isnan(options.humidity_override))
            air.relative_humidity = options.humidity_override;

        const auto& bands = database.band_centres();
        const auto retained = bands_in_range(bands, options.band_limits);
        if (retained == 0)
        {
            input_error("--band_limits ", options.band_limits.low, " to ",
                    options.band_limits.high, " Hz excludes every one of the database's ",
                    bands.size(), " bands (", bands.front(), " to ", bands.back(), " Hz).");
        }

        absorption = absorption_coefficients(air, bands);
        band_gains = a_weighting_gains(bands, options.band_limits);
        while (first_band < bands.size() && band_gains[first_band] == 0.0)
            first_band++;
        last_band = bands.size();
        while (last_band > first_band && band_gains[last_band - 1] == 0.0)
            last_band--;

        log("Propagating spectra: ", retained, " of ", bands.size(), " bands within ",
                options.band_limits.low, " to ", options.band_limits.high, " Hz, at ",
                air.temperature, " K, ", air.pressure, " kPa, ", air.relative_humidity,
                "% RH.");
        // Say what that atmosphere costs at the extremes actually in play, so a
        // wrong humidity is visible in the log rather than only in the answer.
        log("  absorption spans ", absorption.front(), " dB/m at ", bands.front(), " Hz to ",
                absorption.back(), " dB/m at ", bands.back(), " Hz.");

        sound_speed = air.sound_speed();
        if (options.ground_effect)
        {
            ground_model = &options.ground;
            ground_impedance_bands = ground_impedance(options.ground, bands);
            log("  ground effect on: flow resistance ", options.ground.flow_resistance,
                    " kPa s/m^2, receiver ", options.ground.receiver_height, " m up, turbulence ",
                    options.ground.turbulence, ", sound speed ", sound_speed, " m/s.");
        }
    }

    const bool culling = !std::isnan(db_cutoff);
    const double cull_below = db_cutoff - sel_window_db;
    // A culled sample contributes no acoustic energy.  Filling with the
    // threshold instead would be far worse than useless -- thousands of samples
    // at the floor integrate to well above it.
    const double culled_level = -std::numeric_limits<double>::infinity();
    std::vector<double> cull_range_squared(culling ? nt : 0);
    for (std::size_t it = 0; it < cull_range_squared.size(); it++)
    {
        const double range = database.cull_range(source_states[it].blend, cull_below);
        cull_range_squared[it] = range * range;
    }
    std::atomic<std::size_t> culled_samples {0};

    // Find each observer's peak and expand out to the 10-dB-down points, rather
    // than marching every emission time -- see windowed_sound_exposure_level.
    // It is exact, so it is on whenever it applies, which is the plain
    // broadband run.  Everything else marches, for one of two reasons:
    //
    //   the level is transformed after eval, so the bound -- which is on what
    //   eval returns -- would have to be transformed with it:
    //     --spectral     the bound is on the broadband dBA, which is not a
    //                    bound on the band-limited weighted level (the same
    //                    reason --db_cutoff is refused with --spectral)
    //     -a / --ambient
    //
    //   every sample is wanted, not just the window's:
    //     -i / --history writes the whole time history, skipped states included
    //     -D / --dnlfile the day-night window is found on a penalized history,
    //                    which a window found on the raw one does not cover
    //
    // A new option that does either belongs in the matching term below.
    // --db_cutoff stands apart: it already culls on its own terms, and its
    // culling is a deliberate metric change, so leave the two mechanisms
    // separate rather than layering this one on top.
    const bool level_is_transformed = options.spectral || add_ambient;
    const bool needs_full_history = write_time_history || want_dnl;
    DirectionBounds direction_bounds;
    if (!level_is_transformed && !needs_full_history && !culling)
    {
        // direction_bounds() reports itself unusable on an empty trajectory, so
        // nt > 0 does not need testing here as well.
        direction_bounds = database.direction_bounds(source_states);
    }
    const bool expand_from_peak = direction_bounds.usable();
    std::atomic<std::size_t> evaluated_samples {0};

    const auto acoustics_start = std::chrono::steady_clock::now();

    parallel_chunks(observers.size(), workers, [&](std::size_t begin, std::size_t end) {
        // One scratch history per chunk rather than per observer; parallel.hpp
        // deals out several chunks per worker, so this runs a few times a
        // thread, which it measures as noise against the work in the chunk.
        std::vector<double> levels(nt);
        std::vector<double> penalized;
        // Which states of `levels` this observer has filled, stamped by a
        // counter so the array does not have to be cleared per observer.
        std::vector<unsigned> visited(expand_from_peak ? nt : 0, 0u);
        unsigned generation = 0;
        std::size_t culled_here = 0;
        std::size_t evaluated_here = 0;
        for (std::size_t index = begin; index < end; index++)
        {
            const auto& observer = observers[index];
            if (expand_from_peak)
            {
                if (++generation == 0)
                {
                    // Wrapped, so old stamps could be mistaken for this pass.
                    std::fill(visited.begin(), visited.end(), 0u);
                    generation = 1;
                }
                integrated_levels[index] = windowed_sound_exposure_level(database, source_states,
                        direction_bounds, observer, integration_step, levels, visited,
                        generation, evaluated_here);
                continue;
            }
            for (std::size_t it = 0; it < nt; it++)
            {
                if (culling
                        && (observer - source_states[it].position).squaredNorm()
                                > cull_range_squared[it])
                {
                    levels[it] = culled_level;
                    culled_here++;
                    continue;
                }
                double level = options.spectral
                        ? database.spectral_level_at_observer(observer, source_states[it],
                                absorption.data(), band_gains, ground_model,
                                &ground_impedance_bands, sound_speed, first_band, last_band)
                        : database.level_at_observer(observer, source_states[it]);
                if (add_ambient)
                {
                    // Add ambient noise
                    level = 20.0 * std::log10(std::pow(10.0, level / 20.) + ambient_pressure);
                }
                levels[it] = level;
            }
            integrated_levels[index] = sound_exposure_level(levels, integration_step);
            if (want_dnl)
            {
                // Energies add across operations; the logarithm comes at the
                // end, once every trajectory has contributed.  Each observer
                // index belongs to exactly one worker, and the trajectories run
                // one after another, so this needs no synchronisation.
                dnl_energy[index] += std::pow(
                        10.0, 0.1 * day_night_level(levels, integration_step, night, penalized));
            }
            if (write_time_history)
            {
                // Write out time history files
                auto time_history_filename = time_history_prefix + "_i" + std::to_string(index) + ".csv";
                write_time_history_file(levels, integration_step, observer(0), observer(1),
                        time_history_filename);
            }
        }
        culled_samples += culled_here;
        evaluated_samples += evaluated_here;
    });

    if (expand_from_peak)
    {
        const auto total = observers.size() * nt;
        log("Evaluated ", evaluated_samples.load(), " of ", total, " (",
                100.0 * static_cast<double>(evaluated_samples.load())
                        / static_cast<double>(total),
                "%); the rest is provably outside every observer's 10 dB window.");
    }

    if (culling)
    {
        const auto total = observers.size() * nt;
        log("Culled ", culled_samples.load(), " of ", total, " evaluations (",
                100.0 * static_cast<double>(culled_samples.load()) / static_cast<double>(total),
                "%) as provably below ", cull_below, " dB.");
    }

    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - acoustics_start)
                                 .count();
    log("Acoustic calculations completed in ", elapsed / 1000.0, " seconds.");
    write_levels(x_obs, y_obs, integrated_levels, out_file, scattered);
    if (options.latlon)
    {
        write_geodetic_levels(x_obs, y_obs, integrated_levels, *options.geodetic.frame,
                latlon_companion(out_file), scattered);
    }

    if (!metric_file.empty())
    {
        auto casename = std::filesystem::path(path_file).stem().string();
        append_metrics_file(metric_file, casename, integrated_levels);
    }

    if (write_state)
    {
        trajectory.write_state(job.state);
    }

    return trajectory;
}

} // namespace

// Everything after the command line is parsed.  Wrapped by main below, which
// owns the decision about what a failure costs: one trajectory, or the run.
int run(int argc, char* argv[])
{
    // Default to demo case
    std::string log_file = "log/niceops.log";
    // Each of these may be given once, to apply to every trajectory, or once per
    // trajectory.  Several trajectories in one run is what makes the day-night
    // level meaningful: DNL sums the acoustic energy of a whole day's
    // operations, not of one flight.
    std::vector<std::string> database_files = {"databases/S-76D_M3.nod"};
    std::vector<std::string> path_files = {"tracks/S-76D/197421AC.csv"};
    std::vector<std::string> vehicle_files = {"databases/S-76D.json"};
    std::vector<std::string> out_files = {"integrated/197421_LAE.csv"};
    std::vector<std::string> state_files;
    std::string dnl_file = "";
    double utc_offset = 0.0;
    double trim_start_ft = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> grid_bounds;
    double grid_step = 25.0;
    std::string sphere_file = "";
    std::string aam_file = "";
    std::string power_file = "";
    std::string log_level = "info";
    std::vector<double> sphere_condition;
    bool write_state = false;
    bool plotting = false;
    std::string figure_file = "";
    std::string state_figure_file = "";
    // 65 dB is the usual reporting threshold, and 100 clears the peak of a
    // helicopter footprint.  The port inherited 45/90 from c74bab3, whose upper
    // level sits *below* the peak of the demo case, leaving the loudest patch
    // unfilled; the figures in figures/ were made with these.
    double contour_min = 65.0;
    double contour_max = 100.0;
    double contour_step = 5.0;
    bool scattered = false;
    std::string observer_file = "";
    double dt = 0.5;
    // NOTE: D default-initializes `double` to NaN, not 0, so the !isNaN guard
    // below is FALSE by default and D never mixes in an ambient level unless
    // -a is given.  Initialising this to 0.0 added a spurious 0 dB ambient
    // floor to every run, which pinned quiet observers near 0 dB.
    double ambient_level = std::numeric_limits<double>::quiet_NaN();
    bool write_time_history = false;
    std::string time_history_prefix = "";
    std::string metric_file = "";
    bool qsam = true;
    bool self_test = false;
    unsigned thread_count = 0;
    double db_cutoff = std::numeric_limits<double>::quiet_NaN();
    bool spectral = false;
    std::vector<double> band_limits;
    bool ground_effect = false;
    double flow_resistance = flow_resistance::grassy_field;
    double receiver_height = 1.2;
    double turbulence = 0.0;
    double humidity_override = std::numeric_limits<double>::quiet_NaN();
    double temperature_override = std::numeric_limits<double>::quiet_NaN();
    double pressure_override = std::numeric_limits<double>::quiet_NaN();
    // WGS84 input and output.  An unset origin is taken from the first geodetic
    // trajectory; an unset bearing leaves the axes east/north, which is what a
    // map wants -- the shipped AC tracks are the other case, their +x running
    // along the flight track.
    std::vector<double> origin;
    double frame_bearing = std::numeric_limits<double>::quiet_NaN();
    double ground_elevation = 0.0;
    double start_utc = 0.0;
    double path_smoothing = 0.0;
    bool latlon = false;

    CLI::App app{"Noise Informed Community Environment Operation Planning System (NICE-OPS)."};
    app.add_option("-d,--database", database_files,
            "netCDF noise database file. Once, or once per trajectory.");
    app.add_option("-l,--logfile", log_file, "Output log file.");
    app.add_option("-L,--loglevel", log_level, "Log verbosity: trace, info, warn, error, fatal.")
            ->check(CLI::IsMember({"trace", "info", "warn", "error", "fatal"}));
    app.add_option("-t,--trajectoryfile", path_files,
            "CSV trajectory input file. Repeat for several flight paths.");
    app.add_option("-v,--vehicle", vehicle_files,
            "JSON vehicle definition file. Once, or once per trajectory.");
    app.add_option("-o,--output", out_files,
            "Output file name. Once per trajectory; derived from the track name if omitted.");
    app.add_flag("-p,--plotting", plotting, "Enable plotting.");
    app.add_option("-f,--figfile", figure_file, "Filename to save figure.");
    app.add_option("-F,--statefigfile", state_figure_file, "Filename to save state figures.");
    app.add_option("--cmin", contour_min, "Minimum contour level, dB.");
    app.add_option("--cmax", contour_max, "Maximum contour level, dB.");
    app.add_option("--cstep", contour_step, "Contour level step, dB.");
    app.add_option("-s,--step", dt, "Resampling time step for --adapt, in seconds. "
            "The sound exposure level integral uses the trajectory's own sample interval.");
    app.add_option("-a,--ambient", ambient_level, "Ambient noise level, dB");
    app.add_option("-i,--history", time_history_prefix, "Time history file prefix.");
    app.add_option("-m,--metricfile", metric_file, "File to append overall noise metrics.");
    app.add_option("-b,--observer", observer_file, "Observer file name.");
    app.add_flag("-q,--qsam,!--no-qsam", qsam, "Extended QSAM (default on; --no-qsam disables)");
    app.add_option("-S,--state", state_files,
            "Filename to save state data. Once per trajectory.");
    app.add_flag("--self_test", self_test, "Run self tests.");
    app.add_option("-j,--threads", thread_count,
            "Worker threads for the observer loop (0 = one per hardware thread).");
    app.add_option("--trim", trim_start_ft,
            "Crop each trajectory from this x (ft) to its lowest/slowest point. Off by default.");
    app.add_option("--grid", grid_bounds,
            "Observer grid bounds in metres: xmin xmax ymin ymax. Default -1500 1500 -1500 1500.")
            ->expected(4);
    app.add_option("--gridstep", grid_step, "Observer grid spacing in metres (default 25).");
    app.add_option("-D,--dnlfile", dnl_file,
            "Write the day-night level, summed over every trajectory, to this file.");
    app.add_option("--utc_offset", utc_offset,
            "Hours from UTC to local time, for the 22:00-07:00 DNL night window.");
    app.add_option("--powerfile", power_file,
            "Write per-state shaft power, if the database carries any, to this file.");
    app.add_option("--export_sphere", sphere_file,
            "Interpolate a whole noise sphere and write it as a one-group .nod database.");
    app.add_option("--export_aam", aam_file,
            "Interpolate a whole noise sphere and write it in AAM/RNM format (needs spectra).");
    app.add_option("--sphere_condition", sphere_condition,
            "Operating condition for --export_sphere: advance ratio, flight path angle, CT.")
            ->expected(3);
    app.add_option("--db_cutoff", db_cutoff,
            "Skip evaluations that provably fall below this level, dB. Off by default.");
    app.add_flag("--spectral", spectral,
            "Propagate the database's per-band spectra to the ground with ISO 9613-1 "
            "absorption, instead of one broadband level with the database's EAA factor. "
            "Requires a database carrying spectra.");
    app.add_option("--band_limits", band_limits,
            "Frequency range integrated into the metric, Hz: low high. Default 20 20000.")
            ->expected(2);
    app.add_option("--humidity", humidity_override,
            "Relative humidity in percent, overriding the vehicle file.");
    app.add_option("--temperature", temperature_override,
            "Temperature in K, overriding the vehicle file.");
    app.add_option("--pressure", pressure_override,
            "Ambient pressure in kPa, overriding the vehicle file.");
    app.add_flag("--ground_effect", ground_effect,
            "Add excess ground attenuation (Chien-Soroka over a Delany-Bazley impedance) to "
            "the spectral path. Needs --spectral.");
    app.add_option("--flow_resistance", flow_resistance,
            "Specific flow resistance of the ground, kPa s/m^2. 30 snow, 225 grass (default), "
            "650 soil, 3000 dirt, 10000 concrete, 50000 asphalt.");
    app.add_option("--receiver_height", receiver_height,
            "Receiver height above ground for --ground_effect, m (default 1.2).");
    app.add_option("--origin", origin,
            "WGS84 reference for the local frame: latitude longitude, in degrees. "
            "Defaults to the first position of the first geodetic trajectory.")
            ->expected(2);
    app.add_option("--frame_bearing", frame_bearing,
            "Compass bearing, degrees, of the local +x axis. Default: x east, y north. "
            "The shipped AC tracks use 210.");
    app.add_option("--ground_elevation", ground_elevation,
            "Elevation of the observer plane above mean sea level, ft. Subtracted from a "
            "geodetic track's altitude column to give height above the observers (default 0).");
    app.add_option("--start_utc", start_utc,
            "Seconds since UTC midnight at t = 0 of a geodetic track, for the day-night "
            "window. Ignored when the file carries its own utcsec column.");
    app.add_option("--path_smoothing", path_smoothing,
            "How far, in metres RMS, the path fitted through a geodetic track may sit from "
            "the given positions. Default 0, which flies through them exactly.");
    app.add_flag("--latlon", latlon,
            "Also write each level file as latitude, longitude, level. Needs --origin unless "
            "the trajectory is itself geodetic.");
    app.add_option("--turbulence", turbulence,
            "Turbulence parameter for --ground_effect, rad s/sqrt(m); 0 to 16e-4. "
            "Zero gives unrealistically deep interference nulls.");

    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        return app.exit(e);
    }

    if ((!sphere_file.empty() || !aam_file.empty()) && sphere_condition.size() != 3)
    {
        input_error("--export_sphere/--export_aam need "
              "--sphere_condition <advance_ratio> <flight_path_angle> <CT>.");
    }

    // Each per-trajectory list is either one entry to share across every
    // trajectory, or one entry per trajectory.
    const std::size_t path_count = path_files.size();
    if (path_count == 0)
        input_error("No trajectory file given.");
    auto expand_per_path = [path_count](std::vector<std::string>& list, const char* what) {
        if (list.size() == path_count)
            return;
        if (list.size() == 1)
        {
            // Copy first: assign() clears the vector before reading its
            // argument, and list.front() points into the vector.
            const std::string shared = list.front();
            list.assign(path_count, shared);
            return;
        }
        input_error("Got ", list.size(), " ", what, " for ", path_count,
                " trajectories; give one, or one per trajectory.");
    };
    if (out_files.size() == 1 && path_count > 1)
    {
        // A single output name cannot serve several tracks, so name each after
        // its own, next to where the given name pointed.
        const auto directory = std::filesystem::path(out_files.front()).parent_path();
        out_files.clear();
        for (const auto& path : path_files)
        {
            const auto stem = std::filesystem::path(path).stem().string();
            out_files.push_back((directory / (stem + "_LAE.csv")).string());
        }
    }
    expand_per_path(database_files, "database files");
    expand_per_path(vehicle_files, "vehicle files");
    expand_per_path(out_files, "output files");
    if (!state_files.empty())
        expand_per_path(state_files, "state files");

    if (grid_bounds.empty())
    {
        grid_bounds = {-1500.0, 1500.0, -1500.0, 1500.0};
    }
    if (!(grid_step > 0.0))
    {
        input_error("--gridstep must be positive, got ", grid_step, ".");
    }
    if (grid_bounds[1] <= grid_bounds[0] || grid_bounds[3] <= grid_bounds[2])
    {
        input_error("--grid needs xmin < xmax and ymin < ymax.");
    }

    if (!observer_file.empty())
    {
        scattered = true;
    }
    if (!time_history_prefix.empty())
    {
        write_time_history = true;
    }
    if (!state_files.empty())
    {
        write_state = true;
    }

    // Set default logging to file
    Logger::instance().set_file(log_file);
    // The eval loop's trace() calls are gated so they cost nothing at info, but
    // they are there to be turned on when something needs explaining.
    const std::map<std::string, LogLevel> levels {{"trace", LogLevel::trace},
            {"info", LogLevel::info}, {"warn", LogLevel::warn}, {"error", LogLevel::error},
            {"fatal", LogLevel::fatal}};
    Logger::instance().set_level(levels.at(log_level));

    // Settle the local frame before anything is loaded.  Every trajectory in a
    // run shares the observers, so a position in one has to mean the same place
    // as the same position in another; one frame for the whole run is what
    // makes that true.
    if (!(path_smoothing >= 0.0))
        input_error("--path_smoothing must not be negative, got ", path_smoothing, ".");
    std::optional<LocalFrame> frame;
    std::optional<Geodetic> reference;
    if (!origin.empty())
    {
        if (!(std::fabs(origin[0]) <= 90.0) || !(std::fabs(origin[1]) <= 180.0))
        {
            input_error("--origin is not a WGS84 position: latitude ", origin[0],
                    ", longitude ", origin[1], ". Both are in degrees.");
        }
        reference = Geodetic {origin[0], origin[1], 0.0};
    }
    else
    {
        // No origin given: the first geodetic trajectory names one.  A run of
        // AC-format tracks finds none, and needs none -- those are already
        // Cartesian.
        for (const auto& path : path_files)
        {
            if (auto first = first_geodetic_position(path))
            {
                reference = *first;
                log("Taking the local origin from the first position of ", path, ": ",
                        reference->latitude, ", ", reference->longitude, ".");
                break;
            }
        }
    }
    if (reference)
    {
        frame = std::isnan(frame_bearing) ? LocalFrame(*reference)
                                          : LocalFrame::track_aligned(*reference, frame_bearing);
    }
    if (latlon && !frame)
    {
        input_error("--latlon has no frame to invert: give --origin <latitude> <longitude> "
                    "(and usually --frame_bearing) saying where this run's x/y origin sits, "
                    "or use a geodetic trajectory, which carries one.");
    }

    // Start stopwatch
    auto watch_start = std::chrono::steady_clock::now();
    // Each Database owns its own triangulation, so several can be live at once;
    // they are still loaded one at a time and reused across consecutive paths
    // that name the same file, which for the usual case of one database shared
    // by every path means loading once.
    std::string loaded_database = database_files.front();
    auto database = Database(loaded_database);
    auto t1 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()
            - watch_start)
                        .count();
    log("Database loaded in ", t1 / 1000.0, " seconds.");

    if (!sphere_file.empty() || !aam_file.empty())
    {
        // A standalone job: interpolate the database down to one operating
        // condition and write that sphere out.  No trajectory is involved.
        const auto sphere = database.interpolate_sphere(
                sphere_condition[0], sphere_condition[1], sphere_condition[2]);
        if (!sphere_file.empty())
            write_sphere(sphere, sphere_file);
        if (!aam_file.empty())
        {
            // The tip speed turns the non-dimensional advance ratio back into
            // the airspeed the AAM format records.  Prefer the one the database
            // itself records, so an exported sphere cannot pick up the tip speed
            // of whatever unrelated vehicle file happened to be on the command
            // line; fall back to the vehicle file for older databases.
            const double tip_speed = !std::isnan(database.tip_speed())
                    ? database.tip_speed()
                    : configure_vehicle(vehicle_files.front()).tip_speed;
            write_aam_sphere(sphere, aam_file, tip_speed);
        }
        return 0;
    }

    if (self_test)
    {
        // Historical note: this used to fail, and a D-era comment here recorded
        // one sample disagreeing by ~1.7 dB.  The cause was legacy hover
        // directivity (see repair_legacy_hover_directivity); with that repaired
        // the self test passes on every shipped database.
        database.self_test(loaded_database);
    }

    std::vector<double> integrated_levels;
    std::vector<double> x_obs;
    std::vector<double> y_obs;
    std::vector<Vec3> observers;
    if (scattered)
    {
        // Read in scattered observer data
        read_observer_file(observer_file, observers);
        // Set length for levels and observer x, y arrays
        integrated_levels.resize(observers.size());
        x_obs.reserve(observers.size());
        y_obs.reserve(observers.size());
        // Deal observers out to x_obs and y_obs arrays
        for (const auto& observer : observers)
        {
            x_obs.push_back(observer(0));
            y_obs.push_back(observer(1));
        }
    }
    else
    {
        // Set up observer array, mirroring D's iota(-1500.0, 1500.0, 25.0).
        // NOTE: D's iota(begin, end, step) is EXCLUSIVE of end -- it yields
        // ceil((end-begin)/step) values, so this is 120 points spanning
        // -1500..1475 m, NOT 121 points ending at 1500.  Verified against LDC.
        x_obs = iota_exclusive(grid_bounds[0], grid_bounds[1], grid_step);
        y_obs = iota_exclusive(grid_bounds[2], grid_bounds[3], grid_step);
        const std::size_t nx = x_obs.size();
        const std::size_t ny = y_obs.size();
        // Set corresponding length for levels
        integrated_levels.resize(nx * ny);
        observers.reserve(nx * ny);
        for (std::size_t iy = 0; iy < ny; iy++)
        {
            for (std::size_t ix = 0; ix < nx; ix++)
            {
                observers.push_back(Vec3(x_obs[ix], y_obs[iy], 0.0));
            }
        }
    }
    // DNL is a property of the whole set of operations, not of one flight: each
    // trajectory's day-night level is an energy, and they sum.
    const bool want_dnl = !dnl_file.empty();
    std::vector<double> dnl_energy(want_dnl ? observers.size() : 0, 0.0);

    // Outlives the loop so the plotting block below, which only runs for a
    // single trajectory, can still reach it.
    std::optional<Trajectory> trajectory_holder;

    RunOptions options;
    options.qsam = qsam;
    options.dt = dt;
    options.ambient_level = ambient_level;
    options.db_cutoff = db_cutoff;
    options.spectral = spectral;
    if (spectral && !std::isnan(db_cutoff))
    {
        // cull_range bounds the broadband level from the loudest sphere in the
        // blend; it says nothing about the A-weighted, band-limited level the
        // spectral path produces, so the two together could silently discard
        // samples that mattered.
        input_error("--spectral cannot be combined with --db_cutoff: the cutoff's range "
                    "bound is derived from the broadband level, not the band-limited one.");
    }
    if (band_limits.size() == 2)
    {
        options.band_limits.low = band_limits[0];
        options.band_limits.high = band_limits[1];
        if (!(options.band_limits.low < options.band_limits.high))
        {
            input_error("--band_limits needs low < high, got ", options.band_limits.low, " and ",
                    options.band_limits.high, ".");
        }
    }
    if (ground_effect && !spectral)
    {
        input_error("--ground_effect needs --spectral: the ground effect is a function of "
                    "frequency, and the broadband path has no spectrum to apply it to.");
    }
    options.ground_effect = ground_effect;
    options.ground.flow_resistance = flow_resistance;
    options.ground.receiver_height = receiver_height;
    options.ground.turbulence = turbulence;
    options.humidity_override = humidity_override;
    options.temperature_override = temperature_override;
    options.pressure_override = pressure_override;
    options.trim_start_ft = trim_start_ft;
    options.utc_offset = utc_offset;
    options.thread_count = thread_count;
    options.scattered = scattered;
    options.want_dnl = want_dnl;
    options.write_time_history = write_time_history;
    options.write_state = write_state;
    options.time_history_prefix = time_history_prefix;
    options.metric_file = metric_file;
    options.power_file = power_file;
    options.geodetic.frame = frame;
    options.geodetic.ground_elevation_ft = ground_elevation;
    options.geodetic.start_utc = start_utc;
    // -s is the resampling step for both adapt() and a geodetic track that
    // needs resampling; neither ever coarsens a track that is already finer.
    options.geodetic.resample_dt = dt;
    options.geodetic.path_smoothing = path_smoothing;
    options.latlon = latlon;

    std::size_t failed_paths = 0;

    for (std::size_t path_index = 0; path_index < path_count; path_index++)
    {
        try
        {
        if (database_files[path_index] != loaded_database)
        {
            log("Loading database ", database_files[path_index], " for ",
                    path_files[path_index]);
            database = Database(database_files[path_index]);
            // Only once it loaded.  Recording the name first would leave the cache
            // claiming a database that failed to open, so a later trajectory naming
            // the same file would skip the reload and quietly use the previous one.
            loaded_database = database_files[path_index];
        }
            PathJob job;
            job.trajectory = path_files[path_index];
            job.vehicle = vehicle_files[path_index];
            job.output = out_files[path_index];
            job.state = state_files.empty() ? std::string() : state_files[path_index];
            job.count = path_count;
            trajectory_holder = process_trajectory(database, job, options, observers,
                    x_obs, y_obs, integrated_levels, dnl_energy);
        }
        catch (const std::exception& problem)
        {
            failed_paths++;
            error("Skipping ", path_files[path_index], ": ", problem.what());
        }
    }

    if (failed_paths > 0)
    {
        log(failed_paths, " of ", path_count, " trajectories could not be processed.");
    }

    if (want_dnl)
    {
        std::vector<double> dnl(dnl_energy.size());
        for (std::size_t i = 0; i < dnl_energy.size(); i++)
            dnl[i] = 10.0 * std::log10(dnl_energy[i]);
        log("Writing day-night level over ", path_count,
                path_count == 1 ? " operation to " : " operations to ", dnl_file);
        write_levels(x_obs, y_obs, dnl, dnl_file, scattered);
        if (latlon)
            write_geodetic_levels(x_obs, y_obs, dnl, *frame, latlon_companion(dnl_file), scattered);
    }

    const bool want_plots = plotting || !figure_file.empty() || !state_figure_file.empty();
    if (!plt::available() && want_plots)
    {
        // Say so up front rather than running the whole case and drawing
        // nothing at the end.
        input_error("This build has no plotting support. Reconfigure with "
              "-DNICEOPS_PLOTTING=ON, or drop -p/-f/-F.");
    }

    if (path_count > 1 && want_plots)
    {
        log("Plotting is skipped for multi-trajectory runs: one figure file cannot hold ",
                path_count, " flight paths.");
    }
    else if (want_plots)
    {
        auto& trajectory = *trajectory_holder;
        auto x = trajectory.column("x");
        auto y = trajectory.column("y");
        auto z = trajectory.column("z");

        // Convert to ft for plotting
        // TODO add units throughout
        for (auto& v : x_obs)
            v *= 1.0 / ft2m;
        for (auto& v : y_obs)
            v *= 1.0 / ft2m;
        for (auto& v : x)
            v *= 1.0 / ft2m;
        for (auto& v : y)
            v *= 1.0 / ft2m;
        for (auto& v : z)
            v *= 1.0 / ft2m;

        if (plotting || !figure_file.empty())
        {
            // Set contour levels, mirroring D's
            // iota(contour_min, contour_max + contour_step, contour_step).
            // iota's end is exclusive, so this runs contour_min..contour_max
            // inclusive -- one level fewer than an inclusive reading gives.
            const auto clevels = iota_exclusive(contour_min, contour_max + contour_step, contour_step);
            plt::figure();
            if (scattered)
            {
                plt::tricontourf(x_obs, y_obs, integrated_levels, clevels, "YlOrRd");
                plt::plot(x_obs, y_obs, "ko");
            }
            else
            {
                std::vector<std::vector<double>> gridded_levels;
                const std::size_t nx = x_obs.size();
                const std::size_t ny = y_obs.size();
                gridded_levels.resize(ny);
                for (std::size_t iy = 0; iy < ny; iy++)
                {
                    gridded_levels[iy].resize(nx);
                    for (std::size_t ix = 0; ix < nx; ix++)
                    {
                        gridded_levels[iy][ix] = integrated_levels[ix + iy * nx];
                    }
                }
                plt::contourf(x_obs, y_obs, gridded_levels, clevels, "YlOrRd");
            }
            // The flight track over the footprint.  This was commented out in the
            // particular niceops.d the port was derived from, and the port
            // copied the comment; it is live both in that file's own ancestor
            // and on master, and every reference figure in figures/ has the
            // track drawn on it.
            plt::plot(x, y, "k-");
            plt::xlabel("X, ft", 14);
            plt::ylabel("Y, ft", 14);
            plt::axis("equal");
            plt::colorbar();
            plt::tight_layout();

            if (!figure_file.empty())
            {
                plt::savefig(figure_file);
            }
        }
        if (plotting || !state_figure_file.empty())
        {
            auto CT = trajectory.column("CT");
            auto mu = trajectory.column("mu");
            auto gamma = trajectory.column("fpa");
            auto gamma_e = trajectory.column("fpa_e");

            auto p = std::filesystem::path(state_figure_file);
            const std::string ext = p.extension().string();
            const std::string base = p.replace_extension().string();

            plt::figure();
            plt::plot(x, CT);
            plt::xlabel("X, ft", 14);
            plt::ylabel("$C_T$", 14);
            if (!state_figure_file.empty())
            {
                plt::savefig(base + "_CT_" + ext);
            }

            plt::figure();
            plt::plot(x, mu);
            plt::xlabel("X, ft", 14);
            plt::ylabel("$\\mu$", 14);
            if (!state_figure_file.empty())
            {
                plt::savefig(base + "_mu_" + ext);
            }

            plt::figure();
            plt::plot(x, gamma, x, gamma_e);
            plt::xlabel("X, ft", 14);
            plt::ylabel("$\\gamma$", 14);
            if (!state_figure_file.empty())
            {
                plt::savefig(base + "_gamma_" + ext);
            }

            plt::figure();
            plt::plot(x, z);
            plt::xlabel("X, ft", 14);
            plt::ylabel("Z, ft", 14);
            if (!state_figure_file.empty())
            {
                plt::savefig(base + "_z_" + ext);
            }
        }
        if (plotting)
        {
            plt::show();
        }
        plt::clear();
    }
    // Nonzero when any trajectory was skipped, so a batch script notices.
    return failed_paths > 0 ? 1 : 0;
}

int main(int argc, char* argv[])
{
    // An input that condemns the whole run -- an unreadable database named
    // before any trajectory is processed, a malformed observer file, a bad
    // combination of flags -- arrives here as an exception rather than an
    // abort, so it is reported and the process exits cleanly.  Failures that
    // cost only one trajectory are caught inside the loop instead.  That puts
    // the decision about what a failure costs in one place, rather than having
    // each call site guess it by choosing fatal() or input_error().
    try
    {
        return run(argc, argv);
    }
    catch (const std::exception& problem)
    {
        error("Stopping: ", problem.what());
        return 2;
    }
}
