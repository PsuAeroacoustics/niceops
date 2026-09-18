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

    std::size_t failed_paths = 0;

    for (std::size_t path_index = 0; path_index < path_count; path_index++)
    {
    // One unreadable track should cost that track, not the whole batch.  The D
    // original wrapped its per-path work the same way (commit fbb3ab6); the
    // failures that reach here are the ones that condemn a single input --
    // trajectory and vehicle files -- rather than the run.
    try
    {
    const auto& path_file = path_files[path_index];
    const auto& vehicle_file = vehicle_files[path_index];
    const auto& out_file = out_files[path_index];
    if (database_files[path_index] != loaded_database)
    {
        log("Loading database ", database_files[path_index], " for ", path_file);
        database = Database(database_files[path_index]);
        // Only once it loaded.  Recording the name first would leave the cache
        // claiming a database that failed to open, so a later trajectory naming
        // the same file would skip the reload and quietly use the previous one.
        loaded_database = database_files[path_index];
    }

    // Create trajectory from CSV data
    trajectory_holder.emplace(path_file, vehicle_file);
    auto& trajectory = *trajectory_holder;

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
            log("--powerfile given but ", loaded_database,
                    " carries no ShaftPower; skipping ", path_file, ".");
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

    parallel_chunks(observers.size(), workers, [&](std::size_t begin, std::size_t end) {
        // One scratch history per worker rather than per observer.
        std::vector<double> levels(nt);
        std::vector<double> penalized;
        std::size_t culled_here = 0;
        for (std::size_t index = begin; index < end; index++)
        {
            const auto& observer = observers[index];
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
                double level = database.level_at_observer(observer, source_states[it]);
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
    });

    if (culling)
    {
        const auto total = observers.size() * nt;
        log("Culled ", culled_samples.load(), " of ", total, " evaluations (",
                100.0 * static_cast<double>(culled_samples.load()) / static_cast<double>(total),
                "%) as provably below ", cull_below, " dB.");
    }

    auto t2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()
            - watch_start)
                        .count();
    log("Acoustic calculations completed in ", (t2 - t1) / 1000.0, " seconds.");
    write_levels(x_obs, y_obs, integrated_levels, out_file, scattered);

    if (!metric_file.empty())
    {
        auto casename = std::filesystem::path(path_file).stem().string();
        append_metrics_file(metric_file, casename, integrated_levels);
    }

    if (write_state)
    {
        trajectory.write_state(state_files[path_index]);
    }
    }
    catch (const std::exception& problem)
    {
        failed_paths++;
        error("Skipping ", path_files[path_index], ": ", problem.what());
    }
    } // per-trajectory loop

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
