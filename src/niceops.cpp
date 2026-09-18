// Port of niceops.d.  Pipeline: load Database + Trajectory -> for each
// observer and each trajectory state, compute instantaneous dBA via
// Database::level_at_observer -> integrate the time history into a sound
// exposure level (SEL) using the 10-dB-down-from-peak window -> write CSV
// output and optional matplotlib figures.

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>
#include <Eigen/Dense>

#include "plotting.hpp"
#include "sphere_database.hpp"
#include "trajectory.hpp"
#include "util/csv.hpp"
#include "util/logger.hpp"
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
    const auto directory = std::filesystem::path(path).parent_path();
    if (!directory.empty() && !std::filesystem::exists(directory))
    {
        std::filesystem::create_directories(directory);
    }
    std::ofstream file(path, mode);
    if (!file)
        fatal("Unable to open file for writing: ", path);
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

void read_observer_file(const std::string& observer_file, std::vector<Eigen::Vector3d>& observers)
{
    trace("Reading observers from ", observer_file, ".");
    // D reads this with csvReader!(double)(contents, null), which consumes the
    // first record as a header and discards it whether or not it is numeric.
    auto records = read_numeric_csv(observer_file, /*skip_first_record=*/true);
    for (const auto& data : records)
    {
        if (data.size() == 3)
        {
            observers.push_back(Eigen::Vector3d(data[0] * ft2m, data[1] * ft2m, data[2] * ft2m));
        }
        else if (data.size() == 4)
        {
            // Skip channel ID in first column
            observers.push_back(Eigen::Vector3d(data[1] * ft2m, data[2] * ft2m, data[3] * ft2m));
        }
        else
        {
            fatal("Observer file format unreadable.");
        }
    }
}

} // namespace

int main(int argc, char* argv[])
{
    // Default to demo case
    std::string log_file = "log/niceops.log";
    std::string database_file = "databases/S-76D_M3.nod";
    std::string path_file = "tracks/S-76D/197421AC.csv";
    std::string vehicle_file = "databases/S-76D.json";
    std::string out_file = "integrated/197421_LAE.csv";
    std::string state_file = "";
    bool write_state = false;
    bool plotting = false;
    std::string figure_file = "";
    std::string state_figure_file = "";
    double contour_min = 45.0;
    double contour_max = 90.0;
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

    CLI::App app{"Noise Informed Community Environment Operation Planning System (NICE-OPS)."};
    app.add_option("-d,--database", database_file, "netCDF noise database file.");
    app.add_option("-l,--logfile", log_file, "Output log file.");
    app.add_option("-t,--trajectoryfile", path_file, "CSV trajectory input file.");
    app.add_option("-v,--vehicle", vehicle_file, "JSON vehicle definition file.");
    app.add_option("-o,--output", out_file, "Output file name.");
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
    app.add_option("-S,--state", state_file, "Filename to save state data");
    app.add_flag("--self_test", self_test, "Run self tests.");

    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        return app.exit(e);
    }

    if (!observer_file.empty())
    {
        scattered = true;
    }
    if (!time_history_prefix.empty())
    {
        write_time_history = true;
    }
    if (!state_file.empty())
    {
        write_state = true;
    }

    // Set default logging to file
    Logger::instance().set_file(log_file);
    Logger::instance().set_level(LogLevel::info);

    // Start stopwatch
    auto watch_start = std::chrono::steady_clock::now();
    auto database = Database(database_file);
    auto t1 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()
            - watch_start)
                        .count();
    log("Database loaded in ", t1 / 1000.0, " seconds.");

    if (self_test)
    {
        // Historical note: this used to fail, and a D-era comment here recorded
        // one sample disagreeing by ~1.7 dB.  The cause was legacy hover
        // directivity (see repair_legacy_hover_directivity); with that repaired
        // the self test passes on every shipped database.
        database.self_test(database_file);
    }

    // Create trajectory from CSV data
    auto trajectory = Trajectory(path_file, vehicle_file);

    // Trim trajectory (per Mike)
    //immutable double xstart_ft = -6000.;
    //trajectory.watts_trim(xstart_ft);

    // Trajectory adaptation parameters, for the adapt() call below.  Both the
    // call and its arguments are disabled, matching the D original.
    //const int degree = 5;
    //const double smooth = 1e6;
    //trajectory.adapt(dt, degree, smooth);

    // And time series data
    const std::size_t nt = trajectory.states.size();
    std::vector<double> time_history(nt);

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

    std::vector<double> integrated_levels;
    std::vector<double> x_obs;
    std::vector<double> y_obs;
    std::vector<Eigen::Vector3d> observers;
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
        x_obs = iota_exclusive(-1500.0, 1500.0, 25.0);
        y_obs = iota_exclusive(-1500.0, 1500.0, 25.0);
        const std::size_t nx = x_obs.size();
        const std::size_t ny = y_obs.size();
        // Set corresponding length for levels
        integrated_levels.resize(nx * ny);
        observers.reserve(nx * ny);
        for (std::size_t iy = 0; iy < ny; iy++)
        {
            for (std::size_t ix = 0; ix < nx; ix++)
            {
                observers.push_back(Eigen::Vector3d(x_obs[ix], y_obs[iy], 0.0));
            }
        }
    }
    // The source geometry and the blend of database operating conditions depend
    // only on the trajectory state, so build them once here rather than once
    // per (observer, state) pair.  The blend is a CGAL Delaunay locate() walk
    // and was by far the largest single cost in the program; this turns
    // observers*states of them into states of them.
    std::vector<Database::SourceState> source_states;
    source_states.reserve(nt);
    for (const auto& state : trajectory.states)
    {
        source_states.push_back(database.source_state(Eigen::Vector3d(state.x, state.y, state.z),
                Eigen::Vector3d(state.vx, state.vy, state.vz),
                Eigen::Vector3d(state.nx, state.ny, state.nz), state.mu,
                qsam ? state.fpa_e : state.fpa, state.CT));
    }
    // ambient_level stays NaN unless -a was given, which is how D's default-NaN
    // double suppressed the ambient mix-in.  The reference pressure it converts
    // to is loop-invariant; it used to be a pow() per evaluation.
    const bool add_ambient = !std::isnan(ambient_level);
    const double ambient_pressure = add_ambient ? std::pow(10.0, ambient_level / 20.) : 0.0;

    for (std::size_t index = 0; index < observers.size(); index++)
    {
        const auto& observer = observers[index];
        for (std::size_t it = 0; it < nt; it++)
        {
            double level = database.level_at_observer(observer, source_states[it]);
            if (add_ambient)
            {
                // Add ambient noise
                level = 20.0 * std::log10(std::pow(10.0, level / 20.) + ambient_pressure);
            }
            time_history[it] = level;
        }
        integrated_levels[index] = sound_exposure_level(time_history, integration_step);
        if (write_time_history)
        {
            // Write out time history files
            auto time_history_filename = time_history_prefix + "_i" + std::to_string(index) + ".csv";
            write_time_history_file(time_history, integration_step, observer(0), observer(1),
                    time_history_filename);
        }
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
        trajectory.write_state(state_file);
    }

    if (plotting || !figure_file.empty() || !state_figure_file.empty())
    {
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
            //plt::plot(x, y, "k-");
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
    return 0;
}
