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
#include <cstdlib>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>
#include "util/iota.hpp"
#include "vec3.hpp"

#include "plotting.hpp"
#include <complex>
#include <sstream>

#include "atmosphere.hpp"
#include "atmosphere_profile.hpp"
#include "refraction.hpp"
#include "stratified.hpp"
#include "ground_effect.hpp"
#include "plate_receiver.hpp"
#include "propagation_table.hpp"
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
        const std::vector<double>& levels, const std::string& out_file, bool scattered,
        LengthUnit unit)
{
    info("Writing output to ", out_file);
    auto file = open_output(out_file);

    if (scattered)
    {
        info("Scattered output.");
        assert(levels.size() == x_obs.size() && levels.size() == y_obs.size());
        file << "x(" << unit.name << "), y(" << unit.name << "), Level(dB)\n";
        for (std::size_t i = 0; i < levels.size(); i++)
        {
            file << x_obs[i] / unit.to_meters << "," << y_obs[i] / unit.to_meters << ","
                 << levels[i] << "\n";
        }
    }
    else
    {
        info("Gridded output.");
        const std::size_t nx = x_obs.size();
        const std::size_t ny = y_obs.size();
        assert(nx * ny == levels.size());

        // Write header
        file << "Line 1: x(" << unit.name << "), Line 2: y(" << unit.name
             << "),  Lines 3+ levels (x cols y rows)\n";
        // Comma-separated, newline at the end of each row.
        auto write_row = [&file](const std::vector<double>& values, std::size_t first,
                                 std::size_t count, double scale) {
            for (std::size_t i = 0; i < count; i++)
                file << (i ? "," : "") << values[first + i] / scale;
            file << "\n";
        };
        write_row(x_obs, 0, nx, unit.to_meters);
        write_row(y_obs, 0, ny, unit.to_meters);
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
        const std::string& filename, LengthUnit unit)
{
    info("Writing time history to ", filename);
    auto file = open_output(filename);
    file << "Level time history. Next row [dt (s), x (" << unit.name << "), y(" << unit.name
         << ")], remaining rows, dB.\n";
    file << dt << "," << x / unit.to_meters << "," << y / unit.to_meters << "\n";
    for (auto level : time_history)
    {
        file << level << "\n";
    }
}

// One observer's per-band history (--band_history): a row per trajectory
// state, each band's unweighted level at the microphone beside the geometry it
// was reached along, so a measured spectrum can be put against it frame by
// frame.  Bands the run did not evaluate (outside --band_limits) and culled
// states are left empty.
struct BandHistoryRow
{
    double emission_time = 0.0;
    double reception_time = 0.0;
    double range = 0.0;
    double elevation = 0.0;
    double phi = 0.0;
    double theta = 0.0;
};

void write_band_history_file(const std::vector<BandHistoryRow>& rows,
        const std::vector<double>& levels, const std::vector<double>& band_levels,
        const std::vector<double>& band_centers, std::size_t first_band, std::size_t last_band,
        const Vec3& microphone, LengthUnit unit, const std::string& filename)
{
    auto file = open_output(filename);
    file << "# NICE-OPS band history. Microphone x, y, z (" << unit.name << "): "
         << microphone.x / unit.to_meters << ", " << microphone.y / unit.to_meters << ", "
         << microphone.z / unit.to_meters << ". Times s, range " << unit.name
         << ", angles deg; level_dBA is the A-weighted, band-limited level, band_<Hz> each "
            "band's unweighted level at the microphone, dB re 20 uPa.\n";
    file << "emission_time,reception_time,range,elevation,phi,theta,level_dBA";
    for (auto center : band_centers)
        file << ",band_" << center;
    file << "\n";
    file << std::setprecision(10);
    const std::size_t nb = band_centers.size();
    for (std::size_t it = 0; it < rows.size(); it++)
    {
        const auto& row = rows[it];
        const bool evaluated = std::isfinite(levels[it]);
        file << row.emission_time << ",";
        if (evaluated)
        {
            file << row.reception_time << "," << row.range / unit.to_meters << ","
                 << row.elevation << "," << row.phi << "," << row.theta << "," << levels[it];
        }
        else
        {
            file << ",,,,," << levels[it];
        }
        for (std::size_t band = 0; band < nb; band++)
        {
            file << ",";
            if (evaluated && band >= first_band && band < last_band)
                file << band_levels[it * nb + band];
        }
        file << "\n";
    }
}

void append_metrics_file(const std::string& filename, const std::string& casename,
        const std::vector<double>& levels, double threshold = 65.0)
{
    auto file = open_output(filename, std::ios::app);
    double energy = 0;
    double maximum = -std::numeric_limits<double>::infinity();
    std::size_t count = 0;
    // The average is an energy mean, 10 lg of the mean of 10^(L/10): an
    // arithmetic mean of decibels is not an average of noise.
    for (auto level : levels)
    {
        energy += std::pow(10.0, 0.1 * level);
        maximum = std::max(maximum, level);
        if (level >= threshold)
        {
            count++;
        }
    }
    const double average = 10.0 * std::log10(energy / levels.size());
    auto area = double(count) / levels.size();
    file << casename << "," << average << "," << maximum << "," << area << "\n";
}

// Sound exposure level at one observer without marching every emission time.
//
// The integral only ever uses the 10 dB-down duration, which on a shipped case
// is a little over half the track -- the rest is computed and thrown away.  The
// bound does two jobs: it proves the peak is the peak, and it proves a state
// below the threshold, so the scans for the window's first and last states can
// step over most of the track unevaluated.  Everything between those two is
// evaluated, since dips below the threshold are inside the duration.
// DirectionBounds is tight enough (median slack 8.6 dB, inside the 10 dB
// window) that this costs little beyond the window itself.
//
// The result is identical to evaluating every state, not an approximation: a
// skipped state is provably neither the peak nor an end of the window, and
// nothing between the ends is skipped.  The integral is in reception time --
// see reception_factor in sel.hpp.
//
// `geometry_of(it)` gives state `it`'s geometry at this observer: at least the
// squared range and the three components, in the sphere frame, of the
// displacement the level spreads over and looks the sphere up along.  For a
// straight ray that is the displacement itself; for a refracted one it is the
// launch direction scaled by the ray tube's range (refraction.hpp), which is
// what makes the same DirectionBounds bound it.  `level_of(it, geometry)`
// evaluates the state from it, and `bounds` must bound what it returns.
// `range_at(it)` is the range whose rate of change gives reception time.
//
// The bound test needs exactly what an evaluation needs, which is why the
// geometry is formed once and handed on rather than recomputed -- doing it
// twice was measured at 1.48x down to 1.08x.
struct StateGeometry
{
    double squared_range, along, up, across;
};

template <typename GeometryOf, typename LevelOf, typename RangeAt>
double windowed_exposure(std::size_t nt, const DirectionBounds& bounds,
        double sphere_radius_squared, double step, double sound_speed, std::vector<double>& levels,
        std::vector<unsigned>& visited, std::vector<double>& skipped_numerator,
        std::vector<double>& skipped_squared_range, unsigned generation, std::size_t& evaluated,
        GeometryOf&& geometry_of, LevelOf&& level_of, RangeAt&& range_at)
{
    // A coarse sweep seeds the peak search; the stride is not sensitive -- 16
    // and 1024 land within a couple of percent of each other.
    const std::size_t seed_stride = std::max<std::size_t>(1, nt / 128);

    using Geometry = decltype(geometry_of(std::size_t {0}));
    auto evaluate_with = [&](std::size_t it, Geometry&& g) {
        visited[it] = generation;
        evaluated++;
        return levels[it] = level_of(it, std::move(g));
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
    for (std::size_t it = 0; it < nt; it += seed_stride)
        best = std::max(best, evaluate(it));
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
        Geometry g = geometry_of(it);
        // The bound subtracts spherical spreading, which is only a loss outside
        // the sphere; inside it, evaluate rather than trust it.
        if (g.squared_range > sphere_radius_squared)
        {
            const double numerator = bounds.numerator(
                    it, DirectionBounds::bin(g.along, g.up, g.across, g.squared_range));
            if (numerator < best_energy * g.squared_range)
            {
                // Kept for the scans below.  Every state this pass leaves
                // unevaluated is skipped here, so they all have these two by the
                // time the window's ends are looked for; recomputing them there
                // (a full ray solve, under refraction) was 30 percent of the
                // bound work and the reason this path lost to marching.
                skipped_numerator[it] = numerator;
                skipped_squared_range[it] = g.squared_range;
                continue;
            }
        }
        const double level = evaluate_with(it, std::move(g));
        if (level > best)
        {
            best = level;
            best_energy = std::pow(10.0, 0.1 * best);
        }
    }

    // `best` is now provably the peak of the whole history, so the window and
    // the integral are the ordinary ones -- ten_db_down_interval() and
    // sel_over_window() from sel.hpp, reading through `evaluate` instead of a
    // filled array.  Sharing them is what keeps this path's window identical to
    // peak_range's rather than merely intended to be.
    //
    // The window runs from the first to the last state at or above the
    // threshold, so its ends are found by scanning in from both ends of the
    // track.  A state the bound puts provably below the threshold cannot be an
    // end and is skipped unevaluated; every state *between* the ends is then
    // evaluated, dips included, since they are inside the duration.
    const double threshold = best - sel_window_db;
    const double threshold_energy = std::pow(10.0, 0.1 * threshold);
    auto certainly_below = [&](std::size_t it) {
        if (visited[it] == generation)
            return levels[it] < threshold;
        // Not evaluated, so the main pass skipped it outside the sphere, and left
        // its numerator and squared range: the same test as before, on the same
        // two numbers, so the same answer.
        return skipped_numerator[it] < threshold_energy * skipped_squared_range[it];
    };
    const auto index = ten_db_down_interval(nt, threshold, evaluate, certainly_below);
    return sel_over_window(index[0], index[1], step, evaluate, [&](std::size_t it) {
        return reception_factor(it, nt, step, sound_speed, range_at);
    });
}

// The straight-ray form: the displacement from each state to the observer.
template <typename LevelOf>
double windowed_sound_exposure_level(const Database& database,
        const std::vector<Database::SourceState>& states, const DirectionBounds& bounds,
        const Vec3& observer, double step, double sound_speed, std::vector<double>& levels,
        std::vector<unsigned>& visited, std::vector<double>& skipped_numerator,
        std::vector<double>& skipped_squared_range, unsigned generation, std::size_t& evaluated,
        LevelOf&& level_of)
{
    return windowed_exposure(
            states.size(), bounds, database.sphere_radius() * database.sphere_radius(), step,
            sound_speed, levels, visited, skipped_numerator, skipped_squared_range, generation,
            evaluated,
            [&](std::size_t it) {
                const Vec3 delta = observer - states[it].position;
                return StateGeometry {delta.squaredNorm(), states[it].frame.nose.dot(delta),
                        states[it].frame.normal.dot(delta), states[it].frame.starboard.dot(delta)};
            },
            std::forward<LevelOf>(level_of),
            [&](std::size_t it) { return (observer - states[it].position).norm(); });
}

void read_observer_file(const std::string& observer_file, std::vector<Vec3>& observers,
        std::vector<double>& microphone_heights, LengthUnit unit)
{
    trace("Reading observers from ", observer_file, ".");
    // D reads this with csvReader!(double)(contents, null), which consumes the
    // first record as a header and discards it whether or not it is numeric.
    CsvReport report;
    auto records = read_numeric_csv(observer_file, /*skip_first_record=*/true, &report);
    // A row that did not parse used to be dropped without a word, so a typo
    // ("4oo") removed an observer and everything after it moved up a place.
    if (report.dropped > 0)
    {
        input_error("Observer file ", observer_file, ": line ", report.first_dropped_line,
                " is not a row of numbers", report.dropped > 1
                        ? " (and " + std::to_string(report.dropped - 1) + " more are not)"
                        : std::string(),
                ". After the header every line must be x,y,z, id,x,y,z or id,x,y,z,height.");
    }
    // One shape for the whole file: a 3-field row is x,y,z and a 4-field one is
    // id,x,y,z, and the width was decided row by row, so a row cut short became
    // an observer at the wrong coordinates (x,y,z of "4,700,800" is 4 ft east,
    // 700 ft north, 800 ft up).
    for (std::size_t row = 1; row < records.size(); row++)
    {
        if (records[row].size() != records.front().size())
        {
            input_error("Observer file ", observer_file, ": row ", row + 2, " has ",
                    records[row].size(), " fields where the first has ", records.front().size(),
                    ". Every row must have the same columns.");
        }
    }
    for (const auto& data : records)
    {
        if (data.size() == 3)
        {
            observers.push_back(Vec3(data[0] * unit.to_meters, data[1] * unit.to_meters,
                    data[2] * unit.to_meters));
        }
        else if (data.size() == 4)
        {
            // Skip channel ID in first column
            observers.push_back(Vec3(data[1] * unit.to_meters, data[2] * unit.to_meters,
                    data[3] * unit.to_meters));
        }
        else if (data.size() == 5)
        {
            // id, x, y, z, and the microphone's height above the ground there:
            // 0 for a flush microphone (a plate, with --plate_table), more for
            // one on a pole.
            observers.push_back(Vec3(data[1] * unit.to_meters, data[2] * unit.to_meters,
                    data[3] * unit.to_meters));
            if (!(data[4] >= 0.0))
            {
                input_error("Observer file ", observer_file, " gives a microphone height of ",
                        data[4], " ", unit.name, "; it must not be negative.");
            }
            microphone_heights.push_back(data[4] * unit.to_meters);
        }
        else
        {
            input_error("Observer file format unreadable.");
        }
    }
    if (!microphone_heights.empty() && microphone_heights.size() != observers.size())
    {
        input_error("Observer file ", observer_file, " gives a microphone height on some rows "
                "and not others; give it on every row (id,x,y,z,height) or on none.");
    }
    if (observers.empty())
    {
        // Every row failing to parse leaves this empty, and the run then
        // "succeeds" having evaluated nothing and written a lone header line.
        // The first record is consumed as a header whether or not it holds
        // numbers, so a one-line file lands here too.
        input_error("Observer file ", observer_file, " yielded no observers. The first record "
                "is taken as a header and discarded; every row after it must be x,y,z, "
                "id,x,y,z or id,x,y,z,height in ", unit.name, ".");
    }
}

// The observers a run evaluates, and their coordinates as the CSV writers want
// them: for a grid, the two axes; for a scattered set, one entry per observer.
struct ObserverSet
{
    std::vector<Vec3> observers;
    std::vector<double> x;
    std::vector<double> y;
    // Each observer's microphone height above the ground, meters, when the
    // observer file gives one (id,x,y,z,height); empty otherwise, and every
    // observer takes the run's receiver (Receivers in process_trajectory).
    std::vector<double> microphone_heights;
};

ObserverSet scattered_observers(const std::string& observer_file, LengthUnit unit)
{
    ObserverSet set;
    read_observer_file(observer_file, set.observers, set.microphone_heights, unit);
    set.x.reserve(set.observers.size());
    set.y.reserve(set.observers.size());
    for (const auto& observer : set.observers)
    {
        set.x.push_back(observer(0));
        set.y.push_back(observer(1));
    }
    return set;
}

ObserverSet gridded_observers(const std::vector<double>& bounds, double step)
{
    // Mirrors D's iota(-1500.0, 1500.0, 25.0).
    // NOTE: D's iota(begin, end, step) is EXCLUSIVE of end -- it yields
    // ceil((end-begin)/step) values, so the default is 120 points spanning
    // -1500..1475 m, NOT 121 points ending at 1500.  Verified against LDC.
    ObserverSet set;
    set.x = iota_exclusive(bounds[0], bounds[1], step);
    set.y = iota_exclusive(bounds[2], bounds[3], step);
    // Row-major over (y, x), which is the order write_levels reads back.
    set.observers.reserve(set.x.size() * set.y.size());
    for (const double y : set.y)
        for (const double x : set.x)
            set.observers.push_back(Vec3(x, y, 0.0));
    return set;
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
    // dBA with the database's fitted EAA factor.  On by default; a database
    // without spectra falls back to broadband unless spectral_explicit says
    // --spectral was asked for by name.
    bool spectral = true;
    bool spectral_explicit = false;
    BandLimits band_limits;
    // Excess ground attenuation.  On by default, but applied only where the
    // spectral path runs on a database with spectra (process_trajectory decides
    // that per trajectory); ground_effect_given says the user asked for it, so a
    // request that cannot be honored is an error rather than a silent no-op.
    // It is the single most expensive thing in the evaluation.
    bool ground_effect = true;
    bool ground_effect_given = false;
    // A fixed_load_factor database's load-factor law is an extrapolation above
    // the load factors its source supports; it is held at this many g.
    double max_load_factor = 2.0;
    // Warn, rather than refuse, a vehicle file that disagrees with the database
    // by more than 2 percent.
    bool accept_vehicle_mismatch = false;
    GroundModel ground;
    // Evaluate the point-receiver ground effect directly, not from GroundTable.
    bool exact_ground = false;
    // A flush microphone on a ground plate (--plate_table), read once and
    // bound to each trajectory's atmosphere.  Needs ground_effect.
    std::optional<PlateReceiver> plate;
    PlateOrientation plate_orientation;
    // A propagation loss table (--propagation_table), read once for the run:
    // where it covers a pair it replaces the ground model and the ray model
    // (propagation_table.hpp).
    std::shared_ptr<const PropagationTable> propagation_table;
    std::string propagation_table_path;
    // Overrides for the atmosphere read from the vehicle JSON; NaN means "use
    // whatever the file said".
    double humidity_override = std::numeric_limits<double>::quiet_NaN();
    double temperature_override = std::numeric_limits<double>::quiet_NaN();
    double pressure_override = std::numeric_limits<double>::quiet_NaN();
    // A stratified atmosphere (--atmosphere), read once for the whole run.
    // It replaces the vehicle file's single atmosphere for absorption, which
    // is integrated through it (LayeredAbsorption); its surface level stands
    // in for the rest.  Pressure may be missing, in which case each trajectory
    // fills it hydrostatically from its own surface pressure.
    std::optional<AtmosphereProfile> atmosphere;
    // --wind: the profile is built per trajectory from one surface observation
    // instead (synthetic_profile), since its surface air is the trajectory's
    // vehicle file's and its top has to clear the trajectory.  Never set
    // together with `atmosphere`.
    std::optional<SurfaceWind> surface_wind;
    // --refraction linear: curved rays through the profile's equivalent linear
    // effective sound speed, fitted per path (refraction.hpp).  --refraction
    // stratified: through the profile's own layers (stratified.hpp).  Both need
    // the profile.  The shadow zone's attenuation is capped at shadow_floor dB.
    bool refraction = false;
    bool stratified = false;
    // --tunneling, --skip_fill (stratified only): arrivals ray theory leaves
    // out, with their Airy tails per band (layered_wave_arrivals).
    WaveOptions waves;
    double shadow_floor = 20.0;
    // --scattering: the HARMONOISE turbulence strength whose scattered energy
    // fills the shadows (ShadowScattering); 0 keeps shadow_floor instead.
    double scattering = 0.0;
    // Compass bearing of the local +x axis, degrees, for turning the profile's
    // wind into the local frame; NaN means --frame_bearing was not given and
    // +x is east.
    double frame_bearing = std::numeric_limits<double>::quiet_NaN();
    // Which velocity the flight conditions are taken from (--speed_reference).
    // Unset follows each database's own speed_reference attribute.
    std::optional<SpeedReference> speed_reference;
    // Meters, after the --units conversion; NaN means the trim is off.
    double trim_start = 0.0;
    double utc_offset = 0.0;
    bool scattered = false;
    bool want_dnl = false;
    bool want_lamax = false;
    bool want_nax = false;
    double nax_threshold = std::numeric_limits<double>::quiet_NaN();
    bool want_coverage = false;
    bool full_history = false;
    bool write_time_history = false;
    bool write_state = false;
    std::string time_history_prefix;
    // Per-band histories (--band_history); empty for none.
    std::string band_history_prefix;
    std::string metric_file;
    std::string power_file;
    // How WGS84 trajectory input is projected, and the frame it is projected
    // onto.  The frame is settled once for the whole run: the observers are
    // common to every trajectory in it, so their coordinates have to mean the
    // same thing in all of them.
    GeodeticInput geodetic;
    // Write a WGS84 companion beside every level file.
    bool latlon = false;
    // The length unit the command line and the files speak (--units).
    LengthUnit unit;
};

// The files one trajectory is computed from and written to, and where it sits
// in the run.
struct PathJob
{
    std::string trajectory;
    std::string vehicle;
    std::string output;
    std::string state;
    // --write_atmosphere: where the --wind profile goes; empty for nowhere.
    std::string atmosphere_out;
    std::size_t count = 1;
};

// Take the trajectory's flight conditions against the air when its database's
// were: the database's speed_reference attribute, unless --speed_reference
// overrides it.  The wind at each state is the profile's at its height (the
// --atmosphere file's, or the one --wind built), turned into the local frame by
// --frame_bearing.  A run with no wind profile cannot refer anything to the
// air, and treating it as calm would quietly mix the two references, so it is
// refused -- --speed_reference ground says so explicitly.  A --wind observation
// of zero is the exception: calm was stated, not assumed, and the air velocity
// is the ground velocity.  Returns whether the trajectory is now air-referenced.
bool refer_to_database_speed(Trajectory& trajectory, const Database& database,
        const RunOptions& options, const AtmosphereProfile* atmosphere)
{
    const SpeedReference built = database.speed_reference();
    const SpeedReference wanted = options.speed_reference.value_or(built);
    if (wanted != built)
    {
        log("--speed_reference ", to_string(wanted), " overrides the database, whose flight "
                "conditions are referred to the ", to_string(built), ".");
    }
    if (wanted == SpeedReference::ground)
        return false;
    if (options.surface_wind && options.surface_wind->speed == 0.0)
    {
        log("Flight conditions referred to the air, which --wind says is calm: the air velocity is "
            "the ground velocity.");
        return false;
    }
    if (!atmosphere || !atmosphere->has_wind())
    {
        input_error(options.speed_reference ? "--speed_reference air"
                                            : "This database is referred to the air, and",
                " needs the wind at the aircraft: give --atmosphere with wind columns, or --wind, "
                "or --speed_reference ground to take the ground velocity as the airspeed.");
    }
    const CompassToLocal to_local(std::isnan(options.frame_bearing) ? 90.0 : options.frame_bearing);
    const AtmosphereProfile& profile = *atmosphere;
    const double strongest = trajectory.relative_to_air(
            [&](double z) { return to_local(profile.wind_at(z)); });
    log("Flight conditions referred to the air: mu, the path angle and the sphere's nose follow "
        "the velocity less the profile's wind, which reaches ",
            strongest / options.unit.to_meters, " ", options.unit.name, "/s at the aircraft.");
    if (std::isnan(options.frame_bearing))
    {
        // For an AC track, whose x/y are rotated, this is a mistake rather than a
        // remark: the wind is turned by the wrong angle and the footprint is off
        // by several dB with only a line in the log to say why.  A geodetic run
        // has settled its own frame (east/north unless told otherwise), so there
        // it stays a note.
        const char* note = "  the wind directions are compass bearings, turned into the local frame "
                           "taking +x as east. If this run's x/y are rotated -- the shipped AC tracks "
                           "run +x along 210 degrees -- give --frame_bearing.";
        if (options.geodetic.frame)
            log(note);
        else
            warn(note);
    }
    return true;
}

// The --wind profile for one trajectory: the surface air is the trajectory's
// own (its vehicle file's, under --temperature, --humidity and --pressure), and
// the top clears the highest state and microphone.  Written to
// job.atmosphere_out when --write_atmosphere asked for it.
AtmosphereProfile wind_profile(const Trajectory& trajectory, const PathJob& job,
        const RunOptions& options, const std::vector<Vec3>& observers,
        const std::vector<double>& microphone_heights)
{
    const SurfaceWind& wind = *options.surface_wind;
    const auto& vehicle = trajectory.vehicle;
    Atmosphere surface = vehicle.atmosphere;
    if (!std::isnan(options.temperature_override))
        surface.temperature = options.temperature_override;
    if (!std::isnan(options.pressure_override))
        surface.pressure = options.pressure_override;
    if (!std::isnan(options.humidity_override))
        surface.relative_humidity = options.humidity_override;
    // The same refusal the spectral path makes without a profile: the whole
    // humidity profile is carried up from this one value.
    if (std::isnan(options.humidity_override) && !vehicle.humidity_given)
    {
        input_error("--wind builds the humidity profile from the surface humidity, and none is set "
                    "(no --humidity, and ", job.vehicle, " gives none). Give --humidity, or add "
                    "\"relative humidity\" to the vehicle file or its Atmosphere section.");
    }
    if (const auto problem = temperature_problem(surface.temperature); !problem.empty())
        input_error("--wind: the surface temperature from ", job.vehicle, ": ", problem);
    if (const auto problem = pressure_problem(surface.pressure); !problem.empty())
        input_error("--wind: the surface pressure from ", job.vehicle, ": ", problem);
    if (const auto problem = humidity_problem(surface.relative_humidity); !problem.empty())
        input_error("--wind: the surface humidity from ", job.vehicle, ": ", problem);

    // Every height a path reaches.  The microphone heights are taken at their
    // largest, which can only raise the top.
    double highest = 0.0;
    for (const auto& state : trajectory.states)
        highest = std::max(highest, state.z);
    const double receiver = std::max(options.ground.receiver_height,
            microphone_heights.empty()
                    ? 0.0
                    : *std::max_element(microphone_heights.begin(), microphone_heights.end()));
    for (const auto& observer : observers)
        highest = std::max(highest, observer.z + receiver);

    const AtmosphereProfile profile = synthetic_profile(wind, surface, highest);
    const LengthUnit unit = options.unit;
    const std::size_t top = profile.height.size() - 1;
    log("Atmosphere built from --wind: ", wind.speed, " m/s (", wind.speed * 3600.0 / 1852.0,
            " kt) from ", wind.direction, " degrees at ", wind.height / unit.to_meters, " ",
            unit.name, ", roughness length ", wind.roughness_length / unit.to_meters, " ", unit.name,
            ", ", describe(wind.stability), "; ", profile.height.size(), " levels to ",
            profile.top() / unit.to_meters, " ", unit.name, " (", wind.levels_per_decade,
            " per factor of ten from ", wind.first_level() / unit.to_meters, " ", unit.name, ").");
    log("  wind ", wind.speed_at(wind.surface_layer_depth), " m/s at the top of the surface "
            "layer, ", wind.surface_layer_depth / unit.to_meters, " ", unit.name, ", and held above "
            "it; ", profile.temperature.front(), " K, ", profile.pressure.front(), " kPa, ",
            profile.relative_humidity.front(), "% RH at the ground, and ", profile.temperature[top],
            " K, ", profile.pressure[top], " kPa, ", profile.relative_humidity[top],
            "% RH at the top (water vapor held).");
    if (*std::max_element(profile.relative_humidity.begin(), profile.relative_humidity.end())
            >= 100.0)
    {
        warn("--wind: the surface air saturates on its way up to ", profile.top() / unit.to_meters,
                " ", unit.name, "; the relative humidity is held at 100 percent above that.");
    }
    if (!job.atmosphere_out.empty())
    {
        // --humidity 0 is completely dry air, which humidity_problem accepts;
        // read_atmosphere_profile refuses a file whose humidity never exceeds 1
        // percent, taking it for a fraction.  The run is right, but its file
        // will not replay, and that should be said now rather than then.
        if (*std::max_element(profile.relative_humidity.begin(), profile.relative_humidity.end())
                <= 1.0)
        {
            warn("--write_atmosphere: ", job.atmosphere_out, " is completely dry air, and "
                    "--atmosphere refuses a profile whose relative humidity never exceeds 1 percent "
                    "(it takes it for a fraction), so it will not replay as written.");
        }
        write_atmosphere_profile(profile, job.atmosphere_out, unit);
        log("  written to ", job.atmosphere_out, ".");
    }
    return profile;
}

// Compute one trajectory's footprint and write its outputs, returning the
// trajectory itself for the single-path plotting block.  Anything that throws
// in here condemns this trajectory alone; the caller catches and moves on.
Trajectory process_trajectory(const Database& database, const PathJob& job,
        const RunOptions& options, const std::vector<Vec3>& observers,
        const std::vector<double>& microphone_heights, const std::vector<double>& x_obs,
        const std::vector<double>& y_obs, std::vector<double>& integrated_levels,
        std::vector<double>& dnl_energy, std::vector<double>& lamax, std::vector<double>& nax,
        std::vector<std::array<double, 2>>& coverage)
{
    // Create trajectory from CSV data
    Trajectory trajectory(job.trajectory, job.vehicle, options.geodetic, options.unit);
    // Ground effect is on by default, but it acts on per-band levels, so it is
    // applied only where the spectral path runs on a database that has spectra.
    const bool ground_effect = options.ground_effect && options.spectral && database.has_spectra();

    // mu and CT are non-dimensionalized by the vehicle JSON's tip speed, and
    // they index the database's operating conditions -- which were generated at
    // whatever tip speed the database records.  Disagree, and every lookup lands
    // at the wrong condition with nothing downstream able to notice.  Worth
    // checking because the two are easy to state in different units:
    // S-76D_M3.nod was built at 220 m/s, and 600 is that same figure read as
    // ft/s, which puts the tip at Mach 1.7.
    //
    // Weight and rotor radius set CT and mu in the same way and the database
    // records them too, but only the tip speed was ever compared: a vehicle file
    // 10 percent heavy moved the footprint 0.75 dB median with nothing said, and
    // a tip speed 4.6 percent off moved it 1.1 dB.  All three are warnings, on
    // stderr.  The weight is compared at 1 percent: panam builds with g = 9.82
    // and this program uses 9.80665, a 0.14 percent difference that is not a
    // disagreement.
    // Past 2 percent the run is refused, since the footprint would be wrong by
    // more than the grid's own error; between 1 and 2 it warns.  The refusal
    // names the same check as the warning, and --accept_vehicle_mismatch turns
    // it back into a warning.
    const auto vehicle_disagrees = [&](const char* what, const char* unit_text,
                                           double from_file, double from_database) {
        if (std::isnan(from_database) || approx_equal(from_file, from_database))
            return;
        const double relative = std::fabs(from_file - from_database) / std::fabs(from_database);
        if (relative > 0.02 && !options.accept_vehicle_mismatch)
        {
            input_error("Vehicle ", what, " disagrees: ", job.vehicle, " gives ", from_file, " ",
                    unit_text, " but the database was built at ", from_database, " ", unit_text,
                    ", ", relative * 100.0, " percent off. mu and CT would index the wrong "
                    "conditions, so the run is refused. Give the right vehicle file, or "
                    "--accept_vehicle_mismatch to propagate it anyway.");
        }
        warn("Vehicle ", what, " disagrees: ", job.vehicle, " gives ", from_file, " ", unit_text,
                " but the database was built at ", from_database, " ", unit_text,
                ". mu and CT will index the wrong conditions.");
    };
    vehicle_disagrees("main rotor tip speed", "m/s", trajectory.vehicle.tip_speed, database.tip_speed());
    vehicle_disagrees("main rotor radius", "m", trajectory.vehicle.radius, database.rotor_radius());
    vehicle_disagrees("weight", "N", trajectory.vehicle.weight, database.vehicle_weight());
    // The thrust coefficient the vehicle file gives at a load factor of 1,
    // against a fixed_load_factor database's reference.  It folds in the
    // density, which no check above sees.
    if (database.fixed_load_factor())
    {
        const auto& vehicle = trajectory.vehicle;
        const double ct_at_one_g = vehicle.weight
                / (vehicle.rho * vehicle.tip_speed * vehicle.tip_speed * pi * vehicle.radius * vehicle.radius);
        vehicle_disagrees("thrust coefficient at 1 g", "", ct_at_one_g, database.reference_thrust_coefficient());
    }

    // Trim trajectory (per Mike): crop from the given x to whichever comes
    // first, the lowest or the slowest point.  Off unless asked for -- it
    // changes which part of the flight is integrated, and so the levels.
    if (!std::isnan(options.trim_start))
    {
        trajectory.watts_trim(options.trim_start);
    }

    // Resample onto the -s grid, smoothing the path first.  Disabled, matching
    // the D original; --path_smoothing covers the geodetic input that needs it.
    // The arguments are (dt, smoothing, tolerance): a negative tolerance uses
    // `smoothing` as the penalty weight directly, and a non-negative one solves
    // for the weight instead.  The D-era pair here was (degree, smooth) from
    // the vendored B-spline library, whose second and third arguments meant
    // something else entirely -- uncommented as it stood, it asked for a fit
    // allowed to sit 1e6 away from the measured path.
    //trajectory.adapt(options.dt, /*smoothing=*/0.3);

    // Refer mu, the path angle and the sphere's nose to the velocity the
    // database's conditions were taken from.  After the trim, which cuts the
    // approach at its slowest point over the ground, not through the air.
    const std::optional<AtmosphereProfile> synthesized = options.surface_wind
            ? std::optional<AtmosphereProfile>(wind_profile(trajectory, job, options, observers,
                      microphone_heights))
            : std::nullopt;
    // The profile this trajectory propagates through, whichever option gave it.
    const std::optional<AtmosphereProfile>& atmosphere
            = options.surface_wind ? synthesized : options.atmosphere;
    const bool air_referenced = refer_to_database_speed(trajectory, database, options,
            atmosphere ? &*atmosphere : nullptr);

    // And time series data
    const std::size_t nt = trajectory.states.size();

    // The sound exposure level integral weights each sample by the trajectory's
    // own time step, so take that from the states rather than from -s.  Those
    // two coincide only when adapt() resamples onto the -s grid, and adapt() is
    // disabled below; every shipped track is sampled at 0.02 s, so the 0.5 s
    // default would overstate every level by 10*log10(0.5/0.02) = 14 dB.
    double integration_step = options.dt;
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
    // per (observer, state) pair.  The blend is a Delaunay locate() walk
    // and was by far the largest single cost in the program; this turns
    // observers*states of them into states of them.
    // Night classification for the DNL penalty, from the AC format's UTC
    // seconds column shifted into local time.
    std::vector<char> night(options.want_dnl ? nt : 0);
    for (std::size_t it = 0; it < night.size(); it++)
    {
        if (std::isnan(trajectory.states[it].utcsec))
        {
            // A NaN fails both halves of is_night's window test, so a track with
            // no clock was classed as daytime -- or, when the start defaulted to
            // zero, as UTC midnight, i.e. night: +10 dB on the whole footprint,
            // silently.  Either way the answer was made up.
            input_error("-D needs the time of day, and ", job.trajectory, " does not carry one: "
                    "it has no utcsec column. Give --start_utc <seconds since UTC midnight at "
                    "t = 0>.");
        }
        night[it] = is_night(trajectory.states[it].utcsec, options.utc_offset) ? 1 : 0;
    }

    // The sphere's nose: the air velocity, or at hover the heading (nose_direction).
    const auto offset = heading_offset(trajectory.states, options.frame_bearing);
    std::size_t by_heading = 0, slow = 0, load_capped = 0;
    std::size_t mu_clamped = 0, fpa_clamped = 0, ct_clamped = 0, load_scaled = 0;
    std::vector<Database::SourceState> source_states;
    source_states.reserve(nt);
    for (const auto& state : trajectory.states)
    {
        const auto velocity = air_velocity(state);
        const auto nose = nose_direction(state, offset);
        if (std::hypot(velocity[0], velocity[1]) < velocity_nose_speed)
        {
            slow++;
            if (nose != velocity)
                by_heading++;
        }
        // On a fixed_load_factor database the load-factor law is the analytic
        // correction 20 lg(C_T / C_T,ref), an extrapolation above the load
        // factors its source supports.  It is held at --max_load_factor.
        double thrust_coefficient = state.CT;
        if (database.fixed_load_factor())
        {
            const double ceiling = options.max_load_factor * database.reference_thrust_coefficient();
            if (thrust_coefficient > ceiling)
            {
                thrust_coefficient = ceiling;
                load_capped++;
            }
        }
        source_states.push_back(database.source_state(Vec3(state.x, state.y, state.z),
                Vec3(nose[0], nose[1], nose[2]),
                Vec3(state.nx, state.ny, state.nz), state.mu,
                options.qsam ? state.fpa_e : state.fpa, thrust_coefficient));
        const auto& blend = source_states.back().blend;
        mu_clamped += blend.advance_ratio_clamped;
        fpa_clamped += blend.flight_path_angle_clamped;
        ct_clamped += blend.thrust_coefficient_clamped;
        load_scaled += blend.level_correction_db > 0.0;
    }
    log("Source axis coverage over ", nt, " states: advance ratio clamped ", mu_clamped,
            ", flight-path angle clamped ", fpa_clamped, ", thrust coefficient clamped ", ct_clamped,
            "; analytic load-factor correction applied ", load_scaled, ".");
    if (mu_clamped || fpa_clamped || ct_clamped)
        warn("Flight conditions exceed the source database's axis ranges: advance ratio ", mu_clamped,
                ", flight-path angle ", fpa_clamped, ", thrust coefficient ", ct_clamped,
                " of ", nt, " states. These queries use clamped source conditions.");
    if (load_capped > 0)
    {
        warn(load_capped, " of ", nt, " states exceed ", options.max_load_factor, " g on a "
                "fixed_load_factor database, where the load-factor correction is an extrapolation. "
                "Held at ", options.max_load_factor, " g (--max_load_factor).");
    }
    if (slow > 0)
    {
        if (by_heading > 0)
            info(by_heading, " of ", slow, " states under 8 kt turn the sphere's nose from the velocity "
                    "toward the heading (all of the way below 3 kt), turned into the local frame by ",
                    *offset * 180.0 / pi, " degrees.");
        else
            warn(slow, " states under 8 kt orient the sphere by their drift: ",
                    std::isnan(trajectory.states.front().heading)
                            ? "the track carries no heading."
                            : "too few fast states to relate heading to the frame, and no --frame_bearing.");
    }

    if (!options.power_file.empty())
    {
        // Shaft power is a property of the flight condition, so it falls out of
        // the per-state blends with no extra geometry.  It is written per
        // trajectory rather than mixed into the state file, whose columns are
        // pinned by the state golden.
        if (!database.has_shaft_power())
        {
            log("--powerfile given but this database carries no ShaftPower; skipping ",
                    job.trajectory, ".");
        }
        else
        {
            const auto name = job.count == 1
                    ? options.power_file
                    : (std::filesystem::path(options.power_file).parent_path()
                              / (std::filesystem::path(job.trajectory).stem().string() + "_power.csv"))
                              .string();
            auto file = open_output(name);
            info("Writing shaft power to ", name);
            file << "t(s),advance_ratio,flight_path_angle(deg),CT,shaft_power\n";
            for (std::size_t it = 0; it < nt; it++)
            {
                const auto& state = trajectory.states[it];
                file << state.t << "," << state.mu << "," << (options.qsam ? state.fpa_e : state.fpa) << ","
                     << state.CT << "," << database.shaft_power(source_states[it].blend) << "\n";
            }
        }
    }
    // The ambient level stays NaN unless -a was given, which is how D's default-NaN
    // double suppressed the ambient mix-in.  The energy it converts to is
    // loop-invariant; it used to be a pow() per evaluation.  It is mixed in as
    // energy (level_with_ambient), not as the pressure amplitude D summed.
    const bool add_ambient = !std::isnan(options.ambient_level);
    const double ambient_energy = add_ambient ? std::pow(10.0, 0.1 * options.ambient_level) : 0.0;

    // Observers are independent of one another: each reads the database and the
    // source_states array, which are finished and const by now, and writes only
    // its own level.  Nothing in here reaches the Delaunay triangulation -- the
    // blends were all resolved above -- so the loop splits with no locking.
    const unsigned workers = options.thread_count > 0 ? options.thread_count : default_thread_count();
    info("Evaluating ", observers.size(), " observers on ", workers,
            workers == 1 ? " thread." : " threads.");

    if (options.write_time_history)
    {
        // Create the output directory once, here.  Left to the writers it would
        // be a create_directories() race between workers, and a stat per file
        // besides.
        const auto directory = std::filesystem::path(options.time_history_prefix).parent_path();
        if (!directory.empty())
            std::filesystem::create_directories(directory);
    }
    const bool write_band_history = !options.band_history_prefix.empty();
    if (write_band_history)
    {
        const auto directory = std::filesystem::path(options.band_history_prefix).parent_path();
        if (!directory.empty())
            std::filesystem::create_directories(directory);
    }
    // The observer index restarts at 0 for every trajectory, so with several in
    // one run a bare prefix has each trajectory overwrite the files of the one
    // before it -- silently, with exit status 0, and only the last flight's
    // histories survive.  A run of one keeps the plain names; with more, each
    // file is named for its track as well (run() has checked the stems differ).
    const std::string history_infix = job.count == 1
            ? std::string()
            : "_" + std::filesystem::path(job.trajectory).stem().string();

    // Each observer's receiver, and where its level is evaluated: at the
    // microphone, the observer's point raised by the microphone's height above
    // the ground there.  A flush microphone (height 0) is on a plate when
    // --plate_table is given; any other is a point receiver at its height, for
    // the ground model.  An observer file giving heights (id,x,y,z,height) says
    // it per observer, so one run can hold both kinds -- most test arrays are
    // flush microphones with a few on poles.  Otherwise every observer takes
    // the run's: the plate, or --receiver_height with --ground_effect, or none.
    const std::size_t observer_count = observers.size();
    std::vector<Vec3> microphones(observer_count);
    std::vector<double> microphone_height(observer_count, 0.0);
    std::vector<char> on_plate(observer_count, 0);
    // Point receivers, by distinct height: each needs its own ground table.
    std::vector<double> point_heights;
    std::vector<std::size_t> point_of(observer_count, 0);
    {
        const double run_height = options.plate ? 0.0
                : ground_effect || options.propagation_table ? options.ground.receiver_height
                                                                     : 0.0;
        for (std::size_t i = 0; i < observer_count; i++)
        {
            const double height = microphone_heights.empty() ? run_height : microphone_heights[i];
            microphone_height[i] = height;
            microphones[i] = observers[i] + Vec3(0.0, 0.0, height);
            on_plate[i] = options.plate && height == 0.0 ? 1 : 0;
            if (ground_effect && !on_plate[i])
            {
                const auto found = std::find(point_heights.begin(), point_heights.end(), height);
                point_of[i] = static_cast<std::size_t>(found - point_heights.begin());
                if (found == point_heights.end())
                    point_heights.push_back(height);
            }
        }
    }
    const bool any_plate = std::any_of(on_plate.begin(), on_plate.end(), [](char c) { return c != 0; });
    const bool any_point = !point_heights.empty();
    if (ground_effect || !microphone_heights.empty())
    {
        std::vector<std::size_t> at_height(point_heights.size(), 0);
        std::size_t flush_count = 0;
        std::size_t bare_count = 0;
        for (std::size_t i = 0; i < observer_count; i++)
        {
            if (on_plate[i])
                flush_count++;
            else if (ground_effect)
                at_height[point_of[i]]++;
            else
                bare_count++;
        }
        std::ostringstream summary;
        auto item = [&summary](const auto&... parts) {
            if (summary.tellp() > 0)
                summary << ", ";
            (summary << ... << parts);
        };
        if (flush_count > 0)
            item(flush_count, " flush on plates");
        for (std::size_t k = 0; k < point_heights.size(); k++)
            item(at_height[k], " at ", point_heights[k] / options.unit.to_meters, " ", options.unit.name);
        if (bare_count > 0)
            item(bare_count, " with no ground model");
        log("Receivers: ", summary.str(), "; each level is evaluated at its microphone.");
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
    // that margin is what makes the result exact where it matters: for any
    // observer whose peak reaches the cutoff, the window's ends and everything
    // near the peak are at or above cutoff - 10 and are never culled.  The one
    // exception is a dip inside the window (the 10 dB-down duration runs from
    // the first to the last sample within 10 dB of the peak, dips included)
    // deeper than cutoff - 10: it can be culled, and it loses what is at most a
    // tenth of the peak's energy per sample.  On a 101x101 grid over an EC130
    // approach, 11 of the 404 observers reaching a 60 dB cutoff had such dips,
    // and none changed.  Observers that never reach the cutoff are the ones
    // whose level changes, and they are precisely the ones the cutoff says are
    // not of interest.

    // Spectral propagation setup.  Both of these depend only on the atmosphere
    // and the database's band centers, so they are computed once here rather
    // than per observer: the coefficients cost a pow() and two exp() per band.
    std::vector<double> absorption;
    std::vector<double> band_gains;
    // With --atmosphere, absorption integrated through the profile instead
    // (LayeredAbsorption), and its running integral at each state's height,
    // [state][band], so the observer loop never looks a state's up twice.
    std::optional<LayeredAbsorption> layered;
    std::vector<double> state_absorption_integral;
    double z_low = 0.0;
    double z_high = 0.0;
    // When every observer stands at one height -- any grid -- a state's mean
    // coefficient is the same for all of them, so it is formed here once,
    // [state][band], and the loop reads it instead of forming it per pair.
    std::vector<double> state_absorption_mean;
    const GroundModel* ground_model = nullptr;
    std::vector<std::complex<double>> ground_impedance_bands;
    std::optional<PlateReceiver> plate;
    // Per distinct point-receiver height: its ground model (receiver_height
    // set) and, unless --exact_ground, its table.
    std::vector<GroundModel> point_models;
    std::vector<GroundTable> point_tables;
    // The contiguous run of bands inside the metric's limits. Band centers are
    // sorted, so the limits always select a range rather than a scattered set.
    std::size_t first_band = 0;
    std::size_t last_band = 0;
    // --propagation_table: where the run's bands sit in it, and, for a table
    // that already includes absorption, the zero coefficients that keep the
    // run from taking it twice.
    PropagationTable::BandMap table_bands;
    std::vector<double> table_absorption;
    // [azimuth][band]: K for plate microphones read from the table.
    constexpr int table_plate_azimuths = 36;
    std::vector<std::vector<double>> table_plate_offset;

    // The atmosphere is needed on every path now, not only the spectral one:
    // its sound speed turns the emission-time level history into reception time
    // for the exposure integral (see reception_factor in sel.hpp).
    Atmosphere air = trajectory.vehicle.atmosphere;
    if (!std::isnan(options.temperature_override))
        air.temperature = options.temperature_override;
    if (!std::isnan(options.pressure_override))
        air.pressure = options.pressure_override;
    if (!std::isnan(options.humidity_override))
        air.relative_humidity = options.humidity_override;
    // With a profile, its level at the ground stands in for the one
    // atmosphere everything but absorption still takes: the sound speed of the
    // reception-time integral, and the ground impedance.  Both are properties
    // of the air at the surface.  A profile without pressure takes the surface
    // pressure the vehicle file (or --pressure) gives, as `air` now holds.
    std::optional<AtmosphereProfile> profile;
    if (atmosphere)
    {
        profile = atmosphere->has_pressure() ? *atmosphere
                                             : atmosphere->with_hydrostatic_pressure(air.pressure);
        air = profile->at(0.0);
    }
    const double sound_speed = air.sound_speed();

    // Spectral propagation is the default, since per-band ISO 9613-1 absorption
    // is what the EAA factor approximates -- and approximates badly: EAA is one
    // A-weighted attenuation fitted at 1000 m and scaled linearly with range,
    // which under-attenuates at the few hundred meters a footprint mostly sees.
    // Against the 2017 approach measurements (139 runs, ~6200 mic-runs) the
    // broadband path read +1.2 dB high in the median (rms 3.1 dB) and the
    // spectral one -0.2 dB (rms 2.5).  A database without spectra falls back to
    // broadband unless --spectral was asked for by name.
    bool spectral = options.spectral;
    if (spectral && !database.has_spectra())
    {
        if (options.spectral_explicit || options.ground_effect_given || atmosphere
                || options.propagation_table)
        {
            const std::string profile_option = options.surface_wind ? "--wind" : "--atmosphere";
            const std::string asked = options.propagation_table ? "--propagation_table"
                    : options.ground_effect_given               ? "--ground_effect"
                    : atmosphere                                ? profile_option
                                                                : "--spectral";
            input_error(asked
                    + " needs a database carrying per-band spectra; this one has none. Rebuild "
                      "it with spectra, or propagate the broadband level with --no-spectral"
                    + (options.ground_effect_given ? std::string(" and no ground effect.")
                                    : atmosphere ? " and no " + profile_option + "."
                                                 : std::string(".")));
        }
        if (!options.band_history_prefix.empty())
        {
            input_error("--band_history needs a database carrying per-band spectra; this one "
                        "has none.");
        }
        log("Database carries no per-band spectra: propagating the broadband level with its "
            "EAA factor.");
        spectral = false;
    }
    {
        // What this trajectory is propagated with, on stdout, so a footprint
        // states its own configuration rather than leaving it in the log file.
        const char* humidity_from = options.atmosphere ? "the --atmosphere profile's ground level"
                : !std::isnan(options.humidity_override) ? "--humidity"
                : trajectory.vehicle.humidity_given ? "vehicle file" : "default, not given";
        const double receiver = options.ground.receiver_height / options.unit.to_meters;
        std::cout << "Physics: "
                  << (spectral ? "spectral, per-band ISO 9613-1 absorption" : "broadband, database EAA factor")
                  << "; ground "
                  << (ground_effect ? "on, receiver " + std::to_string(receiver) + " " + options.unit.name + " up"
                                    : std::string("off"))
                  << "; air " << air.temperature << " K, " << air.pressure << " kPa, "
                  << air.relative_humidity << " % RH (" << humidity_from << ")"
                  << (profile ? (options.surface_wind ? ", layered profile from --wind" : ", layered profile")
                              : "")
                  << (options.refraction ? "; refraction on" : "") << ".\n";
        // The broadband path carries the EAA factor the database was built
        // with, so the run's air only enters through the spreading and the
        // speed of sound.  A different build atmosphere is an absorption error
        // no later step can see, so say so where the run is stated.
        if (!spectral && database.build_atmosphere())
        {
            const auto& build = *database.build_atmosphere();
            if (std::fabs(air.temperature - build.temperature_kelvin) > 0.5
                    || std::fabs(air.pressure - build.pressure_kpa) > 0.1
                    || std::fabs(air.relative_humidity - build.relative_humidity_percent) > 1.0)
            {
                warn("The broadband EAA factor was built at ", build.temperature_kelvin, " K, ",
                        build.pressure_kpa, " kPa, ", build.relative_humidity_percent,
                        " % RH, and this run's air differs. Its absorption is the build's, not "
                        "this run's: use --spectral for the run's own air, or rebuild the "
                        "database at this atmosphere.");
            }
        }
    }
    if (spectral)
    {

        const auto& bands = database.band_centers();
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

        // Where each value came from, so a default is not mistaken for a
        // measurement: the vehicle file usually gives a temperature and a
        // density, derives the pressure from them, and gives no humidity at all.
        const auto origin = [](bool overridden, const char* option, bool given, const char* otherwise) {
            return std::string(overridden ? option : given ? "vehicle file" : otherwise);
        };
        const auto& vehicle = trajectory.vehicle;
        const std::string temperature_from = origin(!std::isnan(options.temperature_override),
                "--temperature", vehicle.temperature_given, "default");
        const std::string pressure_from = origin(!std::isnan(options.pressure_override), "--pressure",
                vehicle.pressure_given, "from the vehicle file's density and temperature");
        const std::string humidity_from = origin(!std::isnan(options.humidity_override), "--humidity",
                vehicle.humidity_given, "default");
        log("Propagating spectra: ", retained, " of ", bands.size(), " bands within ",
                options.band_limits.low, " to ", options.band_limits.high, " Hz, at ",
                air.temperature, " K (", temperature_from, "), ", air.pressure, " kPa (",
                pressure_from, "), ", air.relative_humidity, "% RH (", humidity_from, ")",
                profile ? ", the ground values" : "", ".");
        // Refused rather than assumed: the spectral path attenuates every band
        // at this humidity, and going from 20 to 50 percent moved a footprint
        // 2 dB, so a silent 20 percent would be a prediction nobody chose.
        if (!profile && humidity_from == "default")
        {
            input_error("Relative humidity is not set (no --humidity, and ", job.vehicle,
                    " gives none), and the spectral path needs it: absorption depends strongly on "
                    "it. Give --humidity, add \"relative humidity\" to the vehicle file or its "
                    "Atmosphere section, or pass --atmosphere.");
        }
        // Say what that atmosphere costs at the extremes actually in play, so a
        // wrong humidity is visible in the log rather than only in the answer.
        log("  absorption spans ", absorption.front(), " dB/m at ", bands.front(), " Hz to ",
                absorption.back(), " dB/m at ", bands.back(), " Hz", profile ? " at the ground" : "",
                ".");

        if (profile)
        {
            // The heights any path spans: every state's and every microphone's,
            // and the ground.
            for (const auto& state : source_states)
            {
                z_low = std::min(z_low, state.position.z);
                z_high = std::max(z_high, state.position.z);
            }
            for (const auto& microphone : microphones)
            {
                z_low = std::min(z_low, microphone.z);
                z_high = std::max(z_high, microphone.z);
            }
            // One level is a uniform atmosphere and reaches every height.  Past
            // the top of several, the profile holds its last level, which is a
            // guess about air nobody measured -- refuse rather than guess.
            if (profile->height.size() > 1 && z_high > profile->top())
            {
                input_error("--atmosphere ", profile->source, " reaches ",
                        profile->top() / options.unit.to_meters, " ", options.unit.name,
                        " above the observers, but paths from ", job.trajectory, " reach ",
                        z_high / options.unit.to_meters, " ", options.unit.name,
                        ". Extend the profile above the highest state and observer.");
            }
            if (profile->height.size() > 1 && z_low < profile->bottom())
            {
                log("  the profile's lowest level, ", profile->bottom() / options.unit.to_meters,
                        " ", options.unit.name, ", is held down to ", z_low / options.unit.to_meters,
                        " ", options.unit.name, ".");
            }
            layered = LayeredAbsorption(*profile, bands, first_band, last_band, z_low, z_high);
            state_absorption_integral.assign(nt * bands.size(), 0.0);
            for (std::size_t it = 0; it < nt; it++)
            {
                layered->integral(source_states[it].position.z,
                        state_absorption_integral.data() + it * bands.size());
            }
            const bool common_height = std::all_of(microphones.begin(), microphones.end(),
                    [&](const Vec3& m) { return m.z == microphones.front().z; });
            if (common_height && !microphones.empty())
            {
                std::vector<double> at_observers(bands.size(), 0.0);
                const double z_observers = microphones.front().z;
                layered->integral(z_observers, at_observers.data());
                state_absorption_mean.assign(nt * bands.size(), 0.0);
                for (std::size_t it = 0; it < nt; it++)
                {
                    layered->mean_coefficients(source_states[it].position.z,
                            state_absorption_integral.data() + it * bands.size(), z_observers,
                            at_observers.data(), state_absorption_mean.data() + it * bands.size());
                }
            }
            std::vector<double> at_top(bands.size(), 0.0);
            layered->coefficients(z_high, at_top.data());
            log("  absorption layered through ", profile->source, " from ",
                    z_low / options.unit.to_meters, " to ", z_high / options.unit.to_meters, " ",
                    options.unit.name, ": at ", bands[last_band - 1], " Hz it is ",
                    absorption[last_band - 1], " dB/m at the ground and ", at_top[last_band - 1],
                    " dB/m at the top.");
        }

        if (ground_effect)
        {
            ground_model = &options.ground;
            ground_impedance_bands = ground_impedance(options.ground, bands, sound_speed);
            log("  ground effect on: impedance model ", static_cast<int>(options.ground.impedance_model),
                    ", flow resistance ", options.ground.flow_resistance,
                    " kPa s/m^2, receiver ", options.ground.receiver_height, " m up, turbulence ",
                    options.ground.turbulence_model == TurbulenceModel::chessell
                            ? options.ground.turbulence : options.ground.gamma_t,
                    options.ground.turbulence_model == TurbulenceModel::chessell ? " (Chessell)"
                                                                                 : " (HARMONOISE gamma_T)",
                    ", roughness ", options.ground.roughness, " m, sound speed ", sound_speed, " m/s.");
            // Farthest any observer is from any state: how far the ground
            // tables must reach.
            auto reach = [&]() {
                double lo_x = 1e300, hi_x = -1e300, lo_y = 1e300, hi_y = -1e300;
                for (const auto& o : observers)
                {
                    lo_x = std::min(lo_x, o.x); hi_x = std::max(hi_x, o.x);
                    lo_y = std::min(lo_y, o.y); hi_y = std::max(hi_y, o.y);
                }
                double r_max = 1.0;
                for (const auto& state : source_states)
                    for (const double cx : {lo_x, hi_x})
                        for (const double cy : {lo_y, hi_y})
                            r_max = std::max(r_max, std::hypot(std::hypot(cx - state.position.x,
                                    cy - state.position.y), state.position.z));
                return r_max * 1.01;
            };
            for (const double height : point_heights)
            {
                GroundModel model = options.ground;
                model.receiver_height = height;
                point_models.push_back(model);
                if (!options.exact_ground)
                {
                    const auto table_start = std::chrono::steady_clock::now();
                    const double r_max = reach();
                    point_tables.emplace_back(model, ground_impedance_bands, sound_speed, bands,
                            first_band, last_band, 1.0, r_max);
                    log("  ground effect for microphones ", height / options.unit.to_meters, " ",
                            options.unit.name, " up tabulated to ", r_max, " m in ",
                            std::chrono::duration<double>(std::chrono::steady_clock::now() - table_start)
                                    .count(),
                            " s (--exact_ground to evaluate it directly).");
                }
            }
            if (any_plate)
            {
                plate = options.plate->bound(bands, options.ground, sound_speed, first_band,
                        last_band);
                if (!options.exact_ground)
                {
                    const auto table_start = std::chrono::steady_clock::now();
                    const double r_max = reach();
                    plate->tabulate(r_max, 161, 40);
                    log("  plate reflection coefficients tabulated to ", r_max, " m in ",
                            std::chrono::duration<double>(std::chrono::steady_clock::now() - table_start)
                                    .count(),
                            " s (--exact_ground to evaluate them directly).");
                }
                log("  observers are flush microphones on a ground plate: ", plate->description(),
                        options.plate_orientation.outboard
                                ? std::string(", microphone outboard of the track.")
                                : ", microphone offset toward "
                                        + std::to_string(options.plate_orientation.azimuth * 180.0 / pi)
                                        + " deg.");
            }
        }

        if (options.propagation_table)
        {
            const PropagationTable& table = *options.propagation_table;
            try
            {
                table_bands = table.match_bands(bands, first_band, last_band);
            }
            catch (const std::runtime_error& e)
            {
                input_error("--propagation_table ", options.propagation_table_path, ": ", e.what(), ".");
            }
            // Receivers it does not reach are refused rather than handed to
            // the fallback wholesale: a table that misses the run's microphones
            // was built for another run.
            std::vector<double> uncovered;
            for (const double height : microphone_height)
                if (!table.covers_receiver(height)
                        && std::find(uncovered.begin(), uncovered.end(), height) == uncovered.end())
                    uncovered.push_back(height);
            if (!uncovered.empty())
            {
                std::ostringstream heights;
                for (std::size_t k = 0; k < uncovered.size(); k++)
                    heights << (k ? ", " : "") << uncovered[k] / options.unit.to_meters;
                input_error("--propagation_table ", options.propagation_table_path, " holds receivers from ",
                        table.axes().receiver_heights.front() / options.unit.to_meters, " to ",
                        table.axes().receiver_heights.back() / options.unit.to_meters, " ", options.unit.name,
                        "; this run's microphones are also at ", heights.str(), " ", options.unit.name,
                        ". Build the table with those receiver heights.");
            }
            if (table.includes_absorption())
                table_absorption.assign(bands.size(), 0.0);
            // A flush microphone on a plate reads the table at its lowest
            // receiver height, plus K, the plate over bare ground there, per band
            // and azimuth -- as ShadowTransition places a plate in the shadow:
            // near grazing K barely depends on elevation, so it is taken at 2
            // degrees, a kilometer out.
            if (plate)
            {
                GroundModel bare = options.ground;
                bare.receiver_height = table.axes().receiver_heights.front();
                const double elevation = 2.0 * pi / 180.0;
                const double h = 1000.0 * std::sin(elevation), d = 1000.0 * std::cos(elevation);
                std::vector<double> bare_db(bands.size(), 0.0), plate_db(bands.size(), 0.0);
                excess_ground_attenuation(bare, ground_impedance_bands, h, d, sound_speed, bands,
                        bare_db.data(), first_band, last_band, false);
                table_plate_offset.assign(table_plate_azimuths, std::vector<double>(bands.size(), 0.0));
                for (int i = 0; i < table_plate_azimuths; i++)
                {
                    plate->level(h, d, 2.0 * pi * i / table_plate_azimuths, plate_db.data(), first_band,
                            last_band);
                    for (std::size_t band = first_band; band < last_band; band++)
                        table_plate_offset[i][band] = plate_db[band] - bare_db[band];
                }
            }
            log("  propagation table ", options.propagation_table_path, ": ", table.method(), "; ",
                    table.axes().bearings.size(), " bearings, source heights ",
                    table.axes().source_heights.front(), " to ", table.axes().source_heights.back(),
                    " m, ranges to ", table.axes().ranges.back(), " m; absorption ",
                    table.includes_absorption() ? "included in it" : "applied by the run", ". Pairs it "
                    "does not cover take ", ground_model ? "the ground model" : "free-field propagation", ".");
            if (!table.isotropic() && std::isnan(options.frame_bearing))
            {
                // A table that depends on bearing, read in a frame that is not
                // east/north, is mis-rotated: up to 3.2 dB on the demo track with
                // only an info line to say so.  A warning for an AC run, whose x/y
                // are rotated; a note where the frame is a geodetic one.
                const char* note = "  the table's bearings are compass bearings, and this run takes +x as "
                                   "east. If its x/y are rotated -- the shipped AC tracks run +x along "
                                   "210 degrees -- give --frame_bearing.";
                if (options.geodetic.frame)
                    log(note);
                else
                    warn(note);
            }
        }
    }

    // --refraction linear: the effective sound speed tabulated for the
    // per-path linear fits, and each state's moments.  Rays end at the
    // microphones, like everything else.
    std::optional<EffectiveSoundSpeed> effective_sound_speed;
    std::vector<EffectiveSoundSpeed::Moments> state_moments;
    ShadowAttenuation shadow_attenuation;
    // Point receivers near a shadow boundary: one table per microphone height;
    // flush plate microphones, one more.
    std::vector<ShadowTransition> shadow_transitions;
    std::optional<ShadowTransition> plate_transition;
    // --refraction stratified: the profile's layers, with the linear model kept
    // alongside for its warm start and for what the layered solve leaves to it.
    std::optional<LayeredSoundSpeed> layered_sound_speed;
    if (options.refraction && spectral && profile)
    {
        const double bearing = std::isnan(options.frame_bearing) ? 90.0 : options.frame_bearing;
        effective_sound_speed = EffectiveSoundSpeed(*profile, bearing, z_low,
                z_high + EffectiveSoundSpeed::min_fit_span);
        if (options.stratified)
            layered_sound_speed = LayeredSoundSpeed(*profile, bearing);
        state_moments.reserve(nt);
        for (const auto& state : source_states)
            state_moments.push_back(effective_sound_speed->moments(state.position.z));
        // The creeping wave runs along the ground: over its admittance when
        // there is a ground model, over rigid ground (the old rate) when not.
        shadow_attenuation = ground_model
                ? ShadowAttenuation(database.band_centers(), ground_impedance_bands, first_band, last_band)
                : ShadowAttenuation(database.band_centers(), first_band, last_band);
        if (ground_model != nullptr)
        {
            for (const double height : point_heights)
            {
                const auto start = std::chrono::steady_clock::now();
                shadow_transitions.emplace_back(database.band_centers(), ground_impedance_bands,
                        height, sound_speed, first_band, last_band, nullptr, nullptr,
                        options.thread_count);
                if (options.scattering > 0.0)
                    shadow_transitions.back().set_scattering(
                            ShadowScattering::for_receiver(options.scattering, height, false));
                log("  shadow-boundary field for microphones ", height / options.unit.to_meters, " ",
                        options.unit.name, " up tabulated from the residue series in ",
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(),
                        " s.");
            }
            if (any_plate)
            {
                const auto start = std::chrono::steady_clock::now();
                plate_transition.emplace(database.band_centers(), ground_impedance_bands,
                        plate->thickness(), sound_speed, first_band, last_band, &*plate, ground_model,
                        options.thread_count);
                if (options.scattering > 0.0)
                    plate_transition->set_scattering(
                            ShadowScattering::for_receiver(options.scattering, plate->thickness(), true));
                log("  shadow-boundary field for plate microphones tabulated in ",
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(),
                        " s.");
            }
        }
        std::string surface_bands;
        for (std::size_t band = first_band; band < last_band; band++)
            if (shadow_attenuation.admits_surface_wave(band))
                surface_bands += (surface_bands.empty() ? "" : ", ") + std::to_string(
                        static_cast<int>(std::lround(database.band_centers()[band]))) + " Hz";
        if (!surface_bands.empty())
            log("  warning: the ground's admittance allows a surface wave in the shadow at ",
                    surface_bands, ", which the creeping-wave loss leaves out. Those bands take "
                    "the rigid ground's rate, and their shadow levels are likely too low.");
        // What the fit makes of the air along and against +x from the
        // highest state, so a profile read wrong shows in the log.
        const auto top = std::max_element(source_states.begin(), source_states.end(),
                [](const auto& a, const auto& b) { return a.position.z < b.position.z; });
        const auto at_top = effective_sound_speed->moments(top->position.z);
        const auto at_ground = effective_sound_speed->moments(0.0);
        const auto along = effective_sound_speed->fit(top->position.z, at_top, 0.0, at_ground,
                1.0, 0.0);
        const auto against = effective_sound_speed->fit(top->position.z, at_top, 0.0, at_ground,
                -1.0, 0.0);
        log(options.stratified
                        ? "Refraction: curved rays through the profile's layers, the effective sound "
                          "speed linear between its levels. The linear fit per path"
                        : "Refraction: curved rays through a linear effective sound speed fitted "
                          "per path.",
                " From the highest state (",
                top->position.z / options.unit.to_meters, " ", options.unit.name,
                ") the gradient is ", along.gradient, " /s toward +x and ", against.gradient,
                " /s toward -x; shadow zones attenuate by the first creeping-wave mode over ",
                shadow_attenuation.rigid() ? "rigid ground" : "the ground's impedance",
                options.scattering > 0.0 ? ", filled by turbulent scattering at gamma_T " : ", capped at ",
                options.scattering > 0.0 ? options.scattering : options.shadow_floor,
                options.scattering > 0.0 ? "." : " dB.");
        if (profile->has_wind() && std::isnan(options.frame_bearing))
        {
            log("  the profile's wind directions are compass bearings, turned into the local "
                "frame taking +x as east. If this run's x/y are rotated -- the shipped AC tracks "
                "run +x along 210 degrees -- give --frame_bearing.");
        }
    }
    else if (profile && profile->has_wind() && !air_referenced)
    {
        log("  the profile's wind is read but not used: it matters only to --refraction and to "
            "a database referred to the air (--speed_reference).");
    }

    // The spectral path's bound (Database::spectral_direction_bounds) serves
    // both --db_cutoff and the peak search below.  With --ground_effect, each
    // band's gain is raised by the most the ground can add to it
    // (ground_gain_bounds), which keeps it a bound with the ground included.
    const bool culling = !std::isnan(options.db_cutoff);
    const bool level_is_transformed = add_ambient;
    if (options.want_coverage)
    {
        if (!database.has_coverage())
            input_error("--coverage_file needs a database that records its coverage, and ",
                    job.trajectory, " is run against one that does not.");
        if (options.refraction || level_is_transformed)
            input_error("--coverage_file reports the straight-ray levels: it does not combine with "
                    "--refraction or -a.");
    }
    const bool needs_full_history = options.write_time_history || options.want_dnl
            || options.full_history || options.want_lamax || options.want_nax
            || write_band_history;
    // The coverage report needs every state's level, so nothing is skipped.
    const bool bounded = !level_is_transformed && !options.want_coverage;
    DirectionBounds direction_bounds;
    // Only the spectral path has a bound worth building.  A broadband evaluation
    // costs about 20 ns, and the bound test costs more than the evaluations it
    // saves: under the outermost-crossing window the search evaluates 94 percent
    // of the pairs on the demo, and ran 4.6 percent slower than marching every
    // state (output identical).  A spectral evaluation is 150 ns, a ground-effect
    // one a microsecond, and there the same search pays.  Broadband --db_cutoff
    // culls on a per-state range (cull_range) and needs no direction table.
    if (bounded && spectral && (culling || !needs_full_history))
    {
        // Both report themselves unusable on an empty trajectory, so nt > 0
        // does not need testing here as well.
        const auto bounds_start = std::chrono::steady_clock::now();
        if (spectral)
        {
            std::vector<double> bound_gains = band_gains;
            // What the path can add above spherical spreading, a power ratio
            // per band.
            std::vector<double> path_gains(bound_gains.size(), 1.0);
            if (ground_model != nullptr)
            {
                // The loudest the ground can make any observer: a plate's, a
                // point receiver's, or the larger of the two when a run holds
                // both.  The point receiver's does not depend on its height.
                std::vector<double> ground_gains(bound_gains.size(), 0.0);
                if (any_plate)
                {
                    ground_gains = plate->gain_bounds(first_band, last_band);
                    if (plate_transition)
                        for (std::size_t band = 0; band < ground_gains.size(); band++)
                            ground_gains[band] = std::max(ground_gains[band],
                                    plate_transition->max_power_ratio(band));
                }
                if (any_point)
                {
                    const auto point_gains = ground_gain_bounds(*ground_model,
                            ground_impedance_bands, first_band, last_band);
                    for (std::size_t band = 0; band < ground_gains.size(); band++)
                        ground_gains[band] = std::max(ground_gains[band], point_gains[band]);
                    // Near a shadow boundary the transition's field stands in
                    // for the ground term, and the bound has to cover it too.
                    for (const auto& transition : shadow_transitions)
                        for (std::size_t band = 0; band < ground_gains.size(); band++)
                            ground_gains[band] = std::max(ground_gains[band],
                                    transition.max_power_ratio(band));
                }
                path_gains = ground_gains;
            }
            // A pair the table covers can gain as much as its largest value;
            // one it does not, as much as the ground model above.
            if (options.propagation_table)
            {
                const auto table_gains = options.propagation_table->max_power_ratio(table_bands);
                for (std::size_t band = 0; band < path_gains.size(); band++)
                    path_gains[band] = std::max(path_gains[band], table_gains[band]);
            }
            for (std::size_t band = 0; band < bound_gains.size(); band++)
                bound_gains[band] *= path_gains[band];
            direction_bounds = database.spectral_direction_bounds(source_states, bound_gains);
        }
        log("Direction bounds built in ",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - bounds_start)
                        .count(),
                " s.");
    }

    const double cull_below = options.db_cutoff - sel_window_db;
    // A culled sample contributes no acoustic energy.  Filling with the
    // threshold instead would be far worse than useless -- thousands of samples
    // at the floor integrate to well above it.
    const double culled_level = -std::numeric_limits<double>::infinity();
    // Broadband culls on a per-state range (cull_range, from the loudest
    // sphere the state blends); spectral culls per direction, on the table
    // above, which is the tighter test.
    const bool cull_by_direction = culling && spectral;
    if (cull_by_direction && !direction_bounds.usable())
    {
        input_error("--db_cutoff: no bound on this run's spectral level (no band inside "
                    "--band_limits?), so nothing could be culled safely.");
    }
    const double cull_energy = std::pow(10.0, 0.1 * cull_below);
    const double sphere_radius_squared = database.sphere_radius() * database.sphere_radius();
    std::vector<double> cull_range_squared(culling && !cull_by_direction ? nt : 0);
    for (std::size_t it = 0; it < cull_range_squared.size(); it++)
    {
        const double range = database.cull_range(source_states[it].blend, cull_below);
        cull_range_squared[it] = range * range;
    }
    auto is_culled = [&](std::size_t it, const Vec3& observer) {
        const Vec3 delta = observer - source_states[it].position;
        const double squared_range = delta.squaredNorm();
        if (!cull_by_direction)
            return squared_range > cull_range_squared[it];
        const auto& frame = source_states[it].frame;
        return squared_range > sphere_radius_squared
                && direction_bounds.numerator(it,
                           DirectionBounds::bin(frame.nose.dot(delta), frame.normal.dot(delta),
                                   frame.starboard.dot(delta), squared_range))
                        < cull_energy * squared_range;
    };
    std::atomic<std::size_t> culled_samples {0};

    // Find each observer's peak and expand out to the 10-dB-down points, rather
    // than marching every emission time -- see windowed_sound_exposure_level.
    // It is exact, so it is on whenever it applies.  Everything else marches,
    // for one of two reasons:
    //
    //   there is no bound on the level being integrated:
    //     -a / --ambient    the level is transformed after eval
    //
    //   every sample is wanted, not just the window's:
    //     -i / --history writes the whole time history, skipped states included
    //     --band_history the same, per band
    //     -D / --dnlfile the day-night window is found on a penalized history,
    //                    which a window found on the raw one does not cover
    //
    // A new option that does either belongs in the matching term above.
    // --db_cutoff stands apart: its culling is a deliberate metric change, so
    // leave the two mechanisms separate rather than layering this one on top.
    // Under --refraction the same bounds hold, taken along each ray's launch
    // direction at its tube's range -- see windowed_exposure.
    // A table beside refraction mixes two geometries per state -- the
    // straight line the table is read along and the ray's launch direction --
    // which the bounds were not built for, so every state is evaluated.
    const bool expand_from_peak = direction_bounds.usable() && !culling && !needs_full_history
            && !(options.propagation_table && options.refraction);

    // Each blend's spectra packed contiguously (Database::pack_spectra), one
    // copy per distinct blend the track visits, so the observer loop's reads
    // stay within one corner's run of bands.  Measured on an approach track
    // (3,520 states, 14,641 observers): 14% faster, and identical to the bit.
    std::vector<Database::PackedSpectra> packed_spectra;
    std::vector<const Database::PackedSpectra*> packed_of_state(spectral ? nt : 0, nullptr);
    if (spectral)
    {
        std::vector<std::array<std::size_t, condition_dims + 1>> packed_indices;
        std::size_t last = 0;
        for (std::size_t it = 0; it < nt; it++)
        {
            const auto& indices = source_states[it].blend.indices;
            std::size_t found = packed_indices.size();
            if (last < packed_indices.size() && packed_indices[last] == indices)
                found = last;
            for (std::size_t b = 0; found == packed_indices.size() && b < packed_indices.size(); b++)
            {
                if (packed_indices[b] == indices)
                    found = b;
            }
            if (found == packed_indices.size())
            {
                packed_indices.push_back(indices);
                packed_spectra.push_back(database.pack_spectra(indices, first_band, last_band));
            }
            last = found;
        }
        // Pointers only once the outer vector has stopped growing.
        for (std::size_t it = 0; it < nt; it++)
        {
            const auto& indices = source_states[it].blend.indices;
            const auto where = std::find(packed_indices.begin(), packed_indices.end(), indices);
            packed_of_state[it] = &packed_spectra[where - packed_indices.begin()];
        }
        log("Packed spectra for ", packed_spectra.size(), " blends (",
                static_cast<double>(packed_spectra.size() * packed_spectra.front().values.size())
                        * sizeof(double) / 1048576.0,
                " MiB).");
    }
    std::atomic<std::size_t> evaluated_samples {0};
    std::atomic<std::size_t> shadowed_samples {0};
    std::atomic<std::size_t> table_samples {0};
    std::atomic<std::size_t> fallback_samples {0};
    // The table's bearings are compass bearings; the run's x/y are turned from
    // east/north by --frame_bearing (CompassToLocal), whose inverse this is.
    const CompassToLocal compass(std::isnan(options.frame_bearing) ? 90.0 : options.frame_bearing);

    const auto acoustics_start = std::chrono::steady_clock::now();

    parallel_chunks(observers.size(), workers, [&](std::size_t begin, std::size_t end) {
        // One scratch history per chunk rather than per observer; parallel.hpp
        // deals out several chunks per worker, so this runs a few times a
        // thread, which it measures as noise against the work in the chunk.
        std::vector<double> levels(nt);
        std::vector<double> ranges(nt);
        std::vector<double> penalized;
        // --band_history: every state's band levels, and the geometry each
        // was reached along.
        const std::size_t history_bands = write_band_history ? database.band_centers().size() : 0;
        std::vector<double> band_levels(nt * history_bands);
        std::vector<BandHistoryRow> band_rows(write_band_history ? nt : 0);
        auto record_path = [&](std::size_t it, const RayPath& path) {
            if (!write_band_history)
                return static_cast<double*>(nullptr);
            auto& row = band_rows[it];
            const auto angles = lookup_angles(source_states[it].frame, path.direction,
                    !database.upper_hemisphere_measured());
            row.range = path.range;
            row.elevation = std::asin(std::clamp(-path.direction.z, -1.0, 1.0)) * 180.0 / pi;
            row.phi = angles[0];
            row.theta = angles[1];
            return band_levels.data() + it * history_bands;
        };
        // Which states of `levels` this observer has filled, stamped by a
        // counter so the array does not have to be cleared per observer.
        std::vector<unsigned> visited(expand_from_peak ? nt : 0, 0u);
        // What the peak search's main pass saw of each state it skipped (see
        // windowed_exposure); valid for a state exactly when it is not visited.
        std::vector<double> skipped_numerator(expand_from_peak ? nt : 0);
        std::vector<double> skipped_squared_range(expand_from_peak ? nt : 0);
        // Refraction through the peak search: each state's reception range,
        // kept from when its ray was formed, and stamped like `visited`.
        const bool windowed_refraction = expand_from_peak && effective_sound_speed.has_value();
        std::vector<double> travel(windowed_refraction ? nt : 0);
        std::vector<unsigned> travel_stamp(windowed_refraction ? nt : 0, 0u);
        // A refracted ray's geometry for the bounds: its launch direction
        // scaled by its tube's range, in the sphere frame, which is what the
        // level spreads over and looks the sphere up along.  Absorption, the
        // shadow's loss and the ground's gain are bounded apart from it.
        auto equivalent_geometry = [&](std::size_t it, const RayPath& path) {
            const Vec3 d = path.direction * path.range;
            const auto& frame = source_states[it].frame;
            return StateGeometry {path.range * path.range, frame.nose.dot(d), frame.normal.dot(d),
                    frame.starboard.dot(d)};
        };
        // With more than one eigenray (Arrival), a geometry whose bound covers
        // their sum: each is bounded by its own direction's numerator over its
        // own squared range, so the largest of those numerators over the
        // harmonic sum of the squared ranges bounds them all.
        auto bound_geometry = [&](std::size_t it, const RefractedRay& ray) {
            const StateGeometry primary = equivalent_geometry(it, ray.path);
            if (ray.extra.empty())
                return primary;
            double inverse = 1.0 / primary.squared_range;
            StateGeometry pick = primary;
            double best = direction_bounds.numerator(it, DirectionBounds::bin(primary.along, primary.up,
                    primary.across, primary.squared_range));
            for (const Arrival& a : ray.extra)
            {
                const StateGeometry g = equivalent_geometry(it, a.path);
                inverse += 1.0 / g.squared_range;
                const double n = direction_bounds.numerator(it,
                        DirectionBounds::bin(g.along, g.up, g.across, g.squared_range));
                if (n > best)
                {
                    best = n;
                    pick = g;
                }
            }
            const double scale = std::sqrt(1.0 / (inverse * pick.squared_range));
            return StateGeometry {1.0 / inverse, pick.along * scale, pick.up * scale, pick.across * scale};
        };
        struct RefractedGeometry
        {
            double squared_range, along, up, across;
            RefractedRay ray;
        };
        // The two halves of a refracted ray, from whichever model the run
        // asked for; the layered one keeps its path's layers in `layer_scratch`.
        std::vector<double> layer_scratch;
        // `transition`: the observer is a point receiver that ShadowTransition
        // serves, so past the layered boundary its lit level is the linear fit's.
        auto form_ray = [&](std::size_t it, const Vec3& observer,
                                const EffectiveSoundSpeed::Moments& at_receiver, bool transition) {
            if (layered_sound_speed)
            {
                return layered_direct_ray(*layered_sound_speed, *effective_sound_speed,
                        source_states[it].position, state_moments[it], observer, at_receiver,
                        layer_scratch, transition, options.waves);
            }
            return direct_ray(*effective_sound_speed, source_states[it].position, state_moments[it],
                    observer, at_receiver);
        };
        auto finish_ray = [&](RefractedRay& ray, GroundGeometry ground, double microphone) {
            if (layered_sound_speed)
                layered_finish_ground(*layered_sound_speed, ray, ground, microphone, layer_scratch);
            else
                finish_ground(ray, ground, microphone);
        };
        unsigned generation = 0;
        std::size_t culled_here = 0;
        std::size_t evaluated_here = 0;
        // Layered absorption: each state's mean coefficients, when the observers
        // share a height; otherwise the running integral at the observer's
        // height, once per observer, and the per-band mean over each path's
        // heights, once per evaluation.  Without a profile every path takes the
        // one run-wide array, as before.
        const std::size_t band_count = layered ? database.band_centers().size() : 0;
        std::vector<double> observer_absorption_integral(band_count, 0.0);
        std::vector<double> path_absorption(band_count, 0.0);
        auto absorption_for = [&](std::size_t it, const Vec3& observer) -> const double* {
            if (!layered)
                return absorption.data();
            if (!state_absorption_mean.empty())
                return state_absorption_mean.data() + it * band_count;
            layered->mean_coefficients(source_states[it].position.z,
                    state_absorption_integral.data() + it * band_count, observer.z,
                    observer_absorption_integral.data(), path_absorption.data());
            return path_absorption.data();
        };
        // Refraction: the moments at each observer's receiver, and the shadow
        // zone's per-band loss.
        std::vector<double> shadow(effective_sound_speed ? database.band_centers().size() : 0, 0.0);
        std::vector<double> table_loss(options.propagation_table ? database.band_centers().size() : 0, 0.0);
        std::size_t tabled_here = 0, fallback_here = 0;
        std::vector<double> shadow_weight(shadow.size(), 1.0);
        std::size_t shadowed_here = 0;
        // --refraction stratified's other eigenrays: each one's loss from its
        // earlier bounces, its band levels, and a ground weight of zero for one
        // that never reaches the ground.
        std::vector<double> arrival_loss(shadow.size(), 0.0);
        std::vector<double> arrival_bands(shadow.size(), 0.0);
        const std::vector<double> no_ground(shadow.size(), 0.0);
        for (std::size_t index = begin; index < end; index++)
        {
            // Evaluated at the microphone, with this observer's receiver.
            const auto& observer = microphones[index];
            if (layered)
                layered->integral(observer.z, observer_absorption_integral.data());
            const bool flush = on_plate[index] != 0;
            const PlateReceiver* plate_receiver = flush ? &*plate : nullptr;
            const GroundModel* ground_here = ground_model == nullptr ? nullptr
                    : flush                                  ? ground_model
                                                             : &point_models[point_of[index]];
            const GroundTable* table_here = ground_model != nullptr && !flush && !point_tables.empty()
                    ? &point_tables[point_of[index]]
                    : nullptr;
            const double plate_mic_azimuth = flush && options.plate_orientation.outboard
                    ? outboard_azimuth(observer, source_states)
                    : options.plate_orientation.azimuth;
            // The straight-ray spectral level: from --propagation_table where
            // it covers the pair, with the table's level in place of the ground
            // model's and spreading and absorption as ever; otherwise the run's
            // own model.
            // The table's level for a pair it covers, into `level`; false where
            // it does not.  A plate microphone reads it at the table's lowest
            // receiver height and adds K at the pair's azimuth in its frame.
            auto table_level = [&](std::size_t it, double* bands_out, double& level) {
                const Vec3& source = source_states[it].position;
                const Vec3 separation = observer - source;
                const auto east_north = compass.east_north({separation.x, separation.y});
                const double east = east_north[0];
                const double north = east_north[1];
                const double height = flush ? options.propagation_table->axes().receiver_heights.front()
                                            : microphone_height[index];
                if (!options.propagation_table->lookup(std::atan2(east, north) * 180.0 / pi, source.z, height,
                            separation.horizontal_norm(), table_bands, table_loss.data()))
                {
                    fallback_here++;
                    return false;
                }
                double k_weight = 0.0;
                std::size_t k_lo = 0, k_hi = 0;
                if (flush && !table_plate_offset.empty())
                {
                    const double azimuth = std::atan2(separation.y, separation.x) - plate_mic_azimuth + 0.5 * pi;
                    double x = std::fmod(azimuth / (2.0 * pi) * table_plate_azimuths, double(table_plate_azimuths));
                    if (x < 0.0)
                        x += table_plate_azimuths;
                    k_lo = std::min(static_cast<std::size_t>(x), std::size_t(table_plate_azimuths - 1));
                    k_hi = (k_lo + 1) % table_plate_azimuths;
                    k_weight = x - static_cast<double>(k_lo);
                }
                for (std::size_t band = first_band; band < last_band; band++)
                {
                    double gain = table_loss[band];
                    if (flush && !table_plate_offset.empty())
                        gain += (1.0 - k_weight) * table_plate_offset[k_lo][band]
                                + k_weight * table_plate_offset[k_hi][band];
                    table_loss[band] = -gain;
                }
                tabled_here++;
                level = database.spectral_level_along(straight_ray(source, observer), source_states[it],
                        table_absorption.empty() ? absorption_for(it, observer) : table_absorption.data(),
                        band_gains, nullptr, nullptr, sound_speed, first_band, last_band, packed_of_state[it],
                        nullptr, 0.0, nullptr, table_loss.data(), nullptr, bands_out);
                return true;
            };
            auto straight_level = [&](std::size_t it, double* bands_out) {
                double level = 0.0;
                if (options.propagation_table && table_level(it, bands_out, level))
                    return level;
                return database.spectral_level_at_observer(observer, source_states[it], absorption_for(it, observer),
                        band_gains, ground_here, &ground_impedance_bands, sound_speed, first_band, last_band,
                        packed_of_state[it], plate_receiver, plate_mic_azimuth, table_here, bands_out);
            };
            const GroundGeometry ray_ground = ground_model == nullptr ? GroundGeometry::none
                    : flush                                           ? GroundGeometry::plate
                                                                      : GroundGeometry::point;
            // Refraction: the table's moments at this microphone, and the level
            // of a state whose direct ray is formed -- its ground geometry
            // (the reflected ray, for a point receiver) is solved only now,
            // once the state has survived any bound.
            const auto at_receiver = effective_sound_speed
                    ? effective_sound_speed->moments(observer.z)
                    : EffectiveSoundSpeed::Moments {};
            const ShadowTransition* transition_here = flush
                    ? (plate_transition ? &*plate_transition : nullptr)
                    : (!shadow_transitions.empty() ? &shadow_transitions[point_of[index]] : nullptr);
            auto refracted_level = [&](std::size_t it, RefractedRay ray) {
                // With a table too, the pairs it covers take it; the rest the ray model.
                if (options.propagation_table)
                {
                    double level = 0.0;
                    double* table_bands_out = write_band_history
                            ? record_path(it, straight_ray(source_states[it].position, observer))
                            : nullptr;
                    if (table_level(it, table_bands_out, level))
                        return add_ambient ? level_with_ambient(level, ambient_energy) : level;
                }
                finish_ray(ray, ray_ground, microphone_height[index]);
                double* bands_out = record_path(it, ray.path);
                const double* extra_loss = nullptr;
                const double* ground_weight = nullptr;
                // The plate's frame, as Database::spectral_level_along takes it.
                const double plate_azimuth = flush ? plate_frame_azimuth(ray.path.separation, plate_mic_azimuth) : 0.0;
                if (transition_here != nullptr
                        && transition_here->apply(ray, shadow_attenuation, options.shadow_floor,
                                first_band, last_band, shadow_weight.data(), shadow.data(),
                                plate_azimuth))
                {
                    // Near the boundary the transition's field replaces the
                    // ground term and the creeping-wave loss.
                    extra_loss = shadow.data();
                    ground_weight = shadow_weight.data();
                    // Past the layered boundary counts, lit linear fit or not;
                    // the fit's ray from handover_fraction of it on does not.
                    if (ray.shadow_depth > 0.0
                            || (ray.from_linear_fit && ray.horizontal > ray.direct.boundary))
                        shadowed_here++;
                }
                else if (ray.shadow_depth > 0.0)
                {
                    shadow_loss(ray, shadow_attenuation, options.shadow_floor, first_band, last_band,
                            shadow.data());
                    extra_loss = shadow.data();
                    shadowed_here++;
                }
                // A primary that bounced before arriving (a duct, past every
                // direct ray) loses what those bounces took.
                if (ray.bounces > 0 && ground_here != nullptr)
                {
                    if (extra_loss != nullptr)
                        std::copy(extra_loss, extra_loss + shadow.size(), arrival_loss.begin());
                    else
                        std::fill(arrival_loss.begin(), arrival_loss.end(), 0.0);
                    add_bounce_loss(ray.ground_grazing, ray.bounces, ground_impedance_bands, first_band,
                            last_band, arrival_loss.data());
                    extra_loss = arrival_loss.data();
                }
                // All arrivals share these endpoints and their mean absorption.
                const double* pair_absorption = absorption_for(it, observer);
                double level = database.spectral_level_along(ray.path, source_states[it],
                        pair_absorption, band_gains, ground_here,
                        &ground_impedance_bands, sound_speed, first_band, last_band,
                        packed_of_state[it], plate_receiver, plate_mic_azimuth, table_here,
                        extra_loss, ground_weight, bands_out);
                if (!ray.extra.empty())
                {
                    // The other eigenrays, each its own level, summed in energy.
                    double energy = std::pow(10.0, 0.1 * level);
                    for (const Arrival& a : ray.extra)
                    {
                        std::fill(arrival_loss.begin(), arrival_loss.end(), 0.0);
                        if (ground_here != nullptr)
                            add_bounce_loss(a.grazing, a.bounces, ground_impedance_bands, first_band,
                                    last_band, arrival_loss.data());
                        if (a.wave.active())
                            add_wave_loss(a.wave, (observer - source_states[it].position).norm(),
                                    a.path.range, database.band_centers(), first_band, last_band,
                                    arrival_loss.data());
                        // One that tunnels reaches the microphone as the field's
                        // tail below its turning point, plate or not.
                        const bool grounded = (a.grazing > 0.0 || flush) && !(a.wave.evanescent > 0.0);
                        const double here = database.spectral_level_along(a.path, source_states[it],
                                pair_absorption, band_gains, ground_here,
                                &ground_impedance_bands, sound_speed, first_band, last_band,
                                packed_of_state[it], plate_receiver, plate_mic_azimuth, table_here,
                                arrival_loss.data(), grounded ? nullptr : no_ground.data(),
                                bands_out != nullptr ? arrival_bands.data() : nullptr);
                        energy += std::pow(10.0, 0.1 * here);
                        if (bands_out != nullptr)
                        {
                            for (std::size_t band = first_band; band < last_band; band++)
                                bands_out[band] = 10.0 * std::log10(std::pow(10.0, 0.1 * bands_out[band])
                                        + std::pow(10.0, 0.1 * arrival_bands[band]));
                        }
                    }
                    level = 10.0 * std::log10(energy);
                }
                if (add_ambient)
                    level = level_with_ambient(level, ambient_energy);
                return level;
            };
            if (expand_from_peak)
            {
                if (++generation == 0)
                {
                    // Wrapped, so old stamps could be mistaken for this pass.
                    std::fill(visited.begin(), visited.end(), 0u);
                    std::fill(travel_stamp.begin(), travel_stamp.end(), 0u);
                    generation = 1;
                }
                if (effective_sound_speed)
                {
                    integrated_levels[index] = windowed_exposure(nt, direction_bounds,
                            sphere_radius_squared, integration_step, sound_speed, levels, visited,
                            skipped_numerator, skipped_squared_range, generation, evaluated_here,
                            [&](std::size_t it) {
                                RefractedGeometry g;
                                g.ray = form_ray(it, observer, at_receiver, transition_here != nullptr);
                                const auto e = bound_geometry(it, g.ray);
                                g.squared_range = e.squared_range;
                                g.along = e.along;
                                g.up = e.up;
                                g.across = e.across;
                                travel[it] = sound_speed * g.ray.travel_time;
                                travel_stamp[it] = generation;
                                return g;
                            },
                            [&](std::size_t it, RefractedGeometry&& g) {
                                return refracted_level(it, std::move(g.ray));
                            },
                            [&](std::size_t it) {
                                if (travel_stamp[it] == generation)
                                    return travel[it];
                                return sound_speed * form_ray(it, observer, at_receiver, transition_here != nullptr).travel_time;
                            });
                }
                else
                {
                    // expand_from_peak holds only when the spectral bound was built.
                    integrated_levels[index] = windowed_sound_exposure_level(database,
                            source_states, direction_bounds, observer, integration_step,
                            sound_speed, levels, visited, skipped_numerator,
                            skipped_squared_range, generation, evaluated_here,
                            [&](std::size_t it, const StateGeometry&) {
                                return straight_level(it, nullptr);
                            });
                }
                continue;
            }
            if (effective_sound_speed)
            {
                for (std::size_t it = 0; it < nt; it++)
                {
                    auto ray = form_ray(it, observer, at_receiver, transition_here != nullptr);
                    // Reception time from the curved ray's travel time, as the
                    // range a straight ray at the surface sound speed would take.
                    ranges[it] = sound_speed * ray.travel_time;
                    if (culling)
                    {
                        const auto g = bound_geometry(it, ray);
                        if (g.squared_range > sphere_radius_squared
                                && direction_bounds.numerator(it,
                                           DirectionBounds::bin(g.along, g.up, g.across, g.squared_range))
                                        < cull_energy * g.squared_range)
                        {
                            levels[it] = culled_level;
                            culled_here++;
                            continue;
                        }
                    }
                    levels[it] = refracted_level(it, std::move(ray));
                }
            }
            else
            {
                for (std::size_t it = 0; it < nt; it++)
                {
                    // Range for the reception-time integral (see reception_factor).
                    ranges[it] = (observer - source_states[it].position).norm();
                    if (culling && is_culled(it, observer))
                    {
                        levels[it] = culled_level;
                        culled_here++;
                        continue;
                    }
                    double* bands_out = write_band_history
                            ? record_path(it, straight_ray(source_states[it].position, observer))
                            : nullptr;
                    double level = spectral ? straight_level(it, bands_out)
                                            : database.level_at_observer(observer, source_states[it]);
                    if (add_ambient)
                    {
                        // Add ambient noise
                        level = level_with_ambient(level, ambient_energy);
                    }
                    levels[it] = level;
                }
            }
            integrated_levels[index] = options.full_history
                    ? full_sound_exposure_level(levels, integration_step, ranges, sound_speed)
                    : sound_exposure_level(levels, integration_step, ranges, sound_speed);
            if (options.want_lamax || options.want_nax)
            {
                // The history is complete here (needs_full_history), so this is
                // the true peak over the whole operation, not a windowed one.
                const double peak = *std::max_element(levels.begin(), levels.end());
                if (options.want_lamax)
                    lamax[index] = std::max(lamax[index], peak);
                // NAx: operations whose own peak reaches the threshold.  Each
                // trajectory is one operation, so the count adds per trajectory.
                if (options.want_nax && peak >= options.nax_threshold)
                    nax[index] += 1.0;
            }
            if (options.want_coverage)
            {
                // The energy-weighted share of this observer's operation that
                // the database did not measure.  Every state's level is in
                // `levels` (the report turns bounding off), so the share is
                // summed over the same states the SEL integrated.
                double energy = 0.0;
                double unmeasured_energy = 0.0;
                for (std::size_t it = 0; it < nt; it++)
                {
                    if (!std::isfinite(levels[it]))
                        continue;
                    const double weight = std::pow(10.0, 0.1 * levels[it]);
                    const Vec3 direction = (observers[index] - source_states[it].position).normalized();
                    const auto angles = lookup_angles(source_states[it].frame, direction,
                            !database.upper_hemisphere_measured());
                    energy += weight;
                    unmeasured_energy += weight
                            * database.unmeasured_share(angles, source_states[it].blend);
                }
                coverage[index][0] += unmeasured_energy;
                coverage[index][1] += energy;
            }
            if (options.want_dnl)
            {
                // Energies add across operations; the logarithm comes at the
                // end, once every trajectory has contributed.  Each observer
                // index belongs to exactly one worker, and the trajectories run
                // one after another, so this needs no synchronization.
                dnl_energy[index] += std::pow(
                        10.0, 0.1 * day_night_level(levels, integration_step, night, penalized, ranges,
                                              sound_speed));
            }
            if (options.write_time_history)
            {
                // Write out time history files
                auto time_history_filename = options.time_history_prefix + history_infix + "_i"
                        + std::to_string(index) + ".csv";
                write_time_history_file(levels, integration_step, observer(0), observer(1),
                        time_history_filename, options.unit);
            }
            if (write_band_history)
            {
                for (std::size_t it = 0; it < nt; it++)
                {
                    band_rows[it].emission_time = trajectory.states[it].t;
                    band_rows[it].reception_time = trajectory.states[it].t + ranges[it] / sound_speed;
                }
                write_band_history_file(band_rows, levels, band_levels, database.band_centers(),
                        first_band, last_band, observer, options.unit,
                        options.band_history_prefix + history_infix + "_i" + std::to_string(index)
                                + ".csv");
            }
        }
        culled_samples += culled_here;
        evaluated_samples += evaluated_here;
        shadowed_samples += shadowed_here;
        table_samples += tabled_here;
        fallback_samples += fallback_here;
    });

    if (options.propagation_table)
    {
        const auto looked_up = table_samples.load() + fallback_samples.load();
        log("The propagation table covered ", table_samples.load(), " of the ", looked_up, " evaluations made (",
                100.0 * static_cast<double>(table_samples.load()) / static_cast<double>(std::max<std::size_t>(looked_up, 1)),
                "%); the rest -- outside its ranges or source heights, or where it has no value, as "
                "steeper than a PE's angle limit -- took ",
                ground_model ? "the ground model." : "free-field propagation.");
    }

    if (effective_sound_speed)
    {
        const auto total = observers.size() * nt;
        const auto evaluations = expand_from_peak ? evaluated_samples.load()
                                                  : total - culled_samples.load();
        log("Refraction put ", shadowed_samples.load(), " of the ", evaluations,
                " evaluations made (",
                100.0 * static_cast<double>(shadowed_samples.load())
                        / static_cast<double>(std::max<std::size_t>(evaluations, 1)),
                "%) in a shadow zone.");
    }

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
    write_levels(x_obs, y_obs, integrated_levels, job.output, options.scattered,
            options.unit);
    if (options.latlon)
    {
        write_geodetic_levels(x_obs, y_obs, integrated_levels, *options.geodetic.frame,
                latlon_companion(job.output), options.scattered);
    }

    if (!options.metric_file.empty())
    {
        auto casename = std::filesystem::path(job.trajectory).stem().string();
        append_metrics_file(options.metric_file, casename, integrated_levels);
    }

    if (options.write_state)
    {
        trajectory.write_state(job.state);
    }

    return trajectory;
}

// Draw the footprint and the state figures for a single-trajectory run.
//
// Lifted out of run() unchanged.  It takes x_obs/y_obs by value because the
// first thing it does is convert them out of meters, and run() still holds the
// metric copies -- in place, that conversion used to run over the caller's
// arrays.
void draw_figures(const Trajectory& trajectory, const std::vector<double>& integrated_levels,
        std::vector<double> x_obs, std::vector<double> y_obs, bool scattered, bool plotting,
        const std::string& figure_file, const std::string& state_figure_file,
        double contour_min, double contour_max, double contour_step, LengthUnit unit)
{
    const std::string x_label = std::string("X, ") + unit.name;
    auto x = trajectory.column("x");
    auto y = trajectory.column("y");
    auto z = trajectory.column("z");

    // Out of meters and into whatever --units selected, which is also what the
    // axes are labeled with.
    for (auto& v : x_obs)
        v *= 1.0 / unit.to_meters;
    for (auto& v : y_obs)
        v *= 1.0 / unit.to_meters;
    for (auto& v : x)
        v *= 1.0 / unit.to_meters;
    for (auto& v : y)
        v *= 1.0 / unit.to_meters;
    for (auto& v : z)
        v *= 1.0 / unit.to_meters;

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
        plt::xlabel(x_label, 14);
        plt::ylabel(std::string("Y, ") + unit.name, 14);
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

        // One figure per state variable against x, all four identical but
        // for the series and the y label.  NOTE: the separator lands before
        // the extension, so these come out as "state_CT_.png" -- the D
        // original's naming, which the reference figures carry.
        auto state_figure = [&](const char* name, const char* label,
                                    const std::vector<double>& series) {
            plt::figure();
            plt::plot(x, series);
            plt::xlabel(x_label, 14);
            plt::ylabel(label, 14);
            if (!state_figure_file.empty())
                plt::savefig(base + "_" + name + "_" + ext);
        };

        state_figure("CT", "$C_T$", CT);
        state_figure("mu", "$\\mu$", mu);
        // The one exception: both flight path angles on the same axes.
        plt::figure();
        plt::plot(x, gamma, x, gamma_e);
        plt::xlabel(x_label, 14);
        plt::ylabel("$\\gamma$", 14);
        if (!state_figure_file.empty())
        {
            plt::savefig(base + "_gamma_" + ext);
        }
        state_figure("z", (std::string("Z, ") + unit.name).c_str(), z);
    }
    if (plotting)
    {
        plt::show();
    }
    plt::clear();
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
    std::string lamax_file = "";
    std::string nax_file = "";
    double nax_threshold = std::numeric_limits<double>::quiet_NaN();
    std::string coverage_file = "";
    bool full_history = false;
    double utc_offset = 0.0;
    double trim_start = std::numeric_limits<double>::quiet_NaN();
    // Lengths on this command line are in whatever --units selects, feet by
    // default, and are converted to meters at one seam below.  A *default* is a
    // fixed physical quantity, not a number in the selected unit, so it is
    // written here in meters and skips that conversion -- otherwise --units
    // would silently change the geometry of a run that named no grid.
    std::string unit_name = "ft";
    std::vector<double> grid_bounds;
    double grid_step = 25.0; // meters
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
    // -a is given.  Initializing this to 0.0 added a spurious 0 dB ambient
    // floor to every run, which pinned quiet observers near 0 dB.
    double ambient_level = std::numeric_limits<double>::quiet_NaN();
    bool write_time_history = false;
    std::string time_history_prefix = "";
    std::string band_history_prefix;
    std::string metric_file = "";
    bool qsam = true;
    bool self_test = false;
    unsigned thread_count = 0;
    double db_cutoff = std::numeric_limits<double>::quiet_NaN();
    bool spectral = true;
    std::vector<double> band_limits;
    bool ground_effect = true;
    double max_load_factor = 2.0;
    bool accept_vehicle_mismatch = false;
    double flow_resistance = flow_resistance::grassy_field;
    double receiver_height = 1.2; // meters
    double turbulence = 0.0;
    std::string impedance_model = "delany_bazley";
    double porosity_rate = 0.0;      // 1/m
    double layer_depth = 0.05;       // m
    std::string turbulence_model = "chessell";
    double gamma_t = 0.0;
    double roughness = 0.0;          // m
    std::string plate_table;
    std::string propagation_table;
    std::string plate_mic_azimuth = "outboard";
    bool exact_ground = false;
    double humidity_override = std::numeric_limits<double>::quiet_NaN();
    double temperature_override = std::numeric_limits<double>::quiet_NaN();
    double pressure_override = std::numeric_limits<double>::quiet_NaN();
    std::string atmosphere_file;
    // --wind and the options that shape the profile built from it.  The
    // lengths are meters until the --units seam, like every default.
    std::string wind_text;
    double roughness_length = SurfaceWind().roughness_length;
    double surface_layer_depth = SurfaceWind().surface_layer_depth;
    std::string stability = "neutral";
    int profile_levels = SurfaceWind().levels_per_decade;
    std::string write_atmosphere;
    std::string refraction = "none";
    double shadow_floor = 20.0;
    bool tunneling = false;
    bool skip_fill = false;
    double scattering = 0.0;
    // Empty follows the database's speed_reference attribute.
    std::string speed_reference;
    // WGS84 input and output.  An unset origin is taken from the first geodetic
    // trajectory; an unset bearing leaves the axes east/north, which is what a
    // map wants -- the shipped AC tracks are the other case, their +x running
    // along the flight track.
    std::vector<double> origin;
    double frame_bearing = std::numeric_limits<double>::quiet_NaN();
    double ground_elevation = 0.0;
    double start_utc = std::numeric_limits<double>::quiet_NaN();
    // Waypoints state intent, not a constraint to be met to the inch: forcing
    // a path exactly through them manufactures curvature, and the load factor
    // is that curvature, so it lands in CT and then in the levels.  Ten feet is
    // loose enough to absorb both coordinate rounding and the corner a planner
    // did not mean literally.
    double path_smoothing = 10.0 * ft2m; // 10 ft, in meters
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
            "Output file name, once per trajectory. Default integrated/197421_LAE.csv; with several "
            "trajectories and one name (or the default), each is named after its track.");
    app.add_flag("-p,--plotting", plotting, "Enable plotting.");
    app.add_option("-f,--figfile", figure_file, "Filename to save figure.");
    app.add_option("-F,--statefigfile", state_figure_file, "Filename to save state figures.");
    app.add_option("--cmin", contour_min, "Minimum contour level, dB.");
    app.add_option("--cmax", contour_max, "Maximum contour level, dB.");
    app.add_option("--cstep", contour_step, "Contour level step, dB.");
    app.add_option("-s,--step", dt, "Resampling time step for --adapt, in seconds. "
            "The sound exposure level integral uses the trajectory's own sample interval.");
    app.add_option("-a,--ambient", ambient_level,
            "Ambient noise level, dB, added to every sample's level in energy.");
    app.add_option("-i,--history", time_history_prefix, "Time history file prefix.");
    app.add_option("--band_history", band_history_prefix,
            "Per-band time history file prefix: one file per observer, a row per trajectory "
            "state with its emission and reception times, geometry, A-weighted level and each "
            "band's unweighted level at the microphone. Needs the spectral path.");
    app.add_option("-m,--metricfile", metric_file, "File to append overall noise metrics.");
    app.add_option("-b,--observer", observer_file, "Observer file name.");
    app.add_flag("-q,--qsam,!--no-qsam", qsam, "Extended QSAM (default on; --no-qsam disables)");
    app.add_option("-S,--state", state_files,
            "Filename to save state data, once per trajectory. With several trajectories and one "
            "name, each is named after its track (<stem>_state.csv).");
    app.add_flag("--self_test", self_test, "Check every channel of the database against what eval() returns for it, then exit.");
    app.add_option("-j,--threads", thread_count,
            "Worker threads for the observer loop (0 = one per hardware thread).");
    auto* trim_option = app.add_option("--trim", trim_start,
            "Crop each trajectory from this x to its lowest/slowest point, in --units. "
            "Off by default.");
    app.add_option("--units", unit_name,
            "Length unit for the command line and for every file carrying a length: "
            "ft (default) or m. Angles stay in degrees, time in seconds, levels in dB, "
            "and the -S state file stays metric.")
            ->check(CLI::IsMember({"ft", "m"}));
    auto* grid_option = app.add_option("--grid", grid_bounds,
            "Observer grid bounds, in --units: xmin xmax ymin ymax. Default +/-1500 m "
            "(+/-4921.26 ft).")
            ->expected(4);
    auto* gridstep_option = app.add_option("--gridstep", grid_step,
            "Observer grid spacing, in --units (default 25 m, = 82.021 ft).");
    app.add_option("-D,--dnlfile", dnl_file,
            "Write the day-night level (DNL: the energy average over a 24-hour day, night "
            "penalized), summed over every trajectory, to this file.");
    app.add_option("--lamax_file", lamax_file,
            "Write each observer's LAmax, the highest A-weighted level over the operation, "
            "dB, to this file.");
    app.add_option("--nax_file", nax_file,
            "Write each observer's NAx, the number of operations whose LAmax reaches "
            "--nax_threshold, to this file.");
    app.add_option("--nax_threshold", nax_threshold,
            "The LAmax, dBA, that an operation must reach to count toward NAx. Required with "
            "--nax_file: there is no default, so the study's threshold is always stated.");
    app.add_option("--coverage_file", coverage_file,
            "Write, per observer, the energy-weighted share of its levels that came from cells "
            "the database did not measure (straight rays only; needs a database with coverage).");
    app.add_flag("--full_history", full_history,
            "Integrate each SEL over the whole history rather than the event's 10 dB-down "
            "window, as full-history tools do. Slower: no window is found.");
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
    auto* spectral_option = app.add_flag("--spectral,!--no-spectral", spectral,
            "Propagate the database's per-band spectra to the ground with ISO 9613-1 "
            "absorption (default on; falls back to broadband when the database carries no "
            "spectra). --no-spectral propagates one broadband level with the database's EAA "
            "factor instead.");
    app.add_option("--band_limits", band_limits,
            "Frequency range integrated into the metric, Hz: low high. Default 20 20000.")
            ->expected(2);
    app.add_option("--max_load_factor", max_load_factor,
            "Largest load factor, in g, a fixed_load_factor database's load-factor correction "
            "is taken to. Above it the correction is held, and the run says so (default 2).")
            ->check(CLI::PositiveNumber);
    app.add_flag("--accept_vehicle_mismatch", accept_vehicle_mismatch,
            "Warn, rather than refuse, a vehicle file that disagrees with the database by more "
            "than 2 percent.");
    app.add_option("--humidity", humidity_override,
            "Relative humidity in percent, overriding the vehicle file.");
    app.add_option("--temperature", temperature_override,
            "Temperature in K, overriding the vehicle file.");
    app.add_option("--pressure", pressure_override,
            "Ambient pressure in kPa, overriding the vehicle file.");
    app.add_option("--atmosphere", atmosphere_file,
            "A stratified atmosphere: CSV of temperature, relative humidity and optionally "
            "pressure against height (see README). Per-band absorption is integrated through "
            "it, and its ground level replaces the vehicle file's atmosphere. Spectral path "
            "only.");
    app.add_option("--wind", wind_text,
            "Build the stratified atmosphere from one surface wind instead of an --atmosphere file: "
            "SPEED,DIRECTION[,HEIGHT]. SPEED names its unit (12kt, 6.2mps or 14mph); DIRECTION is "
            "the compass bearing, degrees, the wind blows FROM (as an --atmosphere file's "
            "wind_direction), turned into the local frame by --frame_bearing; HEIGHT is the "
            "anemometer's, in --units (default 10 m). The wind follows the neutral log law over "
            "--roughness_length up to --surface_layer_depth and is held above; the temperature "
            "lapses dry-adiabatically from the surface air (vehicle file, --temperature, "
            "--humidity, --pressure), pressure is hydrostatic and the water vapor is held. "
            "Not with --atmosphere. Spectral path only.");
    auto* roughness_length_option = app.add_option("--roughness_length", roughness_length,
            "--wind: aerodynamic roughness length z0 of the ground, in --units (default 0.03 m, "
            "= 0.0984 ft: open flat terrain, short grass, airfields; Davenport-Wieringa). "
            "Not --roughness, which is the acoustic ground's RMS height.");
    auto* surface_layer_option = app.add_option("--surface_layer_depth", surface_layer_depth,
            "--wind: height, in --units, up to which the wind follows its law; above it the wind "
            "is held (default 200 m, = 656.2 ft).");
    auto* stability_option = app.add_option("--stability", stability,
            "--wind: neutral (default; the log law, the one checked against measurements), or a "
            "TA Luft class I, II, III/1, III/2, IV, V (Pasquill F to A): a power-law wind and a "
            "temperature gradient in the surface layer. The classes are unvalidated.")
            ->check(CLI::IsMember(stability_names(), CLI::ignore_case));
    auto* profile_levels_option = app.add_option("--profile_levels", profile_levels,
            "--wind: levels per factor of ten in height in the built profile (default 16). Results "
            "do not depend on it; raise it to check that they do not.");
    app.add_option("--write_atmosphere", write_atmosphere,
            "--wind: write the profile built for each trajectory as an --atmosphere CSV (z, T_K, "
            "RH, p_kPa, u_mps, v_mps), so the run can be repeated and inspected. With several "
            "trajectories and one name, each is named after its track "
            "(<track>_atmosphere.csv).");
    app.add_option("--refraction", refraction,
            "Refraction through the --atmosphere (or --wind) profile: none (default); linear -- curved rays "
            "through its effective sound speed (temperature and wind) fitted to a straight line "
            "per path; or stratified -- through the profile's own layers, slower but without "
            "the fit. Both give shadow zones under upward refraction.")
            ->check(CLI::IsMember({"none", "linear", "stratified"}));
    app.add_flag("--tunneling", tunneling,
            "--refraction stratified: also take the rays turned back from aloft that turn again just "
            "above the microphone, which the wave reaches through the layer below them (an Airy "
            "tail per band).  Off by default.");
    app.add_flag("--skip_fill", skip_fill,
            "--refraction stratified: fill the skip zone short of where rays from aloft first land "
            "with the caustic's Airy tail per band.  Off by default.");
    app.add_option("--shadow_floor", shadow_floor,
            "Most a refractive shadow zone attenuates, dB (default 20). Deeper in, energy "
            "scattered by turbulence dominates.");
    app.add_option("--scattering", scattering,
            "Fill refraction shadows with the energy turbulence scatters into them instead of "
            "--shadow_floor: the HARMONOISE strength gamma_T = C_T^2/T0^2 + 22 C_v^2/(3 c0^2) "
            "(0, the default, is off). Needs --refraction and --ground_effect.");
    auto* ground_effect_option = app.add_flag("--ground_effect,!--no-ground_effect", ground_effect,
            "Add excess ground attenuation (Chien-Soroka over a Delany-Bazley impedance) to "
            "the spectral path. On by default wherever the spectral path runs; "
            "--no-ground_effect turns it off. Asking for it with --no-spectral is an error.");
    app.add_option("--flow_resistance", flow_resistance,
            "Specific flow resistance of the ground, kPa s/m^2. 30 snow, 225 grass (default), "
            "650 soil, 3000 dirt, 10000 concrete, 50000 asphalt.");
    auto* receiver_height_option = app.add_option("--receiver_height", receiver_height,
            "Receiver height above ground for --ground_effect, in --units (default 1.2 m, "
            "= 3.937 ft).")
            ->check(CLI::NonNegativeNumber);
    app.add_option("--origin", origin,
            "WGS84 reference for the local frame: latitude longitude, in degrees. "
            "Defaults to the first position of the first geodetic trajectory.")
            ->expected(2);
    app.add_option("--frame_bearing", frame_bearing,
            "Compass bearing, degrees, of the local +x axis. Default: x east, y north. "
            "The shipped AC tracks use 210.");
    app.add_option("--speed_reference", speed_reference,
            "Which velocity sets the flight condition -- mu, the path angle and the sphere's "
            "nose: ground, or air (the ground velocity less the --atmosphere profile's wind at "
            "the aircraft). Default: what the database records it was built with; ground if it "
            "records nothing.")
            ->check(CLI::IsMember({"ground", "air"}));
    auto* ground_elevation_option = app.add_option("--ground_elevation", ground_elevation,
            "Elevation of the observer plane above mean sea level, in --units. Subtracted "
            "from a geodetic track's altitude column to give height above the observers "
            "(default 0).");
    app.add_option("--start_utc", start_utc,
            "Seconds since UTC midnight at t = 0 of a geodetic track, for the day-night "
            "window (-D). Required with -D when the track has no utcsec column; ignored "
            "when it has one.");
    auto* path_smoothing_option = app.add_option("--path_smoothing", path_smoothing,
            "How far, RMS in --units, the path fitted through a geodetic track may sit from "
            "the given positions, per axis. Default 10 ft. Zero flies through them exactly, "
            "which on measured positions makes the load factor track coordinate rounding.");
    app.add_flag("--latlon", latlon,
            "Also write each level file as latitude, longitude, level. Needs --origin unless "
            "the trajectory is itself geodetic.");
    app.add_option("--turbulence", turbulence,
            "Turbulence parameter for --ground_effect, rad s/sqrt(m); 0 to 16e-4. "
            "Zero gives unrealistically deep interference nulls. (--turbulence_model chessell)");
    app.add_option("--impedance_model", impedance_model,
            "Ground impedance model for --ground_effect: delany_bazley (default), miki, "
            "variable_porosity (with --porosity_rate) or delany_bazley_layer (with "
            "--layer_depth). --flow_resistance is the (effective) flow resistivity.")
            ->check(CLI::IsMember({"delany_bazley", "miki", "variable_porosity", "delany_bazley_layer"}));
    app.add_option("--porosity_rate", porosity_rate,
            "variable_porosity: rate of change of porosity with depth, alpha_e, 1/m (a material "
            "parameter, in meters whatever --units says).");
    auto* layer_depth_option = app.add_option("--layer_depth", layer_depth,
            "delany_bazley_layer: depth of the hard-backed layer, in --units (default 0.05 m).");
    app.add_option("--turbulence_model", turbulence_model,
            "chessell (default; --turbulence) or harmonoise (--gamma_t): the HARMONOISE "
            "coherence of the direct and reflected paths.")
            ->check(CLI::IsMember({"chessell", "harmonoise"}));
    app.add_option("--gamma_t", gamma_t,
            "HARMONOISE turbulence strength C_T^2/T0^2 + 22 C_v^2/(3 c0^2); ~1e-6 moderate, "
            "~1e-5 strong. The 2017 Amedee pole data fitted ~3e-5.");
    auto* roughness_option = app.add_option("--roughness", roughness,
            "RMS roughness height of the ground, in --units: the coherent reflection falls by "
            "exp(-2 (k sigma cos theta)^2). Default 0.");
    app.add_flag("--exact_ground", exact_ground,
            "Evaluate the --ground_effect model (or --plate_table's reflection coefficients) at "
            "every observer and state instead of interpolating tables of it built per "
            "trajectory (level histories agree to ~0.001 dB; the tables are 3-9x faster).");
    app.add_option("--plate_table", plate_table,
            "Observers are flush microphones on a rigid ground plate (e.g. GRAS 67AX in a "
            "400 mm GR1425 plate, per SAE ARP 4055): a netCDF table from panam's "
            "axisymmetric_bem.write_netcdf, computed over the same ground as --impedance_model "
            "and --flow_resistance. Replaces the point receiver of --ground_effect, which it needs.");
    app.add_option("--plate_mic_azimuth", plate_mic_azimuth,
            "Direction the microphone is offset on its plate: 'outboard' (default; away from "
            "the track) or degrees in the local frame, from +x toward +y.");
    app.add_option("--propagation_table", propagation_table,
            "A propagation loss table (netCDF; tools/propagation_table.cpp builds one with the "
            "parabolic equation): each band's level re free field by bearing, source and receiver "
            "height and range, for one stratified atmosphere over uniform ground. It replaces the "
            "ground and ray models wherever it covers a pair; elsewhere --ground_effect (or free "
            "field) applies. Spectral path only; not with --refraction or --plate_table.");

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
    std::vector<std::string> atmosphere_out_files;
    if (!write_atmosphere.empty())
    {
        if (path_count == 1)
            atmosphere_out_files.push_back(write_atmosphere);
        else
        {
            // As the state file: one name for several tracks would hold only the last profile.
            const auto directory = std::filesystem::path(write_atmosphere).parent_path();
            for (const auto& path : path_files)
            {
                const auto stem = std::filesystem::path(path).stem().string();
                atmosphere_out_files.push_back((directory / (stem + "_atmosphere.csv")).string());
            }
        }
    }
    if (state_files.size() == 1 && path_count > 1)
    {
        // The same goes for the state file: one name for several tracks would
        // be rewritten by each in turn, leaving only the last one's.
        const auto directory = std::filesystem::path(state_files.front()).parent_path();
        state_files.clear();
        for (const auto& path : path_files)
        {
            const auto stem = std::filesystem::path(path).stem().string();
            state_files.push_back((directory / (stem + "_state.csv")).string());
        }
    }
    expand_per_path(database_files, "database files");
    expand_per_path(vehicle_files, "vehicle files");
    expand_per_path(out_files, "output files");
    if (!state_files.empty())
        expand_per_path(state_files, "state files");
    // Two trajectories writing one file is never what was meant: the later
    // overwrites the earlier and the run still exits 0.  This catches equal
    // track stems (tracks from two directories) and an explicit -o or -S given
    // the same name twice.
    auto require_distinct = [](const std::vector<std::string>& names, const char* what) {
        std::set<std::string> seen;
        for (const auto& name : names)
            if (!seen.insert(std::filesystem::path(name).lexically_normal().string()).second)
                input_error("Two trajectories would write the same ", what, " (", name,
                        "). Name them apart.");
    };
    if (path_count > 1)
    {
        require_distinct(out_files, "output file");
        if (!state_files.empty())
            require_distinct(state_files, "state file");
        if (!atmosphere_out_files.empty())
            require_distinct(atmosphere_out_files, "--write_atmosphere file (their tracks share a file name)");
        if (!time_history_prefix.empty() || !band_history_prefix.empty())
        {
            std::vector<std::string> stems;
            for (const auto& path : path_files)
                stems.push_back(std::filesystem::path(path).stem().string());
            require_distinct(stems, "history file name (their tracks share a file name)");
        }
    }

    if (grid_bounds.empty())
    {
        grid_bounds = {-1500.0, 1500.0, -1500.0, 1500.0}; // meters
    }
    if (!(grid_step > 0.0))
    {
        input_error("--gridstep must be positive, got ", grid_step, ".");
    }
    // Say so up front rather than running the whole case and drawing nothing at
    // the end: this check was written to do that, and sat after the run.
    const bool want_plots = plotting || !figure_file.empty() || !state_figure_file.empty();
    if (want_plots && !plt::available())
    {
        input_error("This build has no plotting support. Reconfigure with "
              "-DNICEOPS_PLOTTING=ON, or drop -p/-f/-F.");
    }
    // The atmosphere overrides are checked where they are read, before anything
    // loads: a bad one used to run the whole case and write a footprint.
    if (!std::isnan(temperature_override))
        if (const auto problem = temperature_problem(temperature_override); !problem.empty())
            input_error("--temperature: ", problem);
    if (!std::isnan(pressure_override))
        if (const auto problem = pressure_problem(pressure_override); !problem.empty())
            input_error("--pressure: ", problem);
    if (!std::isnan(humidity_override))
        if (const auto problem = humidity_problem(humidity_override); !problem.empty())
            input_error("--humidity: ", problem);
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
    // Settled before anything that reads a length: the run's frame is found by
    // reading the first row of a geodetic track, and which altitude spellings
    // that read accepts depends on the unit (alt_m exists only under --units m).
    // This was declared further down and the origin lookup took the default, so
    // a metric track headed alt_m looked unpositioned, the frame was never
    // settled, and --frame_bearing was silently dropped.
    const LengthUnit unit = unit_name == "m" ? meters : feet;
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
            if (auto first = first_geodetic_position(path, unit))
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
        // The writers create with NC_CLOBBER.  `-d X --export_sphere X` replaced a
        // database built over hours in panam with a one-sphere file, and exited 0.
        for (const auto* output : {&sphere_file, &aam_file})
        {
            std::error_code ignored;
            if (output->empty() || !std::filesystem::exists(*output, ignored))
                continue;
            for (const auto& source : database_files)
            {
                if (std::filesystem::equivalent(*output, source, ignored))
                    input_error("--export_sphere/--export_aam would overwrite the database it reads, ", source,
                                ". Give another file name.");
            }
        }
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
        // And stop: the rest of main is the demo case, whose output file the
        // golden tests read, so a self test that ran it too rewrote that file
        // under any test running beside it.
        return 0;
    }

    // The one place the user's units become meters.  Everything above is in
    // what they typed; everything below is metric, as the computation is.
    //
    // Only options the user actually supplied are converted: a default is a
    // physical quantity already written in meters above, and scaling it here
    // would make --units m move the observer grid by 3.28x rather than just
    // renaming its axis.  CLI11's count() is what tells the two apart.
    log("Lengths are in ", unit.name, ".");
    auto to_meters = [&unit](const CLI::Option* option, double& value) {
        if (option->count() > 0)
            value *= unit.to_meters;
    };
    if (grid_option->count() > 0)
    {
        for (double& bound : grid_bounds)
            bound *= unit.to_meters;
    }
    to_meters(gridstep_option, grid_step);
    to_meters(receiver_height_option, receiver_height);
    to_meters(layer_depth_option, layer_depth);
    to_meters(roughness_option, roughness);
    to_meters(path_smoothing_option, path_smoothing);
    to_meters(trim_option, trim_start);
    to_meters(ground_elevation_option, ground_elevation);
    to_meters(roughness_length_option, roughness_length);
    to_meters(surface_layer_option, surface_layer_depth);

    ObserverSet observer_set = scattered
            ? scattered_observers(observer_file, unit)
            : gridded_observers(grid_bounds, grid_step);
    auto& observers = observer_set.observers;
    auto& x_obs = observer_set.x;
    auto& y_obs = observer_set.y;
    std::vector<double> integrated_levels(observers.size());
    // DNL is a property of the whole set of operations, not of one flight: each
    // trajectory's day-night level is an energy, and they sum.
    const bool want_dnl = !dnl_file.empty();
    std::vector<double> dnl_energy(want_dnl ? observers.size() : 0, 0.0);
    const bool want_lamax = !lamax_file.empty();
    std::vector<double> lamax(want_lamax ? observers.size() : 0,
            -std::numeric_limits<double>::infinity());
    const bool want_nax = !nax_file.empty();
    std::vector<double> nax(want_nax ? observers.size() : 0, 0.0);
    const bool want_coverage = !coverage_file.empty();
    std::vector<std::array<double, 2>> coverage(want_coverage ? observers.size() : 0,
            std::array<double, 2> {0.0, 0.0});
    if (want_nax && !std::isfinite(nax_threshold))
        input_error("--nax_file needs --nax_threshold, the LAmax in dBA that counts an operation. "
                "There is no default: give the threshold your study uses.");

    // Outlives the loop so the plotting block below, which only runs for a
    // single trajectory, can still reach it.
    std::optional<Trajectory> trajectory_holder;

    // The flag combinations that cannot be served, checked before anything is
    // built from them.
    const bool spectral_explicit = spectral_option->count() > 0 && spectral;
    const bool ground_explicit = ground_effect_option->count() > 0 && ground_effect;
    if (ground_explicit && !spectral)
    {
        input_error("--ground_effect needs the spectral path: the ground effect is a function "
                    "of frequency, and the broadband path (--no-spectral or --db_cutoff) has no "
                    "spectrum to apply it to.");
    }
    if (!band_history_prefix.empty() && !spectral)
    {
        input_error("--band_history needs the spectral path: the broadband path (--no-spectral "
                    "or --db_cutoff) carries no bands.");
    }
    if (!propagation_table.empty() && !spectral)
        input_error("--propagation_table needs the spectral path.");
    if (!plate_table.empty() && !ground_effect)
    {
        input_error("--plate_table needs --ground_effect: the plate's level depends on the "
                    "ground's reflection coefficient.");
    }
    if (band_limits.size() == 2 && !(band_limits[0] < band_limits[1]))
    {
        input_error("--band_limits needs low < high, got ", band_limits[0], " and ",
                band_limits[1], ".");
    }
    if (!atmosphere_file.empty() && !spectral)
    {
        input_error("--atmosphere needs the spectral path: layered absorption is per band, and "
                    "the broadband path (--no-spectral) carries one level attenuated by the "
                    "database's own EAA factor.");
    }
    if (!wind_text.empty() && !atmosphere_file.empty())
    {
        input_error("--wind builds the atmosphere that --atmosphere would read; give one or the "
                    "other. (--write_atmosphere saves the built one in the file format.)");
    }
    if (!wind_text.empty() && !spectral)
    {
        input_error("--wind needs the spectral path: the profile it builds is integrated per band, "
                    "and the broadband path (--no-spectral) carries one level attenuated by the "
                    "database's own EAA factor.");
    }
    if (wind_text.empty())
    {
        for (const auto& [option, name] : {std::pair {roughness_length_option, "--roughness_length"},
                     std::pair {surface_layer_option, "--surface_layer_depth"},
                     std::pair {stability_option, "--stability"},
                     std::pair {profile_levels_option, "--profile_levels"}})
        {
            if (option->count() > 0)
                input_error(name, " shapes the profile --wind builds, and there is no --wind.");
        }
        if (!write_atmosphere.empty())
            input_error("--write_atmosphere writes the profile --wind builds, and there is no --wind.");
    }
    if (refraction != "none" && atmosphere_file.empty() && wind_text.empty())
    {
        input_error("--refraction ", refraction, " needs --atmosphere or --wind: it refracts "
                    "through the profile's temperature and wind.");
    }
    if (!(shadow_floor >= 0.0))
        input_error("--shadow_floor must not be negative, got ", shadow_floor, ".");
    if (!(scattering >= 0.0))
        input_error("--scattering must not be negative, got ", scattering, ".");
    if (scattering > 0.0 && (refraction == "none" || !ground_effect))
    {
        input_error("--scattering fills the shadows --refraction casts, through the shadow-boundary "
                    "field that --ground_effect's ground model sets up; give both.");
    }
    if (!atmosphere_file.empty()
            && (!std::isnan(temperature_override) || !std::isnan(humidity_override)))
    {
        input_error("--atmosphere gives temperature and humidity at every height; --temperature "
                    "and --humidity would contradict it.");
    }

    // Everything the trajectory loop needs that does not vary between
    // trajectories, named once.  Written as an initializer rather than thirty
    // assignments so a field added to RunOptions and forgotten here reads as
    // its declared default rather than as whatever the struct was left holding.
    RunOptions options {};
    options.qsam = qsam;
    options.dt = dt;
    options.ambient_level = ambient_level;
    options.thread_count = thread_count;
    options.db_cutoff = db_cutoff;
    options.spectral = spectral;
    options.spectral_explicit = spectral_explicit;
    if (band_limits.size() == 2)
        options.band_limits = BandLimits {band_limits[0], band_limits[1]};
    options.ground_effect = ground_effect;
    options.ground_effect_given = ground_explicit;
    options.max_load_factor = max_load_factor;
    options.accept_vehicle_mismatch = accept_vehicle_mismatch;
    options.ground.flow_resistance = flow_resistance;
    options.ground.receiver_height = receiver_height;
    options.ground.turbulence = turbulence;
    options.ground.impedance_model = parse_impedance_model(impedance_model);
    options.ground.porosity_rate = porosity_rate;
    options.ground.layer_depth = layer_depth;
    options.ground.turbulence_model = parse_turbulence_model(turbulence_model);
    options.ground.gamma_t = gamma_t;
    options.ground.roughness = roughness;
    options.exact_ground = exact_ground;
    if (!propagation_table.empty())
    {
        options.propagation_table = std::make_shared<const PropagationTable>(PropagationTable::load(propagation_table));
        options.propagation_table_path = propagation_table;
    }
    if (!plate_table.empty())
    {
        options.plate = PlateReceiver::load(plate_table);
        if (plate_mic_azimuth != "outboard")
        {
            options.plate_orientation.outboard = false;
            try
            {
                options.plate_orientation.azimuth = std::stod(plate_mic_azimuth) * pi / 180.0;
            }
            catch (const std::exception&)
            {
                input_error("--plate_mic_azimuth takes 'outboard' or degrees, not '",
                        plate_mic_azimuth, "'.");
            }
        }
    }
    options.humidity_override = humidity_override;
    options.temperature_override = temperature_override;
    options.pressure_override = pressure_override;
    options.refraction = refraction != "none";
    options.stratified = refraction == "stratified";
    if ((tunneling || skip_fill) && !options.stratified)
        input_error("--tunneling and --skip_fill need --refraction stratified.");
    options.waves.tunneling = tunneling;
    options.waves.skip_fill = skip_fill;
    options.shadow_floor = shadow_floor;
    options.scattering = scattering;
    options.frame_bearing = frame_bearing;
    if (!speed_reference.empty())
        options.speed_reference = parse_speed_reference(speed_reference);
    if (!wind_text.empty())
    {
        SurfaceWind wind = parse_surface_wind(wind_text, unit);
        wind.roughness_length = roughness_length;
        wind.surface_layer_depth = surface_layer_depth;
        wind.stability = parse_stability(stability);
        wind.levels_per_decade = profile_levels;
        if (const auto problem = surface_wind_problem(wind, unit); !problem.empty())
            input_error("--wind: ", problem);
        if (wind.stability != Stability::neutral)
        {
            warn("--stability ", stability, ": ", describe(wind.stability), ". The stability classes "
                    "are unvalidated: on the 2017 windy day an assumed stable class was 0.3 dB rms "
                    "worse than the neutral log law. Prefer neutral unless the class is known.");
        }
        if (wind.height < 10.0 * wind.roughness_length)
        {
            warn("--wind: the anemometer at ", wind.height / unit.to_meters, " ", unit.name,
                    " stands within ten roughness lengths of the ground (z0 ",
                    wind.roughness_length / unit.to_meters, " ", unit.name, "), where the surface "
                    "layer's law does not yet hold; the profile carried up from it is uncertain.");
        }
        if (wind.speed > 0.0 && refraction != "none" && std::isnan(frame_bearing) && !frame)
        {
            warn("--wind: the direction is a compass bearing, turned into the local frame taking "
                 "+x as east. If this run's x/y are rotated -- the shipped AC tracks run +x along "
                 "210 degrees -- give --frame_bearing.");
        }
        options.surface_wind = wind;
    }
    if (!atmosphere_file.empty())
    {
        options.atmosphere = read_atmosphere_profile(atmosphere_file, unit, ground_elevation);
        if (options.atmosphere->has_pressure() && !std::isnan(pressure_override))
        {
            input_error("--atmosphere ", atmosphere_file, " gives pressure at every height; "
                    "--pressure would contradict it. (Without a pressure column, --pressure "
                    "sets the surface pressure the profile is integrated from.)");
        }
    }
    options.trim_start = trim_start;
    options.utc_offset = utc_offset;
    options.scattered = scattered;
    options.want_dnl = want_dnl;
    options.want_lamax = want_lamax;
    options.want_nax = want_nax;
    options.nax_threshold = nax_threshold;
    options.want_coverage = want_coverage;
    options.full_history = full_history;
    options.write_time_history = write_time_history;
    options.write_state = write_state;
    options.time_history_prefix = time_history_prefix;
    options.band_history_prefix = band_history_prefix;
    options.metric_file = metric_file;
    options.power_file = power_file;
    options.latlon = latlon;
    options.unit = unit;
    options.geodetic.frame = frame;
    options.geodetic.ground_elevation = ground_elevation;
    options.geodetic.start_utc = start_utc;
    // -s is the resampling step for both adapt() and a geodetic track that
    // needs resampling; neither ever coarsens a track that is already finer.
    options.geodetic.resample_dt = dt;
    options.geodetic.path_smoothing = path_smoothing;

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
            job.atmosphere_out = atmosphere_out_files.empty() ? std::string()
                                                              : atmosphere_out_files[path_index];
            job.count = path_count;
            trajectory_holder = process_trajectory(database, job, options, observers,
                    observer_set.microphone_heights, x_obs, y_obs, integrated_levels, dnl_energy,
                    lamax, nax, coverage);
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
        // DNL is the energy average over a 24-hour day: the summed event energies
        // over 86 400 s, not their raw sum (which reads 49.37 dB high).
        constexpr double seconds_per_day = 86400.0;
        std::vector<double> dnl(dnl_energy.size());
        for (std::size_t i = 0; i < dnl_energy.size(); i++)
            dnl[i] = 10.0 * std::log10(dnl_energy[i] / seconds_per_day);
        log("Writing day-night level over ", path_count,
                path_count == 1 ? " operation to " : " operations to ", dnl_file);
        write_levels(x_obs, y_obs, dnl, dnl_file, scattered, unit);
        if (latlon)
            write_geodetic_levels(x_obs, y_obs, dnl, *frame, latlon_companion(dnl_file), scattered);
    }

    if (want_lamax)
    {
        log("Writing LAmax over ", path_count,
                path_count == 1 ? " operation to " : " operations to ", lamax_file);
        write_levels(x_obs, y_obs, lamax, lamax_file, scattered, unit);
        if (latlon)
            write_geodetic_levels(x_obs, y_obs, lamax, *frame, latlon_companion(lamax_file), scattered);
    }

    if (want_nax)
    {
        log("Writing NAx (operations with LAmax at or above ", nax_threshold, " dBA) over ",
                path_count, path_count == 1 ? " operation to " : " operations to ", nax_file);
        write_levels(x_obs, y_obs, nax, nax_file, scattered, unit);
        if (latlon)
            write_geodetic_levels(x_obs, y_obs, nax, *frame, latlon_companion(nax_file), scattered);
    }

    if (want_coverage)
    {
        // The energy-weighted fraction of each observer's levels that came from
        // cells the database did not measure.  An observer that heard nothing
        // has no fraction to give, and is written as NaN.
        std::vector<double> share(coverage.size());
        for (std::size_t i = 0; i < coverage.size(); i++)
        {
            share[i] = coverage[i][1] > 0.0 ? coverage[i][0] / coverage[i][1]
                                            : std::numeric_limits<double>::quiet_NaN();
        }
        log("Writing the unmeasured share of the levels over ", path_count,
                path_count == 1 ? " operation to " : " operations to ", coverage_file);
        write_levels(x_obs, y_obs, share, coverage_file, scattered, unit);
        if (latlon)
            write_geodetic_levels(x_obs, y_obs, share, *frame, latlon_companion(coverage_file), scattered);
    }

    if (path_count > 1 && want_plots)
    {
        log("Plotting is skipped for multi-trajectory runs: one figure file cannot hold ",
                path_count, " flight paths.");
    }
    else if (want_plots && trajectory_holder)
    {
        draw_figures(*trajectory_holder, integrated_levels, x_obs, y_obs,
                scattered, plotting, figure_file, state_figure_file,
                contour_min, contour_max, contour_step, unit);
    }
    else if (want_plots)
    {
        // The only trajectory failed, so there is nothing to draw; dereferencing
        // the empty holder here was undefined behavior.
        warn("No figure drawn: the trajectory could not be processed.");
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
