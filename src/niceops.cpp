// Port of niceops.d.  Pipeline: load Database + Trajectory -> for each
// observer and each trajectory state, compute instantaneous dBA via
// Database::level_at_observer -> integrate the time history into a sound
// exposure level (SEL) using the 10-dB-down-from-peak window -> write CSV
// output and optional matplotlib figures.

#include <algorithm>
#include <array>
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
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>
#include "util/iota.hpp"
#include "vec3.hpp"

#include "plotting.hpp"
#include <complex>
#include <sstream>

#include "atmosphere.hpp"
#include "atmosphere_profile.hpp"
#include "contour.hpp"
#include "band_frames.hpp"
#include "pnl.hpp"
#include "refraction.hpp"
#include "stratified.hpp"
#include "ground_effect.hpp"
#include "manifest.hpp"
#include "operations.hpp"
#include "plate_receiver.hpp"
#include "propagation_table.hpp"
#include "raster.hpp"
#include "terrain.hpp"
#include "weighting.hpp"
#include "sphere_database.hpp"
#include "sphere_energy.hpp"
#include "state_condition.hpp"
#include "trajectory.hpp"
#include "util/csv.hpp"
#include "util/logger.hpp"
#include "util/paths.hpp"
#include "util/parallel.hpp"
#include "util/sel.hpp"
#include "util/time_weighting.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

// --manifest: run() fills it in, and main() writes it once run() has returned
// or failed, so the exit status it records is the one the process returns.
struct ManifestRequest
{
    std::string path;
    RunManifest manifest;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
};

std::optional<ManifestRequest>& manifest_request()
{
    static std::optional<ManifestRequest> request;
    return request;
}

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

// The name a database is cached and compared under: its canonical path, so
// two spellings of one file are one database.  A path that cannot be resolved
// keys as given; loading it then fails and names it.
std::string database_key(const std::string& path)
{
    std::error_code ec;
    const auto canonical = std::filesystem::weakly_canonical(path, ec);
    return ec ? path : canonical.string();
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

// One observer's EPNL (--epnl_file), from its 0.5 s reception-time frames.
struct EpnlRow
{
    // False when no frame has a finite PNLT -- the observer heard nothing above
    // the noy floor -- and the metric has no value.
    bool valid = false;
    // The result without its per-record histories, which on a grid would be
    // kept for every observer only to be thrown away.
    EpnlResult result;
    // The tone correction PNLTM carries (PNLTM - PNL at its record: the
    // record's own C, or the five-record mean when the band-sharing check
    // raised it), and the band its record's C came from, Hz; NaN for none.
    double tone = 0.0;
    double tone_band = 0.0;
    double pnltm_time = 0.0;
    double first_time = 0.0;
    double last_time = 0.0;
};

void write_epnl_file(const std::vector<EpnlRow>& rows, const std::vector<Vec3>& microphones,
        LengthUnit unit, const std::string& filename)
{
    info("Writing EPNL to ", filename);
    auto file = open_output(filename);
    file << std::setprecision(10);
    const std::string u = unit.name;
    file << "observer,x(" << u << "),y(" << u << "),z(" << u << "),EPNL(EPNdB),PNLTM(TPNdB),"
         << "duration_correction(dB),tone_correction(dB),tone_band(Hz),bandshare_adjusted,"
         << "PNLTM_time(s),duration_first(s),duration_last(s),clipped\n";
    for (std::size_t i = 0; i < rows.size(); i++)
    {
        const auto& m = microphones[i];
        file << i << "," << m.x / unit.to_meters << "," << m.y / unit.to_meters << ","
             << m.z / unit.to_meters << ",";
        const auto& row = rows[i];
        if (!row.valid)
        {
            // Written as silence rather than left out, so the rows stay one per
            // observer, in order.
            file << "-inf,-inf,nan,nan,nan,0,nan,nan,nan,0\n";
            continue;
        }
        const auto& r = row.result;
        file << r.epnl << "," << r.pnltm << "," << r.duration_correction << "," << row.tone << ","
             << row.tone_band << "," << (r.bandshare_adjusted ? 1 : 0) << "," << row.pnltm_time << ","
             << row.first_time << "," << row.last_time << "," << (r.clipped ? 1 : 0) << "\n";
    }
}

// One line per trajectory: the case, the energy mean and the maximum of its
// SELs, the fraction of observers at or above `threshold`, and on a grid the
// area at or above it, contoured as --contour_file contours it, in `unit`
// squared (NaN for scattered observers, which enclose nothing).
void append_metrics_file(const std::string& filename, const std::string& casename,
        const std::vector<double>& x_obs, const std::vector<double>& y_obs, bool scattered,
        const LengthUnit& unit, const std::vector<double>& levels, double threshold = 65.0)
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
    double enclosed = std::numeric_limits<double>::quiet_NaN();
    if (!scattered)
    {
        enclosed = 0.0;
        for (const auto& polygon : contour_polygons(x_obs, y_obs, levels, threshold))
            enclosed += polygon_area(polygon);
        enclosed /= unit.to_meters * unit.to_meters;
    }
    file << casename << "," << average << "," << maximum << "," << area << "," << enclosed << "\n";
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
    // --blend, --load_factor_correction and --min_load_factor: how a state's
    // condition picks its spheres, and whether the load-factor law also runs
    // below the reference (condition_selection.hpp).  The tip speed is filled
    // per trajectory from its vehicle.
    Database::BlendOptions blend;
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
    // --terrain, read once for the run: the ground models take each path's
    // mean ground line from it (terrain_ray), resampled per trajectory onto the
    // local frame every terrain_step meters with terrain_offset (m) added.
    // Grid observers already stand on it (run()).
    std::shared_ptr<const Terrain> terrain;
    std::string terrain_path;
    double terrain_offset = 0.0;
    double terrain_step = 0.0;
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
    // --shadow_model resolved (stratified only): each band's shadow boundary
    // from the near-ground profile it resolves at shadow_phase radians of
    // vertical phase (resolved_shadow_bands), and a cramped band's field from
    // cramped_ground (ShadowTransition::apply_resolved).  Off: the legacy
    // transition, byte for byte.
    bool resolved_shadow = false;
    double shadow_phase = default_resolved_phase;
    CrampedGround cramped_ground = CrampedGround::series;
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
    // Where the sphere's nose comes from (--nose).  Unset follows each
    // database's own azimuth_reference attribute.
    std::optional<NoseSource> nose;
    // Meters, after the --units conversion; NaN means the trim is off.
    double trim_start = 0.0;
    double utc_offset = 0.0;
    ClockWindow laeq_window;
    bool scattered = false;
    // -D, --lden_file, --cnel_file and the period levels, indexed by DayMetric.
    std::array<bool, day_metric_count> want_day_metric {};
    bool want_lamax = false;
    bool want_nax = false;
    double nax_threshold = std::numeric_limits<double>::quiet_NaN();
    // --ta_file: seconds the weighted level exceeds ta_threshold.
    bool want_ta = false;
    double ta_threshold = std::numeric_limits<double>::quiet_NaN();
    // The time weighting LAmax, NAx and TA read the level through
    // (--lamax_weighting), in reception time.
    LevelWeighting lamax_weighting = LevelWeighting::fast;
    // --csel_file: the C-weighted SEL, from the band levels.
    bool want_csel = false;
    bool want_coverage = false;
    bool want_provenance = false;
    bool full_history = false;
    // --epnl_file's tone correction: tones in bands below their noy threshold
    // excluded (pnl.hpp), unless --no-epnl_masking.
    ToneMasking epnl_masking = ToneMasking::below_noy_floor;
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
    std::string database;
    std::string vehicle;
    std::string output;
    std::string state;
    // --write_atmosphere: where the --wind profile goes; empty for nowhere.
    std::string atmosphere_out;
    // The C-weighted SEL's file (--csel_file); empty for none.
    std::string csel_output;
    // --epnl_file and --band_frames for this trajectory; empty for none.
    std::string epnl;
    std::string band_frames;
    // The run's name for this trajectory, which its per-track files are named
    // after: the track's stem.
    std::string name;
    // How many trajectories the run has (for naming only).
    std::size_t count = 1;
    // How many operations this trajectory stands for: what its contribution
    // to the run totals is multiplied by, after its own SEL is formed.
    double weight = 1.0;
    // Whether it came from an operations file (for the messages).
    bool operations = false;
    // Day, evening and night counts, in place of the track's clock: each day
    // metric then weights the plain SEL by period_count_weight (util/sel.hpp).
    std::optional<std::array<double, 3>> period_counts;
};

// What every trajectory of a run adds to, per observer: each day-evening-night
// metric's penalized energy (summed, normalized to the day at the end), the
// largest weighted level, the NAx count, the time above, and the coverage
// report's two energies.  Each trajectory's share is multiplied by its
// operations' weight, except the largest level, which no count raises.  A
// metric not asked for keeps an empty vector.
struct RunTotals
{
    std::array<std::vector<double>, day_metric_count> day_energy;
    std::vector<double> lamax;
    std::vector<double> nax;
    std::vector<double> time_above;
    std::vector<std::array<double, 2>> coverage;
    // --provenance_file, per observer: the SEL energy, and the shares of it
    // from states with each ConditionFlag (state_condition.hpp), in bit
    // order, and the energy each sphere supplied, sorted by
    // sphere, through the states' blend weights.
    std::vector<std::array<double, 1 + condition_flag_count>> provenance;
    std::vector<std::vector<std::pair<std::uint32_t, double>>> sphere_energy;
};

// --provenance_file: one row per observer, grid or scattered, of where its
// SEL energy came from over every operation.  The shares are of the energy
// the SEL integrated; the top sphere is the one the blend weights gave the
// most of it, with that share and its stored condition.  Spheres are credited
// by their blend weights as the blend gave them, so a sphere with nothing
// measured toward the observer, which the evaluation leaves out, still is.  An observer that
// heard nothing has no shares, and is written as NaN.
void write_provenance(const std::vector<Vec3>& observers, const RunTotals& totals, const Database& database,
        const std::string& out_file, LengthUnit unit)
{
    auto file = open_output(out_file);
    file << "x(" << unit.name << "),y(" << unit.name << "),exposure_db,share_mu_clamped,share_gamma_clamped,"
         << "share_ct_clamped,share_load_factor_held_max,share_load_factor_held_min,top_sphere,"
         << "top_sphere_share,top_sphere_mu,top_sphere_gamma,top_sphere_ct\n";
    file << std::setprecision(10);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i = 0; i < totals.provenance.size(); i++)
    {
        const double x = observers[i].x, y = observers[i].y;
        const auto& p = totals.provenance[i];
        // The exposure shared out: the SEL itself for one operation, the
        // weighted energy sum of the operations' SELs for several.
        file << x / unit.to_meters << "," << y / unit.to_meters << "," << 10.0 * std::log10(p[0]);
        for (int axis = 1; axis < 6; axis++)
            file << "," << (p[0] > 0.0 ? p[axis] / p[0] : nan);
        const auto& list = totals.sphere_energy[i];
        // The blend weights sum to one, so the spheres share out all of it.
        assert(std::abs(std::accumulate(list.begin(), list.end(), 0.0,
                               [](double sum, const auto& entry) { return sum + entry.second; })
                       - p[0])
                <= 1e-9 * p[0]);
        const auto top = std::max_element(list.begin(), list.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });
        if (top == list.end() || !(p[0] > 0.0))
        {
            file << ",nan,nan,nan,nan,nan\n";
            continue;
        }
        const auto& condition = database.condition(top->first);
        file << "," << top->first << "," << top->second / p[0] << "," << condition[0] << ","
             << condition[1] << "," << condition[2] << "\n";
    }
}

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

// The compass bearing of the local +x axis when the run knows it: --frame_bearing,
// or east when a geodetic frame was settled without one (it is then east/north).
// NaN when neither, which for an AC track means its x/y may be turned by any
// angle.
double known_frame_bearing(const RunOptions& options)
{
    if (!std::isnan(options.frame_bearing))
        return options.frame_bearing;
    return options.geodetic.frame ? 90.0 : std::numeric_limits<double>::quiet_NaN();
}

// Where one trajectory's sphere nose comes from, and the offset that turns its
// compass headings into the local frame (see heading_offset).
struct NoseSetup
{
    NoseSource source = NoseSource::velocity;
    std::optional<double> offset;
};

// Settle the nose: --nose, or what the database's azimuth_reference implies
// (heading for "heading", the velocity for "track" or nothing).  --nose against
// a database that records the other is a warning, since the sphere is then read
// a crab angle off the direction it was built along.
//
// Heading mode takes its offset from the frame and nowhere else.  heading_offset's
// fit, the fallback the hover nose uses when the run knows no frame, is the
// circular mean of heading against track: it folds the track's mean crab into
// the offset, and heading mode read through it would point the nose back along
// the mean track, undoing the very thing it is for.  So a run with no frame is
// refused.  A track with no heading at all is refused when --nose heading was
// asked for; when only the database implied it, the nose falls back to the
// velocity with a warning, there being nothing the run could be given instead.
NoseSetup settle_nose(const Trajectory& trajectory, const Database& database,
        const RunOptions& options, const std::string& trajectory_file)
{
    const auto& recorded = database.recorded_azimuth_reference();
    const AzimuthReference built = database.azimuth_reference();
    NoseSetup setup;
    setup.source = built == AzimuthReference::heading ? NoseSource::heading : NoseSource::velocity;
    if (options.nose)
    {
        const bool wants_heading = *options.nose == NoseSource::heading;
        if (recorded && wants_heading != (built == AzimuthReference::heading))
        {
            warn("--nose ", to_string(*options.nose), " overrides the database, whose sphere "
                    "azimuth is measured from the ", to_string(built), " (its azimuth_reference "
                    "attribute): the sphere is read turned by the crab angle from the direction "
                    "it was built along.");
        }
        else if (!recorded && wants_heading)
        {
            log("--nose heading on a database that does not record its azimuth reference, "
                "which is taken to be the track.");
        }
        setup.source = *options.nose;
    }
    const auto& states = trajectory.states;
    const double frame = known_frame_bearing(options);
    if (setup.source != NoseSource::heading)
    {
        setup.offset = heading_offset(states, frame, setup.source == NoseSource::ground_velocity);
        return setup;
    }
    const char* asked = options.nose ? "--nose heading" : "The database's azimuth_reference (heading)";
    const auto with_heading = std::count_if(states.begin(), states.end(),
            [](const State& s) { return std::isfinite(s.heading); });
    if (with_heading == 0)
    {
        if (options.nose)
        {
            input_error("--nose heading needs the nose's compass heading, and ", trajectory_file,
                    " carries none (the AC format's column 9, a WGS84 track's heading column "
                    "and a static operation's heading column are read). Use --nose velocity.");
        }
        warn("The database measures its sphere azimuth from the heading, and ", trajectory_file,
                " carries none (a WGS84 track needs a heading column): the nose follows the "
                "velocity instead, as if the aircraft flew with no crab. Give --nose velocity "
                "to say so.");
        setup.source = NoseSource::velocity;
        setup.offset = heading_offset(states, frame);
        return setup;
    }
    if (static_cast<std::size_t>(with_heading) != states.size())
    {
        const auto first = std::find_if(states.begin(), states.end(),
                [](const State& s) { return !std::isfinite(s.heading); });
        input_error(asked, " needs a heading at every state, and ", trajectory_file, " has ",
                states.size() - static_cast<std::size_t>(with_heading), " without one, the first "
                "at t = ", first->t, " s.");
    }
    if (std::isnan(frame))
    {
        input_error(asked, " turns each state's compass heading into the local frame, and "
                "needs that frame's bearing: give --frame_bearing, the compass bearing of +x "
                "(the shipped AC tracks run +x along 210 degrees)", options.nose ? "" : ", or --nose "
                "velocity", ". It is not inferred from the track, as the hover nose's offset is "
                "when no frame is given: that fit absorbs the mean crab (heading minus track), "
                "and the nose would come back along the mean track.");
    }
    setup.offset = frame * pi / 180.0;
    return setup;
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
    log("Atmosphere built from --wind: ", wind.speed, " m/s (", wind.speed / knot,
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
        RunTotals& totals)
{
    auto& lamax = totals.lamax;
    auto& nax = totals.nax;
    auto& coverage = totals.coverage;
    // Create trajectory from CSV data
    Trajectory trajectory(job.trajectory, job.vehicle, options.geodetic, options.unit);
    // Ground effect is on by default, but it acts on per-band levels, so it is
    // applied only where the spectral path runs on a database that has spectra.
    const bool ground_effect = options.ground_effect && options.spectral && database.has_spectra();

    // The vehicle file against what the database was built with
    // (check_vehicle, state_condition.hpp): refused past 2 percent.
    check_vehicle(trajectory.vehicle, database, job.vehicle, options.accept_vehicle_mismatch);

    // Trim trajectory (per Mike): crop from the given x to whichever comes
    // first, the lowest or the slowest point.  Off unless asked for -- it
    // changes which part of the flight is integrated, and so the levels.
    if (!std::isnan(options.trim_start))
    {
        trajectory.watts_trim(options.trim_start);
    }

    // Over flat ground the observers stand at z = 0, so a state below it is
    // underground: almost always an altitude datum that disagrees with the
    // run's ground.  It ran with exit status 0 and no word.  Checked on the
    // states the run integrates (after --trim), with a tolerance of the
    // smoothing a WGS84 track's fit is allowed (a track landing at the ground
    // dips that far), and at least a meter.  Over --terrain the datum may stand
    // above the ground, and z < 0 can be flown.
    if (!options.terrain)
    {
        const bool geodetic = trajectory.frame().has_value();
        const double tolerance = std::max(1.0, geodetic ? options.geodetic.path_smoothing : 0.0);
        std::size_t below = 0;
        const State* lowest = nullptr;
        for (const auto& state : trajectory.states)
        {
            if (state.z < -tolerance)
                below++;
            if (lowest == nullptr || state.z < lowest->z)
                lowest = &state;
        }
        if (below > 0)
            warn(below, " of ", trajectory.states.size(), " states of ", job.trajectory, " are below the "
                    "ground, the lowest at z = ", lowest->z / options.unit.to_meters, " ", options.unit.name,
                    " (t = ", lowest->t, " s). ",
                    geodetic ? "Check the track's altitude datum against --ground_elevation."
                             : "Check the track's z, which is read as height above the observers' ground.");
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
    // And the sphere's nose, which may be refused, before anything costs time.
    const NoseSetup nose_setup = settle_nose(trajectory, database, options, job.trajectory);

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

    // An exponential time weighting over fewer than two samples per time
    // constant reports the sampling, not the signal: on a 0.5 s history the
    // fast detector (tau 0.125 s) settles fully within every sample and its
    // "maximum" is the largest sample again, under another name.  Refused,
    // since the number would look like a meter reading.
    //
    // The check is on the integration step, the mean interval, because that is
    // the step every sample is weighted by (T times its reception factor), in
    // the time weighting as in the SEL.  Not on the longest interval: the 2017
    // AC tracks are 0.02 s with the odd dropout of up to 0.14 s, which the
    // longest interval would refuse while the history the detector sees is
    // uniform.  A part in 1e9 of slack lets a track at exactly tau/2 pass.
    const bool want_weighted_level = options.want_lamax || options.want_nax || options.want_ta;
    if (want_weighted_level && nt >= 2
            && integration_step > coarsest_weighting_step(options.lamax_weighting) * (1.0 + 1e-9))
    {
        input_error("--lamax_weighting ", to_string(options.lamax_weighting), " (", maximum_name(options.lamax_weighting),
                ", time constant ", weighting_time_constant(options.lamax_weighting), " s) needs the level history "
                "sampled at least every ", coarsest_weighting_step(options.lamax_weighting), " s, half its time "
                "constant, and ", job.trajectory, " is sampled every ", integration_step, " s. A WGS84 track "
                "is refitted onto a finer -s; the 27-column AC format is used as sampled, so it needs a finer "
                "file. Or choose another --lamax_weighting: slow needs 0.5 s, leq0.5 0.25 s, and instantaneous "
                "takes any step.");
    }

    // What this trajectory's energy counts for in each day metric: its
    // operation count, or, with day, evening and night counts, each period's
    // count times the period's weighting (period_count_weight).  It multiplies
    // the exposure after it is formed, never the levels, the window or the
    // bound.  A metric not asked for keeps a weight of 0 (no pass), and so
    // does LAeq, which has no period weighting: period counts with
    // --laeq_file are refused before any job runs.
    std::array<double, day_metric_count> day_weight {};
    for (std::size_t m = 0; m < day_metric_count; m++)
    {
        if (!options.want_day_metric[m])
            continue;
        const auto metric = static_cast<DayMetric>(m);
        if (!job.period_counts)
            day_weight[m] = job.weight;
        else if (metric != DayMetric::laeq)
            day_weight[m] = period_count_weight(metric, (*job.period_counts)[0], (*job.period_counts)[1],
                    (*job.period_counts)[2]);
    }

    // The day-evening-night penalties (-D, --lden_file, --cnel_file), per
    // state, from the track's UTC seconds shifted into local time.  Each is a
    // level added to the sample (0 outside its periods).  Period counts take
    // the clock's place, so they need none.
    std::array<std::vector<double>, day_metric_count> day_penalty;
    for (std::size_t m = 0; m < day_metric_count; m++)
    {
        if (!options.want_day_metric[m] || job.period_counts)
            continue;
        const auto metric = static_cast<DayMetric>(m);
        day_penalty[m].resize(nt);
        for (std::size_t it = 0; it < nt; it++)
        {
            if (std::isnan(trajectory.states[it].utcsec))
            {
                // A NaN fails both halves of a period's window test, so a track
                // with no clock was classed as daytime -- or, when the start
                // defaulted to zero, as UTC midnight, i.e. night: +10 dB on the
                // whole footprint, silently.  Either way the answer was made up.
                input_error(day_metric_option(metric), " needs the time of day, and ", job.trajectory,
                        " does not carry one (a geodetic track without a utcsec column, or a static hold). ",
                        job.operations ? "Give --start_utc, or local_time in the operations file."
                                       : "Give --start_utc <seconds since UTC midnight at t = 0>.");
            }
            day_penalty[m][it] = schedule_penalty(trajectory.states[it].utcsec, options.utc_offset, metric,
                    options.laeq_window);
        }
    }
    // A flight usually lies inside one period, so most metrics need no pass of
    // their own: with no penalty at any state the penalized SEL is the SEL
    // itself, bit for bit (x + 0 == x), and with every state outside a period
    // level's period (-inf) its energy is zero.
    enum class DayPass { integrate, the_sel, none };
    std::array<DayPass, day_metric_count> day_pass {};
    for (std::size_t m = 0; m < day_metric_count; m++)
    {
        if (job.period_counts)
        {
            // Every period's energy is the plain SEL, weighted by its count.
            day_pass[m] = day_weight[m] > 0.0 ? DayPass::the_sel : DayPass::none;
            continue;
        }
        const auto& penalty = day_penalty[m];
        if (std::all_of(penalty.begin(), penalty.end(), [](double p) { return p == 0.0; }))
            day_pass[m] = DayPass::the_sel;
        else if (std::all_of(penalty.begin(), penalty.end(), [](double p) { return p == -std::numeric_limits<double>::infinity(); }))
            day_pass[m] = DayPass::none;
    }

    // How each state's condition picks its spheres (--blend), and the
    // load-factor law's sides (--load_factor_correction).  The tip speed is
    // this trajectory's vehicle's, for a database that records none of its own.
    Database::BlendOptions blend_options = options.blend;
    blend_options.tip_speed = trajectory.vehicle.tip_speed;
    const bool two_sided = blend_options.load_factor_correction == LoadFactorCorrection::two_sided;
    if (two_sided && !database.fixed_load_factor())
    {
        input_error("--load_factor_correction two_sided applies a fixed_load_factor database's law, and the "
                "database given for ", job.trajectory, " is not one: its load factors are spheres, "
                "blended like any other condition.");
    }
    if (blend_options.policy != BlendPolicy::barycentric)
    {
        log("Sphere selection: --blend ", to_string(blend_options.policy),
                blend_options.policy == BlendPolicy::nearest
                        ? " (the nearest stored condition at every state)."
                        : " (AAM v3.1: bracketing in speed, then flight path angle, within 20 kt and 3/5 deg).");
        if (blend_options.policy == BlendPolicy::aam && !database.fixed_load_factor())
        {
            log("--blend aam: AAM has no thrust axis; each state searches the nearest of this "
                    "database's ", database.aam_layer_count(), " thrust-coefficient layers.");
            if (const std::size_t single = database.aam_single_sphere_layers(); single > 0)
            {
                warn("--blend aam: ", single, " of the database's ", database.aam_layer_count(),
                        " thrust-coefficient layers hold a single sphere, which a state nearest that "
                        "C_T takes whatever its speed and path angle. AAM's search needs the spheres "
                        "layered in C_T, not scattered.");
            }
        }
    }
    if (two_sided)
        log("Load-factor correction is two-sided, held at ", blend_options.min_load_factor, " g and ",
                options.max_load_factor, " g.");

    // The source geometry and the blend of database operating conditions depend
    // only on the trajectory state, so build them once here rather than once
    // per (observer, state) pair.  The blend is a Delaunay locate() walk
    // and was by far the largest single cost in the program; this turns
    // observers*states of them into states of them.
    // The sphere's nose: the velocity, or at hover the heading; or under --nose
    // heading the heading throughout (nose_direction, settle_nose).
    const auto& offset = nose_setup.offset;
    const NoseSource nose_source = nose_setup.source;
    std::size_t by_heading = 0, slow = 0, load_capped = 0, load_floored = 0;
    // Heading minus ground track over the states the track is defined at, as a
    // circular mean and an rms about it (heading mode's log).
    double crab_cos = 0.0, crab_sin = 0.0;
    std::vector<double> crabs;
    std::size_t mu_clamped = 0, fpa_clamped = 0, ct_clamped = 0, load_scaled = 0, load_lowered = 0;
    std::array<std::size_t, 4> spheres_blended {};
    std::vector<Database::SourceState> source_states;
    source_states.reserve(nt);
    // Per state, for --provenance_file: its ConditionFlag bits.
    std::vector<unsigned char> state_provenance(nt, 0);
    for (const auto& state : trajectory.states)
    {
        const auto velocity = nose_source == NoseSource::ground_velocity ? ground_velocity(state)
                                                                          : air_velocity(state);
        const auto nose = nose_direction(state, offset, nose_source);
        if (nose_source == NoseSource::heading)
        {
            if (std::hypot(state.vx, state.vy) >= velocity_nose_speed)
            {
                // Compass heading minus compass track: the track's local angle
                // less the nose's, the frame turning the other way.
                const double crab = std::remainder(std::atan2(state.vy, state.vx)
                        - std::atan2(nose[1], nose[0]), 2.0 * pi);
                crab_cos += std::cos(crab);
                crab_sin += std::sin(crab);
                crabs.push_back(crab);
            }
        }
        else if (std::hypot(velocity[0], velocity[1]) < velocity_nose_speed)
        {
            slow++;
            if (nose != velocity)
                by_heading++;
        }
        const auto condition = state_condition(state, database, options.qsam, options.max_load_factor,
                blend_options);
        if (condition.load_capped)
            load_capped++;
        if (condition.load_floored)
            load_floored++;
        source_states.push_back(database.source_state(Vec3(state.x, state.y, state.z),
                Vec3(nose[0], nose[1], nose[2]),
                Vec3(state.nx, state.ny, state.nz), condition.advance_ratio,
                condition.flight_path_angle, condition.thrust_coefficient, blend_options));
        const auto& blend = source_states.back().blend;
        state_provenance[source_states.size() - 1]
                = static_cast<unsigned char>(condition_flags(condition, blend));
        mu_clamped += blend.advance_ratio_clamped;
        fpa_clamped += blend.flight_path_angle_clamped;
        ct_clamped += blend.thrust_coefficient_clamped;
        load_scaled += blend.level_correction_db > 0.0;
        load_lowered += blend.level_correction_db < 0.0;
        if (blend_options.policy == BlendPolicy::aam)
        {
            std::size_t used = 0;
            for (const double weight : blend.weights)
                used += weight > 0.0;
            spheres_blended[std::max<std::size_t>(used, 1) - 1]++;
        }
    }
    if (blend_options.policy == BlendPolicy::aam)
    {
        // How often the gates let AAM interpolate at all: a state on one
        // sphere is either on a stored condition or past every gate.
        log("--blend aam over ", nt, " states: 4 spheres at ", spheres_blended[3], ", 3 at ",
                spheres_blended[2], ", 2 at ", spheres_blended[1], ", 1 at ", spheres_blended[0], ".");
    }
    log("Source axis coverage over ", nt, " states: advance ratio clamped ", mu_clamped,
            ", flight-path angle clamped ", fpa_clamped, ", thrust coefficient clamped ", ct_clamped,
            "; analytic load-factor correction applied ", load_scaled,
            two_sided ? " (raising), " : ".", two_sided ? std::to_string(load_lowered) : std::string(),
            two_sided ? " (lowering)." : "");
    if (load_floored > 0)
    {
        warn(load_floored, " of ", nt, " states are below ", blend_options.min_load_factor, " g on a "
                "fixed_load_factor database; the two-sided load-factor correction is held at ",
                blend_options.min_load_factor, " g (--min_load_factor).");
    }
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
    if (nose_source == NoseSource::heading)
    {
        std::ostringstream crab_text;
        if (!crabs.empty())
        {
            const double mean = std::atan2(crab_sin, crab_cos);
            double sum = 0.0;
            for (const double crab : crabs)
            {
                const double d = std::remainder(crab - mean, 2.0 * pi);
                sum += d * d;
            }
            const double mean_degrees = mean * 180.0 / pi;
            crab_text << " Over the " << crabs.size() << " states above 8 kt it sits "
                      << mean_degrees << " degrees from the ground track (heading minus track, "
                      << "positive with the nose right of the track), rms "
                      << std::sqrt(sum / static_cast<double>(crabs.size())) * 180.0 / pi
                      << " degrees about that.";
            // A helicopter crabs by the wind's share of its airspeed; tens of
            // degrees on average over a whole track is a frame turned the wrong
            // way, not a crosswind.
            if (std::fabs(mean_degrees) > 45.0)
            {
                warn("--nose heading: the heading sits ", mean_degrees, " degrees from the ground "
                        "track on average over ", job.trajectory, ". Is --frame_bearing (",
                        *offset * 180.0 / pi, ") the bearing of this track's +x axis?");
            }
        }
        log("The sphere's nose is the compass heading at all ", nt, " states, turned into the "
            "local frame by a frame bearing of ", *offset * 180.0 / pi, " degrees.", crab_text.str());
    }
    else if (slow > 0)
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
                    : (std::filesystem::path(options.power_file).parent_path() / (job.name + "_power.csv")).string();
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
    // The C-weighted SEL is formed from the same per-band levels the band
    // history writes, so both record them.
    const bool record_bands = write_band_history || options.want_csel;
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
            : "_" + job.name;

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

    // --terrain: the ground under this trajectory's scene, on the run's frame,
    // and the checks that the track and the terrain mean the same place and the
    // same heights.  The observers already stand on it (run()).
    std::optional<TerrainRaster> terrain_raster;
    if (options.terrain)
    {
        const LocalFrame& frame = *options.geodetic.frame;
        const LengthUnit& unit = options.unit;
        const auto& recorded = trajectory.recorded;
        if (!recorded.empty())
        {
            // An AC state file records each state's latitude and longitude beside
            // its x/y, so the run's frame can be checked against the track's own:
            // with a wrong --origin or --frame_bearing the terrain would be read
            // under some other place.  The demo track agrees with its own frame
            // to 0.85 ft rms, the 2017 Amedee tracks with theirs to 1-5 ft.
            const std::size_t stride = std::max<std::size_t>(1, recorded.size() / 500);
            double squares = 0.0;
            std::size_t count = 0;
            for (std::size_t k = 0; k < recorded.size(); k += stride)
            {
                const Vec3 p = frame.to_local(Geodetic {recorded[k].latitude, recorded[k].longitude, 0.0});
                squares += (p.x - recorded[k].x) * (p.x - recorded[k].x) + (p.y - recorded[k].y) * (p.y - recorded[k].y);
                count++;
            }
            const double rms = std::sqrt(squares / static_cast<double>(count));
            constexpr double refused_rms = 30.0; // m
            constexpr double warned_rms = 3.0;   // m
            if (rms > refused_rms)
            {
                input_error("--terrain: ", job.trajectory, "'s own latitudes and longitudes put its states ",
                        rms / unit.to_meters, " ", unit.name, " rms from its x/y in this run's frame (--origin ",
                        frame.reference().latitude, " ", frame.reference().longitude, ", --frame_bearing ",
                        std::isnan(options.frame_bearing) ? std::string("unset, +x east")
                                                          : std::to_string(options.frame_bearing),
                        "), so the terrain would be read under the wrong place. Give the origin and bearing "
                        "its x/y were made with.");
            }
            const std::string frame_note = std::string("  the track's latitudes and longitudes sit ")
                    + std::to_string(rms / unit.to_meters) + " " + unit.name + " rms from its x/y in this run's frame.";
            if (rms > warned_rms)
                warn(frame_note, " Check --origin and --frame_bearing.");
            else
                log(frame_note);
            // And its altitude: z = 0 is, on the track's own datum, at alt - z,
            // which is what --ground_elevation must be for the terrain's heights
            // (put on that datum by --terrain_offset) to land at the right z.
            // Checked only where alt - z is one number along the track.
            std::vector<double> datum(recorded.size());
            for (std::size_t k = 0; k < recorded.size(); k++)
                datum[k] = recorded[k].altitude - recorded[k].z;
            const auto [low, high] = std::minmax_element(datum.begin(), datum.end());
            std::nth_element(datum.begin(), datum.begin() + datum.size() / 2, datum.end());
            const double median = datum[datum.size() / 2];
            constexpr double datum_tolerance = 1.0; // m
            if (*high - *low <= datum_tolerance)
            {
                if (std::fabs(median - options.geodetic.ground_elevation) > datum_tolerance)
                {
                    input_error("--terrain: ", job.trajectory, " puts its z = 0 at altitude ",
                            median / unit.to_meters, " ", unit.name, " (its alt less its z), but "
                            "--ground_elevation is ", options.geodetic.ground_elevation / unit.to_meters, " ",
                            unit.name, ". The terrain's heights would be ",
                            (median - options.geodetic.ground_elevation) / unit.to_meters, " ", unit.name,
                            " off against the track. Give --ground_elevation ", median / unit.to_meters, ".");
                }
                log("  the track's z = 0 is at altitude ", median / unit.to_meters, " ", unit.name,
                        " on its own datum, which --ground_elevation matches.");
            }
            else
            {
                log("  the track's alt less its z varies by ", (*high - *low) / unit.to_meters, " ", unit.name,
                        ", so --ground_elevation is not checked against it.");
            }
        }

        // The raster covers every observer and every state's footprint, with a
        // margin for the rays GroundPlaneCache samples: a step past the far end
        // of a path, and a bearing step to either side of it.
        double x_low = std::numeric_limits<double>::infinity(), x_high = -x_low;
        double y_low = x_low, y_high = -x_low;
        for (const auto& o : observers)
        {
            x_low = std::min(x_low, o.x); x_high = std::max(x_high, o.x);
            y_low = std::min(y_low, o.y); y_high = std::max(y_high, o.y);
        }
        for (const auto& state : source_states)
        {
            x_low = std::min(x_low, state.position.x); x_high = std::max(x_high, state.position.x);
            y_low = std::min(y_low, state.position.y); y_high = std::max(y_high, state.position.y);
        }
        const double step = options.terrain_step;
        const double margin = 2.0 * step
                + std::hypot(x_high - x_low, y_high - y_low) * (2.0 * pi / GroundPlaneCache::bearings);
        const auto raster_start = std::chrono::steady_clock::now();
        terrain_raster.emplace(*options.terrain, frame, x_low - margin, x_high + margin, y_low - margin,
                y_high + margin, step, options.terrain_offset - options.geodetic.ground_elevation);
        log("  terrain sampled onto the local frame: ", terrain_raster->nodes(), " points every ",
                step / unit.to_meters, " ", unit.name, " in ",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - raster_start).count(), " s.");

        // A state under the ground is a datum mistake, not a flight: missing
        // --ground_elevation puts the terrain a whole site elevation above every
        // state.  A few meters is the terrain's own error (3DEP's is 1-2 m) or a
        // landing, and is let through: such a state's height above a path's mean
        // ground is taken as zero where it falls below it (terrain_ray).
        constexpr double refused_depth = 10.0; // m
        std::size_t below = 0;
        double deepest = 0.0;
        double clearance = std::numeric_limits<double>::infinity();
        for (const auto& state : source_states)
        {
            const double depth = terrain_raster->height(state.position.x, state.position.y) - state.position.z;
            clearance = std::min(clearance, -depth);
            if (depth > 0.0)
            {
                below++;
                deepest = std::max(deepest, depth);
            }
        }
        if (deepest > refused_depth)
        {
            input_error("--terrain: ", below, " of ", nt, " states of ", job.trajectory, " are below the "
                    "terrain under them, the deepest by ", deepest / unit.to_meters, " ", unit.name,
                    ". The trajectory and the terrain are on different vertical datums: the ground's z is "
                    "the terrain height + --terrain_offset - --ground_elevation, which must be the "
                    "track's own height of z = 0 on the terrain's datum.");
        }
        if (below > 0)
        {
            warn(below, " of ", nt, " states of ", job.trajectory, " are up to ", deepest / unit.to_meters, " ",
                    unit.name, " below the terrain under them, within the terrain's own error; where a "
                    "source falls below a path's mean ground its height above it is taken as zero.");
        }
        else
        {
            log("  the trajectory clears the terrain under it by ", clearance / unit.to_meters, " ",
                    unit.name, " at least.");
        }
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
    // --csel_file's C-weighting, over the same --band_limits as band_gains.
    std::vector<double> c_gains;
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
    // --epnl_file: where each of the 24 PNL bands sits among the database's.
    std::array<std::size_t, pnl_band_count> pnl_bands {};
    // --propagation_table: where the run's bands sit in it, and, for a table
    // that already includes absorption, the zero coefficients that keep the
    // run from taking it twice.
    PropagationTable::BandMap table_bands;
    std::vector<double> table_absorption;
    // K for plate microphones read from the table, laid out
    // [(elevation * table_plate_azimuths + azimuth) * bands + band], on
    // elevations k * table_plate_step from 0 to 90 degrees.
    constexpr int table_plate_azimuths = 36;
    std::vector<double> table_plate_offset;
    double table_plate_step = 0.0;
    std::size_t table_plate_elevations = 0;

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
        if (options.want_csel)
        {
            // C-weighting is a weighting of the spectrum: a broadband database
            // has one A-weighted level and nothing to reweight.
            input_error("--csel_file needs a database carrying per-band spectra, to C-weight; this "
                        "one has none.");
        }
        if (!job.epnl.empty() || !job.band_frames.empty())
        {
            input_error(!job.epnl.empty() ? "--epnl_file" : "--band_frames",
                    " needs a database carrying per-band spectra; this one has none.");
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
        if (options.want_csel)
            c_gains = c_weighting_gains(bands, options.band_limits);
        while (first_band < bands.size() && band_gains[first_band] == 0.0)
            first_band++;
        last_band = bands.size();
        while (last_band > first_band && band_gains[last_band - 1] == 0.0)
            last_band--;

        // EPNL is defined over the 24 one-third-octave bands 50 Hz to 10 kHz
        // (A36.4), and the noy tables and the tone correction need every one of
        // them: a missing band is not a quiet band, and leaving it out would
        // move the tone correction's background and the noy sum alike.  Bands
        // the database has outside them (10 to 40 Hz) are ignored, as the
        // regulation ignores them.
        if (!job.epnl.empty())
        {
            std::vector<double> missing;
            pnl_bands = pnl_band_indices(bands, missing);
            if (!missing.empty())
            {
                std::ostringstream list;
                for (std::size_t i = 0; i < missing.size(); i++)
                    list << (i ? ", " : "") << missing[i];
                input_error("--epnl_file needs the 24 one-third-octave bands from 50 Hz to 10 kHz "
                            "(14 CFR 36 A36.4); the database (bands ", bands.front(), " to ",
                        bands.back(), " Hz) lacks ", list.str(), " Hz.");
            }
            for (const std::size_t band : pnl_bands)
            {
                if (band < first_band || band >= last_band)
                {
                    input_error("--epnl_file needs the bands from 50 Hz to 10 kHz evaluated, and "
                                "--band_limits ", options.band_limits.low, " to ",
                            options.band_limits.high, " Hz leaves out the ", bands[band],
                            " Hz band. Widen it to take in at least 50 to 10000 Hz.");
                }
            }
            log(options.epnl_masking == ToneMasking::below_noy_floor
                            ? "EPNL: tones in bands below their noy threshold are excluded from the tone "
                              "correction (--no-epnl_masking for the procedure as written)."
                            : "EPNL: the tone correction as written, every band eligible "
                              "(--no-epnl_masking).");
        }

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
                if (terrain_raster)
                {
                    // Over terrain the geometry is taken over each path's mean
                    // ground, which can lie below z = 0: reach down to the lowest
                    // observer, and past it by the observers' relief and the
                    // tallest microphone, which bound how far below a microphone
                    // its foot on a mean ground line can be.
                    double z_low = std::numeric_limits<double>::infinity();
                    double z_high = -z_low;
                    for (const auto& o : observers)
                    {
                        z_low = std::min(z_low, o.z);
                        z_high = std::max(z_high, o.z);
                    }
                    const double tallest = *std::max_element(microphone_height.begin(), microphone_height.end());
                    for (const auto& state : source_states)
                        for (const double cx : {lo_x, hi_x})
                            for (const double cy : {lo_y, hi_y})
                                r_max = std::max(r_max, std::hypot(std::hypot(cx - state.position.x,
                                        cy - state.position.y), state.position.z - z_low)
                                                + (z_high - z_low) + tallest);
                }
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
            // A flush microphone reads it at its lowest receiver height (below),
            // so only point receivers are checked.
            std::vector<double> uncovered;
            for (std::size_t i = 0; i < microphone_height.size(); i++)
            {
                const double height = microphone_height[i];
                if (!on_plate[i] && !table.covers_receiver(height)
                        && std::find(uncovered.begin(), uncovered.end(), height) == uncovered.end())
                    uncovered.push_back(height);
            }
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
            // receiver height, plus K, the plate over bare ground there, per
            // band, azimuth and elevation, on a 1 km sphere as ShadowTransition
            // places a plate.  K was once taken at 2 degrees for every pair, as
            // if it barely depended on elevation: through a flush row that put
            // the demo footprint up to 0.34 dB low against the plate run.  Over
            // a receiver above the ground the bare term's interference fringes
            // are lambda / (2 h) apart in the sine of the elevation, so the grid
            // takes eight steps per fringe at the highest band (clamped), and K
            // is interpolated linearly in elevation.  A raised row is still
            // only an approximation of the plate, growing with the row's height
            // and not from the grid, so it is warned about.
            if (plate)
            {
                GroundModel bare = options.ground;
                const double lowest = table.axes().receiver_heights.front();
                if (lowest > 0.0)
                    warn("Plate microphones read --propagation_table ", options.propagation_table_path,
                            " at its lowest receiver height, ", lowest / options.unit.to_meters, " ",
                            options.unit.name, ", plus the plate's offset over bare ground there, which "
                            "approximates the plate less well the higher the row. Give the table a "
                            "receiver height of 0 for plate microphones.");
                bare.receiver_height = lowest;
                const double shortest_wavelength = sound_speed / bands[last_band - 1];
                const double fringe = lowest > 0.0 ? shortest_wavelength / (2.0 * lowest) : 1.0;
                table_plate_step = std::clamp(fringe / 8.0, 0.05 * pi / 180.0, 0.5 * pi / 180.0);
                table_plate_elevations = static_cast<std::size_t>(std::ceil(0.5 * pi / table_plate_step)) + 1;
                table_plate_step = 0.5 * pi / static_cast<double>(table_plate_elevations - 1);
                table_plate_offset.assign(table_plate_elevations * table_plate_azimuths * bands.size(), 0.0);
                std::vector<double> bare_db(bands.size(), 0.0), plate_db(bands.size(), 0.0);
                for (std::size_t e = 0; e < table_plate_elevations; e++)
                {
                    // Exact grazing has no geometry to speak of: a hundredth of
                    // a degree stands in for it.
                    const double elevation = std::max(static_cast<double>(e) * table_plate_step, 1e-2 * pi / 180.0);
                    const double h = 1000.0 * std::sin(elevation), d = 1000.0 * std::cos(elevation);
                    excess_ground_attenuation(bare, ground_impedance_bands, h, d, sound_speed, bands,
                            bare_db.data(), first_band, last_band, false);
                    for (int i = 0; i < table_plate_azimuths; i++)
                    {
                        plate->level(h, d, 2.0 * pi * i / table_plate_azimuths, plate_db.data(), first_band,
                                last_band);
                        double* k = &table_plate_offset[(e * table_plate_azimuths + i) * bands.size()];
                        for (std::size_t band = first_band; band < last_band; band++)
                            k[band] = plate_db[band] - bare_db[band];
                    }
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
    // --shadow_model resolved: per receiver (the point heights, then the plate).
    std::optional<ResolvedShadowTable> resolved_table;
    // --refraction stratified: the profile's layers, with the linear model kept
    // alongside for its warm start and for what the layered solve leaves to it.
    std::optional<LayeredSoundSpeed> layered_sound_speed;
    if (options.refraction && spectral && profile)
    {
        const double bearing = std::isnan(options.frame_bearing) ? 90.0 : options.frame_bearing;
        effective_sound_speed = EffectiveSoundSpeed(*profile, bearing, z_low,
                z_high + EffectiveSoundSpeed::min_fit_span);
        // The layered solve follows flat air at a path's fastest height only
        // as far as the highest band evaluated resolves it (stratified.hpp).
        if (options.stratified)
            layered_sound_speed = LayeredSoundSpeed(*profile, bearing, database.band_centers()[last_band - 1]);
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
        // --shadow_model resolved: the cramped bands' field on each transition,
        // and every band's boundary per bearing, source height and receiver.
        if (options.resolved_shadow && ground_model != nullptr && layered_sound_speed)
        {
            const auto start = std::chrono::steady_clock::now();
            std::vector<double> receivers;
            std::vector<bool> flush_receivers;
            for (std::size_t k = 0; k < shadow_transitions.size(); k++)
            {
                shadow_transitions[k].set_cramped_ground(options.cramped_ground, &point_models[k],
                        options.thread_count);
                receivers.push_back(point_heights[k]);
                flush_receivers.push_back(point_heights[k] < ShadowTransition::flush_height);
            }
            if (plate_transition)
            {
                plate_transition->set_cramped_ground(options.cramped_ground, nullptr, options.thread_count);
                receivers.push_back(plate->thickness());
                flush_receivers.push_back(true);
            }
            const double cramped_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            double state_low = std::numeric_limits<double>::infinity(), state_high = 0.0;
            for (const auto& state : source_states)
            {
                state_low = std::min(state_low, state.position.z);
                state_high = std::max(state_high, state.position.z);
            }
            resolved_table.emplace(*layered_sound_speed, *effective_sound_speed, database.band_centers(),
                    first_band, last_band, options.shadow_phase, std::max(state_low, 0.0), state_high,
                    receivers, flush_receivers, options.thread_count);
            log("  --shadow_model resolved: each band's shadow boundary from the profile it resolves at ",
                    options.shadow_phase, " rad of vertical phase, tabulated over bearing and source "
                    "height", options.cramped_ground == CrampedGround::series
                            ? ", and the residue series at the source's own height for cramped bands,"
                            : options.cramped_ground == CrampedGround::straight
                            ? ", with cramped bands' lit side on the straight path's ground geometry,"
                            : ",",
                    " in ", std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(),
                    " s (", cramped_seconds, " s of it the cramped bands' tables).");
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
        if (options.stratified)
            log("  a ray along flat air at a path's fastest height is followed only as far as "
                "the highest band evaluated, ", layered_sound_speed->resolution_frequency(),
                    " Hz, resolves it.");
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
    const bool any_day_metric = std::any_of(options.want_day_metric.begin(),
            options.want_day_metric.end(), [](bool b) { return b; });
    // --epnl_file and --band_frames need every state's bands, in reception-time
    // frames.  The EPNL's window is found on PNLT, which the A-weighted bound
    // says nothing about, so the windowed search cannot serve them.
    const bool want_epnl = !job.epnl.empty();
    const bool want_band_frames = !job.band_frames.empty();
    const bool want_frames = want_epnl || want_band_frames;
    if (want_frames && (!spectral || add_ambient))
    {
        // run() refuses both up front; this is the backstop.
        input_error(want_epnl ? "--epnl_file" : "--band_frames",
                " needs the spectral path and no -a: the ambient is an A-weighted level with no "
                "spectrum to add to the bands.");
    }
    const bool needs_full_history = options.write_time_history || any_day_metric
            || options.full_history || options.want_lamax || options.want_nax || options.want_ta
            || options.want_csel || write_band_history || want_frames;
    // run() refuses --provenance_file with -a.
    assert(!(options.want_provenance && level_is_transformed));
    // The coverage report needs every state's level, so nothing is skipped.
    const bool bounded = !level_is_transformed && !options.want_coverage && !options.want_provenance;
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
            // A plate microphone reading the table adds its K, interpolated
            // between nodes and so never above the largest of them.
            if (options.propagation_table)
            {
                auto table_gains = options.propagation_table->max_power_ratio(table_bands);
                if (!table_plate_offset.empty())
                {
                    for (std::size_t band = 0; band < table_gains.size(); band++)
                    {
                        double largest = -std::numeric_limits<double>::infinity();
                        for (std::size_t node = 0; node < table_plate_elevations * table_plate_azimuths; node++)
                            largest = std::max(largest, table_plate_offset[node * table_gains.size() + band]);
                        table_gains[band] *= std::pow(10.0, 0.1 * std::max(largest, 0.0));
                    }
                }
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
    //     --epnl_file / --band_frames  every state's bands, into frames
    //     -D / --dnlfile the day-night window is found on a penalized history,
    //                    which a window found on the raw one does not cover;
    //                    --lden_file and --cnel_file the same
    //     --lamax_file, --nax_file, --ta_file
    //                    the time weighting runs over the whole history
    //     --csel_file    the C-weighted window is not the A-weighted one, and
    //                    every state's bands are needed to find it
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
    std::atomic<std::size_t> terrain_samples {0};
    std::atomic<std::size_t> terrain_raised {0};
    // Which observers had a height raised onto their path's mean ground plane
    // (each observer is one worker's, so no synchronization).
    std::vector<char> raised_at(terrain_raster ? observers.size() : 0, 0);
    // The table's bearings are compass bearings; the run's x/y are turned from
    // east/north by --frame_bearing (CompassToLocal), whose inverse this is.
    const CompassToLocal compass(std::isnan(options.frame_bearing) ? 90.0 : options.frame_bearing);

    // This trajectory's C-weighted SEL per observer (--csel_file).
    std::vector<double> csel_levels(options.want_csel ? observers.size() : 0);
    // Per observer, each written only by the worker that owns the observer.
    std::vector<EpnlRow> epnl_rows(want_epnl ? observers.size() : 0);
    std::vector<ObserverFrames> observer_frames(want_band_frames ? observers.size() : 0);

    const auto acoustics_start = std::chrono::steady_clock::now();

    parallel_chunks(observers.size(), workers, [&](std::size_t begin, std::size_t end) {
        // One scratch history per chunk rather than per observer; parallel.hpp
        // deals out several chunks per worker, so this runs a few times a
        // thread, which it measures as noise against the work in the chunk.
        std::vector<double> levels(nt);
        std::vector<double> sel_weight;
        std::vector<double> ranges(nt);
        std::vector<double> penalized;
        // Each sample's reception-time duration, for the time weighting.
        std::vector<double> durations;
        // --csel_file: the C-weighted history.
        std::vector<double> c_levels(options.want_csel ? nt : 0);
        // --band_history, --csel_file, --epnl_file and --band_frames: every
        // state's band levels, and (for the history) the geometry each was
        // reached along.  The EPNL and frame outputs keep the levels for the
        // frames.
        const bool keep_bands = record_bands || want_frames;
        const std::size_t history_bands = keep_bands ? database.band_centers().size() : 0;
        std::vector<double> band_levels(nt * history_bands);
        std::vector<BandHistoryRow> band_rows(write_band_history ? nt : 0);
        BandFrames frames(want_frames ? last_band - first_band : 0, pnl_record_seconds);
        std::vector<PnlSpectrum> pnl_history;
        auto record_path = [&](std::size_t it, const RayPath& path) {
            if (!keep_bands)
                return static_cast<double*>(nullptr);
            if (!write_band_history)
                return band_levels.data() + it * history_bands;
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
        // --terrain: each observer's mean ground lines (GroundPlaneCache), and
        // how many pairs took one, and how many had a height raised to zero.
        std::optional<GroundPlaneCache> planes;
        if (terrain_raster)
            planes.emplace(*terrain_raster);
        std::size_t terrain_here = 0, raised_here = 0;
        std::vector<double> shadow_weight(shadow.size(), 1.0);
        std::vector<ResolvedBand> resolved_bands(resolved_table ? database.band_centers().size() : 0);
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
            // Over terrain the ground model's geometry is taken over each path's
            // mean ground line, from this observer's foot toward the source's.
            const bool over_terrain = planes.has_value() && ground_here != nullptr;
            if (over_terrain)
                planes->set_observer(observers[index].x, observers[index].y);
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
                // K at the pair's azimuth in the plate's frame and its elevation,
                // bilinear between the four nodes around it.
                const double* k[4] = {nullptr, nullptr, nullptr, nullptr};
                double w[4] = {0.0, 0.0, 0.0, 0.0};
                if (flush && !table_plate_offset.empty())
                {
                    const double azimuth = std::atan2(separation.y, separation.x) - plate_mic_azimuth + 0.5 * pi;
                    double x = std::fmod(azimuth / (2.0 * pi) * table_plate_azimuths, double(table_plate_azimuths));
                    if (x < 0.0)
                        x += table_plate_azimuths;
                    const std::size_t a_lo = std::min(static_cast<std::size_t>(x), std::size_t(table_plate_azimuths - 1));
                    const std::size_t a_hi = (a_lo + 1) % table_plate_azimuths;
                    const double a_weight = x - static_cast<double>(a_lo);
                    const double elevation = std::atan2(std::max(source.z, 0.0), separation.horizontal_norm());
                    const double y = std::min(elevation / table_plate_step, double(table_plate_elevations - 1));
                    const std::size_t e_lo = std::min(static_cast<std::size_t>(y), table_plate_elevations - 2);
                    const double e_weight = y - static_cast<double>(e_lo);
                    const auto node = [&](std::size_t e, std::size_t a) {
                        return &table_plate_offset[(e * table_plate_azimuths + a) * table_loss.size()];
                    };
                    k[0] = node(e_lo, a_lo);
                    k[1] = node(e_lo, a_hi);
                    k[2] = node(e_lo + 1, a_lo);
                    k[3] = node(e_lo + 1, a_hi);
                    w[0] = (1.0 - e_weight) * (1.0 - a_weight);
                    w[1] = (1.0 - e_weight) * a_weight;
                    w[2] = e_weight * (1.0 - a_weight);
                    w[3] = e_weight * a_weight;
                }
                for (std::size_t band = first_band; band < last_band; band++)
                {
                    double gain = table_loss[band];
                    if (k[0] != nullptr)
                        gain += w[0] * k[0][band] + w[1] * k[1][band] + w[2] * k[2][band] + w[3] * k[3][band];
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
                if (over_terrain)
                {
                    const Vec3& source = source_states[it].position;
                    const double dx = source.x - observer.x;
                    const double dy = source.y - observer.y;
                    bool raised = false;
                    const RayPath path = terrain_ray(source, observer,
                            planes->line_to(dx, dy, std::sqrt(dx * dx + dy * dy)), flush, &raised);
                    terrain_here++;
                    raised_here += raised;
                    if (raised)
                        raised_at[index] = 1;
                    return database.spectral_level_along(path, source_states[it], absorption_for(it, observer),
                            band_gains, ground_here, &ground_impedance_bands, sound_speed, first_band, last_band,
                            packed_of_state[it], plate_receiver, plate_mic_azimuth, table_here, nullptr, nullptr,
                            bands_out);
                }
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
                    double* table_bands_out = keep_bands
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
                // --shadow_model resolved: this pair's bands from the table, for
                // a ray that carries a boundary (which apply_resolved needs).
                const bool resolved_here = resolved_table && transition_here != nullptr && ray.curved
                        && ray.direct.boundary_gradient < 0.0 && std::isfinite(ray.direct.boundary);
                if (resolved_here)
                {
                    const double hx = ray.path.separation.x / ray.horizontal;
                    const double hy = ray.path.separation.y / ray.horizontal;
                    resolved_table->lookup(hx, hy, source_states[it].position.z,
                            flush ? shadow_transitions.size() : point_of[index], resolved_bands.data());
                }
                if (transition_here != nullptr
                        && (resolved_here
                                        ? transition_here->apply_resolved(ray, resolved_bands.data(),
                                                  shadow_attenuation, options.shadow_floor, first_band,
                                                  last_band, shadow_weight.data(), shadow.data(), plate_azimuth)
                                        : transition_here->apply(ray, shadow_attenuation, options.shadow_floor,
                                                  first_band, last_band, shadow_weight.data(), shadow.data(),
                                                  plate_azimuth)))
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
                    double* bands_out = keep_bands
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
            // With the SEL, for the reports that share it out (--coverage_file,
            // --provenance_file), what it integrated per state: its 10 dB-down
            // window (the whole history under --full_history), each sample's
            // energy times the step and its share of reception time, and zero
            // outside.  Summing every sample alike counted states the SEL
            // leaves out, and across operations weighted each by its sample
            // count rather than its exposure.
            if (options.want_coverage || options.want_provenance)
            {
                integrated_levels[index] = sound_exposure_level_by_sample(levels, integration_step, ranges,
                        sound_speed, options.full_history, sel_weight);
            }
            else
            {
                integrated_levels[index] = options.full_history
                        ? full_sound_exposure_level(levels, integration_step, ranges, sound_speed)
                        : sound_exposure_level(levels, integration_step, ranges, sound_speed);
            }
            if (options.want_csel)
            {
                // The C-weighted level of every state, from the bands the
                // A-weighted one was summed from, and its SEL by the same rule
                // as SEL's: the 10 dB-down window of the C-weighted history
                // itself (or the whole of it under --full_history).  A state
                // with no A-weighted energy has no band energy either.
                for (std::size_t it = 0; it < nt; it++)
                {
                    c_levels[it] = levels[it] == -std::numeric_limits<double>::infinity()
                            ? levels[it]
                            : weighted_level(band_levels.data() + it * history_bands, c_gains);
                }
                csel_levels[index] = options.full_history
                        ? full_sound_exposure_level(c_levels, integration_step, ranges, sound_speed)
                        : sound_exposure_level(c_levels, integration_step, ranges, sound_speed);
            }
            if (want_weighted_level)
            {
                // The history is complete here (needs_full_history), so the
                // weighting runs over the whole operation, in reception time:
                // each sample lasts what it lasts in the SEL integral.
                reception_durations(nt, integration_step, ranges, sound_speed, durations);
                // One pass gives both; a threshold of +infinity costs nothing.
                const auto weighted = weighted_maximum_and_time_above(levels, durations, options.lamax_weighting,
                        options.want_ta ? options.ta_threshold : std::numeric_limits<double>::infinity());
                if (options.want_lamax)
                    lamax[index] = std::max(lamax[index], weighted.maximum);
                // NAx: operations whose own peak reaches the threshold.  A
                // trajectory stands for job.weight operations, so the count
                // adds by its weight (a fractional count gives a fractional NAx).
                if (options.want_nax && weighted.maximum >= options.nax_threshold)
                    nax[index] += job.weight;
                // Time above: seconds, adding across operations.
                if (options.want_ta)
                    totals.time_above[index] += job.weight * weighted.time_above;
            }
            if (options.want_coverage)
            {
                // The energy-weighted share of this observer's operation that
                // the database did not measure.  Every state's level is in
                // `levels` (the report turns bounding off), and each is weighted
                // as the SEL integrated it.
                double energy = 0.0;
                double unmeasured_energy = 0.0;
                for (std::size_t it = 0; it < nt; it++)
                {
                    const double weight = sel_weight[it];
                    if (weight == 0.0)
                        continue;
                    const Vec3 direction = (observers[index] - source_states[it].position).normalized();
                    const auto angles = lookup_angles(source_states[it].frame, direction,
                            !database.upper_hemisphere_measured());
                    energy += weight;
                    unmeasured_energy += weight
                            * database.unmeasured_share(angles, source_states[it].blend);
                }
                // Scaled at the add, so a count of N sums as N repeats would.
                coverage[index][0] += job.weight * unmeasured_energy;
                coverage[index][1] += job.weight * energy;
            }
            if (options.want_provenance)
            {
                // Each state's energy, to the axes it was clamped on and, by
                // its blend weights, to the spheres it was blended from.  The
                // per-sphere sums are dense per thread, and only the spheres
                // touched are reset and merged into the observer's sparse list.
                thread_local std::vector<double> by_sphere;
                thread_local std::vector<std::uint32_t> touched;
                by_sphere.resize(database.condition_count(), 0.0);
                auto& totals_here = totals.provenance[index];
                for (std::size_t it = 0; it < nt; it++)
                {
                    // The state's SEL energy, times the operations it stands
                    // for, so the totals and the spheres' shares stay one sum.
                    const double weight = sel_weight[it] * job.weight;
                    if (weight == 0.0)
                        continue;
                    totals_here[0] += weight;
                    for (int axis = 0; axis < condition_flag_count; axis++)
                        if (state_provenance[it] & (1u << axis))
                            totals_here[1 + axis] += weight;
                    const auto& blend = source_states[it].blend;
                    for (std::size_t k = 0; k < blend.weights.size(); k++)
                    {
                        if (!(blend.weights[k] > 0.0))
                            continue;
                        const auto sphere = static_cast<std::uint32_t>(blend.indices[k]);
                        if (by_sphere[sphere] == 0.0)
                            touched.push_back(sphere);
                        by_sphere[sphere] += weight * blend.weights[k];
                    }
                }
                std::sort(touched.begin(), touched.end());
                auto& list = totals.sphere_energy[index];
                std::vector<std::pair<std::uint32_t, double>> merged;
                merged.reserve(list.size() + touched.size());
                std::size_t a = 0;
                for (const std::uint32_t sphere : touched)
                {
                    while (a < list.size() && list[a].first < sphere)
                        merged.push_back(list[a++]);
                    double energy = by_sphere[sphere];
                    if (a < list.size() && list[a].first == sphere)
                        energy += list[a++].second;
                    merged.emplace_back(sphere, energy);
                    by_sphere[sphere] = 0.0;
                }
                merged.insert(merged.end(), list.begin() + static_cast<std::ptrdiff_t>(a), list.end());
                list.swap(merged);
                touched.clear();
            }
            for (std::size_t m = 0; m < day_metric_count; m++)
            {
                if (!options.want_day_metric[m] || day_pass[m] == DayPass::none)
                    continue;
                if (day_pass[m] == DayPass::the_sel)
                {
                    totals.day_energy[m][index] += day_weight[m] * std::pow(10.0, 0.1 * integrated_levels[index]);
                    continue;
                }
                // Energies add across operations; the logarithm comes at the
                // end, once every trajectory has contributed.  Each observer
                // index belongs to exactly one worker, and the trajectories run
                // one after another, so this needs no synchronization.  The
                // operation's window is its own penalized 10 dB-down window, or
                // the whole history under --full_history.  The weight
                // multiplies the energy; folded into the exponent it would
                // round differently.
                totals.day_energy[m][index] += day_weight[m] * std::pow(10.0,
                        0.1 * penalized_exposure_level(levels, integration_step, day_penalty[m], penalized,
                                      ranges, sound_speed, options.full_history));
            }
            if (options.write_time_history)
            {
                // Write out time history files
                auto time_history_filename = options.time_history_prefix + history_infix + "_i"
                        + std::to_string(index) + ".csv";
                write_time_history_file(levels, integration_step, observer(0), observer(1),
                        time_history_filename, options.unit);
            }
            if (want_frames)
            {
                // Every state's bands into 0.5 s reception-time frames (see
                // band_frames.hpp), weighted by its share of reception time as
                // the SEL integral weights it.  A culled state (level -inf)
                // carries no energy here either.
                frames.accumulate(
                        nt, [&](std::size_t it) { return std::isfinite(levels[it]); },
                        [&](std::size_t it) { return trajectory.states[it].t + ranges[it] / sound_speed; },
                        [&](std::size_t it) {
                            return integration_step
                                    * reception_factor(it, nt, integration_step, sound_speed,
                                            [&](std::size_t k) { return ranges[k]; });
                        },
                        [&](std::size_t it) { return band_levels.data() + it * history_bands + first_band; });
                if (want_epnl)
                {
                    auto& row = epnl_rows[index];
                    pnl_history.resize(frames.count());
                    for (std::size_t k = 0; k < frames.count(); k++)
                        for (std::size_t b = 0; b < pnl_band_count; b++)
                            pnl_history[k][b] = frames.level(k, pnl_bands[b] - first_band);
                    // Every frame below the noy floor (or none at all): no
                    // finite PNLT, and no EPNL.  Checked here rather than caught,
                    // since a far observer on a grid is an ordinary case.
                    row.valid = false;
                    for (const auto& spectrum : pnl_history)
                        row.valid = row.valid || std::isfinite(perceived_noise_level(spectrum));
                    if (row.valid)
                    {
                        auto& r = row.result;
                        r = effective_perceived_noise_level(pnl_history, pnl_record_seconds, true,
                                options.epnl_masking);
                        row.tone = r.pnltm - r.pnl[r.peak];
                        row.tone_band = r.tone[r.peak] > 0.0 ? pnl_band_centers[r.tone_band[r.peak]]
                                                             : std::numeric_limits<double>::quiet_NaN();
                        row.pnltm_time = frames.center(r.peak);
                        row.first_time = frames.center(r.first);
                        row.last_time = frames.center(r.last);
                        r.pnlt = {};
                        r.pnl = {};
                        r.tone = {};
                        r.tone_band = {};
                    }
                }
                if (want_band_frames)
                {
                    auto& kept = observer_frames[index];
                    kept.first = frames.first();
                    kept.count = frames.count();
                    kept.levels.resize(frames.count() * frames.bands());
                    for (std::size_t k = 0; k < frames.count(); k++)
                        for (std::size_t b = 0; b < frames.bands(); b++)
                            kept.levels[k * frames.bands() + b] = static_cast<float>(frames.level(k, b));
                }
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
        terrain_samples += terrain_here;
        terrain_raised += raised_here;
    });

    if (terrain_raster)
    {
        if (ground_model == nullptr)
            log("--terrain: no ground model runs, so the terrain sets only where the observers stand.");
        else if (terrain_raised.load() == 0)
            log("--terrain: ", terrain_samples.load(), " evaluations took their path's mean ground plane.");
        else
        {
            // A source behind rising ground is heard as if the ground were not
            // there: the level is too high, and the run used to say so only in
            // the log.  A receiver in a dip raises the same flag.
            const auto observers_raised = std::count(raised_at.begin(), raised_at.end(), char(1));
            warn("--terrain: at ", observers_raised, " of ", observers.size(), " observers, paths had a "
                    "source or receiver below their mean ground plane (", terrain_raised.load(), " of ",
                    terrain_samples.load(), " evaluations), taken onto it. A source behind rising ground "
                    "is not screened, so those observers' levels are too high there; a receiver in a dip "
                    "is only approximate.");
        }
    }

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

    if (options.want_csel)
    {
        log("Writing the C-weighted sound exposure level to ", job.csel_output, ".");
        write_levels(x_obs, y_obs, csel_levels, job.csel_output, options.scattered, options.unit);
        if (options.latlon)
        {
            write_geodetic_levels(x_obs, y_obs, csel_levels, *options.geodetic.frame,
                    latlon_companion(job.csel_output), options.scattered);
        }
    }

    if (!options.metric_file.empty())
    {
        append_metrics_file(options.metric_file, job.name, x_obs, y_obs, options.scattered, options.unit,
                integrated_levels);
    }

    if (want_epnl)
        write_epnl_file(epnl_rows, microphones, options.unit, job.epnl);
    if (want_band_frames)
    {
        BandFramesFile meta;
        meta.trajectory = job.trajectory;
        meta.database = job.database;
        meta.frame_length = pnl_record_seconds;
        meta.band_centers.assign(database.band_centers().begin() + static_cast<std::ptrdiff_t>(first_band),
                database.band_centers().begin() + static_cast<std::ptrdiff_t>(last_band));
        meta.length_unit = options.unit.name;
        for (std::size_t i = 0; i < observer_count; i++)
        {
            meta.x.push_back(microphones[i].x / options.unit.to_meters);
            meta.y.push_back(microphones[i].y / options.unit.to_meters);
            meta.z.push_back(microphones[i].z / options.unit.to_meters);
            meta.height.push_back(microphone_height[i] / options.unit.to_meters);
        }
        if (!options.scattered)
        {
            for (const double x : x_obs)
                meta.grid_x.push_back(x / options.unit.to_meters);
            for (const double y : y_obs)
                meta.grid_y.push_back(y / options.unit.to_meters);
        }
        write_band_frames(job.band_frames, meta, observer_frames);
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

// What an operations file adds to a run beyond the per-trajectory lists it
// fills: each row's name, how its failures are reported, its weight, its
// period counts, whether its clock replaces the track's, and its lateral
// offset (in --units; converted at the seam).
struct OperationRun
{
    OperationsFile file;
    std::vector<std::string> names;
    std::vector<std::string> labels;
    std::vector<double> weights;
    std::vector<std::optional<std::array<double, 3>>> periods;
    std::vector<char> clock_replaces_track;
    std::vector<double> offsets;
};

// The per-trajectory lists -t and its companions fill on the command line.
struct PerPathLists
{
    std::vector<std::string>& path_files;
    std::vector<std::string>& database_files;
    std::vector<std::string>& vehicle_files;
    std::vector<std::string>& out_files;
    std::vector<std::string>& atmosphere_files;
    std::vector<std::string>& wind_texts;
    std::vector<double>& humidity_overrides;
    std::vector<double>& temperature_overrides;
    std::vector<double>& pressure_overrides;
    std::vector<double>& start_utcs;
    std::vector<std::string>& state_files;
    std::vector<std::string>& csel_files;
    std::vector<std::string>& epnl_files;
    std::vector<std::string>& band_frames_files;
};

// Reads --operations and fills the per-trajectory lists from it, one entry per
// row, so everything after (the naming, the distinctness checks, the manifest,
// the list-level air checks, the job loop) runs as it does for several -t.
// A per-trajectory option given once fills the rows' empty cells; given more
// than once it is refused, since the file is where per-operation values go.
OperationRun expand_operations(const std::string& operations_file, PerPathLists lists, bool trajectories_given,
        bool database_given, bool vehicle_given, bool output_given, double utc_offset_hours)
{
    if (trajectories_given)
        input_error("--operations names the trajectories; give them in its track column, not with -t.");
    OperationRun run;
    run.file = read_operations(operations_file);
    const auto& rows = run.file.rows;
    const auto where = [&](const Operation& row) {
        return "operation " + row.name + " (" + operations_file + " line " + std::to_string(row.line) + ")";
    };
    const auto at_most_one = [](std::size_t given, const char* option, const char* column) {
        if (given > 1)
            input_error("With --operations, give per-operation ", option, " in the file's ", column, " column.");
    };
    for (const auto& [list, option] : {std::pair{&lists.state_files, "-S"}, {&lists.csel_files, "--csel_file"},
                 {&lists.epnl_files, "--epnl_file"}, {&lists.band_frames_files, "--band_frames"}})
        if (list->size() > 1)
            input_error("With --operations, give ", option, " once: each operation's file is named after its row.");
    if (output_given)
        at_most_one(lists.out_files.size(), "-o", "output");
    if (database_given)
        at_most_one(lists.database_files.size(), "-d", "database");
    if (vehicle_given)
        at_most_one(lists.vehicle_files.size(), "-v", "vehicle");

    // The files each row names, or the one -d/-v for the rows that name none.
    // The demo's S-76D defaults never fill an operations file: a row without a
    // database would otherwise run the demo aircraft without a word.
    const auto file_for = [&](const Operation& row, const std::string& cell, bool given,
                                  const std::vector<std::string>& list, const char* column, const char* option) {
        if (!cell.empty())
            return cell;
        if (!given)
            input_error(where(row), " names no ", column, ": fill the file's ", column, " column, or give ", option,
                    " once for the rows that leave it empty.");
        return list.front();
    };
    std::vector<std::string> databases, vehicles, outputs;
    const auto output_folder = std::filesystem::path(lists.out_files.front()).parent_path();
    for (const auto& row : rows)
    {
        databases.push_back(file_for(row, row.database, database_given, lists.database_files, "database", "-d"));
        vehicles.push_back(file_for(row, row.vehicle, vehicle_given, lists.vehicle_files, "vehicle", "-v"));
        // A row's output cell beats -o; one row keeps -o as given, and
        // several are named after their rows, beside it.
        if (!row.output.empty())
            outputs.push_back(row.output);
        else if (rows.size() == 1 && output_given)
            outputs.push_back(lists.out_files.front());
        else
            outputs.push_back((output_folder / (row.name + "_LAE.csv")).string());
    }

    // The air: a column present must be filled in every row, unless its option
    // was given once, which then fills the empty cells.  The list-level checks
    // (wind with --atmosphere, --atmosphere with --humidity, ...) then see the
    // filled lists, as they would see the command line's.
    const auto fill = [&](auto& list, const char* option, const char* column, bool has_column, auto value_of) {
        at_most_one(list.size(), option, column);
        if (!has_column)
            return;
        std::remove_reference_t<decltype(list)> filled;
        for (const auto& row : rows)
        {
            if (const auto value = value_of(row))
                filled.push_back(*value);
            else if (list.size() == 1)
                filled.push_back(list.front());
            else
                input_error(where(row), ": the ", column, " column is empty. Fill it in every row, or give ", option,
                        " once for the empty cells.");
        }
        list = std::move(filled);
    };
    const auto text = [](const std::string& value) {
        return value.empty() ? std::optional<std::string>() : std::optional<std::string>(value);
    };
    fill(lists.atmosphere_files, "--atmosphere", "atmosphere", run.file.has("atmosphere"),
            [&](const Operation& row) { return text(row.atmosphere); });
    fill(lists.wind_texts, "--wind", "wind_speed and wind_direction", run.file.has("windspeed"),
            [&](const Operation& row) { return text(row.wind_text); });
    fill(lists.humidity_overrides, "--humidity", "humidity", run.file.has("humidity"),
            [](const Operation& row) { return row.humidity; });
    fill(lists.temperature_overrides, "--temperature", "temperature", run.file.has("temperature"),
            [](const Operation& row) { return row.temperature; });
    fill(lists.pressure_overrides, "--pressure", "pressure", run.file.has("pressure"),
            [](const Operation& row) { return row.pressure; });

    // The clock: a row's local_time replaces whatever clock its track carries,
    // dating t = 0; a row without one keeps the track's own clock, or the one
    // --start_utc.
    at_most_one(lists.start_utcs.size(), "--start_utc", "local_time");
    std::vector<double> start_utcs;
    for (const auto& row : rows)
    {
        run.clock_replaces_track.push_back(row.local_time.has_value());
        if (row.local_time)
        {
            double utc = std::fmod(*row.local_time - utc_offset_hours * 3600.0, 86400.0);
            if (utc < 0.0)
                utc += 86400.0;
            start_utcs.push_back(utc);
        }
        else
            start_utcs.push_back(lists.start_utcs.empty() ? std::numeric_limits<double>::quiet_NaN()
                                                          : lists.start_utcs.front());
    }
    if (run.file.has("localtime") || !lists.start_utcs.empty())
        lists.start_utcs = std::move(start_utcs);

    lists.path_files.clear();
    for (const auto& row : rows)
    {
        lists.path_files.push_back(row.track);
        run.names.push_back(row.name);
        run.labels.push_back(where(row));
        run.weights.push_back(operation_count(row));
        run.periods.push_back(row.periods);
        run.offsets.push_back(row.offset.value_or(0.0));
    }
    lists.database_files = std::move(databases);
    lists.vehicle_files = std::move(vehicles);
    lists.out_files = std::move(outputs);
    // The SEL files are single events whatever a row's count: the count
    // weights only the run totals.
    if (rows.size() == 1 && run.weights.front() != 1.0)
        info("The SEL file ", lists.out_files.front(), " is one operation's; its count of ", run.weights.front(),
                " weights only the run totals.");
    return run;
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
    std::string lden_file = "";
    std::string cnel_file = "";
    std::string lamax_file = "";
    std::string lamax_weighting = "fast";
    std::string nax_file = "";
    double nax_threshold = std::numeric_limits<double>::quiet_NaN();
    std::string ta_file = "";
    double ta_threshold = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::string> csel_files;
    std::string coverage_file = "";
    std::string provenance_file;
    std::string provenance_spheres_file;
    bool full_history = false;
    bool epnl_masking = true;
    double utc_offset = 0.0;
    std::string lday_file, levening_file, lnight_file, laeq_file;
    std::vector<double> laeq_window_hours;
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
    std::string contour_file;
    std::string contour_metric = "sel";
    std::string raster_file;
    double dt = 0.5;
    // NOTE: D default-initializes `double` to NaN, not 0, so the !isNaN guard
    // below is FALSE by default and D never mixes in an ambient level unless
    // -a is given.  Initializing this to 0.0 added a spurious 0 dB ambient
    // floor to every run, which pinned quiet observers near 0 dB.
    double ambient_level = std::numeric_limits<double>::quiet_NaN();
    bool write_time_history = false;
    std::string time_history_prefix = "";
    std::string band_history_prefix;
    std::vector<std::string> epnl_files;
    std::vector<std::string> band_frames_files;
    std::string metric_file = "";
    bool qsam = true;
    bool self_test = false;
    unsigned thread_count = 0;
    double db_cutoff = std::numeric_limits<double>::quiet_NaN();
    bool spectral = true;
    std::vector<double> band_limits;
    bool ground_effect = true;
    double max_load_factor = 2.0;
    std::string blend_policy = "barycentric";
    std::string load_factor_correction = "one_sided";
    double min_load_factor = 0.5;
    bool accept_vehicle_mismatch = false;
    DatabaseLoadOptions database_load;
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
    std::string terrain_file;
    double terrain_offset = 0.0;  // meters
    double terrain_step = 0.0;    // meters; 0 follows the terrain's own spacing
    bool observers_on_terrain = false;
    std::string plate_mic_azimuth = "outboard";
    bool exact_ground = false;
    // The air: once for the run, or once per trajectory (empty: not given).
    std::vector<double> humidity_overrides;
    std::vector<double> temperature_overrides;
    std::vector<double> pressure_overrides;
    std::vector<std::string> atmosphere_files;
    // --wind and the options that shape the profile built from it.  The
    // lengths are meters until the --units seam, like every default.
    std::vector<std::string> wind_texts;
    double roughness_length = SurfaceWind().roughness_length;
    double surface_layer_depth = SurfaceWind().surface_layer_depth;
    std::string stability = "neutral";
    int profile_levels = SurfaceWind().levels_per_decade;
    std::string write_atmosphere;
    std::string refraction = "none";
    double shadow_floor = 20.0;
    std::string shadow_model = "legacy";
    double shadow_phase = default_resolved_phase;
    std::string cramped_ground = "series";
    bool tunneling = false;
    bool skip_fill = false;
    double scattering = 0.0;
    // Empty follows the database's speed_reference attribute.
    std::string speed_reference;
    // Empty follows the database's azimuth_reference attribute.
    std::string nose;
    // WGS84 input and output.  An unset origin is taken from the first geodetic
    // trajectory; an unset bearing leaves the axes east/north, which is what a
    // map wants -- the shipped AC tracks are the other case, their +x running
    // along the flight track.
    std::vector<double> origin;
    double frame_bearing = std::numeric_limits<double>::quiet_NaN();
    double ground_elevation = 0.0;
    std::vector<double> start_utcs;
    // Waypoints state intent, not a constraint to be met to the inch: forcing
    // a path exactly through them manufactures curvature, and the load factor
    // is that curvature, so it lands in CT and then in the levels.  Ten feet is
    // loose enough to absorb both coordinate rounding and the corner a planner
    // did not mean literally.
    double path_smoothing = 10.0 * ft2m; // 10 ft, in meters
    bool latlon = false;

    CLI::App app{"Noise Informed Community Environment Operation Planning System (NICE-OPS)."};
    auto* database_option = app.add_option("-d,--database", database_files,
            "netCDF noise database file. Once, or once per trajectory.");
    app.add_option("-l,--logfile", log_file, "Output log file.");
    app.add_option("-L,--loglevel", log_level, "Log verbosity: trace, info, warn, error, fatal.")
            ->check(CLI::IsMember({"trace", "info", "warn", "error", "fatal"}));
    auto* trajectory_option = app.add_option("-t,--trajectoryfile", path_files,
            "CSV trajectory input file. Repeat for several flight paths.");
    std::string operations_file;
    app.add_option("--operations", operations_file,
            "CSV of operations: per row a track, its count, and optionally its database, vehicle, clock, "
            "lateral offset and air. Replaces -t.");
    auto* vehicle_option = app.add_option("-v,--vehicle", vehicle_files,
            "JSON vehicle definition file. Once, or once per trajectory.");
    auto* output_option = app.add_option("-o,--output", out_files,
            "Output file name, once per trajectory. Default integrated/197421_LAE.csv; with several "
            "trajectories and one name (or the default), each is named after its track.");
    app.add_flag("-p,--plotting", plotting, "Enable plotting.");
    app.add_option("-f,--figfile", figure_file, "Filename to save figure.");
    app.add_option("-F,--statefigfile", state_figure_file, "Filename to save state figures.");
    app.add_option("--contour_file", contour_file,
            "Write contour polygons of --contour_metric at the --cmin/--cmax/--cstep levels, with "
            "their areas, as GeoJSON (WGS84). Needs the observer grid and a frame (--origin, or a "
            "geodetic track).");
    app.add_option("--raster_file", raster_file,
            "Write every gridded level this run writes (SEL for one trajectory, the day and period "
            "levels, LAmax, NAx, TA, the unmeasured share) to one CF netCDF, with WGS84 latitude "
            "and longitude when the run has a frame. Needs the observer grid.");
    std::vector<std::string> contour_metrics = {"sel"};
    for (std::size_t m = 0; m < day_metric_count; m++)
        contour_metrics.push_back(day_metric_key(static_cast<DayMetric>(m)));
    contour_metrics.push_back("lamax");
    app.add_option("--contour_metric", contour_metric,
            "What --contour_file contours: sel (default; one trajectory), dnl, lden, cnel, lday, "
            "levening, lnight, laeq or lamax, whose own output must be asked for too.")
            ->check(CLI::IsMember(contour_metrics));
    app.add_option("--cmin", contour_min, "Minimum contour level, dB.");
    app.add_option("--cmax", contour_max, "Maximum contour level, dB.");
    app.add_option("--cstep", contour_step, "Contour level step, dB.");
    app.add_option("-s,--step", dt, "Resampling step, in seconds, for WGS84 tracks that are "
            "uneven or coarser than this; never coarsens a track, and an AC track is used as "
            "sampled. The sound exposure level integral uses the trajectory's own sample "
            "interval.");
    app.add_option("-a,--ambient", ambient_level,
            "Ambient noise level, dB, added to every sample's level in energy.");
    app.add_option("-i,--history", time_history_prefix, "Time history file prefix.");
    app.add_option("--band_history", band_history_prefix,
            "Per-band time history file prefix: one file per observer, a row per trajectory "
            "state with its emission and reception times, geometry, A-weighted level and each "
            "band's unweighted level at the microphone. Needs the spectral path.");
    app.add_option("--epnl_file", epnl_files,
            "Write each observer's EPNL (14 CFR 36 A36.4), PNLTM and duration correction to this "
            "CSV, from one-third-octave levels in 0.5 s reception-time frames. Once, or once per "
            "trajectory; with several trajectories and one name, each is named after its track "
            "(<stem>_EPNL.csv). Needs the spectral path and the 24 bands 50 Hz-10 kHz; refused "
            "with -a.");
    app.add_flag("--epnl_masking,!--no-epnl_masking", epnl_masking,
            "--epnl_file's tone correction: a band below its noy threshold (Table A36-3's SPL(d)) "
            "enters the tone correction at that threshold and carries no tone itself, as tones in "
            "masked bands are excluded from a measured spectrum's (default on; --no-epnl_masking "
            "takes the procedure as written, which reads absorption's plunge in the top bands as "
            "tones).");
    app.add_option("--band_frames", band_frames_files,
            "Write each observer's band levels in 0.5 s reception-time frames to this netCDF file "
            "(level by observer, frame and band, with the coordinates). Once, or once per "
            "trajectory; with several trajectories and one name, each is named after its track "
            "(<stem>_band_frames.nc). Needs the spectral path; refused with -a.");
    app.add_option("-m,--metricfile", metric_file,
            "File to append overall noise metrics: the case name, the energy mean and maximum "
            "level, and the fraction of observers at or above 65 dB.");
    app.add_option("-b,--observer", observer_file, "Observer file name.");
    app.add_flag("-q,--qsam,!--no-qsam", qsam, "Q-SAM flight path angle: the velocity's angle to the tip-path plane, not to the "
            "horizontal (default on; --no-qsam disables)");
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
            "Write the day-night level (DNL: the energy average over a 24-hour day, 22:00-07:00 "
            "night +10 dB), summed over every trajectory, to this file.");
    app.add_option("--lden_file", lden_file,
            "Write the day-evening-night level (Lden, EU Directive 2002/49/EC: 19:00-23:00 "
            "evening +5 dB, 23:00-07:00 night +10 dB, over 24 hours), summed over every "
            "trajectory, to this file.");
    app.add_option("--cnel_file", cnel_file,
            "Write the community noise equivalent level (CNEL, California Title 21: 19:00-22:00 "
            "evening +4.77 dB, 22:00-07:00 night +10 dB, over 24 hours), summed over every "
            "trajectory, to this file.");
    app.add_option("--lday_file", lday_file,
            "Write the day level (Lday, EU Directive 2002/49/EC: the energy heard 07:00-19:00 over "
            "those 12 hours), summed over every trajectory, to this file.");
    app.add_option("--levening_file", levening_file,
            "Write the evening level (Levening: 19:00-23:00 over 4 hours), summed over every "
            "trajectory, to this file.");
    app.add_option("--lnight_file", lnight_file,
            "Write the night level (Lnight: 23:00-07:00 over 8 hours), summed over every "
            "trajectory, to this file.");
    app.add_option("--laeq_file", laeq_file,
            "Write the equivalent level LAeq,T over the local clock window --laeq_window, summed "
            "over every trajectory, to this file.");
    app.add_option("--laeq_window", laeq_window_hours,
            "Local clock window of --laeq_file, hours: START END, END not included; START > END "
            "wraps through midnight (\"7 23\" is LAeq,16h).")
            ->expected(2);
    app.add_option("--lamax_file", lamax_file,
            "Write each observer's LAmax, the highest A-weighted level over every operation, "
            "time-weighted per --lamax_weighting (LAFmax by default), dB, to this file.");
    auto* lamax_weighting_option = app.add_option("--lamax_weighting", lamax_weighting,
            "Time weighting for --lamax_file, --nax_file and --ta_file, in reception time: fast "
            "(IEC 61672 F, 0.125 s; default), slow (S, 1 s), leq0.5 (the largest 0.5 s Leq, ECAC "
            "Doc 32's LASmax surrogate) or instantaneous (the largest sample, what --lamax_file "
            "wrote before). Refused when the history is coarser than half the time constant.")
            ->check(CLI::IsMember({"fast", "slow", "leq0.5", "instantaneous"}));
    app.add_option("--nax_file", nax_file,
            "Write each observer's NAx, the number of operations whose LAmax (per "
            "--lamax_weighting) reaches --nax_threshold, to this file.");
    app.add_option("--nax_threshold", nax_threshold,
            "The LAmax, dBA, that an operation must reach to count toward NAx. Required with "
            "--nax_file: there is no default, so the study's threshold is always stated.");
    app.add_option("--ta_file", ta_file,
            "Write each observer's time above --ta_threshold: seconds of reception time the "
            "level, time-weighted per --lamax_weighting, exceeds it, summed over every "
            "trajectory.");
    app.add_option("--ta_threshold", ta_threshold,
            "The level, dBA, that --ta_file counts the time above. Required with --ta_file: "
            "there is no default.");
    app.add_option("--csel_file", csel_files,
            "Write the C-weighted sound exposure level (LCE) per observer, by SEL's window rule, "
            "to this file: one name, or one per trajectory (with one name and several "
            "trajectories, each writes <track>_LCE.csv beside it). Spectral path only; refused "
            "with -a or --db_cutoff.");
    app.add_option("--provenance_file", provenance_file,
            "Write, per observer, the shares of its SEL energy from states the database clamped "
            "on each axis or held at a load-factor limit, and the sphere that supplied the most.");
    app.add_option("--provenance_spheres", provenance_spheres_file,
            "With --provenance_file, also write every sphere's share of every observer's SEL "
            "energy, and each sphere's stored condition, to this netCDF file.");
    app.add_option("--coverage_file", coverage_file,
            "Write, per observer, the energy-weighted share of its levels that came from cells "
            "the database did not measure (straight rays only; needs a database with coverage).");
    app.add_flag("--full_history", full_history,
            "Integrate each exposure (SEL, C-SEL, and each operation's share of DNL, Lden, "
            "CNEL and the period levels) over the whole history rather than the event's 10 dB-down window, as AAM and "
            "NORAH2 do. Slower: no window is found.");
    app.add_option("--utc_offset", utc_offset,
            "Hours from UTC to local time, for the DNL, Lden and CNEL periods and the period levels.");
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
            "spectra, except that asking for --spectral, or for an option that needs bands, "
            "with such a database is an error). --no-spectral propagates one broadband level "
            "with the database's EAA "
            "factor instead.");
    app.add_option("--band_limits", band_limits,
            "Frequency range integrated into the metric, Hz: low high. Default 20 20000.")
            ->expected(2);
    app.add_option("--max_load_factor", max_load_factor,
            "Largest load factor, in g, a fixed_load_factor database's load-factor correction "
            "is taken to. Above it the correction is held, and the run says so (default 2).")
            ->check(CLI::PositiveNumber);
    app.add_option("--blend", blend_policy,
            "How a state's flight condition picks the database spheres it is evaluated from: "
            "barycentric (default) -- blended in dB over the Delaunay simplex holding it, "
            "nearest outside the hull; nearest -- the one stored condition nearest it, at every "
            "state (the closest a trajectory tool comes to ECAC Doc 32's one hemisphere per "
            "segment); or aam -- AAM v3.1's search, up to four spheres bracketing it in speed "
            "and then flight path angle, within 20 kt and 3 deg descending or 5 deg otherwise, "
            "the nearest sphere where nothing brackets it.")
            ->check(CLI::IsMember({"barycentric", "nearest", "aam"}));
    app.add_option("--load_factor_correction", load_factor_correction,
            "A fixed_load_factor database's load-factor law, 20 lg(CT/CT_ref): one_sided "
            "(default) applies it only above the reference load factor, as panam does; "
            "two_sided below it too, down to --min_load_factor.")
            ->check(CLI::IsMember({"one_sided", "two_sided"}));
    auto* min_load_factor_option = app.add_option("--min_load_factor", min_load_factor,
            "--load_factor_correction two_sided: the smallest load factor, in g, the law is "
            "taken to; below it the correction is held, and the run says so (default 0.5).")
            ->check(CLI::Validator(
                    [](std::string& text) -> std::string {
                        double value = 0.0;
                        try
                        {
                            value = std::stod(text);
                        }
                        catch (const std::exception&)
                        {
                            return "not a number";
                        }
                        return value > 0.0 && value <= 1.0 ? std::string() : std::string("must be in (0, 1]");
                    },
                    "in (0, 1]"));
    app.add_option("--screen_database", database_load.screen_file,
            "Write every database entry the load-time screen finds suspect (a dBA cell far from "
            "its neighbors, a band far above its cell's spectrum) to this CSV.");
    auto* mask_suspect_option = app.add_flag("--mask_suspect_cells", database_load.mask_suspect_cells,
            "Drop the entries the load-time screen finds suspect, as gating drops a direction, "
            "instead of only warning about them.");
    app.add_flag("--accept_vehicle_mismatch", accept_vehicle_mismatch,
            "Warn, rather than refuse, a vehicle file that disagrees with the database by more "
            "than 2 percent.");
    app.add_option("--humidity", humidity_overrides,
            "Relative humidity in percent, overriding the vehicle file. Not with --atmosphere. "
            "Once, or once per trajectory.");
    app.add_option("--temperature", temperature_overrides,
            "Temperature in K, overriding the vehicle file. Not with --atmosphere. Once, or once "
            "per trajectory.");
    app.add_option("--pressure", pressure_overrides,
            "Ambient pressure in kPa, overriding the vehicle file. With --atmosphere, only when "
            "the file has no pressure column; it is then the surface pressure. Once, or once per "
            "trajectory.");
    app.add_option("--atmosphere", atmosphere_files,
            "A stratified atmosphere: CSV of temperature, relative humidity and optionally "
            "pressure against height (see the manual's Stratified atmosphere), with optional wind "
            "columns (wind_speed_mps or wind_speed_kt with wind_direction, or u_mps and v_mps). "
            "Per-band absorption is integrated through "
            "it, and its ground level replaces the vehicle file's atmosphere. Spectral path "
            "only. Once, or once per trajectory.");
    app.add_option("--wind", wind_texts,
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
            "the fit, resolving them only as finely as the highest band evaluated can. Both "
            "give shadow zones under upward refraction.")
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
    app.add_option("--shadow_model", shadow_model,
            "--refraction stratified, with --ground_effect: how a shadow boundary's field is set per "
            "band. legacy (default): each band's boundary interpolated between the layered rays' and "
            "the path's linear fit's. resolved: each band's boundary from the profile it resolves -- "
            "below the height where it gathers --shadow_phase of vertical phase, the profile's tangent "
            "there -- and a band whose boundary is within two diffraction lengths of the source takes "
            "--cramped_ground.")
            ->check(CLI::IsMember({"legacy", "resolved"}));
    auto* shadow_phase_option = app.add_option("--shadow_phase", shadow_phase,
            "--shadow_model resolved: the vertical phase, radians, a band must gather to resolve the "
            "profile near the ground (default 1, the radian the layered rays' flat-floor limit takes).")
            ->check(CLI::PositiveNumber);
    auto* cramped_ground_option = app.add_option("--cramped_ground", cramped_ground,
            "--shadow_model resolved: a cramped band takes the residue series at the source's own "
            "height (series, default); or its lit side the ground model on the straight path "
            "(straight), or on the ray's equivalent geometry (ray).")
            ->check(CLI::IsMember({"series", "straight", "ray"}));
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
    app.add_option("--nose", nose,
            "Where the sphere's nose (its azimuth reference) points: velocity, the velocity "
            "--speed_reference selects; ground_velocity, the ground velocity whatever "
            "--speed_reference says; or heading, the track's compass heading at every speed "
            "(needs --frame_bearing for an AC track, and a heading column). The first two turn "
            "to the heading below 8 kt. Default: what the database's azimuth_reference records; "
            "velocity if it records nothing.")
            ->check(CLI::IsMember({"velocity", "ground_velocity", "heading"}));
    app.add_option("--speed_reference", speed_reference,
            "Which velocity sets the flight condition -- mu, the path angle and (under --nose "
            "velocity) the sphere's nose: ground, or air (the ground velocity less the wind at the "
            "aircraft, from --atmosphere or --wind). Default: what the database records it was built with; ground if it "
            "records nothing.")
            ->check(CLI::IsMember({"ground", "air"}));
    auto* ground_elevation_option = app.add_option("--ground_elevation", ground_elevation,
            "Elevation of the observer plane above mean sea level, in --units. Subtracted "
            "from a geodetic track's altitude column to give height above the observers "
            "(default 0). With --terrain it places the local z = 0 against the terrain's heights; "
            "an AC track's own alt - z must match it.");
    app.add_option("--start_utc", start_utcs,
            "Seconds since UTC midnight at t = 0 of a geodetic track or static hold, for the "
            "day-night periods. Required with -D, --lden_file, --cnel_file or a period level when the "
            "track has no utcsec column; ignored when it has one. Once, or once per trajectory.");
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
            "field) applies. Spectral path only.");
    app.add_option("--terrain", terrain_file,
            "A bare-earth elevation grid (netCDF: latitude, longitude, elevation(latitude, "
            "longitude) in meters, e.g. USGS 3DEP). Grid observers stand on it, and the ground "
            "model and plate take each path's least-squares mean ground plane (CNOSSOS-EU, "
            "NORAH2 / ECAC Doc 32) in place of the plane z = 0. Needs the run's frame (--origin, "
            "and --frame_bearing for rotated x/y, or a geodetic track). No screening or "
            "diffraction by terrain; not with --refraction or --propagation_table.");
    auto* terrain_offset_option = app.add_option("--terrain_offset", terrain_offset,
            "Added to every --terrain height, in --units, to put it on the trajectory's vertical "
            "datum (default 0). At Amedee 3DEP is NAVD88 and the 2017 tracks are ellipsoidal: "
            "-23.6 m (-77.4 ft). The ground's local z is then the height plus this, less "
            "--ground_elevation.");
    auto* terrain_step_option = app.add_option("--terrain_step", terrain_step,
            "Spacing, in --units, of the terrain samples each path's mean ground plane is fitted "
            "to (default: the terrain grid's own spacing, at least 5 m).");
    app.add_flag("--observers_on_terrain", observers_on_terrain,
            "With --terrain and -b: take each observer's ground height from the terrain, "
            "ignoring the file's z. Without it the file's z stands, and a disagreement over 1 m "
            "is warned about.");

    std::string manifest_file;
    app.add_option("--manifest", manifest_file,
            "Write a JSON record of the run to this file: the build, the command line and "
            "resolved options, the SHA-256 of every input and output, each database's build "
            "attributes, and the warnings.");
    auto* config_option = app.set_config("--config", "",
            "Read options from this TOML file, as `name = value` lines (the command line "
            "overrides it).");
    // CLI11 drops a key it does not know by default: a misspelled option would
    // be a run with the default, and exit 0.
    app.allow_config_extras(CLI::config_extras_mode::error);
    app.set_version_flag("--version", build_info().version);

    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        return app.exit(e);
    }

    // --operations fills the per-trajectory lists, one entry per row, before
    // anything reads them.
    std::optional<OperationRun> operations;
    if (!operations_file.empty())
    {
        operations = expand_operations(operations_file,
                PerPathLists {path_files, database_files, vehicle_files, out_files, atmosphere_files, wind_texts,
                        humidity_overrides, temperature_overrides, pressure_overrides, start_utcs, state_files,
                        csel_files, epnl_files, band_frames_files},
                trajectory_option->count() > 0, database_option->count() > 0, vehicle_option->count() > 0,
                output_option->count() > 0, utc_offset);
        // Read, but not run yet: the day, evening and night counts and the
        // lateral offsets are refused here, and only here.
        if (operations->file.has("day") || operations->file.has("evening") || operations->file.has("night"))
            input_error(operations_file, ": day, evening and night counts are not supported yet. Give count, "
                    "with local_time for each operation's clock.");
        for (std::size_t i = 0; i < operations->offsets.size(); i++)
            if (operations->offsets[i] != 0.0)
                input_error(operations->labels[i], ": lateral offsets are not supported yet.");
    }

    if (!manifest_file.empty())
    {
        // Handed to main() only once it passes its checks: a refused one
        // leaves nothing to write.
        ManifestRequest request;
        request.path = manifest_file;
        request.manifest.set_command_line(argc, argv);
        request.manifest.set_options_given(app.config_to_str(false, false));
        request.manifest.set_log_file(log_file);
        if (same_file(manifest_file, log_file))
            input_error("--manifest would overwrite the log file ", log_file, ". Give another file name.");
        if (config_option->count() > 0)
            request.manifest.add_input("config", config_option->as<std::string>());
        request.manifest.add_input("operations", operations_file);
        for (const auto& f : path_files)
            request.manifest.add_input("trajectory", f);
        for (const auto& f : database_files)
            request.manifest.add_database(f);
        for (const auto& f : vehicle_files)
            request.manifest.add_input("vehicle", f);
        request.manifest.add_input("observers", observer_file);
        for (const auto& f : atmosphere_files)
            request.manifest.add_input("atmosphere", f);
        request.manifest.add_input("terrain", terrain_file);
        request.manifest.add_input("plate_table", plate_table);
        request.manifest.add_input("propagation_table", propagation_table);
        request.manifest.refuse_overwriting_an_input(manifest_file);
        request.manifest.start_hashing();
        manifest_request().emplace(std::move(request));
    }

    Database::BlendOptions blend_options;
    blend_options.policy = *parse_blend_policy(blend_policy);
    blend_options.load_factor_correction = *parse_load_factor_correction(load_factor_correction);
    blend_options.min_load_factor = min_load_factor;
    if (min_load_factor_option->count() > 0
            && blend_options.load_factor_correction != LoadFactorCorrection::two_sided)
    {
        input_error("--min_load_factor is the floor of the two-sided load-factor law; it needs "
                "--load_factor_correction two_sided.");
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
    auto expand_per_path = [path_count](auto& list, const char* what) {
        if (list.size() == path_count)
            return;
        if (list.size() == 1)
        {
            // Copy first: assign() clears the vector before reading its
            // argument, and list.front() points into the vector.
            const auto shared = list.front();
            list.assign(path_count, shared);
            return;
        }
        input_error("Got ", list.size(), " ", what, " for ", path_count,
                " trajectories; give one, or one per trajectory.");
    };
    // Each trajectory's name in the run, which its per-track files are named
    // after: the track's stem, or an operations file's row name.
    std::vector<std::string> run_names;
    if (operations)
        run_names = operations->names;
    else
        for (const auto& path : path_files)
            run_names.push_back(std::filesystem::path(path).stem().string());
    if (out_files.size() == 1 && path_count > 1)
    {
        // A single output name cannot serve several tracks, so name each after
        // its own, next to where the given name pointed.
        const auto directory = std::filesystem::path(out_files.front()).parent_path();
        out_files.clear();
        for (const auto& name : run_names)
            out_files.push_back((directory / (name + "_LAE.csv")).string());
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
            for (const auto& name : run_names)
                atmosphere_out_files.push_back((directory / (name + "_atmosphere.csv")).string());
        }
    }
    if (state_files.size() == 1 && path_count > 1)
    {
        // The same goes for the state file: one name for several tracks would
        // be rewritten by each in turn, leaving only the last one's.
        const auto directory = std::filesystem::path(state_files.front()).parent_path();
        state_files.clear();
        for (const auto& name : run_names)
            state_files.push_back((directory / (name + "_state.csv")).string());
    }
    // --csel_file, --epnl_file and --band_frames are per operation, as -S is:
    // one name for several tracks is split into one per track, beside it.
    auto name_per_track = [&](std::vector<std::string>& names, const char* suffix) {
        if (names.size() != 1 || path_count <= 1)
            return;
        const auto directory = std::filesystem::path(names.front()).parent_path();
        names.clear();
        for (const auto& name : run_names)
            names.push_back((directory / (name + suffix)).string());
    };
    name_per_track(csel_files, "_LCE.csv");
    name_per_track(epnl_files, "_EPNL.csv");
    name_per_track(band_frames_files, "_band_frames.nc");
    if (!epnl_files.empty())
        expand_per_path(epnl_files, "--epnl_file names");
    if (!band_frames_files.empty())
        expand_per_path(band_frames_files, "--band_frames names");
    expand_per_path(database_files, "database files");
    // The air, when given: once for every trajectory, or once each.
    for (auto [list, what] : {std::pair{&atmosphere_files, "--atmosphere files"}, {&wind_texts, "--wind values"}})
        if (!list->empty())
            expand_per_path(*list, what);
    for (auto [list, what] : {std::pair{&humidity_overrides, "--humidity values"},
                 {&temperature_overrides, "--temperature values"}, {&pressure_overrides, "--pressure values"},
                 {&start_utcs, "--start_utc values"}})
        if (!list->empty())
            expand_per_path(*list, what);
    if (!provenance_spheres_file.empty() && provenance_file.empty())
        input_error("--provenance_spheres writes beside --provenance_file: give that too.");
    if (!provenance_file.empty() && !std::isnan(ambient_level))
        input_error("--provenance_file attributes the aircraft's energy to the states and spheres "
                "it came from: it does not combine with -a, whose energy came from neither.");
    // Each trajectory's database by its canonical path: two spellings of one
    // file are one database, loaded once.
    std::vector<std::string> database_keys;
    for (const auto& f : database_files)
        database_keys.push_back(database_key(f));
    const std::string front_key = database_keys.front();
    const bool several_databases = std::any_of(database_keys.begin(), database_keys.end(),
            [&](const std::string& key) { return key != front_key; });
    // The provenance names spheres by their index in one database.
    if (!provenance_file.empty() && several_databases)
        input_error("--provenance_file names spheres of one database: run it once per database.");
    // Each database load writes the screen afresh, so a second database would
    // overwrite the first one's list.
    if (!database_load.screen_file.empty() && several_databases)
        input_error("--screen_database lists one database's entries: run it once per database.");
    expand_per_path(vehicle_files, "vehicle files");
    expand_per_path(out_files, "output files");
    if (!state_files.empty())
        expand_per_path(state_files, "state files");
    if (!csel_files.empty())
    {
        expand_per_path(csel_files, "C-SEL files");
        for (std::size_t i = 0; i < path_count; i++)
        {
            if (std::filesystem::path(csel_files[i]).lexically_normal()
                    == std::filesystem::path(out_files[i]).lexically_normal())
                input_error("--csel_file and -o name the same file (", csel_files[i],
                        "); the C-weighted SEL would overwrite the A-weighted one.");
        }
    }
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
        if (!csel_files.empty())
            require_distinct(csel_files, "C-SEL file");
        if (!epnl_files.empty())
            require_distinct(epnl_files, "--epnl_file");
        if (!band_frames_files.empty())
            require_distinct(band_frames_files, "--band_frames file");
        if (!time_history_prefix.empty() || !band_history_prefix.empty())
            require_distinct(run_names, "history file name (their tracks share a file name)");
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
    for (const double temperature : temperature_overrides)
        if (const auto problem = temperature_problem(temperature); !problem.empty())
            input_error("--temperature: ", problem);
    for (const double pressure : pressure_overrides)
        if (const auto problem = pressure_problem(pressure); !problem.empty())
            input_error("--pressure: ", problem);
    for (const double humidity : humidity_overrides)
        if (const auto problem = humidity_problem(humidity); !problem.empty())
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
    // An operations file's lengths (offset, wind_height) are in --units, as the
    // command line's are.  A header naming another unit means a file in it,
    // which nothing converts.
    if (operations)
    {
        const std::set<std::string> spellings = unit.to_meters == 1.0
                ? std::set<std::string> {"m", "meter", "meters", "metre", "metres"}
                : std::set<std::string> {"ft", "feet", "foot"};
        for (const auto& [token, column] : {std::pair {operations->file.offset_unit, "offset"},
                     std::pair {operations->file.wind_height_unit, "wind_height"}})
            if (!token.empty() && spellings.count(token) == 0)
                input_error(operations_file, ": the ", column, " column declares ", token, ", and this run reads "
                        "lengths in ", unit.name, " (--units). Nothing is converted: give the column in ", unit.name,
                        ", or change --units.");
    }
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
    const std::array<const std::string*, day_metric_count> day_files = {&dnl_file, &lden_file, &cnel_file,
            &lday_file, &levening_file, &lnight_file, &laeq_file};
    if (!contour_file.empty())
    {
        if (!frame)
            input_error("--contour_file writes WGS84 polygons and needs the run's frame: give --origin "
                        "<latitude> <longitude> (and usually --frame_bearing), or use a geodetic trajectory.");
        if (!observer_file.empty())
            input_error("--contour_file contours the observer grid; scattered observers (-b) have none.");
        if (!(contour_step > 0.0) || !(contour_max >= contour_min))
            input_error("--contour_file needs --cmin <= --cmax and a positive --cstep.");
        const std::string* need = contour_metric == "lamax" ? &lamax_file : nullptr;
        for (std::size_t m = 0; m < day_metric_count; m++)
            if (contour_metric == day_metric_key(static_cast<DayMetric>(m)))
                need = day_files[m];
        if (need && need->empty())
            input_error("--contour_metric ", contour_metric, " contours that metric's output: ask for it too.");
        if (contour_metric == "sel" && path_files.size() > 1)
            input_error("--contour_metric sel contours one trajectory's SEL; with several, contour a "
                        "metric that sums them (dnl, lden, cnel, a period level or lamax).");
    }
    if (!raster_file.empty() && !observer_file.empty())
        input_error("--raster_file writes the observer grid; scattered observers (-b) have none.");
    if (latlon && !frame)
    {
        input_error("--latlon has no frame to invert: give --origin <latitude> <longitude> "
                    "(and usually --frame_bearing) saying where this run's x/y origin sits, "
                    "or use a geodetic trajectory, which carries one.");
    }

    // Start stopwatch
    auto watch_start = std::chrono::steady_clock::now();
    // Each distinct database is loaded once, on first use, and held until the
    // last trajectory that names it; the front one is loaded here, for the
    // standalone jobs below as well, and held to the end, where the provenance
    // reads it.  Each Database owns its own triangulation, so several can be
    // live at once.  Keyed by database_key().
    std::map<std::string, Database> databases;
    std::map<std::string, std::size_t> last_use;
    for (std::size_t i = 0; i < path_count; i++)
        last_use[database_keys[i]] = i;
    const std::string& front_database = database_files.front();
    // The screen writes its list while the database loads, with a truncating
    // open: `-d X --screen_database X` would replace the database with it.
    // Checked against every database, since any of them could be the one read
    // when it is written.
    if (!database_load.screen_file.empty())
    {
        for (const auto& f : database_files)
            if (same_file(database_load.screen_file, f))
                input_error("--screen_database would overwrite the database it reads, ", f,
                        ". Give another file name.");
    }
    if (self_test && mask_suspect_option->count() > 0)
        input_error("--self_test checks the database as the file holds it; it does not combine "
                "with --mask_suspect_cells.");
    Database& database = databases.try_emplace(front_key, front_database, database_load).first->second;
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
        // --blend and --load_factor_correction apply here as in a run: an
        // exported sphere is what a state at this condition would hear.
        auto export_blend = blend_options;
        if (export_blend.policy == BlendPolicy::aam && std::isnan(database.tip_speed()))
            export_blend.tip_speed = configure_vehicle(vehicle_files.front()).tip_speed;
        if (export_blend.load_factor_correction == LoadFactorCorrection::two_sided
                && !database.fixed_load_factor())
        {
            input_error("--load_factor_correction two_sided applies a fixed_load_factor database's "
                    "law, and ", front_database, " is not one: its load factors are spheres.");
        }
        const auto sphere = database.interpolate_sphere(
                sphere_condition[0], sphere_condition[1], sphere_condition[2], export_blend);
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
        database.self_test(front_database);
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
    to_meters(terrain_offset_option, terrain_offset);
    to_meters(terrain_step_option, terrain_step);

    ObserverSet observer_set = scattered
            ? scattered_observers(observer_file, unit)
            : gridded_observers(grid_bounds, grid_step);
    auto& observers = observer_set.observers;
    auto& x_obs = observer_set.x;
    auto& y_obs = observer_set.y;
    std::vector<double> integrated_levels(observers.size());
    // DNL, Lden and CNEL are properties of the whole set of operations, not of
    // one flight: each trajectory's penalized exposure is an energy, and they sum.
    if (laeq_file.empty() != laeq_window_hours.empty())
        input_error(laeq_file.empty() ? "--laeq_window is for --laeq_file, which is not given."
                                      : "--laeq_file needs --laeq_window START END, local clock hours: there "
                                        "is no default period.");
    if (!laeq_window_hours.empty())
    {
        const double start = laeq_window_hours[0], end = laeq_window_hours[1];
        if (!(start >= 0.0 && start < 24.0) || !(end > 0.0 && end <= 24.0))
            input_error("--laeq_window takes local clock hours, a start in [0, 24) and an end in (0, 24]; ",
                    start, " ", end, " is not one.");
        if (start == end)
            input_error("--laeq_window ", start, " ", end, " is empty; for a whole day give 0 24.");
    }
    std::array<bool, day_metric_count> want_day_metric {};
    RunTotals totals;
    for (std::size_t m = 0; m < day_metric_count; m++)
    {
        want_day_metric[m] = !day_files[m]->empty();
        totals.day_energy[m].assign(want_day_metric[m] ? observers.size() : 0, 0.0);
    }
    const bool want_lamax = !lamax_file.empty();
    totals.lamax.assign(want_lamax ? observers.size() : 0, -std::numeric_limits<double>::infinity());
    const bool want_nax = !nax_file.empty();
    totals.nax.assign(want_nax ? observers.size() : 0, 0.0);
    const bool want_ta = !ta_file.empty();
    totals.time_above.assign(want_ta ? observers.size() : 0, 0.0);
    const bool want_coverage = !coverage_file.empty();
    totals.coverage.assign(want_coverage ? observers.size() : 0, std::array<double, 2> {0.0, 0.0});
    const bool want_provenance = !provenance_file.empty();
    totals.provenance.assign(want_provenance ? observers.size() : 0,
            std::array<double, 1 + condition_flag_count> {});
    totals.sphere_energy.assign(want_provenance ? observers.size() : 0, {});
    const bool want_csel = !csel_files.empty();
    if (want_nax && !std::isfinite(nax_threshold))
        input_error("--nax_file needs --nax_threshold, the LAmax in dBA that counts an operation. "
                "There is no default: give the threshold your study uses.");
    if (want_ta && !std::isfinite(ta_threshold))
        input_error("--ta_file needs --ta_threshold, the level in dBA it counts the time above. "
                "There is no default: give the threshold your study uses.");
    const LevelWeighting level_weighting = parse_level_weighting(lamax_weighting);
    if (lamax_weighting_option->count() > 0 && !want_lamax && !want_nax && !want_ta)
        input_error("--lamax_weighting sets the time weighting of --lamax_file, --nax_file and "
                "--ta_file, and none of them is given.");

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
                    "of frequency, and the broadband path (--no-spectral) has no "
                    "spectrum to apply it to.");
    }
    if (want_csel && !spectral)
    {
        input_error("--csel_file needs the spectral path: C-weighting reweights the band levels, "
                    "and the broadband path (--no-spectral) carries one A-weighted level.");
    }
    if (want_csel && !std::isnan(db_cutoff))
    {
        input_error("--csel_file does not combine with --db_cutoff: the cutoff is proven on the "
                    "A-weighted level, and a state it skips can still carry C-weighted energy, "
                    "which sits mostly below 100 Hz.");
    }
    if (want_csel && !std::isnan(ambient_level))
    {
        input_error("--csel_file does not combine with -a: the ambient is an A-weighted level, "
                    "with no spectrum to C-weight.");
    }
    if (!band_history_prefix.empty() && !spectral)
    {
        input_error("--band_history needs the spectral path: the broadband path (--no-spectral) "
                    "carries no bands.");
    }
    if (!propagation_table.empty() && !spectral)
        input_error("--propagation_table needs the spectral path.");
    for (const auto& [names, option] : {std::pair {&epnl_files, "--epnl_file"},
                 std::pair {&band_frames_files, "--band_frames"}})
    {
        if (names->empty())
            continue;
        if (!spectral)
        {
            input_error(option, " needs the spectral path: the broadband path (--no-spectral) "
                        "carries no bands.");
        }
        if (!std::isnan(ambient_level))
        {
            input_error(option, " does not combine with -a: the ambient is an A-weighted level, "
                        "with no spectrum to add to the bands.");
        }
    }
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
    if (!atmosphere_files.empty() && !spectral)
    {
        input_error("--atmosphere needs the spectral path: layered absorption is per band, and "
                    "the broadband path (--no-spectral) carries one level attenuated by the "
                    "database's own EAA factor.");
    }
    if (!wind_texts.empty() && !atmosphere_files.empty())
    {
        input_error("--wind builds the atmosphere that --atmosphere would read; give one or the "
                    "other. (--write_atmosphere saves the built one in the file format.)");
    }
    if (!wind_texts.empty() && !spectral)
    {
        input_error("--wind needs the spectral path: the profile it builds is integrated per band, "
                    "and the broadband path (--no-spectral) carries one level attenuated by the "
                    "database's own EAA factor.");
    }
    if (wind_texts.empty())
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
    if (refraction != "none" && atmosphere_files.empty() && wind_texts.empty())
    {
        input_error("--refraction ", refraction, " needs --atmosphere or --wind: it refracts "
                    "through the profile's temperature and wind.");
    }
    if (!(shadow_floor >= 0.0))
        input_error("--shadow_floor must not be negative, got ", shadow_floor, ".");
    if (!(scattering >= 0.0))
        input_error("--scattering must not be negative, got ", scattering, ".");
    if (shadow_model == "resolved" && (refraction != "stratified" || !ground_effect))
    {
        input_error("--shadow_model resolved sets the shadow boundaries --refraction stratified casts, "
                    "through the shadow-boundary field --ground_effect's ground model sets up; give both.");
    }
    if (shadow_model != "resolved" && (shadow_phase_option->count() > 0 || cramped_ground_option->count() > 0))
        input_error("--shadow_phase and --cramped_ground shape --shadow_model resolved, which is not given.");
    if (scattering > 0.0 && (refraction == "none" || !ground_effect))
    {
        input_error("--scattering fills the shadows --refraction casts, through the shadow-boundary "
                    "field that --ground_effect's ground model sets up; give both.");
    }
    if (!atmosphere_files.empty() && (!temperature_overrides.empty() || !humidity_overrides.empty()))
    {
        input_error("--atmosphere gives temperature and humidity at every height; --temperature "
                    "and --humidity would contradict it.");
    }

    // --terrain, and what it does to the observers.
    std::shared_ptr<const Terrain> terrain;
    if (terrain_file.empty())
    {
        if (terrain_offset_option->count() > 0 || terrain_step_option->count() > 0 || observers_on_terrain)
            input_error("--terrain_offset, --terrain_step and --observers_on_terrain need --terrain.");
    }
    else
    {
        // Ground-reflected rays, shadow boundaries and the table's receiver
        // heights are all built over the plane z = 0; handing them a tilted
        // plane per path would be a different model, not this one.
        if (refraction != "none")
        {
            input_error("--terrain does not combine with --refraction: the ray models reflect off, "
                        "and cast their shadows onto, the flat ground z = 0. Use --terrain with "
                        "straight rays, or --refraction over flat ground.");
        }
        if (!propagation_table.empty())
        {
            input_error("--terrain does not combine with --propagation_table: the table's source "
                        "and receiver heights are over flat ground.");
        }
        if (observers_on_terrain && !scattered)
            input_error("--observers_on_terrain is for an observer file (-b); grid observers always "
                        "stand on the terrain.");
        if (!frame)
        {
            input_error("--terrain needs to know where this run's x/y lie on the ground: give "
                        "--origin <latitude> <longitude> of the x/y origin (and --frame_bearing for "
                        "rotated axes -- the shipped AC tracks use 210), or use a geodetic "
                        "trajectory, which carries its own.");
        }
        if (terrain_step_option->count() > 0 && !(terrain_step > 0.0 && std::isfinite(terrain_step)))
            input_error("--terrain_step must be a positive length, not ", terrain_step / unit.to_meters, " ",
                    unit.name, ".");
        if (!std::isfinite(terrain_offset))
            input_error("--terrain_offset must be a number, not ", terrain_offset, ".");
        terrain = std::make_shared<const Terrain>(terrain_file);
        if (terrain_step_option->count() == 0)
            terrain_step = std::max(terrain->spacing_meters(), min_terrain_step);
        // The local frame's z is height above --ground_elevation on the
        // trajectory's datum; --terrain_offset puts the grid on that datum.
        const double height_offset = terrain_offset - ground_elevation;
        log("Terrain: ", terrain_file, " (vertical datum ", terrain->vertical_datum(), "; latitudes ",
                terrain->south(), " to ", terrain->north(), ", longitudes ", terrain->west(), " to ",
                terrain->east(), ", spacing ", terrain->spacing_meters(), " m). Ground z = height + "
                "--terrain_offset ", terrain_offset / unit.to_meters, " ", unit.name,
                " - --ground_elevation ", ground_elevation / unit.to_meters, " ", unit.name,
                "; each path's mean ground plane is fitted to samples every ",
                terrain_step / unit.to_meters, " ", unit.name, ".");
        // Each observer's ground, read from the grid itself rather than the
        // per-trajectory raster, which is a resampling of it.
        std::vector<double> ground(observers.size());
        std::size_t uncovered = 0;
        std::size_t first_uncovered = 0;
        for (std::size_t i = 0; i < observers.size(); i++)
        {
            const Geodetic where = frame->to_geodetic(Vec3(observers[i].x, observers[i].y, 0.0));
            ground[i] = terrain->elevation(where.latitude, where.longitude) + height_offset;
            if (!std::isfinite(ground[i]) && uncovered++ == 0)
                first_uncovered = i;
        }
        if (uncovered > 0)
        {
            const Geodetic where = frame->to_geodetic(Vec3(observers[first_uncovered].x,
                    observers[first_uncovered].y, 0.0));
            input_error("--terrain ", terrain_file, " has no height at ", uncovered, " of the ",
                    observers.size(), " observers, the first at x ", observers[first_uncovered].x / unit.to_meters,
                    ", y ", observers[first_uncovered].y / unit.to_meters, " ", unit.name, " (latitude ",
                    where.latitude, ", longitude ", where.longitude, "): outside the grid, which covers "
                    "latitudes ", terrain->south(), " to ", terrain->north(), " and longitudes ",
                    terrain->west(), " to ", terrain->east(), ", or on a missing node. Check --origin "
                    "and --frame_bearing.");
        }
        if (!scattered || observers_on_terrain)
        {
            for (std::size_t i = 0; i < observers.size(); i++)
                observers[i].z = ground[i];
            const auto [low, high] = std::minmax_element(ground.begin(), ground.end());
            log("  observers stand on the terrain, from ", *low / unit.to_meters, " to ",
                    *high / unit.to_meters, " ", unit.name, " in the local frame's z",
                    scattered ? " (the observer file's z is not used)." : ".");
        }
        else
        {
            // The file's z stands, but a surveyed microphone and a terrain grid on
            // one datum agree to the grid's own accuracy -- at Amedee 0.03 m
            // standard deviation once the 23.6 m datum offset is applied -- so a
            // larger disagreement is most likely a datum or frame mistake.
            constexpr double disagreement = 1.0; // m
            std::size_t off = 0;
            std::size_t worst = 0;
            std::vector<double> differences(observers.size());
            for (std::size_t i = 0; i < observers.size(); i++)
            {
                differences[i] = observers[i].z - ground[i];
                if (std::fabs(differences[i]) > disagreement)
                    off++;
                if (std::fabs(differences[i]) > std::fabs(differences[worst]))
                    worst = i;
            }
            std::vector<double> sorted = differences;
            std::nth_element(sorted.begin(), sorted.begin() + sorted.size() / 2, sorted.end());
            const double median = sorted[sorted.size() / 2];
            log("  observers keep the observer file's z; it sits ", median / unit.to_meters, " ",
                    unit.name, " from the terrain in the median, ", differences[worst] / unit.to_meters,
                    " at worst.");
            if (off > 0)
            {
                warn("--terrain and the observer file disagree by more than ", disagreement / unit.to_meters,
                        " ", unit.name, " at ", off, " of ", observers.size(), " observers (median ",
                        median / unit.to_meters, ", worst ", differences[worst] / unit.to_meters, " ",
                        unit.name, " at data row ", worst + 1, "). A uniform offset is a datum mismatch "
                        "(--terrain_offset, --ground_elevation); give --observers_on_terrain to stand "
                        "them on the terrain instead.");
            }
        }
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
    options.blend = blend_options;
    if (blend_options.load_factor_correction == LoadFactorCorrection::two_sided
            && min_load_factor > max_load_factor)
    {
        input_error("--min_load_factor ", min_load_factor, min_load_factor_option->count() > 0 ? "" : " (the default)",
                " is above --max_load_factor ", max_load_factor, ": the two-sided law's floor must not "
                "be above its ceiling.");
    }
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
    options.terrain = terrain;
    options.terrain_path = terrain_file;
    options.terrain_offset = terrain_offset;
    options.terrain_step = terrain_step;
    if (!propagation_table.empty())
    {
        options.propagation_table = std::make_shared<const PropagationTable>(PropagationTable::load(propagation_table));
        options.propagation_table_path = propagation_table;
        // A table pair takes its turbulence and roughness from the table; the
        // run's reach only the pairs the table does not cover (and the plate
        // offset's bare ground).  Either way round, a difference would treat
        // the two kinds of pair differently without a word.
        const TableCoherence table_has = table_coherence(options.propagation_table->provenance());
        const std::string mismatch = coherence_mismatch(ground_coherence(options.ground), table_has);
        if (!mismatch.empty())
        {
            warn("--propagation_table ", propagation_table, ": ", mismatch, ". Build the table with the "
                    "run's settings (niceops_propagation_table --turbulence, --gamma_t, --roughness), or give "
                    "the run the table's.");
        }
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
    // Each trajectory's air (--atmosphere, --wind, --humidity, --temperature,
    // --pressure, --start_utc) goes into its own copy of the options, in the
    // trajectory loop.
    const auto value_for = [](const std::vector<double>& list, std::size_t path) {
        return list.empty() ? std::numeric_limits<double>::quiet_NaN() : list[path];
    };
    options.refraction = refraction != "none";
    options.stratified = refraction == "stratified";
    if ((tunneling || skip_fill) && !options.stratified)
        input_error("--tunneling and --skip_fill need --refraction stratified.");
    options.waves.tunneling = tunneling;
    options.waves.skip_fill = skip_fill;
    options.shadow_floor = shadow_floor;
    options.resolved_shadow = shadow_model == "resolved";
    options.shadow_phase = shadow_phase;
    options.cramped_ground = cramped_ground == "ray" ? CrampedGround::ray
            : cramped_ground == "straight"           ? CrampedGround::straight
                                                     : CrampedGround::series;
    options.scattering = scattering;
    options.frame_bearing = frame_bearing;
    if (!speed_reference.empty())
        options.speed_reference = parse_speed_reference(speed_reference);
    if (!nose.empty())
        options.nose = parse_nose_source(nose);
    // Each distinct --wind and --atmosphere is parsed once, and its warnings
    // given once.
    std::vector<std::optional<SurfaceWind>> winds(path_count);
    std::vector<std::optional<AtmosphereProfile>> profiles(path_count);
    {
        std::map<std::string, SurfaceWind> parsed_winds;
        std::map<std::string, AtmosphereProfile> read_profiles;
        bool warned_bearing = false;
        for (std::size_t path = 0; path < path_count; path++)
        {
            if (!wind_texts.empty())
            {
                const std::string& wind_text = wind_texts[path];
                auto found = parsed_winds.find(wind_text);
                if (found == parsed_winds.end())
                {
                    SurfaceWind wind = parse_surface_wind(wind_text, unit);
                    wind.roughness_length = roughness_length;
                    wind.surface_layer_depth = surface_layer_depth;
                    wind.stability = parse_stability(stability);
                    wind.levels_per_decade = profile_levels;
                    if (const auto problem = surface_wind_problem(wind, unit); !problem.empty())
                        input_error("--wind: ", problem);
                    if (wind.stability != Stability::neutral && parsed_winds.empty())
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
                    // Once, for the first wind that blows: a calm one has no direction.
                    if (wind.speed > 0.0 && refraction != "none" && std::isnan(frame_bearing) && !frame
                            && !warned_bearing)
                    {
                        warned_bearing = true;
                        warn("--wind: the direction is a compass bearing, turned into the local frame taking "
                             "+x as east. If this run's x/y are rotated -- the shipped AC tracks run +x along "
                             "210 degrees -- give --frame_bearing.");
                    }
                    found = parsed_winds.emplace(wind_text, wind).first;
                }
                winds[path] = found->second;
            }
            if (!atmosphere_files.empty())
            {
                const std::string& atmosphere_file = atmosphere_files[path];
                auto found = read_profiles.find(atmosphere_file);
                if (found == read_profiles.end())
                    found = read_profiles.emplace(atmosphere_file,
                            read_atmosphere_profile(atmosphere_file, unit, ground_elevation)).first;
                if (found->second.has_pressure() && !std::isnan(value_for(pressure_overrides, path)))
                {
                    input_error("--atmosphere ", atmosphere_file, " gives pressure at every height; "
                            "--pressure would contradict it. (Without a pressure column, --pressure "
                            "sets the surface pressure the profile is integrated from.)");
                }
                profiles[path] = found->second;
            }
        }
    }
    options.trim_start = trim_start;
    options.utc_offset = utc_offset;
    if (!laeq_window_hours.empty())
        options.laeq_window = {laeq_window_hours[0], laeq_window_hours[1]};
    options.scattered = scattered;
    options.want_day_metric = want_day_metric;
    options.want_lamax = want_lamax;
    options.want_nax = want_nax;
    options.nax_threshold = nax_threshold;
    options.want_ta = want_ta;
    options.ta_threshold = ta_threshold;
    options.lamax_weighting = level_weighting;
    options.want_csel = want_csel;
    options.want_coverage = want_coverage;
    options.want_provenance = want_provenance;
    options.full_history = full_history;
    options.epnl_masking = epnl_masking ? ToneMasking::below_noy_floor : ToneMasking::none;
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
    // -s is the resampling step for both adapt() and a geodetic track that
    // needs resampling; neither ever coarsens a track that is already finer.
    options.geodetic.resample_dt = dt;
    options.geodetic.path_smoothing = path_smoothing;

    std::size_t failed_paths = 0;
    // The operations the run totals stand for: each finished trajectory's
    // weight, summed.
    double operations_run = 0.0;

    for (std::size_t path_index = 0; path_index < path_count; path_index++)
    {
        try
        {
            const std::string& key = database_keys[path_index];
            auto cached = databases.find(key);
            if (cached == databases.end())
            {
                log("Loading database ", database_files[path_index], " for ", path_files[path_index]);
                // try_emplace enters it only once it has loaded.  An entry for a
                // database that failed to open would let a later trajectory naming
                // the same file skip the load and run on nothing.
                cached = databases.try_emplace(key, database_files[path_index], database_load).first;
                if (databases.size() > 1)
                    info(databases.size(), " databases held in memory.");
                if (auto& request = manifest_request())
                    request->manifest.set_database_attributes(database_files[path_index],
                            cached->second.root_attributes());
            }
            PathJob job;
            job.trajectory = path_files[path_index];
            job.database = database_files[path_index];
            job.name = run_names[path_index];
            job.epnl = epnl_files.empty() ? std::string() : epnl_files[path_index];
            job.band_frames = band_frames_files.empty() ? std::string() : band_frames_files[path_index];
            job.vehicle = vehicle_files[path_index];
            job.output = out_files[path_index];
            job.state = state_files.empty() ? std::string() : state_files[path_index];
            job.atmosphere_out = atmosphere_out_files.empty() ? std::string()
                                                              : atmosphere_out_files[path_index];
            job.csel_output = csel_files.empty() ? std::string() : csel_files[path_index];
            job.count = path_count;
            RunOptions path_options = options;
            path_options.atmosphere = profiles[path_index];
            path_options.surface_wind = winds[path_index];
            path_options.humidity_override = value_for(humidity_overrides, path_index);
            path_options.temperature_override = value_for(temperature_overrides, path_index);
            path_options.pressure_override = value_for(pressure_overrides, path_index);
            path_options.geodetic.start_utc = value_for(start_utcs, path_index);
            if (operations)
            {
                job.weight = operations->weights[path_index];
                job.operations = true;
                path_options.geodetic.clock_replaces_track = operations->clock_replaces_track[path_index] != 0;
                if (path_options.geodetic.clock_replaces_track)
                    info("local_time dates ", job.name, " at ", path_options.geodetic.start_utc,
                            " s UTC, in place of any clock its track carries.");
            }
            trajectory_holder = process_trajectory(cached->second, job, path_options, observers,
                    observer_set.microphone_heights, x_obs, y_obs, integrated_levels, totals);
            operations_run += job.weight;
        }
        catch (const std::exception& problem)
        {
            failed_paths++;
            report_again(problem, "Skipping " + (operations ? operations->labels[path_index]
                                                            : path_files[path_index]) + ": ");
        }
        // Done with it once no later trajectory names it; the front one stays.
        const std::string& key = database_keys[path_index];
        if (last_use.at(key) == path_index && key != front_key)
            databases.erase(key);
    }
    // `database` refers to the front entry, which the provenance reads below;
    // erasing it would leave a dangling reference that no output shows.
    assert(databases.count(front_key) == 1 && &databases.find(front_key)->second == &database);

    if (failed_paths > 0)
    {
        log(failed_paths, " of ", path_count, " trajectories could not be processed.");
    }

    // "over 3 operations" in the logs below: the weighted count, with the
    // runs it came from when the two differ.
    const std::size_t runs_finished = path_count - failed_paths;
    const auto operations_phrase = [&](const char* after) {
        std::ostringstream phrase;
        phrase << operations_run << (operations_run == 1.0 ? " operation" : " operations");
        if (operations_run != static_cast<double>(runs_finished))
            phrase << " (" << runs_finished << (runs_finished == 1 ? " run)" : " runs)");
        phrase << after;
        return phrase.str();
    };

    std::array<std::vector<double>, day_metric_count> day_levels;
    for (std::size_t m = 0; m < day_metric_count; m++)
    {
        if (!want_day_metric[m])
            continue;
        // Each is an energy average: over a 24-hour day for the penalized
        // levels, the summed event energies over 86 400 s, not their raw sum
        // (which reads 49.37 dB high); over its own period for a period level.
        const double seconds = averaging_seconds(static_cast<DayMetric>(m), options.laeq_window);
        const auto& energy = totals.day_energy[m];
        const std::string& file = *day_files[m];
        std::vector<double> level(energy.size());
        for (std::size_t i = 0; i < energy.size(); i++)
            level[i] = 10.0 * std::log10(energy[i] / seconds);
        day_levels[m] = level;
        log("Writing the ", day_metric_name(static_cast<DayMetric>(m)), " over ", operations_phrase(""),
                full_history ? ", each integrated over its whole history," : ", each over its own 10 dB-down window,",
                " to ", file);
        write_levels(x_obs, y_obs, level, file, scattered, unit);
        if (latlon)
            write_geodetic_levels(x_obs, y_obs, level, *frame, latlon_companion(file), scattered);
    }

    if (want_lamax)
    {
        log("Writing LAmax (", maximum_name(level_weighting), ") over ", operations_phrase(" to "), lamax_file);
        write_levels(x_obs, y_obs, totals.lamax, lamax_file, scattered, unit);
        if (latlon)
            write_geodetic_levels(x_obs, y_obs, totals.lamax, *frame, latlon_companion(lamax_file), scattered);
    }

    if (want_nax)
    {
        log("Writing NAx (operations with ", maximum_name(level_weighting), " at or above ", nax_threshold,
                " dBA) over ", operations_phrase(" to "), nax_file);
        write_levels(x_obs, y_obs, totals.nax, nax_file, scattered, unit);
        if (latlon)
            write_geodetic_levels(x_obs, y_obs, totals.nax, *frame, latlon_companion(nax_file), scattered);
    }

    if (want_ta)
    {
        log("Writing the time above ", ta_threshold, " dBA (", to_string(level_weighting),
                " weighting, seconds of reception time) over ", operations_phrase(" to "), ta_file);
        write_levels(x_obs, y_obs, totals.time_above, ta_file, scattered, unit);
        if (latlon)
            write_geodetic_levels(x_obs, y_obs, totals.time_above, *frame, latlon_companion(ta_file), scattered);
    }

    if (want_provenance)
    {
        log("Writing where each observer's energy came from, over ", operations_phrase(", to "), provenance_file);
        // The front database: the provenance refuses a run of several.
        assert(!several_databases);
        write_provenance(observers, totals, database, provenance_file, unit);
        if (!provenance_spheres_file.empty())
        {
            SphereEnergyFile meta;
            meta.grid = !scattered;
            for (std::size_t i = 0; i < totals.provenance.size(); i++)
            {
                meta.x.push_back(observers[i].x);
                meta.y.push_back(observers[i].y);
                meta.exposure.push_back(totals.provenance[i][0]);
            }
            for (std::size_t k = 0; k < database.condition_count(); k++)
                meta.conditions.push_back(database.condition(k));
            meta.database = database_files.front();
            log("Writing each sphere's share of each observer's energy to ", provenance_spheres_file);
            write_sphere_energy(provenance_spheres_file, meta, totals.sphere_energy);
        }
    }

    // The energy-weighted fraction of each observer's levels that came from
    // cells the database did not measure.  An observer that heard nothing has
    // no fraction to give, and is written as NaN.
    std::vector<double> unmeasured_share;
    if (want_coverage)
    {
        const auto& coverage = totals.coverage;
        unmeasured_share.resize(coverage.size());
        for (std::size_t i = 0; i < coverage.size(); i++)
        {
            unmeasured_share[i] = coverage[i][1] > 0.0 ? coverage[i][0] / coverage[i][1]
                                                       : std::numeric_limits<double>::quiet_NaN();
        }
        log("Writing the unmeasured share of the levels over ", operations_phrase(" to "), coverage_file);
        write_levels(x_obs, y_obs, unmeasured_share, coverage_file, scattered, unit);
        if (latlon)
        {
            write_geodetic_levels(x_obs, y_obs, unmeasured_share, *frame, latlon_companion(coverage_file),
                    scattered);
        }
    }

    // Every gridded level this run wrote, each as it went to its own file, for
    // --contour_file and --raster_file.  A skipped trajectory writes no -o, and
    // its SEL vector holds only zeros, so the SEL is here only when the one
    // trajectory ran.
    std::vector<RasterField> gridded;
    if (path_count == 1 && failed_paths == 0)
        gridded.push_back({"sel", "sound exposure level (SEL), A-weighted", "dB", &integrated_levels});
    for (std::size_t m = 0; m < day_metric_count; m++)
    {
        const auto metric = static_cast<DayMetric>(m);
        if (want_day_metric[m])
            gridded.push_back({day_metric_key(metric), day_metric_name(metric), "dB", &day_levels[m]});
    }
    if (want_lamax)
        gridded.push_back({"lamax", std::string("LAmax, ") + maximum_name(level_weighting), "dB", &totals.lamax});
    if (want_nax)
    {
        std::ostringstream name;
        name << "operations with " << maximum_name(level_weighting) << " at or above " << nax_threshold << " dB";
        gridded.push_back({"nax", name.str(), "1", &totals.nax});
    }
    if (want_ta)
    {
        std::ostringstream name;
        name << "time above " << ta_threshold << " dB (" << to_string(level_weighting)
             << " weighting, reception time)";
        gridded.push_back({"time_above", name.str(), "s", &totals.time_above});
    }
    if (want_coverage)
    {
        gridded.push_back({"unmeasured_share", "energy share of the level from cells the database did not "
                "measure", "1", &unmeasured_share});
    }

    if (!contour_file.empty())
    {
        // The metric's own output was required up front, and SEL with one
        // trajectory, so only a skipped trajectory's SEL can be missing.
        const auto field = std::find_if(gridded.begin(), gridded.end(),
                [&](const RasterField& f) { return f.name == contour_metric; });
        if (field == gridded.end())
        {
            warn("--contour_file ", contour_file, " is not written: its trajectory was skipped, so there is no "
                    "SEL to contour.");
        }
        else
        {
            // The figures' levels, as D's iota(cmin, cmax + cstep, cstep) runs them.
            const auto levels = iota_exclusive(contour_min, contour_max + contour_step, contour_step);
            const auto areas = write_contour_geojson(contour_file, x_obs, y_obs, *field->values, levels, *frame,
                    contour_metric);
            std::ostringstream summary;
            for (std::size_t k = 0; k < levels.size(); k++)
                summary << (k ? ", " : "") << levels[k] << " dB: " << areas[k] / (unit.to_meters * unit.to_meters)
                        << " " << unit.name << "^2";
            log("Wrote ", contour_metric, " contours to ", contour_file, "; areas at or above ", summary.str(),
                    ".");
        }
    }

    if (!raster_file.empty())
    {
        if (gridded.empty())
        {
            warn("--raster_file ", raster_file, " has no gridded level to hold: SEL goes in it for one "
                    "trajectory, and a day or period level, LAmax, NAx, TA or coverage for several.");
        }
        log("Writing ", gridded.size(), gridded.size() == 1 ? " gridded level" : " gridded levels", " to ",
                raster_file, frame ? ", with WGS84 positions." : " (no frame: local x/y only).");
        write_raster(raster_file, x_obs, y_obs, gridded, frame ? &*frame : nullptr, frame_bearing);
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
    int status = 0;
    try
    {
        status = run(argc, argv);
    }
    catch (const std::exception& problem)
    {
        report_again(problem, "Stopping: ");
        status = 2;
    }
    if (auto& request = manifest_request())
    {
        const double seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - request->start).count();
        try
        {
            request->manifest.write(request->path, status, seconds);
        }
        catch (const std::exception& problem)
        {
            report_again(problem, "Stopping: ");
            status = 2;
        }
    }
    return status;
}
