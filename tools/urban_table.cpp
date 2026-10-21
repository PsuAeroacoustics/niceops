// Build an urban table for niceops --urban_table (src/urban.hpp).
//
//   niceops_urban_table --canyons receivers.csv --database databases/B407_ambient_gated.nod
//       --frame_bearing 210 -o tables/street.nc
//
// Each row of --canyons is one microphone standing in an ideal street canyon,
// by column name, lengths in --units:
//   x, y, z          the observer, as the run's observer file gives it (a
//                    length column may name its unit, x_ft or x_m, which must
//                    then be --units; an id column is ignored)
//   height           the microphone above the ground there (the observer
//                    file's fifth column; the run's microphone is z + height)
//   axis             compass bearing of the street's axis, degrees (either way
//                    along it)
//   width            between the two walls; 0 for an open microphone, with no
//                    street around it (0 dB and no facades; the rest of the
//                    row but its position is not read)
//   offset           from wall A, the wall on the left looking along `axis`
//   height_a,        the walls' tops above the ground
//   height_b
//   absorption_a,    optional: each wall's energy absorption coefficient,
//   absorption_b     every band (default --absorption)
// Or, with --buildings and --receivers, microphones among LoD1 buildings
// (src/building_paths.hpp):
//
//   niceops_urban_table --buildings district.csv --receivers mics.csv --database ... -o tables/district.nc
//
// --buildings may also be GeoJSON (.geojson or .json): footprints in WGS84,
// projected into the run's frame by --origin and --frame_bearing, their
// holes courtyards, their heights in metres from OpenStreetMap's tags
// (src/geojson_buildings.hpp).  As a CSV, --buildings has one row per
// footprint vertex, by column name, lengths in
// --units: building (an id; a building's rows consecutive, its vertices in
// order around it), x, y, top (its roof above the ground), and optionally
// absorption (its walls' energy absorption coefficient, the same on each of
// its rows; default --absorption) and scattering (the share of what its walls
// reflect that they scatter diffusely, likewise; default --facade_scattering).  --receivers has the observer file's x, y,
// z and height; the buildings stand on the ground at each receiver's z.  Each
// receiver's facades within --facade_radius that face it are written for the
// run's first-order images, and left out of its table.
//
// The table is the canyon's under plane-wave incidence (canyon_ratio): the
// specular rays from above and off the floor and the roof edges' diffracted
// rays, summed coherently and band-averaged, but the first-order reflections,
// relative to open ground; the run adds the first-order reflections itself,
// from the two walls written beside the table, each --facade_length long.  The file is SI and records --frame_bearing, the
// compass bearing of the local +x axis, which the run checks against its own.

#include <CLI/CLI.hpp>

#include <algorithm>
#include <chrono>
#include <limits>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "building_paths.hpp"
#include "building_input.hpp"
#include "buildings.hpp"
#include "canyon_wave.hpp"
#include "geodetic.hpp"
#include "geojson_buildings.hpp"
#include "gpu/gpu.hpp"
#include "propagation_table.hpp"
#include "sphere_database.hpp"
#include "urban.hpp"
#include "util/constants.hpp"
#include "util/csv.hpp"
#include "util/parallel.hpp"

using namespace niceops; // NOLINT

int main(int argc, char** argv)
{
    CLI::App app {"Build an urban table for niceops --urban_table."};
    std::string canyons_file;
    std::string buildings_file;
    std::string receivers_file;
    std::size_t max_order = 64;
    std::size_t max_diffracted_order = 6;
    double radius = 150.0;        // m unless given
    double facade_radius = 100.0; // m unless given
    double facade_scattering = 0.0;
    double patch_size = 4.0;      // m unless given
    double diffuse_radius = 100.0; // m unless given
    std::vector<double> origin;
    double level_height = 3.0;
    double default_height = 0.0;
    double orientation_tolerance = 1.0; // degrees
    std::string database_file;
    std::string output;
    std::string units = "ft";
    std::vector<double> bands;
    std::vector<double> elevations;
    int azimuths = 72;
    double frame_bearing = 90.0;
    double absorption = 0.0;
    double floor_absorption = 0.0;
    double facade_length = 10000.0;
    double sound_speed = 343.2;
    bool no_diffraction = false;
    double wave_below = 0.0;
    std::size_t wave_blend = 2;
    double wave_ppw = 15.0;
    int threads = 0;
    bool use_gpu = true;
    app.add_option("--canyons", canyons_file, "Receivers in ideal canyons, one per row (see the header).");
    app.add_option("--buildings", buildings_file,
            "LoD1 buildings, one row per footprint vertex (see the header); with --receivers.");
    app.add_option("--receivers", receivers_file, "With --buildings: the microphones, as the run's observer file.");
    app.add_option("--max_order", max_order,
            "With --buildings: facade reflections on a specular ray, at most (default 64).");
    app.add_option("--max_diffracted_order", max_diffracted_order,
            "With --buildings: facade reflections after a diffraction, at most (default 6).");
    app.add_option("--radius", radius,
            "With --buildings: the edges, and the facades a diffracted ray reflects in, within this of a "
            "microphone, in --units (default 150 m).")
            ->check(CLI::PositiveNumber);
    app.add_option("--facade_scattering", facade_scattering,
            "With --buildings: the share of what a wall reflects that it scatters diffusely, where a building "
            "gives none (default 0, specular).")
            ->check(CLI::Range(0.0, 1.0));
    app.add_option("--patch_size", patch_size,
            "With --buildings and scattering walls: the walls' patches for the diffuse field, at most this "
            "size, in --units (default 4 m).")
            ->check(CLI::PositiveNumber);
    app.add_option("--diffuse_radius", diffuse_radius,
            "With --buildings and scattering walls: the walls within this of a microphone scatter, in --units "
            "(default 100 m).")
            ->check(CLI::PositiveNumber);
    app.add_option("--origin", origin,
            "With GeoJSON --buildings: the run's WGS84 origin, latitude longitude in degrees, to project the "
            "footprints into its frame (with --frame_bearing).")
            ->expected(2);
    app.add_option("--level_height", level_height,
            "With GeoJSON --buildings: metres per storey, for buildings that give only building:levels (default 3).")
            ->check(CLI::PositiveNumber);
    app.add_option("--default_height", default_height,
            "With GeoJSON --buildings: the height, metres, of buildings that give none (default: refuse them).")
            ->check(CLI::PositiveNumber);
    app.add_option("--orientation_tolerance", orientation_tolerance,
            "With --buildings: walls whose orientations agree within this, degrees, reflect alike, so a ray "
            "bounces up a street a little out of parallel as up the canyon it nearly is (default 1).")
            ->check(CLI::NonNegativeNumber);
    app.add_option("--facade_radius", facade_radius,
            "With --buildings: the facades written for the run's first-order images, within this of a "
            "microphone, in --units (default 100 m).")
            ->check(CLI::NonNegativeNumber);
    app.add_option("--database", database_file, "Take the band centers from this database.");
    app.add_option("--bands", bands, "Band centers, Hz, instead of --database's.");
    app.add_option("--azimuths", azimuths, "Azimuths, uniform from 0 over 360 degrees (default 72).")
            ->check(CLI::PositiveNumber);
    app.add_option("--elevations", elevations,
            "Elevations, degrees, ascending (default 0 to 30 by 1.5, then 35 to 90 by 5).");
    app.add_option("--frame_bearing", frame_bearing,
            "Compass bearing of the run's local +x axis, degrees (default 90, +x east).");
    app.add_option("--absorption", absorption, "Walls' energy absorption coefficient where a row gives none "
                                               "(default 0, rigid).")
            ->check(CLI::Range(0.0, 1.0));
    app.add_option("--floor_absorption", floor_absorption,
            "The street's energy absorption coefficient, for the rays the floor reflects (default 0, rigid "
            "pavement).")
            ->check(CLI::Range(0.0, 1.0));
    app.add_option("--facade_length", facade_length,
            "Length of each wall written for the run's first-order images, in --units (default 10 km).")
            ->check(CLI::PositiveNumber);
    app.add_option("--sound_speed", sound_speed, "For the diffracted path's wavelengths, m/s (default 343.2).")
            ->check(CLI::PositiveNumber);
    app.add_flag("--no_diffraction", no_diffraction,
            "Leave out the roof edges' diffraction and sum the specular rays in energy (for comparison).");
    app.add_option("--wave_below", wave_below,
            "Solve each flush-roof canyon by the wave equation (src/canyon_wave) and take its table for the "
            "bands centered at or below this, Hz (default 0, none).")
            ->check(CLI::NonNegativeNumber);
    app.add_option("--wave_blend", wave_blend,
            "Bands above --wave_below over which the wave table fades into the geometric one (default 2).");
    app.add_option("--wave_ppw", wave_ppw, "The wave solver's grid points per wavelength (default 15).")
            ->check(CLI::Range(6.0, 40.0));
    app.add_flag("--gpu,!--no-gpu", use_gpu,
            "March --wave_below's wave equation on the GPU (Metal on macOS), in single precision, two canyons at a "
            "time: on by default wherever there is one. --gpu asks for it, and is refused where there is none or "
            "without --wave_below; --no-gpu keeps the CPU.");
    app.add_option("-j,--threads", threads,
            "Worker threads, for the geometric table's directions and the canyons the wave equation solves "
            "(default: NICEOPS_THREADS when set, else one per hardware thread); the table is the same bits at "
            "any count.");
    app.add_option("--units", units, "Length unit: ft (default) or m.")->check(CLI::IsMember({"ft", "m"}));
    app.add_option("-o,--output", output, "The table file to write.")->required();
    CLI11_PARSE(app, argc, argv);

    try
    {
        const double to_meters = units == "m" ? 1.0 : ft2m;
        if (canyons_file.empty() == buildings_file.empty())
            throw std::runtime_error("give either --canyons or --buildings");
        // Before any work: the GPU marches the wave equation and nothing else,
        // wherever there is one; asked for, it must be there and have one.
        // The flag's count takes --no-gpu too, so asking is its value as well.
        if (use_gpu && app.count("--gpu") > 0)
        {
            if (wave_below <= 0.0)
                throw std::runtime_error("--gpu marches --wave_below's wave equation, and there is none to march");
            gpu::Device::require("--gpu");
        }
        use_gpu = use_gpu && wave_below > 0.0 && gpu::Device::get() != nullptr;
        if (buildings_file.empty() != receivers_file.empty())
            throw std::runtime_error("--receivers goes with --buildings, and --buildings needs it");
        // Defaults are lengths in metres; a value given is in --units.
        if (app.count("--radius") > 0)
            radius *= to_meters;
        if (app.count("--facade_radius") > 0)
            facade_radius *= to_meters;
        if (app.count("--patch_size") > 0)
            patch_size *= to_meters;
        if (app.count("--diffuse_radius") > 0)
            diffuse_radius *= to_meters;
        if (!bands.empty() == !database_file.empty())
            throw std::runtime_error("give either --database or --bands");
        if (!database_file.empty())
        {
            const Database database(database_file);
            if (database.band_centers().empty())
                throw std::runtime_error(database_file + " carries no spectra");
            bands = database.band_centers();
        }
        if (elevations.empty())
        {
            for (int k = 0; k <= 20; k++)
                elevations.push_back(1.5 * k);
            for (int k = 7; k <= 18; k++)
                elevations.push_back(5.0 * k);
        }
        UrbanTable::Axes axes;
        for (int k = 0; k < azimuths; k++)
            axes.azimuths.push_back(360.0 * k / azimuths);
        axes.elevations = elevations;
        axes.bands = bands;

        auto column_of = [&](const LabeledCsv& csv, const std::string& file, const std::string& name, bool required,
                                 bool length) { return length_column(csv, file, units, name, required, length); };
        auto read_rows = [](const std::string& file, const std::string& what) { return read_numeric_rows(file, what); };
        const std::size_t nb = bands.size();

        if (!buildings_file.empty())
        {
            if (wave_below > 0.0)
                throw std::runtime_error("--wave_below solves ideal canyons only, not --buildings");
            BuildingInputOptions input;
            input.units = units;
            input.to_meters = to_meters;
            input.origin = origin;
            input.frame_bearing = frame_bearing;
            input.level_height = level_height;
            input.default_height = default_height;
            input.default_height_given = app.count("--default_height") > 0;
            input.absorption = absorption;
            input.scattering = facade_scattering;
            const BuildingInput district = read_buildings(buildings_file, input);
            const std::vector<Footprint>& footprints = district.footprints;
            const std::vector<double>& building_absorption = district.absorption;
            const std::vector<double>& building_scattering = district.scattering;
            const std::vector<std::string>& ids = district.ids;
            Buildings city;
            try
            {
                city = Buildings(footprints, std::max(orientation_tolerance, 1e-6) * pi / 180.0);
            }
            catch (const std::invalid_argument& e)
            {
                throw std::runtime_error(buildings_file + ": " + e.what() + " (counting buildings from 1 in file order)");
            }
            if (city.shared_sides() > 0)
            {
                std::cerr << "Buildings stand against each other along " << city.shared_sides() / 2
                          << " shared walls." << std::endl;
            }
            // Each facade's reflection and scattering, every band: its
            // building's walls.  The rays, and the run's own reflections,
            // carry the specular share, rho (1 - s).
            std::vector<double> reflection(city.facades().size() * nb), scattering(reflection.size());
            std::vector<double> specular(reflection.size());
            bool scatters = false;
            for (std::size_t f = 0; f < city.facades().size(); f++)
                for (std::size_t b = 0; b < nb; b++)
                {
                    const std::size_t building = city.facade_building(f);
                    reflection[f * nb + b] = 1.0 - building_absorption[building];
                    scattering[f * nb + b] = building_scattering[building];
                    specular[f * nb + b] = reflection[f * nb + b] * (1.0 - scattering[f * nb + b]);
                    scatters = scatters || building_scattering[building] > 0.0;
                }
            DiffuseOptions diffuse_options;
            diffuse_options.radius = diffuse_radius;
            diffuse_options.patch_size = patch_size;
            diffuse_options.max_order = max_order;

            const LabeledCsv rcsv = read_rows(receivers_file, "receivers");
            const std::size_t rx = column_of(rcsv, receivers_file, "x", true, true);
            const std::size_t ry = column_of(rcsv, receivers_file, "y", true, true);
            const std::size_t rz = column_of(rcsv, receivers_file, "z", true, true);
            const std::size_t rh = column_of(rcsv, receivers_file, "height", true, true);
            BuildingPathOptions options;
            options.max_order = max_order;
            options.max_diffracted_order = max_diffracted_order;
            options.radius = radius;
            options.diffraction = !no_diffraction;
            std::vector<UrbanTable::Receiver> receivers;
            std::vector<double> gains;
            std::size_t truncated = 0, unconverged = 0, patches = 0;
            const auto start = std::chrono::steady_clock::now();
            for (std::size_t row = 0; row < rcsv.rows.size(); row++)
            {
                const auto& r = rcsv.rows[row];
                const std::string where = receivers_file + " row " + std::to_string(row + 1) + ": ";
                const double x = r[rx] * to_meters, y = r[ry] * to_meters, floor = r[rz] * to_meters;
                const double height = r[rh] * to_meters;
                if (!(height >= 0.0))
                    throw std::runtime_error(where + "the microphone's height must not be negative");
                // The builder's frame: the local x and y, z above this ground.
                // A microphone on a wall takes the wall's reflection there, and
                // the table carries its first-order reflections too
                // (BuildingReceiver); the tables have no roofs, so one on or
                // over a roof is refused, as is one on an edge.
                const Vec3 local(x, y, height);
                const Buildings::Support on = city.support(local);
                if (on.kind == Buildings::Support::Kind::edge)
                    throw std::runtime_error(where + "the microphone lies on a building's edge, or where two of its "
                                                     "surfaces meet: move it off by a millimetre or more");
                const bool on_wall = on.kind == Buildings::Support::Kind::facade;
                if (!on_wall && city.inside(x, y) < city.size())
                    throw std::runtime_error(where + "the microphone stands inside building "
                            + ids[city.inside(x, y)]
                            + " or over its roof, which the tables do not model (niceops --buildings does)");
                const auto run = on_wall ? std::vector<std::size_t> {} : facing_facades(city, local, facade_radius);
                const BuildingReceiver receiver(city, local, run, options);
                truncated += receiver.truncated();
                std::unique_ptr<DiffuseFacades> diffuse;
                if (scatters)
                {
                    diffuse = std::make_unique<DiffuseFacades>(city, local, reflection, scattering,
                            1.0 - floor_absorption, bands, sound_speed, diffuse_options, tool_thread_count(threads));
                    unconverged += !diffuse->converged();
                    patches += diffuse->patches().size();
                }
                const auto row_gains = building_gains(receiver, axes, frame_bearing, specular,
                        1.0 - floor_absorption, sound_speed, tool_thread_count(threads), diffuse.get());
                gains.insert(gains.end(), row_gains.begin(), row_gains.end());
                UrbanTable::Receiver out;
                out.microphone = Vec3(x, y, floor + height);
                for (std::size_t f : run)
                {
                    Facade facade = city.facades()[f];
                    for (double* z : {&facade.bottom, &facade.top, &facade.top_start, &facade.top_end,
                                 &facade.bottom_start, &facade.bottom_end})
                        *z += floor;
                    out.facades.push_back(facade);
                    out.reflection.insert(out.reflection.end(), specular.begin() + static_cast<std::ptrdiff_t>(f * nb),
                            specular.begin() + static_cast<std::ptrdiff_t>((f + 1) * nb));
                }
                receivers.push_back(std::move(out));
            }
            std::cerr << "Traced " << receivers.size() << " receivers among " << city.size() << " buildings in "
                      << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << " s."
                      << std::endl;
            if (scatters)
                std::cerr << "The walls' diffuse field: " << patches << " patches over the receivers." << std::endl;
            if (unconverged > 0)
                std::cerr << "niceops_urban_table: warning: for " << unconverged
                          << " receivers the diffuse field's solve did not converge; give the walls some absorption."
                          << std::endl;
            if (truncated > 0)
                std::cerr << "niceops_urban_table: warning: for " << truncated
                          << " receivers the diffracted rays' image search stopped at its limit; lower "
                             "--max_diffracted_order or --radius."
                          << std::endl;
            std::map<std::string, std::string> provenance {
                    {"buildings", buildings_file},
                    {"receivers", receivers_file},
                    {"diffraction",
                            no_diffraction ? "none; specular rays in energy" : "UTD at roof and corner edges, coherent"},
                    {"max_order", std::to_string(max_order)},
                    {"orientation_tolerance_degrees", exact_text(orientation_tolerance)},
                    {"max_diffracted_order", std::to_string(max_diffracted_order)},
                    {"radius_m", exact_text(radius)},
                    {"facade_radius_m", exact_text(facade_radius)},
                    {"facade_scattering", scatters ? "per building (see the buildings file), default "
                                    + exact_text(facade_scattering) : "none"},
                    {"sound_speed_mps", exact_text(sound_speed)},
                    {"floor_absorption", exact_text(floor_absorption)},
            };
            if (origin.size() == 2)
                provenance["origin"] = exact_text(origin[0]) + " " + exact_text(origin[1]);
            if (scatters)
            {
                provenance["patch_size_m"] = exact_text(patch_size);
                provenance["diffuse_radius_m"] = exact_text(diffuse_radius);
            }
            const UrbanTable table(axes, std::move(receivers), std::move(gains), frame_bearing,
                    scatters ? "lod1_buildings+diffuse" : "lod1_buildings", provenance);
            table.save(output);
            std::cerr << "Wrote " << output << " (" << table.receivers().size() << " receivers x "
                      << axes.azimuths.size() << " azimuths x " << axes.elevations.size() << " elevations x " << nb
                      << " bands)." << std::endl;
            return 0;
        }

        const LabeledCsv csv = read_rows(canyons_file, "receivers");
        auto column = [&](const std::string& name, bool required, bool length) {
            return column_of(csv, canyons_file, name, required, length);
        };
        const std::size_t cx = column("x", true, true), cy = column("y", true, true), cz = column("z", true, true);
        const std::size_t ch = column("height", true, true), axis_c = column("axis", true, false);
        const std::size_t cw = column("width", true, true), co = column("offset", true, true);
        const std::size_t cha = column("height_a", true, true), chb = column("height_b", true, true);
        const std::size_t caa = column("absorption_a", false, false), cab = column("absorption_b", false, false);

        std::vector<UrbanTable::Receiver> receivers;
        std::vector<double> gains;
        // The wave band: its weight per band, and the canyons to solve.
        const std::vector<double> weight = wave_band_weights(bands, wave_below, wave_blend);
        std::size_t wave_bands = 0;
        for (std::size_t b = 0; b < nb; b++)
            if (weight[b] > 0.0)
                wave_bands = b + 1;
        struct WaveJob
        {
            std::size_t receiver;
            IdealCanyon canyon;
            double axis, alpha_a, alpha_b;
        };
        std::vector<WaveJob> jobs;
        for (std::size_t row = 0; row < csv.rows.size(); row++)
        {
            const auto& r = csv.rows[row];
            const std::string where = canyons_file + " row " + std::to_string(row + 1) + ": ";
            const double floor = r[cz] * to_meters;
            if (r[cw] == 0.0)
            {
                // No street: an open microphone, 0 dB everywhere and no facades.
                if (!(r[ch] >= 0.0))
                    throw std::runtime_error(where + "the microphone's height must not be negative");
                UrbanTable::Receiver open;
                open.microphone = Vec3(r[cx] * to_meters, r[cy] * to_meters, floor + r[ch] * to_meters);
                gains.insert(gains.end(), axes.azimuths.size() * axes.elevations.size() * nb, 0.0);
                receivers.push_back(std::move(open));
                continue;
            }
            IdealCanyon canyon;
            canyon.width = r[cw] * to_meters;
            canyon.offset = r[co] * to_meters;
            canyon.height_a = r[cha] * to_meters;
            canyon.height_b = r[chb] * to_meters;
            canyon.microphone_height = r[ch] * to_meters;
            const double alpha_a = caa == LabeledCsv::npos ? absorption : r[caa];
            const double alpha_b = cab == LabeledCsv::npos ? absorption : r[cab];
            if (!(canyon.width > 0.0 && canyon.offset > 0.0 && canyon.offset < canyon.width))
                throw std::runtime_error(where + "the offset must lie strictly between the walls");
            if (!(canyon.height_a > 0.0 && canyon.height_b > 0.0 && canyon.microphone_height >= 0.0))
                throw std::runtime_error(where + "wall heights must be positive and the microphone's not negative");
            if (!(alpha_a >= 0.0 && alpha_a <= 1.0 && alpha_b >= 0.0 && alpha_b <= 1.0))
                throw std::runtime_error(where + "absorption coefficients must lie within [0, 1]");
            const double axis = r[axis_c];
            if (!std::isfinite(axis))
                throw std::runtime_error(where + "the axis bearing is not finite");

            UrbanTable::Receiver receiver;
            receiver.microphone = Vec3(r[cx] * to_meters, r[cy] * to_meters, floor + canyon.microphone_height);
            receiver.facades = ideal_canyon_facades(canyon, receiver.microphone, floor, axis, frame_bearing,
                    facade_length * to_meters);
            const std::vector<double> reflect_a(nb, 1.0 - alpha_a), reflect_b(nb, 1.0 - alpha_b);
            receiver.reflection = reflect_a;
            receiver.reflection.insert(receiver.reflection.end(), reflect_b.begin(), reflect_b.end());
            const auto row_gains = ideal_canyon_gains(canyon, axis, axes, reflect_a, reflect_b, 1.0 - floor_absorption,
                    sound_speed, !no_diffraction, tool_thread_count(threads));
            gains.insert(gains.end(), row_gains.begin(), row_gains.end());
            if (wave_bands > 0)
            {
                // The wave solver needs flush roofs, and a microphone below them.
                if (std::abs(canyon.height_a - canyon.height_b) > 1e-3)
                    std::cerr << "niceops_urban_table: warning: " << where << "walls of unequal height ("
                              << r[cha] << " and " << r[chb] << " " << units
                              << "); the wave solver needs flush roofs, so its wave band is left geometric."
                              << std::endl;
                else if (canyon.microphone_height < canyon.height_a)
                    jobs.push_back({receivers.size(), canyon, axis, alpha_a, alpha_b});
            }
            receivers.push_back(std::move(receiver));
        }
        // Each distinct canyon once, several at a time, each on one thread.
        double worst_residual = -std::numeric_limits<double>::infinity();
        if (!jobs.empty())
        {
            const auto start = std::chrono::steady_clock::now();
            auto same = [](const WaveJob& a, const WaveJob& b) {
                return a.canyon.width == b.canyon.width && a.canyon.offset == b.canyon.offset
                        && a.canyon.height_a == b.canyon.height_a
                        && a.canyon.microphone_height == b.canyon.microphone_height && a.alpha_a == b.alpha_a
                        && a.alpha_b == b.alpha_b;
            };
            std::vector<std::size_t> solve_of(jobs.size());
            std::vector<std::size_t> distinct;
            for (std::size_t k = 0; k < jobs.size(); k++)
            {
                std::size_t d = 0;
                while (d < distinct.size() && !same(jobs[distinct[d]], jobs[k]))
                    d++;
                if (d == distinct.size())
                    distinct.push_back(k);
                solve_of[k] = d;
            }
            const double top = bands[wave_bands - 1] * std::pow(2.0, 1.0 / 6.0);
            std::vector<ApertureSpectrum> spectra(distinct.size());
            // On the GPU two canyons at a time: each march fills the device, and
            // one canyon's transforms (on the CPU) overlap the next's march,
            // the two sharing the run's threads.
            const unsigned workers = tool_thread_count(threads);
            const unsigned marching = use_gpu ? std::min(2u, workers) : workers;
            parallel_chunks(distinct.size(), marching,
                    [&](std::size_t begin, std::size_t end) {
                for (std::size_t d = begin; d < end; d++)
                {
                    const WaveJob& job = jobs[distinct[d]];
                    CanyonWaveSetup setup;
                    setup.width = job.canyon.width;
                    setup.depth = job.canyon.height_a;
                    setup.offset = job.canyon.offset;
                    setup.microphone_height = job.canyon.microphone_height;
                    setup.max_frequency = top;
                    setup.points_per_wavelength = wave_ppw;
                    setup.sound_speed = sound_speed;
                    setup.wall_absorption = job.alpha_a;
                    setup.wall_absorption_b = job.alpha_b;
                    setup.gpu = use_gpu;
                    if (use_gpu)
                        setup.threads = std::max(1u, workers / marching);
                    spectra[d] = solve_canyon(setup);
                }
            });
            for (const auto& a : spectra)
                worst_residual = std::max(worst_residual, a.residual_db);
            const std::size_t cells = axes.azimuths.size() * axes.elevations.size();
            for (std::size_t k = 0; k < jobs.size(); k++)
            {
                const WaveJob& job = jobs[k];
                const auto wave = canyon_wave_gains(spectra[solve_of[k]], job.axis, axes.azimuths, axes.elevations,
                        bands, wave_bands, 5.0, tool_thread_count(threads));
                std::vector<double> slice(gains.begin() + static_cast<std::ptrdiff_t>(job.receiver * cells * nb),
                        gains.begin() + static_cast<std::ptrdiff_t>((job.receiver + 1) * cells * nb));
                blend_wave_band(slice, receivers[job.receiver].reflection, wave, weight);
                std::copy(slice.begin(), slice.end(),
                        gains.begin() + static_cast<std::ptrdiff_t>(job.receiver * cells * nb));
            }
            std::cerr << "Solved " << distinct.size() << " canyons by the wave equation to " << top << " Hz in "
                      << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
                      << " s; the slowest to decay ended " << -worst_residual << " dB down." << std::endl;
            if (worst_residual > -25.0)
                std::cerr << "niceops_urban_table: warning: a canyon's field had fallen only " << -worst_residual
                          << " dB after 4 s, so its wave band is short of its tail; rigid walls ring longest -- "
                             "give them a little absorption."
                          << std::endl;
        }
        std::map<std::string, std::string> provenance {
                {"canyons", canyons_file},
                {"diffraction", no_diffraction ? "none; specular rays in energy" : "UTD at the roof edges, coherent"},
                {"sound_speed_mps", exact_text(sound_speed)},
                {"facade_length_m", exact_text(facade_length * to_meters)},
                {"floor_absorption", exact_text(floor_absorption)},
        };
        if (!jobs.empty())
        {
            provenance["wave_below_hz"] = exact_text(wave_below);
            provenance["wave_blend_bands"] = std::to_string(wave_blend);
            provenance["wave_points_per_wavelength"] = exact_text(wave_ppw);
            provenance["wave_receivers"] = std::to_string(jobs.size());
        }
        const UrbanTable table(axes, std::move(receivers), std::move(gains), frame_bearing,
                jobs.empty() ? "ideal_canyon" : "ideal_canyon+wave", provenance);
        table.save(output);
        std::cerr << "Wrote " << output << " (" << table.receivers().size() << " receivers x "
                  << axes.azimuths.size() << " azimuths x " << axes.elevations.size() << " elevations x " << nb
                  << " bands)." << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "niceops_urban_table: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
