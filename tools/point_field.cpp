// The field of a point source among LoD1 buildings, at microphones
// (src/point_paths.hpp): each band's level re the source's own at 1 m, the
// rays summed coherently and band-averaged.
//
//   niceops_point_field --buildings district.csv --source 12.5 30 2 --microphones mics.csv
//       --bands 160 200 250 --units m -o levels.csv
//
// --buildings is niceops_urban_table's: a CSV of footprint vertices (building,
// x, y, top, and optionally absorption) or GeoJSON with --origin and
// --frame_bearing (tools/building_input.hpp).  The buildings stand on a flat
// floor at z = 0.  --source is the source's x, y and height above the floor;
// or --sources a CSV of them, x, y and height by name, for a source that
// moves -- each source's image tree is then traced once for every
// microphone, and each microphone's once for a block of sources.
// --microphones has x, y and height columns, by name, and any others are
// carried through.  Lengths are in --units, the output's too.
//
// Each output row is the microphone's x, y and height, then one level per
// band, dB: 10 lg of |p|^2 band-averaged for a source whose free field is
// 1/r.  With --sources each row starts with the source's row in that file
// (from 1), the sources one after another, each over every microphone.
// Open floor at the same points is the floor's ray and the direct ray's; a
// level re free field at the microphone is the output less -20 lg r.
//
// --rays writes every ray too, for drawing (after the source's row, with
// --sources): the microphone's row, the ray's kind, reflections and length,
// its level at each band center (20 lg of its amplitude, no band averaging),
// and its path from the source as "x y z;x y z;..." (a leg off the floor
// unfolded below it).

#include <CLI/CLI.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "building_input.hpp"
#include "buildings.hpp"
#include "point_paths.hpp"
#include "sphere_database.hpp"
#include "util/constants.hpp"
#include "util/parallel.hpp"
#include "util/paths.hpp"

using namespace niceops; // NOLINT

int main(int argc, char** argv)
{
    CLI::App app {"The field of a point source among buildings, at microphones."};
    std::string buildings_file;
    std::string microphones_file;
    std::string sources_file;
    std::string output;
    std::string database_file;
    std::string units = "ft";
    std::vector<double> source;
    std::vector<double> bands;
    std::vector<double> origin;
    double frame_bearing = 90.0;
    double level_height = 3.0;
    double default_height = 0.0;
    double absorption = 0.0;
    double floor_absorption = 0.0;
    double orientation_tolerance = 1.0;
    double sound_speed = 343.2;
    double radius = 150.0;
    std::size_t max_reflections = 4;
    std::size_t max_edge_reflections = 2;
    bool no_diffraction = false;
    bool no_double_diffraction = false;
    bool cascade = false;
    std::string rays_file;
    int threads = 0;
    app.add_option("--buildings", buildings_file, "The district: CSV of footprint vertices, or GeoJSON.")->required();
    auto* one_source = app.add_option("--source", source, "The source's x, y and height above the floor.")->expected(3);
    auto* many_sources = app.add_option("--sources", sources_file,
            "Instead of --source, a CSV of sources: x, y and height columns. Each row is a source position.");
    one_source->excludes(many_sources);
    app.add_option("--microphones", microphones_file, "The microphones: x, y and height columns.")->required();
    app.add_option("--bands", bands, "Band centers, Hz.");
    app.add_option("--database", database_file, "Take the band centers from this database instead.");
    app.add_option("--absorption", absorption, "Walls' energy absorption coefficient where a building gives none "
                                               "(default 0).");
    app.add_option("--floor_absorption", floor_absorption, "The floor's energy absorption coefficient (default 0).");
    app.add_option("--max_reflections", max_reflections,
            "Facade reflections on a ray that diffracts nowhere (default 4).");
    app.add_option("--max_edge_reflections", max_edge_reflections,
            "Facade reflections on each side of a single diffraction (default 2).");
    app.add_option("--radius", radius, "Facades and edges within this horizontal distance of the source or a "
                                       "microphone, in --units (default 150 m).");
    app.add_flag("--no_diffraction", no_diffraction, "Specular rays only.");
    app.add_flag("--no_double_diffraction", no_double_diffraction, "No ray diffracted at two edges.");
    app.add_flag("--cascade", cascade, "Over a shared face, the cascade of two UTD coefficients instead of the "
                                       "uniform form (double_wedge.hpp), for comparison.");
    app.add_option("--origin", origin, "With GeoJSON --buildings: the WGS84 origin, latitude longitude in degrees.")
            ->expected(2);
    app.add_option("--frame_bearing", frame_bearing, "With GeoJSON --buildings: compass bearing of +x (default 90).");
    app.add_option("--level_height", level_height, "With GeoJSON --buildings: metres per storey (default 3).");
    app.add_option("--default_height", default_height, "With GeoJSON --buildings: the height of a building that "
                                                        "gives none, m.");
    app.add_option("--orientation_tolerance", orientation_tolerance,
            "Walls within this many degrees of one orientation reflect alike (default 1).");
    app.add_option("--sound_speed", sound_speed, "m/s (default 343.2).");
    app.add_option("--rays", rays_file, "Also write every ray, with its level per band and its path.");
    app.add_option("-j,--threads", threads, "Microphones evaluated at once. Default the environment's "
                                            "NICEOPS_THREADS, else one per hardware thread.");
    app.add_option("--units", units, "Length unit: ft (default) or m.")->check(CLI::IsMember({"ft", "m"}));
    app.add_option("-o,--output", output, "The levels file to write.")->required();
    CLI11_PARSE(app, argc, argv);

    try
    {
        const double to_meters = units == "m" ? 1.0 : ft2m;
        if (!bands.empty() == !database_file.empty())
            throw std::runtime_error("give either --database or --bands");
        if (!database_file.empty())
        {
            const Database database(database_file);
            if (database.band_centers().empty())
                throw std::runtime_error(database_file + " carries no spectra");
            bands = database.band_centers();
        }
        if (!(absorption >= 0.0 && absorption <= 1.0) || !(floor_absorption >= 0.0 && floor_absorption <= 1.0))
            throw std::runtime_error("absorption coefficients must lie within [0, 1]");
        // Defaults are lengths in metres; a value given is in --units.
        if (app.count("--radius") > 0)
            radius *= to_meters;

        BuildingInputOptions input;
        input.units = units;
        input.to_meters = to_meters;
        input.origin = origin;
        input.frame_bearing = frame_bearing;
        input.level_height = level_height;
        input.default_height = default_height;
        input.default_height_given = app.count("--default_height") > 0;
        input.absorption = absorption;
        const BuildingInput district = read_buildings(buildings_file, input);
        Buildings city;
        try
        {
            city = Buildings(district.footprints, std::max(orientation_tolerance, 1e-6) * pi / 180.0);
        }
        catch (const std::invalid_argument& e)
        {
            throw std::runtime_error(buildings_file + ": " + e.what() + " (counting buildings from 1 in file order)");
        }
        const std::size_t nb = bands.size();
        std::vector<double> reflection(city.facades().size() * nb);
        for (std::size_t f = 0; f < city.facades().size(); f++)
            for (std::size_t b = 0; b < nb; b++)
                reflection[f * nb + b] = 1.0 - district.absorption[city.facade_building(f)];

        // What is wrong with where a point stands, "" for nothing.
        auto misplaced = [&](const Vec3& p) -> std::string {
            if (!(p.z >= 0.0))
                return "stands below the floor";
            const std::size_t b = city.inside(p.x, p.y);
            if (b < city.size() && p.z < city.footprints()[b].top)
                return "stands inside building " + district.ids[b];
            return "";
        };
        // Reading x, y and height columns by name, in --units.
        auto points_of = [&](const LabeledCsv& rows, const std::string& file) {
            const std::size_t x = length_column(rows, file, units, "x", true, true);
            const std::size_t y = length_column(rows, file, units, "y", true, true);
            const std::size_t h = length_column(rows, file, units, "height", true, true);
            std::vector<Vec3> points;
            for (const auto& r : rows.rows)
                points.emplace_back(r[x] * to_meters, r[y] * to_meters, r[h] * to_meters);
            return points;
        };

        if (source.empty() && sources_file.empty())
            throw std::runtime_error("give --source or --sources");
        const std::vector<Vec3> sources = sources_file.empty()
                ? std::vector<Vec3> {Vec3(source[0] * to_meters, source[1] * to_meters, source[2] * to_meters)}
                : points_of(read_numeric_rows(sources_file, "sources"), sources_file);
        for (std::size_t k = 0; k < sources.size(); k++)
            if (const std::string why = misplaced(sources[k]); !why.empty())
                throw std::runtime_error(
                        (sources_file.empty() ? "--source" : sources_file + " row " + std::to_string(k + 1)) + " " + why);

        const LabeledCsv mics = read_numeric_rows(microphones_file, "microphones");
        const std::size_t mx = length_column(mics, microphones_file, units, "x", true, true);
        const std::size_t my = length_column(mics, microphones_file, units, "y", true, true);
        const std::size_t mh = length_column(mics, microphones_file, units, "height", true, true);
        const std::vector<Vec3> points = points_of(mics, microphones_file);
        for (std::size_t row = 0; row < points.size(); row++)
            if (const std::string why = misplaced(points[row]); !why.empty())
                throw std::runtime_error(microphones_file + " row " + std::to_string(row + 1) + ": the microphone " + why);

        PointPathOptions options;
        options.max_reflections = max_reflections;
        options.max_edge_reflections = max_edge_reflections;
        options.diffraction = !no_diffraction;
        options.double_diffraction = !no_double_diffraction;
        options.shared_face = !cascade;
        options.radius = radius;

        // Every edge near some source, and near some microphone: what each
        // end's tree must hold for any pair it joins (PathEnd).
        const std::vector<std::size_t> near_sources = edges_near(city, sources, options.radius),
                                       near_microphones = edges_near(city, points, options.radius);
        const unsigned workers = tool_thread_count(threads);
        const std::size_t ns = sources.size(), nm = points.size();
        const double from_meters = 1.0 / to_meters;
        // For --rays: each band's reflections, [band][facade].
        std::vector<std::vector<double>> band_reflection(rays_file.empty() ? 0 : nb);
        for (std::size_t b = 0; b < band_reflection.size(); b++)
            for (std::size_t f = 0; f < city.facades().size(); f++)
                band_reflection[b].push_back(reflection[f * nb + b]);

        ensure_parent_directory(output);
        std::ofstream out(output);
        if (!out)
            throw std::runtime_error("cannot write " + output);
        out << std::setprecision(10) << (sources_file.empty() ? "" : "source,") << "x,y,height";
        for (double f : bands)
            out << "," << f;
        out << "\n";
        std::ofstream rays;
        if (!rays_file.empty())
        {
            ensure_parent_directory(rays_file);
            rays.open(rays_file);
            if (!rays)
                throw std::runtime_error("cannot write " + rays_file);
            rays << (sources_file.empty() ? "" : "source,") << "microphone,ray,kind,reflections,floor,length";
            for (double f : bands)
                rays << "," << f;
            rays << ",path\n";
        }

        // The sources go a block at a time, so that a long track holds one
        // block's trees and rays rather than the whole flight's.  Within a
        // block the pairs run microphone by microphone, which keeps every
        // worker busy however few microphones there are, and a worker keeps
        // the microphone's end while its pairs last.
        const std::size_t block = std::max<std::size_t>(64, 4 * std::size_t {workers});
        std::size_t total = 0;
        std::vector<std::unique_ptr<PathEnd>> source_ends;
        std::vector<double> levels;
        std::vector<std::size_t> counts;
        std::vector<std::string> ray_rows;
        for (std::size_t first = 0; first < ns; first += block)
        {
            const std::size_t nk = std::min(block, ns - first);
            source_ends.clear();
            source_ends.resize(nk);
            parallel_chunks(nk, workers, [&](std::size_t begin, std::size_t end) {
                for (std::size_t k = begin; k < end; k++)
                    source_ends[k] = std::make_unique<PathEnd>(city, sources[first + k], true, options, near_microphones);
            });
            levels.assign(nk * nm * nb, 0.0);
            counts.assign(nk * nm, 0);
            ray_rows.assign(rays_file.empty() ? 0 : nk * nm, std::string());
            parallel_chunks(nm * nk, workers, [&](std::size_t begin, std::size_t end) {
                std::unique_ptr<PathEnd> microphone;
                std::size_t microphone_row = nm;
                for (std::size_t index = begin; index < end; index++)
                {
                    const std::size_t row = index / nk, k = index % nk;
                    if (row != microphone_row)
                    {
                        microphone = std::make_unique<PathEnd>(city, points[row], false, options, near_sources);
                        microphone_row = row;
                    }
                    const std::size_t pair = k * nm + row;
                    const PointSourceField field(city, *source_ends[k], *microphone, options);
                    field.band_levels(reflection, 1.0 - floor_absorption, bands, sound_speed, &levels[pair * nb]);
                    counts[pair] = field.rays().size();
                    if (rays_file.empty())
                        continue;
                    std::ostringstream text;
                    text << std::setprecision(8);
                    for (std::size_t i = 0; i < field.rays().size(); i++)
                    {
                        const PointRay& ray = field.rays()[i];
                        const char* kind = ray.wedges == 0 ? "specular"
                                : ray.wedges == 1           ? "diffracted"
                                : ray.shared                ? "diffracted twice over a shared face"
                                                            : "diffracted twice";
                        if (!sources_file.empty())
                            text << first + k + 1 << ",";
                        text << row + 1 << "," << i + 1 << "," << kind << "," << ray.facades.size() << ","
                             << ray.floor << "," << ray.length * from_meters;
                        for (std::size_t b = 0; b < nb; b++)
                        {
                            const double kb = 2.0 * pi * bands[b] / sound_speed;
                            const double a = std::abs(
                                    point_ray_amplitude(ray, band_reflection[b], 1.0 - floor_absorption, kb));
                            text << "," << 20.0 * std::log10(std::max(a, 1e-15));
                        }
                        text << ",";
                        for (std::size_t v = 0; v < ray.path.size(); v++)
                            text << (v ? ";" : "") << ray.path[v].x * from_meters << " "
                                 << ray.path[v].y * from_meters << " " << ray.path[v].z * from_meters;
                        text << "\n";
                    }
                    ray_rows[pair] = text.str();
                }
            });

            for (std::size_t k = 0; k < nk; k++)
                for (std::size_t row = 0; row < nm; row++)
                {
                    const auto& r = mics.rows[row];
                    if (!sources_file.empty())
                        out << first + k + 1 << ",";
                    out << r[mx] << "," << r[my] << "," << r[mh];
                    for (std::size_t b = 0; b < nb; b++)
                        out << "," << levels[(k * nm + row) * nb + b];
                    out << "\n";
                }
            for (std::size_t c : counts)
                total += c;
            for (const std::string& text : ray_rows)
                rays << text;
        }
        if (!out || (!rays_file.empty() && !rays))
            throw std::runtime_error("cannot write " + (!out ? output : rays_file));
        if (!sources_file.empty())
            std::cerr << ns << " sources, ";
        std::cerr << nm << " microphones, " << total << " rays" << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "niceops_point_field: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
