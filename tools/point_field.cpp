// The field of a point source among LoD1 buildings, at microphones
// (src/point_paths.hpp): each band's level re the source's own at 1 m, the
// rays summed coherently and band-averaged.
//
//   niceops_point_field --buildings district.csv --source 12.5 30 2 --microphones mics.csv
//       --bands 160 200 250 --humidity 70 --units m -o levels.csv
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
// 1/r, every ray absorbed by the air along its length (--humidity,
// --temperature, --pressure; --no_air_absorption leaves it out).  With --sources each row starts with the source's row in that file
// (from 1), the sources one after another, each over every microphone.
// Open floor at the same points is the floor's ray and the direct ray's; a
// level re free field at the microphone is the output less -20 lg r.
//
// Walls that scatter (--scattering, or a building's own scattering column)
// reflect specularly only rho (1 - s) of what reaches them, and the rest of
// rho diffusely, by Lambert's law: the urban tables' radiosity
// (diffuse_facades.hpp), lit by the source's specular paths to each patch and
// gathered at each microphone, added to the rays' band energy.  The roofs,
// which the rays reflect in, scatter with their building's coefficients.
// Patches are --patch_size at most, over the buildings with a wall within
// --diffuse_radius of a microphone, walls and roofs.  --scattering_spectrum gives s band by band (as
// tools/profile_scattering.py writes it), and a building's coefficient then
// weights it.
//
// --rays writes every ray too, for drawing (after the source's row, with
// --sources): the microphone's row, the ray's kind, reflections and length,
// its level at each band center (20 lg of its amplitude, no band averaging),
// and its path from the source as "x y z;x y z;..." (a leg off the floor
// unfolded below it).

#include <CLI/CLI.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "atmosphere.hpp"
#include "building_input.hpp"
#include "buildings.hpp"
#include "diffuse_facades.hpp"
#include "gpu/band_energy.hpp"
#include "gpu/band_stream.hpp"
#include "gpu/point_rays.hpp"
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
    double scattering = 0.0;
    std::string scattering_spectrum_file;
    double patch_size = 4.0;
    double diffuse_radius = 0.0;
    double floor_absorption = 0.0;
    double orientation_tolerance = 1.0;
    double sound_speed = 343.2;
    double humidity = std::numeric_limits<double>::quiet_NaN();
    double temperature = reference_temperature;
    double pressure = reference_pressure;
    bool no_air_absorption = false;
    bool fade = false;
    bool no_prune = false;
    bool edge_source = false;
    bool fade_incident = false;
    double fade_reach_given = 0.0;
    double radius = 150.0;
    std::size_t max_reflections = 4;
    std::size_t max_edge_reflections = 2;
    bool no_diffraction = false;
    bool no_double_diffraction = false;
    bool no_roofs = false;
    bool exact_shared_faces = false;
    bool cascade = false;
    std::string rays_file;
    int threads = 0;
    bool use_gpu = true;
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
    app.add_option("--scattering", scattering, "Walls' and roofs' scattering coefficient where a building gives none "
                                               "(default 0): the share of what they reflect that they scatter "
                                               "diffusely.");
    app.add_option("--scattering_spectrum", scattering_spectrum_file,
            "A CSV of band_hz and s: the scattering coefficient band by band, times each building's scattering "
            "column or --scattering (default 1 with this), interpolated in log frequency.");
    app.add_option("--patch_size", patch_size, "With scattering walls, the diffuse field's patches' size at most, in "
                                               "--units (default 4 m).");
    app.add_option("--diffuse_radius", diffuse_radius, "With scattering walls, a building with a wall within this "
                                                       "horizontal distance of a microphone scatters, walls and "
                                                       "roof, in --units (default --radius).");
    app.add_option("--floor_absorption", floor_absorption, "The floor's energy absorption coefficient (default 0).");
    app.add_option("--max_reflections", max_reflections,
            "Facade reflections on a ray that diffracts nowhere (default 4).");
    app.add_option("--max_edge_reflections", max_edge_reflections,
            "Facade reflections on each side of a single diffraction (default 2).");
    app.add_option("--radius", radius, "Facades and edges within this horizontal distance of the source or a "
                                       "microphone, in --units (default 150 m).");
    app.add_flag("--no_diffraction", no_diffraction, "Specular rays only.");
    app.add_flag("--no_double_diffraction", no_double_diffraction, "No ray diffracted at two edges.");
    app.add_flag("--no_roofs", no_roofs, "Roofs do not reflect, for comparison.");
    app.add_flag("--no_prune", no_prune,
            "Solve and check every diffracted path rather than ruling most out first: the same rays, slower; for "
            "checking the pruning.");
    app.add_flag("--cascade", cascade, "Over a shared face, the cascade of two UTD coefficients instead of the "
                                       "uniform form (double_wedge.hpp), for comparison.");
    app.add_flag("--exact_shared_faces", exact_shared_faces,
            "Over a shared face, the uniform form's integral at every wavenumber, not only where the cascade cannot "
            "stand in for it: the reference; much slower.");
    app.add_flag("--edge_source", edge_source,
            "Diffraction by the edge-source integral over each edge's whole length (edge_source.hpp), which keeps "
            "the edges' ends that UTD leaves out: every once-diffracted ray, and the twice-diffracted ones near a "
            "transition, an end or a partner's boundary. CPU only; not with --gpu or --fade.");
    app.add_option("--origin", origin, "With GeoJSON --buildings: the WGS84 origin, latitude longitude in degrees.")
            ->expected(2);
    app.add_option("--frame_bearing", frame_bearing, "With GeoJSON --buildings: compass bearing of +x (default 90).");
    app.add_option("--level_height", level_height, "With GeoJSON --buildings: metres per storey (default 3).");
    app.add_option("--default_height", default_height, "With GeoJSON --buildings: the height of a building that "
                                                        "gives none, m.");
    app.add_option("--orientation_tolerance", orientation_tolerance,
            "Walls within this many degrees of one orientation reflect alike (default 1).");
    app.add_option("--sound_speed", sound_speed, "m/s (default 343.2).");
    auto* humidity_option = app.add_option("--humidity", humidity,
            "Relative humidity, percent, for the air's absorption along every ray (ISO 9613-1). Required: the "
            "absorption depends strongly on it, so it is not assumed. Not with --no_air_absorption.");
    auto* temperature_option = app.add_option("--temperature", temperature,
            "Air temperature for the absorption, K (default 293.15).");
    auto* pressure_option = app.add_option("--pressure", pressure,
            "Ambient pressure for the absorption, kPa (default 101.325).");
    app.add_flag("--no_air_absorption", no_air_absorption, "Leave the air's absorption out, instead of --humidity.")
            ->excludes(humidity_option)
            ->excludes(temperature_option)
            ->excludes(pressure_option);
    app.add_flag("--fade", fade,
            "Make the diffracted ray the orders traced leave out at a boundary from the ray it partners, so that "
            "the field does not step there (specular and once-diffracted rays).");
    app.add_option("--fade_reach", fade_reach_given,
            "With --fade, the largest detour by an edge a partner is made for, in --units (default: where the "
            "lowest band's Fresnel parameter reaches 3, at most 0.5 m).");
    app.add_flag("--fade_incident", fade_incident,
            "With --fade, the made partner's coefficient with the incident wave's terms only, for comparison.");
    app.add_option("--rays", rays_file, "Also write every ray, with its level per band and its path.");
    app.add_flag("--gpu,!--no-gpu", use_gpu,
            "Form the rays' amplitudes and sum each band on the GPU (Metal on macOS) while the CPU traces the "
            "next pairs' rays: on by default wherever there is one, except with --edge_source or --fade, which "
            "are the CPU's. --gpu asks for it and is refused where there is none; --no-gpu keeps the CPU.");
    app.add_option("-j,--threads", threads, "Pairs of source and microphone evaluated at once. Default the environment's "
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
        // The air's absorption per metre of path: per sub-band for the levels,
        // at each band's centre for --rays.  Refused rather than assumed.
        std::vector<double> air_absorption, band_air(bands.size(), 0.0);
        if (!no_air_absorption)
        {
            if (std::isnan(humidity))
                throw std::runtime_error("the air's absorption needs --humidity (percent; it depends strongly on "
                                         "it), or --no_air_absorption to leave it out");
            for (const std::string& problem :
                    {temperature_problem(temperature), pressure_problem(pressure), humidity_problem(humidity)})
                if (!problem.empty())
                    throw std::runtime_error(problem);
            Atmosphere air;
            air.temperature = temperature;
            air.pressure = pressure;
            air.relative_humidity = humidity;
            air_absorption = sub_band_absorption(air, bands, sound_speed);
            const std::vector<double> at_centres = absorption_coefficients(air, bands);
            for (std::size_t b = 0; b < bands.size(); b++)
                band_air[b] = at_centres[b] * std::log(10.0) / 20.0;
            std::cerr << "Air: " << temperature << " K, " << pressure << " kPa, " << humidity << "% RH; absorption "
                      << at_centres.front() * 1000.0 << " dB/km at " << bands.front() << " Hz to "
                      << at_centres.back() * 1000.0 << " dB/km at " << bands.back() << " Hz." << std::endl;
        }
        // Defaults are lengths in metres; a value given is in --units.
        if (app.count("--radius") > 0)
            radius *= to_meters;
        if (app.count("--patch_size") > 0)
            patch_size *= to_meters;
        diffuse_radius = app.count("--diffuse_radius") > 0 ? diffuse_radius * to_meters : radius;

        BuildingInputOptions input;
        input.units = units;
        input.to_meters = to_meters;
        input.origin = origin;
        input.frame_bearing = frame_bearing;
        input.level_height = level_height;
        input.default_height = default_height;
        input.default_height_given = app.count("--default_height") > 0;
        input.absorption = absorption;
        // With a spectrum, a building's scattering coefficient weights it: 1, all
        // relief, by default.
        input.scattering = !scattering_spectrum_file.empty() && app.count("--scattering") == 0 ? 1.0 : scattering;
        if (!(input.scattering >= 0.0 && input.scattering <= 1.0))
            throw std::runtime_error("--scattering must lie within [0, 1]");
        std::vector<double> band_scattering(bands.size(), 1.0);
        if (!scattering_spectrum_file.empty())
        {
            const LabeledCsv rows = read_numeric_rows(scattering_spectrum_file, "scattering coefficients");
            const std::size_t fc = length_column(rows, scattering_spectrum_file, units, "band_hz", true, false);
            const std::size_t sc = length_column(rows, scattering_spectrum_file, units, "s", true, false);
            std::vector<std::pair<double, double>> spectrum;
            for (const auto& r : rows.rows)
            {
                if (!(r[fc] > 0.0) || !(r[sc] >= 0.0 && r[sc] <= 1.0))
                    throw std::runtime_error(scattering_spectrum_file + ": each band_hz must be positive and each s "
                                                                        "within [0, 1]");
                spectrum.emplace_back(r[fc], r[sc]);
            }
            std::sort(spectrum.begin(), spectrum.end());
            for (std::size_t i = 1; i < spectrum.size(); i++)
                if (spectrum[i].first == spectrum[i - 1].first)
                    throw std::runtime_error(scattering_spectrum_file + " gives " + exact_text(spectrum[i].first)
                            + " Hz twice");
            for (std::size_t b = 0; b < bands.size(); b++)
            {
                const double f = bands[b];
                if (f < spectrum.front().first / 1.01 || f > spectrum.back().first * 1.01)
                    throw std::runtime_error(scattering_spectrum_file + " covers " + exact_text(spectrum.front().first)
                            + " to " + exact_text(spectrum.back().first) + " Hz, not the band at " + exact_text(f)
                            + " Hz");
                std::size_t i = 1;
                while (i + 1 < spectrum.size() && spectrum[i].first < f)
                    i++;
                if (spectrum.size() == 1)
                    band_scattering[b] = spectrum[0].second;
                else
                {
                    const double t = std::clamp(std::log(f / spectrum[i - 1].first)
                                    / std::log(spectrum[i].first / spectrum[i - 1].first), 0.0, 1.0);
                    band_scattering[b] = spectrum[i - 1].second + t * (spectrum[i].second - spectrum[i - 1].second);
                }
            }
        }
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
        // A wall reflects rho = 1 - alpha of what reaches it and scatters s
        // of that: the rays carry rho (1 - s), and the diffuse field the rest.
        const std::size_t nb = bands.size();
        std::vector<double> rho(city.facades().size() * nb), scatter(rho.size()), reflection(rho.size());
        bool diffuse = false;
        for (std::size_t f = 0; f < city.facades().size(); f++)
            for (std::size_t b = 0; b < nb; b++)
            {
                const std::size_t building = city.facade_building(f);
                rho[f * nb + b] = 1.0 - district.absorption[building];
                scatter[f * nb + b] = district.scattering[building] * band_scattering[b];
                reflection[f * nb + b] = rho[f * nb + b] * (1.0 - scatter[f * nb + b]);
                diffuse = diffuse || scatter[f * nb + b] > 0.0;
            }

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
        options.roofs = !no_roofs;
        options.shared_face = !cascade;
        options.sample_shared_faces = !exact_shared_faces;
        options.prune = !no_prune;
        options.edge_source = edge_source;
        options.edge_source_reach = 100.0 / (2.0 * pi * *std::min_element(bands.begin(), bands.end()) / sound_speed);
        // The GPU wherever there is one and the run is one it takes; asked
        // for (--gpu), it must be.  The flag's count takes --no-gpu too, so
        // asking is its value as well.
        const bool gpu_asked = use_gpu && app.count("--gpu") > 0;
        if (edge_source && (gpu_asked || fade))
            throw std::runtime_error("--edge_source does not run with --gpu or --fade");
        if (fade)
        {
            const double k_min = 2.0 * pi * *std::min_element(bands.begin(), bands.end()) / sound_speed;
            options.fade_reach = app.count("--fade_reach") > 0 ? fade_reach_given * to_meters : fade_reach(k_min);
            options.fade_incident = fade_incident;
            if (gpu_asked)
                throw std::runtime_error("--fade does not run on the GPU yet");
        }
        else if (app.count("--fade_reach") > 0 || fade_incident)
            throw std::runtime_error("--fade_reach and --fade_incident go with --fade");
        options.radius = radius;

        // Every edge near some source, and near some microphone: what each
        // end's tree must hold for any pair it joins (PathEnd).
        const std::vector<std::size_t> near_sources = edges_near(city, sources, options.radius),
                                       near_microphones = edges_near(city, points, options.radius);
        const unsigned workers = tool_thread_count(threads);
        const std::size_t ns = sources.size(), nm = points.size();
        const double from_meters = 1.0 / to_meters;

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
        // A ray of --rays: its row of the file.
        auto ray_text = [&](const PointSourceField& field, std::size_t source_row, std::size_t row,
                                const std::vector<double>& facade_reflection) {
            // Each band's reflections, [band][facade].
            std::vector<std::vector<double>> band_reflection(nb);
            for (std::size_t b = 0; b < nb; b++)
                for (std::size_t f = 0; f < city.facades().size(); f++)
                    band_reflection[b].push_back(facade_reflection[f * nb + b]);
            std::ostringstream text;
            text << std::setprecision(8);
            for (std::size_t i = 0; i < field.rays().size(); i++)
            {
                const PointRay& ray = field.rays()[i];
                const std::string kind = std::string(ray.wedges == 0 ? "specular"
                                                             : ray.wedges == 1 ? "diffracted"
                                                             : ray.shared      ? "diffracted twice over a shared face"
                                                                               : "diffracted twice")
                        + (ray.made_edges.empty() ? "" : " (made partner)");
                if (!sources_file.empty())
                    text << source_row + 1 << ",";
                text << row + 1 << "," << i + 1 << "," << kind << "," << ray.facades.size() << "," << ray.floor << ","
                     << ray.length * from_meters;
                // The edge-source factors at the band centres.
                std::vector<std::complex<double>> factors;
                if (options.edge_source)
                {
                    std::vector<double> centres(nb);
                    for (std::size_t b = 0; b < nb; b++)
                        centres[b] = 2.0 * pi * bands[b] / sound_speed;
                    factors = edge_source_factors(city, ray, centres);
                }
                for (std::size_t b = 0; b < nb; b++)
                {
                    const double kb = 2.0 * pi * bands[b] / sound_speed;
                    const double a = (factors.empty()
                                             ? std::abs(point_ray_amplitude(ray, band_reflection[b],
                                                       1.0 - floor_absorption, kb))
                                             : reflected_amplitude(ray, band_reflection[b], 1.0 - floor_absorption)
                                                     * std::abs(factors[b]))
                            * std::exp(-band_air[b] * ray.length);
                    text << "," << 20.0 * std::log10(std::max(a, 1e-15));
                }
                text << ",";
                for (std::size_t v = 0; v < ray.path.size(); v++)
                    text << (v ? ";" : "") << ray.path[v].x * from_meters << " " << ray.path[v].y * from_meters << " "
                         << ray.path[v].z * from_meters;
                text << "\n";
            }
            return text.str();
        };
        // --gpu: the device, two sets of buffers to alternate between, and
        // the pairs a batch holds -- enough for every core several times over
        // (a GPU batch is a few milliseconds against a pair's tens).
        gpu::Device* device = nullptr;
        const std::size_t gpu_batch = std::max<std::size_t>(64, 4 * std::size_t {workers});
        std::vector<std::array<SubBand, 3>> parts;
        for (double f : bands)
            parts.push_back(sub_bands(f, sound_speed));
        if (use_gpu && !edge_source && !fade)
            device = gpu_asked ? &gpu::Device::require("--gpu") : gpu::Device::get();
        if (device != nullptr)
            std::cerr << "GPU: " << device->name() << std::endl;
        // With scattering walls, each microphone's diffuse field (its patches'
        // weights, diffuse_facades.hpp), for every source to light.
        std::vector<std::unique_ptr<DiffuseFacades>> diffuse_fields(diffuse ? nm : 0);
        if (diffuse)
        {
            DiffuseOptions diffuse_options;
            diffuse_options.radius = diffuse_radius;
            diffuse_options.patch_size = patch_size;
            // The rays reflect in the roofs, so the roofs scatter too.
            diffuse_options.roofs = !no_roofs;
            if (!no_air_absorption)
                diffuse_options.air = band_air;
            for (std::size_t row = 0; row < nm; row++)
            {
                diffuse_fields[row] = std::make_unique<DiffuseFacades>(city, points[row], rho, scatter,
                        1.0 - floor_absorption, bands, sound_speed, diffuse_options, workers);
                diffuse_fields[row]->prepare_point_lighting(max_reflections, workers);
                if (!diffuse_fields[row]->converged())
                    std::cerr << "Microphone " << row + 1 << ": the diffuse field's solve did not converge." << std::endl;
            }
            std::cerr << "Diffuse field: " << diffuse_fields[0]->patches().size() << " patches at microphone 1"
                      << std::endl;
        }
        // Each microphone's walls, and the roofs that take their coefficients:
        // those its diffuse field covers give the rays rho (1 - s), and the
        // rest reflect their whole rho, since nothing returns what they would
        // scatter.
        std::vector<std::vector<double>> microphone_reflection(diffuse ? nm : 0);
        for (std::size_t row = 0; row < microphone_reflection.size(); row++)
        {
            microphone_reflection[row] = rho;
            for (std::size_t f = 0; f < city.facades().size(); f++)
                if (diffuse_fields[row]->scatters(f))
                    for (std::size_t b = 0; b < nb; b++)
                        microphone_reflection[row][f * nb + b] = reflection[f * nb + b];
        }
        auto reflection_for = [&](std::size_t row) -> const std::vector<double>& {
            return diffuse ? microphone_reflection[row] : reflection;
        };
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
            // Traces the pairs [begin, end) of the block, microphone by
            // microphone, and hands each field to `use` with its place; a
            // worker keeps the microphone's end while its pairs last.
            auto trace_pairs = [&](std::size_t begin, std::size_t end, const auto& use) {
                parallel_chunks(end - begin, workers, [&](std::size_t chunk_begin, std::size_t chunk_end) {
                    std::unique_ptr<PathEnd> microphone;
                    std::size_t microphone_row = nm;
                    for (std::size_t index = begin + chunk_begin; index < begin + chunk_end; index++)
                    {
                        const std::size_t row = index / nk, k = index % nk;
                        if (row != microphone_row)
                        {
                            microphone = std::make_unique<PathEnd>(city, points[row], false, options, near_sources);
                            microphone_row = row;
                        }
                        const std::size_t pair = k * nm + row;
                        const PointSourceField field(city, *source_ends[k], *microphone, options);
                        counts[pair] = field.rays().size();
                        if (!rays_file.empty())
                            ray_rows[pair] = ray_text(field, first + k, row, reflection_for(row));
                        use(index - begin, pair, field);
                    }
                });
            };
            if (!device)
                trace_pairs(0, nm * nk, [&](std::size_t, std::size_t pair, const PointSourceField& field) {
                    field.band_levels(reflection_for(pair % nm), 1.0 - floor_absorption, bands, sound_speed,
                            &levels[pair * nb], air_absorption);
                });
            else
            {
                // Workers trace each pair and pack its rays (the GPU forms
                // their amplitudes) without pause; the stream takes them a
                // batch at a time (gpu/band_stream.hpp) and writes their
                // levels.
                gpu::BandStream stream(*device, gpu_batch, [&](std::size_t pair, const double* energy) {
                    for (std::size_t band = 0; band < nb; band++)
                        levels[pair * nb + band] = band_level(bands[band], parts[band], &energy[3 * band]);
                });
                try
                {
                    trace_pairs(0, nm * nk, [&](std::size_t, std::size_t pair, const PointSourceField& field) {
                        gpu::PointRaySet set;
                        gpu::pack_point_rays(field, reflection_for(pair % nm), 1.0 - floor_absorption, bands,
                                sound_speed, set, air_absorption);
                        stream.push(pair, std::move(set));
                    });
                }
                catch (...)
                {
                    // The GPU's own failure, where it is the cause.
                    stream.finish();
                    throw;
                }
                stream.finish();
            }

            // The diffuse energy adds to the rays' in each band.
            if (diffuse)
                parallel_chunks(nk * nm, workers, [&](std::size_t begin, std::size_t end) {
                    std::vector<double> energy(nb);
                    for (std::size_t pair = begin; pair < end; pair++)
                    {
                        const std::size_t k = pair / nm, row = pair % nm;
                        diffuse_fields[row]->point_energy(sources[first + k], max_reflections, band_air, energy.data());
                        for (std::size_t b = 0; b < nb; b++)
                        {
                            double& level = levels[pair * nb + b];
                            level = 10.0 * std::log10(std::pow(10.0, 0.1 * level) + energy[b]);
                        }
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
