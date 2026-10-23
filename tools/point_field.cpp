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
#include "facade_relief.hpp"
#include "grating_paths.hpp"
#include "ground_effect.hpp"
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
    std::string roofs_file;
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
    std::vector<std::string> relief_files;
    std::string relief_scattering = "hybrid";
    std::size_t grating_reflections = 2;
    double patch_size = 4.0;
    double diffuse_radius = 0.0;
    double floor_absorption = 0.0;
    double flow_resistance = std::numeric_limits<double>::quiet_NaN();
    std::string impedance_model = "delany_bazley";
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
    std::size_t max_twice_reflections = 2;
    std::size_t max_between_reflections = 1;
    std::size_t max_thrice_reflections = 0;
    std::size_t max_chain = 8;
    std::string triple_diffraction = "roofs";
    bool no_diffraction = false;
    bool no_double_diffraction = false;
    bool no_roofs = false;
    bool no_slope_facades = false;
    bool no_twice_roofs = false;
    bool no_roof_lighting = false;
    bool no_profile_pe = false;
    bool tilted_step_split = false;
    bool unending_edges = false;
    bool exact_shared_faces = false;
    bool cascade = false;
    std::string rays_file;
    int threads = 0;
    bool use_gpu = true;
    app.add_option("--buildings", buildings_file, "The district: CSV of footprint vertices, or GeoJSON.")->required();
    app.add_option("--roofs", roofs_file, "Pitched roofs (LoD2) for a CSV district: one row per roof plane, building (its id), x, y and z (a "
            "point on the plane) and slope_x and slope_y (its rise along x and y). A building's roof is the "
            "least of its planes, capped at its top; one given none stays flat.");
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
    app.add_option("--relief", relief_files,
            "Reliefs' grating orders (tools/profile_scattering.py --orders), one or more: a building carries the one "
            "its relief column names (from 1, 0 none), or without the column the first wherever it scatters "
            "(--scattering or its column, 1 by default with this), on every wall, each wall's +along its base "
            "direction. Its specular rays take the relief's order 0 at the incidence they meet it at, and the rest goes "
            "by --relief_scattering.");
    app.add_option("--relief_scattering", relief_scattering,
            "With --relief: hybrid (default), the other orders as rays that take one each and what a wall scatters "
            "onto other walls diffusely; orders, the rays alone; or lambert, all diffusely by Lambert's law at the "
            "incidence each lighting path meets the wall at.")
            ->check(CLI::IsMember({"orders", "lambert", "hybrid"}));
    app.add_option("--grating_reflections", grating_reflections,
            "With --relief_scattering orders: facade reflections on a grating ray besides its order (default 2).");
    app.add_option("--patch_size", patch_size, "With scattering walls, the diffuse field's patches' size at most, in "
                                               "--units (default 4 m).");
    app.add_option("--diffuse_radius", diffuse_radius, "With scattering walls, a building with a wall within this "
                                                       "horizontal distance of a microphone scatters, walls and "
                                                       "roof, in --units (default --radius).");
    app.add_option("--floor_absorption", floor_absorption, "The floor's energy absorption coefficient (default 0).");
    app.add_option("--flow_resistance", flow_resistance,
            "A floor of finite impedance instead, its flow resistance in kPa s/m^2: each floor bounce takes the "
            "ground's spherical-wave reflection coefficient, as niceops --buildings --ground_effect does. Not with "
            "--floor_absorption.");
    app.add_option("--impedance_model", impedance_model,
            "With --flow_resistance: delany_bazley (default) or miki, as niceops takes them.")
            ->check(CLI::IsMember({"delany_bazley", "miki"}));
    app.add_option("--max_reflections", max_reflections,
            "Facade reflections on a ray that diffracts nowhere (default 4).");
    app.add_option("--max_edge_reflections", max_edge_reflections,
            "Facade reflections on each side of a single diffraction (default 2).");
    app.add_option("--max_twice_reflections", max_twice_reflections,
            "Facade reflections on each side of a double diffraction (default 2).");
    app.add_option("--max_between_reflections", max_between_reflections,
            "Reflections, in a facade or the floor, between the two edges of a double diffraction, 0 or 1 "
            "(default 1).")
            ->check(CLI::Range(0, 1));
    app.add_option("--radius", radius, "Facades and edges within this horizontal distance of the source or a "
                                       "microphone, in --units (default 150 m).");
    app.add_flag("--no_diffraction", no_diffraction, "Specular rays only.");
    app.add_flag("--no_double_diffraction", no_double_diffraction, "No ray diffracted at two edges or more.");
    app.add_option("--triple_diffraction", triple_diffraction,
            "Rays diffracted at three edges: 'roofs' (default; three roof edges, the rows a street's sound crosses), "
            "'all' (any three, every chain of corners round a block too: far more rays), or 'none'.")
            ->check(CLI::IsMember({"none", "roofs", "all"}));
    app.add_option("--max_thrice_reflections", max_thrice_reflections,
            "Facade reflections on each side of a triple diffraction (default 0).");
    app.add_option("--max_chain", max_chain,
            "Rays diffracted at more than three roof edges, up to this many (default 8, at most 14; 3 for none): "
            "the corners of the taut string over the buildings' cross-section between the source and the "
            "microphone, a mansard's or several rows of roofs from low beside them.")
            ->check(CLI::Range(3, 14));
    app.add_flag("--no_roofs", no_roofs, "Roofs do not reflect, for comparison.");
    app.add_flag("--no_slope_facades", no_slope_facades,
            "No ray off both a slope and a wall, for comparison: no specular ray off walls and a slope, and no "
            "once-diffracted ray with walls and a slope on one side.");
    app.add_flag("--no_roof_lighting", no_roof_lighting,
            "With scattering, no lighting path off a level roof or a slope onto the patches, for comparison.");
    app.add_flag("--no_profile_pe", no_profile_pe,
            "Over a row of flat-roofed blocks with a street, crossed near grazing, the cascades' rays rather than the "
            "row's profile marched by a parabolic equation, for comparison.");
    app.add_flag("--tilted_step_split", tilted_step_split,
            "Across the line through a ridge or roof edge and the edge of the lower roof it steps down to, split the "
            "field's jump either side rather than take it as the rays leave it, stepping where the microphone "
            "crosses the line (opt-in: fitted on one step).");
    app.add_flag("--no_twice_roofs", no_twice_roofs,
            "No level roof's or slope's reflection on a twice-diffracted ray, for comparison.");
    app.add_flag("--unending_edges", unending_edges,
            "UTD's edges without their ends, each taken as its whole infinite line, for comparison.");
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
        // The floor: an energy coefficient, or a ground's admittance at the
        // sub-bands' middles (and, for --rays, at the bands' centres), with
        // what it reflects of a diffuse field for the diffuse facades.
        FloorModel floor(1.0 - floor_absorption), centre_floor = floor;
        double floor_energy = 1.0 - floor_absorption;
        if (app.count("--flow_resistance") > 0)
        {
            if (app.count("--floor_absorption") > 0)
                throw std::runtime_error("--flow_resistance gives the floor an impedance, which --floor_absorption's "
                                         "energy coefficient would replace; give one");
            GroundModel ground;
            ground.flow_resistance = flow_resistance;
            ground.impedance_model = parse_impedance_model(impedance_model);
            floor.admittance = ground_impedance(ground, sub_band_middles(bands, sound_speed), sound_speed);
            centre_floor.admittance = ground_impedance(ground, bands, sound_speed);
            // The diffuse field's floor, each band's at its middle sub-band,
            // as a run takes it.
            floor_energy = 0.0;
            for (std::size_t b = 0; b < bands.size(); b++)
                floor_energy += diffuse_reflection(floor.admittance[3 * b + 1]) / static_cast<double>(bands.size());
            std::cerr << "Floor: " << impedance_model << " at " << flow_resistance << " kPa s/m^2." << std::endl;
        }
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
        input.roofs = roofs_file;
        // With a spectrum, a building's scattering coefficient weights it: 1, all
        // relief, by default.
        input.scattering = (!scattering_spectrum_file.empty() || !relief_files.empty()) && app.count("--scattering") == 0
                ? 1.0
                : scattering;
        if (!relief_files.empty() && !scattering_spectrum_file.empty())
            throw std::runtime_error("give --relief or --scattering_spectrum: each is the walls' scattering");
        std::vector<FacadeRelief> relief_tables;
        for (const std::string& file : relief_files)
            relief_tables.push_back(FacadeRelief::read(file));
        const bool relief = !relief_tables.empty();
        if (!(input.scattering >= 0.0 && input.scattering <= 1.0))
            throw std::runtime_error("--scattering must lie within [0, 1]");
        const std::vector<double> band_scattering = scattering_spectrum_file.empty()
                ? std::vector<double>(bands.size(), 1.0)
                : scattering_at(read_scattering_spectrum(scattering_spectrum_file, units), bands,
                          scattering_spectrum_file);
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
        // A relief's walls carry rho on the rays, which take its order 0 at
        // their incidence; by orders they scatter nothing diffusely, and by
        // Lambert's law their random-incidence s, scaled to each path's.
        const std::size_t nb = bands.size();
        std::vector<double> rho(city.facades().size() * nb), scatter(rho.size()), reflection(rho.size());
        std::vector<const FacadeRelief*> reliefs(city.facades().size(), nullptr);
        // Each relief table's random-incidence s per band, once.
        std::vector<std::vector<double>> random_s(relief_tables.size(), std::vector<double>(nb, 0.0));
        if (relief_scattering != "orders")
            for (std::size_t t = 0; t < relief_tables.size(); t++)
                for (std::size_t b = 0; b < nb; b++)
                    random_s[t][b] = relief_tables[t].random_incidence(bands[b], sound_speed);
        bool diffuse = false;
        for (std::size_t f = 0; f < city.facades().size(); f++)
            for (std::size_t b = 0; b < nb; b++)
            {
                const std::size_t building = city.facade_building(f);
                rho[f * nb + b] = 1.0 - district.absorption[building];
                const int which = relief ? relief_of(district, building, relief_tables.size()) : 0;
                if (which > 0)
                {
                    const FacadeRelief& r = relief_tables[static_cast<std::size_t>(which - 1)];
                    reliefs[f] = &r;
                    scatter[f * nb + b] = random_s[static_cast<std::size_t>(which - 1)][b];
                    reflection[f * nb + b] = rho[f * nb + b];
                }
                else
                {
                    scatter[f * nb + b] = relief ? 0.0 : district.scattering[building] * band_scattering[b];
                    reflection[f * nb + b] = rho[f * nb + b] * (1.0 - scatter[f * nb + b]);
                }
                diffuse = diffuse || scatter[f * nb + b] > 0.0;
            }
        if (relief)
        {
            std::size_t carrying = 0;
            for (const FacadeRelief* r : reliefs)
                carrying += r != nullptr;
            for (const FacadeRelief& r : relief_tables)
                std::cerr << "Relief: " << r.source() << ", period " << r.period() << " m." << std::endl;
            std::cerr << "Reliefs on " << carrying << " walls, scattering by " << relief_scattering << "." << std::endl;
        }

        // What is wrong with where a point stands, "" for nothing.
        auto misplaced = [&](const Vec3& p) { return misplaced_point(city, district.ids, p); };
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
        options.max_twice_reflections = max_twice_reflections;
        options.max_between_reflections = max_between_reflections;
        options.triple_diffraction = triple_diffraction == "none" ? PointPathOptions::Triple::none
                : triple_diffraction == "all"                      ? PointPathOptions::Triple::all
                                                                   : PointPathOptions::Triple::roofs;
        options.max_thrice_reflections = max_thrice_reflections;
        options.max_chain = max_chain;
        options.diffraction = !no_diffraction;
        options.double_diffraction = !no_double_diffraction;
        options.roofs = !no_roofs;
        options.slope_facades = !no_slope_facades;
        options.twice_roofs = !no_twice_roofs;
        options.profile_pe = !no_profile_pe;
        options.tilted_step_split = tilted_step_split;
        options.edge_ends = !unending_edges;
        options.shared_face = !cascade;
        options.sample_shared_faces = !exact_shared_faces;
        options.prune = !no_prune;
        options.edge_source = edge_source;
        if (relief)
        {
            options.reliefs = &reliefs;
            // Grating rays' source side reaches as deep as their own orders.
            if (relief_scattering != "lambert")
                options.grating_reflections = grating_reflections;
        }
        const double lowest_k = 2.0 * pi * *std::min_element(bands.begin(), bands.end()) / sound_speed;
        options.edge_source_reach = 100.0 / lowest_k;
        // UTD's end waves fade out by end_fade_hi radians, at the lowest band last.
        options.end_reach = end_fade_hi / lowest_k;
        // The split's foot fields are taken where its fade at the lowest
        // sub-band's lower edge is above its floor.
        options.tilted_step_k_min = lowest_k / std::pow(2.0, 1.0 / 6.0);
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
        // A run of fewer pairs than threads shares each pair's integrals out.
        if (edge_source && ns * nm < workers)
            options.edge_source_threads = static_cast<unsigned>((workers + ns * nm - 1) / (ns * nm));
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
                                const std::vector<double>& facade_reflection, const std::vector<GratingRay>& gratings) {
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
                const ProfileTerm& part = ray.profile_term;
                const std::string kind = ray.profile
                        ? std::string("profile field") + (part.source_image ? ", off the floor at the source" : "")
                                + (part.microphone_image ? ", off the floor at the microphone" : "")
                        : std::string(ray.wedges == 0 ? "specular"
                                                             : ray.wedges == 1 ? "diffracted"
                                                             : ray.wedges == 3 ? "diffracted three times"
                                                             : ray.wedges > 3  ? "diffracted " + std::to_string(ray.wedges)
                                                                           + " times"
                                                             : ray.shared && ray.wedge[1].grazing
                                                                     ? "diffracted twice over a shared face"
                                                             : ray.shared ? "diffracted twice uniformly over the pair"
                                                                          : "diffracted twice")
                        + (ray.made_edges.empty() ? "" : " (made partner)");
                if (!sources_file.empty())
                    text << source_row + 1 << ",";
                text << row + 1 << "," << i + 1 << "," << kind << "," << ray.facades.size() << "," << ray.floor << ","
                     << ray.length * from_meters;
                // The floor's factor and the edge-source factors at the band
                // centres.
                std::vector<double> centres(nb);
                for (std::size_t b = 0; b < nb; b++)
                    centres[b] = 2.0 * pi * bands[b] / sound_speed;
                std::vector<std::complex<double>> at_floor(nb);
                floor_factors(ray, centre_floor, centres, nullptr, at_floor.data());
                // A relief's specular share where the ray takes one.
                if (meets_relief(ray, field.reliefs()))
                    relief_factors(ray, *field.reliefs(), bands, sound_speed, at_floor.data());
                std::vector<std::complex<double>> factors;
                if (options.edge_source)
                    factors = edge_source_factors(city, ray, centres);
                for (std::size_t b = 0; b < nb; b++)
                {
                    const double kb = 2.0 * pi * bands[b] / sound_speed;
                    const double a = (factors.empty()
                                             ? std::abs(point_ray_amplitude(ray, band_reflection[b],
                                                       centre_floor.reflection, kb))
                                             : reflected_amplitude(ray, band_reflection[b], centre_floor.reflection)
                                                     * std::abs(factors[b]))
                            * std::abs(at_floor[b]) * std::exp(-band_air[b] * ray.length);
                    text << "," << 20.0 * std::log10(std::max(a, 1e-15));
                }
                text << ",";
                for (std::size_t v = 0; v < ray.path.size(); v++)
                    text << (v ? ";" : "") << ray.path[v].x * from_meters << " " << ray.path[v].y * from_meters << " "
                         << ray.path[v].z * from_meters;
                text << "\n";
            }
            // The grating rays: each carries energy in its own band alone.
            for (std::size_t i = 0; i < gratings.size(); i++)
            {
                const GratingRay& ray = gratings[i];
                if (!sources_file.empty())
                    text << source_row + 1 << ",";
                text << row + 1 << "," << field.rays().size() + i + 1 << ",grating order " << ray.order << ","
                     << ray.path.size() - 2 << "," << (ray.source_floor ? 1 : 0) + (ray.microphone_floor ? 1 : 0)
                     << "," << ray.length * from_meters;
                for (std::size_t b = 0; b < nb; b++)
                    text << "," << (b == ray.band ? 10.0 * std::log10(std::max(ray.energy, 1e-30)) : -300.0);
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
            diffuse_options.roof_lighting = !no_roof_lighting;
            if (relief && relief_scattering != "orders")
                diffuse_options.scattering_at = [&](std::size_t f, double sin_psi, double cos_elevation, std::size_t b) {
                    return reliefs[f] != nullptr ? 1.0 - reliefs[f]->specular(sin_psi, bands[b] * cos_elevation, sound_speed)
                                                 : -1.0;
                };
            if (relief && relief_scattering == "hybrid")
                diffuse_options.first_scattering_elsewhere = [&](std::size_t f) { return reliefs[f] != nullptr; };
            if (!no_air_absorption)
                diffuse_options.air = band_air;
            for (std::size_t row = 0; row < nm; row++)
            {
                diffuse_fields[row] = std::make_unique<DiffuseFacades>(city, points[row], rho, scatter,
                        floor_energy, bands, sound_speed, diffuse_options, workers);
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
        // By Lambert's law alone, a relief beyond a microphone's diffuse reach
        // reflects whole, as a plain wall there does.
        std::vector<std::vector<const FacadeRelief*>> microphone_reliefs;
        std::vector<PointPathOptions> microphone_options(nm, options);
        if (relief && relief_scattering == "lambert" && diffuse)
        {
            microphone_reliefs.assign(nm, reliefs);
            for (std::size_t row = 0; row < nm; row++)
            {
                for (std::size_t f = 0; f < reliefs.size(); f++)
                    if (!diffuse_fields[row]->scatters(f))
                        microphone_reliefs[row][f] = nullptr;
                microphone_options[row].reliefs = &microphone_reliefs[row];
            }
        }
        std::size_t total = 0;
        // With a relief's orders as rays, each pair's energy from them per
        // band, and how many.
        const bool gratings = relief && relief_scattering != "lambert";
        GratingInput grating_input;
        grating_input.reliefs = &reliefs;
        grating_input.bands = bands;
        grating_input.sound_speed = sound_speed;
        if (!no_air_absorption)
            grating_input.air = band_air;
        grating_input.max_reflections = grating_reflections;
        // The hybrid's diffuse field carries what a relief scatters onto
        // other walls: its grating rays reach the microphone directly.
        if (relief_scattering == "hybrid")
            grating_input.max_microphone_reflections = 0;
        grating_input.source_reflections = options.source_reflections;
        if (centre_floor.impedance())
            grating_input.floor = [&](double sin_grazing, double leg, std::size_t b) {
                const auto q = spherical_reflection_coefficient(sin_grazing, pi * bands[b] * leg / sound_speed,
                        centre_floor.admittance[b]);
                return std::complex<double>(q[0], q[1]);
            };
        else
            grating_input.floor = [&](double, double, std::size_t) {
                return std::complex<double>(std::sqrt(centre_floor.reflection));
            };
        std::vector<double> grating;
        std::vector<std::size_t> grating_count;
        std::size_t gratings_total = 0;
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
            grating.assign(gratings ? nk * nm * nb : 0, 0.0);
            grating_count.assign(gratings ? nk * nm : 0, 0);
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
                        const PointSourceField field(city, *source_ends[k], *microphone, microphone_options[row]);
                        counts[pair] = field.rays().size();
                        std::vector<GratingRay> found;
                        if (gratings)
                        {
                            GratingInput g = grating_input;
                            g.reflection = &reflection_for(row);
                            grating_rays(city, *source_ends[k], *microphone, g, found);
                            grating_band_energies(found, bands, sound_speed, nullptr, &grating[pair * nb]);
                            grating_count[index - begin] = found.size();
                        }
                        if (!rays_file.empty())
                            ray_rows[pair] = ray_text(field, first + k, row, reflection_for(row), found);
                        use(index - begin, pair, field);
                    }
                });
            };
            if (!device)
                trace_pairs(0, nm * nk, [&](std::size_t, std::size_t pair, const PointSourceField& field) {
                    field.band_levels(reflection_for(pair % nm), floor, bands, sound_speed,
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
                        gpu::pack_point_rays(field, reflection_for(pair % nm), floor, bands,
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

            // The grating rays' energy adds to the rays' in each band.
            for (std::size_t i = 0; i < grating.size(); i++)
                if (grating[i] > 0.0)
                    levels[i] = 10.0 * std::log10(std::pow(10.0, 0.1 * levels[i]) + grating[i]);
            for (std::size_t c : grating_count)
                gratings_total += c;
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
        std::cerr << nm << " microphones, " << total << " rays";
        if (gratings)
            std::cerr << ", " << gratings_total << " grating rays";
        std::cerr << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "niceops_point_field: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
