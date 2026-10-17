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
// --microphones has x, y and height columns, by name, and any others are
// carried through.  Lengths are in --units, the output's too.
//
// Each output row is the microphone's x, y and height, then one level per
// band, dB: 10 lg of |p|^2 band-averaged for a source whose free field is
// 1/r.  Open floor at the same points is the floor's ray and the direct
// ray's; a level re free field at the microphone is the output less
// -20 lg r.

#include <CLI/CLI.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "building_input.hpp"
#include "buildings.hpp"
#include "point_paths.hpp"
#include "sphere_database.hpp"
#include "util/constants.hpp"
#include "util/paths.hpp"

using namespace niceops; // NOLINT

int main(int argc, char** argv)
{
    CLI::App app {"The field of a point source among buildings, at microphones."};
    std::string buildings_file;
    std::string microphones_file;
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
    app.add_option("--buildings", buildings_file, "The district: CSV of footprint vertices, or GeoJSON.")->required();
    app.add_option("--source", source, "The source's x, y and height above the floor.")->required()->expected(3);
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

        const Vec3 s(source[0] * to_meters, source[1] * to_meters, source[2] * to_meters);
        if (!(s.z >= 0.0))
            throw std::runtime_error("--source must stand on or above the floor");
        if (city.inside(s.x, s.y) < city.size() && s.z < city.footprints()[city.inside(s.x, s.y)].top)
            throw std::runtime_error("--source stands inside building " + district.ids[city.inside(s.x, s.y)]);

        const LabeledCsv mics = read_numeric_rows(microphones_file, "microphones");
        const std::size_t mx = length_column(mics, microphones_file, units, "x", true, true);
        const std::size_t my = length_column(mics, microphones_file, units, "y", true, true);
        const std::size_t mh = length_column(mics, microphones_file, units, "height", true, true);

        PointPathOptions options;
        options.max_reflections = max_reflections;
        options.max_edge_reflections = max_edge_reflections;
        options.diffraction = !no_diffraction;
        options.double_diffraction = !no_double_diffraction;
        options.shared_face = !cascade;
        options.radius = radius;

        ensure_parent_directory(output);
        std::ofstream out(output);
        if (!out)
            throw std::runtime_error("cannot write " + output);
        out << std::setprecision(10) << "x,y,height";
        for (double f : bands)
            out << "," << f;
        out << "\n";
        std::vector<double> levels(nb);
        for (std::size_t row = 0; row < mics.rows.size(); row++)
        {
            const auto& r = mics.rows[row];
            const Vec3 m(r[mx] * to_meters, r[my] * to_meters, r[mh] * to_meters);
            const std::string where = microphones_file + " row " + std::to_string(row + 1) + ": ";
            if (!(m.z >= 0.0))
                throw std::runtime_error(where + "the microphone's height must not be negative");
            if (city.inside(m.x, m.y) < city.size() && m.z < city.footprints()[city.inside(m.x, m.y)].top)
                throw std::runtime_error(where + "the microphone stands inside building "
                        + district.ids[city.inside(m.x, m.y)]);
            const PointSourceField field(city, s, m, options);
            field.band_levels(reflection, 1.0 - floor_absorption, bands, sound_speed, levels.data());
            out << r[mx] << "," << r[my] << "," << r[mh];
            for (double level : levels)
                out << "," << level;
            out << "\n";
            std::cerr << "microphone " << row + 1 << ": " << field.rays().size() << " rays" << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "niceops_point_field: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
