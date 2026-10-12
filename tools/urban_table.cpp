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
// The table is the canyon's under plane-wave incidence (canyon_ratio): the
// specular paths beyond the first-order ones and, where a wall blocks the
// direct path, the diffracted path over its top; the run adds the first-order
// reflections itself, from the two walls written beside the table, each
// --facade_length long.  The file is SI and records --frame_bearing, the
// compass bearing of the local +x axis, which the run checks against its own.

#include <CLI/CLI.hpp>

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "propagation_table.hpp"
#include "sphere_database.hpp"
#include "urban.hpp"
#include "util/constants.hpp"
#include "util/csv.hpp"

using namespace niceops; // NOLINT

int main(int argc, char** argv)
{
    CLI::App app {"Build an urban table for niceops --urban_table."};
    std::string canyons_file;
    std::string database_file;
    std::string output;
    std::string units = "ft";
    std::vector<double> bands;
    std::vector<double> elevations;
    int azimuths = 72;
    double frame_bearing = 90.0;
    double absorption = 0.0;
    double facade_length = 10000.0;
    double sound_speed = 343.2;
    bool no_diffraction = false;
    app.add_option("--canyons", canyons_file, "Receivers in ideal canyons, one per row (see the header).")->required();
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
    app.add_option("--facade_length", facade_length,
            "Length of each wall written for the run's first-order images, in --units (default 10 km).")
            ->check(CLI::PositiveNumber);
    app.add_option("--sound_speed", sound_speed, "For the diffracted path's wavelengths, m/s (default 343.2).")
            ->check(CLI::PositiveNumber);
    app.add_flag("--no_diffraction", no_diffraction, "Leave out the diffracted path over a blocking wall.");
    app.add_option("--units", units, "Length unit: ft (default) or m.")->check(CLI::IsMember({"ft", "m"}));
    app.add_option("-o,--output", output, "The table file to write.")->required();
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

        CsvReport report;
        const LabeledCsv csv = read_labeled_csv(canyons_file, &report);
        if (report.dropped > 0)
            throw std::runtime_error(canyons_file + ": line " + std::to_string(report.first_dropped_line)
                    + " is not a row of numbers");
        if (csv.rows.empty())
            throw std::runtime_error(canyons_file + " holds no receivers");
        // A length column may name its unit, as an observer file's do (x_ft,
        // height_m, "width (ft)"): it must then be --units.
        auto column = [&](const std::string& name, bool required, bool length) {
            const std::string ft = name + "_ft", m = name + "_m";
            const std::size_t c = length ? csv.column({name.c_str(), ft.c_str(), m.c_str()}) : csv.column({name.c_str()});
            if (c == LabeledCsv::npos)
            {
                if (required)
                    throw std::runtime_error(canyons_file + " has no '" + name + "' column");
                return c;
            }
            std::string declared = csv.unit_token(c);
            if (declared.empty() && csv.column_matches(c, {ft.c_str()}))
                declared = "ft";
            if (declared.empty() && csv.column_matches(c, {m.c_str()}))
                declared = "m";
            if (length && !declared.empty() && declared != units)
                throw std::runtime_error(canyons_file + ": column '" + csv.header[c] + "' is in " + declared
                        + ", and --units is " + units);
            return c;
        };
        const std::size_t cx = column("x", true, true), cy = column("y", true, true), cz = column("z", true, true);
        const std::size_t ch = column("height", true, true), axis_c = column("axis", true, false);
        const std::size_t cw = column("width", true, true), co = column("offset", true, true);
        const std::size_t cha = column("height_a", true, true), chb = column("height_b", true, true);
        const std::size_t caa = column("absorption_a", false, false), cab = column("absorption_b", false, false);

        const std::size_t nb = bands.size();
        std::vector<UrbanTable::Receiver> receivers;
        std::vector<double> gains;
        for (std::size_t row = 0; row < csv.rows.size(); row++)
        {
            const auto& r = csv.rows[row];
            if (r.size() != csv.header.size())
                throw std::runtime_error(canyons_file + ": row " + std::to_string(row + 1) + " has "
                        + std::to_string(r.size()) + " fields where the header names "
                        + std::to_string(csv.header.size()));
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
            const auto row_gains = ideal_canyon_gains(canyon, axis, axes, reflect_a, reflect_b, sound_speed,
                    !no_diffraction);
            gains.insert(gains.end(), row_gains.begin(), row_gains.end());
            receivers.push_back(std::move(receiver));
        }
        std::map<std::string, std::string> provenance {
                {"canyons", canyons_file},
                {"diffraction", no_diffraction ? "none" : "cnossos_dif placeholder over the blocking wall"},
                {"sound_speed_mps", exact_text(sound_speed)},
                {"facade_length_m", exact_text(facade_length * to_meters)},
        };
        const UrbanTable table(axes, std::move(receivers), std::move(gains), frame_bearing, "ideal_canyon",
                provenance);
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
