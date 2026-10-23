// The parabolic equation's whole field on one vertical plane, for plotting.
//
//   niceops_pe_field --atmosphere profile.csv --toward 312 --source_height 164
//                    --bands 250 500 1000 --max_range 9842 --top 400 > field.csv
//
// For a source in the vertical plane heading along a compass bearing, writes the
// band-averaged level relative to a point source's free field at the straight-
// line distance, on a grid of ranges (--step) and heights (--dz, from the ground
// to --top), as CSV rows band,range,height,level in --units.  The ground,
// turbulence, grid and thread options are niceops_pe_compare's (pe_options.hpp).  The heights are the
// PE's receivers, so asking for a tall grid raises the PE's own domain with it.

#include <CLI/CLI.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "atmosphere_profile.hpp"
#include "ground_effect.hpp"
#include "pe.hpp"
#include "stratified.hpp"
#include "util/constants.hpp"
#include "util/logger.hpp"
#include "util/parallel.hpp"

#include "pe_options.hpp"

using namespace niceops; // NOLINT

// A refused input -- an unreadable profile, a bad grid -- ends the run with its
// reason, not an abort.
int main(int argc, char** argv)
try
{
    CLI::App app {"The parabolic equation's field on one vertical plane."};
    std::string atmosphere_file;
    std::string units = "ft";
    double toward = 0.0;
    double source_height = 0.0;
    std::vector<double> bands = {250.0, 500.0, 1000.0};
    double max_range = 0.0;
    double step = 0.0;
    double top = 0.0;
    double dz_out = 0.0;
    PeOptions pe_options;
    app.add_option("--atmosphere", atmosphere_file, "Atmosphere profile, as for niceops.")->required();
    app.add_option("--units", units, "Length unit: ft (default) or m.")->check(CLI::IsMember({"ft", "m"}));
    app.add_option("--toward", toward, "Compass bearing the sound travels along, degrees.")->required();
    app.add_option("--source_height", source_height, "Source height, in --units.")->required()->check(CLI::NonNegativeNumber);
    app.add_option("--bands", bands, "One-third-octave band centers, Hz (default 250 500 1000).");
    app.add_option("--max_range", max_range, "Farthest range, in --units.")->required()->check(CLI::PositiveNumber);
    app.add_option("--step", step, "Range spacing of the output, in --units (default max_range / 200).")->check(CLI::NonNegativeNumber);
    app.add_option("--top", top, "Highest height written, in --units (default twice the source's).")->check(CLI::NonNegativeNumber);
    app.add_option("--dz", dz_out, "Height spacing of the output, in --units (default top / 150).")->check(CLI::NonNegativeNumber);
    pe_options.add_to(app);
    CLI11_PARSE(app, argc, argv);

    Logger::instance().set_file("log/pe_field.log");
    const unsigned workers = tool_thread_count(pe_options.threads);
    const LengthUnit unit = units == "m" ? meters : feet;
    const double z_s = source_height * unit.to_meters;
    const double r_max = max_range * unit.to_meters;
    const double dr_out = step > 0.0 ? step * unit.to_meters : r_max / 200.0;
    const double z_top = top > 0.0 ? top * unit.to_meters : 2.0 * z_s;
    const double dz = dz_out > 0.0 ? dz_out * unit.to_meters : z_top / 150.0;
    if (!(z_top > 0.0))
    {
        // A source on the ground has no 'twice its height' to default the top to;
        // the height loop below then never ended (dz = 0) and grew to 21 GB in 4 s.
        std::cerr << "niceops_pe_field: --source_height 0 gives no default for --top; "
                     "give --top, the highest height to write.\n";
        return 2;
    }

    auto profile = read_atmosphere_profile(atmosphere_file, unit, 0.0);
    if (!profile.has_pressure())
        profile = profile.with_hydrostatic_pressure(reference_pressure);
    const double sound_speed = profile.at(0.0).sound_speed();
    const double bearing = toward * pi / 180.0;
    const LayeredSoundSpeed layered(profile, 90.0);
    std::vector<double> c(layered.levels());
    layered.along(std::sin(bearing), std::cos(bearing), c.data());

    const auto admittance = ground_impedance(pe_options.ground(), bands, sound_speed);

    std::vector<double> heights;
    for (std::size_t j = 0; static_cast<double>(j) * dz <= z_top + 1e-9; j++)
        heights.push_back(static_cast<double>(j) * dz);
    std::vector<PeSetup> setups(bands.size());
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        PeSetup& setup = setups[b] = pe_options.setup(b, workers);
        setup.source_height = z_s;
        setup.receiver_heights = heights;
        setup.max_range = r_max;
        setup.output_step = dr_out;
        setup.admittance = admittance[b];
        setup.z = layered.heights();
        setup.c = c;
    }
    std::cerr << "PE field: " << bands.size() << " bands x " << pe_options.samples << " frequencies, "
              << heights.size() << " heights on " << workers << (workers == 1 ? " thread..." : " threads...")
              << std::endl;
    const auto pe = solve_pe_bands(setups, bands, pe_options.samples);

    std::cout << "band,range_" << unit.name << ",height_" << unit.name << ",level_db\n" << std::setprecision(6);
    for (std::size_t b = 0; b < bands.size(); b++)
        for (std::size_t q = 0; q < pe[b].range.size(); q++)
            for (std::size_t k = 0; k < heights.size(); k++)
                std::cout << bands[b] << "," << pe[b].range[q] / unit.to_meters << "," << heights[k] / unit.to_meters
                          << "," << pe[b].relative_level[k][q] << "\n";
    return 0;
}
catch (const std::exception& e)
{
    std::cerr << "niceops_pe_field: " << e.what() << std::endl;
    return 1;
}
