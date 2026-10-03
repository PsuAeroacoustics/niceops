// The parabolic equation's whole field on one vertical plane, for plotting.
//
//   niceops_pe_field --atmosphere profile.csv --toward 312 --source_height 164
//                    --bands 250 500 1000 --max_range 9842 --top 400 > field.csv
//
// For a source in the vertical plane heading along a compass bearing, writes the
// band-averaged level relative to a point source's free field at the straight-
// line distance, on a grid of ranges (--step) and heights (--dz, from the ground
// to --top), as CSV rows band,range,height,level in --units.  The ground,
// turbulence and grid options are niceops_pe_compare's.  The heights are the
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

using namespace niceops; // NOLINT

int main(int argc, char** argv)
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
    double flow_resistance = flow_resistance::grassy_field;
    std::string impedance_model = "delany_bazley";
    int samples = 5;
    double points_per_wavelength = 0.0;
    double range_points_per_wavelength = 10.0;
    std::string scheme = "pade";
    double range_step_wavelengths = 8.0;
    double max_angle = 50.0;
    int vertical_order = 4;
    int boundary_order = 4;
    std::string starter = "self";
    double gamma_t = 0.0;
    double outer_scale = 10.0;
    unsigned seed = 1;
    app.add_option("--atmosphere", atmosphere_file, "Atmosphere profile, as for niceops.")->required();
    app.add_option("--units", units, "Length unit: ft (default) or m.")->check(CLI::IsMember({"ft", "m"}));
    app.add_option("--toward", toward, "Compass bearing the sound travels along, degrees.")->required();
    app.add_option("--source_height", source_height, "Source height, in --units.")->required()->check(CLI::NonNegativeNumber);
    app.add_option("--bands", bands, "One-third-octave band centres, Hz (default 250 500 1000).");
    app.add_option("--max_range", max_range, "Farthest range, in --units.")->required()->check(CLI::PositiveNumber);
    app.add_option("--step", step, "Range spacing of the output, in --units (default max_range / 200).")->check(CLI::NonNegativeNumber);
    app.add_option("--top", top, "Highest height written, in --units (default twice the source's).")->check(CLI::NonNegativeNumber);
    app.add_option("--dz", dz_out, "Height spacing of the output, in --units (default top / 150).")->check(CLI::NonNegativeNumber);
    app.add_option("--flow_resistance", flow_resistance, "Ground flow resistance, kPa s/m^2.");
    app.add_option("--impedance_model", impedance_model, "As for niceops.");
    app.add_option("--samples", samples, "PE frequencies per band (default 5).")->check(CLI::PositiveNumber);
    app.add_option("--ppw", points_per_wavelength, "PE grid points per wavelength in height (default 10 where the grid is fourth order "
            "throughout, 20 under --scheme cn or either order 2).")->check(CLI::NonNegativeNumber);
    app.add_option("--scheme", scheme,
            "How the PE marches in range: pade (default), the split-step Pade march, or cn, Crank-Nicolson.")
            ->check(CLI::IsMember({"pade", "cn"}));
    app.add_option("--range_step", range_step_wavelengths,
            "--scheme pade: longest range step, in wavelengths (default 8).")->check(CLI::PositiveNumber);
    app.add_option("--max_angle", max_angle,
            "--scheme pade: elevation angle the march is accurate to, degrees (default 50, at most 60).  Steeper "
            "paths from the source come out as NaN.")->check(CLI::Range(0.1, 60.0));
    app.add_option("--vertical_order", vertical_order,
            "--scheme pade: order of the height derivative on the grid, 4 (default) or 2.")
            ->check(CLI::IsMember({2, 4}));
    app.add_option("--boundary_order", boundary_order,
            "--scheme pade: order of the ground condition and of the source's and receivers' "
            "interpolation, 4 (default) or 2.")
            ->check(CLI::IsMember({2, 4}));
    app.add_option("--starter", starter,
            "--scheme pade: self (default), the wide-angle self-starter, or gaussian.")
            ->check(CLI::IsMember({"self", "gaussian"}));
    app.add_option("--range_ppw", range_points_per_wavelength,
            "--scheme cn: range steps per wavelength (default 10).")->check(CLI::PositiveNumber);
    app.add_option("--gamma_t", gamma_t, "Turbulence in the PE, as for niceops_pe_compare (default 0).");
    app.add_option("--outer_scale", outer_scale, "von Karman outer scale, m (default 10).");
    app.add_option("--seed", seed, "First realisation's seed (default 1).");
    CLI11_PARSE(app, argc, argv);

    Logger::instance().set_file("log/pe_field.log");
    const LengthUnit unit = units == "m" ? metres : feet;
    const double z_s = source_height * unit.to_metres;
    const double r_max = max_range * unit.to_metres;
    const double dr_out = step > 0.0 ? step * unit.to_metres : r_max / 200.0;
    const double z_top = top > 0.0 ? top * unit.to_metres : 2.0 * z_s;
    const double dz = dz_out > 0.0 ? dz_out * unit.to_metres : z_top / 150.0;
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

    GroundModel ground;
    ground.flow_resistance = flow_resistance;
    ground.impedance_model = parse_impedance_model(impedance_model);
    const auto admittance = ground_impedance(ground, bands, sound_speed);

    std::vector<double> heights;
    for (std::size_t j = 0; static_cast<double>(j) * dz <= z_top + 1e-9; j++)
        heights.push_back(static_cast<double>(j) * dz);
    std::vector<PeSetup> setups(bands.size());
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        PeSetup& setup = setups[b];
        setup.source_height = z_s;
        setup.receiver_heights = heights;
        setup.max_range = r_max;
        setup.output_step = dr_out;
        setup.admittance = admittance[b];
        setup.z = layered.heights();
        setup.c = c;
        setup.points_per_wavelength = points_per_wavelength;
        setup.range_points_per_wavelength = range_points_per_wavelength;
        setup.scheme = scheme == "cn" ? PeScheme::crank_nicolson : PeScheme::split_step_pade;
        setup.range_step_wavelengths = range_step_wavelengths;
        setup.max_angle = max_angle;
        setup.vertical_order = vertical_order;
        setup.boundary_order = boundary_order;
        setup.starter = starter == "gaussian" ? PeStarter::gaussian : PeStarter::self;
        setup.gamma_t = gamma_t;
        setup.outer_scale = outer_scale;
        setup.seed = seed + 1000u * static_cast<unsigned>(b);
    }
    std::cerr << "PE field: " << bands.size() << " bands x " << samples << " frequencies, " << heights.size()
              << " heights..." << std::endl;
    const auto pe = solve_pe_bands(setups, bands, samples);

    std::cout << "band,range_" << unit.name << ",height_" << unit.name << ",level_db\n" << std::setprecision(6);
    for (std::size_t b = 0; b < bands.size(); b++)
        for (std::size_t q = 0; q < pe[b].range.size(); q++)
            for (std::size_t k = 0; k < heights.size(); k++)
                std::cout << bands[b] << "," << pe[b].range[q] / unit.to_metres << "," << heights[k] / unit.to_metres
                          << "," << pe[b].relative_level[k][q] << "\n";
    return 0;
}
