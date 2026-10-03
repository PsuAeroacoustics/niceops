// Build a propagation loss table for niceops --propagation_table
// (src/propagation_table.hpp).
//
//   niceops_propagation_table --method pe --atmosphere profile.csv
//       --database databases/B407_fixed_load_factor.nod --band_limits 50 5000
//       --source_heights 33 66 100 200 330 660 1000 --receiver_heights 0 3.937
//       --max_range 16000 --step 33 --bearings 36 -o tables/windy.nc
//
// For each bearing and source height the method gives each band's level
// relative to free field at every range and receiver height:
//   pe            the wide-angle parabolic equation (pe.hpp) through the
//                 profile's effective sound speed along the bearing,
//                 band-averaged over --samples frequencies; NaN steeper than
//                 --max_angle from the source
//   ground_model  straight rays and the spherical-wave ground model niceops
//                 --ground_effect uses (ground_effect.hpp): the same footprint
//                 as --ground_effect, up to interpolation, which is what makes
//                 it the check on the table's plumbing.  It knows no
//                 refraction, so it takes one bearing.
// Neither includes atmospheric absorption; niceops applies its own.  Lengths
// on the command line are in --units (feet by default); the file is SI.

#include <CLI/CLI.hpp>

#include <algorithm>
#include <chrono>
#include <limits>
#include <optional>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "atmosphere_profile.hpp"
#include "ground_effect.hpp"
#include "pe.hpp"
#include "propagation_table.hpp"
#include "sphere_database.hpp"
#include "stratified.hpp"
#include "util/constants.hpp"
#include "util/logger.hpp"
#include "util/parallel.hpp"

using namespace niceops; // NOLINT

namespace {

std::string join(const std::vector<double>& values)
{
    std::ostringstream out;
    for (std::size_t k = 0; k < values.size(); k++)
        out << (k ? " " : "") << values[k];
    return out.str();
}

} // namespace

int main(int argc, char** argv)
{
    CLI::App app {"Build a propagation loss table for niceops --propagation_table."};
    std::string method = "pe";
    std::string atmosphere_file;
    std::string database_file;
    std::string output;
    std::string units = "ft";
    std::vector<double> bands;
    std::vector<double> band_limits;
    std::vector<double> source_heights;
    std::vector<double> receiver_heights;
    double max_range = 0.0;
    double step = 0.0;
    int bearings = 36;
    double flow_resistance = flow_resistance::grassy_field;
    std::string impedance_model = "delany_bazley";
    double porosity_rate = 0.0;
    double layer_depth = 0.05;
    int samples = 5;
    double max_angle = 50.0;
    double points_per_wavelength = 0.0;
    double sound_speed_override = 0.0;
    app.add_option("--method", method, "pe (default) or ground_model.")->check(CLI::IsMember({"pe", "ground_model"}));
    app.add_option("--atmosphere", atmosphere_file,
            "Atmosphere profile, as for niceops; required for pe, and for ground_model gives the sound speed.");
    app.add_option("--sound_speed", sound_speed_override,
            "Sound speed at the ground, m/s, for the ground's impedance and the ground model (default: the "
            "profile's at the ground, else 343.2).  Match the run's: niceops logs it.");
    app.add_option("--database", database_file, "Take the band centers from this database.");
    app.add_option("--bands", bands, "Band centers, Hz, instead of --database's.");
    app.add_option("--band_limits", band_limits, "Keep only bands from LOW to HIGH Hz.")->expected(2);
    app.add_option("--source_heights", source_heights, "Source heights above the ground, in --units.")->required();
    app.add_option("--receiver_heights", receiver_heights, "Receiver heights above the ground, in --units.")
            ->required();
    app.add_option("--max_range", max_range, "Farthest horizontal range, in --units.")->required()->check(CLI::PositiveNumber);
    app.add_option("--step", step, "Range spacing, in --units.")->required()->check(CLI::PositiveNumber);
    app.add_option("--bearings", bearings,
            "Bearings, evenly from 0 over 360 degrees (default 36; 1 when the profile has no wind).");
    app.add_option("--units", units, "Length unit: ft (default) or m.")->check(CLI::IsMember({"ft", "m"}));
    app.add_option("--flow_resistance", flow_resistance, "Ground flow resistance, kPa s/m^2.");
    app.add_option("--impedance_model", impedance_model, "As for niceops.");
    app.add_option("--porosity_rate", porosity_rate, "As for niceops, 1/m.");
    app.add_option("--layer_depth", layer_depth, "As for niceops, m.");
    app.add_option("--samples", samples, "pe: frequencies per band (default 5).")->check(CLI::PositiveNumber);
    app.add_option("--max_angle", max_angle, "pe: angle limit, degrees (default 50).")->check(CLI::Range(0.1, 60.0));
    app.add_option("--ppw", points_per_wavelength, "pe: grid points per wavelength in height (default: the PE's).")->check(CLI::NonNegativeNumber);
    app.add_option("-o,--output", output, "The table to write (netCDF).")->required();
    CLI11_PARSE(app, argc, argv);

    Logger::instance().set_file("log/propagation_table.log");
    try
    {
        const LengthUnit unit = units == "m" ? meters : feet;
        for (double& h : source_heights)
            h *= unit.to_meters;
        for (double& h : receiver_heights)
            h *= unit.to_meters;
        std::sort(source_heights.begin(), source_heights.end());
        std::sort(receiver_heights.begin(), receiver_heights.end());
        const double r_max = max_range * unit.to_meters;
        const double dr = step * unit.to_meters;
        if (!(dr > 0.0) || !(r_max >= dr))
            throw std::runtime_error("--step must be positive and no more than --max_range");

        if (bands.empty())
        {
            if (database_file.empty())
                throw std::runtime_error("give --bands or --database");
            const Database database(database_file);
            if (!database.has_spectra())
                throw std::runtime_error(database_file + " carries no spectra");
            bands = database.band_centers();
        }
        if (!band_limits.empty())
        {
            std::vector<double> kept;
            for (const double f : bands)
                if (f >= band_limits[0] * std::pow(2.0, -1.0 / 6.0) && f <= band_limits[1] * std::pow(2.0, 1.0 / 6.0))
                    kept.push_back(f);
            bands = kept;
        }
        std::sort(bands.begin(), bands.end());
        if (bands.empty())
            throw std::runtime_error("no bands left");

        std::optional<AtmosphereProfile> profile;
        double sound_speed = 343.2;
        if (!atmosphere_file.empty())
        {
            profile = read_atmosphere_profile(atmosphere_file, unit, 0.0);
            if (!profile->has_pressure())
                profile = profile->with_hydrostatic_pressure(reference_pressure);
            sound_speed = profile->at(0.0).sound_speed();
        }
        if (sound_speed_override > 0.0)
            sound_speed = sound_speed_override;
        if (method == "pe" && !profile)
            throw std::runtime_error("--method pe needs --atmosphere");
        if (method == "ground_model" || !profile->has_wind())
            bearings = 1;
        if (bearings < 1)
            throw std::runtime_error("--bearings must be at least 1");

        GroundModel ground;
        ground.flow_resistance = flow_resistance;
        ground.impedance_model = parse_impedance_model(impedance_model);
        ground.porosity_rate = porosity_rate;
        ground.layer_depth = layer_depth;
        const auto admittance = ground_impedance(ground, bands, sound_speed);

        PropagationTable::Axes axes;
        for (int b = 0; b < bearings; b++)
            axes.bearings.push_back(360.0 * b / bearings);
        axes.source_heights = source_heights;
        axes.receiver_heights = receiver_heights;
        for (std::size_t q = 1; q * dr <= r_max + 1e-9; q++)
            axes.ranges.push_back(q * dr);
        axes.bands = bands;
        const std::size_t nb = bands.size(), nr = axes.ranges.size(), nh = receiver_heights.size(),
                          ns = source_heights.size();
        std::vector<double> values(axes.bearings.size() * ns * nh * nr * nb, 0.0);
        auto at = [&](std::size_t b, std::size_t s, std::size_t h, std::size_t r, std::size_t band) -> double& {
            return values[(((b * ns + s) * nh + h) * nr + r) * nb + band];
        };

        std::map<std::string, std::string> provenance;
        provenance["atmosphere"] = atmosphere_file.empty() ? "none" : atmosphere_file;
        provenance["sound_speed_m_s"] = std::to_string(sound_speed);
        provenance["ground"] = impedance_model + ", flow resistance " + std::to_string(flow_resistance)
                + " kPa s/m^2, porosity rate " + std::to_string(porosity_rate) + " /m, layer depth "
                + std::to_string(layer_depth) + " m";
        provenance["command"] = [&] {
            std::string line;
            for (int k = 0; k < argc; k++)
                line += (k ? " " : "") + std::string(argv[k]);
            return line;
        }();

        const auto start = std::chrono::steady_clock::now();
        std::string description;
        if (method == "ground_model")
        {
            description = "straight rays and the spherical-wave ground model (ground_effect.hpp)";
            parallel_chunks(ns * nh, default_thread_count(), [&](std::size_t lo, std::size_t hi) {
                for (std::size_t job = lo; job < hi; job++)
                {
                    const std::size_t s = job / nh, h = job % nh;
                    GroundModel here = ground;
                    here.receiver_height = receiver_heights[h];
                    for (std::size_t r = 0; r < nr; r++)
                    {
                        const auto gain = excess_ground_attenuation(here, source_heights[s], axes.ranges[r],
                                sound_speed, bands);
                        for (std::size_t band = 0; band < nb; band++)
                            at(0, s, h, r, band) = gain[band];
                    }
                }
            });
        }
        else
        {
            std::ostringstream text;
            text << "wide-angle parabolic equation (pe.hpp), split-step Pade, self-starter, fourth order in "
                    "height; "
                 << samples << " frequencies a band; max_angle " << max_angle << " degrees";
            description = text.str();
            const LayeredSoundSpeed layered(*profile, 90.0);
            std::vector<double> c(layered.levels());
            for (std::size_t b = 0; b < axes.bearings.size(); b++)
            {
                const double bearing = axes.bearings[b] * pi / 180.0;
                layered.along(std::sin(bearing), std::cos(bearing), c.data());
                for (std::size_t s = 0; s < ns; s++)
                {
                    std::vector<PeSetup> setups(nb);
                    for (std::size_t band = 0; band < nb; band++)
                    {
                        PeSetup& setup = setups[band];
                        setup.source_height = source_heights[s];
                        setup.receiver_heights = receiver_heights;
                        setup.max_range = r_max;
                        setup.output_step = dr;
                        setup.admittance = admittance[band];
                        setup.z = layered.heights();
                        setup.c = c;
                        setup.max_angle = max_angle;
                        setup.points_per_wavelength = points_per_wavelength;
                    }
                    const auto results = solve_pe_bands(setups, bands, samples);
                    for (std::size_t band = 0; band < nb; band++)
                        for (std::size_t h = 0; h < nh; h++)
                            for (std::size_t r = 0; r < nr; r++)
                                at(b, s, h, r, band) = r < results[band].range.size()
                                        ? results[band].relative_level[h][r]
                                        : std::numeric_limits<double>::quiet_NaN();
                    std::cerr << "bearing " << axes.bearings[b] << ", source " << source_heights[s] << " m: "
                              << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
                              << " s" << std::endl;
                }
            }
        }
        provenance["bands_hz"] = join(bands);
        const PropagationTable table(std::move(axes), std::move(values), description, false, provenance);
        table.save(output);
        std::cerr << "Wrote " << output << " (" << table.axes().bearings.size() << " bearings x " << ns
                  << " source heights x " << nh << " receiver heights x " << nr << " ranges x " << nb << " bands) in "
                  << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << " s."
                  << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "niceops_propagation_table: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
