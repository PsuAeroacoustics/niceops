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
//
// --march receiver (the default) builds the PE table by reciprocity, one march
// per bearing and receiver height with the source at the receiver's height,
// and reads every source height from the field; --march source marches once
// per source height, as the tool always did (pe_table.hpp has why they agree).
// --gamma_t, --outer_scale, --realizations, --seed and --gamma_profile add
// turbulence by phase screens (pe.hpp): a report-only reference, off by
// default.

#include <CLI/CLI.hpp>

#include <algorithm>
#include <chrono>
#include <limits>
#include <optional>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "atmosphere_profile.hpp"
#include "ground_effect.hpp"
#include "pe.hpp"
#include "pe_table.hpp"
#include "propagation_table.hpp"
#include "sphere_database.hpp"
#include "stratified.hpp"
#include "util/constants.hpp"
#include "util/csv.hpp"
#include "util/logger.hpp"
#include "util/parallel.hpp"
#include "util/paths.hpp"

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
    std::string march = "receiver";
    double gamma_t = 0.0;
    double outer_scale = 10.0;
    int realizations = 1;
    unsigned seed = 1;
    std::string gamma_profile;
    std::string validation_csv;
    int validation_samples = 21;
    double validation_max_error = 1.0;
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
    app.add_option("--validation_csv", validation_csv,
            "Check cell midpoints against direct solves and write errors in dB (SI coordinates). "
            "Exit 2 if a covered reference exceeds --validation_max_error or a table gap hides it.");
    app.add_option("--validation_samples", validation_samples,
            "pe validation: frequencies per band in the reference (default 21).")
            ->check(CLI::PositiveNumber);
    app.add_option("--validation_max_error", validation_max_error,
            "Largest permitted absolute midpoint error in dB (default 1).")
            ->check(CLI::NonNegativeNumber);
    app.add_option("--max_angle", max_angle, "pe: angle limit, degrees (default 50).")->check(CLI::Range(0.1, 60.0));
    app.add_option("--ppw", points_per_wavelength, "pe: grid points per wavelength in height (default: the PE's).")->check(CLI::NonNegativeNumber);
    app.add_option("--march", march,
            "pe: receiver (default) marches once per receiver height and reads the source heights by "
            "reciprocity; source marches once per source height.")
            ->check(CLI::IsMember({"receiver", "source"}));
    app.add_option("--gamma_t", gamma_t,
            "pe: turbulence strength (HARMONOISE gamma_t, m^-2/3) for phase screens; 0 (default) is still air.  "
            "With --gamma_profile, the reference the screens are drawn at (default: the profile's largest).")
            ->check(CLI::NonNegativeNumber);
    app.add_option("--outer_scale", outer_scale, "pe: turbulence outer scale, m (default 10).")
            ->check(CLI::PositiveNumber);
    app.add_option("--realizations", realizations,
            "pe: realizations per frequency sample under turbulence, energy-averaged (default 1).")
            ->check(CLI::PositiveNumber);
    app.add_option("--seed", seed, "pe: first seed of the turbulence (default 1).");
    app.add_option("--gamma_profile", gamma_profile,
            "pe: CSV of z (in --units) and gamma (m^-2/3) for a strength that varies with height: each screen "
            "row scaled by sqrt(gamma(z) / gamma_t), an approximation to a field whose statistics change with "
            "height.");
    app.add_option("-o,--output", output, "The table to write (netCDF).")->required();
    CLI11_PARSE(app, argc, argv);

    Logger::instance().set_file("log/propagation_table.log");
    try
    {
        if (!validation_csv.empty())
        {
            if (!std::isfinite(validation_max_error))
                throw std::runtime_error("--validation_max_error must be finite");
            if (std::filesystem::weakly_canonical(validation_csv)
                    == std::filesystem::weakly_canonical(output))
                throw std::runtime_error("--validation_csv must differ from the table output");
        }
        const LengthUnit unit = units == "m" ? meters : feet;
        std::vector<double> turbulence_z, turbulence_gamma;
        if (!gamma_profile.empty())
        {
            for (const auto& row : read_numeric_csv(gamma_profile, false))
            {
                if (row.size() < 2)
                    throw std::runtime_error(gamma_profile + ": each row needs z and gamma");
                turbulence_z.push_back(row[0] * unit.to_meters);
                turbulence_gamma.push_back(row[1]);
            }
            if (turbulence_z.empty())
                throw std::runtime_error(gamma_profile + " holds no rows");
            for (std::size_t k = 0; k < turbulence_z.size(); k++)
                if (!std::isfinite(turbulence_gamma[k]) || turbulence_gamma[k] < 0.0
                        || (k > 0 && !(turbulence_z[k] > turbulence_z[k - 1])))
                    throw std::runtime_error(gamma_profile + ": heights must increase and gamma be finite and >= 0");
            if (gamma_t == 0.0)
                gamma_t = *std::max_element(turbulence_gamma.begin(), turbulence_gamma.end());
            if (!(gamma_t > 0.0))
                throw std::runtime_error(gamma_profile + " is zero everywhere");
        }
        if (method != "pe" && (gamma_t > 0.0 || march != "receiver"))
            throw std::runtime_error("--march and turbulence apply to --method pe only");
        if (gamma_t > 0.0 && !validation_csv.empty())
            throw std::runtime_error("--validation_csv checks a still table; leave out --gamma_t");
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
                 << samples << " frequencies a band; max_angle " << max_angle << " degrees; "
                 << (march == "receiver" ? "one march per receiver height, source heights by reciprocity"
                                         : "one march per source height");
            if (gamma_t > 0.0)
                text << "; turbulence by phase screens, gamma_t " << gamma_t << " m^-2/3"
                     << (turbulence_z.empty() ? "" : " (reference; strength varies with height)") << ", outer scale "
                     << outer_scale << " m, " << realizations << " realizations a frequency, energy mean";
            description = text.str();
            const LayeredSoundSpeed layered(*profile, 90.0);
            PeTableRequest request;
            request.bearings = axes.bearings;
            request.z = layered.heights();
            for (std::size_t b = 0; b < axes.bearings.size(); b++)
            {
                const double bearing = axes.bearings[b] * pi / 180.0;
                std::vector<double> c(layered.levels());
                layered.along(std::sin(bearing), std::cos(bearing), c.data());
                request.c_by_bearing.push_back(std::move(c));
            }
            request.source_heights = source_heights;
            request.receiver_heights = receiver_heights;
            request.max_range = r_max;
            request.step = dr;
            request.bands = bands;
            request.admittance = admittance;
            request.samples = samples;
            request.max_angle = max_angle;
            request.points_per_wavelength = points_per_wavelength;
            request.march = march == "receiver" ? PeTableMarch::receiver : PeTableMarch::source;
            request.gamma_t = gamma_t;
            request.outer_scale = outer_scale;
            request.realizations = realizations;
            request.seed = seed;
            request.turbulence_z = turbulence_z;
            request.turbulence_gamma = turbulence_gamma;
            values = build_pe_table(request, [&](const std::string& what) {
                std::cerr << what << ": "
                          << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << " s"
                          << std::endl;
            });
            provenance["march"] = march;
            if (gamma_t > 0.0)
            {
                provenance["gamma_t"] = std::to_string(gamma_t);
                provenance["outer_scale_m"] = std::to_string(outer_scale);
                provenance["realizations"] = std::to_string(realizations);
                provenance["seed"] = std::to_string(seed);
                if (!gamma_profile.empty())
                    provenance["gamma_profile"] = gamma_profile;
            }
        }
        provenance["bands_hz"] = join(bands);
        const PropagationTable table(std::move(axes), std::move(values), description, false, provenance);
        table.save(output);
        if (!validation_csv.empty())
        {
            ensure_parent_directory(validation_csv);
            std::ofstream report(validation_csv);
            if (!report)
                throw std::runtime_error("cannot write " + validation_csv);
            report << std::setprecision(17)
                   << "bearing_deg,source_height_m,receiver_height_m,range_m,band_hz,table_db,reference_db,error_db,status\n";
            const auto midpoint = [](const std::vector<double>& axis) {
                std::vector<double> points;
                if (axis.size() == 1)
                    points = axis;
                else
                    for (std::size_t i = 0; i + 1 < axis.size(); i++)
                        points.push_back(0.5 * (axis[i] + axis[i + 1]));
                return points;
            };
            const auto sources = midpoint(table.axes().source_heights);
            const auto receivers = midpoint(table.axes().receiver_heights);
            const auto ranges = midpoint(table.axes().ranges);
            const auto band_map = table.match_bands(bands, 0, nb);
            struct Errors { std::size_t count = 0, gaps = 0; double sum = 0.0, square = 0.0, maximum = 0.0; };
            std::vector<Errors> errors(nb);
            std::size_t unavailable = 0;
            std::optional<LayeredSoundSpeed> reference_profile;
            if (method == "pe")
                reference_profile.emplace(*profile, 90.0);
            for (std::size_t b = 0; b < table.axes().bearings.size(); b++)
            {
                const double bearing = table.isotropic() ? 0.0 : table.axes().bearings[b] + 180.0 / bearings;
                std::vector<double> c;
                if (reference_profile)
                {
                    c.resize(reference_profile->levels());
                    reference_profile->along(std::sin(bearing * pi / 180.0), std::cos(bearing * pi / 180.0), c.data());
                }
                for (const double source : sources)
                {
                    std::vector<PeResult> reference;
                    if (method == "pe")
                    {
                        std::vector<PeSetup> setups(nb);
                        for (std::size_t band = 0; band < nb; band++)
                        {
                            auto& setup = setups[band];
                            setup.source_height = source;
                            setup.receiver_heights = receivers;
                            setup.max_range = r_max;
                            setup.output_step = dr / 2.0;
                            setup.admittance = admittance[band];
                            setup.z = reference_profile->heights();
                            setup.c = c;
                            setup.max_angle = max_angle;
                            setup.points_per_wavelength = points_per_wavelength > 0.0 ? 2.0 * points_per_wavelength : 20.0;
                        }
                        reference = solve_pe_bands(setups, bands, validation_samples);
                    }
                    for (std::size_t h = 0; h < receivers.size(); h++)
                        for (std::size_t r = 0; r < ranges.size(); r++)
                        {
                            std::vector<double> interpolated(nb, std::numeric_limits<double>::quiet_NaN());
                            const bool covered = table.lookup(bearing, source, receivers[h], ranges[r], band_map, interpolated.data());
                            GroundModel here = ground;
                            here.receiver_height = receivers[h];
                            const auto exact = method == "ground_model"
                                    ? excess_ground_attenuation(here, source, ranges[r], sound_speed, bands)
                                    : std::vector<double>();
                            for (std::size_t band = 0; band < nb; band++)
                            {
                                // PE output starts at dr/2, so the first range midpoint (1.5 dr) is index 2.
                                const std::size_t q = table.axes().ranges.size() == 1 ? 1 : 2 * r + 2;
                                const double ref = method == "ground_model" ? exact[band]
                                        : q < reference[band].range.size() ? reference[band].relative_level[h][q]
                                                                           : std::numeric_limits<double>::quiet_NaN();
                                const bool valid_reference = std::isfinite(ref);
                                const double error = covered && valid_reference ? interpolated[band] - ref
                                                                                : std::numeric_limits<double>::quiet_NaN();
                                auto& stats = errors[band];
                                if (!valid_reference)
                                    unavailable++;
                                else if (!covered)
                                    stats.gaps++;
                                else
                                {
                                    stats.count++;
                                    stats.sum += error;
                                    stats.square += error * error;
                                    stats.maximum = std::max(stats.maximum, std::abs(error));
                                }
                                report << bearing << ',' << source << ',' << receivers[h] << ',' << ranges[r] << ','
                                       << bands[band] << ',' << (covered ? interpolated[band] : std::numeric_limits<double>::quiet_NaN())
                                       << ',' << ref << ',' << error << ','
                                       << (!valid_reference ? "reference_unavailable" : covered ? "covered" : "table_gap") << '\n';
                            }
                        }
                    std::cerr << "Validated bearing " << bearing << ", source " << source << " m\n";
                }
            }
            report.close();
            if (!report)
                throw std::runtime_error("failed writing " + validation_csv);
            bool failed = false;
            for (std::size_t band = 0; band < nb; band++)
            {
                const auto& stats = errors[band];
                std::cerr << bands[band] << " Hz: " << stats.count << " covered, " << stats.gaps << " gaps; bias "
                          << (stats.count ? stats.sum / stats.count : 0.0) << ", rms "
                          << (stats.count ? std::sqrt(stats.square / stats.count) : 0.0)
                          << ", max " << stats.maximum << " dB\n";
                failed = failed || stats.count == 0 || stats.gaps > 0 || stats.maximum > validation_max_error;
            }
            std::cerr << unavailable << " reference values unavailable; midpoint validation " << (failed ? "FAILED" : "passed")
                      << ". Table retained at " << output << ".\n";
            if (failed)
                return 2;
        }
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
