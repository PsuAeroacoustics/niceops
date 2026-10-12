// Build an --atmosphere profile from measured columns, or fit a site's
// roughness length to wind profiles (src/atmosphere_builder.hpp).
//
//   niceops_build_atmosphere build --column lidar.csv --column string.csv
//       --humidity 62.5 --z0 0.1 --top 5000 -o atmospheres/run.csv
//
//   niceops_build_atmosphere fit_z0 runs/*.csv -o site_z0.json
//
// build writes the profile in the form --atmosphere reads, every value to 17
// digits, reads it back through the run's own reader as a check, and writes a
// provenance sidecar (OUT.provenance.json): the inputs and their digests, the
// rules, z0 and where it came from, which heights were measured and which
// filled, and every warning.  fit_z0 takes one column file per run (the file
// name is the run's label) and writes the fit as JSON.
//
// Lengths on the command line and in the files are in --units (feet by
// default); a speed names its unit (12kt, 6.2mps).

#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "atmosphere.hpp"
#include "atmosphere_builder.hpp"
#include "atmosphere_profile.hpp"
#include "manifest.hpp"
#include "util/constants.hpp"
#include "util/logger.hpp"
#include "util/paths.hpp"
#include "util/sha256.hpp"

using namespace niceops; // NOLINT

namespace {

using json = nlohmann::json;

std::string trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t");
    if (first == std::string::npos)
        return "";
    const auto last = text.find_last_not_of(" \t");
    return text.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string& text, char separator)
{
    std::vector<std::string> fields;
    std::stringstream stream(text);
    std::string field;
    while (std::getline(stream, field, separator))
        fields.push_back(trim(field));
    return fields;
}

double number(const std::string& text, const std::string& option)
{
    std::size_t used = 0;
    double value = 0.0;
    try
    {
        value = std::stod(text, &used);
    }
    catch (const std::exception&)
    {
        used = 0;
    }
    if (used == 0 || used != text.size() || !std::isfinite(value))
        input_error(option, ": '", text, "' is not a number.");
    return value;
}

// A speed with its unit named: 3mps, 6kt, 7mph.  A bare number is refused, as
// --wind refuses one.
double speed_mps(const std::string& text, const std::string& option)
{
    const char* begin = text.c_str();
    char* end = nullptr;
    const double value = std::strtod(begin, &end);
    if (end == begin)
        input_error(option, ": '", text, "' is not a speed.");
    std::string unit = trim(end);
    std::transform(unit.begin(), unit.end(), unit.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (unit == "mps" || unit == "m/s" || unit == "ms")
        return value;
    if (unit == "kt" || unit == "kts" || unit == "knot" || unit == "knots")
        return value * knot;
    if (unit == "mph")
        return value * 1609.344 / 3600.0;
    if (unit.empty())
    {
        input_error(option, ": name the speed's unit, as ", text, "mps or ", text, "kt. Nothing "
                "about the number tells them apart.");
    }
    input_error(option, ": '", unit, "' is not a speed unit; give mps, kt or mph.");
}

json file_entry(const std::string& role, const std::string& path)
{
    json entry {{"role", role}, {"path", path}};
    std::error_code error;
    const auto size = std::filesystem::file_size(path, error);
    if (!error)
        entry["bytes"] = size;
    if (const auto digest = sha256_of_file(path))
        entry["sha256"] = *digest;
    return entry;
}

json build_json(int argc, char** argv)
{
    const BuildInfo& info = build_info();
    std::vector<std::string> command(argv, argv + argc);
    return {{"version", info.version}, {"compiler", info.compiler},
            {"build_type", info.build_type}, {"command_line", command}};
}

void write_json(const json& document, const std::string& path)
{
    ensure_parent_directory(path);
    std::ofstream out(path);
    out << std::setw(2) << document << '\n';
    out.close();
    if (!out)
        input_error("niceops_build_atmosphere: writing ", path, " failed.");
}

const char* fill_name(BuildRules::Fill fill)
{
    switch (fill)
    {
    case BuildRules::Fill::log:
        return "log";
    case BuildRules::Fill::hold:
        return "hold";
    case BuildRules::Fill::none:
        return "none";
    }
    return "";
}

} // namespace

int main(int argc, char** argv)
{
    CLI::App app {"Build an --atmosphere profile from measured columns, or fit a site's roughness "
                  "length to wind profiles."};
    app.require_subcommand(1);
    std::string units = "ft";
    app.add_option("--units", units, "Length unit of the options and the files: ft (default) or m.")
            ->check(CLI::IsMember({"ft", "m"}));

    // --- build -------------------------------------------------------------
    CLI::App* build = app.add_subcommand("build", "Build a profile from measured columns.");
    std::vector<std::string> columns;
    std::string wind_from, temperature_from, humidity_from, pressure_from;
    std::vector<double> missing;
    bool drop_implausible = false;
    double ground_elevation = 0.0;
    double temperature = 0.0, humidity = 0.0, pressure = 0.0;
    std::string surface_wind;
    std::string fill = "log";
    double z0 = 0.0;
    std::string z0_from;
    std::string levels = "geometric";
    int levels_per_decade = 16;
    double top = 0.0;
    std::string above = "hold";
    std::string stability = "neutral";
    double surface_layer_depth = 0.0;
    std::string humidity_aloft = "interpolate";
    bool write_pressure = false;
    std::string calm_levels = "refuse";
    bool clip_humidity = false;
    std::string output, provenance;
    build->add_option("--column", columns,
            "A column file: every variable it carries (z or alt, T_K/T_C, RH or Td_C or q_kgkg, "
            "p_kPa, wind_speed_kt/_mps with wind_direction or u/v).  Repeatable.");
    build->add_option("--wind_from", wind_from, "Take the wind from this file only.");
    build->add_option("--temperature_from", temperature_from, "Take the temperature from this file only.");
    build->add_option("--humidity_from", humidity_from, "Take the humidity from this file only.");
    build->add_option("--pressure_from", pressure_from, "Take the pressure from this file only.");
    build->add_option("--missing", missing, "Cell values that are fill, not data (9999 99999 ...).");
    build->add_flag("--drop_implausible", drop_implausible,
            "Drop implausible cells, counted, instead of refusing the file.");
    build->add_option("--ground_elevation", ground_elevation,
            "The observers' elevation, in --units, taken off an altitude column.");
    build->add_option("--temperature", temperature, "The surface station's temperature, K, at the ground.");
    build->add_option("--humidity", humidity, "The surface station's relative humidity, percent.");
    build->add_option("--pressure", pressure, "The surface station's pressure, kPa.");
    build->add_option("--wind", surface_wind,
            "The surface station's wind, as niceops --wind: SPEED,DIRECTION[,HEIGHT] (12kt,250,33).");
    build->add_option("--fill", fill,
            "Wind below the lowest level: log (default; the --wind law through it), hold, or none.")
            ->check(CLI::IsMember({"log", "hold", "none"}));
    build->add_option("--z0", z0, "Roughness length, in --units.  This or --z0_from is required.");
    build->add_option("--z0_from", z0_from, "Take the roughness length from a fit_z0 result.");
    build->add_option("--levels", levels,
            "geometric (default), or a comma list of levels in --units, to which every measured "
            "level, the ground and the top are added.");
    build->add_option("--levels_per_decade", levels_per_decade,
            "Geometric levels per factor of ten in height (default 16).");
    build->add_option("--top", top, "The profile's top, in --units (required).");
    build->add_option("--above", above,
            "Temperature above its highest measurement: hold (default) or lapse.")
            ->check(CLI::IsMember({"hold", "lapse"}));
    build->add_option("--stability", stability,
            "With a surface temperature alone and --above lapse: the class whose gradient it "
            "lapses with (and whose power law fills below one wind level).");
    build->add_option("--surface_layer_depth", surface_layer_depth,
            "With one wind level: carry the law up to here, in --units, and hold it above.");
    build->add_option("--humidity_aloft", humidity_aloft,
            "Humidity above its highest measurement: interpolate (held) or held_concentration.")
            ->check(CLI::IsMember({"interpolate", "held_concentration"}));
    build->add_flag("--write_pressure", write_pressure,
            "Write a hydrostatic pressure column from --pressure.");
    build->add_option("--calm_levels", calm_levels, "A measured calm level above the ground: refuse "
            "(default) or drop.")->check(CLI::IsMember({"refuse", "drop"}));
    build->add_flag("--clip_humidity", clip_humidity, "Clip relative humidity above 100 to 100.");
    build->add_option("-o,--output", output, "The profile to write.")->required();
    build->add_option("--provenance", provenance,
            "The provenance sidecar (default: the output's name with .provenance.json).");

    // --- fit_z0 ------------------------------------------------------------
    CLI::App* fit_app = app.add_subcommand("fit_z0",
            "Fit one roughness length to wind profiles over near-neutral runs.");
    std::vector<std::string> run_files;
    std::string gates;
    FitRules fit_rules;
    std::string min_speed;
    std::string fit_output, fit_table;
    std::vector<double> fit_missing;
    bool fit_drop_implausible = false;
    double fit_ground_elevation = 0.0;
    fit_app->add_option("runs", run_files, "One column file per run.")->required();
    fit_app->add_option("--gates", gates, "LOW,HIGH: the gate heights fitted, in --units "
            "(default 46,197 ft).");
    fit_app->add_option("--max_richardson", fit_rules.max_richardson,
            "Largest |bulk Richardson number| across the gates (default 0.05).");
    fit_app->add_option("--min_speed", min_speed,
            "Least speed at the highest gate, with its unit (default 3mps).");
    fit_app->add_option("--max_veer", fit_rules.max_veer,
            "Most turning across the gates, degrees (default 20).");
    fit_app->add_flag("--assume_neutral", fit_rules.assume_neutral,
            "Keep runs without the temperature to screen them, with a warning.");
    fit_app->add_option("--seed", fit_rules.seed, "The bootstrap's seed (default 1).");
    fit_app->add_option("--draws", fit_rules.draws, "Bootstrap draws (default 2000).");
    fit_app->add_option("--min_runs", fit_rules.min_runs, "Fewest runs a fit takes (default 5).")
            ->check(CLI::PositiveNumber);
    fit_app->add_option("--missing", fit_missing, "Cell values that are fill, not data.");
    fit_app->add_flag("--drop_implausible", fit_drop_implausible, "Drop implausible cells.");
    fit_app->add_option("--ground_elevation", fit_ground_elevation,
            "The observers' elevation, in --units, taken off an altitude column.");
    fit_app->add_option("-o,--output", fit_output, "The fit, as JSON.")->required();
    fit_app->add_option("--table", fit_table,
            "Also write the fit as a one-row CSV: z0, its interval (in --units), runs, rms.");

    CLI11_PARSE(app, argc, argv);

    Logger::instance().set_file("log/build_atmosphere.log");
    try
    {
        const LengthUnit unit = units == "m" ? meters : feet;

        if (*fit_app)
        {
            if (fit_app->count("--gates"))
            {
                const auto fields = split(gates, ',');
                if (fields.size() != 2)
                    input_error("--gates takes LOW,HIGH in ", unit.name, ", not '", gates, "'.");
                fit_rules.gate_low = number(fields[0], "--gates") * unit.to_meters;
                fit_rules.gate_high = number(fields[1], "--gates") * unit.to_meters;
            }
            if (!min_speed.empty())
                fit_rules.min_speed = speed_mps(min_speed, "--min_speed");
            if (fit_app->count("--ground_elevation"))
                fit_ground_elevation *= unit.to_meters;
            std::vector<ColumnSource> runs;
            std::vector<std::string> labels;
            for (const auto& path : run_files)
            {
                runs.push_back(read_column_source(path, unit, fit_ground_elevation, fit_missing,
                        fit_drop_implausible));
                labels.push_back(std::filesystem::path(path).stem().string());
            }
            const RoughnessFit fit = fit_roughness_length(runs, labels, fit_rules);
            if (fit.at_bound)
            {
                warn("fit_z0: the fit came to an end of its search (", fit.z0, " m): these profiles "
                        "do not follow the log law.");
            }
            json document = build_json(argc, argv);
            document["z0_m"] = fit.z0;
            document["z0_low_m"] = fit.z0_low;
            document["z0_high_m"] = fit.z0_high;
            document["at_bound"] = fit.at_bound;
            document["residual_rms_mps"] = fit.residual_rms;
            document["median_run_z0_m"] = fit.median_run_z0;
            document["runs_used"] = fit.used;
            json per_run = json::object();
            for (std::size_t j = 0; j < fit.used.size(); j++)
                per_run[fit.used[j]] = std::isfinite(fit.per_run_z0[j]) ? json(fit.per_run_z0[j]) : json();
            document["per_run_z0_m"] = per_run;
            json excluded = json::object();
            for (const auto& [label, why] : fit.excluded)
                excluded[label] = why;
            document["runs_excluded"] = excluded;
            document["rules"] = {{"gate_low_m", fit_rules.gate_low},
                    {"gate_high_m", fit_rules.gate_high},
                    {"max_richardson", fit_rules.max_richardson},
                    {"min_speed_mps", fit_rules.min_speed}, {"max_veer_deg", fit_rules.max_veer},
                    {"assume_neutral", fit_rules.assume_neutral}, {"seed", fit_rules.seed},
                    {"draws", fit_rules.draws}, {"min_runs", fit_rules.min_runs}};
            json inputs = json::array();
            for (const auto& path : run_files)
                inputs.push_back(file_entry("run", path));
            document["inputs"] = inputs;
            document["warnings"] = Logger::instance().first_problems();
            write_json(document, fit_output);
            if (!fit_table.empty())
            {
                ensure_parent_directory(fit_table);
                std::ofstream table(fit_table);
                const std::string suffix = unit.to_meters == 1.0 ? "_m" : "_ft";
                table << "z0" << suffix << ",z0_low" << suffix << ",z0_high" << suffix
                      << ",runs,residual_rms_mps\n"
                      << std::setprecision(17) << fit.z0 / unit.to_meters << ','
                      << fit.z0_low / unit.to_meters << ',' << fit.z0_high / unit.to_meters << ','
                      << fit.used.size() << ',' << fit.residual_rms << '\n';
                table.close();
                if (!table)
                    input_error("fit_z0: writing ", fit_table, " failed.");
            }
            std::cout << "z0 = " << fit.z0 / unit.to_meters << " " << unit.name << " (" << fit.z0
                      << " m), 95% [" << fit.z0_low / unit.to_meters << ", "
                      << fit.z0_high / unit.to_meters << "] " << unit.name << "; "
                      << fit.used.size() << " runs used, " << fit.excluded.size()
                      << " excluded; residual rms " << fit.residual_rms << " m/s\n";
            for (const auto& [label, why] : fit.excluded)
                std::cout << "  excluded " << label << ": " << why << "\n";
            return 0;
        }

        // build.  The unit seam: every length given is converted once, here.
        if (build->count("--z0") == 0 && z0_from.empty())
        {
            input_error("niceops_build_atmosphere build: give the roughness length, --z0 or "
                    "--z0_from. It sets the wind below the lowest measurement, nothing measures "
                    "it, and there is no default.");
        }
        if (build->count("--z0") > 0 && !z0_from.empty())
            input_error("niceops_build_atmosphere build: give --z0 or --z0_from, not both.");
        if (build->count("--top") == 0)
        {
            input_error("niceops_build_atmosphere build: give --top. A run refuses a track above "
                    "the profile's top, so the profile has to say where it ends.");
        }
        BuildRules rules;
        json z0_origin;
        if (!z0_from.empty())
        {
            std::ifstream in(z0_from);
            if (!in)
                input_error("--z0_from: cannot open ", z0_from, ".");
            json fit;
            try
            {
                in >> fit;
            }
            catch (const std::exception& e)
            {
                input_error("--z0_from ", z0_from, " is not JSON: ", e.what());
            }
            if (!fit.contains("z0_m") || !fit["z0_m"].is_number())
                input_error("--z0_from ", z0_from, " holds no z0_m: is it a fit_z0 result?");
            rules.roughness_length = fit["z0_m"].get<double>();
            // A fit at an end of its search found no roughness length, only a limit.
            if (fit.contains("at_bound") && fit["at_bound"].is_boolean() && fit["at_bound"].get<bool>())
            {
                warn("--z0_from ", z0_from, ": that fit came to an end of its search (", rules.roughness_length,
                        " m), so its profiles did not follow the log law and this z0 is a limit, not a fit.");
            }
            z0_origin = file_entry("z0_from", z0_from);
            for (const char* key : {"z0_low_m", "z0_high_m", "runs_used"})
            {
                if (fit.contains(key))
                    z0_origin[key] = fit[key];
            }
        }
        else
        {
            rules.roughness_length = z0 * unit.to_meters;
            z0_origin = "given (--z0)";
        }
        rules.fill = fill == "log" ? BuildRules::Fill::log
                : fill == "hold"   ? BuildRules::Fill::hold
                                   : BuildRules::Fill::none;
        rules.top = top * unit.to_meters;
        rules.levels_per_decade = levels_per_decade;
        if (levels != "geometric")
        {
            for (const auto& field : split(levels, ','))
                rules.levels.push_back(number(field, "--levels") * unit.to_meters);
            if (rules.levels.empty())
                input_error("--levels takes geometric or a comma list of heights.");
        }
        rules.above = above == "lapse" ? BuildRules::Above::lapse : BuildRules::Above::hold;
        rules.stability = parse_stability(stability);
        if (build->count("--surface_layer_depth"))
            rules.surface_layer_depth = surface_layer_depth * unit.to_meters;
        rules.humidity_aloft = humidity_aloft == "held_concentration"
                ? BuildRules::HumidityAloft::held_concentration
                : BuildRules::HumidityAloft::interpolate;
        rules.write_pressure = write_pressure;
        rules.drop_calm = calm_levels == "drop";
        rules.clip_humidity = clip_humidity;
        if (build->count("--ground_elevation"))
            ground_elevation *= unit.to_meters;

        SurfaceAir surface;
        if (build->count("--temperature"))
        {
            if (const auto problem = temperature_problem(temperature); !problem.empty())
                input_error("--temperature: ", problem);
            surface.temperature = temperature;
        }
        if (build->count("--humidity"))
        {
            if (const auto problem = humidity_problem(humidity); !problem.empty())
                input_error("--humidity: ", problem);
            surface.humidity = humidity;
        }
        if (build->count("--pressure"))
        {
            if (const auto problem = pressure_problem(pressure); !problem.empty())
                input_error("--pressure: ", problem);
            surface.pressure = pressure;
        }
        if (!surface_wind.empty())
        {
            SurfaceWind wind = parse_surface_wind(surface_wind, unit);
            if (!std::isfinite(wind.speed) || wind.speed < 0.0 || wind.speed > 50.0
                    || !(wind.direction >= 0.0 && wind.direction <= 360.0))
            {
                input_error("--wind: a surface wind of ", wind.speed, " m/s from ", wind.direction,
                        " degrees is not plausible (0 to 50 m/s, from 0 to 360 degrees).");
            }
            surface.wind = wind;
        }

        // Every file once, and each variable from one of them.
        std::map<std::string, ColumnSource> sources;
        json inputs = json::array();
        auto source_of = [&](const std::string& path, const std::string& role) -> ColumnSource& {
            inputs.push_back(file_entry(role, path));
            auto found = sources.find(path);
            if (found == sources.end())
            {
                found = sources.emplace(path, read_column_source(path, unit, ground_elevation,
                                                      missing, drop_implausible))
                                .first;
            }
            return found->second;
        };
        for (const auto& path : columns)
            source_of(path, "column");
        MeasuredColumns measured;
        auto pick = [&](const std::string& named, const char* option, const char* variable,
                            auto member, auto& slot) {
            if (!named.empty())
            {
                const ColumnSource& source = source_of(named, option + 2);
                if (!(source.*member))
                    input_error(option, " ", named, " gives no ", variable, ".");
                slot = source.*member;
                return;
            }
            std::vector<std::string> offering;
            for (const auto& path : columns)
            {
                if (sources.at(path).*member)
                    offering.push_back(path);
            }
            std::sort(offering.begin(), offering.end());
            offering.erase(std::unique(offering.begin(), offering.end()), offering.end());
            if (offering.size() > 1)
            {
                input_error("niceops_build_atmosphere: the ", variable, " comes from both ",
                        offering[0], " and ", offering[1], "; name the one to take with ", option,
                        ".");
            }
            if (offering.size() == 1)
                slot = sources.at(offering.front()).*member;
        };
        pick(temperature_from, "--temperature_from", "temperature", &ColumnSource::temperature,
                measured.temperature);
        pick(humidity_from, "--humidity_from", "humidity", &ColumnSource::humidity,
                measured.humidity);
        pick(pressure_from, "--pressure_from", "pressure", &ColumnSource::pressure,
                measured.pressure);
        pick(wind_from, "--wind_from", "wind", &ColumnSource::wind, measured.wind);

        BuildReport report;
        const AtmosphereProfile profile = build_profile(measured, surface, rules, unit, report);
        write_atmosphere_profile(profile, output, unit);

        // The self-check: the run's reader gives the same profile back, the
        // heights to an ULP of the unit conversion and the values exactly.
        const AtmosphereProfile read = read_atmosphere_profile(output, unit, 0.0);
        bool same = read.height.size() == profile.height.size()
                && read.temperature == profile.temperature
                && read.relative_humidity == profile.relative_humidity
                && read.pressure == profile.pressure && read.wind_east == profile.wind_east
                && read.wind_north == profile.wind_north;
        for (std::size_t i = 0; same && i < profile.height.size(); i++)
            same = std::fabs(read.height[i] - profile.height[i]) <= 1e-15 * profile.height[i];
        if (!same)
            throw std::logic_error(output + " does not read back as the profile written");

        json document = build_json(argc, argv);
        document["output"] = output;
        document["inputs"] = inputs;
        json z0_json = {{"z0_m", rules.roughness_length}, {"origin", z0_origin}};
        document["roughness_length"] = z0_json;
        document["rules"] = {{"units", unit.name}, {"fill", fill_name(rules.fill)},
                {"levels", levels}, {"levels_per_decade", rules.levels_per_decade},
                {"top_m", rules.top}, {"above", above}, {"stability", describe(rules.stability)},
                {"surface_layer_depth_m",
                        rules.surface_layer_depth ? json(*rules.surface_layer_depth) : json()},
                {"humidity_aloft", humidity_aloft}, {"write_pressure", rules.write_pressure},
                {"calm_levels", calm_levels}, {"clip_humidity", rules.clip_humidity},
                {"missing", missing}, {"drop_implausible", drop_implausible},
                {"ground_elevation_m", ground_elevation}};
        json surface_json = json::object();
        if (surface.temperature)
            surface_json["temperature_K"] = *surface.temperature;
        if (surface.humidity)
            surface_json["humidity_pct"] = *surface.humidity;
        if (surface.pressure)
            surface_json["pressure_kPa"] = *surface.pressure;
        if (surface.wind)
        {
            surface_json["wind"] = {{"speed_mps", surface.wind->speed},
                    {"direction_deg", surface.wind->direction},
                    {"height_m", surface.wind->height}};
        }
        document["surface"] = surface_json;
        json sources_json = json::object();
        if (measured.temperature)
            sources_json["temperature"] = measured.temperature->source;
        if (measured.humidity)
            sources_json["humidity"] = measured.humidity->source;
        if (measured.pressure)
            sources_json["pressure"] = measured.pressure->source;
        if (measured.wind)
            sources_json["wind"] = measured.wind->source;
        document["variables_from"] = sources_json;
        json spans = json::object();
        for (const auto& [name, span] : report.measured_span)
            spans[name] = {{"low_m", span.first}, {"high_m", span.second}};
        document["measured_m"] = spans;
        json filled = json::array();
        for (const auto& [low, high] : report.filled_wind)
            filled.push_back({{"low_m", low}, {"high_m", high}});
        document["wind_filled_m"] = filled;
        std::size_t sentinels = 0, implausible = 0;
        for (const auto& [path, source] : sources)
        {
            sentinels += source.sentinels_dropped;
            implausible += source.implausible_dropped;
        }
        document["counts"] = {{"fill_values_dropped", sentinels},
                {"implausible_dropped", implausible},
                {"below_ground_dropped", report.below_ground_dropped},
                {"calm_dropped", report.calm_dropped},
                {"humidity_clipped", report.humidity_clipped}, {"levels", profile.height.size()}};
        document["warnings"] = Logger::instance().first_problems();
        const std::string sidecar = provenance.empty() ? output + ".provenance.json" : provenance;
        write_json(document, sidecar);
        std::cerr << "Wrote " << output << " (" << profile.height.size() << " levels to "
                  << rules.top / unit.to_meters << " " << unit.name << ") and " << sidecar << "."
                  << std::endl;
    }
    catch (const InputError&)
    {
        // Already logged, to the file and stderr.
        return 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "niceops_build_atmosphere: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
