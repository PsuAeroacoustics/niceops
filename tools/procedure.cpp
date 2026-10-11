// Fly a procedure into a track niceops reads (src/procedure.hpp).
//
//   niceops_procedure approach.json -o tracks/approach.csv
//   niceops -t tracks/approach.csv -d databases/B407_fixed_load_factor.nod ...
//
// The procedure is JSON: a start (lat, lon, alt, speed, course) and a list of
// segments, every length and speed naming its unit in its key.  The manual's
// Companion tools section has the format.  The track is WGS84 with velocity and
// heading columns, lengths in --units, its alt on the start's datum.
//
// --check DATABASE --vehicle VEHICLE reads the track back as niceops reads it
// and maps every state to the database's conditions as niceops maps them
// (state_condition, Database::blend_at), reporting the states the database
// clamps on each axis or holds at a load-factor limit, and when: a procedure
// can be tuned to its database before any acoustics are run.

#include <CLI/CLI.hpp>

#include <array>
#include <cstdio>
#include <exception>
#include <string>
#include <utility>
#include <vector>

#include "condition_selection.hpp"
#include "procedure.hpp"
#include "sphere_database.hpp"
#include "state_condition.hpp"
#include "trajectory.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

// The spans of time, [first, last] state times, over which `flagged` holds.
std::vector<std::pair<double, double>> spans(const std::vector<State>& states, const std::vector<unsigned>& flags,
        unsigned flag)
{
    std::vector<std::pair<double, double>> result;
    for (std::size_t i = 0; i < states.size(); i++)
    {
        if (!(flags[i] & flag))
            continue;
        if (i > 0 && (flags[i - 1] & flag))
            result.back().second = states[i].t;
        else
            result.emplace_back(states[i].t, states[i].t);
    }
    return result;
}

// Reads `track` back as niceops would (WGS84, at -s's default, which never
// coarsens the track's own step) and
// reports, per axis, the states `database` clamps.
void check_envelope(const std::string& track, LengthUnit unit, const std::string& database_file,
        const std::string& vehicle_file, double max_load_factor, bool qsam, bool accept_vehicle_mismatch)
{
    GeodeticInput geodetic;
    geodetic.resample_dt = 0.5;
    Trajectory trajectory(track, vehicle_file, geodetic, unit);
    const Database database(database_file);
    check_vehicle(trajectory.vehicle, database, vehicle_file, accept_vehicle_mismatch);
    BlendOptions blend;
    blend.tip_speed = trajectory.vehicle.tip_speed;
    const auto& states = trajectory.states;
    // The flags reported, in order; the --min_load_factor floor is not one.
    const std::pair<unsigned, const char*> reported[] = {
            {advance_ratio_clamped, "advance ratio clamped to the database's range"},
            {flight_path_angle_clamped, "flight path angle clamped to the database's range"},
            {thrust_coefficient_clamped, "thrust coefficient clamped to the database's range"},
            {load_held_at_max, "load factor held at --max_load_factor"}};
    const unsigned any_reported = advance_ratio_clamped | flight_path_angle_clamped | thrust_coefficient_clamped
            | load_held_at_max;
    std::vector<unsigned> flags(states.size());
    std::size_t any = 0;
    for (std::size_t i = 0; i < states.size(); i++)
    {
        const auto condition = state_condition(states[i], database, qsam, max_load_factor, blend);
        const auto at = database.blend_at(condition.advance_ratio, condition.flight_path_angle,
                condition.thrust_coefficient, blend);
        flags[i] = condition_flags(condition, at);
        any += (flags[i] & any_reported) != 0;
    }
    std::printf("Envelope of %s: %zu of %zu states (every %.2f s) clamped or held.\n", database_file.c_str(), any,
            states.size(), states.size() > 1 ? states[1].t - states[0].t : 0.0);
    for (const auto& [flag, name] : reported)
    {
        if (flag == load_held_at_max && !database.fixed_load_factor())
            continue;
        const auto found = spans(states, flags, flag);
        std::size_t count = 0;
        for (const unsigned f : flags)
            count += (f & flag) != 0;
        std::printf("  %s: %zu states", name, count);
        for (std::size_t k = 0; k < found.size(); k++)
            std::printf("%s t = %.1f-%.1f s", k == 0 ? "," : ";", found[k].first, found[k].second);
        std::printf("\n");
    }

}

} // namespace

int main(int argc, char* argv[])
{
    CLI::App app {"Fly a procedure into a WGS84 track for niceops."};
    std::string input, output;
    std::string unit_name = "ft";
    app.add_option("procedure", input, "Procedure file (JSON).")->required();
    app.add_option("-o,--output", output, "Track to write (CSV).")->required();
    app.add_option("--units", unit_name, "Length unit of the track: ft (default) or m.")
            ->check(CLI::IsMember({"ft", "m"}));
    std::string check_database, vehicle;
    double max_load_factor = 2.0;
    bool qsam = true;
    bool accept_vehicle_mismatch = false;
    app.add_option("--check", check_database,
            "Report the states this database (.nod) would clamp on each axis, or hold at "
            "--max_load_factor, as niceops maps them. Needs --vehicle.");
    app.add_option("--vehicle", vehicle, "Vehicle file (JSON) for --check.");
    app.add_option("--max_load_factor", max_load_factor,
            "For --check: niceops's --max_load_factor, in g (default 2).")
            ->check(CLI::PositiveNumber);
    app.add_flag("--qsam,!--no-qsam", qsam, "For --check: niceops's --qsam (default on).");
    app.add_flag("--accept_vehicle_mismatch", accept_vehicle_mismatch,
            "For --check: map the states through a vehicle file that disagrees with the database, "
            "as niceops's --accept_vehicle_mismatch.");
    CLI11_PARSE(app, argc, argv);
    // The database and track readers narrate at info; the tool's own report
    // is on stdout.
    Logger::instance().set_level(LogLevel::warn);
    try
    {
        const auto procedure = read_procedure(input);
        const auto states = fly(procedure);
        if (!check_database.empty() && vehicle.empty())
            input_error("--check maps the states through the vehicle's rotor: give --vehicle.");
        const LengthUnit unit = unit_name == "m" ? meters : feet;
        write_procedure_track(output, states, unit);
        std::printf("Wrote %s: %zu states over %.1f s.\n", output.c_str(), states.size(),
                states.empty() ? 0.0 : states.back().t);
        if (!check_database.empty())
            check_envelope(output, unit, check_database, vehicle, max_load_factor, qsam, accept_vehicle_mismatch);
    }
    catch (const std::exception& problem)
    {
        report_again(problem, "Stopping: ");
        return 2;
    }
    return 0;
}
