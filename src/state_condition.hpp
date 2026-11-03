#pragma once

// The operating condition a trajectory state looks a database up at: its
// advance ratio, its flight path angle (Q-SAM's, against the tip-path plane,
// or the plain one) and its thrust coefficient, held at --max_load_factor on a
// fixed_load_factor database.  One function, so that niceops and the
// procedure tool's envelope check (niceops_procedure --check) cannot map a
// state differently.

#include <optional>
#include <string>

#include "condition_selection.hpp"

namespace niceops {

class Database;
struct State;
struct dimensional_data;

// Which air density C_T is formed with (--condition_density): the vehicle
// file's, which every database has been built and run with, or the run
// atmosphere's at each state's height.  A database records which its C_T axis
// used (root attribute thrust_coefficient_density, vehicle or run).
enum class ConditionDensity
{
    vehicle,
    atmosphere,
};

// "vehicle" or "atmosphere"; nullopt for anything else.
std::optional<ConditionDensity> parse_condition_density(const std::string& name);
const char* to_string(ConditionDensity density);

struct StateCondition
{
    double advance_ratio = 0.0;
    double flight_path_angle = 0.0;
    double thrust_coefficient = 0.0;
    // On a fixed_load_factor database: C_T held at max_load_factor times the
    // reference, and (two-sided) below the --min_load_factor floor, which the
    // blend holds.
    bool load_capped = false;
    bool load_floored = false;
};

// density_ratio is the vehicle file's density over the one C_T is to be
// formed with (--condition_density atmosphere: the run's air at the state); the
// trajectory's C_T, formed with the vehicle file's, is scaled by it.  1, the
// default, leaves C_T untouched.
StateCondition state_condition(const State& state, const Database& database, bool qsam,
        double max_load_factor, const BlendOptions& blend, double density_ratio = 1.0);

// The ways a state's lookup is not its own condition, as bits: niceops's
// --provenance_file shares and niceops_procedure --check's report.
enum ConditionFlag : unsigned
{
    advance_ratio_clamped = 1u << 0,
    flight_path_angle_clamped = 1u << 1,
    thrust_coefficient_clamped = 1u << 2,
    load_held_at_max = 1u << 3, // --max_load_factor
    load_held_at_min = 1u << 4, // --min_load_factor
};
inline constexpr int condition_flag_count = 5;

// A state's flags, from its condition and the blend the database gave it
// (Database::Blend).
template <typename Blend>
unsigned condition_flags(const StateCondition& condition, const Blend& blend)
{
    return (blend.advance_ratio_clamped ? advance_ratio_clamped : 0u)
            | (blend.flight_path_angle_clamped ? flight_path_angle_clamped : 0u)
            | (blend.thrust_coefficient_clamped ? thrust_coefficient_clamped : 0u)
            | (condition.load_capped ? load_held_at_max : 0u) | (condition.load_floored ? load_held_at_min : 0u);
}

// The vehicle file's tip speed, rotor radius, weight and (on a
// fixed_load_factor database) C_T at 1 g against what the database records:
// warns past 1 percent, and past 2 throws input_error unless accept_mismatch
// (--accept_vehicle_mismatch).  mu and C_T index the database's conditions
// through the vehicle file, so a wrong one maps every state wrongly.
void check_vehicle(const dimensional_data& vehicle, const Database& database, const std::string& vehicle_path,
        bool accept_mismatch);

} // namespace niceops
