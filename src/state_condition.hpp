#pragma once

// The operating condition a trajectory state looks a database up at: its
// advance ratio, its flight path angle (Q-SAM's, against the tip-path plane,
// or the plain one) and its thrust coefficient, held at --max_load_factor on a
// fixed_load_factor database.  One function, so that niceops and the
// procedure tool's envelope check (niceops_procedure --check) cannot map a
// state differently.

#include <string>

#include "condition_selection.hpp"

namespace niceops {

class Database;
struct State;
struct dimensional_data;

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

StateCondition state_condition(const State& state, const Database& database, bool qsam,
        double max_load_factor, const BlendOptions& blend);

// The vehicle file's tip speed, rotor radius, weight and (on a
// fixed_load_factor database) C_T at 1 g against what the database records:
// warns past 1 percent, and past 2 throws input_error unless accept_mismatch
// (--accept_vehicle_mismatch).  mu and C_T index the database's conditions
// through the vehicle file, so a wrong one maps every state wrongly.
void check_vehicle(const dimensional_data& vehicle, const Database& database, const std::string& vehicle_path,
        bool accept_mismatch);

} // namespace niceops
