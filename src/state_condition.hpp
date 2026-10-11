#pragma once

// The operating condition a trajectory state looks a database up at: its
// advance ratio, its flight path angle (Q-SAM's, against the tip-path plane,
// or the plain one) and its thrust coefficient, held at --max_load_factor on a
// fixed_load_factor database.  One function, so that niceops and the
// procedure tool's envelope check (niceops_procedure --check) cannot map a
// state differently.

#include "condition_selection.hpp"

namespace niceops {

class Database;
struct State;

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

} // namespace niceops
