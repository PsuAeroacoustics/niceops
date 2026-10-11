#include "state_condition.hpp"

#include "sphere_database.hpp"
#include "trajectory.hpp"

namespace niceops {

StateCondition state_condition(const State& state, const Database& database, bool qsam,
        double max_load_factor, const BlendOptions& blend)
{
    StateCondition condition;
    condition.advance_ratio = state.mu;
    condition.flight_path_angle = qsam ? state.fpa_e : state.fpa;
    // On a fixed_load_factor database the load-factor law is the analytic
    // correction 20 lg(C_T / C_T,ref), an extrapolation above the load
    // factors its source supports.  It is held at --max_load_factor, and
    // two-sided at --min_load_factor below (blend_at() applies that floor).
    condition.thrust_coefficient = state.CT;
    if (database.fixed_load_factor())
    {
        const double ceiling = max_load_factor * database.reference_thrust_coefficient();
        if (condition.thrust_coefficient > ceiling)
        {
            condition.thrust_coefficient = ceiling;
            condition.load_capped = true;
        }
        if (blend.load_factor_correction == LoadFactorCorrection::two_sided
                && condition.thrust_coefficient < blend.min_load_factor * database.reference_thrust_coefficient())
            condition.load_floored = true;
    }
    return condition;
}

} // namespace niceops
