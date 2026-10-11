#include "state_condition.hpp"

#include <cmath>
#include <string>

#include "sphere_database.hpp"
#include "trajectory.hpp"
#include "util/constants.hpp"
#include "util/logger.hpp"

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

void check_vehicle(const dimensional_data& vehicle, const Database& database, const std::string& vehicle_path,
        bool accept_mismatch)
{
    // mu and CT are non-dimensionalized by the vehicle JSON's tip speed, and
    // they index the database's operating conditions -- which were generated at
    // whatever tip speed the database records.  Disagree, and every lookup lands
    // at the wrong condition with nothing downstream able to notice.  Worth
    // checking because the two are easy to state in different units:
    // S-76D_M3.nod was built at 220 m/s, and 600 is that same figure read as
    // ft/s, which puts the tip at Mach 1.7.
    //
    // Weight and rotor radius set CT and mu in the same way and the database
    // records them too, but only the tip speed was ever compared: a vehicle file
    // 10 percent heavy moved the footprint 0.75 dB median with nothing said, and
    // a tip speed 4.6 percent off moved it 1.1 dB.  All three are checked,
    // on stderr.  The weight is compared at 1 percent: panam builds with
    // g = 9.82 and this program uses 9.80665, a 0.14 percent difference that is
    // not a disagreement.
    // Past 2 percent the run is refused, since the footprint would be wrong by
    // more than the grid's own error; between 1 and 2 it warns.  The refusal
    // names the same check as the warning, and --accept_vehicle_mismatch turns
    // it back into a warning.
    const auto vehicle_disagrees = [&](const char* what, const char* unit_text,
                                           double from_file, double from_database) {
        if (std::isnan(from_database) || approx_equal(from_file, from_database))
            return;
        const double relative = std::fabs(from_file - from_database) / std::fabs(from_database);
        if (relative > 0.02 && !accept_mismatch)
        {
            input_error("Vehicle ", what, " disagrees: ", vehicle_path, " gives ", from_file, " ",
                    unit_text, " but the database was built at ", from_database, " ", unit_text,
                    ", ", relative * 100.0, " percent off. mu and CT would index the wrong "
                    "conditions, so the run is refused. Give the right vehicle file, or "
                    "--accept_vehicle_mismatch to propagate it anyway.");
        }
        warn("Vehicle ", what, " disagrees: ", vehicle_path, " gives ", from_file, " ", unit_text,
                " but the database was built at ", from_database, " ", unit_text,
                ". mu and CT will index the wrong conditions.");
    };
    vehicle_disagrees("main rotor tip speed", "m/s", vehicle.tip_speed, database.tip_speed());
    vehicle_disagrees("main rotor radius", "m", vehicle.radius, database.rotor_radius());
    vehicle_disagrees("weight", "N", vehicle.weight, database.vehicle_weight());
    // The thrust coefficient the vehicle file gives at a load factor of 1,
    // against a fixed_load_factor database's reference.  It folds in the
    // density, which no check above sees.
    if (database.fixed_load_factor())
    {
        const double ct_at_one_g = vehicle.weight
                / (vehicle.rho * vehicle.tip_speed * vehicle.tip_speed * pi * vehicle.radius * vehicle.radius);
        vehicle_disagrees("thrust coefficient at 1 g", "", ct_at_one_g, database.reference_thrust_coefficient());
    }
}

} // namespace niceops
