#include "state_condition.hpp"

#include <cmath>
#include <string>

#include "atmosphere.hpp"
#include "sphere_database.hpp"
#include "trajectory.hpp"
#include "util/constants.hpp"
#include "util/logger.hpp"

namespace niceops {

std::optional<ConditionDensity> parse_condition_density(const std::string& name)
{
    if (name == "vehicle")
        return ConditionDensity::vehicle;
    if (name == "atmosphere")
        return ConditionDensity::atmosphere;
    return std::nullopt;
}

const char* to_string(ConditionDensity density)
{
    return density == ConditionDensity::atmosphere ? "atmosphere" : "vehicle";
}

StateCondition state_condition(const State& state, const Database& database, bool qsam,
        double max_load_factor, const BlendOptions& blend, double density_ratio)
{
    StateCondition condition;
    condition.advance_ratio = state.mu;
    condition.flight_path_angle = qsam ? state.fpa_e : state.fpa;
    // On a fixed_load_factor database the load-factor law is the analytic
    // correction 20 lg(C_T / C_T,ref), an extrapolation above the load
    // factors its source supports.  It is held at --max_load_factor, and
    // two-sided at --min_load_factor below (blend_at() applies that floor).
    // C_T = n W / (rho V_tip^2 pi R^2) is inverse in the density, so another
    // density is a ratio on the trajectory's.  Not multiplied at all when the
    // ratio is 1, so the vehicle file's C_T reaches the blend as it was.
    condition.thrust_coefficient = density_ratio == 1.0 ? state.CT : state.CT * density_ratio;
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

StateCondition multirotor_condition(const State& state, const Database& database, double weight_newtons,
        double density, double sound_speed)
{
    StateCondition condition;
    double load_factor = state.LF;
    if (!(load_factor >= multirotor_minimum_load_factor))
    {
        load_factor = multirotor_minimum_load_factor;
        condition.load_floored = true;
    }
    const double weight = weight_newtons * load_factor;
    const double rotors = database.rotor_count();
    const double radius = database.multirotor_rotor_radius();
    const double disk = pi * radius * radius;
    const double induced = std::sqrt(weight / (2.0 * density * rotors * disk));
    const double hover_tip_speed = std::sqrt(weight / (rotors * density * disk * database.hover_thrust_coefficient()));
    const auto u = air_velocity(state);
    const double norm = std::sqrt(state.nx * state.nx + state.ny * state.ny + state.nz * state.nz);
    const double n[3] = {state.nx / norm, state.ny / norm, state.nz / norm};
    const double axial = u[0] * n[0] + u[1] * n[1] + u[2] * n[2];
    const double in_plane[3] = {u[0] - axial * n[0], u[1] - axial * n[1], u[2] - axial * n[2]};
    condition.advance_ratio = std::sqrt(in_plane[0] * in_plane[0] + in_plane[1] * in_plane[1]
                                      + in_plane[2] * in_plane[2])
            / induced;
    condition.flight_path_angle = axial / induced;
    condition.thrust_coefficient = hover_tip_speed / sound_speed;
    condition.hover_rotor_speed = hover_tip_speed / radius;
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
    // the S-76D databases were built at 220 m/s, and 600 is that same figure read as
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
    if (database.multirotor())
    {
        // A multirotor's axes are normalized by the state's own weight and air
        // (momentum similarity), so neither the weight nor a tip speed indexes
        // its conditions; the rotor radius sets every one of them.
        vehicle_disagrees("rotor radius", "m", vehicle.radius, database.multirotor_rotor_radius());
        return;
    }
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

std::optional<dimensional_data> database_vehicle(const Database& database, const std::string& name)
{
    const double radius = database.rotor_radius(), tip_speed = database.tip_speed(),
                 weight = database.vehicle_weight(), rho = database.air_density();
    if (!(std::isfinite(radius) && std::isfinite(tip_speed) && std::isfinite(weight) && std::isfinite(rho)))
        return std::nullopt;
    // What it records has to be a vehicle, as a vehicle file's numbers do: a zero
    // or negative value would make every CT and mu infinite or NaN.
    if (!(radius > 0.0 && tip_speed > 0.0 && weight > 0.0 && rho > 0.0))
    {
        input_error("Database ", name, " records its vehicle at the root, but not one that can fly: rotor radius ",
                radius, " m, tip speed ", tip_speed, " m/s, weight ", weight, " N, air density ", rho,
                " kg/m^3. Give a vehicle file (-v).");
    }
    dimensional_data vehicle;
    vehicle.radius = radius;
    vehicle.tip_speed = tip_speed;
    vehicle.weight = weight;
    vehicle.rho = rho;
    // The temperature and pressure the database was built in go with its
    // density, as a vehicle file's Atmosphere section does.  Its humidity does
    // not: that is the sphere-build day's, and the spectral path refuses to
    // assume a humidity nobody chose for this prediction (--humidity, or a
    // vehicle file that gives one).
    if (const auto& built = database.build_atmosphere())
    {
        vehicle.atmosphere.temperature = built->temperature_kelvin;
        vehicle.atmosphere.pressure = built->pressure_kpa;
        vehicle.temperature_given = vehicle.pressure_given = true;
    }
    else
    {
        // As configure_vehicle without a pressure: the one the density implies.
        constexpr double dry_air_gas_constant = 287.05; // J/(kg K)
        vehicle.atmosphere.pressure = vehicle.rho * dry_air_gas_constant * vehicle.atmosphere.temperature / 1000.0;
    }
    // As configure_vehicle checks a file's air (a density in the wrong unit
    // shows up in the pressure it implies).
    for (const auto& problem :
            {temperature_problem(vehicle.atmosphere.temperature), pressure_problem(vehicle.atmosphere.pressure)})
    {
        if (!problem.empty())
            input_error("Database ", name, ", the vehicle it records at the root: ", problem);
    }
    // A density in the wrong unit (g/m^3, or lb/ft^3) would put every CT off
    // by orders of magnitude; air anywhere a rotorcraft flies is 0.3 to 2 kg/m^3.
    // (PANAM records the vehicle file's density, which need not be the build
    // atmosphere's, so the two are not compared.)
    if (!(vehicle.rho >= 0.3 && vehicle.rho <= 2.0))
        input_error("Database ", name, " records an air density of ", vehicle.rho,
                " kg/m^3 at the root, which is not air (0.3 to 2 kg/m^3); is it in another unit? "
                "Give a vehicle file (-v).");
    return vehicle;
}

VehicleChoice choose_vehicle(const Database& database, const std::string& database_name,
        const std::string& vehicle_file, bool vehicle_given)
{
    if (!vehicle_given)
    {
        if (auto vehicle = database_vehicle(database, database_name))
            return {*vehicle, database_name + " (its root vehicle)", true};
        if (vehicle_file.empty())
            input_error("No vehicle file was given (-v, --vehicle, or an operations file's vehicle column), and ",
                    database_name, " records no vehicle at its root: give one.");
    }
    return {configure_vehicle(vehicle_file), vehicle_file, false};
}

} // namespace niceops
