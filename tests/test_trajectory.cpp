// New test (no D original): exercise the Trajectory CSV/JSON parsing path
// with the demo case.  CTest runs with the source root as the working
// directory, so the default relative paths resolve here.

#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "trajectory.hpp"

using namespace niceops;

TEST_CASE("Trajectory parses the demo case", "[trajectory]")
{
    Trajectory trajectory("tracks/S-76D/197421AC.csv", "databases/S-76D.json");
    REQUIRE(trajectory.states.size() > 0);
    for (const auto& state : trajectory.states)
    {
        REQUIRE(std::isfinite(state.CT));
        REQUIRE(std::isfinite(state.mu));
        REQUIRE(std::isfinite(state.fpa));
        REQUIRE(std::isfinite(state.fpa_e));
    }
    auto x = trajectory.column("x");
    REQUIRE(x.size() == trajectory.states.size());
}

TEST_CASE("Hover states do not produce NaN flight path angles", "[trajectory]")
{
    // A stationary aircraft makes v zero, and the flight path angles used to be
    // asin of a division by it -- NaN, which then travelled through mu into the
    // sphere lookup and out into every level.  The database carries hover
    // spheres at mu = 0, so this regime is meant to be reachable.
    dimensional_data vehicle;
    vehicle.tip_speed = 200.0;
    vehicle.radius = 6.7;
    vehicle.weight = 5000.0 * 9.81;
    vehicle.rho = 1.225;

    State hover;
    hover.vx = 0.0;
    hover.vy = 0.0;
    hover.vz = 0.0;
    hover.nx = 0.0;
    hover.ny = 0.0;
    hover.nz = 1.0;
    derive_flight_parameters(hover, vehicle);
    REQUIRE(std::isfinite(hover.fpa));
    REQUIRE(std::isfinite(hover.fpa_e));
    REQUIRE(hover.fpa == 0.0);
    REQUIRE(hover.fpa_e == 0.0);
    REQUIRE(hover.mu == 0.0);
    REQUIRE(std::isfinite(hover.CT));

    // Forward flight is untouched: a climb along +x with the rotor level has
    // its flight path angle set by the vertical component alone.
    State climbing;
    climbing.vx = 50.0;
    climbing.vy = 0.0;
    climbing.vz = 50.0;
    climbing.nx = 0.0;
    climbing.ny = 0.0;
    climbing.nz = 1.0;
    derive_flight_parameters(climbing, vehicle);
    REQUIRE(std::fabs(climbing.fpa - 45.0) < 1e-9);
}
