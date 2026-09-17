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
