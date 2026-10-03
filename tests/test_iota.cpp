#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>

#include "util/constants.hpp"
#include "util/iota.hpp"

using namespace niceops; // NOLINT

TEST_CASE("iota_exclusive excludes its end value", "[iota]")
{
    // D's iota(-1500, 1500, 25) is 120 values ending at 1475.
    const auto values = iota_exclusive(-1500.0, 1500.0, 25.0);
    REQUIRE(values.size() == 120);
    REQUIRE(values.front() == -1500.0);
    REQUIRE(values.back() == 1475.0);
    // A span that is not a whole number of steps keeps its partial last step.
    REQUIRE(iota_exclusive(0.0, 10.0, 4.0).size() == 3);
    REQUIRE(iota_exclusive(0.0, 10.0, 3.0).size() == 4);
    REQUIRE(iota_exclusive(5.0, 5.0, 1.0).empty());
    REQUIRE(iota_exclusive(5.0, 1.0, 1.0).empty());
}

TEST_CASE("a grid given in feet has the count of the same grid in meters", "[iota]")
{
    // run() converts each bound and the step to meters separately.  For many
    // exact multiples the quotient then lands a ULP above the integer, and a
    // bare ceil() gave one row and column too many: --grid -1500 1500 --gridstep
    // 100 was 31 values in feet and 30 in meters.
    std::size_t checked = 0;
    for (int step = 1; step <= 400; step++)
    {
        for (int cells : {1, 2, 3, 10, 30, 41, 60, 120})
        {
            for (int begin : {-6000, -4921, -3000, -1500, -100, 0, 250})
            {
                const double span = static_cast<double>(step) * cells;
                const double end = begin + span;
                const auto meters = iota_exclusive(begin * ft2m, end * ft2m, step * ft2m);
                REQUIRE(meters.size() == static_cast<std::size_t>(cells));
                checked++;
            }
        }
    }
    REQUIRE(checked > 20000);
}
