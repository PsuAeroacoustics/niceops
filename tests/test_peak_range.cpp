// Ported from the unittest block in niceops.d (peak_range).

#include <catch2/catch_test_macros.hpp>

#include <vector>

#include "util/sel.hpp"

using namespace niceops;

TEST_CASE("peak_range", "[sel]")
{
    // Single peak
    std::vector<double> levels = {0.0, 10.0, 20.0, 25.0, 20.0, 10.0, 0.0};
    auto index = peak_range(levels);
    REQUIRE(index[0] == 2);
    REQUIRE(index[1] == 4);

    // Single peak, exactly at the threshold level
    levels = {0.0, 15.0, 20.0, 25.0, 20.0, 15.0, 0.0};
    index = peak_range(levels);
    REQUIRE(index[0] == 1);
    REQUIRE(index[1] == levels.size() - 2);

    // Soft peak, no 10 dB down point
    levels = {20.0, 20.0, 20.0, 25.0, 20.0, 20.0, 20.0};
    index = peak_range(levels);
    REQUIRE(index[0] == 0);
    REQUIRE(index[1] == levels.size() - 1);

    // Two peaks
    levels = {20.0, 10.0, 20.0, 25.0, 20.0, 10.0, 20.0};
    index = peak_range(levels);
    REQUIRE(index[0] == 2);
    REQUIRE(index[1] == 4);
}
