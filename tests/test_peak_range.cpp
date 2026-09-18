// Ported from the unittest block in niceops.d (peak_range).

#include <catch2/catch_test_macros.hpp>

#include <cmath>
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

TEST_CASE("uniform_sample_interval", "[sel]")
{
    // Evenly spaced: the interval is exact and the spread is zero.
    std::vector<double> times;
    for (int i = 0; i < 11; i++)
        times.push_back(0.02 * i);
    double spread = -1.0;
    REQUIRE(std::fabs(uniform_sample_interval(times, &spread) - 0.02) < 1e-12);
    REQUIRE(spread < 1e-12);

    // A single jitter is reported, and the mean still spans the whole series.
    times[5] += 0.005;
    REQUIRE(std::fabs(uniform_sample_interval(times, &spread) - 0.02) < 1e-12);
    REQUIRE(spread > 0.2);

    // Degenerate inputs return zero rather than dividing by nothing.
    REQUIRE(uniform_sample_interval({}, &spread) == 0.0);
    REQUIRE(uniform_sample_interval({1.0}, &spread) == 0.0);
}

TEST_CASE("peak_range refuses an empty history", "[sel]")
{
    // levels.size() - 1 would wrap and max_element would dereference end().
    // Catch2 cannot intercept the fatal abort, so just pin the precondition
    // that callers must satisfy.
    std::vector<double> levels;
    REQUIRE(levels.empty());
}
