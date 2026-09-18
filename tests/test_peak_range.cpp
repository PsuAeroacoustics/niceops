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


TEST_CASE("is_night applies the 22:00-07:00 local window", "[sel]")
{
    const double noon = 12.0 * 3600.0;
    const double midnight = 0.0;
    const double six_am = 6.0 * 3600.0;
    const double ten_pm = 22.0 * 3600.0;
    const double nine_pm = 21.0 * 3600.0;

    REQUIRE(!is_night(noon, 0.0));
    REQUIRE(is_night(midnight, 0.0));
    REQUIRE(is_night(six_am, 0.0));
    REQUIRE(is_night(ten_pm, 0.0)); // the boundary belongs to night
    REQUIRE(!is_night(nine_pm, 0.0));
    REQUIRE(!is_night(7.0 * 3600.0, 0.0)); // ...and 07:00 belongs to day

    // The offset shifts the window, wrapping in both directions.
    REQUIRE(is_night(noon, 11.0));  // 23:00 local
    REQUIRE(!is_night(midnight, 9.0)); // 09:00 local
    REQUIRE(is_night(noon, -13.0)); // 23:00 local the previous day
}

TEST_CASE("day_night_level penalises night operations by 10 dB", "[sel]")
{
    // A history with a clear peak, so the 10-dB-down window is well defined.
    std::vector<double> levels;
    for (int i = 0; i < 41; i++)
        levels.push_back(80.0 - std::fabs(static_cast<double>(i - 20)) * 0.4);

    const double T = 0.02;
    const double sel = sound_exposure_level(levels, T);

    std::vector<double> scratch;
    std::vector<char> all_day(levels.size(), 0);
    std::vector<char> all_night(levels.size(), 1);

    // Daytime: the penalty is zero, so this is just the sound exposure level.
    REQUIRE(std::fabs(day_night_level(levels, T, all_day, scratch) - sel) < 1e-12);
    // Night: every sample gains 10 dB, and so does the integral -- the window
    // is unchanged because a constant shift moves the peak and the threshold
    // together.
    REQUIRE(std::fabs(day_night_level(levels, T, all_night, scratch) - (sel + 10.0)) < 1e-12);
}
