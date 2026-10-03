// Ported from the unittest block in niceops.d (peak_range).

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
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

    // Two peaks: the interval runs from the first to the last sample within
    // 10 dB of the peak, dips included (14 CFR 36 A36.4.5).  The D original's
    // walk out from the peak stopped at the dips and gave [2, 4].
    levels = {20.0, 10.0, 20.0, 25.0, 20.0, 10.0, 20.0};
    index = peak_range(levels);
    REQUIRE(index[0] == 0);
    REQUIRE(index[1] == levels.size() - 1);

    // Samples outside the interval are not counted however the dips fall.
    levels = {0.0, 16.0, 0.0, 25.0, 0.0, 0.0, 0.0};
    index = peak_range(levels);
    REQUIRE(index[0] == 1);
    REQUIRE(index[1] == 3);
}

TEST_CASE("sound_exposure_level integrates in reception time", "[sel]")
{
    // A steady level for 11 samples, from a source receding at a constant
    // rate: each emission step arrives stretched by 1 + v/c.
    const double T = 0.5, c = 340.0, v = 34.0;
    std::vector<double> levels(11, 80.0), ranges;
    for (int i = 0; i < 11; i++)
        ranges.push_back(100.0 + v * T * i);
    const double emission = sound_exposure_level(levels, T);
    const double reception = sound_exposure_level(levels, T, ranges, c);
    REQUIRE(std::fabs(emission - (80.0 + 10.0 * std::log10(11 * T))) < 1e-9);
    REQUIRE(std::fabs(reception - emission - 10.0 * std::log10(1.0 + v / c)) < 1e-9);

    // Approaching, compressed by 1 - v/c.
    std::reverse(ranges.begin(), ranges.end());
    REQUIRE(std::fabs(sound_exposure_level(levels, T, ranges, c) - emission
                    - 10.0 * std::log10(1.0 - v / c)) < 1e-9);

    // Mismatched or absent ranges fall back to emission time.
    REQUIRE(sound_exposure_level(levels, T, {1.0, 2.0}, c) == emission);
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

TEST_CASE("day_night_level penalizes night operations by 10 dB", "[sel]")
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

TEST_CASE("an ambient level is added in energy", "[sel][ambient]")
{
    // Aircraft and ambient noise are uncorrelated, so their energies add.  The D
    // original summed pressure amplitudes, which put two equal levels 6.02 dB up
    // and a level 10 dB below the ambient 2.4 dB up, where energy gives 3.01 and 0.41.
    const double ambient = 60.0;
    const double energy = std::pow(10.0, 0.1 * ambient);
    REQUIRE(level_with_ambient(60.0, energy) - 60.0 == Catch::Approx(10.0 * std::log10(2.0)).margin(1e-12));
    REQUIRE(level_with_ambient(50.0, energy) - 60.0 == Catch::Approx(10.0 * std::log10(1.1)).margin(1e-12));
    REQUIRE(level_with_ambient(80.0, energy) - 80.0 == Catch::Approx(10.0 * std::log10(1.01)).margin(1e-12));
    // Never below either part, and a culled state (-inf) leaves the ambient alone.
    for (double level : {-30.0, 0.0, 40.0, 59.9, 60.1, 90.0})
        REQUIRE(level_with_ambient(level, energy) > std::max(level, ambient));
    REQUIRE(level_with_ambient(-std::numeric_limits<double>::infinity(), energy)
            == Catch::Approx(ambient).margin(1e-12));
    // A negligible ambient changes nothing at the golden's tolerance.
    REQUIRE(level_with_ambient(75.0, std::pow(10.0, -20.0)) == Catch::Approx(75.0).epsilon(1e-12));
}
