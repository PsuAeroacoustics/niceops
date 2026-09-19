// A-weighting and band-limited integration.
//
// The reference values are panam's own dBAw(), extracted from
// flight_acoustics.py and evaluated directly rather than retyped, so this
// checks the port rather than checking the formula against itself.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <vector>

#include "weighting.hpp"

using namespace niceops;

namespace {

struct Reference
{
    double frequency;
    double weight;
};

void check_close(double got, double expected, double tolerance)
{
    REQUIRE(std::fabs(got - expected) <= tolerance);
}

} // namespace

TEST_CASE("A-weighting matches panam", "[weighting]")
{
    const std::vector<Reference> reference = {
        {     10.0, -7.043036976900915e+01},
        {     50.0, -3.027161818696229e+01},
        {    100.0, -1.914277845996175e+01},
        {    200.0, -1.084611992692391e+01},
        {    500.0, -3.247677543915642e+00},
        {   1000.0, -1.487616528087614e-06},
        {   2000.0, +1.201400801961097e+00},
        {   4000.0, +9.632997918558588e-01},
        {   8000.0, -1.147126570799637e+00},
        {  10000.0, -2.491786529170330e+00},
        {  20000.0, -9.346910546198494e+00},
    };

    for (const auto& r : reference)
        check_close(a_weighting(r.frequency), r.weight, 1e-12);

    // The curve's defining property: 0 dB at 1 kHz.
    check_close(a_weighting(1000.0), 0.0, 1e-5);

    // Non-positive frequencies are floored rather than allowed to go to -inf,
    // so a band centre of zero cannot poison a sum.
    REQUIRE(a_weighting(0.0) == -300.0);
    REQUIRE(a_weighting(-5.0) == -300.0);
}

TEST_CASE("band limits exclude bands outright", "[weighting]")
{
    const std::vector<double> bands = {10.0, 50.0, 500.0, 5000.0, 30000.0};

    BandLimits wide;
    REQUIRE(bands_in_range(bands, wide) == 3); // 50, 500, 5000: 10 and 30k are outside 20..20k

    BandLimits narrow;
    narrow.low = 100.0;
    narrow.high = 1000.0;
    REQUIRE(bands_in_range(bands, narrow) == 1); // 500 only

    const auto gains = a_weighting_gains(bands, narrow);
    REQUIRE(gains.size() == bands.size());
    // Excluded bands are exactly zero, not merely small: an excluded band must
    // contribute nothing whatever its level.
    REQUIRE(gains[0] == 0.0);
    REQUIRE(gains[1] == 0.0);
    REQUIRE(gains[2] > 0.0);
    REQUIRE(gains[3] == 0.0);
    REQUIRE(gains[4] == 0.0);

    // Inclusive at the edges.
    BandLimits edge;
    edge.low = 50.0;
    edge.high = 5000.0;
    REQUIRE(bands_in_range(bands, edge) == 3);
}

TEST_CASE("weighted_level integrates bands on an energy basis", "[weighting]")
{
    // A single band at 1 kHz comes back as the level plus the A-weight there.
    // That weight is not identically zero: the pole-zero form gives -1.49e-6 dB
    // at 1 kHz rather than the 0 dB the curve is defined to have, and carrying
    // it explicitly is what keeps the tolerances below meaningful.
    const double at_1k = a_weighting(1000.0);
    const std::vector<double> one_band = {1000.0};
    const auto gains = a_weighting_gains(one_band, BandLimits{});
    const double level = 87.3;
    check_close(weighted_level(&level, gains), 87.3 + at_1k, 1e-12);

    // Two equal bands at 1 kHz sum to +3.01 dB, the energy sum rather than the
    // arithmetic one.
    const std::vector<double> two_bands = {1000.0, 1000.0};
    const auto two_gains = a_weighting_gains(two_bands, BandLimits{});
    const double levels[2] = {87.3, 87.3};
    check_close(weighted_level(levels, two_gains), 87.3 + at_1k + 10.0 * std::log10(2.0), 1e-12);

    // A band well below the A-weighting's reach barely moves the total: 50 Hz
    // is 30 dB down, so an equal-level band there adds under 0.005 dB.
    const std::vector<double> split = {1000.0, 50.0};
    const auto split_gains = a_weighting_gains(split, BandLimits{});
    const double split_levels[2] = {87.3, 87.3};
    const double got = weighted_level(split_levels, split_gains);
    REQUIRE(got > 87.3);
    REQUIRE(got < 87.31);

    // Everything excluded is silence, and silence is -inf rather than 0 dB --
    // 0 dB would read as a real (quiet) level and survive an energy sum.
    BandLimits nothing;
    nothing.low = 1e5;
    nothing.high = 2e5;
    const auto empty_gains = a_weighting_gains(split, nothing);
    REQUIRE(weighted_level(split_levels, empty_gains)
            == -std::numeric_limits<double>::infinity());

    REQUIRE(weighted_level(nullptr, {}) == -std::numeric_limits<double>::infinity());
}
