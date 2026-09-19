// Repair of sphere azimuth labels in pre-versioned databases.
//
// Source spheres cover the lower half of the roll circle, phi in [-90, 90] with
// phi = 0 straight down.  The upper half is the reflection through the
// horizontal plane, phi -> 180 - phi.  Builders before database_version 1
// completed the sphere with `-phi` instead, which is a left/right flip *within*
// the lower half: because the source range is symmetric about zero it maps
// [-90, 90] onto itself, so the mirrored half wears lower-hemisphere labels and
// the upper surface is missing.  Only the labels are wrong, so this is
// recoverable: phi_correct = wrap(180 + phi_stored).

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "sphere_database.hpp"

using namespace niceops;

namespace {

// A sphere as a pre-versioned builder wrote it: [source, -source].
std::vector<double> broken_sphere()
{
    std::vector<double> source;
    for (int p = -90; p <= 90; p += 10)
    {
        source.push_back(static_cast<double>(p));
    }
    std::vector<double> all = source;
    for (auto phi : source)
    {
        all.push_back(-phi);
    }
    return all;
}

} // namespace

TEST_CASE("legacy hover directivity is repaired", "[legacy]")
{
    auto phis = broken_sphere();
    const std::size_t half = phis.size() / 2;

    REQUIRE(repair_legacy_hover_directivity(phis));

    // The source half is untouched.
    for (std::size_t c = 0; c < half; c++)
    {
        REQUIRE(phis[c] == -90.0 + 10.0 * static_cast<double>(c));
    }
    // The mirrored half now carries the reflection of its source row.
    for (std::size_t c = 0; c < half; c++)
    {
        const double source = -90.0 + 10.0 * static_cast<double>(c);
        double expected = std::fmod(180.0 - source + 180.0, 360.0);
        if (expected < 0.0)
        {
            expected += 360.0;
        }
        expected -= 180.0;
        REQUIRE(std::fabs(phis[half + c] - expected) < 1e-9);
    }
}

TEST_CASE("repair completes the azimuth circle", "[legacy]")
{
    auto phis = broken_sphere();
    REQUIRE(repair_legacy_hover_directivity(phis));

    std::vector<double> sorted = phis;
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());

    // 19 source azimuths + 19 reflections, with +/-90 shared by both halves.
    REQUIRE(sorted.size() == 36);
    REQUIRE(sorted.front() == -180.0);
    REQUIRE(sorted.back() == 170.0);
}

TEST_CASE("the horizontal plane reflects onto itself", "[legacy]")
{
    // +/-90 lie in the plane shared by both halves, so they must come back
    // unchanged rather than being moved somewhere new.
    std::vector<double> phis = {-90.0, 90.0, 90.0, -90.0};
    REQUIRE(repair_legacy_hover_directivity(phis));
    REQUIRE(std::fabs(phis[2] - (-90.0)) < 1e-9);
    REQUIRE(std::fabs(phis[3] - 90.0) < 1e-9);
}

TEST_CASE("a correctly completed sphere is left alone", "[legacy]")
{
    // Spans the full circle, so it already has an upper surface.
    std::vector<double> phis;
    for (int p = -180; p <= 170; p += 10)
    {
        phis.push_back(static_cast<double>(p));
    }
    const auto before = phis;
    REQUIRE_FALSE(repair_legacy_hover_directivity(phis));
    REQUIRE(phis == before);
}

TEST_CASE("a sphere that does not match the broken shape is left alone", "[legacy]")
{
    // Confined to [-90, 90] but the halves are not mirrored, so the structural
    // assumption does not hold and the data must not be rewritten.
    std::vector<double> phis = {-90.0, -80.0, -70.0, -60.0};
    const auto before = phis;
    REQUIRE_FALSE(repair_legacy_hover_directivity(phis));
    REQUIRE(phis == before);

    // Odd channel count cannot be [source, mirror].
    std::vector<double> odd = {-90.0, 0.0, 90.0};
    const auto odd_before = odd;
    REQUIRE_FALSE(repair_legacy_hover_directivity(odd));
    REQUIRE(odd == odd_before);
}
