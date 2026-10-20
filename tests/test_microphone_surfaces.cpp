#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "building_input.hpp"
#include "buildings.hpp"
#include "ground_effect.hpp"
#include "point_paths.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT

// A source or a microphone may stand on a surface of the district -- a
// microphone flush on a roof or a wall, a loudspeaker on a roof -- as well as
// on the floor or in the air.  On a wall or a level roof the surface reflects
// at the point: every ray arriving there has a partner off it at the point,
// of one length (PointRay::source_surface).

namespace {

Footprint box(double x0, double y0, double x1, double y1, double top)
{
    Footprint f;
    f.x = {x0, x1, x1, x0};
    f.y = {y0, y0, y1, y1};
    f.top = top;
    return f;
}

const std::vector<double> bands {125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0};
const double c = 343.2;

std::vector<double> levels(const Buildings& city, const Vec3& source, const Vec3& microphone, double rho = 1.0,
        const FloorModel& floor = FloorModel(1.0))
{
    const PointSourceField field(city, source, microphone);
    std::vector<double> out(bands.size());
    field.band_levels(std::vector<double>(city.facades().size() * bands.size(), rho), floor, bands, c, out.data());
    return out;
}

} // namespace

// What a point lies on: a wall, a level roof (a flat one, or a mansard's deck),
// a slope, an edge (a roof's edge, a corner, a ridge), or nothing -- the air
// and the floor alike.
TEST_CASE("a point lies on a wall, a roof, a slope or an edge", "[microphone_surfaces]")
{
    Footprint mansard = box(30.0, 0.0, 42.0, 20.0, 9.0);
    mansard.roof = {{1.0, 0.0, -24.0}, {0.4, 0.0, -5.1}, {-0.4, 0.0, 28.5}, {-1.0, 0.0, 54.0}, {0.0, 0.0, 9.0}};
    const Buildings city({box(0.0, 0.0, 20.0, 20.0, 10.0), mansard});
    using Kind = Buildings::Support::Kind;
    REQUIRE(city.support(Vec3(10.0, 10.0, 10.0)).kind == Kind::roof);
    REQUIRE(city.support(Vec3(10.0, 10.0, 10.0 + 5e-7)).kind == Kind::roof);
    REQUIRE(city.support(Vec3(10.0, 10.0, 10.0 + 5e-6)).kind == Kind::none);
    REQUIRE(city.support(Vec3(20.0, 10.0, 5.0)).kind == Kind::facade);
    REQUIRE(city.support(Vec3(20.0, 10.0, 0.0)).kind == Kind::facade);
    REQUIRE(city.support(Vec3(20.0, 10.0, 10.0)).kind == Kind::edge);
    REQUIRE(city.support(Vec3(20.0, 20.0, 5.0)).kind == Kind::edge);
    REQUIRE(city.support(Vec3(25.0, 10.0, 0.0)).kind == Kind::none);
    REQUIRE(city.support(Vec3(25.0, 10.0, 3.0)).kind == Kind::none);
    // The mansard: its lower slope, its deck, and the break between its slopes.
    REQUIRE(city.support(Vec3(31.0, 10.0, 7.0)).kind == Kind::slope);
    REQUIRE(city.support(Vec3(36.0, 10.0, 9.0)).kind == Kind::roof);
    REQUIRE(city.support(Vec3(31.5, 10.0, 7.5)).kind == Kind::edge);
}

// A point on a wall or a level roof hears what one a hair off it hears, where
// the partners are traced: over the roof's edge, round the wall's corner, in
// the building's shadow, overhead; and a source on a roof is heard as one a
// hair above it.  Rigid surfaces: a diffracted ray's wedge holds its face's
// reflection rigidly, and so does the end on that face.  Within 0.1 dB (0.062
// measured, at 125 Hz in the building's shadow: the backscatter off an edge
// back along its face, which a hair off the face carries and the end on the
// face leaves out); before, a microphone on a roof had no doubling, and none
// at all behind the roof's edge, and one on a wall none either.
TEST_CASE("an end on a roof or a wall is one a hair off it", "[microphone_surfaces]")
{
    const Buildings city({box(0.0, 0.0, 20.0, 20.0, 10.0)});
    struct Case
    {
        Vec3 source, on, off;
    };
    const double hair = 1e-5;
    const std::vector<Case> cases {
        // On the roof: overhead, behind its edge, obliquely behind it.
        {{10.0, 10.0, 60.0}, {10.0, 10.0, 10.0}, {10.0, 10.0, 10.0 + hair}},
        {{-30.0, 10.0, 6.0}, {5.0, 10.0, 10.0}, {5.0, 10.0, 10.0 + hair}},
        {{-30.0, 25.0, 6.0}, {5.0, 10.0, 10.0}, {5.0, 10.0, 10.0 + hair}},
        {{10.0, -30.0, 3.0}, {5.0, 10.0, 10.0}, {5.0, 10.0, 10.0 + hair}},
        // On a wall: facing it, round its corner, behind the building, oblique.
        {{60.0, 10.0, 5.0}, {20.0, 10.0, 5.0}, {20.0 + hair, 10.0, 5.0}},
        {{60.0, 40.0, 5.0}, {20.0, 10.0, 5.0}, {20.0 + hair, 10.0, 5.0}},
        {{-30.0, 10.0, 4.0}, {20.0, 10.0, 5.0}, {20.0 + hair, 10.0, 5.0}},
        {{10.0, -30.0, 12.0}, {20.0, 10.0, 5.0}, {20.0 + hair, 10.0, 5.0}},
    };
    double worst = 0.0;
    for (const Case& k : cases)
    {
        INFO("source " << k.source.x << ", " << k.source.y << ", " << k.source.z << "; microphone " << k.on.x << ", "
                       << k.on.y << ", " << k.on.z);
        const auto on = levels(city, k.source, k.on), off = levels(city, k.source, k.off);
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            REQUIRE(std::isfinite(on[b]));
            worst = std::max(worst, std::abs(on[b] - off[b]));
        }
    }
    INFO("worst " << worst << " dB");
    REQUIRE(worst < 0.1);
    // A source on the roof, heard low behind the building.
    const Vec3 low(40.0, 10.0, 1.5);
    const auto on = levels(city, Vec3(15.0, 10.0, 10.0), low), off = levels(city, Vec3(15.0, 10.0, 10.0 + hair), low);
    for (std::size_t b = 0; b < bands.size(); b++)
        REQUIRE(on[b] == Catch::Approx(off[b]).margin(0.01));
}

// An absorbing roof reflects at a microphone on it as it does a hair above:
// overhead, (1 + sqrt(rho))^2 the free field's energy, where the ray off the
// roof is traced above it.
TEST_CASE("an absorbing roof doubles a microphone on it by its own reflection", "[microphone_surfaces]")
{
    const Buildings city({box(0.0, 0.0, 200.0, 200.0, 10.0)});
    const Vec3 source(100.0, 100.0, 60.0), on(100.0, 100.0, 10.0), off(100.0, 100.0, 10.0 + 1e-5);
    const double rho = 0.6;
    const auto at = levels(city, source, on, rho), above = levels(city, source, off, rho);
    const double expected = -20.0 * std::log10(50.0) + 20.0 * std::log10(1.0 + std::sqrt(rho));
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        REQUIRE(at[b] == Catch::Approx(above[b]).margin(0.01));
        REQUIRE(at[b] == Catch::Approx(expected).margin(0.05));
    }
}

// A plate on a roof: each ray takes P_d + R P_r, R the roof's amplitude
// reflection, where the floor's plate takes P_d and its partner off the floor
// Q P_r.  With P_d and P_r the same at every angle the level moves by
// 20 lg |P_d + R P_r| / (1 + R) from the bare roof's.
TEST_CASE("a plate on a roof takes the roof's reflection", "[microphone_surfaces]")
{
    const Buildings city({box(0.0, 0.0, 200.0, 200.0, 10.0)});
    const Vec3 source(100.0, 100.0, 60.0), microphone(110.0, 95.0, 10.0);
    const double rho = 0.81, r = std::sqrt(rho);
    std::vector<double> middles;
    for (double band : bands)
        for (const SubBand& s : sub_bands(band, c))
            middles.push_back(0.5 * (s.g1 + s.g2));
    FloorModel ground;
    ground.admittance = ground_impedance(GroundModel {}, middles, c);
    const std::complex<double> pd(0.8, 0.1), pr(0.6, -0.2);
    FloorModel plate = ground;
    plate.plate = [&](double, double, std::complex<double>* direct, std::complex<double>* image) {
        std::fill_n(direct, middles.size(), pd);
        std::fill_n(image, middles.size(), pr);
    };
    plate.surface.assign(middles.size(), r);
    const auto bare = levels(city, source, microphone, rho, ground);
    const auto plated = levels(city, source, microphone, rho, plate);
    const double expected = 20.0 * std::log10(std::abs(pd + r * pr) / (1.0 + r));
    for (std::size_t b = 0; b < bands.size(); b++)
        REQUIRE(plated[b] - bare[b] == Catch::Approx(expected).margin(1e-9));
}

// Where no one surface's reflection stands for it -- a pitched roof's slope,
// which reflects nothing in the model, or an edge -- an end is refused, by
// the field and by the tools' placement check; above the slope, below the
// ridge, it is in the air and fine.
TEST_CASE("an end on a slope or an edge is refused", "[microphone_surfaces]")
{
    Footprint gable = box(0.0, 0.0, 10.0, 20.0, 9.0);
    gable.roof = {{0.6, 0.0, 6.0}, {-0.6, 0.0, 12.0}};
    const Buildings city({gable});
    const std::vector<std::string> ids {"1"};
    const Vec3 source(-20.0, 10.0, 2.0);
    const Vec3 slope(2.0, 10.0, 7.2), edge(0.0, 10.0, 6.0), above(2.0, 10.0, 7.7);
    REQUIRE_THROWS_AS(PointSourceField(city, source, slope), std::invalid_argument);
    REQUIRE_THROWS_AS(PointSourceField(city, source, edge), std::invalid_argument);
    REQUIRE_NOTHROW(PointSourceField(city, source, above));
    REQUIRE_THAT(misplaced_point(city, ids, slope), Catch::Matchers::ContainsSubstring("sloping roof"));
    REQUIRE_THAT(misplaced_point(city, ids, edge), Catch::Matchers::ContainsSubstring("edge"));
    REQUIRE(misplaced_point(city, ids, above).empty());
    REQUIRE_THAT(misplaced_point(city, ids, Vec3(2.0, 10.0, 5.0)), Catch::Matchers::ContainsSubstring("inside"));
}
