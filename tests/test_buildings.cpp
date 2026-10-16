#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "buildings.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT

namespace {

Footprint box(double x0, double y0, double x1, double y1, double top)
{
    Footprint f;
    f.x = {x0, x1, x1, x0};
    f.y = {y0, y0, y1, y1};
    f.top = top;
    return f;
}

} // namespace

TEST_CASE("footprints are oriented, and their facades face out", "[buildings]")
{
    Footprint clockwise = box(0.0, 0.0, 10.0, 20.0, 15.0);
    std::reverse(clockwise.x.begin(), clockwise.x.end());
    std::reverse(clockwise.y.begin(), clockwise.y.end());
    const Buildings city({clockwise});
    REQUIRE(city.facades().size() == 4);
    const Vec3 outside(-5.0, 10.0, 2.0), outside_far(30.0, 10.0, 2.0);
    int facing = 0;
    for (const auto& facade : city.facades())
    {
        // Every point of the box is behind every facade.
        REQUIRE(facade_side(facade, Vec3(5.0, 10.0, 2.0)) < 0.0);
        REQUIRE(facade.bottom == 0.0);
        REQUIRE(facade.top == 15.0);
        if (facade_side(facade, outside) > 0.0 || facade_side(facade, outside_far) > 0.0)
            facing++;
    }
    REQUIRE(facing >= 2);
    // Four roof edges and four convex corners, every one a right angle.
    REQUIRE(city.edges().size() == 8);
    for (const auto& edge : city.edges())
        REQUIRE(edge.n == Catch::Approx(1.5).epsilon(1e-12));
}

TEST_CASE("an edge's angles run from its first face through the air", "[buildings]")
{
    const Buildings city({box(0.0, 0.0, 10.0, 10.0, 20.0)});
    for (const auto& edge : city.edges())
    {
        REQUIRE(edge.angle(edge.face0) == Catch::Approx(0.0).margin(1e-12));
        REQUIRE(edge.angle(edge.air) == Catch::Approx(pi / 2.0).epsilon(1e-12));
        // Straight away from the building is in the air, in the middle of it.
        Vec3 mid = edge.start + 0.5 * edge.length * edge.direction;
        const Vec3 center(5.0, 5.0, mid.z);
        Vec3 away = mid - center;
        away = away - away.dot(edge.direction) * edge.direction;
        if (!edge.vertical)
            away = away + Vec3(0.0, 0.0, away.horizontal_norm()); // up and out
        const double theta = edge.angle(away.normalized());
        REQUIRE(theta > 0.0);
        REQUIRE(theta < edge.n * pi);
        REQUIRE(theta == Catch::Approx(edge.n * pi / 2.0).epsilon(1e-9));
    }
}

TEST_CASE("an L-shaped footprint diffracts at its five convex corners only", "[buildings]")
{
    Footprint l;
    l.x = {0.0, 20.0, 20.0, 10.0, 10.0, 0.0};
    l.y = {0.0, 0.0, 10.0, 10.0, 20.0, 20.0};
    l.top = 12.0;
    const Buildings city({l});
    int corners = 0;
    for (const auto& edge : city.edges())
        if (edge.vertical)
        {
            corners++;
            REQUIRE(edge.n == Catch::Approx(1.5).epsilon(1e-12));
            REQUIRE(!(edge.start.x == 10.0 && edge.start.y == 10.0)); // the reentrant corner
        }
    REQUIRE(corners == 5);
    // A triangle's corners are sharper: n = 2 - angle / pi.
    Footprint triangle;
    triangle.x = {0.0, 10.0, 0.0};
    triangle.y = {0.0, 0.0, 10.0};
    triangle.top = 5.0;
    const Buildings t({triangle});
    std::vector<double> n;
    for (const auto& edge : t.edges())
        if (edge.vertical)
            n.push_back(edge.n);
    std::sort(n.begin(), n.end());
    REQUIRE(n.size() == 3);
    REQUIRE(n[0] == Catch::Approx(1.5).epsilon(1e-12));
    REQUIRE(n[1] == Catch::Approx(1.75).epsilon(1e-12));
    REQUIRE(n[2] == Catch::Approx(1.75).epsilon(1e-12));
}

TEST_CASE("blocking: through, over, touching and grazing", "[buildings]")
{
    const Buildings city({box(0.0, 0.0, 10.0, 10.0, 20.0)});
    // Straight through, below the roof.
    REQUIRE(city.blocked(Vec3(-5.0, 5.0, 5.0), Vec3(15.0, 5.0, 5.0)));
    // Over the roof.
    REQUIRE(!city.blocked(Vec3(-5.0, 5.0, 25.0), Vec3(15.0, 5.0, 25.0)));
    // In through a wall and out through the roof.
    REQUIRE(city.blocked(Vec3(-5.0, 5.0, 5.0), Vec3(15.0, 5.0, 35.0)));
    // Clipping a corner of the footprint by a metre.
    REQUIRE(city.blocked(Vec3(-1.0, 0.0, 5.0), Vec3(1.0, 2.0, 5.0)));
    // Ending on a wall, as a reflection does.
    REQUIRE(!city.blocked(Vec3(-5.0, 5.0, 5.0), Vec3(0.0, 5.0, 5.0)));
    // Starting on a roof edge and leaving into the air.
    REQUIRE(!city.blocked(Vec3(0.0, 5.0, 20.0), Vec3(-10.0, 5.0, 2.0)));
    // ... and into the building.
    REQUIRE(city.blocked(Vec3(0.0, 5.0, 20.0), Vec3(10.0, 5.0, 2.0)));
    // Grazing along a wall, and across the roof's plane.
    REQUIRE(!city.blocked(Vec3(0.0, -5.0, 5.0), Vec3(0.0, 15.0, 5.0)));
    REQUIRE(!city.blocked(Vec3(-5.0, 5.0, 20.0), Vec3(15.0, 5.0, 20.0)));
    // Through a vertex exactly.
    REQUIRE(!city.blocked(Vec3(-5.0, 5.0, 5.0), Vec3(5.0, -5.0, 5.0)));
    // Below the floor only in its mirror world.
    REQUIRE(!city.blocked(Vec3(-5.0, 5.0, -5.0), Vec3(15.0, 5.0, -5.0)));
    REQUIRE(city.blocked(Vec3(-5.0, 5.0, -5.0), Vec3(15.0, 5.0, -5.0), true));
    REQUIRE(!city.blocked(Vec3(-5.0, 5.0, -25.0), Vec3(15.0, 5.0, -25.0), true));
}

TEST_CASE("escaping to infinity", "[buildings]")
{
    const Buildings city({box(0.0, 0.0, 10.0, 10.0, 20.0), box(30.0, 0.0, 40.0, 10.0, 50.0)});
    const Vec3 street(20.0, 5.0, 1.5);
    // Up the street's axis, at any elevation.
    REQUIRE(city.escapes(street, Vec3(0.0, 1.0, 0.0)));
    REQUIRE(city.escapes(street, Vec3(0.0, 1.0, 1.0).normalized()));
    // Toward the lower building: blocked below the line over its roof edge.
    const double clear = std::atan2(20.0 - 1.5, 10.0);
    REQUIRE(!city.escapes(street, Vec3(-std::cos(clear - 0.01), 0.0, std::sin(clear - 0.01))));
    REQUIRE(city.escapes(street, Vec3(-std::cos(clear + 0.01), 0.0, std::sin(clear + 0.01))));
    // Straight up from the street, and from the roof.
    REQUIRE(city.escapes(street, Vec3(0.0, 0.0, 1.0)));
    REQUIRE(city.escapes(Vec3(5.0, 5.0, 20.0), Vec3(0.0, 0.0, 1.0)));
    // From the floor's image the mirrored buildings block too: the ray has
    // to climb from 1.5 m below the floor to clear the same roof.
    const Vec3 image(20.0, 5.0, -1.5);
    const double image_clear = std::atan2(20.0 + 1.5, 10.0);
    REQUIRE(!city.escapes(image, Vec3(-std::cos(image_clear - 0.01), 0.0, std::sin(image_clear - 0.01)), true));
    REQUIRE(city.escapes(image, Vec3(-std::cos(image_clear + 0.01), 0.0, std::sin(image_clear + 0.01)), true));
    REQUIRE(!city.escapes(image, Vec3(-1.0, 0.0, 0.05).normalized(), true));
}

TEST_CASE("bad footprints are refused", "[buildings]")
{
    Footprint two;
    two.x = {0.0, 1.0};
    two.y = {0.0, 1.0};
    two.top = 5.0;
    REQUIRE_THROWS_AS(Buildings({two}), std::invalid_argument);
    Footprint bow;
    bow.x = {0.0, 10.0, 0.0, 10.0};
    bow.y = {0.0, 10.0, 10.0, 0.0};
    bow.top = 5.0;
    REQUIRE_THROWS_AS(Buildings({bow}), std::invalid_argument);
    Footprint flat = box(0.0, 0.0, 1.0, 1.0, 0.0);
    REQUIRE_THROWS_AS(Buildings({flat}), std::invalid_argument);
}
