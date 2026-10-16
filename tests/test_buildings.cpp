#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
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

TEST_CASE("the plan grid finds what trying everything finds", "[buildings]")
{
    // A random district: rectangles and skewed quadrilaterals of random
    // heights, some touching.
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<Footprint> footprints;
    for (int i = 0; i < 60; i++)
    {
        const double cx = 400.0 * unit(rng), cy = 400.0 * unit(rng);
        const double a = 5.0 + 30.0 * unit(rng), b = 5.0 + 30.0 * unit(rng), turn = unit(rng) * pi;
        Footprint f;
        for (int k = 0; k < 4; k++)
        {
            const double px = (k == 1 || k == 2 ? a : -a) * (0.7 + 0.3 * unit(rng));
            const double py = (k >= 2 ? b : -b) * (0.7 + 0.3 * unit(rng));
            f.x.push_back(cx + px * std::cos(turn) - py * std::sin(turn));
            f.y.push_back(cy + px * std::sin(turn) + py * std::cos(turn));
        }
        f.top = 5.0 + 45.0 * unit(rng);
        footprints.push_back(f);
    }
    const Buildings city(footprints);
    std::vector<Buildings> alone;
    for (const auto& f : city.footprints())
        alone.emplace_back(std::vector<Footprint> {f});
    auto random_point = [&](double z_low, double z_high) {
        return Vec3(-50.0 + 500.0 * unit(rng), -50.0 + 500.0 * unit(rng), z_low + (z_high - z_low) * unit(rng));
    };
    std::size_t blocked = 0, hits = 0;
    for (int trial = 0; trial < 4000; trial++)
    {
        const bool below = trial % 3 == 0;
        const Vec3 a = random_point(below ? -40.0 : 0.0, 60.0), b = random_point(0.0, 60.0);
        bool any = false;
        for (const auto& one : alone)
            any = any || one.blocked(a, b, below);
        REQUIRE(city.blocked(a, b, below) == any);
        blocked += any;

        // A rising ray, against every facade tried in turn.
        Vec3 dir = Vec3(unit(rng) - 0.5, unit(rng) - 0.5, 0.3 * unit(rng));
        dir.normalize();
        Vec3 hit;
        const std::size_t found = city.first_hit(a, dir, below, hit);
        std::size_t best = city.facades().size();
        double best_t = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < city.facades().size(); i++)
        {
            const Facade& f = city.facades()[i];
            const double side = facade_side(f, a);
            const double length = std::hypot(f.x1 - f.x0, f.y1 - f.y0);
            const double nx = -(f.y1 - f.y0) / length, ny = (f.x1 - f.x0) / length;
            const double approach = dir.x * nx + dir.y * ny;
            if (side <= Buildings::contact || approach >= 0.0)
                continue;
            const double t = side / -approach;
            const Vec3 p = a + t * dir;
            const double along = ((p.x - f.x0) * (f.x1 - f.x0) + (p.y - f.y0) * (f.y1 - f.y0)) / length;
            const double bottom = below ? -f.top : f.bottom;
            if (t < best_t && along >= -1e-9 && along <= length + 1e-9 && p.z >= bottom - 1e-9 && p.z <= f.top + 1e-9)
            {
                best = i;
                best_t = t;
            }
        }
        REQUIRE(found == best);
        hits += found < city.facades().size();
    }
    // Both kinds of answer were exercised.
    REQUIRE(blocked > 400);
    REQUIRE(blocked < 3600);
    REQUIRE(hits > 400);
}
