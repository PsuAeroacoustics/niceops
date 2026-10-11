#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include "contour.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

struct Field
{
    std::vector<double> x, y, v;
};

// f(x, y) on a square grid from -extent to +extent in steps of `step`.
Field sample(double extent, double step, const std::function<double(double, double)>& f)
{
    Field field;
    for (double c = -extent; c <= extent + 1e-9; c += step)
    {
        field.x.push_back(c);
        field.y.push_back(c);
    }
    for (const double y : field.y)
        for (const double x : field.x)
            field.v.push_back(f(x, y));
    return field;
}

constexpr double pi = 3.14159265358979323846;

} // namespace

TEST_CASE("a cone's contour is a circle of the right area, counterclockwise", "[contour]")
{
    // 100 - r: the region at or above 90 dB is the disc r <= 10.
    const auto field = sample(30.0, 0.25, [](double x, double y) { return 100.0 - std::hypot(x, y); });
    const auto polygons = contour_polygons(field.x, field.y, field.v, 90.0);
    REQUIRE(polygons.size() == 1);
    CHECK(polygons[0].holes.empty());
    CHECK(ring_area(polygons[0].outer) > 0.0);
    // Linear interpolation inscribes the circle: the deficit is second order in the step.
    CHECK(std::fabs(polygon_area(polygons[0]) - pi * 100.0) < 0.05);
    for (const auto& point : polygons[0].outer)
        CHECK(std::fabs(std::hypot(point[0], point[1]) - 10.0) < 0.01);
}

TEST_CASE("a ring of high ground is one polygon with one clockwise hole", "[contour]")
{
    // -(r - 20)^2: at or above -25 is the annulus 15 <= r <= 25.
    const auto field = sample(40.0, 0.25, [](double x, double y) {
        const double r = std::hypot(x, y);
        return -(r - 20.0) * (r - 20.0);
    });
    const auto polygons = contour_polygons(field.x, field.y, field.v, -25.0);
    REQUIRE(polygons.size() == 1);
    REQUIRE(polygons[0].holes.size() == 1);
    CHECK(ring_area(polygons[0].holes[0]) < 0.0);
    CHECK(std::fabs(polygon_area(polygons[0]) - pi * (625.0 - 225.0)) < 0.5);
}

TEST_CASE("two peaks give two polygons, and a level above both gives none", "[contour]")
{
    const auto field = sample(30.0, 0.5, [](double x, double y) {
        return std::max(100.0 - std::hypot(x - 12.0, y), 95.0 - std::hypot(x + 12.0, y));
    });
    CHECK(contour_polygons(field.x, field.y, field.v, 90.0).size() == 2);
    CHECK(contour_polygons(field.x, field.y, field.v, 101.0).empty());
}

TEST_CASE("a region reaching the grid's edge is closed along it", "[contour]")
{
    // Everything at or above the level: the polygon is the grid itself.
    const auto field = sample(10.0, 1.0, [](double, double) { return 80.0; });
    const auto polygons = contour_polygons(field.x, field.y, field.v, 70.0);
    REQUIRE(polygons.size() == 1);
    CHECK(std::fabs(polygon_area(polygons[0]) - 400.0) < 1e-9);
    // A half plane: x >= 0 of the 20 x 20 grid.
    const auto half = sample(10.0, 1.0, [](double x, double) { return x; });
    const auto halves = contour_polygons(half.x, half.y, half.v, 0.0);
    REQUIRE(halves.size() == 1);
    CHECK(std::fabs(polygon_area(halves[0]) - 200.0) < 1e-9);
}

TEST_CASE("culled nodes are below every level", "[contour]")
{
    // A disc of finite levels in a sea of -inf, as a --db_cutoff footprint.
    const auto field = sample(20.0, 0.25, [](double x, double y) {
        const double r = std::hypot(x, y);
        return r <= 10.0 ? 100.0 - r : -std::numeric_limits<double>::infinity();
    });
    const auto polygons = contour_polygons(field.x, field.y, field.v, 95.0);
    REQUIRE(polygons.size() == 1);
    CHECK(std::fabs(polygon_area(polygons[0]) - pi * 25.0) < 0.1);
}

TEST_CASE("a saddle is decided by its center", "[contour]")
{
    // Two high corners on a diagonal of one cell.  Center high: one polygon
    // joining them; center low: two.
    const std::vector<double> x = {0.0, 1.0}, y = {0.0, 1.0};
    CHECK(contour_polygons(x, y, {10.0, 0.0, 0.0, 10.0}, 4.0).size() == 1); // mean 5 >= 4
    CHECK(contour_polygons(x, y, {10.0, 0.0, 0.0, 10.0}, 6.0).size() == 2); // mean 5 < 6
}

TEST_CASE("an island standing in a hole is its own polygon, and the hole stays with its outer ring",
        "[contour]")
{
    // High ground with a ring of low nodes around a high center node: an outer
    // ring with one hole, and an island inside the hole.  The hole's vertex
    // mean is the center, inside the island, which once took the hole.
    const std::vector<double> x = {0, 1, 2, 3, 4, 5, 6}, y = x;
    std::vector<double> v;
    for (int iy = 0; iy < 7; iy++)
        for (int ix = 0; ix < 7; ix++)
            v.push_back(std::max(std::abs(ix - 3), std::abs(iy - 3)) == 1 ? 0.0 : 9.0);
    const auto polygons = contour_polygons(x, y, v, 5.0);
    REQUIRE(polygons.size() == 2);
    for (const auto& polygon : polygons)
        CHECK(polygon_area(polygon) > 0.0);
    const auto& big = ring_area(polygons[0].outer) > ring_area(polygons[1].outer) ? polygons[0] : polygons[1];
    const auto& island = &big == &polygons[0] ? polygons[1] : polygons[0];
    CHECK(big.holes.size() == 1);
    CHECK(island.holes.empty());
}

TEST_CASE("a U-shaped hole is kept", "[contour]")
{
    // High ground with a U of low nodes in it: the U's vertex mean lies in
    // the high ground between its arms, inside no hole, and the hole was dropped.
    std::vector<double> x, y, v;
    for (int i = 0; i < 11; i++)
    {
        x.push_back(i);
        y.push_back(i);
    }
    for (int iy = 0; iy < 11; iy++)
        for (int ix = 0; ix < 11; ix++)
        {
            const bool arm = (ix == 3 || ix == 7) && iy >= 3 && iy <= 8;
            const bool base = iy == 3 && ix >= 3 && ix <= 7;
            v.push_back(arm || base ? 0.0 : 9.0);
        }
    const auto polygons = contour_polygons(x, y, v, 5.0);
    REQUIRE(polygons.size() == 1);
    REQUIRE(polygons[0].holes.size() == 1);
    // 100 for the grid less the U: the area is that of a field linear on each cell's triangles' worth, so check its sign and size only.
    CHECK(polygon_area(polygons[0]) < 100.0 - 5.0);
}

TEST_CASE("levels met exactly give simple rings", "[contour]")
{
    // Nodes equal to the level put crossings on the nodes; the ring once ran
    // out along the edge and back, and met itself.
    const std::vector<double> x = {0, 1, 2, 3}, y = {0, 1, 2};
    const auto polygons = contour_polygons(x, y, {0, 0, 0, 0, 0, 3, 0, 0, 0, 2, 2, 2}, 2.0);
    for (const auto& polygon : polygons)
    {
        for (const Ring* ring = &polygon.outer; ring != nullptr; ring = nullptr)
            for (std::size_t i = 0; i < ring->size(); i++)
                for (std::size_t j = i + 1; j < ring->size(); j++)
                    CHECK((*ring)[i] != (*ring)[j]);
    }
}

TEST_CASE("random integer fields give simple rings, positive polygons and placed holes", "[contour]")
{
    // Integer values put many nodes exactly on the level.
    unsigned state = 12345;
    const auto next = [&state]() {
        state = state * 1103515245u + 12345u;
        return (state >> 16) & 0x7fff;
    };
    std::size_t holes = 0;
    for (int trial = 0; trial < 300; trial++)
    {
        const int nx = 3 + static_cast<int>(next() % 8), ny = 3 + static_cast<int>(next() % 8);
        std::vector<double> x, y, v;
        for (int i = 0; i < nx; i++)
            x.push_back(i);
        for (int i = 0; i < ny; i++)
            y.push_back(i);
        for (int i = 0; i < nx * ny; i++)
            v.push_back(static_cast<double>(next() % 5));
        for (const auto& polygon : contour_polygons(x, y, v, 2.0))
        {
            CHECK(polygon_area(polygon) > 0.0);
            std::vector<const Ring*> rings = {&polygon.outer};
            for (const auto& hole : polygon.holes)
            {
                rings.push_back(&hole);
                holes++;
                CHECK(std::fabs(ring_area(hole)) < ring_area(polygon.outer));
            }
            for (const Ring* ring : rings)
            {
                REQUIRE(ring->size() >= 3);
                for (std::size_t i = 0; i < ring->size(); i++)
                    for (std::size_t j = i + 1; j < ring->size(); j++)
                        CHECK((*ring)[i] != (*ring)[j]);
            }
        }
    }
    CHECK(holes > 0);
}
