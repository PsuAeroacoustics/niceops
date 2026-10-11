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
