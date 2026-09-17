// Ported from the unittest blocks in sphere_database.d.  The D tests use
// std.math.approxEqual (relative epsilon 1e-6), so the port uses the C++
// approx_equal() from sphere_database.hpp to keep identical pass/fail
// semantics rather than Catch::Approx.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "sphere_database.hpp"

using namespace niceops;

TEST_CASE("periodic_bilinear_interp", "[sphere_database]")
{
    // Test bilinear interpolation
    double x1 = -1.0, x2 = 1.0, y1 = -1.0, y2 = 1.0;
    double q11 = 0.0, q22 = 1.0, q12 = 0.5, q21 = 0.5;
    auto midpoint = periodic_bilinear_interp(q11, q12, q21, q22, x1, x2, y1, y2, 0.5, 0.5);
    REQUIRE(midpoint <= 1.0);
    REQUIRE(midpoint >= -1.0);

    x1 = 1.0;
    x2 = 1.0;
    y1 = -1.0;
    y2 = 1.0;
    q11 = 0.0;
    q22 = 1.0;
    q12 = 1.0;
    q21 = 0.0;
    midpoint = periodic_bilinear_interp(q11, q12, q21, q22, x1, x2, y1, y2, 0.0, 0.0);
    REQUIRE(approx_equal(midpoint, 0.5));

    x1 = 50.0;
    x2 = -40.0;
    y1 = 85;
    y2 = -95.0;
    q11 = 0.0;
    q22 = 1.0;
    q12 = 1.0;
    q21 = 0.0;
    midpoint = periodic_bilinear_interp(q11, q12, q21, q22, x1, x2, y1, y2, 5.0, -5.0);
    REQUIRE(approx_equal(midpoint, 0.5));
}

TEST_CASE("bracket", "[sphere_database]")
{
    // iota(-90.0, 180.0, 5.0) in D is inclusive of the stop value
    std::vector<double> test_range;
    for (int i = -90; i <= 180; i += 5)
    {
        test_range.push_back(static_cast<double>(i));
    }
    double target = 10.0;
    auto brack = bracket(test_range, target);
    // A matching value in the range
    REQUIRE((test_range[brack[0]] <= target && test_range[brack[1]] >= target));
    // A bracketed value
    target = 2.345;
    brack = bracket(test_range, target);
    REQUIRE((test_range[brack[0]] <= target && test_range[brack[1]] >= target));
    // at lower limit
    target = 0.0;
    brack = bracket(test_range, target);
    REQUIRE((test_range[brack[0]] <= target && test_range[brack[1]] >= target));
    // wrap around upper limit
    target = 180.0;
    brack = bracket(test_range, target);
    REQUIRE((test_range[brack[0]] <= target && brack[1] == 0));
}

TEST_CASE("clamp", "[sphere_database]")
{
    REQUIRE(clamp01(0.5) == 0.5);
    REQUIRE(clamp01(-10.) == 0.);
    REQUIRE(clamp01(2.3) == 1.0);
}

TEST_CASE("geometry helpers", "[sphere_database]")
{
    REQUIRE(dot({1, 2, 3}, {-1, 2, 1}) == 6);
    REQUIRE(approx_equal(norm({1., 2., 3.}), std::sqrt(1. * 1 + 2. * 2 + 3. * 3)));
    // The D original tested a 4-element array; the C++ port is 3D only
    REQUIRE(approx_equal(norm(normalize({1., 2., 3.})), 1.0));

    auto x = cross({1.0, 2.0, 3.0}, {3.0, 1.0, 2.0});
    REQUIRE(approx_equal(x(0), 1.0));
    REQUIRE(approx_equal(x(1), 7.0));
    REQUIRE(approx_equal(x(2), -5.0));

    REQUIRE(acosc(0.5) == std::acos(0.5));
    REQUIRE(acosc(10.0) == std::acos(1.0));
    REQUIRE(acosc(-10.0) == std::acos(-1.0));
}

TEST_CASE("lookup_angles", "[sphere_database]")
{
    const double r = 123.456;
    Eigen::Vector3d source(0, 0, r);
    Eigen::Vector3d observer(r, 0, 0);
    Eigen::Vector3d velocity1(1, 0, 0);
    Eigen::Vector3d velocity2(0, 1, 0);
    Eigen::Vector3d velocity3(0, -1, 0);
    Eigen::Vector3d velocity4(-1, 0, 0);
    const double s2 = std::sqrt(2.);
    Eigen::Vector3d normal1(0, 0, 1);
    Eigen::Vector3d normal2(s2, 0, s2);
    Eigen::Vector3d normal3(-s2, 0, s2);
    Eigen::Vector3d normal4(0, s2, s2);
    Eigen::Vector3d normal5(0, -s2, s2);

    // Check various headings only
    auto angles = lookup_angles(source, velocity1, normal1, observer);
    REQUIRE((approx_equal(angles[0], 0.0) && approx_equal(angles[1], 45.0)));
    angles = lookup_angles(source, velocity2, normal1, observer);
    REQUIRE((approx_equal(angles[0], 45.0) && approx_equal(angles[1], 90.0)));
    angles = lookup_angles(source, velocity3, normal1, observer);
    REQUIRE((approx_equal(angles[0], -45.0) && approx_equal(angles[1], 90.0)));
    angles = lookup_angles(source, velocity4, normal1, observer);
    REQUIRE((approx_equal(angles[0], 0.0) && approx_equal(angles[1], 135.0)));

    // Check various TPP tilt / heading combinations
    angles = lookup_angles(source, velocity1, normal2, observer);
    REQUIRE(approx_equal(angles[1], 0.0));
    angles = lookup_angles(source, velocity2, normal2, observer);
    REQUIRE((approx_equal(angles[0], 90.0) && approx_equal(angles[1], 90.0)));
    angles = lookup_angles(source, velocity3, normal2, observer);
    REQUIRE((approx_equal(angles[0], -90.0) && approx_equal(angles[1], 90.0)));
    angles = lookup_angles(source, velocity1, normal3, observer);
    REQUIRE((approx_equal(angles[0], 0.0) && approx_equal(angles[1], 90.0)));
    angles = lookup_angles(source, velocity2, normal3, observer);
    REQUIRE((approx_equal(angles[0], 0.0) && approx_equal(angles[1], 90.0)));
    angles = lookup_angles(source, velocity3, normal3, observer);
    REQUIRE((approx_equal(angles[0], 0.0) && approx_equal(angles[1], 90.0)));
    angles = lookup_angles(source, velocity2, normal2, observer);
    REQUIRE((approx_equal(angles[0], 90.0) && approx_equal(angles[1], 90.0)));
    angles = lookup_angles(source, velocity1, normal4, observer);
    REQUIRE((approx_equal(angles[0], -45.0) && approx_equal(angles[1], 45.0)));
    angles = lookup_angles(source, velocity1, normal5, observer);
    REQUIRE((approx_equal(angles[0], 45.0) && approx_equal(angles[1], 45.0)));
}
