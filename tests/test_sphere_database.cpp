// Ported from the unittest blocks in sphere_database.d.  The D tests use
// std.math.approxEqual (relative epsilon 1e-6), so the port uses the C++
// approx_equal() from sphere_database.hpp to keep identical pass/fail
// semantics rather than Catch::Approx.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <limits>
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

TEST_CASE("GridAxis::bracket matches the binary search", "[sphere_database]")
{
    // GridAxis short-circuits the binary search when the axis is uniform, and
    // the whole optimisation rests on it giving an identical answer.  Sweep
    // both a uniform axis (the shortcut) and a deliberately irregular one (the
    // fall-back), including values on the knots, between them, and off both
    // ends.
    std::vector<double> uniform;
    for (int i = -180; i < 180; i += 10)
        uniform.push_back(static_cast<double>(i));

    std::vector<double> irregular = {-180.0, -179.0, -100.0, -12.5, 0.0, 3.0, 61.25, 170.0};

    for (const auto& values : {uniform, irregular})
    {
        GridAxis axis;
        axis.reset(values);
        for (int step = -2000; step <= 2000; step++)
        {
            const double val = 0.1 * static_cast<double>(step);
            const auto expected = bracket(values, val);
            const auto got = axis.bracket(val);
            REQUIRE(got[0] == expected[0]);
            REQUIRE(got[1] == expected[1]);
        }
        // Exactly on every knot, where an off-by-one would be easiest to miss.
        for (auto val : values)
        {
            const auto expected = bracket(values, val);
            const auto got = axis.bracket(val);
            REQUIRE(got[0] == expected[0]);
            REQUIRE(got[1] == expected[1]);
        }
    }
}

TEST_CASE("clamp", "[sphere_database]")
{
    REQUIRE(clamp01(0.5) == 0.5);
    REQUIRE(clamp01(-10.) == 0.);
    REQUIRE(clamp01(2.3) == 1.0);
}

TEST_CASE("acosc clamps its argument", "[sphere_database]")
{
    // dot/norm/normalize/cross are one-line vector helpers and are not tested
    // here; testing them would test vec3.hpp, not this file.  acosc is
    // real: it clamps before acos so round-off past the ends is not NaN.
    //
    // Do not assert bit-identity against std::acos on a literal.  acosc is a
    // cross-TU call that reaches libm at run time, while std::acos(0.5) here
    // is folded by the compiler at build time -- GCC folds via MPFR, correctly
    // rounded, and the two answers differ by an ULP.  The contract is the
    // clamping, not agreement with whichever libm the build happens to link,
    // so compare acosc against itself at the clamp points and allow a few ULP
    // in between.
    REQUIRE(std::fabs(acosc(0.5) - std::acos(0.5))
            <= 4.0 * std::numeric_limits<double>::epsilon());
    REQUIRE(acosc(10.0) == acosc(1.0));
    REQUIRE(acosc(-10.0) == acosc(-1.0));

    // The reason it exists: a dot product of unit vectors can land a hair
    // outside [-1, 1], where a bare acos returns NaN.
    REQUIRE(!std::isnan(acosc(1.0 + 1e-15)));
    REQUIRE(!std::isnan(acosc(-1.0 - 1e-15)));
}

TEST_CASE("lookup_angles", "[sphere_database]")
{
    const double r = 123.456;
    Vec3 source(0, 0, r);
    Vec3 observer(r, 0, 0);
    Vec3 velocity1(1, 0, 0);
    Vec3 velocity2(0, 1, 0);
    Vec3 velocity3(0, -1, 0);
    Vec3 velocity4(-1, 0, 0);
    const double s2 = std::sqrt(2.);
    Vec3 normal1(0, 0, 1);
    Vec3 normal2(s2, 0, s2);
    Vec3 normal3(-s2, 0, s2);
    Vec3 normal4(0, s2, s2);
    Vec3 normal5(0, -s2, s2);

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

TEST_CASE("Shaft power round-trips through a written sphere", "[sphere_database]")
{
    // No shipped database carries ShaftPower, so build one that does: take a
    // sphere out of the demo database, attach a power, write it, and read it
    // back.  This exercises the whole optional path -- write, detect, load,
    // blend -- which is otherwise unreachable with the data in this repository.
    Database source("databases/S-76D_M3.nod");
    REQUIRE(!source.has_shaft_power());

    auto sphere = source.interpolate_sphere(0.15, -4.0, 0.006);
    REQUIRE(std::isnan(sphere.shaft_power));
    const double power = 1234.5;
    sphere.shaft_power = power;

    const auto path = (std::filesystem::temp_directory_path() / "niceops_power_test.nod").string();
    write_sphere(sphere, path);

    Database reloaded(path);
    REQUIRE(reloaded.has_shaft_power());
    // One sphere, so every blend lands on it whatever the weights.
    const auto blend = reloaded.blend_at(0.15, -4.0, 0.006);
    REQUIRE(std::fabs(reloaded.shaft_power(blend) - power) < 1e-9);

    // And the levels still survive the round trip.
    REQUIRE(reloaded.sphere_radius() == sphere.radius);
    std::filesystem::remove(path);
}

TEST_CASE("Spectra survive interpolation and re-export", "[sphere_database]")
{
    // Spectra are optional, and none of the shipped databases has any.
    Database plain("databases/S-76D_M3.nod");
    REQUIRE(!plain.has_spectra());
    auto plain_sphere = plain.interpolate_sphere(0.15, -4.0, 0.006);
    REQUIRE(plain_sphere.frequencies.empty());
    REQUIRE(plain_sphere.amplitudes.empty());

    // Build one that does have them, by attaching a synthetic spectrum to an
    // interpolated sphere and writing it out.  A real spectral database comes
    // from panam; this keeps the read/blend/write path covered without one.
    const std::size_t bands = 4;
    plain_sphere.frequencies = {100.0, 200.0, 400.0, 800.0};
    plain_sphere.amplitudes.resize(plain_sphere.cells.size() * bands);
    for (std::size_t cell = 0; cell < plain_sphere.cells.size(); cell++)
    {
        for (std::size_t band = 0; band < bands; band++)
        {
            // Something cell- and band-dependent, so a transposed or shifted
            // read shows up rather than cancelling out.
            plain_sphere.amplitudes[cell * bands + band]
                    = 50.0 + 0.01 * static_cast<double>(cell) + static_cast<double>(band);
        }
    }

    const auto path = (std::filesystem::temp_directory_path() / "niceops_spectra_test.nod").string();
    write_sphere(plain_sphere, path);

    Database reloaded(path);
    REQUIRE(reloaded.has_spectra());
    auto round_tripped = reloaded.interpolate_sphere(0.15, -4.0, 0.006);
    REQUIRE(round_tripped.frequencies == plain_sphere.frequencies);
    REQUIRE(round_tripped.amplitudes.size() == plain_sphere.amplitudes.size());
    // One sphere, so the blend is a copy: every band at every cell must come
    // back exactly as written.
    double worst = 0.0;
    for (std::size_t i = 0; i < round_tripped.amplitudes.size(); i++)
        worst = std::max(worst, std::fabs(round_tripped.amplitudes[i] - plain_sphere.amplitudes[i]));
    REQUIRE(worst < 1e-9);

    std::filesystem::remove(path);
}
