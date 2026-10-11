// Ported from the unittest blocks in sphere_database.d.  The D tests use
// std.math.approxEqual (relative epsilon 1e-6), so the port uses the C++
// approx_equal() from sphere_database.hpp to keep identical pass/fail
// semantics rather than Catch::Approx.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <vector>

#include <netcdf.h>

#include "sphere_database.hpp"
#include "weighting.hpp"
#include "util/logger.hpp"

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
    // the whole optimization rests on it giving an identical answer.  Sweep
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

TEST_CASE("GridAxis::locate matches the bilinear cell it replaces", "[sphere_database]")
{
    // eval() takes its interpolation weights from locate() rather than from
    // periodic_bilinear_interp(), which is what makes the hot path arithmetic
    // instead of four knot lookups and a divide.  The two have to agree, so
    // interpolate a known ramp both ways and compare: locate() gives the weight
    // on its `hi` corner, and periodic_bilinear_interp() applied to the same
    // pair has to land on the same value.
    //
    // Sweeps a uniform axis (locate's arithmetic path), an irregular one (the
    // bracket fall-back) and the wrap-around cell at either end.
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
            const auto cell = axis.locate(val);
            const auto expected_pair = axis.bracket(val);
            REQUIRE(cell.lo == expected_pair[0]);
            REQUIRE(cell.hi == expected_pair[1]);

            // Corner values chosen so the interpolated answer is the weight
            // itself: 0 on `lo` and 1 on `hi`, held constant along the other
            // axis so only this one's weighting shows.
            const double reference = periodic_bilinear_interp(0.0, 0.0, 1.0, 1.0,
                    values[cell.lo], values[cell.hi], 0.0, 1.0, val, 0.5);
            REQUIRE(std::fabs(cell.fraction - reference) <= 1e-12);
        }
        for (auto val : values)
        {
            const auto cell = axis.locate(val);
            const double reference = periodic_bilinear_interp(0.0, 0.0, 1.0, 1.0,
                    values[cell.lo], values[cell.hi], 0.0, 1.0, val, 0.5);
            REQUIRE(std::fabs(cell.fraction - reference) <= 1e-12);
        }
    }
}

TEST_CASE("DirectionBounds really bounds the level", "[sphere_database]")
{
    // The observer loop skips a state when this bound puts it below the best
    // level found so far, so the whole peak search is only correct if the bound
    // is never exceeded.  Sweep a spread of source states and observers and
    // check it directly, in the energy form the hot path compares in.
    if (!std::filesystem::exists("databases/S-76D_M3.nod"))
        return;
    Database db("databases/S-76D_M3.nod");

    std::vector<Database::SourceState> states;
    for (double mu : {0.0, 0.08, 0.15, 0.22})
        for (double gamma : {-6.0, 0.0, 6.0})
            for (double ct : {0.005, 0.007})
            {
                // A few attitudes, so the sphere frame is not always level.
                states.push_back(db.source_state(Vec3(0.0, 0.0, 150.0), Vec3(60.0, 10.0, 2.0),
                        Vec3(0.1, -0.05, 1.0), mu, gamma, ct));
                states.push_back(db.source_state(Vec3(20.0, -40.0, 300.0), Vec3(-30.0, 55.0, -3.0),
                        Vec3(-0.08, 0.12, 1.0), mu, gamma, ct));
            }

    const auto bounds = db.direction_bounds(states);
    REQUIRE(bounds.usable());

    const double radius = db.sphere_radius();
    std::size_t checked = 0;
    for (std::size_t s = 0; s < states.size(); s++)
    {
        for (double x = -2000.0; x <= 2000.0; x += 137.0)
        {
            for (double y = -2000.0; y <= 2000.0; y += 149.0)
            {
                const Vec3 observer(x, y, 0.0);
                const Vec3 delta = observer - states[s].position;
                const double squared_range = delta.squaredNorm();
                if (squared_range <= radius * radius)
                    continue; // inside the sphere, where the bound does not apply
                const int bin = DirectionBounds::bin(states[s].frame.nose.dot(delta),
                        states[s].frame.normal.dot(delta),
                        states[s].frame.starboard.dot(delta), squared_range);
                const double level = db.level_at_observer(observer, states[s]);
                // bound >= level, in the energy form: numerator >= 10^(L/10) r^2.
                REQUIRE(bounds.numerator(s, bin)
                        >= std::pow(10.0, 0.1 * level) * squared_range * (1.0 - 1e-12));
                checked++;
            }
        }
    }
    REQUIRE(checked > 10000);
}

TEST_CASE("DirectionBounds still bounds the level at the nose and tail axes", "[sphere_database][bounds]")
{
    // A cell is filed under the bin of its center direction.  At theta = 0 or 180
    // that direction has no azimuth (up = across = 0), so every phi of the pole
    // row landed in one azimuth bin and the bins either side of it never saw the
    // pole's value -- while a query a degree off the axis, in one of those
    // bins, weights the pole cells at up to 0.9.  The bound was then exceeded by
    // up to 0.8 dB, and the windowed peak search's exactness rests on it never
    // being.  The ordinary sweep above puts only a handful of observers within 2
    // degrees of an axis; this one aims at them, over the whole operating
    // envelope (the violation needs a pole value that is the local maximum of its
    // bin, which not every blend has).
    if (!std::filesystem::exists("databases/S-76D_M3.nod"))
        return;
    Database db("databases/S-76D_M3.nod");

    std::vector<Database::SourceState> states;
    for (double mu = 0.0; mu <= 0.30001; mu += 0.02)
        for (double gamma = -12.0; gamma <= 12.0001; gamma += 3.0)
            for (double ct = 0.003; ct <= 0.0101; ct += 0.0014)
                states.push_back(db.source_state(Vec3(0.0, 0.0, 150.0), Vec3(60.0, 10.0, 2.0),
                        Vec3(0.1, -0.05, 1.0), mu, gamma, ct));
    const auto bounds = db.direction_bounds(states);
    REQUIRE(bounds.usable());

    const double radius = db.sphere_radius();
    std::size_t checked = 0, exceeded = 0;
    double worst_db = 0.0;
    for (std::size_t s = 0; s < states.size(); s++)
    {
        const auto& frame = states[s].frame;
        for (double theta : {0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0})
        {
            for (double sign : {1.0, -1.0}) // the nose axis, then the tail axis
            {
                for (double phi = -180.0; phi < 180.0; phi += 2.5)
                {
                    for (double range : {150.0, 400.0, 1200.0})
                    {
                    const double t = theta * pi / 180.0, p = phi * pi / 180.0;
                    Vec3 direction = sign * std::cos(t) * frame.nose
                            + std::sin(t) * std::cos(p) * frame.normal
                            + std::sin(t) * std::sin(p) * frame.starboard;
                    direction = direction.normalized();
                    const Vec3 observer = states[s].position + range * direction;
                    const Vec3 delta = observer - states[s].position;
                    const double squared_range = delta.squaredNorm();
                    if (squared_range <= radius * radius)
                        continue;
                    const int bin = DirectionBounds::bin(frame.nose.dot(delta),
                            frame.normal.dot(delta), frame.starboard.dot(delta), squared_range);
                    const double level = db.level_at_observer(observer, states[s]);
                    const double needed = std::pow(10.0, 0.1 * level) * squared_range * (1.0 - 1e-12);
                    checked++;
                    if (bounds.numerator(s, bin) < needed)
                    {
                        exceeded++;
                        worst_db = std::max(worst_db, 10.0 * std::log10(needed / bounds.numerator(s, bin)));
                    }
                    }
                }
            }
        }
    }
    INFO(exceeded << " of " << checked << " directions exceed the bound, worst " << worst_db << " dB");
    REQUIRE(checked > 5000000);
    REQUIRE(exceeded == 0);
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

TEST_CASE("speed_reference is read from the database root", "[sphere_database]")
{
    // Every database so far predates the attribute and was built against the
    // ground, so absent has to mean ground.  panam writes it with
    // netCDF4-python, as NC_CHAR text.
    Database source("databases/S-76D_M3.nod");
    REQUIRE(source.speed_reference() == SpeedReference::ground);
    const auto sphere = source.interpolate_sphere(0.15, -4.0, 0.006);
    const auto path = (std::filesystem::temp_directory_path() / "niceops_speed_reference.nod").string();

    const auto with_attribute = [&](const char* value) {
        write_sphere(sphere, path);
        int ncid = 0;
        REQUIRE(nc_open(path.c_str(), NC_WRITE, &ncid) == NC_NOERR);
        REQUIRE(nc_redef(ncid) == NC_NOERR);
        REQUIRE(nc_put_att_text(ncid, NC_GLOBAL, "speed_reference", std::strlen(value), value)
                == NC_NOERR);
        REQUIRE(nc_close(ncid) == NC_NOERR);
    };

    write_sphere(sphere, path);
    REQUIRE(Database(path).speed_reference() == SpeedReference::ground);
    with_attribute("air");
    REQUIRE(Database(path).speed_reference() == SpeedReference::air);
    with_attribute("ground");
    REQUIRE(Database(path).speed_reference() == SpeedReference::ground);
    // Anything else is refused rather than guessed at.
    with_attribute("inertial");
    REQUIRE_THROWS(Database(path));
    std::filesystem::remove(path);

    REQUIRE(parse_speed_reference("air") == SpeedReference::air);
    REQUIRE(!parse_speed_reference("Air").has_value());
    REQUIRE(std::string(to_string(SpeedReference::ground)) == "ground");
}

TEST_CASE("azimuth_reference is read from the database root and survives export", "[sphere_database]")
{
    // Absent is the track, as every database built before the attribute was;
    // the Database still says that it was absent, so a run can tell an
    // explicit --nose against the file from one against a default.
    Database source("databases/S-76D_M3.nod");
    REQUIRE(source.azimuth_reference() == AzimuthReference::track);
    REQUIRE_FALSE(source.recorded_azimuth_reference().has_value());
    auto sphere = source.interpolate_sphere(0.15, -4.0, 0.006);
    REQUIRE_FALSE(sphere.azimuth_reference.has_value());
    const auto path = (std::filesystem::temp_directory_path() / "niceops_azimuth_reference.nod").string();

    const auto with_attribute = [&](const char* value) {
        write_sphere(sphere, path);
        int ncid = 0;
        REQUIRE(nc_open(path.c_str(), NC_WRITE, &ncid) == NC_NOERR);
        REQUIRE(nc_redef(ncid) == NC_NOERR);
        REQUIRE(nc_put_att_text(ncid, NC_GLOBAL, "azimuth_reference", std::strlen(value), value)
                == NC_NOERR);
        REQUIRE(nc_close(ncid) == NC_NOERR);
    };

    // An export of a database that records nothing records nothing.
    write_sphere(sphere, path);
    REQUIRE_FALSE(Database(path).recorded_azimuth_reference().has_value());
    with_attribute("heading");
    {
        const Database reloaded(path);
        REQUIRE(reloaded.azimuth_reference() == AzimuthReference::heading);
        REQUIRE(reloaded.recorded_azimuth_reference() == AzimuthReference::heading);
        // And re-exported, it carries the attribute on.
        const auto again = reloaded.interpolate_sphere(0.15, -4.0, 0.006);
        REQUIRE(again.azimuth_reference == AzimuthReference::heading);
        const auto copy = (std::filesystem::temp_directory_path() / "niceops_azimuth_copy.nod").string();
        write_sphere(again, copy);
        REQUIRE(Database(copy).recorded_azimuth_reference() == AzimuthReference::heading);
        std::filesystem::remove(copy);
    }
    with_attribute("track");
    REQUIRE(Database(path).recorded_azimuth_reference() == AzimuthReference::track);
    // Anything else is refused rather than guessed at, padding aside.
    with_attribute("heading  ");
    REQUIRE(Database(path).azimuth_reference() == AzimuthReference::heading);
    with_attribute("nose");
    REQUIRE_THROWS_WITH(Database(path), Catch::Matchers::ContainsSubstring("azimuth_reference is 'nose'"));
    std::filesystem::remove(path);

    REQUIRE(parse_azimuth_reference("heading") == AzimuthReference::heading);
    REQUIRE(parse_azimuth_reference("track") == AzimuthReference::track);
    REQUIRE(!parse_azimuth_reference("Heading").has_value());
    REQUIRE(std::string(to_string(AzimuthReference::heading)) == "heading");
    REQUIRE(std::string(to_string(AzimuthReference::track)) == "track");
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
            // read shows up rather than canceling out.
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

TEST_CASE("a database with no condition groups is refused, not crashed", "[sphere_database]")
{
    // Every field the constructor fills comes out of the loop over the groups.
    // With none, the loop never ran and the tail of the constructor read an
    // uninitialized sphere radius and dereferenced max_element's end() on the
    // empty condition lists.  That segfaulted, with nothing on stderr, after a
    // log that looked like a successful load.
    const auto path = (std::filesystem::temp_directory_path() / "niceops_no_groups.nod").string();

    int file = 0;
    REQUIRE(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file) == NC_NOERR);
    int scalar_dim = 0;
    REQUIRE(nc_def_dim(file, "scalar", 1, &scalar_dim) == NC_NOERR);
    int same_grid = 0;
    REQUIRE(nc_def_var(file, "same_grid", NC_SHORT, 1, &scalar_dim, &same_grid) == NC_NOERR);
    int version = 0;
    REQUIRE(nc_def_var(file, "database_version", NC_INT, 1, &scalar_dim, &version) == NC_NOERR);
    REQUIRE(nc_enddef(file) == NC_NOERR);
    const short one = 1;
    const int one_int = 1;
    const std::size_t start = 0, count = 1;
    REQUIRE(nc_put_vara_short(file, same_grid, &start, &count, &one) == NC_NOERR);
    REQUIRE(nc_put_vara_int(file, version, &start, &count, &one_int) == NC_NOERR);
    REQUIRE(nc_close(file) == NC_NOERR);

    REQUIRE_THROWS_AS(Database(path), std::runtime_error);
    std::filesystem::remove(path);
}

namespace {

// Writes a minimal .nod with one group per (advance_ratio, flight_path_angle,
// thrust_coefficient, dBA) condition, each a uniform 2x2 (phi, theta) sphere
// (dBA and EAA=0 the same at every channel) -- enough for a real GridAxis on
// both axes without needing to hand-verify any phi/theta interpolation
// arithmetic, which is not what these tests are about. radius is always
// sphere_radius_m, so eval() at range == radius has spreading_loss == 0 and
// (EAA == 0) excess_attenuation == 0, leaving ground_level == sphere_level +
// level_correction_db exactly -- the cleanest point to query panam's
// load-factor scaling law against.
std::string write_fixed_load_factor_test_database(const std::string& name,
        const std::vector<std::array<double, 4>>& conditions, bool fixed_load_factor,
        double sphere_radius_m = 30.48, nc_type storage = NC_DOUBLE, double sphere_1_phi_shift = 0.0,
        double first_channel_eaa = 0.0)
{
    const auto path = (std::filesystem::temp_directory_path() / name).string();
    int file = 0;
    REQUIRE(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file) == NC_NOERR);
    int scalar_dim = 0;
    REQUIRE(nc_def_dim(file, "scalar", 1, &scalar_dim) == NC_NOERR);
    int same_grid_var = 0, version_var = 0, flag_var = 0;
    REQUIRE(nc_def_var(file, "same_grid", NC_SHORT, 1, &scalar_dim, &same_grid_var) == NC_NOERR);
    REQUIRE(nc_def_var(file, "database_version", NC_INT, 1, &scalar_dim, &version_var) == NC_NOERR);
    REQUIRE(nc_def_var(file, "fixed_load_factor", NC_SHORT, 1, &scalar_dim, &flag_var) == NC_NOERR);

    struct GroupVars
    {
        int group, radius, scale, mu, gamma, ct, phi, theta, dba, eaa;
    };
    // Full-sphere coverage (every phi x every theta), or direction_bounds'
    // padding logic legitimately leaves bins this fixture never touches at
    // -infinity -- correct (a direction with no nearby grid data cannot be
    // bounded), but wrong for a test that wants to sweep observers broadly.
    // A real AAM sphere always covers the whole range for the same reason.
    const std::vector<double> phis = {-90.0, -30.0, 30.0, 90.0};
    const std::vector<double> thetas = {0.0, 60.0, 120.0, 180.0};
    const std::size_t n_channels = phis.size() * thetas.size();

    std::vector<GroupVars> groups;
    for (std::size_t i = 0; i < conditions.size(); i++)
    {
        int group = 0;
        REQUIRE(nc_def_grp(file, ("sphere" + std::to_string(i)).c_str(), &group) == NC_NOERR);
        int radii_dim = 0, condition_dim = 0, channel_dim = 0;
        REQUIRE(nc_def_dim(group, "radii", 1, &radii_dim) == NC_NOERR);
        REQUIRE(nc_def_dim(group, "condition", 1, &condition_dim) == NC_NOERR);
        REQUIRE(nc_def_dim(group, "channels", n_channels, &channel_dim) == NC_NOERR);
        GroupVars v {group, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        REQUIRE(nc_def_var(group, "radius", storage, 1, &radii_dim, &v.radius) == NC_NOERR);
        REQUIRE(nc_def_var(group, "rotor_scale", storage, 1, &condition_dim, &v.scale) == NC_NOERR);
        REQUIRE(nc_def_var(group, "advance_ratio", storage, 1, &condition_dim, &v.mu) == NC_NOERR);
        REQUIRE(nc_def_var(group, "flight_path_angle", storage, 1, &condition_dim, &v.gamma)
                == NC_NOERR);
        REQUIRE(nc_def_var(group, "thrust_coefficient", storage, 1, &condition_dim, &v.ct)
                == NC_NOERR);
        REQUIRE(nc_def_var(group, "phi", storage, 1, &channel_dim, &v.phi) == NC_NOERR);
        REQUIRE(nc_def_var(group, "theta", storage, 1, &channel_dim, &v.theta) == NC_NOERR);
        REQUIRE(nc_def_var(group, "dBA", storage, 1, &channel_dim, &v.dba) == NC_NOERR);
        REQUIRE(nc_def_var(group, "EAA", storage, 1, &channel_dim, &v.eaa) == NC_NOERR);
        groups.push_back(v);
    }
    REQUIRE(nc_enddef(file) == NC_NOERR);

    const short one = 1;
    const short flag = fixed_load_factor ? 1 : 0;
    const int version = 1;
    REQUIRE(nc_put_var_short(file, same_grid_var, &one) == NC_NOERR);
    REQUIRE(nc_put_var_int(file, version_var, &version) == NC_NOERR);
    REQUIRE(nc_put_var_short(file, flag_var, &flag) == NC_NOERR);

    const double scale = 1.0;
    std::vector<double> phi_channel, theta_channel;
    for (double p : phis)
        for (double t : thetas)
        {
            phi_channel.push_back(p);
            theta_channel.push_back(t);
        }
    for (std::size_t i = 0; i < conditions.size(); i++)
    {
        const auto& v = groups[i];
        // Uniform across every channel: these tests are about the
        // load-factor correction, not phi/theta interpolation, so the base
        // level is the same everywhere on a given condition's sphere.
        std::vector<double> eaa(n_channels, 0.0);
        if (i == 0)
            eaa[0] = first_channel_eaa;
        const std::vector<double> dba(n_channels, conditions[i][3]);
        REQUIRE(nc_put_var_double(v.group, v.radius, &sphere_radius_m) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.scale, &scale) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.mu, &conditions[i][0]) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.gamma, &conditions[i][1]) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.ct, &conditions[i][2]) == NC_NOERR);
        std::vector<double> labels = phi_channel;
        if (i == 1)
            for (double& label : labels)
                label += sphere_1_phi_shift;
        REQUIRE(nc_put_var_double(v.group, v.phi, labels.data()) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.theta, theta_channel.data()) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.dba, dba.data()) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.eaa, eaa.data()) == NC_NOERR);
    }
    REQUIRE(nc_close(file) == NC_NOERR);
    return path;
}

// Three conditions at a shared reference CT, forming a non-degenerate
// triangle in (mu, gamma) so blend_at() has a real 2D locate to walk rather
// than the single-point nearest-condition fallback.
constexpr double test_ct_reference = 0.005;
const std::vector<std::array<double, 4>> fixed_load_factor_conditions = {
    {0.0, 0.0, test_ct_reference, 100.0},
    {1.0, 0.0, test_ct_reference, 110.0},
    {0.0, 1.0, test_ct_reference, 120.0},
};

// Writes a minimal .nod with one group per (advance_ratio, flight_path_angle,
// thrust_coefficient, dBA) condition, every condition sharing the same 4x4
// (phi, theta) sphere and, when with_spectrum, the same 2 frequency bands --
// the case shared_grid_and_frequency is meant to deduplicate.  When
// shared_grid_and_frequency is true, phi/theta/frequency are written once at
// the root and every group's dBA/EAA/amplitude resolve them by inheriting the
// root's dimensions (a plain netCDF4 group scoping rule), exactly as panam's
// writer does; when false, every group carries its own copies, so a test can
// build both forms of the same conditions and require identical eval()
// results between them.
//
// Two options for the blend and layout tests: channel_dba(condition, channel),
// when given, sets each channel's dBA (use -inf for a direction ambient gating
// emptied); reverse_phi_rows lays the channels -- and so the amplitude rows --
// out with phi descending, the way panam's unsorted completion order arrives,
// and writes each amplitude as 200 + that row's phi + 10 * band so a test can
// see which direction every spectrum landed in.  bands and cell_amplitude(condition,
// cell, band) replace the two bands and their levels.
std::string write_shared_grid_test_database(const std::string& name,
        const std::vector<std::array<double, 4>>& conditions, bool shared_grid_and_frequency,
        bool with_spectrum, double sphere_radius_m = 30.48,
        std::function<double(std::size_t, std::size_t)> channel_dba = {},
        bool reverse_phi_rows = false, std::vector<double> bands = {500.0, 2000.0},
        std::function<double(std::size_t, std::size_t, std::size_t)> cell_amplitude = {})
{
    const auto path = (std::filesystem::temp_directory_path() / name).string();
    int file = 0;
    REQUIRE(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file) == NC_NOERR);
    int scalar_dim = 0;
    REQUIRE(nc_def_dim(file, "scalar", 1, &scalar_dim) == NC_NOERR);
    int same_grid_var = 0, version_var = 0, shared_var = 0, fixed_lf_var = 0;
    REQUIRE(nc_def_var(file, "same_grid", NC_SHORT, 1, &scalar_dim, &same_grid_var) == NC_NOERR);
    // These fixtures give every sphere one thrust coefficient, which is what a
    // fixed_load_factor database is.  Written without the flag, the thrust axis
    // normalized 0/0 and the triangulation was built on NaN -- the loader now
    // refuses that, so the flag says what the fixture is.
    REQUIRE(nc_def_var(file, "fixed_load_factor", NC_SHORT, 1, &scalar_dim, &fixed_lf_var) == NC_NOERR);
    REQUIRE(nc_def_var(file, "database_version", NC_INT, 1, &scalar_dim, &version_var) == NC_NOERR);
    REQUIRE(nc_def_var(file, "shared_grid_and_frequency", NC_SHORT, 1, &scalar_dim, &shared_var)
            == NC_NOERR);

    const std::vector<double> phis = {-90.0, -30.0, 30.0, 90.0};
    const std::vector<double> thetas = {0.0, 60.0, 120.0, 180.0};
    const std::size_t n_channels = phis.size() * thetas.size();
    std::vector<double> phi_channel, theta_channel;
    std::vector<double> row_phis = phis;
    if (reverse_phi_rows)
        std::reverse(row_phis.begin(), row_phis.end());
    for (double p : row_phis)
        for (double t : thetas)
        {
            phi_channel.push_back(p);
            theta_channel.push_back(t);
        }

    int root_channels_dim = 0, root_phi_var = 0, root_theta_var = 0;
    int root_frequency_dim = 0, root_frequency_var = 0, root_phi_dim = 0, root_theta_dim = 0;
    if (shared_grid_and_frequency)
    {
        REQUIRE(nc_def_dim(file, "channels", n_channels, &root_channels_dim) == NC_NOERR);
        REQUIRE(nc_def_var(file, "phi", NC_DOUBLE, 1, &root_channels_dim, &root_phi_var)
                == NC_NOERR);
        REQUIRE(nc_def_var(file, "theta", NC_DOUBLE, 1, &root_channels_dim, &root_theta_var)
                == NC_NOERR);
        if (with_spectrum)
        {
            REQUIRE(nc_def_dim(file, "frequency", bands.size(), &root_frequency_dim)
                    == NC_NOERR);
            REQUIRE(nc_def_var(file, "frequency", NC_DOUBLE, 1, &root_frequency_dim,
                    &root_frequency_var) == NC_NOERR);
            REQUIRE(nc_def_dim(file, "PHI", phis.size(), &root_phi_dim) == NC_NOERR);
            REQUIRE(nc_def_dim(file, "THETA", thetas.size(), &root_theta_dim) == NC_NOERR);
        }
    }

    struct GroupVars
    {
        int group, radius, scale, mu, gamma, ct, phi, theta, dba, eaa, frequency, amplitude;
    };
    std::vector<GroupVars> groups;
    for (std::size_t i = 0; i < conditions.size(); i++)
    {
        int group = 0;
        REQUIRE(nc_def_grp(file, ("sphere" + std::to_string(i)).c_str(), &group) == NC_NOERR);
        int radii_dim = 0, condition_dim = 0;
        REQUIRE(nc_def_dim(group, "radii", 1, &radii_dim) == NC_NOERR);
        REQUIRE(nc_def_dim(group, "condition", 1, &condition_dim) == NC_NOERR);
        GroupVars v {group, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        REQUIRE(nc_def_var(group, "radius", NC_DOUBLE, 1, &radii_dim, &v.radius) == NC_NOERR);
        REQUIRE(nc_def_var(group, "rotor_scale", NC_DOUBLE, 1, &condition_dim, &v.scale) == NC_NOERR);
        REQUIRE(nc_def_var(group, "advance_ratio", NC_DOUBLE, 1, &condition_dim, &v.mu) == NC_NOERR);
        REQUIRE(nc_def_var(group, "flight_path_angle", NC_DOUBLE, 1, &condition_dim, &v.gamma)
                == NC_NOERR);
        REQUIRE(nc_def_var(group, "thrust_coefficient", NC_DOUBLE, 1, &condition_dim, &v.ct)
                == NC_NOERR);
        int channel_dim = root_channels_dim;
        if (!shared_grid_and_frequency)
        {
            REQUIRE(nc_def_dim(group, "channels", n_channels, &channel_dim) == NC_NOERR);
            REQUIRE(nc_def_var(group, "phi", NC_DOUBLE, 1, &channel_dim, &v.phi) == NC_NOERR);
            REQUIRE(nc_def_var(group, "theta", NC_DOUBLE, 1, &channel_dim, &v.theta) == NC_NOERR);
        }
        REQUIRE(nc_def_var(group, "dBA", NC_DOUBLE, 1, &channel_dim, &v.dba) == NC_NOERR);
        REQUIRE(nc_def_var(group, "EAA", NC_DOUBLE, 1, &channel_dim, &v.eaa) == NC_NOERR);
        if (with_spectrum)
        {
            int frequency_dim = root_frequency_dim, phi_dim = root_phi_dim,
                theta_dim = root_theta_dim;
            if (!shared_grid_and_frequency)
            {
                REQUIRE(nc_def_dim(group, "frequency", bands.size(), &frequency_dim) == NC_NOERR);
                REQUIRE(nc_def_dim(group, "PHI", phis.size(), &phi_dim) == NC_NOERR);
                REQUIRE(nc_def_dim(group, "THETA", thetas.size(), &theta_dim) == NC_NOERR);
                REQUIRE(nc_def_var(group, "frequency", NC_DOUBLE, 1, &frequency_dim,
                        &v.frequency) == NC_NOERR);
            }
            const int amplitude_dims[3] = {phi_dim, theta_dim, frequency_dim};
            REQUIRE(nc_def_var(group, "amplitude", NC_DOUBLE, 3, amplitude_dims, &v.amplitude)
                    == NC_NOERR);
        }
        groups.push_back(v);
    }
    REQUIRE(nc_enddef(file) == NC_NOERR);

    const short one = 1;
    const short shared_flag = shared_grid_and_frequency ? 1 : 0;
    const int version = 1;
    REQUIRE(nc_put_var_short(file, same_grid_var, &one) == NC_NOERR);
    REQUIRE(nc_put_var_int(file, version_var, &version) == NC_NOERR);
    REQUIRE(nc_put_var_short(file, shared_var, &shared_flag) == NC_NOERR);
    {
        bool constant_ct = conditions.size() > 1;
        for (const auto& condition : conditions)
            constant_ct = constant_ct && condition[2] == conditions.front()[2];
        const short fixed_flag = constant_ct ? 1 : 0;
        REQUIRE(nc_put_var_short(file, fixed_lf_var, &fixed_flag) == NC_NOERR);
    }
    if (shared_grid_and_frequency)
    {
        REQUIRE(nc_put_var_double(file, root_phi_var, phi_channel.data()) == NC_NOERR);
        REQUIRE(nc_put_var_double(file, root_theta_var, theta_channel.data()) == NC_NOERR);
        if (with_spectrum)
            REQUIRE(nc_put_var_double(file, root_frequency_var, bands.data()) == NC_NOERR);
    }

    const double scale = 1.0;
    for (std::size_t i = 0; i < conditions.size(); i++)
    {
        const auto& v = groups[i];
        const std::vector<double> eaa(n_channels, 0.0);
        std::vector<double> dba(n_channels, conditions[i][3]);
        if (channel_dba)
            for (std::size_t c = 0; c < n_channels; c++)
                dba[c] = channel_dba(i, c);
        REQUIRE(nc_put_var_double(v.group, v.radius, &sphere_radius_m) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.scale, &scale) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.mu, &conditions[i][0]) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.gamma, &conditions[i][1]) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.ct, &conditions[i][2]) == NC_NOERR);
        if (!shared_grid_and_frequency)
        {
            REQUIRE(nc_put_var_double(v.group, v.phi, phi_channel.data()) == NC_NOERR);
            REQUIRE(nc_put_var_double(v.group, v.theta, theta_channel.data()) == NC_NOERR);
        }
        REQUIRE(nc_put_var_double(v.group, v.dba, dba.data()) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.eaa, eaa.data()) == NC_NOERR);
        if (with_spectrum)
        {
            if (!shared_grid_and_frequency)
                REQUIRE(nc_put_var_double(v.group, v.frequency, bands.data()) == NC_NOERR);
            // Distinguishable per condition and per band, uniform over phi/theta
            // (this test is about grid/frequency dedup, not spatial variation).
            std::vector<double> amplitude(phis.size() * thetas.size() * bands.size());
            for (std::size_t cell = 0; cell < phis.size() * thetas.size(); cell++)
                for (std::size_t band = 0; band < bands.size(); band++)
                    amplitude[cell * bands.size() + band] = cell_amplitude
                            ? cell_amplitude(i, cell, band)
                            : reverse_phi_rows
                            ? 200.0 + phi_channel[cell] + 10.0 * static_cast<double>(band)
                            : conditions[i][3] + 10.0 * static_cast<double>(band);
            REQUIRE(nc_put_var_double(v.group, v.amplitude, amplitude.data()) == NC_NOERR);
        }
    }
    REQUIRE(nc_close(file) == NC_NOERR);
    return path;
}

// A single condition on a 4x4 sphere with one channel's dBA at -inf, the way
// ambient gating legitimately empties a direction.  self_test() re-reads
// this raw value from the file to compare against eval(), so it has to know
// eval() can only ever answer no_energy_floor_db for it, never -inf itself.
std::string write_test_database_with_silent_direction(const std::string& name)
{
    const auto path = (std::filesystem::temp_directory_path() / name).string();
    int file = 0;
    REQUIRE(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file) == NC_NOERR);
    int scalar_dim = 0;
    REQUIRE(nc_def_dim(file, "scalar", 1, &scalar_dim) == NC_NOERR);
    int same_grid_var = 0;
    REQUIRE(nc_def_var(file, "same_grid", NC_SHORT, 1, &scalar_dim, &same_grid_var) == NC_NOERR);

    int group = 0;
    REQUIRE(nc_def_grp(file, "sphere0", &group) == NC_NOERR);
    int radii_dim = 0, condition_dim = 0, channel_dim = 0;
    REQUIRE(nc_def_dim(group, "radii", 1, &radii_dim) == NC_NOERR);
    REQUIRE(nc_def_dim(group, "condition", 1, &condition_dim) == NC_NOERR);
    const std::vector<double> phis = {-90.0, -30.0, 30.0, 90.0};
    const std::vector<double> thetas = {0.0, 60.0, 120.0, 180.0};
    const std::size_t n_channels = phis.size() * thetas.size();
    REQUIRE(nc_def_dim(group, "channels", n_channels, &channel_dim) == NC_NOERR);
    int radius_var = 0, scale_var = 0, mu_var = 0, gamma_var = 0, ct_var = 0, phi_var = 0,
        theta_var = 0, dba_var = 0, eaa_var = 0;
    REQUIRE(nc_def_var(group, "radius", NC_DOUBLE, 1, &radii_dim, &radius_var) == NC_NOERR);
    REQUIRE(nc_def_var(group, "rotor_scale", NC_DOUBLE, 1, &condition_dim, &scale_var) == NC_NOERR);
    REQUIRE(nc_def_var(group, "advance_ratio", NC_DOUBLE, 1, &condition_dim, &mu_var) == NC_NOERR);
    REQUIRE(nc_def_var(group, "flight_path_angle", NC_DOUBLE, 1, &condition_dim, &gamma_var)
            == NC_NOERR);
    REQUIRE(nc_def_var(group, "thrust_coefficient", NC_DOUBLE, 1, &condition_dim, &ct_var)
            == NC_NOERR);
    REQUIRE(nc_def_var(group, "phi", NC_DOUBLE, 1, &channel_dim, &phi_var) == NC_NOERR);
    REQUIRE(nc_def_var(group, "theta", NC_DOUBLE, 1, &channel_dim, &theta_var) == NC_NOERR);
    REQUIRE(nc_def_var(group, "dBA", NC_DOUBLE, 1, &channel_dim, &dba_var) == NC_NOERR);
    REQUIRE(nc_def_var(group, "EAA", NC_DOUBLE, 1, &channel_dim, &eaa_var) == NC_NOERR);
    REQUIRE(nc_enddef(file) == NC_NOERR);

    const short one = 1;
    REQUIRE(nc_put_var_short(file, same_grid_var, &one) == NC_NOERR);
    const double radius = 30.48, scale = 1.0, mu = 0.0, gamma = 0.0, ct = 0.005;
    REQUIRE(nc_put_var_double(group, radius_var, &radius) == NC_NOERR);
    REQUIRE(nc_put_var_double(group, scale_var, &scale) == NC_NOERR);
    REQUIRE(nc_put_var_double(group, mu_var, &mu) == NC_NOERR);
    REQUIRE(nc_put_var_double(group, gamma_var, &gamma) == NC_NOERR);
    REQUIRE(nc_put_var_double(group, ct_var, &ct) == NC_NOERR);
    std::vector<double> phi_channel, theta_channel;
    for (double p : phis)
        for (double t : thetas)
        {
            phi_channel.push_back(p);
            theta_channel.push_back(t);
        }
    REQUIRE(nc_put_var_double(group, phi_var, phi_channel.data()) == NC_NOERR);
    REQUIRE(nc_put_var_double(group, theta_var, theta_channel.data()) == NC_NOERR);
    std::vector<double> dba(n_channels, 100.0);
    dba[0] = -std::numeric_limits<double>::infinity();
    const std::vector<double> eaa(n_channels, 0.0);
    REQUIRE(nc_put_var_double(group, dba_var, dba.data()) == NC_NOERR);
    REQUIRE(nc_put_var_double(group, eaa_var, eaa.data()) == NC_NOERR);
    REQUIRE(nc_close(file) == NC_NOERR);
    return path;
}

} // namespace

TEST_CASE("fixed_load_factor: query at or below the reference applies no correction",
        "[sphere_database]")
{
    const auto path = write_fixed_load_factor_test_database(
            "niceops_fixed_lf_no_correction.nod", fixed_load_factor_conditions, true);
    Database db(path);

    // At (mu=0, gamma=0) the blend lands on condition 0 alone -- "a sphere
    // exported at an operating condition the database already holds comes
    // back bit-for-bit as that sphere's own data" (interpolate_sphere's own
    // comment) -- so this is dBA=100.0 exactly, with no interpolation
    // arithmetic to account for.
    for (double ct : {test_ct_reference * 0.2, test_ct_reference * 0.5, test_ct_reference})
    {
        const auto blend = db.blend_at(0.0, 0.0, ct);
        REQUIRE(blend.level_correction_db == 0.0);
        // range == radius: spreading_loss and (EAA == 0) excess_attenuation
        // both vanish, so eval() is exactly the stored level.
        REQUIRE(approx_equal(db.eval(db.sphere_radius(), 0.0, 90.0, blend), 100.0));
    }
    std::filesystem::remove(path);
}

TEST_CASE("fixed_load_factor: query above the reference scales by 20*log10(load factor)",
        "[sphere_database]")
{
    const auto path = write_fixed_load_factor_test_database(
            "niceops_fixed_lf_correction.nod", fixed_load_factor_conditions, true);
    Database db(path);

    for (double load_factor : {1.1, 2.0, 5.0, 50.0}) // no upper bound, by design
    {
        const double ct = test_ct_reference * load_factor;
        const auto blend = db.blend_at(0.0, 0.0, ct);
        const double expected = 20.0 * std::log10(load_factor);
        REQUIRE(approx_equal(blend.level_correction_db, expected));
        REQUIRE(approx_equal(db.eval(db.sphere_radius(), 0.0, 90.0, blend), 100.0 + expected));
    }
    std::filesystem::remove(path);
}

TEST_CASE("source coverage distinguishes clamping from analytic load scaling", "[sphere_database]")
{
    const auto path = write_fixed_load_factor_test_database(
            "niceops_source_coverage.nod", fixed_load_factor_conditions, true);
    Database db(path);
    const auto inside = db.blend_at(0.0, 0.0, test_ct_reference);
    CHECK_FALSE(inside.advance_ratio_clamped);
    CHECK_FALSE(inside.flight_path_angle_clamped);
    CHECK_FALSE(inside.thrust_coefficient_clamped);
    const auto outside = db.blend_at(-0.1, 1.1, 2.0 * test_ct_reference);
    CHECK(outside.advance_ratio_clamped);
    CHECK(outside.flight_path_angle_clamped);
    CHECK_FALSE(outside.thrust_coefficient_clamped);
    CHECK(outside.level_correction_db > 0.0);
    const auto boundary = db.blend_at(0.0, 1.0, 2.0 * test_ct_reference);
    CHECK_FALSE(boundary.advance_ratio_clamped);
    CHECK_FALSE(boundary.flight_path_angle_clamped);
    CHECK(db.eval(db.sphere_radius(), 0.0, 90.0, outside)
            == db.eval(db.sphere_radius(), 0.0, 90.0, boundary));
    std::filesystem::remove(path);
}

TEST_CASE("fixed_load_factor: cull_range accounts for the correction", "[sphere_database]")
{
    const auto path = write_fixed_load_factor_test_database(
            "niceops_fixed_lf_cull.nod", fixed_load_factor_conditions, true);
    Database db(path);

    // A cutoff comfortably above the stored 100 dBA but below the corrected
    // level at 3x load factor (100 + 20*log10(3) ~= 109.5): cull_range must
    // not treat this state as culled at any finite range, or eval() at a
    // range where it is still genuinely audible would be skipped.
    const auto blend = db.blend_at(0.0, 0.0, test_ct_reference * 3.0);
    const double cutoff = 105.0;
    const double range = db.cull_range(blend, cutoff);
    REQUIRE(std::isfinite(range));

    // Sweep ranges out to (and past) the bound and check the invariant
    // cull_range exists to guarantee: eval() must be below cutoff at or
    // beyond it.  (EAA == 0 in this fixture, so level is non-increasing in
    // range and one boundary check would do, but sweeping is what the real
    // DirectionBounds test above does and catches an off-by-factor error
    // the same way.)
    for (double r = range; r <= range * 1.5; r += db.sphere_radius())
        REQUIRE(db.eval(r, 0.0, 90.0, blend) <= cutoff + 1e-9);
    // And it must not be trivially large: just below it, the level this
    // fixture actually has is still above cutoff, so a bound of 0 (or
    // unbounded) would not have been wrong, but would have been useless.
    REQUIRE(db.eval(range * 0.9, 0.0, 90.0, blend) > cutoff);
    std::filesystem::remove(path);
}

TEST_CASE("fixed_load_factor: direction_bounds accounts for varying load factor "
        "at the same (mu, gamma)",
        "[sphere_database]")
{
    // Two states share (mu, gamma) = (0, 0) -- the same 2D blend indices --
    // but query different load factors, which the 2D locate cannot see: this
    // is exactly the case the indices-alone grouping in direction_bounds
    // would have merged, applying one state's correction to both.
    const auto path = write_fixed_load_factor_test_database(
            "niceops_fixed_lf_bounds.nod", fixed_load_factor_conditions, true);
    Database db(path);

    std::vector<Database::SourceState> states;
    for (double load_factor : {1.0, 4.0})
    {
        states.push_back(db.source_state(Vec3(0.0, 0.0, 100.0), Vec3(50.0, 0.0, 0.0),
                Vec3(0.0, 0.0, 1.0), 0.0, 0.0, test_ct_reference * load_factor));
    }
    const auto bounds = db.direction_bounds(states);
    REQUIRE(bounds.usable());

    const double radius = db.sphere_radius();
    for (std::size_t s = 0; s < states.size(); s++)
    {
        for (double x = -300.0; x <= 300.0; x += 37.0)
        {
            for (double y = -300.0; y <= 300.0; y += 41.0)
            {
                const Vec3 observer(x, y, 0.0);
                const Vec3 delta = observer - states[s].position;
                const double squared_range = delta.squaredNorm();
                if (squared_range <= radius * radius)
                    continue;
                const int bin = DirectionBounds::bin(states[s].frame.nose.dot(delta),
                        states[s].frame.normal.dot(delta),
                        states[s].frame.starboard.dot(delta), squared_range);
                const double level = db.level_at_observer(observer, states[s]);
                REQUIRE(bounds.numerator(s, bin)
                        >= std::pow(10.0, 0.1 * level) * squared_range * (1.0 - 1e-9));
            }
        }
    }
    std::filesystem::remove(path);
}

TEST_CASE("fixed_load_factor requires thrust_coefficient to actually be constant",
        "[sphere_database]")
{
    // The flag is trusted for the triangulation (dropping CT from it entirely)
    // but not blindly: a file claiming fixed_load_factor with varying CT would
    // otherwise silently stop scaling correctly rather than fail loudly.
    auto conditions = fixed_load_factor_conditions;
    conditions[1][2] = test_ct_reference * 2.0; // now inconsistent with the others
    const auto path = write_fixed_load_factor_test_database(
            "niceops_fixed_lf_inconsistent.nod", conditions, true);

    REQUIRE_THROWS_AS(Database(path), std::runtime_error);
    std::filesystem::remove(path);
}

TEST_CASE("shared_grid_and_frequency reproduces the per-group form exactly, "
        "levels and spectra alike",
        "[sphere_database]")
{
    const std::vector<std::array<double, 4>> conditions = {
        {0.0, 0.0, 0.005, 100.0},
        {1.0, 0.0, 0.005, 110.0},
        {0.0, 1.0, 0.005, 120.0},
    };
    const auto shared_path = write_shared_grid_test_database(
            "niceops_shared_grid.nod", conditions, true, true);
    const auto per_group_path = write_shared_grid_test_database(
            "niceops_per_group_grid.nod", conditions, false, true);
    Database shared_db(shared_path);
    Database per_group_db(per_group_path);

    for (const auto& c : conditions)
    {
        const auto shared_sphere = shared_db.interpolate_sphere(c[0], c[1], c[2]);
        const auto per_group_sphere = per_group_db.interpolate_sphere(c[0], c[1], c[2]);

        REQUIRE(shared_sphere.phis == per_group_sphere.phis);
        REQUIRE(shared_sphere.thetas == per_group_sphere.thetas);
        REQUIRE(shared_sphere.frequencies == per_group_sphere.frequencies);
        REQUIRE(shared_sphere.cells.size() == per_group_sphere.cells.size());
        for (std::size_t cell = 0; cell < shared_sphere.cells.size(); cell++)
        {
            CHECK(approx_equal(shared_sphere.cells[cell].dBA, per_group_sphere.cells[cell].dBA));
            CHECK(approx_equal(shared_sphere.cells[cell].EAA, per_group_sphere.cells[cell].EAA));
        }
        REQUIRE(shared_sphere.amplitudes.size() == per_group_sphere.amplitudes.size());
        for (std::size_t k = 0; k < shared_sphere.amplitudes.size(); k++)
        {
            CHECK(approx_equal(shared_sphere.amplitudes[k], per_group_sphere.amplitudes[k]));
        }
    }
    std::filesystem::remove(shared_path);
    std::filesystem::remove(per_group_path);
}

TEST_CASE("shared_grid_and_frequency: absent is the same as false (legacy files)",
        "[sphere_database]")
{
    // Every file every other test in this suite loads predates this flag, so
    // this is really just documenting the default rather than adding new
    // coverage: get_int_scalar's own default is what makes it work.
    // Three spheres spanning mu and gamma: two at one gamma left that axis
    // constant, which a database of several spheres is refused for.
    const std::vector<std::array<double, 4>> conditions = {
        {0.0, 0.0, 0.005, 100.0},
        {1.0, 0.0, 0.005, 110.0},
        {0.0, 1.0, 0.005, 120.0},
    };
    const auto path = write_shared_grid_test_database(
            "niceops_no_shared_flag.nod", conditions, false, false);
    Database db(path); // must not throw
    const auto sphere = db.interpolate_sphere(0.0, 0.0, 0.005);
    REQUIRE(approx_equal(sphere.cells[0].dBA, 100.0));
    std::filesystem::remove(path);
}

TEST_CASE("self_test does not flag a direction ambient gating legitimately emptied",
        "[sphere_database]")
{
    // A silent direction is stored as dBA = -inf, and the constructor floors
    // that to no_energy_floor_db so a zero barycentric weight against it
    // cannot produce a NaN (see the member comment). self_test() re-reads
    // the raw file value directly rather than going through that floor, so
    // without the same floor applied to its own comparison it always
    // disagreed with eval() -- found on Be407's real fixed_load_factor
    // database (3399 of 238428 samples), reproduced here on one channel.
    const auto path = write_test_database_with_silent_direction(
            "niceops_self_test_silent_direction.nod");
    Database db(path);
    db.self_test(path); // must not throw
    std::filesystem::remove(path);
}


// ---------------------------------------------------------------------------
// Empty directions drop out of the blends rather than silencing them.

namespace {
const double minus_inf = -std::numeric_limits<double>::infinity();
const std::vector<std::array<double, 4>> uniform_conditions = {
    {0.0, 0.0, 0.005, 100.0},
    {1.0, 0.0, 0.005, 100.0},
    {0.0, 1.0, 0.005, 100.0},
};
// Channel 5 is (phi -30, theta 60): row 1, column 1 of the 4x4 test sphere.
constexpr std::size_t gated_channel = 5;
} // namespace

TEST_CASE("a direction empty in one condition takes the others' level, not silence",
        "[sphere_database][empty]")
{
    // Blending in dB is linear, so before this a third of the weight on the
    // -1000 dB floor put the blended level near -270 dB: a measured direction
    // became silent because one neighboring condition had not measured it.
    const auto path = write_shared_grid_test_database("niceops_empty_one_condition.nod",
            uniform_conditions, false, false, 30.48,
            [](std::size_t condition, std::size_t channel) {
                return condition == 0 && channel == gated_channel ? minus_inf : 100.0;
            });
    Database db(path);
    const auto sphere = db.interpolate_sphere(0.3, 0.3, 0.005); // every condition weighted
    REQUIRE(approx_equal(sphere.cells[gated_channel].dBA, 100.0));
    REQUIRE(approx_equal(sphere.cells[0].dBA, 100.0));
    std::filesystem::remove(path);
}

TEST_CASE("a direction empty in every condition stays empty", "[sphere_database][empty]")
{
    const auto path = write_shared_grid_test_database("niceops_empty_all_conditions.nod",
            uniform_conditions, false, false, 30.48, [](std::size_t, std::size_t channel) {
                return channel == gated_channel ? minus_inf : 100.0;
            });
    Database db(path);
    const auto sphere = db.interpolate_sphere(0.3, 0.3, 0.005);
    REQUIRE(sphere.cells[gated_channel].dBA <= Database::no_energy_floor_db + 1.0);
    std::filesystem::remove(path);
}

TEST_CASE("an empty grid corner drops out of the bilinear step", "[sphere_database][empty]")
{
    // Every condition empty at one grid point, so the empty value reaches the
    // corner stage.  A direction in a cell that touches it must read the same
    // as one in a fully measured cell, the field being 100 dB all round.
    const auto path = write_shared_grid_test_database("niceops_empty_corner.nod",
            uniform_conditions, false, false, 30.48, [](std::size_t, std::size_t channel) {
                return channel == gated_channel ? minus_inf : 100.0;
            });
    Database db(path);
    const double range = 100.0;
    const double touching = db.eval(range, -45.0, 40.0, 0.3, 0.3, 0.005); // cell with (-30, 60)
    const double clear = db.eval(range, 60.0, 150.0, 0.3, 0.3, 0.005);    // (30..90, 120..180)
    REQUIRE(std::isfinite(touching));
    REQUIRE(approx_equal(touching, clear));
    std::filesystem::remove(path);
}

TEST_CASE("nothing empty: the blend is exactly what it always was", "[sphere_database][empty]")
{
    // The empty-aware path only runs when an entry is empty.  Fully measured,
    // querying a condition the database holds returns that condition's value.
    const std::vector<std::array<double, 4>> conditions = {
        {0.0, 0.0, 0.005, 90.0},
        {1.0, 0.0, 0.005, 100.0},
        {0.0, 1.0, 0.005, 110.0},
    };
    const auto path = write_shared_grid_test_database("niceops_empty_none.nod", conditions, false, false);
    Database db(path);
    const auto sphere = db.interpolate_sphere(1.0, 0.0, 0.005);
    REQUIRE(sphere.cells[gated_channel].dBA == 100.0);
    std::filesystem::remove(path);
}

TEST_CASE("spectra land in the direction their channel names, whatever the row order",
        "[sphere_database][spectrum]")
{
    // panam wrote amplitude rows in completion order (lower hemisphere, then
    // the mirrored upper surface), not sorted; read as sorted, every spectrum
    // sat in the wrong direction -- an exported Be407 sphere came back with its
    // lateral asymmetry reversed.  The channels say which row is which.
    for (bool shared : {false, true})
    {
        const auto path = write_shared_grid_test_database(
                shared ? "niceops_reversed_rows_shared.nod" : "niceops_reversed_rows.nod",
                uniform_conditions, shared, true, 30.48, {}, true);
        Database db(path);
        const auto sphere = db.interpolate_sphere(0.0, 0.0, 0.005);
        const std::size_t n_bands = sphere.frequencies.size();
        REQUIRE(n_bands == 2);
        for (std::size_t m = 0; m < sphere.phis.size(); m++)
            for (std::size_t n = 0; n < sphere.thetas.size(); n++)
                for (std::size_t b = 0; b < n_bands; b++)
                    REQUIRE(approx_equal(
                            sphere.amplitudes[(m * sphere.thetas.size() + n) * n_bands + b],
                            200.0 + sphere.phis[m] + 10.0 * static_cast<double>(b)));
        std::filesystem::remove(path);
    }
}

namespace {

// A database of `conditions` as the shared-grid writer lays it out, loaded and
// removed again; returns what loading it said was wrong, or empty.
std::string load_problem_of(const std::string& name,
        const std::vector<std::array<double, 4>>& conditions)
{
    const auto path = write_shared_grid_test_database(name, conditions, false, false);
    std::string problem;
    try
    {
        Database db(path);
    }
    catch (const std::exception& e)
    {
        problem = e.what();
    }
    std::filesystem::remove(path);
    return problem;
}

} // namespace

TEST_CASE("a database whose conditions do not vary is refused, not NaN", "[sphere_database][validation]")
{
    // (x - min) / (max - min) is 0/0 on an axis that does not vary, and the
    // triangulation was built on NaN: every query then returned one arbitrary
    // sphere (114.0 dB at four conditions whose true levels were 100 to 112),
    // silently.  The writer flags equal thrust coefficients as fixed_load_factor,
    // as a real build of reference spheres is.
    const std::vector<std::array<double, 4>> level_flight_only = {
        {0.0, 0.0, 0.004, 100.0}, {0.1, 0.0, 0.004, 105.0}, {0.2, 0.0, 0.004, 110.0},
        {0.0, 0.0, 0.008, 104.0}, {0.1, 0.0, 0.008, 109.0}, {0.2, 0.0, 0.008, 114.0}};
    const auto problem = load_problem_of("niceops_constant_gamma.nod", level_flight_only);
    REQUIRE(problem.find("flight path angle is 0 in all 6 spheres") != std::string::npos);

    const std::vector<std::array<double, 4>> constant_mu = {
        {0.0, 0.0, 0.005, 100.0}, {0.0, 1.0, 0.005, 110.0}, {0.0, 2.0, 0.005, 120.0}};
    REQUIRE(load_problem_of("niceops_constant_mu.nod", constant_mu).find("advance ratio is 0")
            != std::string::npos);

    // A varying set loads.
    const std::vector<std::array<double, 4>> fine = {
        {0.0, 0.0, 0.005, 100.0}, {1.0, 0.0, 0.005, 110.0}, {0.0, 1.0, 0.005, 120.0}};
    REQUIRE(load_problem_of("niceops_varying.nod", fine).empty());
}

TEST_CASE("two spheres at one operating condition are refused", "[sphere_database][validation]")
{
    // The second left a zero-volume simplex in the triangulation: the load
    // aborted with "could not locate a simplex", or every query walked the whole
    // mesh and fell through to the nearest sphere.
    const std::vector<std::array<double, 4>> duplicated = {
        {0.0, 0.0, 0.005, 100.0}, {1.0, 0.0, 0.005, 110.0}, {0.0, 1.0, 0.005, 120.0},
        {1.0, 1.0, 0.005, 130.0}, {1.0, 0.0, 0.005, 111.0}};
    const auto problem = load_problem_of("niceops_duplicate_condition.nod", duplicated);
    REQUIRE(problem.find("Spheres 1 and 4 are at the same operating condition")
            != std::string::npos);
}

TEST_CASE("inside the sphere the level is the surface's", "[sphere_database]")
{
    // Spherical spreading extrapolated inward grows as 1/r and diverges at r = 0:
    // an observer exactly at a source position got +inf in the CSV, an observer
    // within the 30 m sphere up to 4.4 dB more than at its surface, and with the
    // asserts live the r = 0 case aborted the run.  The sphere is the source's own
    // volume, so spreading and absorption start at its surface and the level inside
    // is the surface's.
    if (!std::filesystem::exists("databases/S-76D_M3.nod"))
        return;
    Database db("databases/S-76D_M3.nod");
    const double radius = db.sphere_radius();
    const auto blend = db.blend_at(0.1, 0.0, 0.006);
    const double at_surface = db.eval(radius, 20.0, 60.0, blend);
    for (double fraction : {0.0, 1e-12, 1e-6, 0.01, 0.1, 0.5, 0.99, 1.0})
    {
        const double level = db.eval(fraction * radius, 20.0, 60.0, blend);
        INFO("range " << fraction << " radii");
        REQUIRE(std::isfinite(level));
        REQUIRE(level == at_surface);
    }
    // Beyond the surface the level still falls with range, as it always did.
    double previous = at_surface;
    for (double fraction : {1.01, 1.5, 2.0, 10.0, 100.0})
    {
        const double level = db.eval(fraction * radius, 20.0, 60.0, blend);
        INFO("range " << fraction << " radii");
        REQUIRE(level < previous);
        previous = level;
    }
}

TEST_CASE("a database stored as float32 loads as what it holds", "[sphere_database][validation]")
{
    // nc_get_var copies the file's external type unconverted, so a float32 .nod
    // was read as garbage into doubles: a sphere radius of 0, a maximum level of
    // 3e14, NaN footprints, and an exit status of 0.
    const auto path = write_fixed_load_factor_test_database(
            "niceops_float32.nod", fixed_load_factor_conditions, true, 30.0, NC_FLOAT);
    const Database database(path);
    CHECK(database.sphere_radius() == 30.0);
    // The conditions' own levels, from three different spheres.
    for (const auto& condition : fixed_load_factor_conditions)
    {
        const auto blend = database.blend_at(condition[0], condition[1], condition[2]);
        CHECK(database.eval(30.0, 0.0, 90.0, blend) == Catch::Approx(condition[3]).margin(1e-9));
    }
}

TEST_CASE("a sphere whose channels miss the grid is reported", "[sphere_database][validation]")
{
    // The axes come from sphere 0; a sphere whose phi labels are two degrees off
    // matches none of them, so its channels are dropped and its cells stay 0 dBA.
    // Nothing used to say so.
    const auto path = (std::filesystem::temp_directory_path() / "niceops_coverage.log").string();
    std::filesystem::remove(path);
    Logger::instance().set_file(path);
    const auto database_path = write_fixed_load_factor_test_database(
            "niceops_offgrid.nod", fixed_load_factor_conditions, true, 30.0, NC_DOUBLE, 2.0);
    const Database database(database_path);
    const Database aligned(write_fixed_load_factor_test_database(
            "niceops_ongrid.nod", fixed_load_factor_conditions, true, 30.0));
    Logger::instance().close_file();
    std::ifstream in(path);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    CHECK(text.find("sphere 1 leaves 16 of 16 (phi, theta) cells without a channel") != std::string::npos);
    CHECK(text.find("drops 16 of 16 channels") != std::string::npos);
    // The aligned twin says nothing of the sort: one warning, not two.
    CHECK(text.find("without a channel") == text.rfind("without a channel"));
}

namespace {
constexpr double netcdf_fill_value = 9.9692099683868690e+36;
} // namespace

TEST_CASE("a group holding the netCDF fill value is refused", "[sphere_database][validation][screen]")
{
    // A group panam never finished writing holds the fill value everywhere.  It
    // is finite, so the non-finite floor passed it, and it read as a 1e37 dB
    // level: an infinite SEL wherever a state reached that condition.
    const auto path = write_shared_grid_test_database("niceops_fill_value.nod", uniform_conditions,
            false, false, 30.48,
            [](std::size_t condition, std::size_t) { return condition == 2 ? netcdf_fill_value : 100.0; });
    std::string problem;
    try
    {
        Database db(path);
    }
    catch (const std::exception& e)
    {
        problem = e.what();
    }
    std::filesystem::remove(path);
    REQUIRE(problem.find("sphere 2 holds the netCDF fill value in 16 of 16 dBA entries") != std::string::npos);
}

TEST_CASE("a dBA cell far from the rest of its sphere is reported, and dropped when asked",
        "[sphere_database][screen]")
{
    // The S-76D M4 and fixed_load_factor files carry channels up to 1450 dBA
    // beside 90-100 dBA neighbors, which the blends took as real.
    const auto path = write_shared_grid_test_database("niceops_spike.nod", uniform_conditions, false,
            false, 30.48, [](std::size_t condition, std::size_t channel) {
                return condition == 0 && channel == gated_channel ? 1450.0 : 100.0;
            });
    const auto clean_path = write_shared_grid_test_database("niceops_spike_clean.nod",
            uniform_conditions, false, false);
    const Database kept(path);
    DatabaseLoadOptions mask;
    mask.mask_suspect_cells = true;
    const Database masked(path, mask);
    const Database clean(clean_path);
    CHECK(kept.suspect_cells() == 1);
    CHECK(kept.suspect_bands() == 0);
    CHECK(masked.suspect_cells() == 1);
    CHECK(clean.suspect_cells() == 0);
    // Kept, the spike is blended in; masked, it reads as a gated direction,
    // which takes the other conditions' level.
    CHECK(kept.interpolate_sphere(0.3, 0.3, 0.005).cells[gated_channel].dBA > 200.0);
    CHECK(approx_equal(masked.interpolate_sphere(0.3, 0.3, 0.005).cells[gated_channel].dBA, 100.0));
    // The per-sphere maximum leaves with it, or the cull and the bounds would
    // keep every observer in range of a level no direction holds.
    const double cutoff = 90.0;
    const double clean_range = clean.cull_range(clean.blend_at(0.0, 0.0, 0.005), cutoff);
    CHECK(masked.cull_range(masked.blend_at(0.0, 0.0, 0.005), cutoff) == clean_range);
    CHECK(kept.cull_range(kept.blend_at(0.0, 0.0, 0.005), cutoff) > clean_range);
    std::filesystem::remove(path);
    std::filesystem::remove(clean_path);
}

TEST_CASE("a near-silent cell is reported, and a steep but sound sphere is not",
        "[sphere_database][screen]")
{
    // B407_M3's written zeros and the September files' -70 dBA cells beside
    // 84 dBA neighbors are as wrong as a spike, the other way.
    const auto zero = write_shared_grid_test_database("niceops_zero_cell.nod", uniform_conditions,
            false, false, 30.48, [](std::size_t condition, std::size_t channel) {
                return condition == 1 && channel == gated_channel ? 0.0 : 100.0;
            });
    // 12 dB per row, 36 across the sphere: steeper than most real directivity,
    // well inside the margin.
    const auto steep = write_shared_grid_test_database("niceops_steep.nod", uniform_conditions,
            false, false, 30.48, [](std::size_t, std::size_t channel) {
                return 60.0 + 12.0 * static_cast<double>(channel / 4);
            });
    CHECK(Database(zero).suspect_cells() == 1);
    CHECK(Database(steep).suspect_cells() == 0);
    std::filesystem::remove(zero);
    std::filesystem::remove(steep);
}

TEST_CASE("a sound cell beside a strip of zeros is not flagged with it", "[sphere_database][screen]")
{
    // Channels 7 to 11 are zeros (row 2, and the end of row 1).  Against its
    // grid neighbors, half of them zeros, channel 6 read 50 dB off and the
    // whole of row 3 was cut off behind the strip: next to B407_M3's zero
    // blocks, 224 sound 88 dBA cells were flagged that way.  Against the
    // sphere's median, only the zeros are.
    const auto path = write_shared_grid_test_database("niceops_zero_strip.nod", uniform_conditions,
            false, false, 30.48, [](std::size_t condition, std::size_t channel) {
                return condition == 0 && channel >= 7 && channel <= 11 ? 0.0 : 100.0;
            });
    const auto list = (std::filesystem::temp_directory_path() / "niceops_zero_strip_screen.csv").string();
    DatabaseLoadOptions load;
    load.screen_file = list;
    CHECK(Database(path, load).suspect_cells() == 5);
    std::ifstream in(list);
    std::string line;
    std::getline(in, line);
    std::getline(in, line);
    std::size_t rows = 0;
    while (std::getline(in, line))
    {
        rows++;
        CHECK(line.find(",,0,100,0") != std::string::npos); // a zero, against the median 100
    }
    CHECK(rows == 5);
    std::filesystem::remove(path);
    std::filesystem::remove(list);
}

TEST_CASE("a dropped dBA cell drops its spectrum too", "[sphere_database][screen][spectrum]")
{
    // The spectral path reads only the bands.  A spectrum lowered as a whole has
    // no band out of line with its own cell, so only the dBA screen sees it,
    // and masking the dBA alone left the default path blending it.
    const std::vector<double> bands = {500.0, 1000.0, 2000.0, 4000.0, 8000.0, 10000.0};
    const auto path = write_shared_grid_test_database("niceops_quiet_cell.nod", uniform_conditions,
            false, true, 30.48,
            [](std::size_t condition, std::size_t channel) {
                return condition == 1 && channel == gated_channel ? 10.0 : 100.0;
            },
            false, bands, [](std::size_t condition, std::size_t cell, std::size_t) {
                return condition == 1 && cell == gated_channel ? 5.0 : 90.0;
            });
    DatabaseLoadOptions mask;
    mask.mask_suspect_cells = true;
    const Database masked(path, mask);
    CHECK(masked.suspect_cells() == 1);
    CHECK(masked.suspect_bands() == 0);
    const auto sphere = masked.interpolate_sphere(1.0, 0.0, 0.005);
    for (std::size_t band = 0; band < bands.size(); band++)
        CHECK(sphere.amplitudes[gated_channel * bands.size() + band] <= Database::no_energy_floor_db + 1.0);
    CHECK(sphere.amplitudes[(gated_channel + 1) * bands.size()] == 90.0);
    std::filesystem::remove(path);
}

TEST_CASE("a band far above its cell's spectrum is reported, listed and dropped when asked",
        "[sphere_database][screen][spectrum]")
{
    // NASA's R44 spheres hold 135-170 dB at 8 and 10 kHz in cells whose other
    // bands are 64-82 dB.  A band 40 dB above the rest is inside the margin.
    const std::vector<double> bands = {500.0, 1000.0, 2000.0, 4000.0, 8000.0, 10000.0};
    const auto path = write_shared_grid_test_database("niceops_hot_band.nod", uniform_conditions,
            false, true, 30.48, {}, false, bands,
            [](std::size_t condition, std::size_t cell, std::size_t band) {
                if (condition == 1 && cell == gated_channel && band == 5)
                    return 170.0;
                if (condition == 1 && cell == 6 && band == 0)
                    return 120.0;
                return 80.0;
            });
    const auto list = (std::filesystem::temp_directory_path() / "niceops_hot_band_screen.csv").string();
    std::filesystem::remove(list);
    const Database kept(path);
    CHECK(kept.suspect_bands() == 1);
    CHECK(kept.suspect_cells() == 0);
    CHECK(kept.interpolate_sphere(1.0, 0.0, 0.005).amplitudes[gated_channel * bands.size() + 5] == 170.0);
    DatabaseLoadOptions mask;
    mask.mask_suspect_cells = true;
    mask.screen_file = list;
    const Database masked(path, mask);
    const auto sphere = masked.interpolate_sphere(1.0, 0.0, 0.005);
    CHECK(sphere.amplitudes[gated_channel * bands.size() + 5] <= Database::no_energy_floor_db + 1.0);
    CHECK(sphere.amplitudes[gated_channel * bands.size() + 4] == 80.0);
    CHECK(sphere.amplitudes[6 * bands.size()] == 120.0);
    std::ifstream in(list);
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line);)
        lines.push_back(line);
    REQUIRE(lines.size() == 3);
    CHECK(lines[1] == "sphere,advance_ratio,flight_path_angle,thrust_coefficient,phi,theta,band_hz,"
                      "level_db,reference_db,masked");
    CHECK(lines[2] == "1,1,0,0.005,-30,60,10000,170,80,1");
    std::filesystem::remove(path);
    std::filesystem::remove(list);
}

TEST_CASE("a frame with no velocity is still a frame", "[sphere_database][validation]")
{
    // Zero velocity, or velocity along the thrust axis, left the nose zero and with it
    // the starboard vector: every observer was looked up at (phi 0, theta 90), the same
    // point of the sphere from every side, and the answer jumped when the velocity
    // became a millimeter a second.
    for (const Vec3 velocity : {Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, -2.0), Vec3(0.0, 0.0, 2.0)})
    {
        const Vec3 normal(0.0, 0.0, 1.0);
        const auto frame = sphere_frame(velocity, normal);
        CHECK(frame.nose.norm() == Catch::Approx(1.0).margin(1e-12));
        CHECK(frame.starboard.norm() == Catch::Approx(1.0).margin(1e-12));
        CHECK(frame.nose.dot(frame.normal) == Catch::Approx(0.0).margin(1e-12));
        CHECK(frame.nose.dot(frame.starboard) == Catch::Approx(0.0).margin(1e-12));
        // And directions from different sides are looked up at different places.
        const auto ahead = lookup_angles(frame, Vec3(0.6, 0.0, -0.8));
        const auto beside = lookup_angles(frame, Vec3(0.0, 0.6, -0.8));
        CHECK(!(ahead[0] == beside[0] && ahead[1] == beside[1]));
    }
    // A velocity with a component in the plane is unchanged: the nose is that
    // component's direction.
    const auto moving = sphere_frame(Vec3(3.0, 4.0, 5.0), Vec3(0.0, 0.0, 1.0));
    CHECK(moving.nose.x == Catch::Approx(0.6).margin(1e-12));
    CHECK(moving.nose.y == Catch::Approx(0.8).margin(1e-12));
}

TEST_CASE("an exported sphere reloads as the sphere it was, cell for cell", "[sphere_database][export]")
{
    // Every level and EAA of the exported sphere, and what the root says about the
    // vehicle and its speed reference, comes back exactly.  Swapped dBA and EAA, a
    // dropped load-factor correction or a transposed grid all passed the whole
    // suite before this.
    Database source("databases/S-76D_M3.nod");
    const auto sphere = source.interpolate_sphere(0.15, -4.0, 0.006);
    const auto path = (std::filesystem::temp_directory_path() / "niceops_round_trip.nod").string();
    write_sphere(sphere, path);
    const Database reloaded(path);
    const auto again = reloaded.interpolate_sphere(0.15, -4.0, 0.006);
    REQUIRE(again.phis == sphere.phis);
    REQUIRE(again.thetas == sphere.thetas);
    REQUIRE(again.cells.size() == sphere.cells.size());
    double eaa_sum = 0.0;
    for (std::size_t c = 0; c < sphere.cells.size(); c++)
    {
        REQUIRE(again.cells[c].dBA == sphere.cells[c].dBA);
        REQUIRE(again.cells[c].EAA == sphere.cells[c].EAA);
        eaa_sum += std::fabs(sphere.cells[c].EAA);
    }
    // The test means something only if the two differ: dBA and EAA are not the same field.
    REQUIRE(eaa_sum > 0.0);
    REQUIRE(sphere.cells[0].dBA != sphere.cells[0].EAA);
    // What the root records.
    CHECK(reloaded.speed_reference() == source.speed_reference());
    if (!std::isnan(source.tip_speed()))
        CHECK(reloaded.tip_speed() == source.tip_speed());
    std::filesystem::remove(path);
}

TEST_CASE("an exported sphere keeps its speed reference and rotor", "[sphere_database][export]")
{
    // The loader decides how to read a run's flight conditions from the root's
    // speed_reference, and checks the vehicle file against the rotor it records; an
    // export that dropped both reloaded as ground-referenced with no rotor, and the
    // run neither refused nor referred to the air.
    const auto source_path = write_fixed_load_factor_test_database(
            "niceops_air_source.nod", fixed_load_factor_conditions, true);
    int ncid = 0;
    REQUIRE(nc_open(source_path.c_str(), NC_WRITE, &ncid) == NC_NOERR);
    REQUIRE(nc_redef(ncid) == NC_NOERR);
    REQUIRE(nc_put_att_text(ncid, NC_GLOBAL, "speed_reference", 3, "air") == NC_NOERR);
    int scalar_dim = 0;
    REQUIRE(nc_inq_dimid(ncid, "scalar", &scalar_dim) == NC_NOERR);
    int tip_var = 0, radius_var = 0, weight_var = 0;
    REQUIRE(nc_def_var(ncid, "main_rotor_tip_speed_meters_per_sec", NC_DOUBLE, 1, &scalar_dim, &tip_var) == NC_NOERR);
    REQUIRE(nc_def_var(ncid, "main_rotor_radius_meters", NC_DOUBLE, 1, &scalar_dim, &radius_var) == NC_NOERR);
    REQUIRE(nc_def_var(ncid, "vehicle_weight_newtons", NC_DOUBLE, 1, &scalar_dim, &weight_var) == NC_NOERR);
    REQUIRE(nc_enddef(ncid) == NC_NOERR);
    const double tip = 213.4, radius = 6.7, weight = 49033.2;
    REQUIRE(nc_put_var_double(ncid, tip_var, &tip) == NC_NOERR);
    REQUIRE(nc_put_var_double(ncid, radius_var, &radius) == NC_NOERR);
    REQUIRE(nc_put_var_double(ncid, weight_var, &weight) == NC_NOERR);
    REQUIRE(nc_close(ncid) == NC_NOERR);

    const Database source(source_path);
    REQUIRE(source.speed_reference() == SpeedReference::air);
    // Above the reference thrust coefficient, so the load-factor correction is in
    // the exported levels as well.
    const auto sphere = source.interpolate_sphere(0.2, 0.2, 2.0 * test_ct_reference);
    const auto path = (std::filesystem::temp_directory_path() / "niceops_air_export.nod").string();
    write_sphere(sphere, path);
    const Database reloaded(path);
    CHECK(reloaded.speed_reference() == SpeedReference::air);
    CHECK(reloaded.tip_speed() == tip);
    CHECK(reloaded.rotor_radius() == radius);
    CHECK(reloaded.vehicle_weight() == weight);
    // The correction survived: 20 log10 of the load factor (2), on levels of 100-120 dB.
    const auto again = reloaded.interpolate_sphere(0.2, 0.2, 2.0 * test_ct_reference);
    for (std::size_t c = 0; c < sphere.cells.size(); c++)
        REQUIRE(again.cells[c].dBA == Catch::Approx(sphere.cells[c].dBA).margin(1e-9));
    const auto uncorrected = source.interpolate_sphere(0.2, 0.2, test_ct_reference);
    CHECK(sphere.cells[0].dBA - uncorrected.cells[0].dBA == Catch::Approx(20.0 * std::log10(2.0)).margin(1e-9));
    std::filesystem::remove(path);
    std::filesystem::remove(source_path);
}

TEST_CASE("an AAM export lays out its axes and amplitudes as the format says", "[sphere_database][export]")
{
    // AMPLITUDE(PHI, THETA, FREQUENCY), the axes in degrees and hertz, RADIUS in feet
    // and SPEED in knots from the advance ratio and the tip speed.  Nothing read an
    // AAM file back before this.
    Database source("databases/S-76D_M3.nod");
    auto sphere = source.interpolate_sphere(0.15, -4.0, 0.006);
    const std::size_t bands = 3;
    sphere.frequencies = {125.0, 250.0, 500.0};
    sphere.amplitudes.resize(sphere.cells.size() * bands);
    for (std::size_t cell = 0; cell < sphere.cells.size(); cell++)
        for (std::size_t band = 0; band < bands; band++)
            sphere.amplitudes[cell * bands + band] = 40.0 + 0.5 * static_cast<double>(cell % 97) + 10.0 * band;
    const auto path = (std::filesystem::temp_directory_path() / "niceops_aam_export.nc").string();
    const double tip_speed = 213.0;
    write_aam_sphere(sphere, path, tip_speed);

    int ncid = 0;
    REQUIRE(nc_open(path.c_str(), NC_NOWRITE, &ncid) == NC_NOERR);
    auto length = [&](const char* name) {
        int dim = 0;
        std::size_t n = 0;
        REQUIRE(nc_inq_dimid(ncid, name, &dim) == NC_NOERR);
        REQUIRE(nc_inq_dimlen(ncid, dim, &n) == NC_NOERR);
        return n;
    };
    REQUIRE(length("PHI") == sphere.phis.size());
    REQUIRE(length("THETA") == sphere.thetas.size());
    REQUIRE(length("FREQUENCY") == bands);
    auto read_floats = [&](const char* name, std::size_t n) {
        int var = 0;
        REQUIRE(nc_inq_varid(ncid, name, &var) == NC_NOERR);
        std::vector<float> values(n);
        REQUIRE(nc_get_var_float(ncid, var, values.data()) == NC_NOERR);
        return values;
    };
    const auto phi = read_floats("PHI", sphere.phis.size());
    const auto theta = read_floats("THETA", sphere.thetas.size());
    const auto frequency = read_floats("FREQUENCY", bands);
    for (std::size_t i = 0; i < phi.size(); i++)
        CHECK(phi[i] == static_cast<float>(sphere.phis[i]));
    for (std::size_t i = 0; i < theta.size(); i++)
        CHECK(theta[i] == static_cast<float>(sphere.thetas[i]));
    CHECK(frequency[1] == 250.0f);
    // Element order: phi slowest, then theta, then frequency.
    const auto amplitude = read_floats("AMPLITUDE", sphere.cells.size() * bands);
    for (std::size_t m = 0; m < phi.size(); m += 5)
        for (std::size_t n = 0; n < theta.size(); n += 3)
            for (std::size_t f = 0; f < bands; f++)
                REQUIRE(amplitude[(m * theta.size() + n) * bands + f]
                        == static_cast<float>(sphere.amplitudes[(m * theta.size() + n) * bands + f]));
    CHECK(read_floats("RADIUS", 1)[0] == Catch::Approx(sphere.radius / 0.3048).epsilon(1e-6));
    CHECK(read_floats("SPEED", 1)[0] == Catch::Approx(0.15 * tip_speed / knot).epsilon(1e-6));
    CHECK(read_floats("FLIGHT_PATH_ANGLE", 1)[0] == Catch::Approx(-4.0).margin(1e-6));
    REQUIRE(nc_close(ncid) == NC_NOERR);
    std::filesystem::remove(path);
}

TEST_CASE("a negative excess-attenuation factor turns the bounds off", "[sphere_database][bounds]")
{
    // cull_range() and direction_bounds() are upper bounds only because they drop
    // the excess-attenuation term, which is a loss only if no EAA is negative.  The
    // guard is what keeps the windowed run exact on a database that has one; nothing
    // built a database with one, so removing the guard passed every test.
    std::vector<Database::SourceState> states;
    for (const double first_eaa : {0.0, -20.0})
    {
        const auto path = write_fixed_load_factor_test_database("niceops_negative_eaa.nod",
                fixed_load_factor_conditions, true, 30.48, NC_DOUBLE, 0.0, first_eaa);
        const Database database(path);
        states.clear();
        states.push_back(database.source_state(Vec3(0.0, 0.0, 100.0), Vec3(50.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0),
                0.0, 0.0, test_ct_reference));
        const auto blend = database.blend_at(0.0, 0.0, test_ct_reference);
        INFO("first channel's EAA " << first_eaa);
        if (first_eaa >= 0.0)
        {
            CHECK(database.direction_bounds(states).usable());
            CHECK(std::isfinite(database.cull_range(blend, 60.0)));
        }
        else
        {
            CHECK_FALSE(database.direction_bounds(states).usable());
            CHECK(std::isinf(database.cull_range(blend, 60.0)));
        }
        std::filesystem::remove(path);
    }
}

TEST_CASE("phi reads the upper half above the plane, and folds only when asked",
        "[sphere_database][upper_hemisphere]")
{
    // Directions built from the frame's own axes, so that up and across are
    // exactly what the lookup sees.  psi sweeps the circle about the nose.
    const auto frame = sphere_frame(Vec3(1.0, 0.0, 0.0), Vec3(0.0, 0.0, 1.0));
    const double to_degrees = 180.0 / pi;
    for (const double psi_degrees : {-150.0, -60.0, -20.0, 0.0, 45.0, 170.0})
    {
        const double psi = psi_degrees * pi / 180.0;
        const Vec3 direction = 0.6 * frame.nose + (0.8 * std::cos(psi)) * frame.normal
                + (0.8 * std::sin(psi)) * frame.starboard;
        const double up = frame.normal.dot(direction);
        const double across = frame.starboard.dot(direction);
        INFO("psi = " << psi_degrees);
        const auto full = lookup_angles(frame, direction);
        const auto folded = lookup_angles(frame, direction, true);
        // The full circle: atan2 of (across, -up), the stored layout's own
        // convention (phi = 0 straight down, 180 straight up).
        CHECK(full[0] == Catch::Approx(std::atan2(across, -up) * to_degrees).margin(1e-9));
        // The fold: the mirror image, atan(across / |up|), as before.
        CHECK(folded[0] == Catch::Approx(std::atan(across / std::fabs(up)) * to_degrees).margin(1e-9));
        CHECK(full[1] == folded[1]);
        if (up < 0.0)
            CHECK(full[0] == folded[0]);
        else
            CHECK(std::fabs(full[0]) > 90.0);
    }
}

TEST_CASE("GridAxis::locate crosses the seam of a periodic axis", "[sphere_database][upper_hemisphere]")
{
    // The full-circle phi axis runs -180 to 170 in 10-degree steps, so the one
    // cell that is not a cell of the axis itself is the seam, 170 to 180.
    std::vector<double> phi;
    for (int i = -180; i < 180; i += 10)
        phi.push_back(static_cast<double>(i));
    GridAxis axis;
    axis.reset(phi);
    axis.set_period(360.0);

    // Every value on the circle, seam included, weighs its two corners
    // exactly as the periodic reference does.
    for (int step = -720; step <= 720; step++)
    {
        const double val = 0.25 * static_cast<double>(step);
        const auto cell = axis.locate(val);
        CHECK(cell.fraction >= 0.0);
        CHECK(cell.fraction <= 1.0);
        // 180 is the same direction as -180, so the reference reads it there,
        // as locate() does.
        const double x = val >= 180.0 ? val - 360.0 : val;
        const double reference = periodic_bilinear_interp(0.0, 0.0, 1.0, 1.0,
                phi[cell.lo], phi[cell.hi], 0.0, 1.0, x, 0.5, 360.0);
        REQUIRE(std::fabs(cell.fraction - reference) <= 1e-12);
    }

    // Halfway across the seam sits halfway between the 170 knot and the -180
    // knot one period on, so the weight on the wrapped corner is one half.
    const auto seam = axis.locate(175.0);
    CHECK(seam.lo == phi.size() - 1);
    CHECK(seam.hi == 0);
    CHECK(seam.fraction == Catch::Approx(0.5));
    // 180 is the same direction as -180, the first knot.
    const auto pole = axis.locate(180.0);
    CHECK(phi[pole.lo] + pole.fraction * (phi[pole.hi] - phi[pole.lo]) == Catch::Approx(-180.0).margin(1e-9));
}

TEST_CASE("The build atmosphere and rotor speed are read when the file records them",
        "[sphere_database][build]")
{
    // Older databases carry neither, and must say nothing rather than guess.
    const Database plain("databases/S-76D_M3.nod");
    REQUIRE(!plain.build_atmosphere().has_value());
    REQUIRE(!plain.rotor_rpm().has_value());

    const auto sphere = plain.interpolate_sphere(0.15, -4.0, 0.006);
    const auto path = (std::filesystem::temp_directory_path() / "niceops_build_attributes.nod").string();
    write_sphere(sphere, path);
    int ncid = 0;
    REQUIRE(nc_open(path.c_str(), NC_WRITE, &ncid) == NC_NOERR);
    REQUIRE(nc_redef(ncid) == NC_NOERR);
    const double temperature = 293.15;
    const double pressure = 101.325;
    const double humidity = 20.0;
    const double rpm = 393.0;
    REQUIRE(nc_put_att_double(ncid, NC_GLOBAL, "build_temperature_K", NC_DOUBLE, 1, &temperature)
            == NC_NOERR);
    REQUIRE(nc_put_att_double(ncid, NC_GLOBAL, "build_pressure_kPa", NC_DOUBLE, 1, &pressure)
            == NC_NOERR);
    REQUIRE(nc_put_att_double(ncid, NC_GLOBAL, "build_relative_humidity_percent", NC_DOUBLE, 1,
            &humidity) == NC_NOERR);
    REQUIRE(nc_put_att_double(ncid, NC_GLOBAL, "rotor_rpm", NC_DOUBLE, 1, &rpm) == NC_NOERR);
    REQUIRE(nc_close(ncid) == NC_NOERR);

    const Database recorded(path);
    REQUIRE(recorded.build_atmosphere().has_value());
    CHECK(recorded.build_atmosphere()->temperature_kelvin == Catch::Approx(293.15));
    CHECK(recorded.build_atmosphere()->pressure_kpa == Catch::Approx(101.325));
    CHECK(recorded.build_atmosphere()->relative_humidity_percent == Catch::Approx(20.0));
    REQUIRE(recorded.rotor_rpm().has_value());
    CHECK(*recorded.rotor_rpm() == Catch::Approx(393.0));
    std::filesystem::remove(path);
}

TEST_CASE("The unmeasured share follows the channels the builder did not measure",
        "[sphere_database][coverage]")
{
    // The shipped databases record no coverage, so the share is undefined.
    const Database plain("databases/S-76D_M3.nod");
    REQUIRE(!plain.has_coverage());

    // One sphere written out, with a coverage variable along its channels:
    // measured where phi is negative, unmeasured elsewhere.
    const auto sphere = plain.interpolate_sphere(0.15, -4.0, 0.006);
    const auto path = (std::filesystem::temp_directory_path() / "niceops_coverage.nod").string();
    write_sphere(sphere, path);
    int ncid = 0;
    REQUIRE(nc_open(path.c_str(), NC_WRITE, &ncid) == NC_NOERR);
    int group = 0;
    REQUIRE(nc_inq_grp_ncid(ncid, "sphere0", &group) == NC_NOERR);
    int channels_dim = 0;
    REQUIRE(nc_inq_dimid(group, "channels", &channels_dim) == NC_NOERR);
    std::size_t channels = 0;
    REQUIRE(nc_inq_dimlen(group, channels_dim, &channels) == NC_NOERR);
    int phi_var = 0;
    REQUIRE(nc_inq_varid(group, "phi", &phi_var) == NC_NOERR);
    std::vector<double> phi(channels);
    REQUIRE(nc_get_var_double(group, phi_var, phi.data()) == NC_NOERR);
    std::vector<int> coverage(channels);
    for (std::size_t c = 0; c < channels; c++)
        coverage[c] = phi[c] < 0.0 ? 1 : 0;
    REQUIRE(nc_redef(ncid) == NC_NOERR);
    int coverage_var = 0;
    REQUIRE(nc_def_var(group, "coverage", NC_INT, 1, &channels_dim, &coverage_var) == NC_NOERR);
    REQUIRE(nc_enddef(ncid) == NC_NOERR);
    REQUIRE(nc_put_var_int(group, coverage_var, coverage.data()) == NC_NOERR);
    REQUIRE(nc_close(ncid) == NC_NOERR);

    const Database recorded(path);
    REQUIRE(recorded.has_coverage());
    const auto state = recorded.source_state(Vec3(0.0, 0.0, 150.0), Vec3(60.0, 10.0, 2.0),
            Vec3(0.1, -0.05, 1.0), 0.1, 0.0, 0.006);
    // Both corner columns measured: nothing unmeasured.
    CHECK(recorded.unmeasured_share({-45.0, 60.0}, state.blend) == Catch::Approx(0.0).margin(1e-12));
    // Both unmeasured: all of it.
    CHECK(recorded.unmeasured_share({45.0, 60.0}, state.blend) == Catch::Approx(1.0).margin(1e-12));
    // Halfway between a measured column (-10) and an unmeasured one (0): half.
    CHECK(recorded.unmeasured_share({-5.0, 60.0}, state.blend) == Catch::Approx(0.5).margin(1e-9));
    std::filesystem::remove(path);
}

TEST_CASE("Grid-shaped coverage counts an unmeasured band by its A-weighted energy",
        "[sphere_database][coverage]")
{
    // A spectral sphere, with coverage laid out as its amplitude is:
    // (PHI, THETA, frequency), every band at the same level.  Above the plane
    // one band is unmeasured, so those cells are unmeasured by that band's
    // share of the A-weighted energy: little for 100 Hz, most for 800 Hz.
    const Database plain("databases/S-76D_M3.nod");
    auto sphere = plain.interpolate_sphere(0.15, -4.0, 0.006);
    const std::size_t bands = 4;
    sphere.frequencies = {100.0, 200.0, 400.0, 800.0};
    sphere.amplitudes.assign(sphere.cells.size() * bands, 1.0);
    const auto gains = a_weighting_gains(sphere.frequencies, BandLimits {});
    const double total = gains[0] + gains[1] + gains[2] + gains[3];
    for (const std::size_t unmeasured_band : {std::size_t {0}, std::size_t {3}})
    {
        const auto path = (std::filesystem::temp_directory_path() / "niceops_grid_coverage.nod").string();
        write_sphere(sphere, path);

        int ncid = 0;
        REQUIRE(nc_open(path.c_str(), NC_WRITE, &ncid) == NC_NOERR);
        int group = 0;
        REQUIRE(nc_inq_grp_ncid(ncid, "sphere0", &group) == NC_NOERR);
        // Dimensions for the coverage, in the group: the sphere's own axes, and
        // the frequency dimension write_sphere already defined.
        REQUIRE(nc_redef(ncid) == NC_NOERR);
        int dims[3] = {0, 0, 0};
        REQUIRE(nc_def_dim(group, "coverage_phi", sphere.phis.size(), &dims[0]) == NC_NOERR);
        REQUIRE(nc_def_dim(group, "coverage_theta", sphere.thetas.size(), &dims[1]) == NC_NOERR);
        REQUIRE(nc_inq_dimid(group, "frequency", &dims[2]) == NC_NOERR);
        const std::size_t n_phi = sphere.phis.size();
        const std::size_t n_theta = sphere.thetas.size();
        const auto& phi = sphere.phis;
        std::vector<int> coverage(n_phi * n_theta * bands);
        for (std::size_t p = 0; p < n_phi; p++)
            for (std::size_t t = 0; t < n_theta; t++)
                for (std::size_t b = 0; b < bands; b++)
                    coverage[(p * n_theta + t) * bands + b] = (phi[p] < 0.0 || b != unmeasured_band) ? 1 : 0;
        int coverage_var = 0;
        REQUIRE(nc_def_var(group, "coverage", NC_INT, 3, dims, &coverage_var) == NC_NOERR);
        REQUIRE(nc_enddef(ncid) == NC_NOERR);
        REQUIRE(nc_put_var_int(group, coverage_var, coverage.data()) == NC_NOERR);
        REQUIRE(nc_close(ncid) == NC_NOERR);

        const Database recorded(path);
        REQUIRE(recorded.has_spectra());
        REQUIRE(recorded.has_coverage());
        const auto state = recorded.source_state(Vec3(0.0, 0.0, 150.0), Vec3(60.0, 10.0, 2.0),
                Vec3(0.1, -0.05, 1.0), 0.1, 0.0, 0.006);
        const double share = gains[unmeasured_band] / total;
        INFO("unmeasured band " << sphere.frequencies[unmeasured_band] << " Hz, share " << share);
        CHECK(recorded.unmeasured_share({-45.0, 60.0}, state.blend) == Catch::Approx(0.0).margin(1e-12));
        CHECK(recorded.unmeasured_share({45.0, 60.0}, state.blend) == Catch::Approx(share).margin(1e-6));
        CHECK(recorded.unmeasured_share({-5.0, 60.0}, state.blend) == Catch::Approx(0.5 * share).margin(1e-6));
        std::filesystem::remove(path);
    }
    // The two bands bracket the cases: 100 Hz carries about 1%, 800 Hz about two thirds.
    CHECK(gains[0] / total < 0.02);
    CHECK(gains[3] / total > 0.6);
}

// --blend and --load_factor_correction through Database::blend_at.  The
// fixtures give each sphere one constant level, so the level at the sphere's
// surface (no spreading, EAA 0) says which spheres a blend took and by how
// much.  The search rules themselves are swept on synthetic condition sets in
// test_condition_selection.cpp.
namespace {

double surface_level(const Database& db, const Database::Blend& blend)
{
    return db.eval(db.sphere_radius(), 0.0, 90.0, blend);
}

// A tip speed of 100 kt, so an advance ratio is knots / 100.
constexpr double tip_speed_100_kt = 100.0 * knot;

} // namespace

TEST_CASE("--blend nearest takes the stored condition nearest in the normalized space",
        "[sphere_database][blend]")
{
    // mu spans 0.3 and gamma 20 degrees: in raw units gamma would decide every
    // distance, which is why the search runs where the triangulation does.
    const std::vector<std::array<double, 4>> conditions = {
        {0.0, -10.0, test_ct_reference, 100.0},
        {0.3, -10.0, test_ct_reference, 110.0},
        {0.0, 10.0, test_ct_reference, 120.0},
    };
    const auto path = write_fixed_load_factor_test_database("niceops_blend_nearest.nod", conditions, true);
    const Database db(path);
    Database::BlendOptions nearest;
    nearest.policy = BlendPolicy::nearest;

    // (0.2, 1): normalized (0.67, 0.55), nearest sphere 1; in raw units
    // sphere 2 (gamma 9 degrees away against 11) would have been.
    const auto blend = db.blend_at(0.2, 1.0, test_ct_reference, nearest);
    REQUIRE(blend.weights[0] == 1.0);
    for (std::size_t k = 1; k < blend.weights.size(); k++)
    {
        REQUIRE(blend.weights[k] == 0.0);
        REQUIRE(blend.indices[k] == blend.indices[0]);
    }
    REQUIRE(surface_level(db, blend) == 110.0);
    REQUIRE(surface_level(db, db.blend_at(0.05, -8.0, test_ct_reference, nearest)) == 100.0);
    REQUIRE(surface_level(db, db.blend_at(0.02, 6.0, test_ct_reference, nearest)) == 120.0);
    // Outside the envelope it is still the nearest, from the clamped query.
    REQUIRE(surface_level(db, db.blend_at(0.9, -30.0, test_ct_reference, nearest)) == 110.0);

    // The default is unchanged: inside the triangle, a blend, not a pick.
    // (0.1, -5) is normalized (1/3, 1/4); sphere 0 is nearest.
    REQUIRE(surface_level(db, db.blend_at(0.1, -5.0, test_ct_reference, nearest)) == 100.0);
    REQUIRE(surface_level(db, db.blend_at(0.1, -5.0, test_ct_reference))
            == Catch::Approx(100.0 * 5.0 / 12.0 + 110.0 / 3.0 + 120.0 / 4.0));
    std::filesystem::remove(path);
}

TEST_CASE("--blend aam brackets speed in each path-angle row, then path angle", "[sphere_database][blend]")
{
    // Rows at -6 and -3 degrees with 60 and 75 kt in each; levels that are not
    // linear in (speed, angle), so AAM's bilinear blend and the Delaunay
    // blend differ.
    const std::vector<std::array<double, 4>> conditions = {
        {0.60, -6.0, test_ct_reference, 100.0},
        {0.75, -6.0, test_ct_reference, 106.0},
        {0.60, -3.0, test_ct_reference, 110.0},
        {0.75, -3.0, test_ct_reference, 120.0},
        {0.60, 4.0, test_ct_reference, 130.0},
    };
    const auto path = write_fixed_load_factor_test_database("niceops_blend_aam.nod", conditions, true);
    const Database db(path);
    REQUIRE(db.aam_layer_count() == 1);
    Database::BlendOptions aam;
    aam.policy = BlendPolicy::aam;
    aam.tip_speed = tip_speed_100_kt;

    // 65 kt at -4.5 degrees: a third of the way across each row's speeds, half
    // way across the rows (3 degrees apart, inside the approach gate).
    const double expected = 0.5 * (100.0 + 6.0 / 3.0) + 0.5 * (110.0 + 10.0 / 3.0);
    const auto blend = db.blend_at(0.65, -4.5, test_ct_reference, aam);
    REQUIRE(surface_level(db, blend) == Catch::Approx(expected).epsilon(1e-14));
    REQUIRE(surface_level(db, db.blend_at(0.65, -4.5, test_ct_reference)) != Catch::Approx(expected));

    // Climbing at 0.5 degrees: the -3/+4 rows are 7 degrees apart, past even
    // the departure gate, so the nearer row (-3) alone.
    REQUIRE(surface_level(db, db.blend_at(0.60, 0.5, test_ct_reference, aam)) == 110.0);
    // 110 kt: past every speed, the fastest in the row.
    REQUIRE(surface_level(db, db.blend_at(1.10, -6.0, test_ct_reference, aam)) == 106.0);

    // The 20 kt gate is converted with the tip speed: at 100 kt, mu 0.60 and
    // 0.75 are 15 kt apart and interpolated; at 150 kt they are 22.5 kt
    // apart, past the gate, and a state between them takes the nearer alone.
    aam.tip_speed = tip_speed_100_kt * 1.5;
    REQUIRE(surface_level(db, db.blend_at(0.65, -6.0, test_ct_reference, aam)) == 100.0);

    // A database with no tip speed of its own needs the vehicle's.
    aam.tip_speed = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS(db.blend_at(0.65, -4.5, test_ct_reference, aam));
    std::filesystem::remove(path);
}

TEST_CASE("--blend aam searches the nearest thrust layer of a load-factor database",
        "[sphere_database][blend]")
{
    // Not fixed_load_factor: two thrust coefficients, each a layer with its
    // own spheres.  AAM has no thrust axis, so a state takes the nearer
    // layer and searches it alone.
    const double low = 0.004, high = 0.008;
    const std::vector<std::array<double, 4>> conditions = {
        {0.5, -6.0, low, 100.0},
        {0.7, -6.0, low, 102.0},
        {0.5, 0.0, low, 104.0},
        {0.7, 0.0, low, 106.0},
        {0.5, -6.0, high, 110.0},
        {0.7, -6.0, high, 112.0},
        {0.5, 0.0, high, 114.0},
        {0.7, 0.0, high, 116.0},
    };
    const auto path = write_fixed_load_factor_test_database("niceops_blend_aam_layers.nod", conditions, false);
    const Database db(path);
    REQUIRE_FALSE(db.fixed_load_factor());
    REQUIRE(db.aam_layer_count() == 2);
    REQUIRE(db.aam_single_sphere_layers() == 0);
    Database::BlendOptions aam;
    aam.policy = BlendPolicy::aam;
    aam.tip_speed = tip_speed_100_kt;
    REQUIRE(surface_level(db, db.blend_at(0.5, -6.0, 0.0055, aam)) == 100.0);
    REQUIRE(surface_level(db, db.blend_at(0.5, -6.0, 0.0065, aam)) == 110.0);
    // Exactly between: the higher thrust, as speed ties take the faster.
    REQUIRE(surface_level(db, db.blend_at(0.5, -6.0, 0.006, aam)) == 110.0);
    // Within the layer it is the ordinary search: 60 kt, a row pair 6 degrees
    // apart (past the 3 degree approach gate), so the nearer row, two speeds.
    REQUIRE(surface_level(db, db.blend_at(0.6, -2.0, 0.0079, aam)) == Catch::Approx(115.0));
    std::filesystem::remove(path);
}

TEST_CASE("--blend aam counts the thrust layers that hold one sphere", "[sphere_database][blend]")
{
    // C_T scattered sphere by sphere rather than layered: every layer is one
    // sphere, and the run warns that AAM's search has nothing to bracket.
    const std::vector<std::array<double, 4>> conditions = {
        {0.5, -6.0, 0.0040, 100.0},
        {0.7, -6.0, 0.0050, 102.0},
        {0.5, 0.0, 0.0060, 104.0},
        {0.7, 0.0, 0.0075, 106.0},
    };
    const auto path = write_fixed_load_factor_test_database("niceops_blend_aam_scattered.nod", conditions, false);
    const Database db(path);
    REQUIRE(db.aam_layer_count() == 4);
    REQUIRE(db.aam_single_sphere_layers() == 4);
    std::filesystem::remove(path);
}

TEST_CASE("two-sided load-factor correction follows the law below the reference, to a floor",
        "[sphere_database][blend]")
{
    const auto path = write_fixed_load_factor_test_database(
            "niceops_fixed_lf_two_sided.nod", fixed_load_factor_conditions, true);
    const Database db(path);
    Database::BlendOptions two_sided;
    two_sided.load_factor_correction = LoadFactorCorrection::two_sided;
    two_sided.min_load_factor = 0.5;

    for (double load_factor : {0.5, 0.8, 1.0, 1.3, 3.0})
    {
        const auto blend = db.blend_at(0.0, 0.0, test_ct_reference * load_factor, two_sided);
        REQUIRE(blend.level_correction_db == Catch::Approx(20.0 * std::log10(load_factor)).margin(1e-12));
        REQUIRE(surface_level(db, blend) == Catch::Approx(100.0 + 20.0 * std::log10(load_factor)));
        // One-sided, the same states: nothing below the reference, the same law above it.
        const auto one_sided = db.blend_at(0.0, 0.0, test_ct_reference * load_factor);
        REQUIRE(one_sided.level_correction_db
                == (load_factor > 1.0 ? blend.level_correction_db : 0.0));
    }
    // Below the floor the correction is held there, through zero and below
    // (a state at negative g has a negative C_T, and no logarithm).
    for (double load_factor : {0.49, 0.2, 0.0, -0.5})
    {
        const auto blend = db.blend_at(0.0, 0.0, test_ct_reference * load_factor, two_sided);
        REQUIRE(blend.level_correction_db == Catch::Approx(20.0 * std::log10(0.5)));
    }
    two_sided.min_load_factor = 0.25;
    REQUIRE(db.blend_at(0.0, 0.0, test_ct_reference * 0.3, two_sided).level_correction_db
            == Catch::Approx(20.0 * std::log10(0.3)));
    REQUIRE(db.blend_at(0.0, 0.0, test_ct_reference * 0.1, two_sided).level_correction_db
            == Catch::Approx(20.0 * std::log10(0.25)));

    // An exported sphere carries it too.
    const auto sphere = db.interpolate_sphere(0.0, 0.0, test_ct_reference * 0.5,
            Database::BlendOptions {BlendPolicy::barycentric, LoadFactorCorrection::two_sided, 0.5});
    for (const auto& cell : sphere.cells)
        if (cell.dBA > -1000.0)
            REQUIRE(cell.dBA == Catch::Approx(100.0 + 20.0 * std::log10(0.5)));
    std::filesystem::remove(path);
}

TEST_CASE("direction bounds hold for two-sided corrections and the other policies",
        "[sphere_database][blend][bounds]")
{
    const std::vector<std::array<double, 4>> conditions = {
        {0.60, -6.0, test_ct_reference, 100.0},
        {0.75, -6.0, test_ct_reference, 106.0},
        {0.60, -3.0, test_ct_reference, 110.0},
        {0.75, -3.0, test_ct_reference, 120.0},
        {0.60, 4.0, test_ct_reference, 130.0},
    };
    const auto path = write_fixed_load_factor_test_database("niceops_blend_bounds.nod", conditions, true);
    const Database db(path);
    for (const auto policy : {BlendPolicy::barycentric, BlendPolicy::nearest, BlendPolicy::aam})
    {
        Database::BlendOptions options;
        options.policy = policy;
        options.tip_speed = tip_speed_100_kt;
        options.load_factor_correction = LoadFactorCorrection::two_sided;
        std::vector<Database::SourceState> states;
        // Two states at one condition and different load factors, and a few
        // more around the envelope.
        for (const auto& c : std::vector<std::array<double, 3>> {
                     {0.65, -4.5, 0.6}, {0.65, -4.5, 1.0}, {0.62, 0.5, 0.3}, {0.9, -9.0, 1.7}})
        {
            states.push_back(db.source_state(Vec3(0.0, 0.0, 100.0), Vec3(50.0, 0.0, -5.0),
                    Vec3(0.05, 0.0, 1.0), c[0], c[1], test_ct_reference * c[2], options));
        }
        const auto bounds = db.direction_bounds(states);
        REQUIRE(bounds.usable());
        const double radius = db.sphere_radius();
        for (std::size_t s = 0; s < states.size(); s++)
        {
            for (double x = -300.0; x <= 300.0; x += 37.0)
            {
                for (double y = -300.0; y <= 300.0; y += 41.0)
                {
                    const Vec3 observer(x, y, 0.0);
                    const Vec3 delta = observer - states[s].position;
                    const double squared_range = delta.squaredNorm();
                    if (squared_range <= radius * radius)
                        continue;
                    const int bin = DirectionBounds::bin(states[s].frame.nose.dot(delta),
                            states[s].frame.normal.dot(delta),
                            states[s].frame.starboard.dot(delta), squared_range);
                    const double level = db.level_at_observer(observer, states[s]);
                    INFO("policy " << to_string(policy) << ", state " << s);
                    REQUIRE(bounds.numerator(s, bin)
                            >= std::pow(10.0, 0.1 * level) * squared_range * (1.0 - 1e-9));
                }
            }
        }
        // cull_range stays a bound under a negative correction.
        const double cutoff = 95.0;
        const double range = db.cull_range(states[0].blend, cutoff);
        REQUIRE(std::isfinite(range));
        REQUIRE(db.eval(range, 0.0, 90.0, states[0].blend) <= cutoff + 1e-9);
    }
    std::filesystem::remove(path);
}
