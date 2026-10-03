// Ported from the unittest blocks in sphere_database.d.  The D tests use
// std.math.approxEqual (relative epsilon 1e-6), so the port uses the C++
// approx_equal() from sphere_database.hpp to keep identical pass/fail
// semantics rather than Catch::Approx.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

#include <netcdf.h>

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
        double sphere_radius_m = 30.48)
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
        REQUIRE(nc_def_var(group, "radius", NC_DOUBLE, 1, &radii_dim, &v.radius) == NC_NOERR);
        REQUIRE(nc_def_var(group, "rotor_scale", NC_DOUBLE, 1, &condition_dim, &v.scale) == NC_NOERR);
        REQUIRE(nc_def_var(group, "advance_ratio", NC_DOUBLE, 1, &condition_dim, &v.mu) == NC_NOERR);
        REQUIRE(nc_def_var(group, "flight_path_angle", NC_DOUBLE, 1, &condition_dim, &v.gamma)
                == NC_NOERR);
        REQUIRE(nc_def_var(group, "thrust_coefficient", NC_DOUBLE, 1, &condition_dim, &v.ct)
                == NC_NOERR);
        REQUIRE(nc_def_var(group, "phi", NC_DOUBLE, 1, &channel_dim, &v.phi) == NC_NOERR);
        REQUIRE(nc_def_var(group, "theta", NC_DOUBLE, 1, &channel_dim, &v.theta) == NC_NOERR);
        REQUIRE(nc_def_var(group, "dBA", NC_DOUBLE, 1, &channel_dim, &v.dba) == NC_NOERR);
        REQUIRE(nc_def_var(group, "EAA", NC_DOUBLE, 1, &channel_dim, &v.eaa) == NC_NOERR);
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
        const std::vector<double> eaa(n_channels, 0.0);
        const std::vector<double> dba(n_channels, conditions[i][3]);
        REQUIRE(nc_put_var_double(v.group, v.radius, &sphere_radius_m) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.scale, &scale) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.mu, &conditions[i][0]) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.gamma, &conditions[i][1]) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.ct, &conditions[i][2]) == NC_NOERR);
        REQUIRE(nc_put_var_double(v.group, v.phi, phi_channel.data()) == NC_NOERR);
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
// see which direction every spectrum landed in.
std::string write_shared_grid_test_database(const std::string& name,
        const std::vector<std::array<double, 4>>& conditions, bool shared_grid_and_frequency,
        bool with_spectrum, double sphere_radius_m = 30.48,
        std::function<double(std::size_t, std::size_t)> channel_dba = {},
        bool reverse_phi_rows = false)
{
    const auto path = (std::filesystem::temp_directory_path() / name).string();
    int file = 0;
    REQUIRE(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file) == NC_NOERR);
    int scalar_dim = 0;
    REQUIRE(nc_def_dim(file, "scalar", 1, &scalar_dim) == NC_NOERR);
    int same_grid_var = 0, version_var = 0, shared_var = 0;
    REQUIRE(nc_def_var(file, "same_grid", NC_SHORT, 1, &scalar_dim, &same_grid_var) == NC_NOERR);
    REQUIRE(nc_def_var(file, "database_version", NC_INT, 1, &scalar_dim, &version_var) == NC_NOERR);
    REQUIRE(nc_def_var(file, "shared_grid_and_frequency", NC_SHORT, 1, &scalar_dim, &shared_var)
            == NC_NOERR);

    const std::vector<double> phis = {-90.0, -30.0, 30.0, 90.0};
    const std::vector<double> thetas = {0.0, 60.0, 120.0, 180.0};
    const std::size_t n_channels = phis.size() * thetas.size();
    const std::vector<double> bands = {500.0, 2000.0};
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
                    amplitude[cell * bands.size() + band] = reverse_phi_rows
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
    const std::vector<std::array<double, 4>> conditions = {
        {0.0, 0.0, 0.005, 100.0},
        {1.0, 0.0, 0.005, 110.0},
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
