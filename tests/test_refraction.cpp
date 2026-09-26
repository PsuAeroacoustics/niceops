// Refraction by a linear effective sound-speed profile (refraction.hpp).
//
// The closed forms are checked against rays traced numerically -- a Runge-Kutta
// integration of the ray equations for c(z) = c0 + a z, which knows nothing of
// circles -- for position, travel time, arc length, and (by finite differences
// of neighbouring rays) the ray tube behind the spreading.  Then the limits:
// straight rays as the gradient vanishes, the launch-angle shifts and shadow
// boundaries tabulated in the propagation roadmap, and continuity across the
// shadow boundary.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "refraction.hpp"
#include "util/constants.hpp"

using namespace niceops;

namespace {

constexpr double c0 = 340.0;

struct Traced
{
    double z = 0.0; // height where the ray reaches x = r
    double time = 0.0;
    double length = 0.0;
    double arrival = 0.0; // elevation there
    double grazing = -1.0; // angle at the ground bounce, if any
};

// Trace a ray from (0, z_s) at elevation `launch` until it reaches horizontal
// distance r, by RK4 in arc length on (x, z, theta, t):
//   dx/ds = cos theta, dz/ds = sin theta, dtheta/ds = -a cos theta / c,
//   dt/ds = 1/c.
// With `reflect`, the ray bounces off z = 0 (the first time it gets there).
Traced trace(double a, double z_s, double launch, double r, bool reflect, double sound_speed = c0)
{
    struct State
    {
        double x, z, theta, t;
    };
    auto derivative = [a, sound_speed](const State& s) {
        const double c = sound_speed + a * s.z;
        return State {std::cos(s.theta), std::sin(s.theta), -a * std::cos(s.theta) / c, 1.0 / c};
    };
    auto step = [&](const State& s, double h) {
        auto add = [](const State& p, const State& d, double k) {
            return State {p.x + k * d.x, p.z + k * d.z, p.theta + k * d.theta, p.t + k * d.t};
        };
        const State k1 = derivative(s);
        const State k2 = derivative(add(s, k1, 0.5 * h));
        const State k3 = derivative(add(s, k2, 0.5 * h));
        const State k4 = derivative(add(s, k3, h));
        return State {s.x + h / 6 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x),
                s.z + h / 6 * (k1.z + 2 * k2.z + 2 * k3.z + k4.z),
                s.theta + h / 6 * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta),
                s.t + h / 6 * (k1.t + 2 * k2.t + 2 * k3.t + k4.t)};
    };
    Traced out;
    State s {0.0, z_s, launch, 0.0};
    double length = 0.0;
    bool bounced = false;
    const double h0 = 0.05;
    while (s.x < r)
    {
        double h = h0;
        State next = step(s, h);
        if (reflect && !bounced && next.z < 0.0)
        {
            // Bisect the step onto the ground, then mirror the angle.
            double lo = 0.0, hi = h;
            for (int i = 0; i < 60; i++)
            {
                const double mid = 0.5 * (lo + hi);
                (step(s, mid).z < 0.0 ? hi : lo) = mid;
            }
            s = step(s, lo);
            length += lo;
            out.grazing = -s.theta;
            s.theta = -s.theta;
            s.z = 0.0;
            bounced = true;
            continue;
        }
        if (next.x >= r)
        {
            double lo = 0.0, hi = h;
            for (int i = 0; i < 60; i++)
            {
                const double mid = 0.5 * (lo + hi);
                (step(s, mid).x < r ? lo : hi) = mid;
            }
            const State end = step(s, hi);
            out.z = end.z;
            out.time = end.t;
            out.length = length + hi;
            out.arrival = end.theta;
            return out;
        }
        s = next;
        length += h;
    }
    return out;
}

double chord_angle(double z_s, double z_r, double r)
{
    return std::atan2(z_r - z_s, r);
}

} // namespace

TEST_CASE("direct_eigenray matches a numerically traced ray", "[refraction]")
{
    for (const double a : {0.05, 0.01, -0.01, -0.05})
    {
        for (const double r : {300.0, 1000.0, 1500.0})
        {
            const double z_s = 150.0;
            const double z_r = 1.2;
            const LinearProfile profile {c0, a};
            const auto ray = direct_eigenray(profile, z_s, z_r, r);
            if (!ray.lit)
                continue;
            const auto traced = trace(a, z_s, ray.launch, r, false);
            INFO("a = " << a << ", r = " << r);
            CHECK(traced.z == Catch::Approx(z_r).margin(1e-6));
            CHECK(traced.time == Catch::Approx(ray.travel_time).epsilon(1e-9));
            CHECK(traced.length == Catch::Approx(ray.path_length).epsilon(1e-9));
            CHECK(traced.arrival == Catch::Approx(ray.arrival).margin(1e-8));
        }
    }
}

TEST_CASE("the spreading is the traced ray tube's", "[refraction]")
{
    // r_eq^2 = r |dr/dtheta| sin|arrival| / cos(launch), with dr/dtheta from
    // two neighbouring traced rays and r where each crosses z_r.
    const double z_s = 150.0;
    const double z_r = 1.2;
    for (const double a : {0.05, -0.02})
    {
        const double r = 1200.0;
        const LinearProfile profile {c0, a};
        const auto ray = direct_eigenray(profile, z_s, z_r, r);
        REQUIRE(ray.lit);
        // Where a ray launched at theta crosses z_r: solve on r by bisection,
        // tracing to each trial distance.
        auto landing = [&](double theta) {
            double lo = 0.5 * r, hi = 2.0 * r;
            for (int i = 0; i < 60; i++)
            {
                const double mid = 0.5 * (lo + hi);
                (trace(a, z_s, theta, mid, false).z > z_r ? lo : hi) = mid;
            }
            return 0.5 * (lo + hi);
        };
        const double delta = 1e-5;
        const double drdtheta = (landing(ray.launch + delta) - landing(ray.launch - delta)) / (2 * delta);
        const double tube = r * std::fabs(drdtheta) * std::sin(std::fabs(ray.arrival)) / std::cos(ray.launch);
        INFO("a = " << a);
        CHECK(std::sqrt(tube) == Catch::Approx(ray.spreading_range).epsilon(1e-5));
    }
}

TEST_CASE("direct_eigenray is the straight ray as the gradient vanishes", "[refraction]")
{
    const double z_s = 300.0, z_r = 1.2, r = 2500.0;
    const double slant = std::hypot(r, z_s - z_r);
    for (const double a : {1e-8, -1e-8})
    {
        const auto ray = direct_eigenray({c0, a}, z_s, z_r, r);
        CHECK(ray.lit);
        CHECK(ray.launch == Catch::Approx(chord_angle(z_s, z_r, r)).margin(1e-7));
        CHECK(ray.spreading_range == Catch::Approx(slant).epsilon(1e-7));
        CHECK(ray.path_length == Catch::Approx(slant).epsilon(1e-9));
        CHECK(ray.travel_time == Catch::Approx(slant / c0).epsilon(1e-5));
    }
}

TEST_CASE("launch-angle shifts reproduce the roadmap's table", "[refraction]")
{
    // Source at 150 m, receiver at 1.2 m, downward refraction: the shift from
    // the chord is about d |dc/dz| / (2 c0).  The table was computed with
    // c0 = 343 m/s.
    struct Row
    {
        double gradient, distance, shift_degrees;
    };
    for (const Row row : {Row {0.01, 2000.0, 1.6663}, Row {0.05, 2000.0, 8.2046},
                 Row {0.02, 4000.0, 6.6227}, Row {0.005, 500.0, 0.2086}})
    {
        const auto ray = direct_eigenray({343.0, row.gradient}, 150.0, 1.2, row.distance);
        const double shift = (ray.launch - chord_angle(150.0, 1.2, row.distance)) * 180.0 / pi;
        CHECK(shift == Catch::Approx(row.shift_degrees).margin(1e-4));
    }
}

TEST_CASE("the shadow boundary is where the grazing ray reaches the receiver", "[refraction]")
{
    const double z_s = 150.0, z_r = 1.2, a = -0.02;
    const LinearProfile profile {343.0, a};
    const double boundary = shadow_boundary(profile, z_s, z_r);
    // Exactly sqrt(z (2R - z)) out from the grazing point on each side, R = c0/|a|;
    // the roadmap's estimate sqrt(2R)(sqrt(h_s) + sqrt(h_r)) drops the z^2,
    // which is z/4R, 0.2% here.
    const double R = 343.0 / std::fabs(a);
    CHECK(boundary
            == Catch::Approx(std::sqrt(z_s * (2 * R - z_s)) + std::sqrt(z_r * (2 * R - z_r))).epsilon(1e-12));
    CHECK(boundary == Catch::Approx(std::sqrt(2 * R) * (std::sqrt(z_s) + std::sqrt(z_r))).epsilon(3e-3));
    CHECK(shadow_boundary({343.0, 0.02}, z_s, z_r) == std::numeric_limits<double>::infinity());

    const auto inside = direct_eigenray(profile, z_s, z_r, boundary * (1 - 1e-9));
    const auto outside = direct_eigenray(profile, z_s, z_r, boundary * (1 + 1e-9));
    CHECK(inside.lit);
    CHECK_FALSE(outside.lit);
    CHECK(outside.shadow_depth == Catch::Approx(boundary * 1e-9).epsilon(1e-3));
    // Continuous across it: nothing the level reads jumps.
    CHECK(outside.launch == Catch::Approx(inside.launch).margin(1e-6));
    CHECK(outside.spreading_range == Catch::Approx(inside.spreading_range).epsilon(1e-6));
    CHECK(outside.path_length == Catch::Approx(inside.path_length).epsilon(1e-6));
    CHECK(outside.travel_time == Catch::Approx(inside.travel_time).epsilon(1e-6));
    // The limiting ray touches the ground: traced, it just reaches z = 0.
    const auto grazing = direct_eigenray(profile, z_s, z_r, boundary);
    const double t_s = std::sqrt(z_s * (2 * R - z_s));
    CHECK(std::fabs(trace(a, z_s, grazing.launch, t_s, false, 343.0).z) < 1e-4);
    CHECK_FALSE(reflected_eigenray(profile, z_s, z_r, boundary * 1.01).exists);
}

TEST_CASE("reflected_eigenray matches a numerically traced bounce", "[refraction]")
{
    for (const double a : {0.05, 0.01, -0.01, -0.03})
    {
        for (const double r : {200.0, 800.0, 1500.0})
        {
            const double z_s = 150.0, z_r = 1.2;
            const LinearProfile profile {c0, a};
            const auto ray = reflected_eigenray(profile, z_s, z_r, r);
            if (!ray.exists)
                continue;
            const auto traced = trace(a, z_s, ray.launch, r, true);
            INFO("a = " << a << ", r = " << r);
            REQUIRE(traced.grazing > 0.0);
            CHECK(traced.z == Catch::Approx(z_r).margin(1e-6));
            CHECK(traced.grazing == Catch::Approx(ray.grazing).margin(1e-8));
            CHECK(traced.length == Catch::Approx(ray.path_length).epsilon(1e-9));
        }
    }
}

TEST_CASE("reflected_eigenray is the image ray as the gradient vanishes", "[refraction]")
{
    const double z_s = 150.0, z_r = 1.2, r = 900.0;
    const auto ray = reflected_eigenray({c0, 1e-8}, z_s, z_r, r);
    REQUIRE(ray.exists);
    CHECK(ray.grazing == Catch::Approx(std::atan2(z_s + z_r, r)).margin(1e-7));
    CHECK(ray.path_length == Catch::Approx(std::hypot(r, z_s + z_r)).epsilon(1e-7));
    // Downward refraction at long range launches the reflected ray upward,
    // over its top: the other branch of the range.
    const auto far = reflected_eigenray({c0, 0.05}, 20.0, 1.2, 3000.0);
    REQUIRE(far.exists);
    CHECK(far.launch > 0.0);
    const auto traced = trace(0.05, 20.0, far.launch, 3000.0, true);
    CHECK(traced.z == Catch::Approx(1.2).margin(1e-6));
}

TEST_CASE("EffectiveSoundSpeed fits a linear profile exactly", "[refraction]")
{
    // Still, isothermal air with a wind toward the east growing linearly: the
    // effective sound speed is exactly linear along any bearing.
    AtmosphereProfile profile;
    profile.height = {0.0, 400.0};
    profile.temperature = {288.15, 288.15};
    profile.relative_humidity = {50.0, 50.0};
    profile.pressure = {101.325, 97.0};
    profile.wind_east = {2.0, 10.0}; // shear 0.02 /s
    profile.wind_north = {0.0, 0.0};
    const double c = profile.at(0.0).sound_speed();

    SECTION("x east")
    {
        const EffectiveSoundSpeed air(profile, 90.0, 0.0, 400.0);
        const auto top = air.moments(300.0);
        const auto bottom = air.moments(1.2);
        const auto downwind = air.fit(300.0, top, 1.2, bottom, 1.0, 0.0);
        CHECK(downwind.gradient == Catch::Approx(0.02).epsilon(1e-10));
        CHECK(downwind.sound_speed == Catch::Approx(c + 2.0).epsilon(1e-12));
        const auto upwind = air.fit(300.0, top, 1.2, bottom, -1.0, 0.0);
        CHECK(upwind.gradient == Catch::Approx(-0.02).epsilon(1e-10));
        const auto crosswind = air.fit(300.0, top, 1.2, bottom, 0.0, 1.0);
        CHECK(std::fabs(crosswind.gradient) < 1e-12);
    }

    SECTION("x along 210 degrees: an east wind is (-1/2, sqrt(3)/2) in the local frame")
    {
        const EffectiveSoundSpeed air(profile, 210.0, 0.0, 400.0);
        const auto w = air.wind(100.0);
        CHECK(w[0] == Catch::Approx(-0.5 * 4.0).epsilon(1e-12));
        CHECK(w[1] == Catch::Approx(std::sqrt(3.0) / 2.0 * 4.0).epsilon(1e-12));
    }

    SECTION("a short span is widened to min_fit_span")
    {
        const EffectiveSoundSpeed air(profile, 90.0, 0.0, 400.0);
        const auto fit = air.fit(1.0, air.moments(1.0), 1.5, air.moments(1.5), 1.0, 0.0);
        CHECK(fit.gradient == Catch::Approx(0.02).epsilon(1e-9));
    }
}

TEST_CASE("refracted_ray in uniform air is the straight ray", "[refraction]")
{
    AtmosphereProfile profile;
    profile.height = {0.0};
    profile.temperature = {290.0};
    profile.relative_humidity = {40.0};
    profile.pressure = {100.0};
    const EffectiveSoundSpeed air(profile, 90.0, 0.0, 600.0);
    const Vec3 source(-400.0, 250.0, 300.0);
    const Vec3 observer(700.0, -150.0, 0.0);
    const auto ray = refracted_ray(air, source, air.moments(300.0), observer, 0.0, air.moments(0.0),
            GroundGeometry::point);
    const auto straight = straight_ray(source, observer);
    CHECK(ray.path.range == straight.range);
    CHECK(ray.path.path_length == straight.path_length);
    CHECK(ray.path.direction.x == straight.direction.x);
    CHECK(ray.path.direction.z == straight.direction.z);
    CHECK(ray.path.source_height == straight.source_height);
    CHECK(ray.path.ground_distance == straight.ground_distance);
    CHECK(ray.shadow_depth == 0.0);
    CHECK(ray.travel_time == Catch::Approx(straight.range / profile.at(0.0).sound_speed()).epsilon(1e-12));
}

TEST_CASE("refracted_ray hands the ground model the reflected ray's geometry", "[refraction]")
{
    // A tiny gradient: the equivalent geometry is the real one.
    AtmosphereProfile profile;
    profile.height = {0.0, 600.0};
    profile.temperature = {288.15, 288.15};
    profile.relative_humidity = {50.0, 50.0};
    profile.pressure = {101.325, 94.0};
    profile.wind_east = {0.0, 1e-4};
    profile.wind_north = {0.0, 0.0};
    const EffectiveSoundSpeed air(profile, 90.0, 0.0, 700.0);
    const Vec3 source(0.0, 0.0, 200.0);
    const Vec3 observer(900.0, 0.0, 0.0);
    const double h_r = 1.2;
    const auto ray = refracted_ray(air, source, air.moments(200.0), observer, h_r, air.moments(h_r),
            GroundGeometry::point);
    CHECK(ray.path.source_height == Catch::Approx(200.0).epsilon(1e-6));
    CHECK(ray.path.ground_distance == Catch::Approx(900.0).epsilon(1e-6));
}

TEST_CASE("shadow_loss grows with depth and frequency and stops at the floor", "[refraction]")
{
    const std::vector<double> bands = {100.0, 1000.0, 8000.0};
    const auto rates = shadow_attenuation_rates(bands, 0, bands.size());
    // (pi f)^(1/3): a decade in frequency is 10^(1/3) in rate.
    CHECK(rates[1] / rates[0] == Catch::Approx(std::cbrt(10.0)).epsilon(1e-12));
    RefractedRay ray;
    ray.profile = {343.0, -0.02};
    ray.shadow_depth = 100.0;
    std::vector<double> loss(3, 0.0);
    shadow_loss(ray, rates, 20.0, 0, 3, loss.data());
    // At 1 kHz and 0.02 /s the first mode decays at about 24 dB/km.
    CHECK(loss[1] == Catch::Approx(2.4111).epsilon(1e-4));
    ray.shadow_depth = 5000.0;
    shadow_loss(ray, rates, 20.0, 0, 3, loss.data());
    CHECK(loss[2] == 20.0);
}

TEST_CASE("read_atmosphere_profile reads wind", "[refraction]")
{
    struct TempCsv
    {
        std::string path;
        explicit TempCsv(const std::string& contents)
        {
            char pattern[] = "/tmp/niceops_wind_XXXXXX";
            const int fd = ::mkstemp(pattern);
            REQUIRE(fd != -1);
            ::close(fd);
            path = pattern;
            std::ofstream(path) << contents;
        }
        ~TempCsv()
        {
            std::remove(path.c_str());
        }
    };

    // From the west at 10 m/s blows toward the east; from the north toward the
    // south.
    TempCsv polar("z, T_K, RH, wind_speed_mps, wind_direction\n0, 288, 50, 10, 270\n100, 287, 50, 10, 0\n");
    const auto a = read_atmosphere_profile(polar.path, metres, 0.0);
    REQUIRE(a.has_wind());
    CHECK(a.wind_east[0] == Catch::Approx(10.0));
    CHECK(std::fabs(a.wind_north[0]) < 1e-12);
    CHECK(std::fabs(a.wind_east[1]) < 1e-12);
    CHECK(a.wind_north[1] == Catch::Approx(-10.0));

    TempCsv knots("z, T_K, RH, u_kt, v_kt\n0, 288, 50, 10, -5\n");
    const auto b = read_atmosphere_profile(knots.path, metres, 0.0);
    CHECK(b.wind_east[0] == Catch::Approx(10.0 * 1852.0 / 3600.0));
    CHECK(b.wind_north[0] == Catch::Approx(-5.0 * 1852.0 / 3600.0));

    TempCsv bare("z, T_K, RH, wind_speed, wind_direction\n0, 288, 50, 10, 270\n");
    CHECK_THROWS_AS(read_atmosphere_profile(bare.path, metres, 0.0), std::runtime_error);
    TempCsv half("z, T_K, RH, wind_speed_mps\n0, 288, 50, 10\n");
    CHECK_THROWS_AS(read_atmosphere_profile(half.path, metres, 0.0), std::runtime_error);
}
