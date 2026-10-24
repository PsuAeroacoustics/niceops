// Refraction by a linear effective sound-speed profile (refraction.hpp).
//
// The closed forms are checked against rays traced numerically -- a Runge-Kutta
// integration of the ray equations for c(z) = c0 + a z, which knows nothing of
// circles -- for position, travel time, arc length, and (by finite differences
// of neighboring rays) the ray tube behind the spreading.  Then the limits:
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

TEST_CASE("the spreading is the traced ray tube's", "[extended][refraction]")
{
    // r_eq^2 = r |dr/dtheta| sin|arrival| / cos(launch), with dr/dtheta from
    // two neighboring traced rays and r where each crosses z_r.
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
    const auto ray = refracted_ray(air, source, air.moments(300.0), observer, air.moments(0.0),
            GroundGeometry::point, 0.0);
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
    const auto ray = refracted_ray(air, source, air.moments(200.0), observer + Vec3(0.0, 0.0, h_r),
            air.moments(h_r), GroundGeometry::point, h_r);
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
    const ShadowAttenuation rigid(bands, 0, bands.size());
    shadow_loss(ray, rigid, 20.0, 0, 3, loss.data());
    // At 1 kHz and 0.02 /s the first mode decays at about 24 dB/km.
    CHECK(loss[1] == Catch::Approx(2.4111).epsilon(1e-4));
    ray.shadow_depth = 5000.0;
    shadow_loss(ray, rigid, 20.0, 0, 3, loss.data());
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
    const auto a = read_atmosphere_profile(polar.path, meters, 0.0);
    REQUIRE(a.has_wind());
    CHECK(a.wind_east[0] == Catch::Approx(10.0));
    CHECK(std::fabs(a.wind_north[0]) < 1e-12);
    CHECK(std::fabs(a.wind_east[1]) < 1e-12);
    CHECK(a.wind_north[1] == Catch::Approx(-10.0));

    TempCsv knots("z, T_K, RH, u_kt, v_kt\n0, 288, 50, 10, -5\n");
    const auto b = read_atmosphere_profile(knots.path, meters, 0.0);
    CHECK(b.wind_east[0] == Catch::Approx(10.0 * knot));
    CHECK(b.wind_north[0] == Catch::Approx(-5.0 * knot));

    TempCsv bare("z, T_K, RH, wind_speed, wind_direction\n0, 288, 50, 10, 270\n");
    CHECK_THROWS_AS(read_atmosphere_profile(bare.path, meters, 0.0), std::runtime_error);
    TempCsv half("z, T_K, RH, wind_speed_mps\n0, 288, 50, 10\n");
    CHECK_THROWS_AS(read_atmosphere_profile(half.path, meters, 0.0), std::runtime_error);
}

TEST_CASE("reflected_eigenray reaches a receiver above the source, or reports none", "[refraction]")
{
    // With sound bending down and the receiver above the source, launches steeper
    // than a window about the horizontal never climb to the receiver's height.  The
    // solve used to accept a zero of the clamped continuation there: a = 0.01,
    // z_s = 2, z_r = 3, r = 800 m returned a launch whose traced ray ended 1.1 m
    // short, with a grazing angle 19% low.  Every ray it reports must land.
    int found = 0, reported_none = 0;
    for (const double a : {0.003, 0.01, 0.03, 0.1})
    {
        for (const double z_s : {2.0, 5.0, 20.0})
        {
            for (const double z_r : {3.0, 8.0, 40.0})
            {
                if (!(z_r > z_s))
                    continue;
                for (double r = 60.0; r <= 3000.0; r += 140.0)
                {
                    const auto ray = reflected_eigenray({c0, a}, z_s, z_r, r);
                    INFO("a = " << a << ", z_s = " << z_s << ", z_r = " << z_r << ", r = " << r);
                    if (!ray.exists)
                    {
                        reported_none++;
                        continue;
                    }
                    const auto traced = trace(a, z_s, ray.launch, r, true);
                    REQUIRE(traced.grazing > 0.0);
                    CHECK(traced.z == Catch::Approx(z_r).margin(1e-5));
                    CHECK(traced.grazing == Catch::Approx(ray.grazing).margin(1e-7));
                    CHECK(traced.length == Catch::Approx(ray.path_length).epsilon(1e-8));
                    found++;
                }
            }
        }
    }
    CHECK(found > 200);
    // The case that returned a non-ray has none to return: at 800 m the only launches
    // that put the traced ray at 3 m are the ones that rise past it, turn over (its
    // apex is 3.3 m) and come back down -- a second crossing, arriving from above,
    // which this ray (the first, ascending one) is not.  At 500 m there is a first.
    CHECK_FALSE(reflected_eigenray({c0, 0.01}, 2.0, 3.0, 800.0).exists);
    const auto near = reflected_eigenray({c0, 0.01}, 2.0, 3.0, 500.0);
    REQUIRE(near.exists);
    CHECK(trace(0.01, 2.0, near.launch, 500.0, true).z == Catch::Approx(3.0).margin(1e-5));
    INFO(reported_none << " pairs had no reflected ray");
}

TEST_CASE("the fit projects the wind onto the path's own direction, both components", "[refraction]")
{
    // Every other wind test has one component zero or a profile that is linear in
    // both, so a flipped sign on the northward projection passed them all.  Still,
    // isothermal air with an east wind growing 0.025 /s and a north wind growing
    // 0.015 /s; the frame's +x is east and +y north (bearing 90).
    AtmosphereProfile profile;
    profile.height = {0.0, 400.0};
    profile.temperature = {288.15, 288.15};
    profile.relative_humidity = {50.0, 50.0};
    profile.pressure = {101.325, 97.0};
    profile.wind_east = {0.0, 10.0};
    profile.wind_north = {0.0, 6.0};
    const EffectiveSoundSpeed air(profile, 90.0, 0.0, 400.0);
    const auto top = air.moments(300.0);
    const auto bottom = air.moments(1.2);
    struct Direction
    {
        double hx, hy, gradient;
    };
    for (const Direction d : {Direction {0.6, 0.8, 0.6 * 0.025 + 0.8 * 0.015},
                 Direction {0.6, -0.8, 0.6 * 0.025 - 0.8 * 0.015},
                 Direction {-0.8, 0.6, -0.8 * 0.025 + 0.6 * 0.015},
                 Direction {0.0, 1.0, 0.015}, Direction {-0.6, -0.8, -0.6 * 0.025 - 0.8 * 0.015}})
    {
        INFO("direction (" << d.hx << ", " << d.hy << ")");
        CHECK(air.fit(300.0, top, 1.2, bottom, d.hx, d.hy).gradient == Catch::Approx(d.gradient).epsilon(1e-9));
    }
    // And in a frame turned 210 degrees the same wind is projected through that turn:
    // a path along the compass bearing 120 (toward the south-east) is local (cos 30, -sin 30)
    // when +x is bearing 150... easiest to check against the wind itself.
    const EffectiveSoundSpeed turned(profile, 210.0, 0.0, 400.0);
    const auto w = turned.wind(400.0);
    const double along = 0.8 * w[0] + 0.6 * w[1];
    CHECK(turned.fit(300.0, turned.moments(300.0), 1.2, turned.moments(1.2), 0.8, 0.6).gradient
            == Catch::Approx(along / 400.0).epsilon(1e-9));
}

TEST_CASE("a short fit span is widened upward, past a kink, not downward", "[refraction]")
{
    // A level path spans too little height to fit a gradient through, so the span is
    // widened to min_fit_span: upward, to where the sound will actually be going.  On a
    // linear profile any window gives the same answer, so this one has no wind below
    // 10 m and a shear above it.
    AtmosphereProfile profile;
    profile.height = {0.0, 10.0, 400.0};
    profile.temperature = {288.15, 288.15, 288.15};
    profile.relative_humidity = {50.0, 50.0, 50.0};
    profile.pressure = {101.325, 101.2, 97.0};
    profile.wind_east = {0.0, 0.0, 7.8};
    profile.wind_north = {0.0, 0.0, 0.0};
    const EffectiveSoundSpeed air(profile, 90.0, 0.0, 400.0);
    // The least-squares slope of the wind over [5, 15], by quadrature.
    const double z1 = 5.0, z2 = 5.0 + EffectiveSoundSpeed::min_fit_span;
    auto wind = [](double z) { return z <= 10.0 ? 0.0 : 0.02 * (z - 10.0); };
    const int n = 200000;
    double mean = 0.0, moment = 0.0;
    const double z_mid = 0.5 * (z1 + z2);
    for (int i = 0; i < n; i++)
    {
        const double z = z1 + (z2 - z1) * (i + 0.5) / n;
        mean += wind(z) / n;
        moment += wind(z) * (z - z_mid) / n;
    }
    const double expected = moment / ((z2 - z1) * (z2 - z1) / 12.0);
    REQUIRE(expected > 0.005);
    const auto fit = air.fit(5.0, air.moments(5.0), 5.5, air.moments(5.5), 1.0, 0.0);
    CHECK(fit.gradient == Catch::Approx(expected).epsilon(1e-4));
    // Source and receiver in either order are the same path.
    const auto reversed = air.fit(5.5, air.moments(5.5), 5.0, air.moments(5.0), 1.0, 0.0);
    CHECK(reversed.gradient == Catch::Approx(fit.gradient).epsilon(1e-12));
}

TEST_CASE("a flush microphone sees the direct ray arrive at its own angle", "[refraction]")
{
    // finish_ground's plate branch hands the ground model the direct ray's arrival
    // geometry: a source height and ground distance whose slant has the arrival's
    // elevation and the ray's path length.  Swapping sin and cos in it passed every
    // test, none of which asked for GroundGeometry::plate on a curved ray.
    AtmosphereProfile profile;
    profile.height = {0.0, 600.0};
    profile.temperature = {288.15, 288.15};
    profile.relative_humidity = {50.0, 50.0};
    profile.pressure = {101.325, 94.0};
    profile.wind_east = {0.0, 12.0}; // 0.02 /s toward the east: downwind bends the sound down
    profile.wind_north = {0.0, 0.0};
    const EffectiveSoundSpeed air(profile, 90.0, 0.0, 700.0);
    for (const double r : {400.0, 900.0, 1500.0})
    {
        const Vec3 source(0.0, 0.0, 120.0);
        const Vec3 observer(r, 0.0, 0.0);
        const auto ray = refracted_ray(air, source, air.moments(120.0), observer, air.moments(0.0),
                GroundGeometry::plate, 0.0);
        REQUIRE(ray.curved);
        REQUIRE(ray.direct.lit);
        // The numerically traced ray's elevation where it arrives.
        const auto traced = trace(ray.profile.gradient, 120.0, ray.direct.launch, r, false, ray.profile.sound_speed);
        INFO("r = " << r);
        CHECK(traced.z == Catch::Approx(0.0).margin(1e-6));
        CHECK(std::atan2(ray.path.source_height, ray.path.ground_distance)
                == Catch::Approx(std::fabs(traced.arrival)).margin(1e-7));
        CHECK(std::hypot(ray.path.source_height, ray.path.ground_distance)
                == Catch::Approx(ray.direct.path_length).epsilon(1e-9));
        // The ray is steeper than the straight chord's 120/r: it was bent down.
        CHECK(ray.path.source_height / ray.path.ground_distance > 120.0 / r);
    }
}
