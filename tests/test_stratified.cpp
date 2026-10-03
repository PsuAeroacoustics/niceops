// Refraction through a layered profile (stratified.hpp).
//
// Two independent checks.  On a two-level profile -- a linear one -- the layered
// model must reproduce the linear model's closed forms.  On genuinely layered
// profiles its rays must match Runge-Kutta traces of the ray equations through
// the same piecewise-linear sound speed, including rays that climb, turn aloft
// and come back down into a near-ground shadow.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include "atmosphere_profile.hpp"
#include "refraction.hpp"
#include "stratified.hpp"
#include "util/constants.hpp"

using namespace niceops;

namespace {

struct Traced
{
    double z = 0.0;
    double time = 0.0;
    double length = 0.0;
    double grazing = -1.0;
};

// RK4 in arc length through c(z) given by `layers`, from (0, z_s) at elevation
// `launch` until horizontal distance r; bouncing off z = 0 up to `reflections`
// times (grazing is the last bounce's).
Traced trace(const Layers& layers, double z_s, double launch, double r, int reflections, double h = 0.02)
{
    auto c = [&](double z) { return layers.at(z); };
    // The layer's own gradient: a finite difference here rounds to ~1e-8 /s,
    // which over a kilometre of ray is a tenth of a millimetre.
    // Exactly on a level -- where the trace lands on purpose -- the layer is
    // the one the ray is heading into.
    auto dcdz = [&](double z, bool down) {
        const std::size_t n = layers.levels;
        if (z < layers.z[0] || z > layers.z[n - 1] || (down ? z <= layers.z[0] : z >= layers.z[n - 1]))
            return 0.0;
        std::size_t i = 0;
        while (i + 2 < n && (down ? layers.z[i + 1] < z : layers.z[i + 1] <= z))
            i++;
        return (layers.c[i + 1] - layers.c[i]) / (layers.z[i + 1] - layers.z[i]);
    };
    struct S
    {
        double x, z, theta, t;
    };
    // One gradient for all four stages: the layer the step starts in.  A step
    // ending on a level would otherwise take its last stage's gradient from
    // the layer beyond, which is a first-order error per crossing.
    auto f = [&](const S& s, double g) {
        const double cc = c(s.z);
        return S {std::cos(s.theta), std::sin(s.theta), -g * std::cos(s.theta) / cc, 1.0 / cc};
    };
    auto step = [&](const S& s, double h) {
        auto add = [](const S& a, const S& d, double k) {
            return S {a.x + k * d.x, a.z + k * d.z, a.theta + k * d.theta, a.t + k * d.t};
        };
        const double g = dcdz(s.z, std::sin(s.theta) < 0.0);
        const S k1 = f(s, g), k2 = f(add(s, k1, h / 2), g), k3 = f(add(s, k2, h / 2), g),
                k4 = f(add(s, k3, h), g);
        return S {s.x + h / 6 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x),
                s.z + h / 6 * (k1.z + 2 * k2.z + 2 * k3.z + k4.z),
                s.theta + h / 6 * (k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta),
                s.t + h / 6 * (k1.t + 2 * k2.t + 2 * k3.t + k4.t)};
    };
    Traced out;
    S s {0.0, z_s, launch, 0.0};
    double length = 0.0;
    int bounced = 0;
    // A step that straddles a level would integrate across the gradient's
    // jump there, and RK4 drops to first order; land on the level instead.
    auto crosses_level = [&](double z0, double z1) {
        for (std::size_t k = 0; k < layers.levels; k++)
        {
            const double level = layers.z[k];
            if ((z0 - level) * (z1 - level) < 0.0)
                return level;
        }
        return std::numeric_limits<double>::quiet_NaN();
    };
    while (true)
    {
        S next = step(s, h);
        const double level = crosses_level(s.z, next.z);
        if (!std::isnan(level) && next.x < r && !(bounced < reflections && next.z < 0.0))
        {
            double lo = 0, hi = h;
            for (int i = 0; i < 60; i++)
            {
                const double mid = 0.5 * (lo + hi);
                ((step(s, mid).z - level) * (s.z - level) > 0.0 ? lo : hi) = mid;
            }
            s = step(s, hi);
            s.z = level;
            length += hi;
            continue;
        }
        if (bounced < reflections && next.z < 0.0)
        {
            double lo = 0, hi = h;
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
            bounced++;
            continue;
        }
        if (next.x >= r)
        {
            double lo = 0, hi = h;
            for (int i = 0; i < 60; i++)
            {
                const double mid = 0.5 * (lo + hi);
                (step(s, mid).x < r ? lo : hi) = mid;
            }
            const S end = step(s, hi);
            out.z = end.z;
            out.time = end.t;
            out.length = length + hi;
            return out;
        }
        s = next;
        length += h;
    }
}

} // namespace

TEST_CASE("layered eigenrays reproduce the linear model on a linear profile", "[stratified]")
{
    for (const double a : {0.05, 0.01, -0.01, -0.03})
    {
        const double c0 = 340.0;
        // Tall enough that nothing here turns above the top.
        const std::vector<double> z = {0.0, 5000.0};
        const std::vector<double> c = {c0, c0 + a * 5000.0};
        const Layers layers {z.data(), c.data(), 2};
        const LinearProfile linear {c0, a};
        for (const double r : {300.0, 1200.0, 2500.0})
        {
            const double z_s = 150.0, z_r = 1.2;
            INFO("a = " << a << ", r = " << r);
            const auto expected = direct_eigenray(linear, z_s, z_r, r);
            const auto got = layered_direct_eigenray(layers, z_s, z_r, r, expected.launch + 0.01);
            REQUIRE(got.lit == expected.lit);
            CHECK(got.launch == Catch::Approx(expected.launch).margin(1e-9));
            CHECK(got.path_length == Catch::Approx(expected.path_length).epsilon(1e-9));
            CHECK(got.travel_time == Catch::Approx(expected.travel_time).epsilon(1e-9));
            CHECK(got.shadow_depth == Catch::Approx(expected.shadow_depth).margin(1e-6));
            if (expected.lit)
            {
                CHECK(got.arrival == Catch::Approx(expected.arrival).margin(1e-9));
                CHECK(got.spreading_range == Catch::Approx(expected.spreading_range).epsilon(1e-5));
            }
            const auto reflected_expected = reflected_eigenray(linear, z_s, z_r, r);
            const auto reflected = layered_reflected_eigenray(layers, z_s, z_r, r, -0.1);
            REQUIRE(reflected.exists == reflected_expected.exists);
            if (reflected.exists)
            {
                CHECK(reflected.grazing == Catch::Approx(reflected_expected.grazing).margin(1e-9));
                CHECK(reflected.path_length
                        == Catch::Approx(reflected_expected.path_length).epsilon(1e-9));
            }
        }
    }
}

TEST_CASE("layered eigenrays match traced rays through a layered profile", "[stratified]")
{
    // Faster, slower, faster: a near-ground gradient, a lid, and wind aloft.
    const std::vector<double> z = {0.0, 10.0, 60.0, 250.0, 900.0};
    const std::vector<double> c = {340.0, 343.0, 344.5, 343.5, 350.0};
    const Layers layers {z.data(), c.data(), z.size()};
    const double z_s = 180.0, z_r = 1.2;
    for (const double r : {250.0, 900.0, 1800.0})
    {
        INFO("r = " << r);
        const auto ray = layered_direct_eigenray(layers, z_s, z_r, r, std::atan2(z_r - z_s, r));
        if (ray.lit)
        {
            const auto traced = trace(layers, z_s, ray.launch, r, false, 0.01);
            CHECK(traced.z == Catch::Approx(z_r).margin(1e-6));
            CHECK(traced.time == Catch::Approx(ray.travel_time).epsilon(1e-9));
            CHECK(traced.length == Catch::Approx(ray.path_length).epsilon(1e-9));
        }
        const auto reflected = layered_reflected_eigenray(layers, z_s, z_r, r, -0.2);
        if (reflected.exists)
        {
            const auto traced = trace(layers, z_s, reflected.launch, r, true, 0.01);
            REQUIRE(traced.grazing > 0.0);
            CHECK(traced.z == Catch::Approx(z_r).margin(1e-6));
            CHECK(traced.grazing == Catch::Approx(reflected.grazing).margin(1e-8));
            CHECK(traced.length == Catch::Approx(reflected.path_length).epsilon(1e-9));
        }
    }
}

TEST_CASE("rays that turn aloft reach into a near-ground shadow", "[stratified]")
{
    // Sound speed falls with height near the ground -- upward refraction, a
    // shadow -- and rises fast aloft, which turns steep rays back down.
    const std::vector<double> z = {0.0, 30.0, 300.0, 1500.0};
    const std::vector<double> c = {344.0, 340.0, 340.0, 380.0};
    const Layers layers {z.data(), c.data(), z.size()};
    const double z_s = 150.0, z_r = 1.2;

    // Outward from the source: branch A's rays, launched down; then a shadow,
    // whose limiting ray is A's last; then rays launched up that turn aloft.
    enum class Seen { a, gap, b };
    Seen seen = Seen::a;
    double first_b = 0.0;
    for (double x = 100.0; x < 30000.0 && first_b == 0.0; x += 100.0)
    {
        const auto ray = layered_direct_eigenray(layers, z_s, z_r, x, -0.1);
        if (seen == Seen::a && !ray.lit)
        {
            seen = Seen::gap;
            CHECK(ray.launch < 0.0);
            CHECK(ray.shadow_depth > 0.0);
        }
        else if (seen == Seen::a)
        {
            CHECK(ray.launch < 0.0);
        }
        else if (seen == Seen::gap && ray.lit)
        {
            seen = Seen::b;
            first_b = x;
        }
    }
    REQUIRE(seen == Seen::b);
    for (const double r : {first_b, first_b + 2000.0})
    {
        INFO("r = " << r);
        const auto far = layered_direct_eigenray(layers, z_s, z_r, r, 0.3);
        REQUIRE(far.lit);
        CHECK(far.launch > 0.0);
        const auto traced = trace(layers, z_s, far.launch, r, false, 0.01);
        CHECK(traced.z == Catch::Approx(z_r).margin(1e-5));
        CHECK(traced.time == Catch::Approx(far.travel_time).epsilon(1e-8));
        CHECK(traced.length == Catch::Approx(far.path_length).epsilon(1e-8));
    }
}

TEST_CASE("the layered spreading is the traced ray tube's", "[stratified]")
{
    const std::vector<double> z = {0.0, 10.0, 60.0, 250.0, 900.0};
    const std::vector<double> c = {340.0, 343.0, 344.5, 343.5, 350.0};
    const Layers layers {z.data(), c.data(), z.size()};
    const double z_s = 180.0, z_r = 1.2, r = 900.0;
    const auto ray = layered_direct_eigenray(layers, z_s, z_r, r, -0.2);
    REQUIRE(ray.lit);
    auto landing = [&](double theta) {
        double lo = 0.5 * r, hi = 2.0 * r;
        for (int i = 0; i < 50; i++)
        {
            const double mid = 0.5 * (lo + hi);
            (trace(layers, z_s, theta, mid, false).z > z_r ? lo : hi) = mid;
        }
        return 0.5 * (lo + hi);
    };
    const double delta = 2e-5;
    const double slope = (landing(ray.launch + delta) - landing(ray.launch - delta)) / (2 * delta);
    const double tube = std::sqrt(r * std::fabs(slope) * std::sin(std::fabs(ray.arrival)) / std::cos(ray.launch));
    CHECK(tube == Catch::Approx(ray.spreading_range).epsilon(1e-3));
}

TEST_CASE("layered_direct_ray in uniform air is the straight ray", "[stratified]")
{
    AtmosphereProfile profile;
    profile.height = {0.0, 200.0, 600.0};
    profile.temperature = {290.0, 290.0, 290.0};
    profile.relative_humidity = {40.0, 40.0, 40.0};
    profile.pressure = {100.0, 97.7, 93.2};
    const LayeredSoundSpeed air(profile, 90.0);
    const EffectiveSoundSpeed linear(profile, 90.0, 0.0, 700.0);
    std::vector<double> scratch;
    const Vec3 source(-400.0, 250.0, 300.0);
    const Vec3 receiver(700.0, -150.0, 0.0);
    const auto ray = layered_direct_ray(air, linear, source, linear.moments(300.0), receiver,
            linear.moments(0.0), scratch);
    const auto straight = straight_ray(source, receiver);
    CHECK(ray.path.range == straight.range);
    CHECK(ray.path.direction.z == straight.direction.z);
    CHECK(ray.shadow_depth == 0.0);
}



TEST_CASE("the log-law upwind shadow starts where Snell's law puts it", "[stratified]")
{
    // With no wind at the ground the air there is the fastest upwind, so the
    // limiting ray leaves the source at cos(theta) = c(z_s) / c(0) whatever the
    // profile's shape, and the boundary is the integral of cot(theta(z)) down
    // to the ground: about 520 m for a source at 100 m in
    // atmospheres/example_windy.csv (5.2 m/s at 10 m).  One line fitted over
    // the path puts it more than twice as far out -- the reason for layers.
    const auto profile = read_atmosphere_profile("atmospheres/example_windy.csv", feet, 0.0)
                                 .with_hydrostatic_pressure(101.325);
    const LayeredSoundSpeed air(profile, 90.0);
    std::vector<double> c(air.levels());
    air.along(-1.0, 0.0, c.data()); // the wind blows east; propagate west
    const Layers layers {air.heights().data(), c.data(), air.levels()};
    const double z_s = 100.0;
    // The same integral by quadrature, straight from the layers.
    const double c_ground = layers.at(0.0);
    double boundary = 0.0;
    const int n = 200000;
    for (int i = 0; i < n; i++)
    {
        const double z = z_s * (i + 0.5) / n;
        const double cosine = layers.at(z) / c_ground;
        boundary += (z_s / n) * cosine / std::sqrt(1.0 - cosine * cosine);
    }
    CHECK(boundary == Catch::Approx(520.0).epsilon(0.02));
    CHECK(layered_direct_eigenray(layers, z_s, 0.0, 0.99 * boundary, -0.2).lit);
    CHECK_FALSE(layered_direct_eigenray(layers, z_s, 0.0, 1.01 * boundary, -0.2).lit);
}

TEST_CASE("a receiver at the shadow boundary has a finite ray tube", "[stratified]")
{
    // Where the lit root closes on the end of its family the tube's window in
    // launch angle has no width, its slope is 0/0, and the ray's spreading range
    // came out NaN with the ray still marked lit.  That NaN went through the
    // level to a NaN footprint cell, written with exit status 0, within about 1e-8
    // of the boundary in range.
    const auto profile = read_atmosphere_profile("atmospheres/example_windy.csv", feet, 0.0)
                                 .with_hydrostatic_pressure(101.325);
    const LayeredSoundSpeed air(profile, 90.0);
    std::size_t lit = 0, boundaries = 0;
    for (int bearing = 0; bearing < 12; bearing++)
    {
        const double angle = bearing * 30.0 * 3.14159265358979323846 / 180.0;
        std::vector<double> c(air.levels());
        air.along(std::cos(angle), std::sin(angle), c.data());
        const Layers layers {air.heights().data(), c.data(), air.levels()};
        for (double z_s : {20.0, 30.48, 60.0, 100.0, 150.0})
        {
            for (double z_r : {0.0, 1.2})
            {
                // The solver's own boundary: far enough out to be in the shadow.
                const auto far = layered_direct_eigenray(layers, z_s, z_r, 30000.0, -0.2);
                if (!std::isfinite(far.boundary) || !(far.boundary > 10.0))
                    continue;
                boundaries++;
                for (double eps : {1e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7})
                {
                    const auto ray = layered_direct_eigenray(layers, z_s, z_r,
                            far.boundary * (1.0 - eps), -0.2);
                    if (!ray.lit)
                        continue;
                    lit++;
                    INFO("bearing " << bearing * 30 << " z_s " << z_s << " z_r " << z_r << " eps " << eps);
                    REQUIRE(std::isfinite(ray.spreading_range));
                    REQUIRE(ray.spreading_range > 0.0);
                }
            }
        }
    }
    REQUIRE(boundaries > 20);
    REQUIRE(lit > 100);
}

namespace {

// A calm morning inversion capped near 107 m, day 286's (run 286276) without its
// wind: c(T) at each level, T in Celsius.
struct Duct
{
    std::vector<double> z, c;
    Duct()
    {
        const double table[][2] = {{0, -0.033}, {3.048, 0.100}, {9.144, 0.400}, {15.24, 0.767}, {30.48, 1.333},
                {45.72, 2.211}, {60.96, 2.556}, {76.2, 2.594}, {91.44, 3.233}, {106.07, 3.423},
                {106.68, 3.417}, {121.92, 3.506}, {137.16, 3.572}, {152.4, 3.539}, {1524.0, 3.539}};
        for (const auto& row : table)
        {
            z.push_back(row[0]);
            c.push_back(331.3 * std::sqrt(1.0 + row[1] / 273.15));
        }
    }
    Layers layers() const
    {
        return Layers {z.data(), c.data(), z.size()};
    }
};

} // namespace

TEST_CASE("every arrival in a surface duct lands where layered_arrivals says", "[stratified]")
{
    // Past the reach of every direct ray, only rays that bounced arrive; nearer,
    // they arrive beside the direct one.  Each must be a real ray: traced with
    // its bounces off z = 0, it reaches the receiver's height at r, in the time
    // and over the length given, meeting the ground at the angle given.
    const Duct duct;
    const Layers layers = duct.layers();
    const double z_s = 30.48;
    const double z_r = 1.2;
    int bounced_arrivals = 0;
    for (const double r : {1800.0, 2600.0, 3400.0})
    {
        const auto direct = layered_direct_eigenray(layers, z_s, z_r, r, -0.05);
        const auto arrivals = layered_arrivals(layers, z_s, z_r, r, direct);
        for (const auto& a : arrivals)
        {
            CAPTURE(r, a.bounces, a.launch);
            const auto traced = trace(layers, z_s, a.launch, r, a.bounces, 0.01);
            CHECK(traced.z == Catch::Approx(z_r).margin(1e-5));
            CHECK(traced.time == Catch::Approx(a.travel_time).epsilon(1e-8));
            CHECK(traced.length == Catch::Approx(a.path_length).epsilon(1e-8));
            if (a.bounces > 0)
            {
                CHECK(traced.grazing == Catch::Approx(a.grazing).margin(1e-7));
                bounced_arrivals++;
            }
            CHECK(a.spreading_range > 0.0);
        }
    }
    CHECK(bounced_arrivals > 0);
}

TEST_CASE("the direct solve reports no ray that does not land", "[stratified]")
{
    // A family whose end ray never reaches the receiver brackets falsely: that
    // end stands in as an infinity of the right sign.  Through this duct (day
    // 286's run 286276, calm) the iteration closed on the duct's top from 2.4 km
    // out and reported it lit, its tube 575 times the range, 55 dB low.  Lit has
    // to mean landed.
    auto profile = read_atmosphere_profile("atmospheres/example_inversion.csv", feet, 0.0);
    profile = profile.with_hydrostatic_pressure(101.325);
    const LayeredSoundSpeed air(profile, 90.0);
    std::vector<double> c(air.levels());
    air.along(-1.0, 0.0, c.data());
    const Layers layers {air.heights().data(), c.data(), air.levels()};
    const double z_s = 100.0 * ft2m;
    const double z_r = 1.2;
    int dark = 0;
    for (double r = 600.0; r <= 4000.0; r += 200.0)
    {
        const auto ray = layered_direct_eigenray(layers, z_s, z_r, r, -0.05);
        CAPTURE(r, ray.launch);
        if (ray.lit)
        {
            const auto traced = trace(layers, z_s, ray.launch, r, 0, 0.02);
            CHECK(traced.z == Catch::Approx(z_r).margin(1e-5));
        }
        else
        {
            dark++;
        }
    }
    // Past the single-hop reach the direct solve finds nothing; what arrives
    // there is layered_arrivals'.
    CHECK(dark > 0);
}
