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
    // which over a kilometer of ray is a tenth of a millimeter.
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

namespace {

// The level re free field of a pair's arrivals summed in energy, the primary
// counted as free field if it is lit.  The floor stands for "nothing arrives".
double arrivals_level(const Layers& layers, double z_s, double z_r, double r, int samples)
{
    const auto direct = layered_direct_eigenray(layers, z_s, z_r, r, std::atan2(z_r - z_s, r));
    const auto arrivals = layered_arrivals(layers, z_s, z_r, r, direct, false, samples);
    const double chord = std::hypot(r, z_s - z_r);
    double energy = direct.lit ? 1.0 : 0.0;
    for (const auto& a : arrivals)
        energy += std::pow(chord / a.spreading_range, 2);
    return 10.0 * std::log10(std::max(energy, 1e-6));
}

} // namespace

TEST_CASE("the arrival search has converged in its sample count", "[stratified]")
{
    // A ray that bounces m times travels X_0(p) + m L(p), which turns through m
    // times as many wiggles per sample interval as a direct one, so a fixed 32
    // samples lost multi-hop arrivals, whole rays that arrive nowhere else: over
    // the pairs below, 10 of 385 came out more than 1 dB from the answer at 512
    // samples (11% of a wider sweep, with the grid not yet bunched toward the
    // shallow end as it now is).  The count now grows with the hops the path can
    // hold, and 3 of the 385 are left.  Sources of 10-60 m, ranges to 8 km:
    const Duct duct;
    const Layers layers = duct.layers();
    int pairs = 0, off_by_1db = 0;
    for (const double z_s : {10.0, 20.0, 30.48, 45.0, 60.0})
    {
        for (double r = 400.0; r <= 8000.0; r += 100.0)
        {
            const double got = arrivals_level(layers, z_s, 1.2, r, 0);
            const double converged = arrivals_level(layers, z_s, 1.2, r, 512);
            pairs++;
            off_by_1db += std::fabs(got - converged) > 1.0;
        }
    }
    INFO(off_by_1db << " of " << pairs << " pairs differ by more than 1 dB");
    CHECK(off_by_1db <= 0.015 * pairs);
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

namespace {

// Day 286's last profile toward 270 (run 286288): fastest at the ground, where
// the log-law headwind assumed below the lidar's first gate takes 0.2 m/s off in
// the lowest foot, slowest near 230 ft, and back to nearly the ground's speed at
// the inversion's top near 330 ft.  Heights in feet.
struct Elevated
{
    std::vector<double> z, c;
    Elevated()
    {
        const double ft[] = {0, 1, 3, 10, 20, 30, 46, 70, 100, 120, 160, 230, 280, 330, 400, 500};
        const double speed[] = {334.125, 333.918, 333.835, 333.796, 333.703, 333.648, 333.490, 333.296,
                333.367, 333.409, 333.259, 333.017, 333.369, 333.921, 333.669, 333.516};
        for (std::size_t k = 0; k < std::size(ft); k++)
        {
            z.push_back(ft[k] * ft2m);
            c.push_back(speed[k]);
        }
    }
    Layers layers() const { return {z.data(), c.data(), z.size()}; }
};

} // namespace

TEST_CASE("tunneled arrivals carry the evanescent integral below their turning point", "[stratified]")
{
    // A 100 ft source 12,000 ft out: the duct's return, whose rays turn back up
    // a few feet above a flush microphone.  Each one's integral is checked
    // against a brute-force quadrature of sqrt(p^2 - 1/c^2) where it is real.
    const Elevated profile;
    const Layers layers = profile.layers();
    const double z_s = 100.0 * ft2m, z_r = 0.0, r = 12000.0 * ft2m;
    CHECK(layered_wave_arrivals(layers, z_s, z_r, r, {}).empty());
    const auto arrivals = layered_wave_arrivals(layers, z_s, z_r, r, {true, false});
    REQUIRE(!arrivals.empty());
    const double c_s = layers.at(z_s);
    int checked = 0;
    for (const Arrival& a : arrivals)
    {
        if (!(a.wave.evanescent > 0.0) || a.wave.caustic)
            continue;
        const double p = std::cos(a.launch) / c_s;
        // The turning point: going down from the source, where c first reaches 1/p.
        double z_t = z_s;
        while (z_t > 0.0 && layers.at(z_t) < 1.0 / p)
            z_t -= 1e-5;
        REQUIRE(z_t > z_r);
        constexpr int n = 200000;
        double sum = 0.0;
        for (int k = 0; k < n; k++)
        {
            const double z = z_r + (k + 0.5) * (z_t - z_r) / n;
            const double c = layers.at(z);
            sum += std::sqrt(std::max(0.0, p * p - 1.0 / (c * c)));
        }
        sum *= (z_t - z_r) / n;
        CAPTURE(a.launch, z_t);
        CHECK(a.wave.evanescent == Catch::Approx(sum).epsilon(1e-3));
        CHECK(a.grazing == 0.0);
        checked++;
    }
    CHECK(checked > 0);
}

TEST_CASE("a wave arrival's loss is never negative", "[stratified]")
{
    // Never negative, so the peak search's bound on the tube bounds the band.
    // Not monotonic in frequency at small depth -- the energy near a turning
    // point grows as omega^(1/3) while a shallow tail hardly decays -- but a
    // deep tail costs the top band more than the bottom one.
    const std::vector<double> bands = {50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250,
            1600, 2000, 2500, 3150, 4000, 5000};
    for (const double evanescent : {1e-5, 1e-4, 1e-3})
    {
        WaveCorrection wave;
        wave.evanescent = evanescent;
        std::vector<double> loss(bands.size(), 0.0);
        add_wave_loss(wave, 1000.0, 1000.0, bands, 0, bands.size(), loss.data());
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            CAPTURE(evanescent, bands[b]);
            CHECK(loss[b] >= 0.0);
        }
        if (evanescent >= 1e-4)
            CHECK(loss.back() > loss.front());
    }
    // A caustic's dark side, and its lit side, where the band averages the fringes.
    for (const double depth : {-1e-3, 0.0, 1e-3, 1e-2})
    {
        WaveCorrection wave;
        wave.caustic = true;
        wave.caustic_depth = depth;
        std::vector<double> loss(bands.size(), 0.0);
        add_wave_loss(wave, 1000.0, 1000.0, bands, 0, bands.size(), loss.data());
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            CAPTURE(depth, bands[b]);
            CHECK(loss[b] >= 0.0);
        }
        if (depth > 0.0)
            CHECK(loss.back() > loss.front());
    }
    // A tube focused past the cap is held to it.
    WaveCorrection wave;
    wave.evanescent = 1e-9;
    std::vector<double> loss(1, 0.0);
    const std::vector<double> one = {wave_reference_frequency};
    add_wave_loss(wave, 1000.0, 10.0, one, 0, 1, loss.data());
    CHECK(loss[0] == Catch::Approx(40.0 - max_focusing_gain_db).margin(1e-6));
}

TEST_CASE("a deep duct retains more than 24 arrivals", "[stratified]")
{
    const Duct duct;
    const Layers layers = duct.layers();
    std::size_t most = 0;
    for (double source : {1.0, 2.0, 5.0, 10.0})
        for (double range : {4000.0, 8000.0, 16000.0, 32000.0})
        {
            const auto primary = layered_direct_eigenray(layers, source, 1.2, range,
                    std::atan2(1.2 - source, range));
            const auto arrivals = layered_arrivals(layers, source, 1.2, range, primary, false, 512);
            most = std::max(most, arrivals.size());
            for (const auto& ray : arrivals)
            {
                REQUIRE(std::isfinite(ray.spreading_range));
                REQUIRE(ray.spreading_range > 0.0);
            }
        }
    INFO("most arrivals retained: " << most);
    REQUIRE(most > 24);
}

namespace {

// An upwind log law like day 284's: a 7 kt headwind at 6.6 ft over z0 = 0.003
// ft, and the dry-adiabatic lapse, as effective sound speed along the path.
// With `floor` > 0, the air below `floor` meters is calm and isothermal, as a
// log law written with zero wind below its roughness length leaves it; with
// `extra` levels inside that floor too, describing the same air.
struct Upwind
{
    std::vector<double> z, c;
    Upwind(double floor = 0.0, int extra = 0, double headwind = 7.0 * knot)
    {
        const double ft[] = {0, 1, 3, 6.6, 10, 20, 33, 46, 80, 120, 200, 300, 400, 500, 650, 1000};
        auto speed = [&](double height) {
            const double t = 283.15 - 0.0098 * height;
            const double z0 = 0.003 * ft2m, ref = 6.6 * ft2m;
            const double u = height > z0 ? headwind * std::log(height / z0) / std::log(ref / z0) : 0.0;
            return 20.05 * std::sqrt(t) - u;
        };
        z.push_back(0.0);
        c.push_back(speed(0.0));
        for (int k = 1; k <= extra; k++)
        {
            z.push_back(floor * k / (extra + 1.0));
            c.push_back(c[0]);
        }
        if (floor > 0.0)
        {
            z.push_back(floor);
            c.push_back(c[0]);
        }
        for (double h : ft)
        {
            if (h * ft2m > floor)
            {
                z.push_back(h * ft2m);
                c.push_back(speed(h * ft2m));
            }
        }
    }
    Layers layers(double frequency = default_resolution_frequency) const
    {
        return Layers {z.data(), c.data(), z.size(), frequency};
    }
};

} // namespace

TEST_CASE("a calm floor thinner than a wavelength leaves the shadow boundary where it was", "[stratified]")
{
    // A log law written with no wind below its roughness length, on day 284's
    // wind, leaves a calm, isothermal floor of 0.1 ft: as fast as the ground,
    // so the limiting ray of the upwind shadow ran along it.  The downward
    // family's last ray (1e-12 short of the ground's invariant) crossed it at a
    // microradian, which put the boundary 20 km out, a figure set by the inset
    // alone; flush microphones on it came out 22-28 dB low.  No band sees a
    // floor that thin.  Plain ray theory is resolution_frequency 0.
    const Upwind open;
    for (double floor : {0.001, 0.1 * ft2m})
    {
        const Upwind calm(floor);
        for (double z_s : {30.0, 100.0})
        {
            for (double z_r : {0.0, 0.5 * floor})
            {
                CAPTURE(floor, z_s, z_r);
                const double b0 = layered_direct_eigenray(open.layers(), z_s, z_r, 1e6, -0.2).boundary;
                const double b = layered_direct_eigenray(calm.layers(), z_s, z_r, 1e6, -0.2).boundary;
                REQUIRE(std::isfinite(b0));
                REQUIRE(std::isfinite(b));
                // The floor may stretch the boundary by the cap's slack over
                // what ray theory travels in a resolved turning, (1.5 - 4/3)
                // k h^2: 1.4 m of 186 at 10 kHz, the most.
                CHECK(std::fabs(b - b0) < 0.012 * b0);
                for (double f : {0.8, 0.95, 1.05, 1.2})
                {
                    CAPTURE(f);
                    const auto with = layered_direct_eigenray(calm.layers(), z_s, z_r, f * b0, -0.2);
                    const auto without = layered_direct_eigenray(open.layers(), z_s, z_r, f * b0, -0.2);
                    CHECK(with.lit == without.lit);
                    // Short of the hand-over (handover_fraction), where the
                    // layered ray's own spreading is what a run takes.
                    if (f < handover_fraction)
                        CHECK(std::fabs(20.0 * std::log10(with.spreading_range / without.spreading_range)) < 0.25);
                }
                if (floor > 0.01 && z_r == 0.0)
                    CHECK(layered_direct_eigenray(calm.layers(0.0), z_s, z_r, 1e6, -0.2).boundary > 20.0 * b0);
            }
        }
    }
}

namespace {

// Every field the layered solve reports, compared exactly.
bool same(const DirectEigenray& a, const DirectEigenray& b)
{
    auto eq = [](double x, double y) { return x == y || (std::isnan(x) && std::isnan(y)); };
    return a.lit == b.lit && eq(a.launch, b.launch) && eq(a.arrival, b.arrival) && eq(a.spreading_range, b.spreading_range)
            && eq(a.path_length, b.path_length) && eq(a.travel_time, b.travel_time) && eq(a.boundary, b.boundary)
            && eq(a.shadow_depth, b.shadow_depth);
}

} // namespace

TEST_CASE("the resolution limit leaves every resolved turning alone", "[stratified]")
{
    // In a linear gradient a ray grazing the fastest air travels (4/3) k h^2
    // across the heights within one radian of it, 8/9 of the cap, so the limit
    // must never bind there.  A level a millimeter off the fastest air, on the
    // same line, keeps the walk from deciding that by its shortcut; it has to
    // integrate.  With a cap of k h^2 every one of these would bind.
    for (double gradient : {-5.0, -0.5, -0.05, 0.05, 0.5})
    {
        const double c0 = 340.0;
        const std::vector<double> z {0.0, 1e-3, 0.4, 300.0};
        std::vector<double> c;
        for (double h : z)
            c.push_back(c0 + gradient * h);
        for (double frequency : {1000.0, 4000.0, 10000.0})
        {
            for (double z_s : {30.0, 100.0})
            {
                for (double z_r : {0.0, 1e-3, 1.2})
                {
                    const Layers resolved {z.data(), c.data(), z.size(), frequency};
                    const Layers plain {z.data(), c.data(), z.size(), 0.0};
                    for (double r : {50.0, 200.0, 800.0, 3000.0})
                    {
                        CAPTURE(gradient, frequency, z_s, z_r, r);
                        CHECK(same(layered_direct_eigenray(resolved, z_s, z_r, r, -0.2),
                                layered_direct_eigenray(plain, z_s, z_r, r, -0.2)));
                    }
                }
            }
        }
    }
    // Nor on a log law without a floor, whose lowest layers are steep.
    const Upwind open;
    for (double z_s : {10.0, 30.0, 100.0})
        for (double z_r : {0.0, 0.02, 1.2})
            for (double r : {100.0, 300.0, 600.0, 2000.0})
            {
                CAPTURE(z_s, z_r, r);
                CHECK(same(layered_direct_eigenray(open.layers(), z_s, z_r, r, -0.2),
                        layered_direct_eigenray(open.layers(0.0), z_s, z_r, r, -0.2)));
            }
}

TEST_CASE("a path shorter than the turning region keeps its rays", "[stratified]")
{
    // The turning region belongs to the air, not to the path.  Cut off at the
    // path's ends, a low source's region was as short as the path and the cap
    // bound on smooth gradients: under a plain lapse a 1.5 m source over a
    // 1.2 m microphone had its boundary brought in from 789 m to 34 m at
    // 10 kHz, and at 200 Hz sources to 10 m lost theirs in any light wind.
    // Under a tailwind the source's own height is the fastest air and the
    // downward family's end is the horizontal launch, where the families
    // launched up to turn aloft take over; capping it opened a gap between
    // them, receivers in the lit zone with no ray.  None of these resolve
    // anything differently from plain ray theory.
    for (double wind : {-1.0, -0.1, 0.0, 0.03, 0.3, 3.6})
    {
        const Upwind air(0.0, 0, wind);
        for (double frequency : {200.0, 1000.0, 10000.0})
            for (double z_s : {1.5, 3.0, 5.0, 10.0})
                for (double z_r : {0.0, 1.2})
                    for (double r : {30.0, 100.0, 300.0, 800.0, 1e7})
                    {
                        CAPTURE(wind, frequency, z_s, z_r, r);
                        CHECK(same(layered_direct_eigenray(air.layers(frequency), z_s, z_r, r, -0.05),
                                layered_direct_eigenray(air.layers(0.0), z_s, z_r, r, -0.05)));
                    }
    }
}

TEST_CASE("a flat floor's boundary does not depend on the levels that describe it", "[stratified]")
{
    // The turning region is a set of heights, found by integrating the vertical
    // phase, not a count of layers: the same air described by more levels
    // gives the same boundary.
    for (double floor : {0.1 * ft2m, 1.0})
    {
        const Upwind few(floor), many(floor, 4);
        for (double z_s : {30.0, 100.0})
        {
            for (double z_r : {0.0, 0.25 * floor})
            {
                CAPTURE(floor, z_s, z_r);
                const double a = layered_direct_eigenray(few.layers(), z_s, z_r, 1e6, -0.2).boundary;
                const double b = layered_direct_eigenray(many.layers(), z_s, z_r, 1e6, -0.2).boundary;
                REQUIRE(std::isfinite(a));
                CHECK(b == Catch::Approx(a).epsilon(1e-6));
            }
        }
    }
}

TEST_CASE("a thick flat floor is followed as far as the highest band resolves it", "[stratified]")
{
    // A meter of calm, isothermal air is ten wavelengths thick at 10 kHz and a
    // third of one at 100 Hz.  A wave confined to it stays there for about
    // k h^2, so the higher the run's highest band, the farther the rays along
    // it are followed -- never to plain ray theory's 700 km.
    const Upwind calm(1.0), open;
    for (double z_s : {30.0, 100.0})
    {
        CAPTURE(z_s);
        const double b0 = layered_direct_eigenray(open.layers(), z_s, 0.0, 1e7, -0.2).boundary;
        double previous = b0;
        for (double frequency : {100.0, 1000.0, 10000.0})
        {
            CAPTURE(frequency);
            const double b = layered_direct_eigenray(calm.layers(frequency), z_s, 0.0, 1e7, -0.2).boundary;
            CHECK(b > previous);
            // At most the cap over the floor and a turning region above it no
            // taller than the floor again.
            const double k = 2.0 * pi * frequency / calm.c[0];
            CHECK(b - b0 < 1.5 * k * 4.0 + 100.0);
            previous = b;
        }
        CHECK(layered_direct_eigenray(calm.layers(0.0), z_s, 0.0, 1e7, -0.2).boundary > 1e5);
    }
}

TEST_CASE("a run's resolution frequency reaches every path's layers", "[stratified]")
{
    // niceops hands LayeredSoundSpeed its highest band, and each path's Layers
    // (path_layers) must carry it to resolved_pass, or a run quietly takes the
    // 10 kHz default whatever its bands.  Through layered_direct_ray, on
    // atmospheres/example_windy_calm_floor.csv: plain ray theory (0) still runs
    // the upwind boundary along the floor, the default puts it back where the
    // profile without the floor has it.
    const auto read = [](const char* path) {
        return read_atmosphere_profile(path, feet, 0.0).with_hydrostatic_pressure(101.325);
    };
    const auto floored = read("atmospheres/example_windy_calm_floor.csv");
    const auto open = read("atmospheres/example_windy.csv");
    const EffectiveSoundSpeed linear(open, 90.0, 0.0, open.top() + 20.0);
    constexpr double z_s = 100.0, r = 400.0;
    const Vec3 source(0.0, 0.0, z_s);
    const Vec3 receiver(-r, 0.0, 0.0); // the wind blows east; propagate west
    std::vector<double> scratch;
    auto boundary = [&](const LayeredSoundSpeed& air) {
        const auto ray = layered_direct_ray(air, linear, source, linear.moments(z_s), receiver,
                linear.moments(0.0), scratch);
        return ray.direct.boundary;
    };
    REQUIRE(LayeredSoundSpeed(floored, 90.0).resolution_frequency() == default_resolution_frequency);
    REQUIRE(LayeredSoundSpeed(floored, 90.0, 0.0).resolution_frequency() == 0.0);
    const double b0 = boundary(LayeredSoundSpeed(open, 90.0));
    const double plain = boundary(LayeredSoundSpeed(floored, 90.0, 0.0));
    const double resolved = boundary(LayeredSoundSpeed(floored, 90.0));
    const double low = boundary(LayeredSoundSpeed(floored, 90.0, 200.0));
    INFO("boundary without the floor " << b0 << " m; with it, plain " << plain << ", 10 kHz "
                                       << resolved << ", 200 Hz " << low);
    REQUIRE(std::isfinite(b0));
    CHECK(plain > 20.0 * b0);
    CHECK(std::fabs(resolved - b0) < 0.025 * b0);
    CHECK(low != resolved);
}
