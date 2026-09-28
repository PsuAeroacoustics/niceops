// The field through a shadow boundary (ShadowTransition), and the zeros past
// |q| = 6 that building it needs.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "airy.hpp"
#include "exact_solutions.hpp"
#include "ground_effect.hpp"
#include "refraction.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT

namespace {

using complex = std::complex<double>;

const complex w = std::polar(1.0, 2.0 * pi / 3.0);
const complex soft = 1.0 / complex(5.96, 2.46); // the benchmark's 1 kHz ground

// The ray model's level re free field for a point receiver: spreading, the
// ground term with the grazing ray's straight equivalent geometry, and -- if
// the transition takes the ray -- its weight and loss in place of the
// creeping-wave loss.
double ray_level(const LinearProfile& profile, double z_s, double z_r, double r, double f,
        complex beta, const ShadowTransition* transition)
{
    RefractedRay ray;
    ray.curved = true;
    ray.profile = profile;
    ray.horizontal = r;
    ray.source_z = z_s;
    ray.direct = direct_eigenray(profile, z_s, z_r, r);
    ray.shadow_depth = ray.direct.shadow_depth;
    const auto reflected = ray.direct.lit ? reflected_eigenray(profile, z_s, z_r, r) : ReflectedEigenray {};
    const double grazing = reflected.exists ? reflected.grazing : 0.0;
    const double length = reflected.exists ? reflected.path_length : ray.direct.path_length;
    GroundModel ground;
    ground.receiver_height = z_r;
    double ground_db = 0.0;
    excess_ground_attenuation(ground, {beta}, std::max(length * std::sin(grazing) - z_r, 0.0),
            length * std::cos(grazing), profile.sound_speed, {f}, &ground_db, 0, 1, false);
    const double spreading = 20.0 * std::log10(std::hypot(r, z_s - z_r) / ray.direct.spreading_range);
    const ShadowAttenuation attenuation({f}, {beta}, 0, 1);
    double weight = 1.0, loss = 0.0;
    if (transition == nullptr || !transition->apply(ray, attenuation, 20.0, 0, 1, &weight, &loss))
        loss = ray.shadow_depth > 0.0 ? std::min(attenuation.rate(0, profile) * ray.shadow_depth, 20.0) : 0.0;
    return spreading + weight * ground_db - loss;
}

} // namespace

TEST_CASE("past |q| = 6 the seeded zeros carry on from the followed ones", "[shadow]")
{
    for (const complex beta : {soft, 1.0 / complex(12.81, 11.62), complex(0.3, -0.02)})
    {
        const complex d = complex(0.0, -1.0) * w * w * beta / std::abs(beta);
        const auto inside = impedance_airy_zeros(5.999 * d, 30);
        const auto outside = impedance_airy_zeros(6.001 * d, 30);
        for (std::size_t n = 0; n < 30; n++)
        {
            INFO("beta " << beta << ", zero " << n + 1);
            CHECK(std::abs(inside[n] - outside[n]) < 1e-2);
        }
        // And far out, each is a zero of Ai' + q Ai.
        for (const double size : {30.0, 1e3, 1e5})
            for (const complex t : impedance_airy_zeros(size * d, 50))
            {
                complex ai, aip;
                airy(t, ai, aip);
                CHECK(std::abs(aip + size * d * ai) <= 1e-8 * (std::abs(aip) + std::abs(size * d * ai)));
            }
    }
    // A downward problem past 6 has a surface wave, which seeding cannot place.
    CHECK_THROWS(impedance_airy_zeros(complex(5.0, 5.0), 5));
}

TEST_CASE("the residue series runs on past |q| = 6", "[shadow]")
{
    // The same geometry at two frequencies either side of |q| = 6: the field
    // moves smoothly from one to the other.
    constexpr double c0 = 343.23, a = 0.03;
    std::vector<double> levels;
    std::vector<double> sizes;
    for (const double f : {548.0, 553.0, 558.0})
    {
        ExactProblem p;
        p.frequency = f;
        p.sound_speed = c0;
        p.admittance = soft;
        p.source_height = 30.0;
        p.receiver_height = 1.2;
        const double k0 = 2.0 * pi * f / c0;
        const double q = std::cbrt(c0 / a / (2.0 * k0 * k0)) * k0 * std::abs(soft);
        sizes.push_back(q);
        INFO(f << " Hz: |q| = " << q);
        const auto field = upward_pressure(p, c0 / a, {1600.0});
        REQUIRE(std::isfinite(field[0].real()));
        levels.push_back(20.0 * std::log10(std::abs(field[0]) * 1600.0));
    }
    CHECK(sizes.front() < 6.0);
    CHECK(sizes.back() > 6.0);
    CHECK(std::fabs(levels[1] - 0.5 * (levels[0] + levels[2])) < 0.3);
}

TEST_CASE("through the boundary the transition follows the exact field", "[shadow]")
{
    // Tabulated with the source 20 l up, applied to sources at 30 and 100 m:
    // within 1.5 dB of the residue series from xi = -1.5 to 2, where the ray
    // model's ground term was 10-25 dB off over soft ground (measured 0.1-1 dB;
    // tests/test_shadow_attenuation.cpp covers the rate beyond).
    constexpr double c0 = 343.23, z_r = 1.2;
    const std::vector<double> bands = {250.0, 1000.0};
    for (const complex beta : {soft, complex(1e-9, 0.0)})
    {
        for (const double a : {0.1, 0.03})
            for (std::size_t b = 0; b < bands.size(); b++)
            {
                const double f = bands[b];
                const double k0 = 2.0 * pi * f / c0;
                const double s = std::cbrt(pi * f / a);
                const double diffraction = 2.0 * s * s / k0;
                const ShadowTransition single({f}, {beta}, z_r, c0, 0, 1);
                for (const double z_s : {30.0, 100.0})
                {
                    const LinearProfile profile {c0, -a};
                    const double r_b = shadow_boundary(profile, z_s, z_r);
                    const double exact_r_b = std::sqrt(2.0 * c0 / a * z_s) + std::sqrt(2.0 * c0 / a * z_r);
                    // Past V(0) - 20 dB the model holds at its floor by design,
                    // so compare only above it.
                    double at_boundary = std::numeric_limits<double>::quiet_NaN();
                    for (const double xi : {0.0, -1.5, -1.0, -0.5, 0.5, 1.0, 2.0})
                    {
                        ExactProblem p;
                        p.frequency = f;
                        p.sound_speed = c0;
                        p.admittance = beta;
                        p.source_height = z_s;
                        p.receiver_height = z_r;
                        const double r_exact = exact_r_b + xi * diffraction;
                        const auto field = upward_pressure(p, c0 / a, {r_exact});
                        if (!std::isfinite(field[0].real()))
                            continue;
                        const double exact = 20.0 * std::log10(std::abs(field[0]) * std::hypot(r_exact, z_s - z_r));
                        if (xi == 0.0)
                            at_boundary = exact;
                        if (exact < at_boundary - 18.0)
                            continue;
                        const double model = ray_level(profile, z_s, z_r, r_b + xi * diffraction, f, beta, &single);
                        INFO((beta == soft ? "soft" : "rigid") << ", " << f << " Hz, a " << a << ", z_s " << z_s
                                << ", xi " << xi << ": exact " << exact << ", transition " << model);
                        CHECK(std::fabs(model - exact) < 1.5);
                    }
                }
            }
    }
}

TEST_CASE("the transition leaves lit rays alone and keeps to its floor and cap", "[shadow]")
{
    constexpr double c0 = 343.0, z_r = 1.2, a = 0.05;
    const std::vector<double> bands = {125.0, 500.0, 2000.0};
    const std::vector<complex> beta(bands.size(), soft);
    const ShadowTransition transition(bands, beta, z_r, c0, 0, bands.size());
    const ShadowAttenuation attenuation(bands, beta, 0, bands.size());
    const LinearProfile profile {c0, -a};
    constexpr double z_s = 80.0;
    const double r_b = shadow_boundary(profile, z_s, z_r);
    auto make = [&](double r) {
        RefractedRay ray;
        ray.curved = true;
        ray.profile = profile;
        ray.horizontal = r;
        ray.source_z = z_s;
        ray.direct = direct_eigenray(profile, z_s, z_r, r);
        ray.shadow_depth = ray.direct.shadow_depth;
        return ray;
    };
    std::vector<double> weight(bands.size()), loss(bands.size());

    // More than two diffraction lengths short of the boundary in every band
    // (the lowest band's is the longest): untouched.
    const double k_lo = 2.0 * pi * bands[0] / c0;
    const double s_lo = std::cbrt(pi * bands[0] / a);
    CHECK_FALSE(transition.apply(make(r_b - 2.5 * 2.0 * s_lo * s_lo / k_lo), attenuation, 20.0, 0,
            bands.size(), weight.data(), loss.data()));
    // A straight ray, or one refracting downward, has no boundary.
    RefractedRay down = make(0.5 * r_b);
    down.profile.gradient = a;
    down.direct = direct_eigenray(down.profile, z_s, z_r, 0.5 * r_b);
    CHECK_FALSE(transition.apply(down, attenuation, 20.0, 0, bands.size(), weight.data(), loss.data()));

    // Deep in the shadow: all V, no ground term, and never below V(0) - floor.
    for (const double floor : {10.0, 20.0})
    {
        REQUIRE(transition.apply(make(r_b + 5000.0), attenuation, floor, 0, bands.size(), weight.data(),
                loss.data()));
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            const double s = std::cbrt(pi * bands[b] / a);
            INFO(bands[b] << " Hz, floor " << floor);
            CHECK(weight[b] == 0.0);
            CHECK(-loss[b] == Catch::Approx(transition.level(b, 0.0, s) - floor).margin(0.05));
        }
    }
    // The table never goes above pressure doubling.
    for (std::size_t b = 0; b < bands.size(); b++)
        CHECK(transition.max_power_ratio(b) <= 4.0 * (1.0 + 1e-12));
}
