// The field through a shadow boundary (ShadowTransition), and the zeros past
// |q| = 6 that building it needs.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "airy.hpp"
#include "atmosphere_profile.hpp"
#include "exact_solutions.hpp"
#include "ground_effect.hpp"
#include "plate_receiver.hpp"
#include "refraction.hpp"
#include "stratified.hpp"
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
    return spreading + (weight > 0.0 ? weight * ground_db : 0.0) - loss;
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

TEST_CASE("each band's shadow boundary lies between the layered one and the fit's", "[shadow]")
{
    // example_windy.csv's upwind geometry from 100 m: the layered boundary is
    // 535 m, the linear fit's 1620 m.
    constexpr double r_s = 535.0, r_l = 1620.0, z_s = 100.0, z_r = 1.2, c0 = 340.0;
    double previous = r_l;
    for (const double f : {63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0})
    {
        const double b = ShadowTransition::band_boundary(r_s, r_l, z_s, z_r, f, c0);
        INFO(f << " Hz: " << b << " m");
        CHECK(b > r_s);
        CHECK(b < previous);
        previous = b;
    }
    // Toward the fit's as the diffraction length outgrows the gap, toward the
    // layered one as it shrinks.
    CHECK(ShadowTransition::band_boundary(r_s, r_l, z_s, z_r, 1.0, c0) > r_s + 0.8 * (r_l - r_s));
    CHECK(ShadowTransition::band_boundary(r_s, r_l, z_s, z_r, 1e5, c0) < r_s + 0.15 * (r_l - r_s));
    // One boundary is itself.
    CHECK(ShadowTransition::band_boundary(r_s, r_s, z_s, z_r, 500.0, c0) == r_s);
    CHECK(ShadowTransition::band_boundary(r_s, std::numeric_limits<double>::infinity(), z_s, z_r,
                  500.0, c0) == r_s);
}

TEST_CASE("through a log-law shadow, stratified with the boundary law follows the PE", "[shadow]")
{
    // niceops_pe_compare's stratified column, formed the way it forms it, against
    // the PE's levels from that tool on atmospheres/example_windy.csv: upwind
    // (toward 270), source 100 m, receiver 1.2 m, grass.  Over the shadow edge,
    // from the PE's last level above -3 dB to -25 dB, the model is within 2.7 dB
    // rms at 500 Hz and 5.5 dB at 1 kHz; before the boundary law and the field
    // near the boundary it was 30 and 34 dB off.
    auto profile = read_atmosphere_profile("atmospheres/example_windy.csv", feet, 0.0);
    profile = profile.with_hydrostatic_pressure(reference_pressure);
    const double c0 = profile.at(0.0).sound_speed();
    constexpr double z_s = 100.0, z_r = 1.2;
    const double bearing = 270.0 * pi / 180.0;
    const double hx = std::sin(bearing), hy = std::cos(bearing);
    const LayeredSoundSpeed layered(profile, 90.0);
    const EffectiveSoundSpeed linear(profile, 90.0, 0.0,
            std::max(profile.top(), z_s) + 2.0 * EffectiveSoundSpeed::min_fit_span);
    struct Band
    {
        double frequency, first_range, tolerance;
        std::vector<double> pe;
    };
    const std::vector<Band> bands = {
        {500.0, 500.0, 4.0,
                {0.49, 0.92, 1.41, 1.61, 1.38, 0.65, -0.59, -2.37, -4.75, -7.81, -11.7, -16.44, -21.65}},
        {1000.0, 500.0, 7.0, {4.64, 6.01, 6.27, 5.63, 4.22, 1.98, -1.25, -5.35, -10.19, -15.89, -22.19}},
    };
    GroundModel ground;
    ground.receiver_height = z_r;
    for (const Band& band : bands)
    {
        const std::vector<double> f = {band.frequency};
        const auto admittance = ground_impedance(ground, f, c0);
        const ShadowAttenuation attenuation(f, admittance, 0, 1);
        const ShadowTransition transition(f, admittance, z_r, c0, 0, 1);
        const Vec3 source(0.0, 0.0, z_s);
        std::vector<double> scratch;
        double sum = 0.0;
        for (std::size_t i = 0; i < band.pe.size(); i++)
        {
            const double r = band.first_range + 40.0 * static_cast<double>(i);
            const Vec3 receiver(r * hx, r * hy, z_r);
            RefractedRay ray = layered_direct_ray(layered, linear, source, linear.moments(z_s), receiver,
                    linear.moments(z_r), scratch, true);
            layered_finish_ground(layered, ray, GroundGeometry::point, z_r, scratch);
            const double ground_db = excess_ground_attenuation(ground, ray.path.source_height,
                    ray.path.ground_distance, c0, f)[0];
            const double spreading = 20.0 * std::log10(std::hypot(r, z_s - z_r) / ray.path.range);
            double weight = 1.0, loss = 0.0;
            REQUIRE(transition.apply(ray, attenuation, 20.0, 0, 1, &weight, &loss));
            const double level = spreading + (weight > 0.0 ? weight * ground_db : 0.0) - loss;
            INFO(band.frequency << " Hz, " << r << " m: PE " << band.pe[i] << ", model " << level);
            sum += (level - band.pe[i]) * (level - band.pe[i]);
        }
        const double rms = std::sqrt(sum / static_cast<double>(band.pe.size()));
        INFO(band.frequency << " Hz: rms " << rms << " dB");
        CHECK(rms < band.tolerance);
    }
}

TEST_CASE("a microphone just under the plane is solved for on it", "[shadow]")
{
    // Ground boards sit a few centimetres below the survey datum.  The layers
    // hold their lowest speed below z = 0 and the direct-ray solve knows no
    // ground, so it used to send a ray grazing along under the plane to such a
    // microphone: 33 dB low on day 284 of the 2017 test.  Now it is the pair on
    // the plane.
    auto profile = read_atmosphere_profile("atmospheres/example_windy.csv", feet, 0.0);
    profile = profile.with_hydrostatic_pressure(reference_pressure);
    const LayeredSoundSpeed layered(profile, 90.0);
    const EffectiveSoundSpeed linear(profile, 90.0, -1.0, profile.top() + 20.0);
    constexpr double z_s = 60.0;
    const Vec3 source(0.0, 0.0, z_s);
    std::vector<double> scratch;
    for (const double bearing : {270.0, 90.0})
        for (const double r : {300.0, 800.0, 1500.0})
        {
            const double b = bearing * pi / 180.0;
            auto ray_to = [&](double z) {
                const Vec3 receiver(r * std::sin(b), r * std::cos(b), z);
                return layered_direct_ray(layered, linear, source, linear.moments(z_s), receiver,
                        linear.moments(z), scratch);
            };
            const auto below = ray_to(-0.027);
            const auto on = ray_to(0.0);
            INFO("toward " << bearing << ", " << r << " m: spreading " << below.path.range << " against "
                           << on.path.range);
            CHECK(below.direct.lit == on.direct.lit);
            CHECK(below.path.range == Catch::Approx(on.path.range).epsilon(1e-6));
            CHECK(below.shadow_depth == Catch::Approx(on.shadow_depth).margin(1e-6));
        }
}

TEST_CASE("a receiver millimetres up still spreads finitely in the shadow", "[shadow]")
{
    // 3 mm up, the dipping family covers a sliver of heights that deep, and its
    // tube used to come out 0 / 0 past the boundary.
    auto profile = read_atmosphere_profile("atmospheres/example_windy.csv", feet, 0.0);
    profile = profile.with_hydrostatic_pressure(reference_pressure);
    const LayeredSoundSpeed layered(profile, 90.0);
    const EffectiveSoundSpeed linear(profile, 90.0, 0.0, profile.top() + 20.0);
    constexpr double z_s = 100.0, z_r = 0.003;
    const Vec3 source(0.0, 0.0, z_s);
    std::vector<double> scratch;
    const double b = 270.0 * pi / 180.0;
    for (const double r : {1000.0, 1500.0, 2500.0})
    {
        const Vec3 receiver(r * std::sin(b), r * std::cos(b), z_r);
        const auto ray = layered_direct_ray(layered, linear, source, linear.moments(z_s), receiver,
                linear.moments(z_r), scratch);
        INFO(r << " m: spreading " << ray.path.range << ", shadow " << ray.shadow_depth);
        REQUIRE(ray.shadow_depth > 0.0);
        CHECK(std::isfinite(ray.path.range));
        CHECK(ray.path.range > 0.5 * r);
        CHECK(ray.path.range < 2.0 * r);
    }
}

TEST_CASE("a plate's field near the boundary is the ground's plus the plate over it", "[shadow]")
{
    // K, the GR1425 over bare Amedee ground near grazing: about nothing to
    // 1 kHz, a few dB at 4 kHz, most of pressure doubling's worth by 8 kHz.
    const std::vector<double> bands = {10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200,
            250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000,
            10000};
    GroundModel ground;
    ground.impedance_model = ImpedanceModel::variable_porosity;
    ground.flow_resistance = 200.0;
    const double c0 = 343.0;
    const auto plate = PlateReceiver::load("plates/gr1425_vp200.nc").bound(bands, ground, c0);
    const auto admittance = ground_impedance(ground, bands, c0);
    struct Expect
    {
        std::size_t band;
        double low, high;
    };
    for (const Expect e : {Expect {11, -1.0, 1.0}, Expect {20, -1.0, 1.0}, Expect {26, 1.0, 5.0},
                 Expect {29, 4.0, 11.0}})
    {
        const ShadowTransition transition(bands, admittance, plate.thickness(), c0, e.band, e.band + 1,
                &plate, &ground);
        const ShadowTransition point(bands, admittance, plate.thickness(), c0, e.band, e.band + 1);
        for (int i = 0; i < 12; i++)
        {
            const double azimuth = 2.0 * pi * i / 12.0;
            INFO(bands[e.band] << " Hz, azimuth " << 30 * i << ": K " << transition.plate_offset(e.band, azimuth));
            CHECK(transition.plate_offset(e.band, azimuth) > e.low);
            CHECK(transition.plate_offset(e.band, azimuth) < e.high);
            CHECK(point.plate_offset(e.band, azimuth) == 0.0);
        }
        // The bound covers V + K.
        CHECK(transition.max_power_ratio(e.band) >= point.max_power_ratio(e.band)
                        * std::pow(10.0, std::min(e.low, 0.0) / 10.0));
    }
}
