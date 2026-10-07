// --shadow_model resolved: each band's shadow boundary from the profile it
// resolves (resolved_height, tangent_below, resolved_shadow_bands and their
// table), and the cramped bands' residue series (ShadowTransition).

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "atmosphere_profile.hpp"
#include "exact_solutions.hpp"
#include "ground_effect.hpp"
#include "refraction.hpp"
#include "stratified.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT

namespace {

// A linear profile c0 + g z (g < 0 upward refraction) as layers.
struct LinearLayers
{
    std::vector<double> z, c;
    Layers layers() const { return Layers {z.data(), c.data(), z.size(), 0.0}; }
};

LinearLayers linear_layers(double c0, double g, double top, int levels)
{
    LinearLayers out;
    for (int k = 0; k < levels; k++)
    {
        const double z = top * k / (levels - 1);
        out.z.push_back(z);
        out.c.push_back(c0 + g * z);
    }
    return out;
}

// A file for read_atmosphere_profile: still air whose temperature falls
// linearly so that the effective sound speed's gradient is `gradient` /s at
// the ground (the frozen synthetic set's lin profiles, 2026-10-06).
std::string lapse_profile(double gradient)
{
    const auto path = std::filesystem::temp_directory_path()
            / ("niceops_lapse_" + std::to_string(static_cast<int>(-gradient * 1e5)) + ".csv");
    std::ofstream out(path);
    out << "z_ft,T_C,RH\n";
    const double t0 = 20.0;
    const double dtdz = gradient / (20.05 / (2.0 * std::sqrt(t0 + 273.15)));
    std::vector<double> levels = {0, 1, 3, 10, 20, 30, 40, 46};
    for (int z = 50; z <= 1000; z += 10)
        levels.push_back(z);
    char line[96];
    for (const double z : levels)
    {
        std::snprintf(line, sizeof line, "%g,%.5f,60.0\n", z, t0 + dtdz * z * 0.3048);
        out << line;
    }
    return path.string();
}

} // namespace

TEST_CASE("on a linear profile a band resolves it at the analytic height", "[resolved]")
{
    // k int_0^h sqrt(2 |a| z / c0) dz = phase gives
    // h = (3 phase / (2 k sqrt(2 |a| / c0)))^(2/3), 1.31 creeping-wave layers at a radian.
    constexpr double c0 = 340.0;
    for (const double g : {-0.01, -0.04})
        for (const double f : {125.0, 1000.0, 4000.0})
            for (const double phase : {0.5, 1.0, 2.0})
            {
                const auto profile = linear_layers(c0, g, 300.0, 31);
                const double k = 2.0 * pi * f / c0;
                const double expected = std::pow(3.0 * phase / (2.0 * k * std::sqrt(2.0 * -g / c0)), 2.0 / 3.0);
                INFO(g << " /s, " << f << " Hz, " << phase << " rad");
                CHECK(resolved_height(profile.layers(), f, phase) == Catch::Approx(expected).epsilon(1e-9));
                if (phase == 1.0)
                {
                    const double l = std::cbrt(c0 / (2.0 * -g * k * k));
                    CHECK(expected / l == Catch::Approx(std::cbrt(2.25)).epsilon(1e-12));
                }
            }
    // Downward refraction, or still air: nothing to resolve.
    const auto downward = linear_layers(c0, 0.02, 300.0, 31);
    CHECK(std::isinf(resolved_height(downward.layers(), 1000.0, 1.0)));
    // Faster than the ground only in a thin layer, then slower again: the band
    // gathers the layer's phase and no more.
    const std::vector<double> z = {0.0, 0.3, 1.0, 50.0};
    const std::vector<double> c = {340.0, 339.8, 340.5, 341.0};
    const Layers thin {z.data(), c.data(), z.size(), 0.0};
    CHECK(std::isinf(resolved_height(thin, 500.0, 1.0)));
    CHECK(resolved_height(thin, 50000.0, 1.0) < 1.0);
}

TEST_CASE("the tangent profile is the profile's line carried down below h", "[resolved]")
{
    // A log-law-like profile: steep near the ground, flatter above.
    const std::vector<double> z = {0.0, 0.3, 1.0, 3.0, 10.0, 30.0};
    const std::vector<double> c = {340.0, 339.6, 339.3, 339.0, 338.7, 338.5};
    const Layers layers {z.data(), c.data(), z.size(), 0.0};
    std::vector<double> tz, tc;
    double gradient = 0.0;
    const double h = 2.0; // in the layer from 1 to 3 m
    const Layers tangent = tangent_below(layers, h, tz, tc, gradient);
    const double slope = (339.0 - 339.3) / 2.0;
    CHECK(gradient == Catch::Approx(slope));
    // Below h, the line through (h, c(h)) at that slope; at and above h the
    // profile itself.
    for (const double height : {0.0, 0.2, 0.7, 1.5})
        CHECK(tangent.at(height) == Catch::Approx(layers.at(h) + slope * (height - h)));
    for (const double height : {2.0, 2.5, 3.0, 7.0, 20.0, 40.0})
        CHECK(tangent.at(height) == Catch::Approx(layers.at(height)));
    // So the ground is slower than the profile's own: the shadow boundary
    // moves out.
    CHECK(layered_shadow_boundary(tangent, 6.0, 0.003) > layered_shadow_boundary(layers, 6.0, 0.003));
    // Every h in one layer gives the same profile.
    std::vector<double> uz, uc;
    double g2 = 0.0;
    const Layers other = tangent_below(layers, 1.2, uz, uc, g2);
    CHECK(g2 == gradient);
    CHECK(layered_shadow_boundary(other, 6.0, 0.003)
            == Catch::Approx(layered_shadow_boundary(tangent, 6.0, 0.003)).epsilon(1e-9));
}

TEST_CASE("on a linear profile every band's resolved boundary is the profile's own", "[resolved]")
{
    constexpr double c0 = 340.0, g = -0.02, z_s = 20.0, z_r = 1.2;
    const auto profile = linear_layers(c0, g, 400.0, 81);
    const LinearProfile fit {c0, g};
    const std::vector<double> f = {125.0, 500.0, 2000.0};
    std::vector<ResolvedBand> bands(f.size());
    resolved_shadow_bands(profile.layers(), fit, z_s, z_r, f, 0, f.size(), 1.0, false, bands.data());
    const double r = shadow_boundary(fit, z_s, z_r);
    for (std::size_t b = 0; b < f.size(); b++)
    {
        INFO(f[b] << " Hz");
        CHECK(bands[b].layered_boundary == Catch::Approx(r).epsilon(2e-3));
        CHECK(bands[b].linear_boundary == Catch::Approx(r).epsilon(1e-9));
        CHECK(bands[b].gradient == Catch::Approx(-g).epsilon(1e-9));
        // r_s(f) and r_l(f) equal: band_boundary has nothing to interpolate.
        CHECK(bands[b].boundary == Catch::Approx(r).epsilon(2e-3));
    }
}

TEST_CASE("the resolved table gives a node's bands back", "[resolved]")
{
    auto profile = read_atmosphere_profile("atmospheres/example_windy.csv", feet, 0.0);
    profile = profile.with_hydrostatic_pressure(reference_pressure);
    const LayeredSoundSpeed layered(profile, 90.0);
    const EffectiveSoundSpeed linear(profile, 90.0, 0.0, 400.0);
    const std::vector<double> f = {125.0, 500.0, 2000.0};
    const ResolvedShadowTable table(layered, linear, f, 0, f.size(), 1.0, 0.0, 120.0, {1.2, 0.003},
            {false, true}, 2);
    REQUIRE_FALSE(table.empty());
    std::vector<double> scratch(layered.levels());
    for (const double degrees : {0.0, 90.0, 180.0, 270.0})
        for (const double z_s : {6.25, 25.0, 100.0}) // nodes: sqrt(z) on quarters
            for (std::size_t q = 0; q < 2; q++)
            {
                const double z_r = q == 0 ? 1.2 : 0.003;
                const double hx = std::cos(degrees * pi / 180.0), hy = std::sin(degrees * pi / 180.0);
                layered.along(hx, hy, scratch.data());
                const Layers layers {layered.heights().data(), scratch.data(), layered.levels(),
                        layered.resolution_frequency()};
                std::vector<ResolvedBand> direct(f.size()), looked(f.size());
                resolved_shadow_bands(layers, linear.fit(z_s, linear.moments(z_s), z_r, linear.moments(z_r), hx, hy),
                        z_s, z_r, f, 0, f.size(), 1.0, q == 1, direct.data());
                table.lookup(hx, hy, z_s, q, looked.data());
                for (std::size_t b = 0; b < f.size(); b++)
                {
                    INFO(degrees << " deg, z_s " << z_s << ", z_r " << z_r << ", " << f[b] << " Hz");
                    if (std::isinf(direct[b].boundary))
                        CHECK(std::isinf(looked[b].boundary));
                    else
                        CHECK(looked[b].boundary == Catch::Approx(direct[b].boundary).epsilon(1e-6));
                }
            }
}

TEST_CASE("the cramped table is the residue series at the source's height", "[resolved]")
{
    constexpr double c0 = 343.0, z_r = 0.003;
    const std::vector<double> f = {500.0};
    GroundModel ground;
    ground.impedance_model = ImpedanceModel::variable_porosity;
    ground.flow_resistance = 200.0;
    ground.receiver_height = z_r;
    const auto admittance = ground_impedance(ground, f, c0);
    ShadowTransition transition(f, admittance, z_r, c0, 0, 1, nullptr, nullptr, 2);
    transition.set_cramped_ground(CrampedGround::series, &ground, 2);
    // On a node in s and sqrt(y_s), at xi nodes: the series itself.
    const double s = ShadowTransition::s_min * std::pow(10.0, 8.0 / ShadowTransition::s_per_decade);
    const double k0 = 2.0 * pi * f[0] / c0;
    const double a = pi * f[0] / (s * s * s);
    const double l = s / k0, diffraction = 2.0 * s * s / k0;
    const double y_s = 0.8 * 0.8; // sqrt(y_s) = 0.8, a node
    const double z_s = y_s * l;
    const double boundary = shadow_boundary({c0, -a}, z_s, z_r);
    ExactProblem problem;
    problem.frequency = f[0];
    problem.sound_speed = c0;
    problem.admittance = admittance[0];
    problem.source_height = z_s;
    problem.receiver_height = z_r;
    for (const double xi : {-0.4, 0.0, 1.0, 3.0})
    {
        const double r = boundary + xi * diffraction;
        const auto field = upward_pressure(problem, c0 / a, {r});
        const double expected = 20.0 * std::log10(std::abs(field[0]) * std::hypot(r, z_s - z_r));
        INFO("xi " << xi);
        CHECK(transition.cramped_level(0, xi, s, y_s) == Catch::Approx(expected).margin(1e-3));
    }
    // Past the cramped heights there is no value.
    CHECK(std::isnan(transition.cramped_level(0, 0.0, s, 5.0)));
}

TEST_CASE("a cramped source over a linear lapse follows the PE with --shadow_model resolved", "[resolved]")
{
    // The frozen synthetic set's lin020 case (2026-10-06): still air, -0.02 /s,
    // source 20 ft, flush receiver (3 mm), variable porosity 200, against
    // niceops_pe_compare's PE (split-step Pade, self-starter, five frequencies
    // a band) from 800 to 2000 ft.  The boundary is 1534 ft, within two
    // diffraction lengths of the source at 500 Hz and 1 kHz, so the legacy
    // transition keeps the ray model to it -- up to 17 dB under the PE as the
    // ray's equivalent ground geometry collapses -- where the residue series
    // at the source's own height is within 0.3 dB.
    auto profile = read_atmosphere_profile(lapse_profile(-0.02), feet, 0.0);
    profile = profile.with_hydrostatic_pressure(reference_pressure);
    const double c0 = profile.at(0.0).sound_speed();
    const double z_s = 20.0 * 0.3048, z_r = 0.01 * 0.3048;
    const LayeredSoundSpeed layered(profile, 90.0);
    const EffectiveSoundSpeed linear(profile, 90.0, 0.0,
            std::max(profile.top(), z_s) + 2.0 * EffectiveSoundSpeed::min_fit_span);
    struct Band
    {
        double frequency;
        std::vector<double> pe;
    };
    const std::vector<Band> bands = {
        {500.0, {-7.58, -8.79, -9.96, -11.08, -12.18, -13.25, -14.30, -15.34, -16.37, -17.39, -18.41, -19.42, -20.43}},
        {1000.0, {-10.62, -12.00, -13.35, -14.67, -15.98, -17.28, -18.57, -19.86, -21.16, -22.46, -23.76, -25.08, -26.40}},
    };
    GroundModel ground;
    ground.impedance_model = ImpedanceModel::variable_porosity;
    ground.flow_resistance = 200.0;
    ground.receiver_height = z_r;
    const double hx = 1.0, hy = 0.0;
    std::vector<double> c_along(layered.levels());
    layered.along(hx, hy, c_along.data());
    const Layers path {layered.heights().data(), c_along.data(), layered.levels(), layered.resolution_frequency()};
    for (const Band& band : bands)
    {
        const std::vector<double> f = {band.frequency};
        const auto admittance = ground_impedance(ground, f, c0);
        const ShadowAttenuation attenuation(f, admittance, 0, 1);
        ShadowTransition transition(f, admittance, z_r, c0, 0, 1, nullptr, nullptr, 2);
        transition.set_cramped_ground(CrampedGround::series, &ground, 2);
        std::vector<ResolvedBand> resolved(1);
        resolved_shadow_bands(path, linear.fit(z_s, linear.moments(z_s), z_r, linear.moments(z_r), hx, hy), z_s,
                z_r, f, 0, 1, 1.0, true, resolved.data());
        const Vec3 source(0.0, 0.0, z_s);
        std::vector<double> scratch;
        double worst = 0.0, worst_legacy = 0.0;
        for (std::size_t i = 0; i < band.pe.size(); i++)
        {
            const double r = (800.0 + 100.0 * static_cast<double>(i)) * 0.3048;
            const Vec3 receiver(r * hx, r * hy, z_r);
            RefractedRay ray = layered_direct_ray(layered, linear, source, linear.moments(z_s), receiver,
                    linear.moments(z_r), scratch, true);
            layered_finish_ground(layered, ray, GroundGeometry::point, z_r, scratch);
            const double ground_db = excess_ground_attenuation(ground, ray.path.source_height,
                    ray.path.ground_distance, c0, f)[0];
            const double spreading = 20.0 * std::log10(std::hypot(r, z_s - z_r) / ray.path.range);
            double weight = 1.0, loss = 0.0;
            REQUIRE(transition.apply_resolved(ray, resolved.data(), attenuation, 20.0, 0, 1, &weight, &loss));
            const double level = spreading + (weight > 0.0 ? weight * ground_db : 0.0) - loss;
            double legacy_weight = 1.0, legacy_loss = 0.0;
            if (!transition.apply(ray, attenuation, 20.0, 0, 1, &legacy_weight, &legacy_loss))
                legacy_loss = 0.0;
            const double legacy = spreading + (legacy_weight > 0.0 ? legacy_weight * ground_db : 0.0) - legacy_loss;
            INFO(band.frequency << " Hz, " << r / 0.3048 << " ft: PE " << band.pe[i] << ", resolved " << level
                                << ", legacy " << legacy);
            CHECK(std::fabs(level - band.pe[i]) < 0.5);
            worst = std::max(worst, std::fabs(level - band.pe[i]));
            worst_legacy = std::max(worst_legacy, std::fabs(legacy - band.pe[i]));
        }
        INFO(band.frequency << " Hz: worst " << worst << " dB, legacy " << worst_legacy << " dB");
        CHECK(worst < 0.5);
        CHECK(worst_legacy > 5.0); // the case is the defect, not a lit path
    }
}

TEST_CASE("resolved levels stay inside the peak search's bound", "[resolved][bounds]")
{
    // Whatever a band takes -- the ray, the straight path's ground term, the
    // table or the cramped series -- its gain over the ray's spreading is at most
    // the transition's largest power ratio or the ground's own bound.
    constexpr double c0 = 343.0, z_r = 1.2;
    const std::vector<double> f = {250.0, 1000.0};
    GroundModel ground;
    ground.receiver_height = z_r;
    const auto admittance = ground_impedance(ground, f, c0);
    const auto ground_bound = ground_gain_bounds(ground, admittance, 0, f.size());
    const ShadowAttenuation attenuation(f, admittance, 0, f.size());
    for (const CrampedGround mode : {CrampedGround::series, CrampedGround::straight, CrampedGround::ray})
    {
        ShadowTransition transition(f, admittance, z_r, c0, 0, f.size(), nullptr, nullptr, 2);
        transition.set_cramped_ground(mode, &ground, 2);
        for (const double g : {-0.01, -0.05})
            for (const double z_s : {2.0, 6.0, 30.0})
                for (double r = 20.0; r < 3000.0; r *= 1.25)
                {
                    const LinearProfile profile {c0, g};
                    RefractedRay ray;
                    ray.curved = true;
                    ray.profile = profile;
                    ray.horizontal = r;
                    ray.source_z = z_s;
                    ray.receiver_z = z_r;
                    ray.direct = direct_eigenray(profile, z_s, z_r, r);
                    ray.shadow_depth = ray.direct.shadow_depth;
                    ray.path.range = ray.direct.spreading_range;
                    ray.path.source_height = z_s;
                    ray.path.ground_distance = r;
                    std::vector<ResolvedBand> bands(f.size());
                    for (std::size_t b = 0; b < f.size(); b++)
                        bands[b].boundary = ray.direct.boundary * (b == 0 ? 1.3 : 1.0);
                    std::vector<double> weight(f.size()), loss(f.size());
                    REQUIRE(transition.apply_resolved(ray, bands.data(), attenuation, 20.0, 0, f.size(),
                            weight.data(), loss.data()));
                    for (std::size_t b = 0; b < f.size(); b++)
                    {
                        // Ground term at most the ground's bound, weighted; loss at
                        // least minus the larger of the two.
                        const double cap = 10.0
                                * std::log10(std::max(ground_bound[b], transition.max_power_ratio(b)));
                        INFO("mode " << static_cast<int>(mode) << ", g " << g << ", z_s " << z_s << ", r " << r
                                     << ", band " << b << ": weight " << weight[b] << ", loss " << loss[b]);
                        CHECK(weight[b] * 10.0 * std::log10(ground_bound[b]) - loss[b] <= cap + 1e-9);
                    }
                }
    }
}
