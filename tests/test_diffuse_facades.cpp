#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "building_paths.hpp"
#include "buildings.hpp"
#include "diffuse_facades.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT

namespace {

Footprint box(double x0, double y0, double x1, double y1, double top)
{
    Footprint f;
    f.x = {x0, x1, x1, x0};
    f.y = {y0, y0, y1, y1};
    f.top = top;
    return f;
}

Vec3 direction(double azimuth_degrees, double elevation_degrees)
{
    const double a = azimuth_degrees * pi / 180.0, e = elevation_degrees * pi / 180.0;
    return {std::cos(e) * std::cos(a), std::cos(e) * std::sin(a), std::sin(e)};
}

// The solid angle of a rectangle seen from a point on the normal through one
// of its corners, at distance d, its sides a and b.
double corner_solid_angle(double a, double b, double d)
{
    return std::atan(a * b / (d * std::sqrt(d * d + a * a + b * b)));
}

} // namespace

TEST_CASE("a lit wall's diffuse field is rho s cos(incidence) times its solid angle over pi", "[diffuse]")
{
    // The face x = 0 of a building 4 km wide and 2 km tall, the microphone
    // 200 m in front of its middle; no floor reflection, only that face in
    // reach.
    const Buildings city({box(0.0, -2000.0, 10.0, 2000.0, 2000.0)});
    const Vec3 microphone(-200.0, 0.0, 1000.0);
    const std::vector<double> bands {500.0};
    const std::vector<double> reflection(city.facades().size(), 0.9), scattering(city.facades().size(), 0.4);
    DiffuseOptions options;
    options.radius = 205.0;
    options.patch_size = 20.0;
    const DiffuseFacades diffuse(city, microphone, reflection, scattering, 0.0, bands, 343.2, options, 4);
    REQUIRE(diffuse.patches().size() == 200 * 100);
    REQUIRE(diffuse.converged());
    const double omega = 4.0 * corner_solid_angle(2000.0, 1000.0, 200.0);
    for (double azimuth : {180.0, 150.0, 215.0})
        for (double elevation : {5.0, 30.0, 60.0})
        {
            const Vec3 u = direction(azimuth, elevation);
            double e = 0.0;
            diffuse.energy(u, &e);
            const double incidence = -u.x;
            INFO("azimuth " << azimuth << " elevation " << elevation);
            REQUIRE(e == Catch::Approx(0.9 * 0.4 * incidence * omega / pi).epsilon(2e-3));
        }
    // From behind the face, nothing.
    double e = 1.0;
    diffuse.energy(direction(0.0, 20.0), &e);
    REQUIRE(e == 0.0);
}

TEST_CASE("the reverse solve is the forward one", "[diffuse]")
{
    // A street between two buildings, its walls exchanging energy through the
    // air and off the floor.
    const Buildings city({box(-40.0, -20.0, 40.0, 0.0, 22.0), box(-35.0, 15.0, 35.0, 35.0, 30.0)});
    const Vec3 microphone(3.0, 6.0, 1.5);
    const std::vector<double> bands {125.0, 1000.0, 4000.0};
    std::vector<double> reflection(city.facades().size() * bands.size()), scattering(reflection.size());
    for (std::size_t i = 0; i < reflection.size(); i++)
    {
        reflection[i] = 0.7 + 0.05 * static_cast<double>(i % 5);
        scattering[i] = 0.1 + 0.1 * static_cast<double>(i % 3);
    }
    DiffuseOptions options;
    options.patch_size = 5.0;
    const DiffuseFacades diffuse(city, microphone, reflection, scattering, 0.9, bands, 343.2, options, 4);
    REQUIRE(diffuse.converged());
    std::vector<double> reverse(bands.size()), forward(bands.size());
    for (double azimuth : {80.0, 100.0, 260.0, 300.0})
        for (double elevation : {10.0, 40.0, 75.0})
        {
            const Vec3 u = direction(azimuth, elevation);
            diffuse.energy(u, reverse.data());
            diffuse.energy_forward(u, forward.data());
            for (std::size_t b = 0; b < bands.size(); b++)
            {
                INFO("azimuth " << azimuth << " elevation " << elevation << " band " << bands[b]);
                REQUIRE(reverse[b] > 0.0);
                REQUIRE(reverse[b] == Catch::Approx(forward[b]).epsilon(1e-8));
            }
        }
}

TEST_CASE("form factors are reciprocal and bounded", "[diffuse]")
{
    const Buildings city({box(-30.0, -20.0, 30.0, 0.0, 20.0), box(-30.0, 12.0, 30.0, 32.0, 20.0)});
    DiffuseOptions options;
    options.patch_size = 4.0;
    const std::vector<double> bands {1000.0};
    const std::vector<double> reflection(city.facades().size(), 0.8), scattering(city.facades().size(), 0.5);
    const DiffuseFacades diffuse(city, Vec3(0.0, 6.0, 1.5), reflection, scattering, 1.0, bands, 343.2, options, 4);
    const auto& patches = diffuse.patches();
    std::size_t pairs = 0;
    for (std::size_t p = 0; p < patches.size(); p++)
    {
        double sum = 0.0;
        for (const auto& [q, f] : diffuse.form_factors(p))
        {
            sum += f;
            // A_p F_pq = A_q F_qp, the point form being exact both ways
            // where neither is subdivided.
            const auto& back = diffuse.form_factors(q);
            const auto it = std::find_if(back.begin(), back.end(), [&](const auto& e) { return e.first == p; });
            REQUIRE(it != back.end());
            const double distance = (patches[p].center - patches[q].center).norm();
            if (distance > 3.0 * 4.0 + 1.0)
            {
                REQUIRE(patches[p].area * f == Catch::Approx(patches[q].area * it->second).epsilon(1e-12));
                pairs++;
            }
        }
        // What leaves a patch for the opposite wall and its floor image is at
        // most what a hemisphere and its image could receive.
        REQUIRE(sum <= 2.0);
    }
    REQUIRE(pairs > 1000);
}

TEST_CASE("the diffuse energy adds to the table's band energy", "[diffuse]")
{
    const Buildings city({box(-40.0, -20.0, 40.0, 0.0, 22.0), box(-35.0, 15.0, 35.0, 35.0, 30.0)});
    const Vec3 microphone(3.0, 6.0, 1.5);
    const std::vector<double> bands {250.0, 2000.0};
    const std::size_t nf = city.facades().size();
    const std::vector<double> rho(nf * bands.size(), 0.9), s(nf * bands.size(), 0.3);
    std::vector<double> specular(nf * bands.size());
    for (std::size_t i = 0; i < specular.size(); i++)
        specular[i] = rho[i] * (1.0 - s[i]);
    const DiffuseFacades diffuse(city, microphone, rho, s, 1.0, bands, 343.2, DiffuseOptions {}, 2);
    for (bool diffraction : {false, true})
    {
        BuildingPathOptions options;
        options.diffraction = diffraction;
        const BuildingReceiver receiver(city, microphone, facing_facades(city, microphone, 100.0), options);
        const Vec3 u = direction(95.0, 25.0);
        std::vector<double> with(bands.size()), without(bands.size()), scattered(bands.size()), open(bands.size());
        receiver.ratio(u, specular, 1.0, bands, 343.2, with.data(), &diffuse);
        receiver.ratio(u, specular, 1.0, bands, 343.2, without.data());
        diffuse.energy(u, scattered.data());
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            REQUIRE(scattered[b] > 0.0);
            // ratio = (rays + diffuse) / open, so the difference is the
            // diffuse energy over open ground's: 1 + rho_f in energy, or its
            // two rays band-averaged over the three sub-bands.
            double open = 2.0;
            if (diffraction)
            {
                const double half = std::pow(2.0, 1.0 / 6.0), f1 = bands[b] / half, f2 = bands[b] * half;
                double sum = 0.0;
                for (int part = 0; part < 3; part++)
                {
                    const double g1 = f1 * std::pow(f2 / f1, part / 3.0), g2 = f1 * std::pow(f2 / f1, (part + 1) / 3.0);
                    sum += band_averaged_energy({1.0, 1.0}, {0.0, 2.0 * microphone.z * u.z}, 2.0 * pi * g1 / 343.2,
                                   2.0 * pi * g2 / 343.2)
                            * (g2 - g1);
                }
                open = sum / (f2 - f1);
            }
            REQUIRE(with[b] - without[b] == Catch::Approx(scattered[b] / open).epsilon(1e-9));
        }
    }
}

TEST_CASE("every specular reflection lights the facades, carrying rho (1 - s)", "[diffuse]")
{
    // Two walls 20 m apart and 30 m tall, a wave across them at 20 degrees: it
    // climbs 7.279 m a crossing, and a patch is lit through as many
    // reflections as it takes to climb out, along the wave only if that
    // ends on the right wall.  No floor.
    const Buildings city({box(-20.0, -20.0, 20.0, 0.0, 30.0), box(-20.0, 20.0, 20.0, 40.0, 30.0)});
    const std::vector<double> bands {1000.0};
    std::vector<double> reflection(city.facades().size(), 0.9), scattering(city.facades().size(), 0.0);
    std::size_t wall_a = 0, wall_b = 0;
    for (std::size_t f = 0; f < city.facades().size(); f++)
    {
        const Facade& facade = city.facades()[f];
        if (facade.y0 == 0.0 && facade.y1 == 0.0)
        {
            wall_a = f;
            reflection[f] = 0.8;
            scattering[f] = 0.25;
        }
        if (facade.y0 == 20.0 && facade.y1 == 20.0)
        {
            wall_b = f;
            reflection[f] = 0.9;
            scattering[f] = 0.5;
        }
    }
    const double a = 0.8 * 0.75, b = 0.9 * 0.5; // each wall's specular share
    DiffuseOptions options;
    options.patch_size = 2.0;
    const DiffuseFacades diffuse(city, Vec3(0.0, 10.0, 1.5), reflection, scattering, 0.0, bands, 343.2, options, 4);
    const double elevation = 20.0 * pi / 180.0;
    const Vec3 u(0.0, std::cos(elevation), std::sin(elevation));
    std::vector<double> e;
    diffuse.irradiance(u, e);
    std::size_t checked = 0;
    for (std::size_t p = 0; p < diffuse.patches().size(); p++)
    {
        const auto& patch = diffuse.patches()[p];
        if (patch.facade != wall_a && patch.facade != wall_b)
            continue;
        // Back along the wave from wall A, against it from wall B, bouncing
        // between them: the crossing is 20 m less the patches' millimetre.
        const bool on_a = patch.facade == wall_a;
        const double climb = (20.0 - 2e-3) * std::tan(elevation);
        double z = patch.center.z, carried = std::cos(elevation);
        int hits = 0;
        while (z + climb <= 30.0)
        {
            z += climb;
            hits++;
            const bool at_b = (hits % 2 == 1) == on_a;
            carried *= at_b ? b : a;
        }
        // From A the ray leaves along the wave after an even number of
        // reflections, from B after an odd one.
        const double expected = (hits % 2 == 0) == on_a ? carried : 0.0;
        INFO("wall " << (on_a ? "A" : "B") << " at " << patch.center.z << " m, " << hits << " reflections");
        REQUIRE(e[p] == Catch::Approx(expected).margin(1e-12));
        checked++;
    }
    REQUIRE(checked == 2 * 20 * 15);
}

TEST_CASE("a point source lights a wall's diffuse field as Lambert's law says", "[diffuse]")
{
    // A wall 10 m square, the source and the microphone together 1 m in front
    // of its middle, nothing else lit and no floor.  The microphone gathers
    //   rho s d^2 / pi  integral over the wall of dA / r^6,
    // which over an infinite wall is rho s / (2 d^2); here it is taken over the
    // square itself, by the midpoint rule.
    const double d = 1.0;
    const Buildings city({box(0.0, -5.0, 10.0, 5.0, 10.0)});
    const Vec3 point(-d, 0.0, 5.0);
    const std::vector<double> bands {1000.0};
    const std::vector<double> rho(city.facades().size(), 0.9), s(city.facades().size(), 0.5);
    DiffuseOptions options;
    options.radius = 50.0;
    options.patch_size = 0.2;
    const DiffuseFacades diffuse(city, point, rho, s, 0.0, bands, 343.2, options, 4);
    double energy = 0.0;
    diffuse.point_energy(point, 4, {}, &energy, 4);
    double integral = 0.0;
    const int n = 1000;
    const double step = 10.0 / n;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            const double y = -5.0 + (i + 0.5) * step, z = -5.0 + (j + 0.5) * step;
            integral += step * step / std::pow(d * d + y * y + z * z, 3);
        }
    const double exact = 0.9 * 0.5 * d * d / pi * integral;
    REQUIRE(exact == Catch::Approx(0.9 * 0.5 / (2.0 * d * d)).epsilon(2e-3));
    REQUIRE(energy == Catch::Approx(exact).epsilon(5e-3));
    // The air's absorption along the lighting paths only lowers it.
    double absorbed = 0.0;
    diffuse.point_energy(point, 4, {0.05}, &absorbed, 4);
    REQUIRE(absorbed < energy);
    REQUIRE(absorbed > energy * std::exp(-2.0 * 0.05 * 8.0));
}
