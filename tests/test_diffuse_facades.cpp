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

// The same off a roof (DiffuseOptions::roofs): a low block 20 m square, the
// source and the microphone together d above the middle of its roof.  Its
// walls are hidden from the microphone, and the floor reflects nothing a roof
// sends up, so the roof alone is gathered: rho s d^2 / pi times the integral
// of dA / r^6 over the roof.  Without roofs the field gathers nothing.
TEST_CASE("a point source lights a roof's diffuse field as Lambert's law says", "[diffuse]")
{
    const double d = 2.0, top = 6.0;
    const Buildings city({box(-10.0, -10.0, 10.0, 10.0, top)});
    const Vec3 point(0.0, 0.0, top + d);
    const std::vector<double> bands {1000.0};
    const std::vector<double> rho(city.facades().size(), 0.9), s(city.facades().size(), 0.5);
    DiffuseOptions options;
    options.radius = 50.0;
    options.patch_size = 0.4;
    double walls_only = 0.0;
    DiffuseFacades(city, point, rho, s, 1.0, bands, 343.2, options, 4).point_energy(point, 2, {}, &walls_only, 4);
    REQUIRE(walls_only == 0.0);
    options.roofs = true;
    const DiffuseFacades diffuse(city, point, rho, s, 1.0, bands, 343.2, options, 4);
    std::size_t roof = 0;
    double area = 0.0;
    for (const auto& p : diffuse.patches())
        if (p.roof)
        {
            roof++;
            area += p.area;
            REQUIRE(p.normal.z == 1.0);
            REQUIRE(p.center.z == Catch::Approx(top + 1e-3));
        }
    REQUIRE(roof == 50 * 50);
    REQUIRE(area == Catch::Approx(400.0));
    double energy = 0.0;
    diffuse.point_energy(point, 2, {}, &energy, 4);
    double integral = 0.0;
    const int n = 1000;
    const double step = 20.0 / n;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            const double x = -10.0 + (i + 0.5) * step, y = -10.0 + (j + 0.5) * step;
            integral += step * step / std::pow(d * d + x * x + y * y, 3);
        }
    REQUIRE(energy == Catch::Approx(0.9 * 0.5 * d * d / pi * integral).epsilon(5e-3));

    // An L-shaped roof is patched over the L alone, and a building whose
    // coefficient facade is out of reach keeps its roof whole.
    Footprint ell;
    ell.x = {0.0, 12.0, 12.0, 4.0, 4.0, 0.0};
    ell.y = {0.0, 0.0, 4.0, 4.0, 12.0, 12.0};
    ell.top = 5.0;
    const Buildings two({ell, box(500.0, 0.0, 510.0, 10.0, 5.0)});
    const std::vector<double> rho2(two.facades().size(), 0.9), s2(two.facades().size(), 0.5);
    options.patch_size = 1.0;
    const DiffuseFacades near(two, Vec3(-5.0, -5.0, 1.0), rho2, s2, 1.0, bands, 343.2, options, 2);
    double ell_area = 0.0;
    for (const auto& p : near.patches())
        if (p.roof)
        {
            REQUIRE(p.facade == two.roof_facade(0));
            ell_area += p.area;
        }
    REQUIRE(ell_area == Catch::Approx(80.0));
    for (std::size_t f = 0; f < two.facades().size(); f++)
        REQUIRE(near.scatters(f) == (two.facade_building(f) == 0));

    // With the roofs a building scatters as a whole: a long block whose far
    // wall is out of reach scatters there too, and only then.
    const Buildings long_block({box(0.0, 0.0, 200.0, 10.0, 5.0)});
    const std::vector<double> rho4(long_block.facades().size(), 0.9), s4(long_block.facades().size(), 0.5);
    const Vec3 west(-5.0, 5.0, 1.0);
    const DiffuseFacades whole(long_block, west, rho4, s4, 1.0, bands, 343.2, options, 2);
    DiffuseOptions walls = options;
    walls.roofs = false;
    const DiffuseFacades walls_only_field(long_block, west, rho4, s4, 1.0, bands, 343.2, walls, 2);
    std::size_t out_of_reach = 0;
    for (std::size_t f = 0; f < long_block.facades().size(); f++)
    {
        REQUIRE(whole.scatters(f));
        out_of_reach += !walls_only_field.scatters(f);
    }
    REQUIRE(out_of_reach == 1);

    // A roof sends nothing by way of the floor: a low block's roof and the
    // tower beside it exchange the same energy over a reflecting floor as over
    // none, where the walls facing each other across the floor do not.
    const Buildings pair({box(0.0, 0.0, 20.0, 20.0, 5.0), box(25.0, 0.0, 35.0, 20.0, 30.0)});
    const std::vector<double> rho3(pair.facades().size(), 0.9), s3(pair.facades().size(), 0.5);
    options.patch_size = 2.0;
    const DiffuseFacades bare(pair, Vec3(22.0, 10.0, 1.0), rho3, s3, 0.0, bands, 343.2, options, 2);
    const DiffuseFacades floored(pair, Vec3(22.0, 10.0, 1.0), rho3, s3, 1.0, bands, 343.2, options, 2);
    REQUIRE(bare.patches().size() == floored.patches().size());
    double from_roofs = 0.0, from_walls_bare = 0.0, from_walls_floored = 0.0;
    for (std::size_t p = 0; p < bare.patches().size(); p++)
    {
        if (bare.patches()[p].roof)
            continue;
        auto sum = [&](const DiffuseFacades& field, bool roofs) {
            double total = 0.0;
            for (const auto& [q, f] : field.form_factors(p))
                if (field.patches()[q].roof == roofs)
                    total += f;
            return total;
        };
        REQUIRE(sum(floored, true) == Catch::Approx(sum(bare, true)).epsilon(1e-12));
        from_roofs += sum(bare, true);
        from_walls_bare += sum(bare, false);
        from_walls_floored += sum(floored, false);
    }
    REQUIRE(from_roofs > 0.0);
    REQUIRE(from_walls_floored > from_walls_bare);
}

// A slope scatters as a level roof does, its patches tilted with it: a shed
// 20 m square whose roof rises 0.3 a metre along x, the source and the
// microphone together d along the slope's normal above its middle, no floor.
// In the plane the surface is a 20 sqrt(1.09) by 20 m rectangle centred under
// them, so the microphone gathers rho s d^2 / pi times the integral of
// dA / (d^2 + u^2 + v^2)^3 over it.  The slope's rays take rho (1 - s), so
// the share s must come back here, or a slope's reflection loses it.
TEST_CASE("a point source lights a slope's diffuse field as Lambert's law says", "[diffuse]")
{
    const double d = 2.0;
    Footprint shed = box(-10.0, -10.0, 10.0, 10.0, 20.0);
    shed.roof = {{0.3, 0.0, 6.0}};
    const Buildings city({shed});
    const Vec3 normal = shed.roof[0].normal();
    const Vec3 point = Vec3(0.0, 0.0, 6.0) + d * normal;
    const std::vector<double> bands {1000.0};
    const std::vector<double> rho(city.facades().size(), 0.9), s(city.facades().size(), 0.5);
    DiffuseOptions options;
    options.radius = 50.0;
    options.patch_size = 0.4;
    options.roofs = true;
    const DiffuseFacades diffuse(city, point, rho, s, 0.0, bands, 343.2, options, 4);
    double area = 0.0;
    std::size_t slope = 0;
    for (const auto& p : diffuse.patches())
        if (p.roof)
        {
            slope++;
            area += p.area;
            REQUIRE(p.normal.x == Catch::Approx(normal.x));
            REQUIRE(p.normal.z == Catch::Approx(normal.z));
            REQUIRE(shed.roof[0].above(p.center) == Catch::Approx(1e-3));
        }
    REQUIRE(slope == 50 * 50);
    const double width = 20.0 * std::sqrt(1.09);
    REQUIRE(area == Catch::Approx(width * 20.0));
    double energy = 0.0;
    diffuse.point_energy(point, 2, {}, &energy, 4);
    double integral = 0.0;
    const int n = 1000;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            const double u = -0.5 * width + (i + 0.5) * width / n, v = -10.0 + (j + 0.5) * 20.0 / n;
            integral += (width / n) * (20.0 / n) / std::pow(d * d + u * u + v * v, 3);
        }
    REQUIRE(energy == Catch::Approx(0.9 * 0.5 * d * d / pi * integral).epsilon(5e-3));
}

// A roof's patches cover it once: each sample of a cell lies on one surface,
// so a row of samples on a ridge, which is on both planes, counts for the
// first only.  A gable 4 m square whose ridge runs where a row of samples
// falls (1.5 m in, slopes of 1 either side): its patches' area is the roof's,
// 16 sqrt(2) m^2, where counting the ridge row on both slopes gave a quarter
// of a cell's samples twice.
TEST_CASE("a roof's patches cover a ridge once", "[diffuse]")
{
    Footprint gable = box(0.0, 0.0, 4.0, 4.0, 20.0);
    gable.roof = {{1.0, 0.0, 6.0}, {-1.0, 0.0, 9.0}};
    const Buildings city({gable});
    const std::vector<double> rho(city.facades().size(), 0.9), s(city.facades().size(), 0.5);
    DiffuseOptions options;
    options.radius = 50.0;
    options.patch_size = 4.0;
    options.roofs = true;
    const DiffuseFacades diffuse(city, Vec3(-5.0, 2.0, 12.0), rho, s, 0.0, {1000.0}, 343.2, options);
    double area = 0.0;
    for (const auto& p : diffuse.patches())
        if (p.roof)
            area += p.area;
    REQUIRE(area == Catch::Approx(16.0 * std::sqrt(2.0)).epsilon(1e-12));
}

// A wall lit off a roof (DiffuseOptions::roofs): a low block, flat or with a
// slope rising toward a tall one, and the source over the low roof.  Each of
// the tall block's wall patches over the roof is lit straight from the source
// and by the source's mirror in the roof's plane where the line from the
// mirror to the patch crosses the roof: cos(incidence) / L^2 each, the second
// times the roof's rho (1 - s).  Nothing else reaches the wall, and the floor
// reflects nothing.  The roof counts among the path's reflections: at order
// 0 the wall is lit straight alone.
TEST_CASE("a point source lights a wall off a roof, by its mirror in the roof", "[diffuse]")
{
    for (const bool slope : {false, true})
    {
        Footprint low = box(-100.0, -100.0, 0.0, 100.0, slope ? 30.0 : 10.0);
        if (slope)
            low.roof = {{0.1, 0.0, 15.0}};
        const Buildings city({low, box(0.0, -100.0, 10.0, 100.0, 60.0)});
        const Vec3 source(-50.0, 0.0, slope ? 40.0 : 30.0);
        const std::vector<double> bands {1000.0};
        const std::vector<double> rho(city.facades().size(), 0.8), s(city.facades().size(), 0.3);
        DiffuseOptions options;
        options.radius = 200.0;
        options.patch_size = 4.0;
        options.roofs = true;
        const DiffuseFacades diffuse(city, source, rho, s, 0.0, bands, 343.2, options, 4);
        // The source's mirror in the roof's plane, and that plane's height
        // above a point.
        const RoofPlane plane = slope ? low.roof.front() : RoofPlane {0.0, 0.0, 10.0};
        const Vec3 mirror = source - (2.0 * plane.above(source)) * plane.normal();
        std::vector<double> e, straight;
        diffuse.point_irradiance(source, 2, {}, e, 4);
        diffuse.point_irradiance(source, 0, {}, straight, 4);
        std::size_t checked = 0, off_roof = 0;
        for (std::size_t p = 0; p < diffuse.patches().size(); p++)
        {
            const auto& patch = diffuse.patches()[p];
            if (patch.roof || !(patch.normal.x < -0.99) || std::abs(patch.center.x) > 0.01 || e[p] == 0.0)
                continue;
            const Vec3 n = patch.normal;
            const Vec3 to_source = source - patch.center, to_mirror = mirror - patch.center;
            double expected = n.dot(to_source) / std::pow(to_source.norm(), 3);
            // Where the line to the mirror crosses the roof's plane: on the
            // low block's roof, and the patch above the plane.
            const double d = plane.above(patch.center), dm = -plane.above(mirror);
            const Vec3 bounce = mirror + (dm / (dm + d)) * (patch.center - mirror);
            if (d > 0.0 && bounce.x > -100.0 && bounce.x < 0.0 && std::abs(bounce.y) < 100.0)
            {
                expected += 0.8 * (1.0 - 0.3) * n.dot(to_mirror) / std::pow(to_mirror.norm(), 3);
                off_roof++;
            }
            INFO("slope " << slope << ", patch at " << patch.center.y << ", " << patch.center.z);
            CHECK(e[p] == Catch::Approx(expected).epsilon(1e-9));
            CHECK(straight[p] == Catch::Approx(n.dot(to_source) / std::pow(to_source.norm(), 3)).epsilon(1e-9));
            checked++;
        }
        CHECK(checked > 100);
        CHECK(off_roof > 100);
    }
}
