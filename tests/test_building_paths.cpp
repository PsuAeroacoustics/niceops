#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include "building_paths.hpp"
#include "buildings.hpp"
#include "urban.hpp"
#include "util/constants.hpp"
#include "wedge.hpp"

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

// An ideal canyon built from two long buildings along x: wall A the face of
// the first at y = 0, wall B the second's at y = width.
Buildings canyon_buildings(const IdealCanyon& c)
{
    const double long_way = 20000.0;
    return Buildings({box(-long_way, -20.0, long_way, 0.0, c.height_a),
            box(-long_way, c.width, long_way, c.width + 20.0, c.height_b)});
}

// The field at the microphone at wavenumber k, every facade rigid and the
// floor absorbing, from the rays `keep` accepts.
template <typename Keep>
std::complex<double> field(const std::vector<BuildingRay>& rays, std::size_t facades, double k, Keep keep)
{
    const std::vector<double> rigid(facades, 1.0);
    std::complex<double> sum = 0.0;
    for (const auto& ray : rays)
        if (keep(ray))
            sum += ray_amplitude(ray, rigid, 0.0, k) * std::polar(1.0, k * ray.length);
    return sum;
}

} // namespace

TEST_CASE("open ground is the ratio's unit", "[building_paths]")
{
    const Buildings none;
    BuildingReceiver receiver(none, Vec3(3.0, -4.0, 1.5), {}, BuildingPathOptions {});
    const std::vector<double> bands {125.0, 500.0, 2000.0};
    std::vector<double> out(bands.size());
    for (double elevation : {0.5, 10.0, 45.0, 90.0})
    {
        receiver.ratio(direction(30.0, elevation), {}, 0.8, bands, 343.2, out.data());
        for (double v : out)
            REQUIRE(v == Catch::Approx(1.0).epsilon(1e-12));
    }
}

TEST_CASE("long buildings reproduce the ideal canyon's specular rays exactly", "[building_paths]")
{
    const std::vector<IdealCanyon> canyons {{30.0, 10.0, 25.0, 40.0, 1.5}, {12.0, 6.0, 30.0, 30.0, 4.0},
            {50.0, 45.0, 15.0, 60.0, 0.0}};
    const std::vector<double> bands {250.0, 1000.0};
    for (const auto& c : canyons)
    {
        const Buildings city = canyon_buildings(c);
        const Vec3 microphone(0.0, c.offset, c.microphone_height);
        BuildingPathOptions options;
        options.diffraction = false;
        options.max_order = 400;
        const auto run = facing_facades(city, microphone, options.radius);
        REQUIRE(run.size() == 2);
        // Wall A (the face at y = 0) reflects 0.9 and 0.6, wall B 0.8 and 0.3.
        std::vector<double> reflection(city.facades().size() * bands.size(), 1.0);
        const double ra[] {0.9, 0.6}, rb[] {0.8, 0.3};
        for (std::size_t f = 0; f < city.facades().size(); f++)
        {
            const Facade& facade = city.facades()[f];
            if (facade.y0 == 0.0 && facade.y1 == 0.0)
                for (std::size_t b = 0; b < 2; b++)
                    reflection[f * 2 + b] = ra[b];
            if (facade.y0 == c.width && facade.y1 == c.width)
                for (std::size_t b = 0; b < 2; b++)
                    reflection[f * 2 + b] = rb[b];
        }
        const BuildingReceiver receiver(city, microphone, run, options);
        std::vector<double> general(bands.size()), ideal(bands.size());
        for (double azimuth = 0.0; azimuth < 360.0; azimuth += 7.5)
            for (double elevation : {1.0, 3.0, 8.0, 20.0, 45.0, 75.0})
            {
                const Vec3 u = direction(azimuth, elevation);
                if (std::abs(u.y) < 1e-6)
                    continue;
                receiver.ratio(u, reflection, 0.7, bands, 343.2, general.data());
                canyon_ratio(c, u.y, u.z, ra, rb, 0.7, bands, 343.2, false, ideal.data());
                for (std::size_t b = 0; b < bands.size(); b++)
                {
                    INFO("canyon " << c.width << " azimuth " << azimuth << " elevation " << elevation);
                    REQUIRE(general[b] == Catch::Approx(ideal[b]).epsilon(1e-12).margin(1e-28));
                }
            }
    }
}

TEST_CASE("long buildings reproduce the ideal canyon's coherent field", "[building_paths]")
{
    const std::vector<IdealCanyon> canyons {{30.0, 10.0, 25.0, 40.0, 1.5}, {16.0, 8.0, 30.0, 30.0, 0.0}};
    const std::vector<double> bands {125.0, 500.0, 1600.0};
    for (const auto& c : canyons)
    {
        const Buildings city = canyon_buildings(c);
        const Vec3 microphone(0.0, c.offset, c.microphone_height);
        BuildingPathOptions options;
        options.max_order = 400;
        options.max_diffracted_order = 40;
        const auto run = facing_facades(city, microphone, options.radius);
        std::vector<double> reflection(city.facades().size() * bands.size(), 0.9);
        const std::vector<double> ra(bands.size(), 0.9), rb(bands.size(), 0.9);
        const BuildingReceiver receiver(city, microphone, run, options);
        // The two inner faces, each mirrored up to 40 times in turn.
        REQUIRE(receiver.images().size() == 81);
        std::vector<double> general(bands.size()), ideal(bands.size());
        // Away from the axis, where the diffraction point lies within the
        // 40 km of edge.
        for (double azimuth : {20.0, 45.0, 80.0, 100.0, 150.0, 200.0, 260.0, 300.0, 335.0})
            for (double elevation : {2.0, 6.0, 15.0, 30.0, 60.0, 85.0})
            {
                const Vec3 u = direction(azimuth, elevation);
                receiver.ratio(u, reflection, 0.8, bands, 343.2, general.data());
                canyon_ratio(c, u.y, u.z, ra.data(), rb.data(), 0.8, bands, 343.2, true, ideal.data());
                for (std::size_t b = 0; b < bands.size(); b++)
                {
                    INFO("canyon " << c.width << " azimuth " << azimuth << " elevation " << elevation << " band "
                                   << bands[b]);
                    REQUIRE(10.0 * std::log10(general[b]) == Catch::Approx(10.0 * std::log10(ideal[b])).margin(1e-9));
                }
            }
    }
}

TEST_CASE("an oblique wave at a long roof edge is the 2-D wedge at k sin(beta)", "[building_paths]")
{
    // A roof edge a kilometre up, so that the floor and the wall's foot are
    // out of reach: the edge along the face x = 0 of a building filling x > 0.
    const double top = 1000.0;
    const Buildings city({box(0.0, -50000.0, 100000.0, 50000.0, top)});
    std::size_t edge = city.edges().size();
    for (std::size_t e = 0; e < city.edges().size(); e++)
        if (!city.edges()[e].vertical && city.edges()[e].start.x == 0.0)
            edge = e;
    REQUIRE(edge < city.edges().size());
    BuildingPathOptions options;
    const double k = 2.0 * pi * 500.0 / 343.2;
    for (const Vec3 offset : {Vec3(-3.0, 0.0, -2.0), Vec3(-0.5, 0.0, -6.0), Vec3(-8.0, 0.0, -0.3)})
    {
        const Vec3 microphone = Vec3(0.0, 0.0, top) + offset;
        const BuildingReceiver receiver(city, microphone, {}, options);
        const double r = std::hypot(offset.x, offset.z);
        const double theta = city.edges()[edge].angle(offset);
        for (double azimuth : {0.0, 30.0, 100.0, 150.0, 180.0, 210.0, 300.0})
            for (double elevation : {3.0, 15.0, 40.0, 70.0})
            {
                const Vec3 u = direction(azimuth, elevation);
                const auto rays = receiver.rays(u);
                const auto general = field(rays, city.facades().size(), k,
                        [&](const BuildingRay& ray) { return !ray.diffracted || ray.edge == edge; });
                const double sin_beta = std::sqrt(1.0 - u.y * u.y);
                const double theta0 = city.edges()[edge].angle(Vec3(u.x, 0.0, u.z));
                const auto exact = hard_wedge_field(1.5, k * sin_beta * r, theta, theta0);
                INFO("offset " << offset << " azimuth " << azimuth << " elevation " << elevation);
                REQUIRE(std::abs(general) == Catch::Approx(std::abs(exact)).epsilon(1e-9));
            }
    }
}

TEST_CASE("an oblique wave at a vertical corner is the 2-D wedge at k sin(beta)", "[building_paths]")
{
    // The corner at the origin of a building filling x > 0, y > 0, the
    // microphone half way up it.
    const Buildings city({box(0.0, 0.0, 100000.0, 100000.0, 2000.0)});
    std::size_t edge = city.edges().size();
    for (std::size_t e = 0; e < city.edges().size(); e++)
        if (city.edges()[e].vertical && city.edges()[e].start.x == 0.0 && city.edges()[e].start.y == 0.0)
            edge = e;
    REQUIRE(edge < city.edges().size());
    REQUIRE(city.edges()[edge].n == Catch::Approx(1.5));
    BuildingPathOptions options;
    const double k = 2.0 * pi * 1000.0 / 343.2;
    for (const Vec3 offset : {Vec3(-3.0, -2.0, 0.0), Vec3(-0.4, -7.0, 0.0), Vec3(-6.0, 0.5, 0.0)})
    {
        const Vec3 microphone = Vec3(0.0, 0.0, 1000.0) + offset;
        const BuildingReceiver receiver(city, microphone, {}, options);
        const double r = std::hypot(offset.x, offset.y);
        const double theta = city.edges()[edge].angle(offset);
        for (double azimuth : {190.0, 200.0, 225.0, 260.0, 300.0, 120.0, 30.0})
            for (double elevation : {0.1, 10.0, 35.0, 60.0})
            {
                const Vec3 u = direction(azimuth, elevation);
                const auto rays = receiver.rays(u);
                const auto general = field(rays, city.facades().size(), k,
                        [&](const BuildingRay& ray) { return !ray.diffracted || ray.edge == edge; });
                const Vec3 incident(u.x, u.y, 0.0);
                const double theta0 = city.edges()[edge].angle(incident);
                if (!(theta0 > 0.0 && theta0 < 1.5 * pi))
                    continue; // from inside the building
                const double sin_beta = std::sqrt(1.0 - u.z * u.z);
                const auto exact = hard_wedge_field(1.5, k * sin_beta * r, theta, theta0);
                INFO("offset " << offset << " azimuth " << azimuth << " elevation " << elevation);
                REQUIRE(std::abs(general) == Catch::Approx(std::abs(exact)).epsilon(1e-9));
            }
    }
}

TEST_CASE("the run's first-order reflections are left out of the table", "[building_paths]")
{
    const Buildings city({box(0.0, -100.0, 20.0, 100.0, 30.0)});
    const Vec3 microphone(-10.0, 0.0, 1.5);
    const auto run = facing_facades(city, microphone, 150.0);
    REQUIRE(run.size() == 1);
    BuildingPathOptions options;
    options.diffraction = false;
    const Vec3 u = direction(170.0, 20.0); // from the street side, so the wall reflects it
    const BuildingReceiver with(city, microphone, run, options), without(city, microphone, {}, options);
    auto count = [](const std::vector<BuildingRay>& rays, bool run_evaluates) {
        std::size_t n = 0;
        for (const auto& r : rays)
            n += r.run_evaluates == run_evaluates;
        return n;
    };
    const auto a = with.rays(u), b = without.rays(u);
    // Direct and floor rays, and the wall's reflection with its floor partner.
    REQUIRE(a.size() == 4);
    REQUIRE(count(a, true) == 2);
    REQUIRE(count(b, true) == 0);
    std::vector<double> ra(1), rb(1);
    const std::vector<double> reflection(city.facades().size(), 0.5);
    with.ratio(u, reflection, 1.0, {1000.0}, 343.2, ra.data());
    without.ratio(u, reflection, 1.0, {1000.0}, 343.2, rb.data());
    REQUIRE(ra[0] == Catch::Approx(1.0).epsilon(1e-12));
    REQUIRE(rb[0] == Catch::Approx(1.5).epsilon(1e-12));
}

TEST_CASE("local_direction inverts compass_direction", "[building_paths]")
{
    for (double frame : {0.0, 90.0, 210.0})
        for (double azimuth : {0.0, 45.0, 123.0, 359.0})
            for (double elevation : {0.0, 12.5, 89.0})
            {
                const Vec3 u = local_direction(azimuth, elevation, frame);
                REQUIRE(u.norm() == Catch::Approx(1.0).epsilon(1e-14));
                const CompassDirection back = compass_direction(u, frame);
                // Modulo a turn: a bearing of 0 can come back as 360 by round-off.
                REQUIRE(std::remainder(back.azimuth - azimuth, 360.0) == Catch::Approx(0.0).margin(1e-9));
                REQUIRE(back.elevation == Catch::Approx(elevation).margin(1e-9));
            }
}

TEST_CASE("building_gains over long buildings is ideal_canyon_gains", "[building_paths]")
{
    // Wall A at y = 0 and wall B at y = width: looking west along the axis,
    // A is on the left.  The local +x is east.
    const IdealCanyon c {24.0, 9.0, 20.0, 35.0, 1.2};
    const Buildings city = canyon_buildings(c);
    const Vec3 microphone(0.0, c.offset, c.microphone_height);
    BuildingPathOptions options;
    options.max_order = 400;
    options.max_diffracted_order = 40;
    const BuildingReceiver receiver(city, microphone, facing_facades(city, microphone, options.radius), options);
    UrbanTable::Axes axes;
    axes.azimuths = {0.0, 30.0, 60.0, 120.0, 150.0, 180.0, 210.0, 240.0, 300.0, 330.0};
    axes.elevations = {3.0, 10.0, 30.0, 60.0};
    axes.bands = {200.0, 800.0};
    const std::vector<double> reflection(city.facades().size() * axes.bands.size(), 0.85);
    const std::vector<double> walls(axes.bands.size(), 0.85);
    const auto general = building_gains(receiver, axes, 90.0, reflection, 0.9, 343.2, 2, nullptr, 0.0);
    const auto ideal = ideal_canyon_gains(c, 270.0, axes, walls, walls, 0.9, 343.2, true, 1);
    REQUIRE(general.size() == ideal.size());
    for (std::size_t i = 0; i < general.size(); i++)
        REQUIRE(general[i] == Catch::Approx(ideal[i]).margin(1e-9));
}

TEST_CASE("the beams find every diffracted ray the unpruned search does", "[building_paths]")
{
    // A district of blocks of different heights, the microphone in a street.
    std::vector<Footprint> blocks;
    const double heights[] {12.0, 31.0, 18.0, 26.0, 40.0, 15.0, 22.0, 35.0, 28.0};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            blocks.push_back(box(i * 70.0, j * 70.0, i * 70.0 + 50.0, j * 70.0 + 50.0, heights[3 * i + j]));
    // An irregular one too, at an angle.
    Footprint skew;
    skew.x = {230.0, 280.0, 290.0, 240.0};
    skew.y = {20.0, 35.0, 85.0, 70.0};
    skew.top = 24.0;
    blocks.push_back(skew);
    const Buildings city(blocks);
    for (const Vec3 microphone : {Vec3(60.0, 95.0, 1.5), Vec3(130.0, 60.0, 0.0), Vec3(215.0, 30.0, 4.0)})
    {
        BuildingPathOptions pruned, full;
        pruned.max_diffracted_order = full.max_diffracted_order = 3;
        full.beams = false;
        const auto run = facing_facades(city, microphone, 100.0);
        const BuildingReceiver a(city, microphone, run, pruned), b(city, microphone, run, full);
        REQUIRE(a.images().size() < b.images().size());
        std::size_t rays = 0;
        for (double azimuth = 0.0; azimuth < 360.0; azimuth += 22.5)
            for (double elevation : {2.0, 9.0, 25.0, 50.0})
            {
                const Vec3 u = direction(azimuth, elevation);
                auto key = [](const std::vector<BuildingRay>& all) {
                    std::vector<std::pair<double, std::size_t>> k;
                    for (const auto& r : all)
                        k.emplace_back(r.length, r.facades.size() * 4 + r.floor * 2 + r.diffracted);
                    std::sort(k.begin(), k.end());
                    return k;
                };
                const auto ka = key(a.rays(u)), kb = key(b.rays(u));
                INFO("microphone " << microphone << " azimuth " << azimuth << " elevation " << elevation);
                REQUIRE(ka.size() == kb.size());
                for (std::size_t i = 0; i < ka.size(); i++)
                {
                    REQUIRE(ka[i].first == Catch::Approx(kb[i].first).epsilon(1e-12));
                    REQUIRE(ka[i].second == kb[i].second);
                }
                rays += ka.size();
            }
        REQUIRE(rays > 100);
    }
}

TEST_CASE("a street half a degree out of parallel keeps its long bounce chains", "[building_paths]")
{
    // Wall A along y = 0, wall B from y = 20 at x = -200 to y = 20 + 400 tan(0.5 deg) at x = 200, both
    // 30 m tall; a wave across the street at 3 degrees has to bounce some 27 times to climb out to
    // a microphone 1.5 m up.  Exact orientations make every one of those chains a composition of
    // dozens of different mirrors, and none is a candidate; within a degree the street is a canyon.
    Footprint b;
    const double rise = 400.0 * std::tan(0.5 * pi / 180.0);
    b.x = {-200.0, 200.0, 200.0, -200.0};
    b.y = {20.0, 20.0 + rise, 40.0 + rise, 40.0};
    b.top = 30.0;
    const Vec3 microphone(0.0, 10.0, 1.5);
    const Vec3 u = direction(90.0, 3.0);
    std::size_t longest_exact = 0, longest_tolerant = 0;
    for (const double tolerance : {1e-5, 1.0 * pi / 180.0})
    {
        const Buildings city({box(-200.0, -20.0, 200.0, 0.0, 30.0), b}, tolerance);
        BuildingPathOptions options;
        options.diffraction = false;
        const BuildingReceiver receiver(city, microphone, {}, options);
        std::size_t longest = 0;
        for (const auto& ray : receiver.rays(u))
            longest = std::max(longest, ray.facades.size());
        (tolerance < 1e-3 ? longest_exact : longest_tolerant) = longest;
    }
    INFO("longest chain: exact " << longest_exact << ", within a degree " << longest_tolerant);
    REQUIRE(longest_exact < 5);
    REQUIRE(longest_tolerant > 20);
}
