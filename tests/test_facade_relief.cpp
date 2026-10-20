#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "buildings.hpp"
#include "diffuse_facades.hpp"
#include "facade_relief.hpp"
#include "grating_paths.hpp"
#include "point_paths.hpp"
#include "urban.hpp"
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

// A relief table with the same shares at every frequency and incidence.
std::string write_relief(const std::string& name, double period, const std::map<int, double>& shares)
{
    const std::string path = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") + "/" + name;
    std::ofstream out(path);
    out << "frequency_hz,incidence_deg,order,share,period_m\n";
    for (double f : {50.0, 20000.0})
        for (int psi = -85; psi <= 85; psi += 5)
            for (const auto& [m, s] : shares)
                out << f << "," << psi << "," << m << "," << s << "," << period << "\n";
    return path;
}

// The facade of a building whose front `point` stands before.
std::size_t facing(const Buildings& city, const Vec3& point, double x)
{
    for (std::size_t f = 0; f < city.facades().size(); f++)
    {
        const Facade& g = city.facades()[f];
        if (std::abs(g.x0 - x) < 1e-12 && std::abs(g.x1 - x) < 1e-12 && facade_side(g, point) > 0.0)
            return f;
    }
    return city.facades().size();
}

} // namespace

// The shares are the table's, renormalized over the orders that propagate:
// sin(theta_m) = -sin(psi) + m lambda / period within (-1, 1).
TEST_CASE("a relief's shares are the table's over the orders that propagate", "[facade_relief]")
{
    const FacadeRelief relief = FacadeRelief::read(write_relief("niceops_relief_a.csv", 1.0,
            {{-2, 0.1}, {-1, 0.2}, {0, 0.4}, {1, 0.2}, {2, 0.1}}));
    REQUIRE(relief.period() == 1.0);
    std::vector<std::pair<int, double>> shares;
    // At normal incidence and lambda / period = 0.4, all five propagate.
    relief.orders(0.0, 343.2 / 0.4, 343.2, shares);
    REQUIRE(shares.size() == 5);
    REQUIRE(relief.specular(0.0, 343.2 / 0.4, 343.2) == Catch::Approx(0.4).epsilon(1e-12));
    // At lambda / period = 0.6 the second orders (1.2) do not: 0.4 of 0.8.
    relief.orders(0.0, 343.2 / 0.6, 343.2, shares);
    REQUIRE(shares.size() == 3);
    REQUIRE(relief.specular(0.0, 343.2 / 0.6, 343.2) == Catch::Approx(0.5).epsilon(1e-12));
    // From sin(psi) = 0.5 at 0.6: -0.5 + m 0.6 for m = -2..2 is -1.7, -1.1, -0.5, 0.1, 0.7.
    relief.orders(0.5, 343.2 / 0.6, 343.2, shares);
    REQUIRE(shares.size() == 3);
    REQUIRE(shares.front().first == 0);
    REQUIRE(relief.specular(0.5, 343.2 / 0.6, 343.2) == Catch::Approx(0.4 / 0.7).epsilon(1e-12));
    // Below every order's cut-on the specular takes all.
    REQUIRE(relief.specular(0.0, 100.0, 343.2) == 1.0);
    REQUIRE(relief.random_incidence(100.0, 343.2) == Catch::Approx(0.0).margin(1e-12));
}

// The generalized law: on the plane, the outgoing direction's component along
// the facade is the incident one's plus the kick, and its vertical component
// the same.  Without a kick it is the specular point, and the tube spreads as
// from the image, 1 / L^2.  Source and microphone level with each other, the
// tube is the plan's: L per radian up times the width dS/dphi across, the
// prototype's closed form for one reflection.
TEST_CASE("a grating point obeys the generalized law of reflection", "[facade_relief]")
{
    Facade wall;
    wall.x0 = 0.0;
    wall.y0 = -50.0;
    wall.x1 = 0.0;
    wall.y1 = 50.0;
    wall.top = 30.0;
    // The front is the left of (0, -50) -> (0, 50): -x.
    const Vec3 a(-10.0, -5.0, 3.0), b(-15.0, 8.0, 1.5);
    REQUIRE(facade_side(wall, a) > 0.0);
    const Vec3 along(0.0, 1.0, 0.0);
    double x, z;
    REQUIRE(grating_point(wall, a, b, 0.0, x, z));
    const Vec3 image = mirror_in(wall, b);
    REQUIRE(grating_spreading(wall, a, Vec3(0.0, -50.0 + x, z), b, 0.0)
            == Catch::Approx(1.0 / (image - a).squaredNorm()).epsilon(1e-6));
    for (double kick : {-1.2, -0.5, 0.3, 0.9})
    {
        REQUIRE(grating_point(wall, a, b, kick, x, z));
        const Vec3 p(0.0, -50.0 + x, z);
        const Vec3 u1 = (p - a) / (p - a).norm(), u2 = (b - p) / (b - p).norm();
        REQUIRE(along.dot(u2) - along.dot(u1) == Catch::Approx(kick).margin(1e-10));
        REQUIRE(u2.z == Catch::Approx(u1.z).margin(1e-10));
    }
    // Level ends: the plan's tube.
    const Vec3 s(-4.0, -6.0, 2.0), m(-7.0, 9.0, 2.0);
    for (double kick : {-0.6, 0.0, 0.4})
    {
        REQUIRE(grating_point(wall, s, m, kick, x, z));
        const Vec3 p(0.0, -50.0 + x, z);
        REQUIRE(z == Catch::Approx(2.0).margin(1e-9));
        const double t1 = along.dot((p - s) / (p - s).norm()), t2 = t1 + kick;
        const double dx1 = 4.0, dx2 = 7.0;
        const double length = (p - s).norm() + (m - p).norm();
        const double dy = dx1 / std::pow(1.0 - t1 * t1, 1.5) + dx2 / std::pow(1.0 - t2 * t2, 1.5);
        const double width = dy * std::sqrt(1.0 - t1 * t1) * std::sqrt(1.0 - t2 * t2);
        REQUIRE(grating_spreading(wall, s, p, m, kick) == Catch::Approx(1.0 / (length * width)).epsilon(1e-5));
    }
}

// A relief's specular share is the specular ray's: with its order 0 taking a
// quarter (the other order, its period long, always propagating) the field is
// the plain wall's with that wall reflecting a quarter of the energy.
TEST_CASE("a relief's specular share scales the specular ray", "[facade_relief]")
{
    const Buildings city({box(0.0, -50.0, 1.0, 50.0, 30.0)});
    const Vec3 s(-10.0, -5.0, 3.0), m(-15.0, 8.0, 1.5);
    const std::size_t wall = facing(city, s, 0.0);
    REQUIRE(wall < city.facades().size());
    const FacadeRelief relief = FacadeRelief::read(write_relief("niceops_relief_b.csv", 1000.0, {{0, 0.25}, {1, 0.75}}));
    std::vector<const FacadeRelief*> reliefs(city.facades().size(), nullptr);
    reliefs[wall] = &relief;
    PointPathOptions with;
    with.reliefs = &reliefs;
    const std::vector<double> bands {250.0, 1000.0, 4000.0};
    std::vector<double> plain(city.facades().size() * bands.size(), 1.0), quarter = plain;
    for (std::size_t b = 0; b < bands.size(); b++)
        quarter[wall * bands.size() + b] = 0.25;
    std::vector<double> got(bands.size()), want(bands.size());
    PointSourceField(city, s, m, with).band_levels(plain, 1.0, bands, 343.2, got.data());
    PointSourceField(city, s, m).band_levels(quarter, 1.0, bands, 343.2, want.data());
    for (std::size_t b = 0; b < bands.size(); b++)
        REQUIRE(got[b] == Catch::Approx(want[b]).margin(1e-9));
}

// Off one wall, each grating ray meets it within its rectangle by the
// generalized law, leaves the source along its first leg, and carries its
// order's share at the incidence and frequency it meets the wall at times its
// tube's spreading; every order whose point lies on the wall is found, with
// the floor at either end or not.
TEST_CASE("grating rays leave on the incident cone in their orders' shares", "[facade_relief]")
{
    const Buildings city({box(0.0, -50.0, 1.0, 50.0, 30.0)});
    const Vec3 s(-10.0, -5.0, 3.0), m(-15.0, 8.0, 1.5);
    const std::size_t wall = facing(city, s, 0.0);
    const FacadeRelief relief = FacadeRelief::read(write_relief("niceops_relief_c.csv", 1.0,
            {{-2, 0.1}, {-1, 0.2}, {0, 0.4}, {1, 0.2}, {2, 0.1}}));
    std::vector<const FacadeRelief*> reliefs(city.facades().size(), nullptr);
    reliefs[wall] = &relief;
    PointPathOptions options;
    options.reliefs = &reliefs;
    const PathEnd source(city, s, true, options, edges_near(city, {m}, options.radius));
    const PathEnd microphone(city, m, false, options, edges_near(city, {s}, options.radius));
    GratingInput input;
    input.reliefs = &reliefs;
    input.bands = {1000.0};
    const std::vector<double> rho(city.facades().size(), 1.0);
    input.reflection = &rho;
    std::vector<GratingRay> rays;
    grating_rays(city, source, microphone, input, rays);
    const Facade& f = city.facades()[wall];
    const Vec3 along = Vec3(f.x1 - f.x0, f.y1 - f.y0, 0.0) / std::hypot(f.x1 - f.x0, f.y1 - f.y0);
    std::map<int, int> direct;
    for (const GratingRay& ray : rays)
    {
        REQUIRE(ray.order != 0);
        REQUIRE(ray.path.size() == 3);
        const Vec3& p = ray.path[1];
        const Vec3 pa = ray.source_floor ? Vec3(s.x, s.y, -s.z) : s, pb = ray.microphone_floor ? Vec3(m.x, m.y, -m.z) : m;
        const Vec3 u1 = (p - pa) / (p - pa).norm(), u2 = (pb - p) / (pb - p).norm();
        const double kick = ray.order * 343.2 / 1000.0;
        REQUIRE(along.dot(u2) - along.dot(u1) == Catch::Approx(kick).margin(1e-9));
        REQUIRE(u2.z == Catch::Approx(u1.z).margin(1e-9));
        REQUIRE(ray.length == Catch::Approx((p - pa).norm() + (pb - p).norm()).epsilon(1e-12));
        const Vec3 launch = ray.source_floor ? Vec3(u1.x, u1.y, -u1.z) : u1;
        REQUIRE((ray.launch - launch).norm() < 1e-12);
        std::vector<std::pair<int, double>> shares;
        const double horizontal = std::hypot(u1.x, u1.y);
        relief.orders(-along.dot(u1) / horizontal, 1000.0 * horizontal, 343.2, shares);
        double share = 0.0;
        for (const auto& [order, v] : shares)
            if (order == ray.order)
                share = v;
        REQUIRE(ray.energy == Catch::Approx(share * grating_spreading(f, pa, p, pb, kick)).epsilon(1e-12));
        if (!ray.source_floor && !ray.microphone_floor)
            direct[ray.order]++;
    }
    // lambda / period = 0.3432: orders -2..2 propagate wherever they meet the
    // wall, and each meets it once.
    REQUIRE(direct.size() == 4);
    for (const auto& [order, count] : direct)
        REQUIRE(count == 1);
    INFO(rays.size() << " grating rays");
    REQUIRE(rays.size() >= 8);
}

// Only a ray that diffracts nowhere takes a relief's specular share: one that
// diffracts reflects whole, since what it would scatter reaches nothing in the
// model.  A block beside the relief wall sends diffracted rays onto it.
TEST_CASE("a ray that diffracts reflects in a relief whole", "[facade_relief]")
{
    const Buildings city({box(0.0, -50.0, 1.0, 50.0, 30.0), box(-12.0, 0.0, -8.0, 4.0, 6.0)});
    const Vec3 s(-10.0, -5.0, 3.0), m(-11.0, 9.0, 1.5);
    const std::size_t wall = facing(city, s, 0.0);
    const FacadeRelief relief = FacadeRelief::read(write_relief("niceops_relief_d.csv", 1000.0, {{0, 0.25}, {1, 0.75}}));
    std::vector<const FacadeRelief*> reliefs(city.facades().size(), nullptr);
    reliefs[wall] = &relief;
    PointPathOptions options;
    options.reliefs = &reliefs;
    const PointSourceField field(city, s, m, options);
    std::size_t specular = 0, diffracted = 0;
    for (const PointRay& ray : field.rays())
    {
        const bool on_wall = std::find(ray.facades.begin(), ray.facades.end(), wall) != ray.facades.end();
        REQUIRE(ray.hits.size() == ray.facades.size());
        REQUIRE(meets_relief(ray, &reliefs) == (on_wall && ray.wedges == 0));
        specular += on_wall && ray.wedges == 0;
        diffracted += on_wall && ray.wedges > 0;
    }
    REQUIRE(specular > 0);
    REQUIRE(diffracted > 0);
}

// The diffuse field's incidence-dependent scattering, given the facades' own
// coefficient at every incidence, changes nothing; and a lone wall whose first
// scattering goes some other way sends the microphone nothing, its patches
// seeing no other.
TEST_CASE("the diffuse field takes s at each path's incidence and can leave the first scattering out",
        "[facade_relief]")
{
    const Buildings city({box(0.0, -10.0, 1.0, 10.0, 8.0)});
    const Vec3 s(-6.0, -3.0, 2.0), m(-5.0, 4.0, 1.5);
    const std::vector<double> bands {500.0, 2000.0};
    const std::size_t nf = city.facades().size();
    const std::vector<double> rho(nf * bands.size(), 0.9), scatter(nf * bands.size(), 0.4);
    DiffuseOptions plain;
    plain.patch_size = 1.0;
    DiffuseOptions at = plain;
    at.scattering_at = [](std::size_t, double, double, std::size_t) { return 0.4; };
    std::vector<double> a(bands.size()), b(bands.size());
    DiffuseFacades(city, m, rho, scatter, 1.0, bands, 343.2, plain).point_energy(s, 2, {}, a.data());
    DiffuseFacades(city, m, rho, scatter, 1.0, bands, 343.2, at).point_energy(s, 2, {}, b.data());
    for (std::size_t i = 0; i < bands.size(); i++)
    {
        REQUIRE(a[i] > 0.0);
        REQUIRE(b[i] == Catch::Approx(a[i]).epsilon(1e-12));
    }
    DiffuseOptions elsewhere = plain;
    elsewhere.first_scattering_elsewhere = [](std::size_t) { return true; };
    DiffuseFacades(city, m, rho, scatter, 1.0, bands, 343.2, elsewhere).point_energy(s, 2, {}, b.data());
    for (std::size_t i = 0; i < bands.size(); i++)
        REQUIRE(b[i] == Catch::Approx(0.0).margin(1e-15));
}
