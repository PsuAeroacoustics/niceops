#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "building_input.hpp"
#include "buildings.hpp"
#include "ground_effect.hpp"
#include "point_paths.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT

// A source or a microphone may stand on a surface of the district -- a
// microphone flush on a roof or a wall, a loudspeaker on a roof -- as well as
// on the floor or in the air.  On a wall or a level roof the surface reflects
// at the point: every ray arriving there has a partner off it at the point,
// of one length (PointRay::source_surface).

namespace {

Footprint box(double x0, double y0, double x1, double y1, double top)
{
    Footprint f;
    f.x = {x0, x1, x1, x0};
    f.y = {y0, y0, y1, y1};
    f.top = top;
    return f;
}

const std::vector<double> bands {125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0};
const double c = 343.2;

std::vector<double> levels(const Buildings& city, const Vec3& source, const Vec3& microphone, double rho = 1.0,
        const FloorModel& floor = FloorModel(1.0))
{
    const PointSourceField field(city, source, microphone);
    std::vector<double> out(bands.size());
    field.band_levels(std::vector<double>(city.facades().size() * bands.size(), rho), floor, bands, c, out.data());
    return out;
}

} // namespace

// What a point lies on: a wall, a level roof (a flat one, or a mansard's deck),
// a slope, an edge (a roof's edge, a corner, a ridge), or nothing -- the air
// and the floor alike.
TEST_CASE("a point lies on a wall, a roof, a slope or an edge", "[microphone_surfaces]")
{
    Footprint mansard = box(30.0, 0.0, 42.0, 20.0, 9.0);
    mansard.roof = {{1.0, 0.0, -24.0}, {0.4, 0.0, -5.1}, {-0.4, 0.0, 28.5}, {-1.0, 0.0, 54.0}, {0.0, 0.0, 9.0}};
    const Buildings city({box(0.0, 0.0, 20.0, 20.0, 10.0), mansard});
    using Kind = Buildings::Support::Kind;
    REQUIRE(city.support(Vec3(10.0, 10.0, 10.0)).kind == Kind::roof);
    REQUIRE(city.support(Vec3(10.0, 10.0, 10.0 + 5e-7)).kind == Kind::roof);
    REQUIRE(city.support(Vec3(10.0, 10.0, 10.0 + 5e-6)).kind == Kind::none);
    REQUIRE(city.support(Vec3(20.0, 10.0, 5.0)).kind == Kind::facade);
    REQUIRE(city.support(Vec3(20.0, 10.0, 0.0)).kind == Kind::facade);
    REQUIRE(city.support(Vec3(20.0, 10.0, 10.0)).kind == Kind::edge);
    REQUIRE(city.support(Vec3(20.0, 20.0, 5.0)).kind == Kind::edge);
    REQUIRE(city.support(Vec3(25.0, 10.0, 0.0)).kind == Kind::none);
    REQUIRE(city.support(Vec3(25.0, 10.0, 3.0)).kind == Kind::none);
    // The mansard: its lower slope, its deck, and the break between its slopes.
    REQUIRE(city.support(Vec3(31.0, 10.0, 7.0)).kind == Kind::slope);
    REQUIRE(city.support(Vec3(36.0, 10.0, 9.0)).kind == Kind::roof);
    REQUIRE(city.support(Vec3(31.5, 10.0, 7.5)).kind == Kind::edge);
}

// A point on a wall or a level roof hears what one a hair off it hears, where
// the partners are traced: over the roof's edge, round the wall's corner, in
// the building's shadow, overhead; and a source on a roof is heard as one a
// hair above it.  Rigid surfaces: a diffracted ray's wedge holds its face's
// reflection rigidly, and so does the end on that face.  Within 0.1 dB (0.062
// measured, at 125 Hz in the building's shadow: the backscatter off an edge
// back along its face, which a hair off the face carries and the end on the
// face leaves out); before, a microphone on a roof had no doubling, and none
// at all behind the roof's edge, and one on a wall none either.
TEST_CASE("an end on a roof or a wall is one a hair off it", "[microphone_surfaces]")
{
    const Buildings city({box(0.0, 0.0, 20.0, 20.0, 10.0)});
    struct Case
    {
        Vec3 source, on, off;
    };
    const double hair = 1e-5;
    const std::vector<Case> cases {
        // On the roof: overhead, behind its edge, obliquely behind it.
        {{10.0, 10.0, 60.0}, {10.0, 10.0, 10.0}, {10.0, 10.0, 10.0 + hair}},
        {{-30.0, 10.0, 6.0}, {5.0, 10.0, 10.0}, {5.0, 10.0, 10.0 + hair}},
        {{-30.0, 25.0, 6.0}, {5.0, 10.0, 10.0}, {5.0, 10.0, 10.0 + hair}},
        {{10.0, -30.0, 3.0}, {5.0, 10.0, 10.0}, {5.0, 10.0, 10.0 + hair}},
        // On a wall: facing it, round its corner, behind the building, oblique.
        {{60.0, 10.0, 5.0}, {20.0, 10.0, 5.0}, {20.0 + hair, 10.0, 5.0}},
        {{60.0, 40.0, 5.0}, {20.0, 10.0, 5.0}, {20.0 + hair, 10.0, 5.0}},
        {{-30.0, 10.0, 4.0}, {20.0, 10.0, 5.0}, {20.0 + hair, 10.0, 5.0}},
        {{10.0, -30.0, 12.0}, {20.0, 10.0, 5.0}, {20.0 + hair, 10.0, 5.0}},
    };
    double worst = 0.0;
    for (const Case& k : cases)
    {
        INFO("source " << k.source.x << ", " << k.source.y << ", " << k.source.z << "; microphone " << k.on.x << ", "
                       << k.on.y << ", " << k.on.z);
        const auto on = levels(city, k.source, k.on), off = levels(city, k.source, k.off);
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            REQUIRE(std::isfinite(on[b]));
            worst = std::max(worst, std::abs(on[b] - off[b]));
        }
    }
    INFO("worst " << worst << " dB");
    REQUIRE(worst < 0.1);
    // A source on the roof, heard low behind the building.
    const Vec3 low(40.0, 10.0, 1.5);
    const auto on = levels(city, Vec3(15.0, 10.0, 10.0), low), off = levels(city, Vec3(15.0, 10.0, 10.0 + hair), low);
    for (std::size_t b = 0; b < bands.size(); b++)
        REQUIRE(on[b] == Catch::Approx(off[b]).margin(0.01));
}

// An absorbing roof reflects at a microphone on it as it does a hair above:
// overhead, (1 + sqrt(rho))^2 the free field's energy, where the ray off the
// roof is traced above it.
TEST_CASE("an absorbing roof doubles a microphone on it by its own reflection", "[microphone_surfaces]")
{
    const Buildings city({box(0.0, 0.0, 200.0, 200.0, 10.0)});
    const Vec3 source(100.0, 100.0, 60.0), on(100.0, 100.0, 10.0), off(100.0, 100.0, 10.0 + 1e-5);
    const double rho = 0.6;
    const auto at = levels(city, source, on, rho), above = levels(city, source, off, rho);
    const double expected = -20.0 * std::log10(50.0) + 20.0 * std::log10(1.0 + std::sqrt(rho));
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        REQUIRE(at[b] == Catch::Approx(above[b]).margin(0.01));
        REQUIRE(at[b] == Catch::Approx(expected).margin(0.05));
    }
}

// A plate on a roof: each ray takes P_d + R P_r, R the roof's amplitude
// reflection, where the floor's plate takes P_d and its partner off the floor
// Q P_r.  With P_d and P_r the same at every angle the level moves by
// 20 lg |P_d + R P_r| / (1 + R) from the bare roof's.
TEST_CASE("a plate on a roof takes the roof's reflection", "[microphone_surfaces]")
{
    const Buildings city({box(0.0, 0.0, 200.0, 200.0, 10.0)});
    const Vec3 source(100.0, 100.0, 60.0), microphone(110.0, 95.0, 10.0);
    const double rho = 0.81, r = std::sqrt(rho);
    const std::vector<double> middles = sub_band_middles(bands, c);
    FloorModel ground;
    ground.admittance = ground_impedance(GroundModel {}, middles, c);
    const std::complex<double> pd(0.8, 0.1), pr(0.6, -0.2);
    FloorModel plate = ground;
    plate.plate = [&](double, double, std::complex<double>* direct, std::complex<double>* image) {
        std::fill_n(direct, middles.size(), pd);
        std::fill_n(image, middles.size(), pr);
    };
    plate.surface.assign(middles.size(), r);
    const auto bare = levels(city, source, microphone, rho, ground);
    const auto plated = levels(city, source, microphone, rho, plate);
    const double expected = 20.0 * std::log10(std::abs(pd + r * pr) / (1.0 + r));
    for (std::size_t b = 0; b < bands.size(); b++)
        REQUIRE(plated[b] - bare[b] == Catch::Approx(expected).margin(1e-9));
}

// On a slope a point takes the slope's reflection there, as on a level roof:
// it hears what a point a hair above the slope hears, the slope's rays
// traced -- in front of it, overhead, from below the eave, beyond the ridge
// -- and so does a source on it.  Rigid surfaces, as above.  Within 0.2 dB
// (0.155 measured, at 125 Hz behind the ridge: the ridge's wave down the
// slope, off the eave and back up along it, which a hair off the face carries
// and the end on the face leaves out, as on a wall).
TEST_CASE("an end on a slope is one a hair above it", "[microphone_surfaces]")
{
    Footprint gable = box(0.0, 0.0, 10.0, 20.0, 9.0);
    gable.roof = {{0.6, 0.0, 6.0}, {-0.6, 0.0, 12.0}};
    const Buildings city({gable});
    const double scale = std::sqrt(1.36), hair = 1e-5;
    const Vec3 n(-0.6 / scale, 0.0, 1.0 / scale);
    const Vec3 on(2.0, 10.0, 7.2), off = on + hair * n;
    REQUIRE(city.support(on).kind == Buildings::Support::Kind::slope);
    double worst = 0.0;
    for (const Vec3& source : {Vec3(2.0, 10.0, 40.0), Vec3(-30.0, 10.0, 20.0), Vec3(-30.0, 10.0, 4.0),
                 Vec3(30.0, 10.0, 12.0), Vec3(2.0, -30.0, 15.0)})
    {
        INFO("source " << source.x << ", " << source.y << ", " << source.z);
        const auto a = levels(city, source, on), b = levels(city, source, off);
        for (std::size_t i = 0; i < bands.size(); i++)
        {
            REQUIRE(std::isfinite(a[i]));
            worst = std::max(worst, std::abs(a[i] - b[i]));
        }
        // And the source on the slope.
        const auto c = levels(city, on, source), d = levels(city, off, source);
        for (std::size_t i = 0; i < bands.size(); i++)
            worst = std::max(worst, std::abs(c[i] - d[i]));
    }
    INFO("worst " << worst << " dB");
    REQUIRE(worst < 0.2);
}

// Where no one surface's reflection stands for it -- an edge: an eave, a
// ridge, a corner -- an end is refused, by the field and by the tools'
// placement check; on a slope it is fine, and below the roof it is inside.
TEST_CASE("an end on an edge is refused", "[microphone_surfaces]")
{
    Footprint gable = box(0.0, 0.0, 10.0, 20.0, 9.0);
    gable.roof = {{0.6, 0.0, 6.0}, {-0.6, 0.0, 12.0}};
    const Buildings city({gable});
    const std::vector<std::string> ids {"1"};
    const Vec3 source(-20.0, 10.0, 2.0);
    const Vec3 slope(2.0, 10.0, 7.2), edge(0.0, 10.0, 6.0), ridge(5.0, 10.0, 9.0), above(2.0, 10.0, 7.7);
    REQUIRE_NOTHROW(PointSourceField(city, source, slope));
    REQUIRE_THROWS_AS(PointSourceField(city, source, edge), std::invalid_argument);
    REQUIRE_THROWS_AS(PointSourceField(city, source, ridge), std::invalid_argument);
    REQUIRE_NOTHROW(PointSourceField(city, source, above));
    REQUIRE(misplaced_point(city, ids, slope).empty());
    REQUIRE_THAT(misplaced_point(city, ids, edge), Catch::Matchers::ContainsSubstring("edge"));
    REQUIRE(misplaced_point(city, ids, above).empty());
    REQUIRE_THAT(misplaced_point(city, ids, Vec3(2.0, 10.0, 5.0)), Catch::Matchers::ContainsSubstring("inside"));
}

// A pitched roof's slope reflects as a mirror in its plane: the ray off it
// meets the plane where that plane is the roof, at equal angles, and is as
// long as the line from the source's image to the microphone.
TEST_CASE("a slope reflects a specular ray", "[microphone_surfaces]")
{
    Footprint shed = box(0.0, 0.0, 200.0, 200.0, 60.0);
    shed.roof = {{0.2, 0.0, 10.0}};
    const Buildings city({shed});
    const Vec3 s(60.0, 100.0, 40.0), m(120.0, 110.0, 45.0);
    const PointSourceField field(city, s, m);
    const double scale = std::sqrt(1.04);
    const Vec3 n(-0.2 / scale, 0.0, 1.0 / scale);
    const double ds = (s.z - (0.2 * s.x + 10.0)) / scale;
    const Vec3 image = s - (2.0 * ds) * n;
    std::size_t found = 0;
    for (const PointRay& ray : field.rays())
    {
        if (ray.roof_facades.empty() || ray.wedges > 0)
            continue;
        found++;
        REQUIRE(ray.path.size() == 3);
        const Vec3& p = ray.path[1];
        REQUIRE(p.z == Catch::Approx(0.2 * p.x + 10.0).margin(1e-9));
        REQUIRE(ray.length == Catch::Approx((m - image).norm()).epsilon(1e-12));
        REQUIRE(ray.length == Catch::Approx((p - s).norm() + (m - p).norm()).epsilon(1e-12));
        const Vec3 in = (p - s) / (p - s).norm(), out = (m - p) / (m - p).norm();
        REQUIRE(-in.dot(n) == Catch::Approx(out.dot(n)).epsilon(1e-12));
    }
    REQUIRE(found == 1);
}

// Where a slope's reflection meets its edge -- the eave below it, the ridge
// above -- the field is continuous: the edge's coefficient carries half the
// reflection either side of its boundary, which needs the reflection there.
// Before slopes reflected it stepped by the reflection's whole level.
TEST_CASE("a slope's reflection ends where its edges' boundaries are", "[microphone_surfaces]")
{
    Footprint gable = box(0.0, 0.0, 10.0, 20.0, 9.0);
    gable.roof = {{0.6, 0.0, 6.0}, {-0.6, 0.0, 12.0}};
    const Buildings city({gable});
    const Vec3 s(-20.0, 10.0, 20.0);
    const double scale = std::sqrt(1.36);
    const Vec3 n(-0.6 / scale, 0.0, 1.0 / scale);
    const std::vector<double> bands {250.0, 1000.0, 4000.0};
    for (const Vec3& at : {Vec3(0.0, 10.0, 6.0), Vec3(5.0, 10.0, 9.0)})
    {
        const Vec3 d = (at - s) / (at - s).norm();
        const Vec3 r = d - (2.0 * d.dot(n)) * n;
        // Across the boundary: normal to the edge's line (along y) and to r.
        const Vec3 across = Vec3(0.0, 1.0, 0.0).cross(r);
        const Vec3 m = at + 15.0 * r;
        const auto here = [&](double side) {
            const PointSourceField field(city, s, m + side * across);
            std::vector<double> out(bands.size());
            field.band_levels(std::vector<double>(city.facades().size() * bands.size(), 1.0), FloorModel(1.0), bands,
                    343.2, out.data());
            return out;
        };
        const auto a = here(1e-6), b = here(-1e-6);
        for (std::size_t i = 0; i < bands.size(); i++)
        {
            INFO("edge at x = " << at.x << ", " << bands[i] << " Hz");
            REQUIRE(a[i] == Catch::Approx(b[i]).margin(0.02));
        }
    }
}

// A slope's rays as the floor models read them: the reflection point on the
// plane at equal angles, a floor bounce on z = 0 at either end or not, each
// bounce recording the whole unfolded path as its leg (the slope a mirror,
// not a new source) and its grazing sine at the floor, and the legs summing
// to the ray's length.  An end on the floor is its own bounce: no zero-length
// leg is left in the path, so its launch and arrival are the slope's, which a
// directional source and a plate read.
TEST_CASE("a slope's rays carry their floor bounces and real directions", "[microphone_surfaces]")
{
    Footprint gable = box(0.0, 0.0, 10.0, 20.0, 9.0);
    gable.roof = {{0.6, 0.0, 6.0}, {-0.6, 0.0, 12.0}};
    const Buildings city({gable});
    std::size_t rays = 0, with_floor = 0, on_floor = 0;
    // A source high to the west, over the ridge onto the east slope, which
    // sends it east and a little down to microphones about 40 m out; and the
    // same turned round, the source on the floor.
    const std::vector<std::pair<Vec3, Vec3>> pairs {{{-20.0, 10.0, 40.0}, {42.0, 10.0, 0.0}},
            {{-20.0, 10.0, 40.0}, {40.0, 12.0, 1.5}}, {{-18.0, 7.0, 36.0}, {45.0, 9.0, 0.0}},
            {{42.0, 10.0, 0.0}, {-20.0, 10.0, 40.0}}, {{40.0, 12.0, 1.5}, {-20.0, 10.0, 40.0}}};
    for (const auto& [source, microphone] : pairs)
        {
            const PointSourceField field(city, source, microphone);
            for (const PointRay& ray : field.rays())
            {
                if (ray.roof_facades.empty() || ray.wedges > 0)
                    continue;
                rays++;
                REQUIRE(ray.path.front().x == source.x);
                REQUIRE(ray.path.back().x == microphone.x);
                double length = 0.0;
                for (std::size_t k = 0; k + 1 < ray.path.size(); k++)
                {
                    const double leg = (ray.path[k + 1] - ray.path[k]).norm();
                    REQUIRE(leg > 1e-6);
                    length += leg;
                }
                REQUIRE(length == Catch::Approx(ray.length).epsilon(1e-9));
                // The slope's point, and the mirror law there.
                std::size_t at = 0;
                for (std::size_t k = 1; k + 1 < ray.path.size(); k++)
                    if (ray.path[k].z > 1e-9)
                        at = k;
                REQUIRE(at > 0);
                const Vec3& p = ray.path[at];
                const RoofPlane& plane = p.x < 5.0 ? gable.roof[0] : gable.roof[1];
                REQUIRE(plane.above(p) == Catch::Approx(0.0).margin(1e-9));
                const Vec3 in = (p - ray.path[at - 1]) / (p - ray.path[at - 1]).norm();
                const Vec3 out = (ray.path[at + 1] - p) / (ray.path[at + 1] - p).norm();
                REQUIRE(-in.dot(plane.normal()) == Catch::Approx(out.dot(plane.normal())).epsilon(1e-9));
                for (int i = 0; i < ray.floor; i++)
                {
                    const FloorBounce& b = ray.floor_bounce[static_cast<std::size_t>(i)];
                    REQUIRE(b.leg == Catch::Approx(ray.length).epsilon(1e-12));
                    // The leg between the slope and the floor's side end.
                    const Vec3 d = b.microphone ? ray.path[at + 1] - p : ray.path[at - 1] - p;
                    REQUIRE(b.sin_grazing == Catch::Approx(std::abs(d.z) / d.norm()).epsilon(1e-9));
                }
                with_floor += ray.floor > 0;
                on_floor += source.z == 0.0 || microphone.z == 0.0;
                const Vec3 launch = point_ray_launch(ray);
                REQUIRE(launch.norm() == Catch::Approx(1.0));
                double elevation = 0.0, azimuth = 0.0;
                point_ray_arrival(ray, elevation, azimuth);
                const Vec3 last = ray.path[ray.path.size() - 1] - ray.path[ray.path.size() - 2];
                REQUIRE(elevation == Catch::Approx(std::atan2(std::abs(last.z), std::hypot(last.x, last.y))));
            }
        }
    INFO(rays << " slope rays, " << with_floor << " with the floor, " << on_floor << " with an end on it");
    REQUIRE(with_floor > 0);
    REQUIRE(on_floor > 0);
}

// Where another building's edge cuts off a slope's reflection -- the wave on
// its way to the slope, or on its way from it to the microphone -- the field
// is continuous: the edge's once-diffracted ray with the slope's reflection on
// that side is the partner that carries half of it either side of the
// boundary.  A box's top edge is put on the reflected path, and the
// microphone moved 0.1 mm either side of the boundary it casts (a
// micrometre is within the visibility test's contact).
// Without those partners the level stepped by the reflection's.
TEST_CASE("an edge that cuts off a slope's reflection leaves the field continuous", "[microphone_surfaces]")
{
    Footprint gable = box(0.0, 0.0, 10.0, 20.0, 9.0);
    gable.roof = {{0.6, 0.0, 6.0}, {-0.6, 0.0, 12.0}};
    const RoofPlane& east = gable.roof[1];
    const Vec3 s(-20.0, 10.0, 40.0), p(7.5, 10.0, east.at(7.5, 10.0));
    const Vec3 d = (p - s) / (p - s).norm();
    const Vec3 r = d - (2.0 * d.dot(east.normal())) * east.normal();
    const Vec3 m = p + 30.0 * r;
    const Vec3 across = Vec3(0.0, 1.0, 0.0).cross(r) / Vec3(0.0, 1.0, 0.0).cross(r).norm();
    // A box under the incoming leg, its east top edge on it, the leg grazing
    // over it; then one under the outgoing leg.
    const Vec3 before = s + 0.6 * (p - s), after = p + 0.5 * (m - p);
    for (const Vec3& edge : {before, after})
    {
        const Buildings city({gable, box(edge.x - 1.5, -20.0, edge.x, 40.0, edge.z)});
        // The slope's reflection is there on one side and not on the other.
        const auto slope_rays = [&](double side) {
            std::size_t n = 0;
            const PointSourceField field(city, s, m + side * across);
            for (const PointRay& ray : field.rays())
                n += !ray.roof_facades.empty() && ray.wedges == 0 && ray.floor == 0;
            return n;
        };
        REQUIRE(slope_rays(1e-4) != slope_rays(-1e-4));
        double worst = 0.0;
        const auto at = [&](double side) { return levels(city, s, m + side * across); };
        const auto a = at(1e-4), b = at(-1e-4);
        for (std::size_t i = 0; i < bands.size(); i++)
        {
            REQUIRE(std::isfinite(a[i]));
            worst = std::max(worst, std::abs(a[i] - b[i]));
        }
        INFO("the edge at x = " << edge.x << ": worst " << worst << " dB");
        REQUIRE(worst < 0.05);
    }
}
