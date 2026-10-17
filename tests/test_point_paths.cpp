#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include "building_paths.hpp"
#include "buildings.hpp"
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

const std::vector<double> bands {125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0};

std::vector<double> levels(const Buildings& city, const Vec3& source, const Vec3& microphone, double rho, double floor,
        const PointPathOptions& options = {})
{
    const PointSourceField field(city, source, microphone, options);
    const std::vector<double> reflection(city.facades().size() * bands.size(), rho);
    std::vector<double> out(bands.size());
    field.band_levels(reflection, floor, bands, 343.2, out.data());
    return out;
}

} // namespace

TEST_CASE("over open ground a point source's field is the direct ray and the floor's", "[point_paths]")
{
    const Buildings none;
    const Vec3 s(0.0, 0.0, 2.0), m(30.0, 4.0, 1.2);
    const PointSourceField field(none, s, m);
    REQUIRE(field.rays().size() == 2);
    const double r1 = (m - s).norm(), r2 = (Vec3(m.x, m.y, -m.z) - s).norm();
    const auto out = levels(none, s, m, 1.0, 0.8);
    const double half = std::pow(2.0, 1.0 / 6.0);
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        const double f1 = bands[b] / half, f2 = bands[b] * half;
        double energy = 0.0;
        for (int part = 0; part < 3; part++)
        {
            const double g1 = f1 * std::pow(f2 / f1, part / 3.0), g2 = f1 * std::pow(f2 / f1, (part + 1) / 3.0);
            energy += band_averaged_energy({1.0 / r1, std::sqrt(0.8) / r2}, {r1, r2}, 2.0 * pi * g1 / 343.2,
                              2.0 * pi * g2 / 343.2)
                    * (g2 - g1);
        }
        REQUIRE(out[b] == Catch::Approx(10.0 * std::log10(energy / (f2 - f1))).margin(1e-9));
    }
}

TEST_CASE("the point-source field is reciprocal", "[point_paths]")
{
    // A small district: rays reflect, diffract at roof edges and corners, at
    // edges lit only by a reflection, and at two edges in turn.
    const Buildings city({box(0.0, 0.0, 6.0, 2.4, 7.8), box(8.0, 4.0, 10.4, 16.0, 7.8), box(-2.0, 8.0, 4.0, 10.4, 5.2),
            box(0.0, 14.0, 2.4, 20.0, 7.8)});
    const std::vector<std::pair<Vec3, Vec3>> pairs {{Vec3(3.0, 12.0, 2.0), Vec3(4.0, -3.0, 1.2)},
            {Vec3(6.0, 6.0, 2.0), Vec3(13.0, 10.0, 1.2)}, {Vec3(-4.0, 4.0, 1.5), Vec3(6.0, 18.0, 1.2)},
            {Vec3(-5.0, 4.0, 2.0), Vec3(6.0, 6.0, 1.2)}};
    std::size_t wedges[3] {0, 0, 0}, shared = 0;
    for (const auto& [s, m] : pairs)
    {
        const auto forward = levels(city, s, m, 0.9, 1.0), back = levels(city, m, s, 0.9, 1.0);
        for (std::size_t b = 0; b < bands.size(); b++)
            REQUIRE(forward[b] == Catch::Approx(back[b]).margin(1e-9));
        const PointSourceField field(city, s, m);
        for (const auto& ray : field.rays())
        {
            wedges[ray.wedges]++;
            shared += ray.shared;
            // The path drawn is the ray: it runs from the source to the
            // microphone through its points, and off the floor's mirror it is
            // as long as the ray.
            REQUIRE(ray.path.size() >= 2);
            REQUIRE((ray.path.front() - s).norm() == 0.0);
            REQUIRE((ray.path.back() - m).norm() == 0.0);
            REQUIRE(ray.path.size() == 2 + ray.facades.size() + static_cast<std::size_t>(ray.wedges));
            if (ray.floor == 0)
            {
                double length = 0.0;
                for (std::size_t v = 1; v < ray.path.size(); v++)
                    length += (ray.path[v] - ray.path[v - 1]).norm();
                REQUIRE(length == Catch::Approx(ray.length).epsilon(1e-9));
            }
        }
    }
    // Every kind was exercised, the uniform form over a shared face too.
    REQUIRE(wedges[0] > 0);
    REQUIRE(wedges[1] > 0);
    REQUIRE(wedges[2] > 0);
    REQUIRE(shared > 0);
}

TEST_CASE("a far point source is the plane wave", "[point_paths]")
{
    // One building, the microphone near it; the source a thousand kilometres
    // off in the direction u.  Scaled by the range, the point source's field
    // is the plane-wave builder's, every ray of it.
    const Buildings city({box(0.0, 0.0, 30.0, 12.0, 20.0)});
    const Vec3 m(-6.0, 5.0, 1.5);
    const double range = 1e6;
    BuildingPathOptions plane;
    plane.max_diffracted_order = 2;
    const BuildingReceiver receiver(city, m, {}, plane);
    PointPathOptions point;
    point.double_diffraction = false;
    point.source_reflections = false;
    point.max_edge_reflections = 2;
    point.radius = 2e6;
    const std::vector<double> reflection(city.facades().size() * bands.size(), 0.85);
    for (const Vec3 u : {local_direction(250.0, 10.0, 90.0), local_direction(80.0, 30.0, 90.0),
                 local_direction(120.0, 5.0, 90.0)})
    {
        std::vector<double> ratio(bands.size());
        receiver.ratio(u, reflection, 0.9, bands, 343.2, ratio.data());
        const auto far = levels(city, m + range * u, m, 0.85, 0.9, point);
        // The plane wave's total is its ratio times open ground's band energy.
        const auto open = levels(Buildings {}, m + range * u, m, 0.85, 0.9, point);
        for (std::size_t b = 0; b < bands.size(); b++)
            REQUIRE(far[b] - open[b] == Catch::Approx(10.0 * std::log10(ratio[b])).margin(2e-3));
    }
}

TEST_CASE("the field is continuous where a thick building's roof plane hides its far edge", "[point_paths]")
{
    // A long building 10 m thick and 8 m high, the source low before it, the
    // microphone behind it crossing the roof plane: above, the near edge's ray
    // grazes the roof and reaches it; below, only the ray diffracted at both
    // edges does.  The uniform form over the roof (double_wedge.hpp) carries
    // the one into the other.
    const Buildings city({box(-200.0, 0.0, 200.0, 10.0, 8.0)});
    const Vec3 s(0.0, -15.0, 2.0);
    PointPathOptions options;
    options.max_reflections = 0;
    options.max_edge_reflections = 0;
    for (const double y : {14.0, 25.0})
    {
        const auto above = levels(city, s, Vec3(0.0, y, 8.0 + 1e-3), 1.0, 0.0, options);
        const auto under = levels(city, s, Vec3(0.0, y, 8.0 - 1e-3), 1.0, 0.0, options);
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            INFO("microphone " << y << " m behind, band " << bands[b]);
            REQUIRE(std::abs(above[b] - under[b]) < 0.5);
        }
    }
}

TEST_CASE("over a long block the point-source field is the boundary-element solution's", "[point_paths]")
{
    // A block 2.44 m thick and 5.18 m high -- two stacked shipping containers,
    // DigiCity's (FIDELITY) -- on rigid ground, long enough to be two-
    // dimensional, the source and the microphone in one plane across it.
    // There every ray is a 2-D ray times 1/sqrt(its length), and the field in
    // that plane is the line source's: the exact values are a boundary-element
    // solution for one (tools/thick_barrier_bem.py, 60 elements a wavelength,
    // converged to 0.05 dB here), in dB re the free field.  Sources low, near
    // grazing and above the roof; microphones deep in the shadow, at the roof
    // plane and above it.  The uniform form is within 1 dB everywhere; the
    // cascade misses by up to 10.
    struct Case
    {
        double f, source_y, source_z, behind, z, exact;
    };
    const std::vector<Case> cases {
        {125, -25, 4.5, 5, 1.2, -17.01},
        {125, -25, 4.5, 5, 4.8, -7.37},
        {125, -25, 4.5, 5, 5.26, -7.20},
        {125, -25, 4.5, 5, 6.5, -0.48},
        {125, -25, 4.5, 25, 1.2, -6.88},
        {125, -25, 4.5, 25, 4.8, -5.72},
        {125, -25, 4.5, 25, 5.26, -4.29},
        {125, -25, 4.5, 25, 6.5, -2.38},
        {125, -15, 6.5, 5, 1.2, -19.05},
        {125, -15, 6.5, 5, 4.8, -9.04},
        {125, -15, 6.5, 25, 1.2, -6.60},
        {125, -15, 6.5, 25, 4.8, -4.40},
        {125, -10, 2, 5, 1.2, -19.70},
        {125, -10, 2, 5, 4.8, -11.46},
        {125, -10, 2, 5, 5.26, -12.13},
        {125, -10, 2, 5, 6.5, -5.35},
        {125, -10, 2, 25, 1.2, -10.42},
        {125, -10, 2, 25, 4.8, -10.90},
        {125, -10, 2, 25, 5.26, -9.35},
        {125, -10, 2, 25, 6.5, -7.42},
        {250, -25, 4.5, 5, 1.2, -18.15},
        {250, -25, 4.5, 5, 4.8, -11.96},
        {250, -25, 4.5, 5, 5.26, -7.76},
        {250, -25, 4.5, 5, 6.5, -5.11},
        {250, -25, 4.5, 25, 1.2, -17.58},
        {250, -25, 4.5, 25, 4.8, -12.95},
        {250, -25, 4.5, 25, 5.26, -10.00},
        {250, -25, 4.5, 25, 6.5, -5.81},
        {250, -15, 6.5, 5, 1.2, -10.90},
        {250, -15, 6.5, 5, 4.8, -4.96},
        {250, -15, 6.5, 25, 1.2, -9.85},
        {250, -15, 6.5, 25, 4.8, -5.93},
        {250, -10, 2, 5, 1.2, -20.69},
        {250, -10, 2, 5, 4.8, -16.45},
        {250, -10, 2, 5, 5.26, -11.72},
        {250, -10, 2, 5, 6.5, -9.54},
        {250, -10, 2, 25, 1.2, -20.80},
        {250, -10, 2, 25, 4.8, -20.37},
        {250, -10, 2, 25, 5.26, -17.72},
        {250, -10, 2, 25, 6.5, -11.11},
        {500, -25, 4.5, 5, 1.2, -32.26},
        {500, -25, 4.5, 5, 4.8, -8.09},
        {500, -25, 4.5, 5, 5.26, -6.71},
        {500, -25, 4.5, 5, 6.5, 1.20},
        {500, -25, 4.5, 25, 1.2, -14.40},
        {500, -25, 4.5, 25, 4.8, -7.17},
        {500, -25, 4.5, 25, 5.26, -6.10},
        {500, -25, 4.5, 25, 6.5, -3.73},
        {500, -15, 6.5, 5, 1.2, -26.98},
        {500, -15, 6.5, 5, 4.8, -2.92},
        {500, -15, 6.5, 25, 1.2, -6.93},
        {500, -15, 6.5, 25, 4.8, -0.34},
        {500, -10, 2, 5, 1.2, -37.68},
        {500, -10, 2, 5, 4.8, -15.04},
        {500, -10, 2, 5, 5.26, -14.99},
        {500, -10, 2, 5, 6.5, -8.12},
        {500, -10, 2, 25, 1.2, -19.90},
        {500, -10, 2, 25, 4.8, -15.47},
        {500, -10, 2, 25, 5.26, -14.75},
        {500, -10, 2, 25, 6.5, -14.35},
        {1000, -25, 4.5, 5, 1.2, -19.87},
        {1000, -25, 4.5, 5, 4.8, -10.00},
        {1000, -25, 4.5, 5, 5.26, -7.61},
        {1000, -25, 4.5, 5, 6.5, -1.61},
        {1000, -25, 4.5, 25, 1.2, -24.09},
        {1000, -25, 4.5, 25, 4.8, -10.02},
        {1000, -25, 4.5, 25, 5.26, -8.09},
        {1000, -25, 4.5, 25, 6.5, -5.26},
        {1000, -15, 6.5, 5, 1.2, -15.40},
        {1000, -15, 6.5, 5, 4.8, -4.28},
        {1000, -15, 6.5, 25, 1.2, -16.63},
        {1000, -15, 6.5, 25, 4.8, -1.34},
        {1000, -10, 2, 5, 1.2, -25.89},
        {1000, -10, 2, 5, 4.8, -17.77},
        {1000, -10, 2, 5, 5.26, -16.37},
        {1000, -10, 2, 5, 6.5, -9.36},
        {1000, -10, 2, 25, 1.2, -31.45},
        {1000, -10, 2, 25, 4.8, -20.50},
        {1000, -10, 2, 25, 5.26, -17.71},
        {1000, -10, 2, 25, 6.5, -18.45},
    };
    const double width = 2.44, height = 5.18;
    Footprint block;
    block.x = {-2000.0, 2000.0, 2000.0, -2000.0};
    block.y = {0.0, 0.0, width, width};
    block.top = height;
    const Buildings city({block});
    const std::vector<double> reflection(city.facades().size(), 1.0);
    double worst[2] {0.0, 0.0}, square = 0.0;
    for (const Case& c : cases)
    {
        const double k = 2.0 * pi * c.f / 343.2;
        const Vec3 s(0.0, c.source_y, c.source_z), m(0.0, width + c.behind, c.z);
        const double free = 1.0 / std::sqrt((m - s).norm());
        for (int uniform = 0; uniform < 2; uniform++)
        {
            PointPathOptions options;
            options.shared_face = uniform == 1;
            const PointSourceField field(city, s, m, options);
            std::complex<double> sum(0.0, 0.0);
            for (const PointRay& ray : field.rays())
                sum += point_ray_amplitude(ray, reflection, 1.0, k) * std::sqrt(ray.length)
                        * std::polar(1.0, k * ray.length);
            const double error = 20.0 * std::log10(std::abs(sum) / free) - c.exact;
            worst[uniform] = std::max(worst[uniform], std::abs(error));
            if (uniform == 1)
            {
                square += error * error;
                INFO(c.f << " Hz, source " << c.source_y << " m, " << c.source_z << " m up, microphone " << c.behind
                          << " m behind, " << c.z << " m up");
                REQUIRE(std::abs(error) < 1.0);
            }
        }
    }
    REQUIRE(std::sqrt(square / cases.size()) < 0.3);
    REQUIRE(worst[0] > 5.0);
}
