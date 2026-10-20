#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "building_input.hpp"
#include "buildings.hpp"
#include "point_paths.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT

namespace {

// A gable house 10 m by 20 m, its ridge along y at x = 5: eaves 6 m up, the
// ridge 10 m, the roof's planes rising 0.8 a metre.
Footprint gable_house()
{
    Footprint f;
    f.x = {0.0, 10.0, 10.0, 0.0};
    f.y = {0.0, 0.0, 20.0, 20.0};
    f.roof = {{0.8, 0.0, 6.0}, {-0.8, 0.0, 14.0}};
    return f;
}

} // namespace

// The house's top is its ridge; its long walls stand to the eaves, its gable
// walls in two stretches up to the ridge and down; it diffracts at four
// corners as tall as the eaves, at its eaves, at its four rakes, and at its
// ridge, each wedge the angle its faces make: the eaves 90 degrees plus the
// pitch inside, the rakes square, the ridge 180 less twice the pitch.
TEST_CASE("a gable roof's walls, eaves, rakes and ridge", "[pitched_roofs]")
{
    const Buildings city({gable_house()});
    REQUIRE(city.footprints()[0].top == Catch::Approx(10.0).epsilon(1e-12));
    std::size_t rectangles = 0, sloped = 0;
    for (const Facade& f : city.facades())
    {
        if (!(f.top_start == f.top_start) || f.top_start == f.top_end)
        {
            rectangles++;
            REQUIRE(facade_top(f, Vec3(f.x0, f.y0, 0.0)) == Catch::Approx(6.0).epsilon(1e-12));
            continue;
        }
        sloped++;
        // A gable stretch runs from an eave (6) to the ridge (10).
        REQUIRE(std::min(f.top_start, f.top_end) == Catch::Approx(6.0).epsilon(1e-12));
        REQUIRE(std::max(f.top_start, f.top_end) == Catch::Approx(10.0).epsilon(1e-12));
        const Vec3 middle(0.5 * (f.x0 + f.x1), 0.5 * (f.y0 + f.y1), 0.0);
        REQUIRE(facade_top(f, middle) == Catch::Approx(8.0).epsilon(1e-12));
    }
    REQUIRE(rectangles == 2);
    REQUIRE(sloped == 4);
    const double pitch = std::atan(0.8);
    std::size_t corners = 0, eaves = 0, rakes = 0, ridges = 0;
    for (const Edge& e : city.edges())
    {
        if (e.vertical)
        {
            corners++;
            REQUIRE(e.length == Catch::Approx(6.0).epsilon(1e-12));
            REQUIRE(e.n == Catch::Approx(1.5).epsilon(1e-12));
            continue;
        }
        const Vec3 end = e.start + e.length * e.direction;
        if (std::abs(e.direction.z) < 1e-12 && std::abs(e.start.z - 10.0) < 1e-9)
        {
            ridges++;
            REQUIRE(e.length == Catch::Approx(20.0).epsilon(1e-12));
            REQUIRE(e.start.x == Catch::Approx(5.0).epsilon(1e-12));
            REQUIRE(e.n == Catch::Approx(1.0 + 2.0 * pitch / pi).epsilon(1e-12));
        }
        else if (std::abs(e.direction.z) < 1e-12)
        {
            eaves++;
            REQUIRE(e.start.z == Catch::Approx(6.0).epsilon(1e-12));
            REQUIRE(e.n == Catch::Approx(2.0 - (0.5 * pi + pitch) / pi).epsilon(1e-12));
        }
        else
        {
            rakes++;
            REQUIRE(std::abs(std::abs(end.z - e.start.z) - 4.0) < 1e-9);
            REQUIRE(e.n == Catch::Approx(1.5).epsilon(1e-12));
        }
        // The air lies above the roof, off both faces.
        REQUIRE(e.air.z > 0.0);
    }
    REQUIRE(corners == 4);
    REQUIRE(eaves == 2);
    REQUIRE(rakes == 4);
    REQUIRE(ridges == 1);
}

// Blocked where the segment passes through the house: inside its footprint
// between the floor and the roof over that point, against sampling the
// segment finely; below the floor, the house and its image, from minus the
// roof to the roof.  A point is buried under the roof and not above it.
TEST_CASE("a pitched house blocks what passes under its roof", "[pitched_roofs]")
{
    const Footprint house = gable_house();
    const Buildings city({house});
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> x(-5.0, 15.0), y(-5.0, 25.0), z(-12.0, 14.0);
    std::size_t blocked = 0, open = 0;
    for (int trial = 0; trial < 4000; trial++)
    {
        const Vec3 a(x(rng), y(rng), z(rng)), b(x(rng), y(rng), z(rng));
        for (const bool below : {false, true})
        {
            if (!below && (a.z < 0.0 || b.z < 0.0))
                continue;
            // Sampled: inside somewhere by a margin.
            bool deep = false, near = false;
            for (int k = 1; k < 4000; k++)
            {
                const Vec3 p = a + (k / 4000.0) * (b - a);
                const bool in_plan = p.x > 0.0 && p.x < 10.0 && p.y > 0.0 && p.y < 20.0;
                const double h = roof_height(house, p.x, p.y);
                const double low = below ? -h : 0.0;
                const double margin = std::min({p.x, 10.0 - p.x, p.y, 20.0 - p.y, p.z - low, h - p.z});
                if (in_plan && margin > 0.02)
                    deep = true;
                if (in_plan && margin > -0.02)
                    near = true;
            }
            if (deep)
            {
                REQUIRE(city.blocked(a, b, below));
                blocked++;
            }
            else if (!near)
            {
                REQUIRE(!city.blocked(a, b, below));
                open++;
            }
        }
    }
    REQUIRE(blocked > 500);
    REQUIRE(open > 500);
    REQUIRE(city.buried(Vec3(3.0, 10.0, 8.0)));
    REQUIRE(!city.buried(Vec3(1.0, 10.0, 7.0)));
    REQUIRE(!city.buried(Vec3(5.0, 10.0, 10.5)));
}

// A shed roof (one plane) and a hip roof (four planes): their tops, and a
// hip roof's ridge and four hips, all convex.
TEST_CASE("shed and hip roofs", "[pitched_roofs]")
{
    Footprint shed;
    shed.x = {0.0, 8.0, 8.0, 0.0};
    shed.y = {0.0, 0.0, 6.0, 6.0};
    shed.roof = {{0.5, 0.0, 3.0}};
    const Buildings sheds({shed});
    REQUIRE(sheds.footprints()[0].top == Catch::Approx(7.0).epsilon(1e-12));
    for (const Edge& e : sheds.edges())
        REQUIRE(e.n > 1.0);

    Footprint hip;
    hip.x = {0.0, 12.0, 12.0, 0.0};
    hip.y = {0.0, 0.0, 8.0, 8.0};
    // Every side's plane rising at 1 in 2 from 5 m: the ridge at 7 m from
    // x = 4 to 8.
    hip.roof = {{0.5, 0.0, 5.0}, {-0.5, 0.0, 11.0}, {0.0, 0.5, 5.0}, {0.0, -0.5, 9.0}};
    const Buildings hips({hip});
    REQUIRE(hips.footprints()[0].top == Catch::Approx(7.0).epsilon(1e-12));
    std::size_t level = 0, sloping = 0;
    for (const Edge& e : hips.edges())
    {
        if (e.vertical)
            continue;
        const bool interior = e.start.x > 1e-9 && e.start.x < 12.0 - 1e-9 && e.start.y > 1e-9 && e.start.y < 8.0 - 1e-9;
        const Vec3 end = e.start + e.length * e.direction;
        const bool interior_end = end.x > 1e-9 && end.x < 12.0 - 1e-9 && end.y > 1e-9 && end.y < 8.0 - 1e-9;
        if (!interior && !interior_end)
            continue; // an eave
        if (std::abs(e.direction.z) < 1e-12)
        {
            level++;
            REQUIRE(e.length == Catch::Approx(4.0).epsilon(1e-9));
            REQUIRE(e.start.z == Catch::Approx(7.0).epsilon(1e-12));
        }
        else
            sloping++;
        REQUIRE(e.n > 1.0);
        REQUIRE(e.n < 1.5);
    }
    REQUIRE(level == 1);
    REQUIRE(sloping == 4);
}

namespace {

// A long gable house's boundary-element solutions (tools/thick_barrier_bem.py
// --width 8 --height 6 --ridge 9, 90 elements a wavelength): frequency, the
// source's y and height, the microphone's distance behind the house and
// height, and the level, dB re the free field.
struct GableCase
{
    double f, source_y, source_z, behind, z, exact;
};

const std::vector<GableCase> gable_cases {
        {125, -25, 4.5, 2, 1.2, -10.62},
        {125, -25, 4.5, 2, 4, -9.75},
        {125, -25, 4.5, 2, 7, -3.29},
        {125, -25, 4.5, 5, 1.2, -11.88},
        {125, -25, 4.5, 5, 4, -7.12},
        {125, -25, 4.5, 5, 7, -3.50},
        {125, -25, 4.5, 10, 1.2, -21.32},
        {125, -25, 4.5, 10, 4, -8.16},
        {125, -25, 4.5, 10, 7, -3.21},
        {125, -10, 2, 2, 1.2, -15.32},
        {125, -10, 2, 2, 4, -14.43},
        {125, -10, 2, 2, 7, -8.07},
        {125, -10, 2, 5, 1.2, -16.35},
        {125, -10, 2, 5, 4, -11.66},
        {125, -10, 2, 5, 7, -8.36},
        {125, -10, 2, 10, 1.2, -25.83},
        {125, -10, 2, 10, 4, -12.66},
        {125, -10, 2, 10, 7, -8.07},
        {250, -25, 4.5, 2, 1.2, -20.96},
        {250, -25, 4.5, 2, 4, -17.49},
        {250, -25, 4.5, 2, 7, -5.80},
        {250, -25, 4.5, 5, 1.2, -19.21},
        {250, -25, 4.5, 5, 4, -11.83},
        {250, -25, 4.5, 5, 7, -5.47},
        {250, -25, 4.5, 10, 1.2, -7.31},
        {250, -25, 4.5, 10, 4, -6.73},
        {250, -25, 4.5, 10, 7, -7.55},
        {250, -10, 2, 2, 1.2, -30.52},
        {250, -10, 2, 2, 4, -27.01},
        {250, -10, 2, 2, 7, -15.49},
        {250, -10, 2, 5, 1.2, -28.53},
        {250, -10, 2, 5, 4, -21.22},
        {250, -10, 2, 5, 7, -15.13},
        {250, -10, 2, 10, 1.2, -16.36},
        {250, -10, 2, 10, 4, -15.93},
        {250, -10, 2, 10, 7, -17.39},
        {500, -25, 4.5, 2, 1.2, -26.39},
        {500, -25, 4.5, 2, 4, -18.12},
        {500, -25, 4.5, 2, 7, -9.28},
        {500, -25, 4.5, 5, 1.2, -20.61},
        {500, -25, 4.5, 5, 4, -12.55},
        {500, -25, 4.5, 5, 7, -9.38},
        {500, -25, 4.5, 10, 1.2, -10.57},
        {500, -25, 4.5, 10, 4, -11.62},
        {500, -25, 4.5, 10, 7, -10.26},
        {500, -10, 2, 2, 1.2, -29.32},
        {500, -10, 2, 2, 4, -21.01},
        {500, -10, 2, 2, 7, -12.45},
        {500, -10, 2, 5, 1.2, -23.28},
        {500, -10, 2, 5, 4, -15.23},
        {500, -10, 2, 5, 7, -12.67},
        {500, -10, 2, 10, 1.2, -12.94},
        {500, -10, 2, 10, 4, -14.27},
        {500, -10, 2, 10, 7, -13.69},
        {1000, -25, 4.5, 2, 1.2, -37.83},
        {1000, -25, 4.5, 2, 4, -28.30},
        {1000, -25, 4.5, 2, 7, -15.56},
        {1000, -25, 4.5, 5, 1.2, -31.21},
        {1000, -25, 4.5, 5, 4, -16.96},
        {1000, -25, 4.5, 5, 7, -15.74},
        {1000, -25, 4.5, 10, 1.2, -16.85},
        {1000, -25, 4.5, 10, 4, -21.18},
        {1000, -25, 4.5, 10, 7, -15.58},
        {1000, -10, 2, 2, 4, -29.76},
        {1000, -10, 2, 2, 7, -17.31},
        {1000, -10, 2, 5, 1.2, -32.46},
        {1000, -10, 2, 5, 4, -18.18},
        {1000, -10, 2, 5, 7, -17.69},
        {1000, -10, 2, 10, 1.2, -17.74},
        {1000, -10, 2, 10, 4, -22.43},
        {1000, -10, 2, 10, 7, -17.82},
};

// The level, dB re the free field, of a long house `width` across: every
// ray times the root of its length, the line source's equivalent.
double gable_level(const Buildings& city, const GableCase& c, const PointPathOptions& options, double width = 8.0)
{
    const double k = 2.0 * pi * c.f / 343.2;
    const Vec3 s(0.0, c.source_y, c.source_z), m(0.0, width + c.behind, c.z);
    const std::vector<double> reflection(city.facades().size(), 1.0);
    const PointSourceField field(city, s, m, options);
    std::complex<double> sum(0.0, 0.0);
    for (const PointRay& ray : field.rays())
    {
        std::complex<double> a;
        if (options.edge_source)
        {
            const std::vector<std::complex<double>> f = edge_source_factors(city, ray, {k});
            a = f.empty() ? point_ray_amplitude(ray, reflection, 1.0, k) : reflected_amplitude(ray, reflection, 1.0) * f[0];
        }
        else
            a = point_ray_amplitude(ray, reflection, 1.0, k);
        sum += a * std::sqrt(ray.length) * std::polar(1.0, k * ray.length);
    }
    return 20.0 * std::log10(std::abs(sum) * std::sqrt((m - s).norm()));
}

} // namespace

// A gable house 8 m across, its eaves 6 m up and its ridge 9 m, long enough
// to be two-dimensional, on rigid ground; the source and the microphone in
// one plane across it, as for the long block (test_point_paths.cpp): every
// ray times the root of its length is the line source's, and the exact
// values a boundary-element solution for one (90 elements a wavelength,
// converged to 0.15 dB; the one case that was not is left out).  Sound
// crosses the roof by the near eave, the ridge and the far eave, three edges
// whose wedges are the roof's own.  UTD is within 0.74 dB, rms 0.17, and the
// edge-source integral within 0.53, rms 0.12; a flat roof at the eaves, at
// the mean height or at the ridge -- what LoD1 would give it -- is 6.2, 8.9
// and 11.0 dB rms off.
TEST_CASE("over a long gable house the point-source field is the boundary-element solution's", "[pitched_roofs]")
{
    Footprint house;
    house.x = {-2000.0, 2000.0, 2000.0, -2000.0};
    house.y = {0.0, 0.0, 8.0, 8.0};
    house.roof = {{0.0, 0.75, 6.0}, {0.0, -0.75, 12.0}};
    const Buildings city({house});
    Footprint box = house;
    box.roof.clear();
    box.top = 7.5;
    const Buildings flat({box});
    double square[3] {0.0, 0.0, 0.0};
    for (const GableCase& c : gable_cases)
    {
        INFO(c.f << " Hz, source " << c.source_y << " m, " << c.source_z << " m up, microphone " << c.behind
                  << " m behind, " << c.z << " m up");
        PointPathOptions options;
        const double utd = gable_level(city, c, options) - c.exact;
        REQUIRE(std::abs(utd) < 1.0);
        options.edge_source = true;
        const double edge_source = gable_level(city, c, options) - c.exact;
        REQUIRE(std::abs(edge_source) < 0.7);
        const double lod1 = gable_level(flat, c, PointPathOptions()) - c.exact;
        square[0] += utd * utd;
        square[1] += edge_source * edge_source;
        square[2] += lod1 * lod1;
    }
    const double n = static_cast<double>(gable_cases.size());
    REQUIRE(std::sqrt(square[0] / n) < 0.25);
    REQUIRE(std::sqrt(square[1] / n) < 0.18);
    REQUIRE(std::sqrt(square[2] / n) > 5.0);
}

namespace {

std::string write_file(const std::string& name, const std::string& text)
{
    const auto path = (std::filesystem::temp_directory_path() / name).string();
    std::ofstream(path) << text;
    return path;
}

} // namespace

// A district in feet with a gable house, a hip-roofed house capped flat (a
// deck at 22 ft, below its planes' 22.5), a shed whose top is above its
// plane, and a flat block.  Lengths convert, slopes do not, and a plane is
// the same plane wherever its point is given on it.
TEST_CASE("roof planes come from a CSV beside the district", "[pitched_roofs]")
{
    const std::string district = write_file("niceops_roofs_district.csv",
            "building,x,y,top\n"
            "1,0,0,40\n1,30,0,40\n1,30,60,40\n1,0,60,40\n"
            "2,50,0,22\n2,90,0,22\n2,90,30,22\n2,50,30,22\n"
            "3,100,0,30\n3,120,0,30\n3,120,20,30\n3,100,20,30\n"
            "4,0,100,15\n4,10,100,15\n4,10,110,15\n4,0,110,15\n");
    // The gable rises 0.8 from 18 ft at its eaves, its second plane given by
    // a point on its far eave; the hip rises 0.5 from 15 ft on every side.
    const std::string roofs = write_file("niceops_roofs_planes.csv",
            "building,x,y,z,slope_x,slope_y\n"
            "1,0,0,18,0.8,0\n1,30,45,18,-0.8,0\n"
            "2,50,0,15,0.5,0\n2,90,0,15,-0.5,0\n2,50,0,15,0,0.5\n2,50,30,15,0,-0.5\n"
            "3,100,0,10,0,0.5\n");
    BuildingInputOptions options;
    options.roofs = roofs;
    const BuildingInput in = read_buildings(district, options);
    REQUIRE(in.footprints[0].roof.size() == 2);
    REQUIRE(roof_height(in.footprints[0], 15.0 * ft2m, 10.0 * ft2m) == Catch::Approx(30.0 * ft2m).epsilon(1e-12));
    REQUIRE(roof_height(in.footprints[0], 30.0 * ft2m, 10.0 * ft2m) == Catch::Approx(18.0 * ft2m).epsilon(1e-12));
    REQUIRE(in.footprints[0].roof[1].a == -0.8);
    // The hip's deck: a fifth, level plane at its top.
    REQUIRE(in.footprints[1].roof.size() == 5);
    REQUIRE(roof_height(in.footprints[1], 70.0 * ft2m, 15.0 * ft2m) == Catch::Approx(22.0 * ft2m).epsilon(1e-12));
    REQUIRE(roof_height(in.footprints[1], 52.0 * ft2m, 15.0 * ft2m) == Catch::Approx(16.0 * ft2m).epsilon(1e-12));
    REQUIRE(in.footprints[2].roof.size() == 1);
    REQUIRE(in.footprints[3].roof.empty());
    const Buildings city(in.footprints);
    REQUIRE(city.footprints()[0].top == Catch::Approx(30.0 * ft2m).epsilon(1e-12));
    REQUIRE(city.footprints()[1].top == Catch::Approx(22.0 * ft2m).epsilon(1e-12));
    REQUIRE(city.footprints()[2].top == Catch::Approx(20.0 * ft2m).epsilon(1e-12));
    REQUIRE(city.footprints()[3].top == Catch::Approx(15.0 * ft2m).epsilon(1e-12));

    options.roofs = write_file("niceops_roofs_stray.csv", "building,x,y,z,slope_x,slope_y\n7,0,0,10,0.5,0\n");
    REQUIRE_THROWS_WITH(read_buildings(district, options), Catch::Matchers::ContainsSubstring("building 7 is not"));
    options.roofs = write_file("niceops_roofs_metres.csv", "building,x,y,z_m,slope_x,slope_y\n1,0,0,6,0.5,0\n");
    REQUIRE_THROWS_WITH(read_buildings(district, options), Catch::Matchers::ContainsSubstring("is in m"));
}

namespace {

// A long mansard house's boundary-element solutions (tools/thick_barrier_bem.py
// --width 12 --height 6 --roof '1.5,7.5;5.25,9;6.75,9;10.5,7.5', 90 elements
// a wavelength, the four cases not converged to 0.15 dB left out), as
// gable_cases.
const std::vector<GableCase> mansard_cases {
        {125, -30, 2, 2, 0.5, -26.60},
        {125, -30, 2, 2, 1.5, -19.67},
        {125, -30, 2, 2, 4, -19.89},
        {125, -30, 2, 6, 0.5, -20.78},
        {125, -30, 2, 6, 1.5, -18.28},
        {125, -30, 2, 6, 4, -15.05},
        {125, -30, 2, 15, 0.5, -14.15},
        {125, -30, 2, 15, 1.5, -28.20},
        {125, -30, 2, 15, 4, -14.75},
        {125, -15, 0.5, 2, 0.5, -22.75},
        {125, -15, 0.5, 2, 1.5, -15.82},
        {125, -15, 0.5, 2, 4, -16.03},
        {125, -15, 0.5, 6, 0.5, -16.78},
        {125, -15, 0.5, 6, 1.5, -14.27},
        {125, -15, 0.5, 6, 4, -11.07},
        {125, -15, 0.5, 15, 0.5, -9.96},
        {125, -15, 0.5, 15, 1.5, -24.25},
        {125, -15, 0.5, 15, 4, -10.65},
        {250, -30, 2, 2, 0.5, -26.35},
        {250, -30, 2, 2, 1.5, -21.55},
        {250, -30, 2, 2, 4, -24.12},
        {250, -30, 2, 6, 0.5, -34.43},
        {250, -30, 2, 6, 1.5, -25.02},
        {250, -30, 2, 6, 4, -14.13},
        {250, -30, 2, 15, 0.5, -16.47},
        {250, -30, 2, 15, 1.5, -11.44},
        {250, -30, 2, 15, 4, -16.58},
        {250, -15, 0.5, 2, 0.5, -31.51},
        {250, -15, 0.5, 2, 1.5, -26.70},
        {250, -15, 0.5, 2, 4, -29.28},
        {250, -15, 0.5, 6, 0.5, -39.41},
        {250, -15, 0.5, 6, 1.5, -29.98},
        {250, -15, 0.5, 6, 4, -19.11},
        {250, -15, 0.5, 15, 0.5, -21.16},
        {250, -15, 0.5, 15, 1.5, -16.11},
        {250, -15, 0.5, 15, 4, -21.43},
        {500, -30, 2, 2, 1.5, -40.28},
        {500, -30, 2, 2, 4, -36.09},
        {500, -30, 2, 6, 0.5, -34.20},
        {500, -30, 2, 6, 1.5, -34.13},
        {500, -30, 2, 6, 4, -28.97},
        {500, -30, 2, 15, 0.5, -34.63},
        {500, -30, 2, 15, 1.5, -26.64},
        {500, -30, 2, 15, 4, -24.98},
        {500, -15, 0.5, 2, 0.5, -47.29},
        {500, -15, 0.5, 2, 1.5, -38.89},
        {500, -15, 0.5, 2, 4, -34.64},
        {500, -15, 0.5, 6, 0.5, -32.52},
        {500, -15, 0.5, 6, 1.5, -32.60},
        {500, -15, 0.5, 6, 4, -27.46},
        {500, -15, 0.5, 15, 0.5, -32.75},
        {500, -15, 0.5, 15, 1.5, -24.80},
        {500, -15, 0.5, 15, 4, -23.35},
        {1000, -30, 2, 2, 0.5, -36.76},
        {1000, -30, 2, 2, 1.5, -37.41},
        {1000, -30, 2, 2, 4, -31.29},
        {1000, -30, 2, 6, 0.5, -26.75},
        {1000, -30, 2, 6, 1.5, -27.66},
        {1000, -30, 2, 6, 4, -20.27},
        {1000, -30, 2, 15, 0.5, -18.75},
        {1000, -30, 2, 15, 1.5, -21.23},
        {1000, -30, 2, 15, 4, -15.78},
        {1000, -15, 0.5, 6, 0.5, -34.31},
        {1000, -15, 0.5, 6, 1.5, -35.06},
        {1000, -15, 0.5, 6, 4, -27.72},
        {1000, -15, 0.5, 15, 0.5, -25.93},
        {1000, -15, 0.5, 15, 1.5, -28.41},
        {1000, -15, 0.5, 15, 4, -22.97},
};

} // namespace

// A long mansard house 12 m across, its eaves 6 m up and its deck 9 m (slopes
// of 1 and 0.4 each side), against its boundary-element solution as the
// gable is.  From a source low beside it to a microphone low and close
// behind, the path bends at four edges -- the deck's two, the far break and
// the far eave -- and only a chain (PointPathOptions::max_chain) carries it:
// without one there is no ray at all.  With it those levels are within 1.5
// dB below 4 m (at 4 m the microphone sits on the far lower slope's plane,
// at its eave's boundary).  Over every case the rays are 3.9 dB rms: the
// cascade over a mansard's nearly flat wedges (n 1.1), each pair sharing the
// face the ray runs along, is that far from the exact field with three edges
// as with four, where a gable's is 0.17 (THEORY 6.5).
TEST_CASE("over a long mansard house the rays need four edges", "[pitched_roofs]")
{
    Footprint house;
    house.x = {-2000.0, 2000.0, 2000.0, -2000.0};
    house.y = {0.0, 0.0, 12.0, 12.0};
    house.top = 9.0;
    house.roof = {{0.0, 1.0, 6.0}, {0.0, 0.4, 6.9}, {0.0, -0.4, 11.7}, {0.0, -1.0, 18.0}};
    const Buildings city({house});
    PointPathOptions none;
    none.max_chain = 3;
    double square = 0.0;
    int chained = 0;
    for (const GableCase& c : mansard_cases)
    {
        INFO(c.f << " Hz, source " << c.source_y << " m, " << c.source_z << " m up, microphone " << c.behind
                  << " m behind, " << c.z << " m up");
        const double error = gable_level(city, c, PointPathOptions(), 12.0) - c.exact;
        square += error * error;
        const Vec3 s(0.0, c.source_y, c.source_z), m(0.0, 12.0 + c.behind, c.z);
        const PointSourceField field(city, s, m, PointPathOptions());
        int most = 0;
        for (const PointRay& ray : field.rays())
            most = std::max(most, ray.wedges);
        const bool low_and_close = c.source_z < 1.0 && c.behind < 3.0;
        if (low_and_close)
        {
            REQUIRE(most == 4);
            REQUIRE(PointSourceField(city, s, m, none).rays().empty());
            if (c.z < 2.0)
                REQUIRE(std::abs(error) < 1.5);
            chained++;
        }
        else if (c.source_z > 1.0)
            REQUIRE(most == 3);
    }
    REQUIRE(chained == 9);
    REQUIRE(std::sqrt(square / mansard_cases.size()) < 4.2);
}

// Through three edges a chain's ray is the thrice-diffracted one, term for
// term: here over a hip roof's eave, hip and ridge, none parallel, at points
// that need not be stationary.
TEST_CASE("a chain of three edges is the thrice-diffracted ray", "[pitched_roofs]")
{
    Footprint hip;
    hip.x = {0.0, 12.0, 12.0, 0.0};
    hip.y = {0.0, 0.0, 8.0, 8.0};
    hip.roof = {{0.5, 0.0, 5.0}, {-0.5, 0.0, 11.0}, {0.0, 0.5, 5.0}, {0.0, -0.5, 9.0}};
    const Buildings city({hip});
    std::vector<const Edge*> g;
    for (const Edge& e : city.edges())
        if (!e.vertical && g.size() < 3
                && (g.empty() || std::abs(std::abs(e.direction.dot(g.back()->direction)) - 1.0) > 1e-6))
            g.push_back(&e);
    REQUIRE(g.size() == 3);
    const std::vector<Vec3> q {g[0]->start + 0.4 * g[0]->length * g[0]->direction,
            g[1]->start + 0.6 * g[1]->length * g[1]->direction, g[2]->start + 0.5 * g[2]->length * g[2]->direction};
    const Vec3 pa(-8.0, -3.0, 1.5), pb(19.0, 12.0, 2.0);
    const PointRay chain = chain_diffracted_ray(g, pa, pb, q);
    const PointRay thrice = thrice_diffracted_ray(*g[0], *g[1], *g[2], pa, pb, q[0], q[1], q[2]);
    REQUIRE(chain.wedges == 3);
    REQUIRE(chain.length == Catch::Approx(thrice.length).epsilon(1e-14));
    REQUIRE(chain.spreading == Catch::Approx(thrice.spreading).epsilon(1e-12));
    for (int w = 0; w < 3; w++)
    {
        const WedgeTerm &a = chain.wedge[w], &b = thrice.wedge[w];
        REQUIRE(a.angles.n == b.angles.n);
        for (int k = 0; k < 4; k++)
        {
            REQUIRE(a.angles.c[k] == Catch::Approx(b.angles.c[k]).margin(1e-12));
            REQUIRE(a.angles.cot[k] == Catch::Approx(b.angles.cot[k]).epsilon(1e-10).margin(1e-10));
        }
        REQUIRE(a.sin_beta == Catch::Approx(b.sin_beta).epsilon(1e-12));
        REQUIRE(a.distance == Catch::Approx(b.distance).epsilon(1e-12));
        REQUIRE(a.grazing == b.grazing);
    }
}

// Across a hip roof below its ridge, from low beside it: the far eave, two
// hips and the near eave, the hips oblique to the eaves, and the path over
// them stationary in 3-D -- Keller's law at every edge -- however the solver
// starts.
TEST_CASE("a chain over a hip roof keeps Keller's law at every edge", "[pitched_roofs]")
{
    Footprint hip;
    hip.x = {0.0, 12.0, 12.0, 0.0};
    hip.y = {0.0, 0.0, 8.0, 8.0};
    hip.roof = {{0.5, 0.0, 5.0}, {-0.5, 0.0, 11.0}, {0.0, 0.5, 5.0}, {0.0, -0.5, 9.0}};
    const Buildings city({hip});
    const Vec3 s(-2.0, 2.6, 1.0), m(14.0, 3.4, 1.0);
    const PointSourceField field(city, s, m, PointPathOptions());
    const PointRay* chain = nullptr;
    for (const PointRay& ray : field.rays())
        if (ray.wedges == 4 && ray.floor == 0)
            chain = &ray;
    REQUIRE(chain != nullptr);
    REQUIRE(chain->path.size() == 6);
    std::vector<const Edge*> g;
    int oblique = 0;
    for (std::size_t i = 1; i <= 4; i++)
    {
        const Vec3& p = chain->path[i];
        for (const Edge& e : city.edges())
        {
            const double t = (p - e.start).dot(e.direction);
            if (!e.vertical && t > 0.0 && t < e.length && (e.start + t * e.direction - p).norm() < 1e-9)
                g.push_back(&e);
        }
        REQUIRE(g.size() == i);
        oblique += std::abs(g.back()->direction.z) > 0.1;
        const Vec3 in = (p - chain->path[i - 1]) / (p - chain->path[i - 1]).norm();
        const Vec3 out = (chain->path[i + 1] - p) / (chain->path[i + 1] - p).norm();
        REQUIRE(std::abs(in.dot(g.back()->direction) - out.dot(g.back()->direction)) < 1e-9);
    }
    REQUIRE(oblique == 2);
    // From each edge's middle, the solver finds the same points.
    std::vector<double> t;
    for (const Edge* e : g)
        t.push_back(0.5 * e->length);
    REQUIRE(chain_points(g, s, m, t));
    for (std::size_t i = 0; i < 4; i++)
        REQUIRE((g[i]->start + t[i] * g[i]->direction - chain->path[i + 1]).norm() < 1e-8);
}

// The pruning's visibility proof takes a convex pitched house as its prism cut
// by each roof plane (mirrored in a facade and in the floor as the footprint
// is): pruned, the search finds the rays the unpruned one does, ray for ray,
// among gable, hip, shed and mansard houses, from low and from above, the
// microphones in their shadows and in the street.
TEST_CASE("pruning finds every ray among pitched houses", "[pitched_roofs]")
{
    auto house = [](double x0, double y0, double w, double d, std::vector<RoofPlane> roof, double top) {
        Footprint f;
        f.x = {x0, x0 + w, x0 + w, x0};
        f.y = {y0, y0, y0 + d, y0 + d};
        // The planes given over the house's own corner.
        for (RoofPlane& p : roof)
            p.c -= p.a * x0 + p.b * y0;
        f.roof = std::move(roof);
        f.top = top;
        return f;
    };
    const Buildings city({house(0.0, 0.0, 10.0, 16.0, {{0.7, 0.0, 6.0}, {-0.7, 0.0, 13.0}}, 20.0),
            house(18.0, 0.0, 12.0, 8.0, {{0.5, 0.0, 5.0}, {-0.5, 0.0, 11.0}, {0.0, 0.5, 5.0}, {0.0, -0.5, 9.0}}, 20.0),
            house(0.0, 24.0, 8.0, 6.0, {{0.0, 0.4, 4.0}}, 20.0),
            house(18.0, 18.0, 12.0, 14.0, {{1.0, 0.0, 6.0}, {0.4, 0.0, 6.9}, {-0.4, 0.0, 11.7}, {-1.0, 0.0, 18.0}}, 9.0)});
    const std::vector<std::pair<Vec3, Vec3>> pairs {{{-12.0, 8.0, 1.0}, {14.0, 8.0, 1.2}},
            {{-12.0, 8.0, 1.0}, {34.0, 25.0, 1.2}}, {{24.0, -10.0, 1.5}, {24.0, 12.0, 1.2}},
            {{14.0, 40.0, 30.0}, {14.0, -6.0, 1.2}}, {{-6.0, 27.0, 2.0}, {36.0, 4.0, 1.2}},
            {{12.0, 27.0, 1.2}, {34.0, 25.0, 1.2}}};
    std::size_t diffracted = 0;
    for (const auto& [s, m] : pairs)
    {
        INFO("source " << s.x << ", " << s.y << ", " << s.z << "; microphone " << m.x << ", " << m.y << ", " << m.z);
        PointPathOptions all;
        all.prune = false;
        const PointSourceField pruned(city, s, m);
        const PointSourceField every(city, s, m, all);
        REQUIRE(pruned.rays().size() == every.rays().size());
        for (std::size_t i = 0; i < every.rays().size(); i++)
        {
            REQUIRE(pruned.rays()[i].length == every.rays()[i].length);
            REQUIRE(pruned.rays()[i].wedges == every.rays()[i].wedges);
            diffracted += every.rays()[i].wedges >= 2;
        }
    }
    REQUIRE(diffracted > 100);
}
