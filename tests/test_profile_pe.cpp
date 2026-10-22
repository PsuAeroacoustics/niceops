#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "buildings.hpp"
#include "ground_effect.hpp"
#include "point_paths.hpp"
#include "profile_pe.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT

namespace {

// A boundary-element solution's value (tools/thick_barrier_bem.py --blocks):
// frequency, the source's u and height, the receiver's u and height, and the
// level, dB re the free field.
struct ProfileCase
{
    double f, su, sz, ru, rz, exact;
};

// Two blocks 10 m high, a 10 m street between (45 elements a wavelength).
const std::vector<ProfileCase>& street_cases()
{
    static const std::vector<ProfileCase> cases {
        {125, -100, 12, 32, 1.2, -11.66}, {125, -100, 12, 35, 1.2, -12.56}, {125, -100, 12, 40, 1.2, -17.97},
        {125, -50, 15, 32, 1.2, -11.75}, {125, -50, 15, 35, 1.2, -12.58}, {125, -50, 15, 40, 1.2, -17.84},
        {250, -100, 12, 32, 1.2, -20.04}, {250, -100, 12, 35, 1.2, -29.22}, {250, -100, 12, 40, 1.2, -17.17},
        {250, -50, 15, 32, 1.2, -18.81}, {250, -50, 15, 35, 1.2, -27.91}, {250, -50, 15, 40, 1.2, -15.75},
        {500, -100, 12, 32, 1.2, -34.19}, {500, -100, 12, 35, 1.2, -21.78}, {500, -100, 12, 40, 1.2, -38.06},
        {500, -50, 15, 32, 1.2, -28.44}, {500, -50, 15, 35, 1.2, -15.98}, {500, -50, 15, 40, 1.2, -32.14},
    };
    return cases;
}

// One row's cases against the profile field: the worst error and the rms.
void score(const std::vector<ProfileBlock>& blocks, const std::vector<ProfileCase>& cases, double worst_bound,
        double rms_bound)
{
    double worst = 0.0, square = 0.0;
    for (const ProfileCase& c : cases)
    {
        ProfileProblem p;
        p.blocks = blocks;
        p.su = c.su;
        p.sz = c.sz;
        p.ru = c.ru;
        p.rz = c.rz;
        const double k = 2.0 * pi * c.f / 343.2;
        const double r = std::hypot(c.ru - c.su, c.rz - c.sz);
        const double error = 20.0 * std::log10(std::abs(profile_field(p, k)) * std::sqrt(r)) - c.exact;
        INFO(c.f << " Hz, source (" << c.su << ", " << c.sz << "), receiver (" << c.ru << ", " << c.rz << ")");
        REQUIRE(std::abs(error) < worst_bound);
        worst = std::max(worst, std::abs(error));
        square += error * error;
    }
    INFO("worst " << worst << " dB");
    REQUIRE(std::sqrt(square / static_cast<double>(cases.size())) < rms_bound);
}

} // namespace

// Rows of flat blocks 10 m across, on rigid ground, line sources low and far
// before them, receivers 1.2 m up in the shadow behind (FIDELITY, "Over a row
// of roofs").  The march drops the street's reverberation, which leaves it
// 1-2.5 dB low at 125 Hz and less above; the UTD cascades it replaces in the
// point model are 3.4 dB rms over one street and 6.7 over two.
TEST_CASE("the profile field over a street is the boundary-element solution's", "[profile_pe]")
{
    // Two blocks 10 m high, a 10 m street between (45 elements a wavelength).
    score({{0, 10, 10}, {20, 30, 10}}, street_cases(), 2.6, 1.7);
    // The same at 1 kHz (20 elements a wavelength), where the cascades are
    // still 3.5 dB rms out.
    const std::vector<ProfileCase> kilohertz {
            {1000, -100, 12, 32, 1.2, -23.56}, {1000, -100, 12, 35, 1.2, -25.09}, {1000, -100, 12, 40, 1.2, -21.50},
            {1000, -50, 15, 32, 1.2, -18.25}, {1000, -50, 15, 35, 1.2, -19.76}, {1000, -50, 15, 40, 1.2, -15.96},
    };
    score({{0, 10, 10}, {20, 30, 10}}, kilohertz, 0.6, 0.35);
    // A source higher, 27 degrees over the first roof (30 elements).
    const std::vector<ProfileCase> steep {
            {125, -50, 25, 32, 1.2, -8.12}, {125, -50, 25, 35, 1.2, -8.89}, {125, -50, 25, 40, 1.2, -13.98},
            {250, -50, 25, 32, 1.2, -15.03}, {250, -50, 25, 35, 1.2, -23.99}, {250, -50, 25, 40, 1.2, -11.80},
            {500, -50, 25, 32, 1.2, -25.57}, {500, -50, 25, 35, 1.2, -13.11}, {500, -50, 25, 40, 1.2, -29.08},
    };
    score({{0, 10, 10}, {20, 30, 10}}, steep, 1.8, 1.2);
}

TEST_CASE("the profile field over a step and over two streets is the boundary-element solution's", "[profile_pe]")
{
    // The second block 8 m high (45 elements a wavelength).
    const std::vector<ProfileCase> step {
            {125, -100, 12, 32, 1.2, -13.87}, {125, -100, 12, 35, 1.2, -15.43}, {125, -100, 12, 40, 1.2, -25.12},
            {125, -50, 15, 32, 1.2, -13.55}, {125, -50, 15, 35, 1.2, -15.02}, {125, -50, 15, 40, 1.2, -24.46},
            {250, -100, 12, 32, 1.2, -24.56}, {250, -100, 12, 35, 1.2, -38.42}, {250, -100, 12, 40, 1.2, -18.18},
            {250, -50, 15, 32, 1.2, -22.38}, {250, -50, 15, 35, 1.2, -36.14}, {250, -50, 15, 40, 1.2, -15.79},
            {500, -100, 12, 32, 1.2, -36.67}, {500, -100, 12, 35, 1.2, -26.86}, {500, -100, 12, 40, 1.2, -26.87},
            {500, -50, 15, 32, 1.2, -29.49}, {500, -50, 15, 35, 1.2, -19.61}, {500, -50, 15, 40, 1.2, -19.50},
    };
    score({{0, 10, 10}, {20, 30, 8}}, step, 1.7, 1.2);
    // Three blocks, two streets (30 elements a wavelength).
    const std::vector<ProfileCase> two_streets {
            {125, -100, 12, 52, 1.2, -12.39}, {125, -100, 12, 55, 1.2, -13.28}, {125, -100, 12, 60, 1.2, -18.74},
            {125, -50, 15, 52, 1.2, -12.02}, {125, -50, 15, 55, 1.2, -12.87}, {125, -50, 15, 60, 1.2, -18.23},
            {125, -50, 25, 52, 1.2, -7.84}, {125, -50, 25, 55, 1.2, -8.67}, {125, -50, 25, 60, 1.2, -13.90},
            {250, -100, 12, 52, 1.2, -21.09}, {250, -100, 12, 55, 1.2, -30.34}, {250, -100, 12, 60, 1.2, -18.26},
            {250, -50, 15, 52, 1.2, -19.26}, {250, -50, 15, 55, 1.2, -28.46}, {250, -50, 15, 60, 1.2, -16.31},
            {250, -50, 25, 52, 1.2, -16.40}, {250, -50, 25, 55, 1.2, -25.49}, {250, -50, 25, 60, 1.2, -13.30},
            {500, -100, 12, 52, 1.2, -34.37}, {500, -100, 12, 55, 1.2, -22.01}, {500, -100, 12, 60, 1.2, -38.59},
            {500, -50, 15, 52, 1.2, -28.35}, {500, -50, 15, 55, 1.2, -15.95}, {500, -50, 15, 60, 1.2, -32.43},
            {500, -50, 25, 52, 1.2, -26.29}, {500, -50, 25, 55, 1.2, -13.95}, {500, -50, 25, 60, 1.2, -30.10},
    };
    score({{0, 10, 10}, {20, 30, 10}, {40, 50, 10}}, two_streets, 3.7, 2.2);
}

// One block: the march's start and finish alone, against the same BEM (30
// elements a wavelength).  The start drops what meets the front face below
// the roof, and the finish takes the back corner by UTD.
TEST_CASE("the profile field over one block is the boundary-element solution's", "[profile_pe]")
{
    const std::vector<ProfileCase> block {
            {125, -100, 12, 12, 1.2, -10.39}, {125, -100, 12, 15, 1.2, -11.26}, {125, -100, 12, 20, 1.2, -16.69},
            {125, -50, 15, 12, 1.2, -11.26}, {125, -50, 15, 15, 1.2, -11.99}, {125, -50, 15, 20, 1.2, -17.13},
            {250, -100, 12, 12, 1.2, -18.36}, {250, -100, 12, 15, 1.2, -27.53}, {250, -100, 12, 20, 1.2, -15.47},
            {250, -50, 15, 12, 1.2, -18.24}, {250, -50, 15, 15, 1.2, -27.32}, {250, -50, 15, 20, 1.2, -15.02},
            {500, -100, 12, 12, 1.2, -35.06}, {500, -100, 12, 15, 1.2, -22.46}, {500, -100, 12, 20, 1.2, -38.66},
            {500, -50, 15, 12, 1.2, -29.67}, {500, -50, 15, 15, 1.2, -17.11}, {500, -50, 15, 20, 1.2, -33.07},
    };
    score({{0, 10, 10}}, block, 0.7, 0.45);
}

TEST_CASE("the profile field refuses a row it does not take", "[profile_pe]")
{
    ProfileProblem p;
    p.blocks = {{0, 10, 10}, {20, 30, 10}};
    p.su = -50;
    p.sz = 15;
    p.ru = 35;
    p.rz = 1.2;
    REQUIRE(profile_applies(p));
    // Touching blocks are a row; overlapping ones, a source over a roof, a
    // receiver before the last block's back or below the floor are not.
    p.blocks[1].u0 = 10;
    REQUIRE(profile_applies(p));
    p.blocks[1].u0 = 9;
    REQUIRE_FALSE(profile_applies(p));
    p.blocks[1].u0 = 20;
    p.su = 5;
    REQUIRE_FALSE(profile_applies(p));
    p.su = -50;
    p.ru = 29;
    REQUIRE_FALSE(profile_applies(p));
    p.ru = 35;
    p.rz = -0.1;
    REQUIRE_FALSE(profile_applies(p));
    REQUIRE_THROWS_AS(profile_field(p, 2.0), std::invalid_argument);
    p.rz = 1.2;
    p.blocks.clear();
    REQUIRE_FALSE(profile_applies(p));
}

// The finish's corner Green's function steps at every Huygens point's shadow
// boundary, which its diffraction cancels; the field behind the row is
// smooth as the receiver crosses them, every 5 cm up through the roof plane
// and above.
TEST_CASE("the profile field is smooth across the corner's shadow boundaries", "[profile_pe]")
{
    ProfileProblem p;
    p.blocks = {{0, 10, 10}, {20, 30, 10}};
    p.su = -50;
    p.sz = 12;
    p.ru = 34;
    const double k = 2.0 * pi * 250.0 / 343.2;
    std::vector<std::complex<double>> field;
    for (int i = 0; i <= 300; i++)
    {
        p.rz = 1.0 + 0.05 * i;
        field.push_back(profile_field(p, k));
    }
    double scale = 0.0;
    for (const auto& f : field)
        scale = std::max(scale, std::abs(f));
    double worst = 0.0;
    for (std::size_t i = 1; i + 1 < field.size(); i++)
        worst = std::max(worst, std::abs(field[i + 1] - 2.0 * field[i] + field[i - 1]) / scale);
    REQUIRE(worst < 0.01);
}

// Each part goes as exp(i k l), l its own path's length (profile_path),
// which is what lets the point model carry it as a ray of that length: its
// phase's slope in k is l, to a fraction of a wavelength.  A source high and
// far, one low and near, and a receiver higher behind; 250 Hz less 1% puts a
// Huygens node exactly on the corner's shadow boundary from the second.
TEST_CASE("each part of the profile field goes as its path's length", "[profile_pe]")
{
    const std::vector<std::array<double, 3>> ends {{12, 34, 1.2}, {2, 36, 1.2}, {15, 40, 4}};
    for (const auto& [sz, ru, rz] : ends)
        for (const double f : {250.0, 500.0, 1000.0})
        {
            ProfileProblem p;
            p.blocks = {{0, 10, 10}, {20, 30, 10}};
            p.su = -60;
            p.sz = sz;
            p.ru = ru;
            p.rz = rz;
            const double k = 2.0 * pi * f / 343.2, dk = 0.01 * k;
            const ProfileParts lo = profile_parts(p, k - dk), hi = profile_parts(p, k + dk);
            for (int a = 0; a < 2; a++)
                for (int b = 0; b < 2; b++)
                {
                    const double l = profile_path(p, a, b);
                    // The phase's change, unwrapped by the length's.
                    const double change = std::arg(hi.part[a][b] / lo.part[a][b]);
                    const double turns = std::round((l * 2.0 * dk - change) / (2.0 * pi));
                    const double delay = (change + 2.0 * pi * turns) / (2.0 * dk);
                    INFO(f << " Hz, source " << sz << " m up, receiver (" << ru << ", " << rz << "), part " << a << b
                           << ": delay " << delay << " m, path " << l << " m");
                    REQUIRE(std::abs(delay - l) < 0.3);
                    // Nowhere near the 30 dB a node on a shadow boundary gave.
                    REQUIRE(std::abs(20.0 * std::log10(std::abs(hi.part[a][b] / lo.part[a][b]))) < 3.0);
                }
        }
}

TEST_CASE("a part's path is the taut string over the corners", "[profile_pe]")
{
    ProfileProblem p;
    p.blocks = {{0, 10, 10}, {20, 30, 8}};
    p.su = -60;
    p.sz = 2;
    p.ru = 36;
    p.rz = 1.2;
    std::vector<std::array<double, 2>> points;
    // Over the first roof's two corners and the second's far one; the
    // second's near corner lies under the string.
    const double l = profile_path(p, false, false, &points);
    REQUIRE(points.size() == 5);
    REQUIRE(points[1] == std::array<double, 2> {0, 10});
    REQUIRE(points[2] == std::array<double, 2> {10, 10});
    REQUIRE(points[3] == std::array<double, 2> {30, 8});
    REQUIRE(l == Catch::Approx(std::hypot(60, 8) + 10 + std::hypot(20, 2) + std::hypot(6, 6.8)));
    // From the images, below the floor.
    REQUIRE(profile_path(p, true, true) == Catch::Approx(std::hypot(60, 12) + 10 + std::hypot(20, 2) + std::hypot(6, 9.2)));
}

namespace {

Footprint box(double x0, double y0, double x1, double y1, double top)
{
    Footprint f;
    f.x = {x0, x1, x1, x0};
    f.y = {y0, y0, y1, y1};
    f.top = top;
    return f;
}

// The point field's level re free field at wavenumber k, rigid facades and
// floor, and what its rays are.
struct PointLevel
{
    double level = 0.0;
    std::complex<double> field;
    std::size_t profiles = 0, rays = 0;
    double share = 0.0; // the profile ray's
};

PointLevel point_level(const Buildings& city, const Vec3& s, const Vec3& m, double k, bool profile, bool line = false)
{
    PointPathOptions options;
    options.profile_pe = profile;
    const PointSourceField field(city, s, m, options);
    const std::vector<double> reflection(city.facades().size(), 1.0);
    PointLevel out;
    for (const PointRay& ray : field.rays())
    {
        // As a line source's, each ray times the root of its length (the
        // long-block BEM test in test_point_paths.cpp).
        out.field += point_ray_amplitude(ray, reflection, 1.0, k) * std::polar(1.0, k * ray.length)
                * (line ? std::sqrt(ray.length) : 1.0);
        out.profiles += ray.profile;
        out.rays++;
        if (ray.profile)
            out.share = ray.share;
    }
    const double free = line ? 1.0 / std::sqrt((m - s).norm()) : 1.0 / (m - s).norm();
    out.level = 20.0 * std::log10(std::abs(out.field) / free);
    return out;
}

} // namespace

// Two long blocks and the street between, the source and the microphone in
// one plane across them: the point model's profile ray reproduces the 2-D
// field, as BEM has it, where the cascades it replaces are 3.4 dB rms out.
TEST_CASE("over a street the point-source field takes the profile field", "[profile_pe][point_paths]")
{
    const Buildings city({box(0, -4000, 10, 4000, 10), box(20, -4000, 30, 4000, 10)});
    double square[2] {0.0, 0.0};
    for (const ProfileCase& c : street_cases())
    {
        const double k = 2.0 * pi * c.f / 343.2;
        const Vec3 s(c.su, 0.0, c.sz), m(c.ru, 0.0, c.rz);
        for (int profile = 0; profile < 2; profile++)
        {
            const PointLevel l = point_level(city, s, m, k, profile == 1, true);
            const double error = l.level - c.exact;
            square[profile] += error * error;
            INFO(c.f << " Hz, source (" << c.su << ", " << c.sz << "), microphone " << c.ru << " m");
            if (profile == 1)
            {
                // The row's rays, the street's walls' among them, all went
                // into its four parts'.
                REQUIRE(l.profiles == 4);
                REQUIRE(l.rays == 4);
                REQUIRE(std::abs(error) < 2.6);
            }
        }
    }
    const double n = static_cast<double>(street_cases().size());
    REQUIRE(std::sqrt(square[1] / n) < 1.7);
    REQUIRE(std::sqrt(square[0] / n) > 3.0);
}

// Swapped ends march the other way round and weigh the two the other way:
// the same field, near grazing at both ends (both marches weighed), obliquely
// across the row, and with a street's floor and a facade outside the row.
TEST_CASE("the profile ray is reciprocal", "[profile_pe][point_paths]")
{
    const Buildings city({box(0, -400, 10, 400, 10), box(20, -400, 30, 400, 10), box(45, 30, 60, 50, 14)});
    const std::vector<std::pair<Vec3, Vec3>> pairs {{Vec3(-20, 0, 12), Vec3(40, 0, 11.5)},
            {Vec3(-60, -20, 12), Vec3(40, 25, 1.2)}, {Vec3(-30, 5, 3), Vec3(42, 35, 1.5)}};
    std::size_t profiles = 0;
    for (const auto& [s, m] : pairs)
        for (const double f : {125.0, 500.0})
        {
            const double k = 2.0 * pi * f / 343.2;
            const PointLevel forward = point_level(city, s, m, k, true), back = point_level(city, m, s, k, true);
            INFO(f << " Hz, source (" << s.x << ", " << s.y << ", " << s.z << ")");
            REQUIRE(forward.profiles == back.profiles);
            REQUIRE(forward.level == Catch::Approx(back.level).margin(1e-6));
            profiles += forward.profiles;
        }
    REQUIRE(profiles == 24);
}

// The profile ray is weighed in and the rays it replaces out by the ends'
// elevations, the ends' angles off the vertical over the outer corners and
// the line's distance from the blocks' ends: smooth, where the cascades step
// 3.9 dB with the microphone in the roofs' plane (the chains of three and
// four edges over the row; a twice-diffracted ray over a roof keeps its level
// through its plane, double_wedge.cpp).  Each sweep's largest second
// difference against its largest field, as in the 2-D test, with and without.
TEST_CASE("the profile ray comes in and out smoothly", "[profile_pe][point_paths]")
{
    const double k = 2.0 * pi * 250.0 / 343.2;
    struct Sweep
    {
        const char* name;
        double length; // the blocks' half-length
        Vec3 s0, s1, m0, m1;
    };
    const std::vector<Sweep> sweeps {
            // Up through the roofs' plane 2 m behind the row, and over the
            // last roof's corner.
            {"height", 4000, {-100, 0, 12}, {-100, 0, 12}, {32, 0, 1.2}, {32, 0, 19.2}},
            {"over the corner", 4000, {-100, 0, 12}, {-100, 0, 12}, {38, 0, 12.5}, {28, 0, 12.5}},
            // The source's elevation through 25 to 40 degrees.
            {"elevation", 4000, {-50, 0, 26}, {-50, 0, 62}, {35, 0, 1.2}, {35, 0, 1.2}},
            // The line out past the blocks' ends.
            {"ends", 60, {-60, 0, 12}, {-60, 30, 12}, {40, 30, 1.2}, {40, 60, 1.2}},
            // The microphone's elevation through the source's, the march's
            // two ways round.
            {"lean", 4000, {-20, 0, 12}, {-20, 0, 12}, {40, 0, 10.2}, {40, 0, 13.2}},
    };
    for (const Sweep& w : sweeps)
    {
        const Buildings city({box(0, -w.length, 10, w.length, 10), box(20, -w.length, 30, w.length, 10)});
        double worst[2] {0.0, 0.0};
        std::size_t profiles = 0;
        std::vector<double> share;
        for (int profile = 0; profile < 2; profile++)
        {
            std::vector<std::complex<double>> field;
            const int n = 240;
            for (int i = 0; i <= n; i++)
            {
                const double t = static_cast<double>(i) / n;
                const Vec3 s = w.s0 + (w.s1 - w.s0) * t, m = w.m0 + (w.m1 - w.m0) * t;
                const PointLevel l = point_level(city, s, m, k, profile == 1);
                field.push_back(l.field * (m - s).norm());
                profiles += l.profiles;
                if (profile == 1)
                    share.push_back(l.share);
            }
            double scale = 0.0;
            for (const auto& f : field)
                scale = std::max(scale, std::abs(f));
            for (std::size_t i = 1; i + 1 < field.size(); i++)
                worst[profile] = std::max(worst[profile], std::abs(field[i + 1] - 2.0 * field[i] + field[i - 1]) / scale);
        }
        // The weight itself, which changes in every sweep but the last, has
        // no step.
        double weight = 0.0;
        for (std::size_t i = 1; i + 1 < share.size(); i++)
            weight = std::max(weight, std::abs(share[i + 1] - 2.0 * share[i] + share[i - 1]));
        INFO(w.name << ": cascades " << worst[0] << ", profile " << worst[1] << ", weight " << weight);
        REQUIRE(profiles > 0);
        REQUIRE(weight < 0.05);
        // Nothing the cascades do not do already...
        REQUIRE(worst[1] < std::max(1.05 * worst[0], 0.05));
        // ...and their step in the roofs' plane gone.
        if (std::string(w.name) == "height")
        {
            REQUIRE(worst[0] > 0.08);
            REQUIRE(worst[1] < 0.02);
        }
    }
}

// A row's field is four rays, its parts, each with its own length, so that a
// band's sum averages their interference as it does the traced rays off the
// floor: its band levels are the field's energy averaged over each sub-band
// (sampled finely).  One ray carrying the whole field at each sub-band's
// middle took the interference of the ends' images there, 3.3 dB high at
// 2 kHz and up to 10 dB in a sub-band.
TEST_CASE("a row's band levels average its parts' interference", "[profile_pe][point_paths]")
{
    const Buildings city({box(0, -4000, 10, 4000, 10), box(20, -4000, 30, 4000, 10)});
    const PointSourceField field(city, Vec3(-100, 0, 12), Vec3(34, 0, 1.2));
    const std::vector<double> bands {1000.0, 2000.0};
    const double c = 343.2;
    const std::vector<double> reflection(city.facades().size() * bands.size(), 1.0);
    std::vector<double> levels(bands.size());
    field.band_levels(reflection, 1.0, bands, c, levels.data());
    const std::vector<double> rigid(city.facades().size(), 1.0);
    auto at = [&](double k) {
        std::complex<double> sum(0.0, 0.0);
        for (const PointRay& ray : field.rays())
            sum += point_ray_amplitude(ray, rigid, 1.0, k) * std::polar(1.0, k * ray.length);
        return sum;
    };
    REQUIRE(field.rays().size() == 4);
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        const auto parts = sub_bands(bands[b], c);
        double sampled[3], middle[3];
        for (int part = 0; part < 3; part++)
        {
            const int n = 16;
            sampled[part] = 0.0;
            for (int i = 0; i < n; i++)
                sampled[part] += std::norm(at(parts[part].k1 + (i + 0.5) / n * (parts[part].k2 - parts[part].k1))) / n;
            middle[part] = std::norm(at(0.5 * (parts[part].k1 + parts[part].k2)));
        }
        const double exact = band_level(bands[b], parts, sampled), frozen = band_level(bands[b], parts, middle);
        INFO(bands[b] << " Hz: band_levels " << levels[b] << ", sampled " << exact << ", at the middles " << frozen);
        REQUIRE(std::abs(levels[b] - exact) < 0.1);
        if (bands[b] == 2000.0)
            REQUIRE(frozen - exact > 2.0);
    }
}

// The parts' bounces off the floor at the ends are a ray's: a floor of
// finite impedance takes them by its spherical-wave coefficient, which as
// its admittance vanishes is the rigid floor's.
TEST_CASE("a row's profile rays take the floor's model at their bounces", "[profile_pe][point_paths]")
{
    const Buildings city({box(0, -4000, 10, 4000, 10), box(20, -4000, 30, 4000, 10)});
    const PointSourceField field(city, Vec3(-100, 0, 12), Vec3(34, 0, 1.2));
    const std::vector<double> bands {250.0, 500.0};
    const double c = 343.2;
    const std::vector<double> reflection(city.facades().size() * bands.size(), 1.0);
    std::vector<double> rigid(bands.size()), hard(bands.size()), grass(bands.size());
    field.band_levels(reflection, 1.0, bands, c, rigid.data());
    FloorModel floor;
    floor.admittance.assign(3 * bands.size(), std::complex<double>(1e-9, -1e-9));
    field.band_levels(reflection, floor, bands, c, hard.data());
    floor.admittance = ground_impedance(GroundModel {}, sub_band_middles(bands, c), c);
    field.band_levels(reflection, floor, bands, c, grass.data());
    int bounces = 0;
    for (const PointRay& ray : field.rays())
    {
        REQUIRE(ray.floor == ray.profile_term.source_image + ray.profile_term.microphone_image);
        for (int b = 0; b < ray.floor; b++)
        {
            REQUIRE(ray.floor_bounce[b].sin_grazing > 0.0);
            REQUIRE(ray.floor_bounce[b].leg > 0.0);
            bounces++;
        }
    }
    REQUIRE(bounces == 4);
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        INFO(bands[b] << " Hz: rigid " << rigid[b] << ", hard " << hard[b] << ", grass " << grass[b]);
        REQUIRE(std::abs(hard[b] - rigid[b]) < 0.01);
        REQUIRE(std::abs(grass[b] - rigid[b]) > 0.1);
    }
}

// The domain above the roofs grows with the march's length at high frequency
// (profile_pe.cpp): a layer a fixed number of wavelengths thick reflects the
// waves meeting it near grazing, and with 15 a street's parts were 1.1 dB
// off at 3.15 kHz and 3 dB at 8 kHz.  Against the same march in a domain 250
// wavelengths high (and its Huygens line as high), from a source 100 m off
// and from a microphone 4 m before the row.
TEST_CASE("the profile field is converged in its domain at high frequency", "[profile_pe]")
{
    struct Reference
    {
        ProfileProblem problem;
        double f;
        double level[2][2];
    };
    ProfileProblem far, near;
    far.blocks = {{100, 110, 10}, {120, 130, 10}};
    far.sz = 12;
    far.ru = 134;
    far.rz = 1.2;
    near.blocks = {{4, 14, 10}, {24, 34, 10}};
    near.sz = 1.2;
    near.ru = 134;
    near.rz = 12;
    const std::vector<Reference> references {
            {far, 3150, {{-53.497, -54.661}, {-78.040, -79.220}}},
            {far, 8000, {{-55.598, -56.761}, {-86.146, -87.318}}},
            {near, 3150, {{-54.713, -78.951}, {-55.912, -80.150}}},
            {near, 8000, {{-56.823, -87.057}, {-58.022, -88.256}}},
    };
    for (const Reference& r : references)
    {
        const ProfileParts parts = profile_parts(r.problem, 2.0 * pi * r.f / 343.2);
        for (int a = 0; a < 2; a++)
            for (int b = 0; b < 2; b++)
            {
                INFO(r.f << " Hz, source " << r.problem.sz << " m up, part " << a << b);
                REQUIRE(20.0 * std::log10(std::abs(parts.part[a][b])) == Catch::Approx(r.level[a][b]).margin(0.12));
            }
    }
}

// The spacing puts the roofs on nodes (aligned_spacing): a roof between two
// was taken at the nearer, and the level stepped by up to 0.2 dB as the
// frequency moved it across.  Two roof heights, 480 to 520 Hz in half-hertz
// steps: the level's largest second difference was 0.33 dB.
TEST_CASE("the profile field is smooth in frequency", "[profile_pe]")
{
    ProfileProblem p;
    p.blocks = {{100, 110, 10}, {120, 130, 8}};
    p.sz = 12;
    p.ru = 134;
    p.rz = 1.2;
    std::vector<double> level;
    for (int i = 0; i <= 80; i++)
        level.push_back(20.0 * std::log10(std::abs(profile_parts(p, 2.0 * pi * (480.0 + 0.5 * i) / 343.2).part[0][0])));
    double worst = 0.0;
    for (std::size_t i = 1; i + 1 < level.size(); i++)
        worst = std::max(worst, std::abs(level[i + 1] - 2.0 * level[i] + level[i - 1]));
    REQUIRE(worst < 0.03);
}
