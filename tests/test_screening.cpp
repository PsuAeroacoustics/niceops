#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <random>
#include <vector>

#include "faddeeva.hpp"
#include "screening.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace)

namespace {

// The upper hull by brute force: a point is a vertex when it lies strictly
// above the segment between every pair of points on either side of it.
std::vector<std::size_t> brute_force_hull(const std::vector<ScreenPoint>& p)
{
    std::vector<std::size_t> hull;
    for (std::size_t i = 0; i < p.size(); i++)
    {
        bool vertex = true;
        for (std::size_t j = 0; j < i && vertex; j++)
            for (std::size_t k = i + 1; k < p.size() && vertex; k++)
                vertex = (p[k].s - p[j].s) * (p[i].z - p[j].z) - (p[k].z - p[j].z) * (p[i].s - p[j].s) > 0.0;
        if (vertex)
            hull.push_back(i);
    }
    return hull;
}

std::vector<std::size_t> hull_of(const std::vector<ScreenPoint>& p)
{
    std::vector<std::size_t> hull(p.size());
    hull.resize(upper_hull(p.data(), p.size(), hull.data()));
    return hull;
}

// Sommerfeld's F(a) = exp(-i pi/4)/sqrt(pi) * integral from -inf to a of
// exp(i t^2) dt = erfc(-exp(-i pi/4) a) / 2, through the Faddeeva function:
// with u = exp(-i pi/4) |a|, erfc(u) = exp(i a^2) w(exp(i pi/4) |a|), whose
// argument is in the first quadrant.
std::complex<double> sommerfeld_f(double a)
{
    const std::complex<double> tail = 0.5 * std::exp(std::complex<double>(0.0, a * a))
            * faddeeva(std::polar(std::fabs(a), pi / 4.0));
    return a < 0.0 ? tail : 1.0 - tail;
}

// The exact field behind a rigid half-plane lit by a plane wave (Sommerfeld
// 1896), relative to the incident wave: the half-plane along theta = 0, the
// wave arriving from direction theta0, observed at (r, theta), k the
// wavenumber.  Returned as an insertion loss, dB.
double half_plane_loss(double kr, double theta, double theta0)
{
    const auto term = [&](double angle) {
        return std::exp(std::complex<double>(0.0, -kr * std::cos(angle)))
                * sommerfeld_f(std::sqrt(2.0 * kr) * std::cos(0.5 * angle));
    };
    return -20.0 * std::log10(std::abs(term(theta - theta0) + term(theta + theta0)));
}

} // namespace

TEST_CASE("the upper hull is the brute-force hull", "[screening]")
{
    std::mt19937 random(3);
    std::uniform_real_distribution<double> height(-5.0, 25.0);
    std::uniform_int_distribution<int> level(0, 4);
    for (int trial = 0; trial < 400; trial++)
    {
        const std::size_t n = 2 + trial % 40;
        std::vector<ScreenPoint> p(n);
        for (std::size_t k = 0; k < n; k++)
        {
            p[k].s = 10.0 * static_cast<double>(k);
            // Whole heights on a few levels give plateaus and collinear runs;
            // the rest are random.
            p[k].z = trial % 3 == 0 ? static_cast<double>(level(random))
                    : trial % 3 == 1 ? 2.0 * static_cast<double>(k % 5) : height(random);
        }
        if (trial % 7 == 0)
            p.back().z = 60.0; // a raised source
        REQUIRE(hull_of(p) == brute_force_hull(p));
    }
    // Collinear: only the ends.
    const std::vector<ScreenPoint> line {{0.0, 1.0}, {10.0, 3.5}, {20.0, 6.0}, {30.0, 8.5}};
    CHECK(hull_of(line) == std::vector<std::size_t> {0, 3});
    // Two points (a path shorter than a step): the ends.
    CHECK(hull_of({{0.0, 1.2}, {4.0, 100.0}}) == std::vector<std::size_t> {0, 1});
}

TEST_CASE("the rubber band finds the edges, the path difference and e", "[screening]")
{
    std::vector<ScreenPoint> points, edges;
    std::vector<std::size_t> hull;
    Screen screen;
    // A single ridge: 10 m high at s = 50, the microphone 1 m up, the source
    // 2 m up 100 m out.
    std::vector<double> heights(11, 0.0);
    heights[5] = 10.0;
    REQUIRE(rubber_band(heights, 9, 10.0, 100.0, 1.0, 2.0, points, hull, edges, screen));
    CHECK(screen.edges == 1);
    CHECK(screen.receiver_edge == 5);
    CHECK(screen.source_edge == 5);
    CHECK(screen.e == 0.0);
    CHECK(screen.delta
            == Catch::Approx(std::hypot(50.0, 9.0) + std::hypot(50.0, 8.0) - std::hypot(100.0, 1.0)).epsilon(1e-14));
    // Two ridges: e is the distance between them.
    heights[8] = 10.0;
    REQUIRE(rubber_band(heights, 9, 10.0, 100.0, 1.0, 2.0, points, hull, edges, screen));
    CHECK(screen.edges == 2);
    CHECK(screen.receiver_edge == 5);
    CHECK(screen.source_edge == 8);
    CHECK(screen.e == 30.0);
    CHECK(screen.delta
            == Catch::Approx(std::hypot(50.0, 9.0) + 30.0 + std::hypot(20.0, 8.0) - std::hypot(100.0, 1.0))
                       .epsilon(1e-14));
    // Nothing above the line of sight: no band.
    std::fill(heights.begin(), heights.end(), 0.0);
    CHECK(!rubber_band(heights, 9, 10.0, 100.0, 1.0, 2.0, points, hull, edges, screen));
}

TEST_CASE("Delta_dif follows the Directive's formula", "[screening]")
{
    // One edge: 10 lg(3 + 40 delta / lambda).
    CHECK(cnossos_dif(0.1, 0.0, 1, 1.0) == Catch::Approx(10.0 * std::log10(7.0)).epsilon(1e-15));
    // Toward grazing, 10 lg 3.
    CHECK(cnossos_dif(1e-15, 0.0, 1, 1.0) == Catch::Approx(10.0 * std::log10(3.0)).epsilon(1e-12));
    // C'' between 1 (edges close together) and 3 (far apart).
    const double lambda = 0.5;
    const double r = 5.0 * lambda / 4.0;
    CHECK(cnossos_dif(0.1, 4.0, 2, lambda)
            == Catch::Approx(10.0 * std::log10(3.0 + 40.0 / lambda * (1.0 + r * r) / (1.0 / 3.0 + r * r) * 0.1))
                       .epsilon(1e-15));
    CHECK(cnossos_dif(0.1, 1e6, 2, lambda) == Catch::Approx(10.0 * std::log10(3.0 + 3.0 * 40.0 / lambda * 0.1)).epsilon(1e-9));
    // Edges within 0.3 m are one, and so is a single edge whatever e says.
    CHECK(cnossos_dif(0.1, 0.3, 2, lambda) == cnossos_dif(0.1, 0.0, 1, lambda));
    CHECK(cnossos_dif(0.1, 50.0, 1, lambda) == cnossos_dif(0.1, 0.0, 1, lambda));
    // Below 40 C'' delta / lambda = -2 there is no diffraction term; at -2 it is 0.
    CHECK(cnossos_dif(-0.06, 0.0, 1, 1.0) == 0.0);
    CHECK(cnossos_dif(-0.05, 0.0, 1, 1.0) == Catch::Approx(0.0).margin(1e-14));
}

TEST_CASE("the screened loss: the cap, the side terms and the F clamp", "[screening]")
{
    const std::vector<double> bands {63.0, 250.0, 1000.0, 4000.0};
    const double c = 340.0;
    std::vector<double> loss(bands.size());
    Screen screen;
    screen.edges = 1;
    screen.delta = 0.4;

    SECTION("without ground, Delta_dif alone, capped at 25 dB")
    {
        ScreenTally tally;
        cnossos_screen_loss(screen, bands, c, nullptr, nullptr, 0, bands.size(), loss.data(), &tally);
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            const double expected = std::min(25.0, 10.0 * std::log10(3.0 + 40.0 * bands[b] / c * 0.4));
            CHECK(loss[b] == Catch::Approx(expected).epsilon(1e-12));
        }
        // 4 kHz: 10 lg(3 + 188) = 22.8, under the cap; 40 m of path difference is over it.
        CHECK(!tally.capped);
        screen.delta = 40.0;
        cnossos_screen_loss(screen, bands, c, nullptr, nullptr, 0, bands.size(), loss.data(), &tally);
        CHECK(tally.capped);
        CHECK(loss[3] == 25.0);
    }

    SECTION("each side's ground, weighted by its image's extra diffraction")
    {
        screen.delta_source_image = 0.9;
        screen.delta_receiver_image = 0.5;
        const std::vector<double> gs {4.0, 0.5, 2.0, 0.1}, gr {3.0, 1.0, 0.2, 4.0};
        ScreenTally tally;
        cnossos_screen_loss(screen, bands, c, gs.data(), gr.data(), 0, bands.size(), loss.data(), &tally);
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            const double lambda = c / bands[b];
            const double dif = cnossos_dif(0.4, 0.0, 1, lambda);
            const double fs = std::pow(10.0, -(cnossos_dif(0.9, 0.0, 1, lambda) - dif) / 20.0);
            const double fr = std::pow(10.0, -(cnossos_dif(0.5, 0.0, 1, lambda) - dif) / 20.0);
            const double expected = std::min(dif, 25.0) - 20.0 * std::log10(1.0 + (std::sqrt(gs[b]) - 1.0) * fs)
                    - 20.0 * std::log10(1.0 + (std::sqrt(gr[b]) - 1.0) * fr);
            CHECK(loss[b] == Catch::Approx(expected).epsilon(1e-12));
        }
        CHECK(!tally.source_clamped);
        CHECK(!tally.receiver_clamped);
        // The images exactly as deep as the real path: F = 1, and the sides
        // are the plain ground gains.
        screen.delta_source_image = screen.delta_receiver_image = 0.4;
        cnossos_screen_loss(screen, bands, c, gs.data(), gr.data(), 0, bands.size(), loss.data(), &tally);
        for (std::size_t b = 0; b < bands.size(); b++)
            CHECK(loss[b]
                    == Catch::Approx(cnossos_dif(0.4, 0.0, 1, c / bands[b]) - 10.0 * std::log10(gs[b])
                                     - 10.0 * std::log10(gr[b]))
                               .epsilon(1e-12));
    }

    SECTION("an image shallower than the real path holds F at 1")
    {
        screen.delta_source_image = 0.1;
        screen.delta_receiver_image = -1.0; // under the -2 cutoff: no Delta_dif at all
        const std::vector<double> g(bands.size(), 0.25);
        ScreenTally tally;
        cnossos_screen_loss(screen, bands, c, g.data(), g.data(), 0, bands.size(), loss.data(), &tally);
        CHECK(tally.source_clamped);
        CHECK(tally.receiver_clamped);
        for (std::size_t b = 0; b < bands.size(); b++)
            CHECK(loss[b] == Catch::Approx(cnossos_dif(0.4, 0.0, 1, c / bands[b]) - 20.0 * std::log10(0.25)).epsilon(1e-12));
    }
}

TEST_CASE("screen_sides mirrors in each side's line, and a plate keeps R", "[screening]")
{
    std::vector<ScreenPoint> points, edges;
    std::vector<std::size_t> hull;
    Screen screen;
    std::vector<double> heights(11, 0.0);
    heights[5] = 20.0;
    REQUIRE(rubber_band(heights, 9, 10.0, 100.0, 1.2, 30.0, points, hull, edges, screen));
    // Level lines through the ground: images straight down, and the side
    // geometry is the plain one.
    screen_sides(screen, edges, GroundLine {0.0, 0.0}, GroundLine {0.0, 0.0}, false);
    CHECK(screen.delta_source_image
            == Catch::Approx(std::hypot(50.0, 18.8) + std::hypot(50.0, 50.0) - std::hypot(100.0, 31.2)).epsilon(1e-13));
    CHECK(screen.delta_receiver_image
            == Catch::Approx(std::hypot(50.0, 21.2) + std::hypot(50.0, 10.0) - std::hypot(100.0, 31.2)).epsilon(1e-13));
    CHECK(screen.source_side.source_height == 30.0);
    CHECK(screen.source_side.receiver_height == 20.0);
    CHECK(screen.source_side.ground_distance == 50.0);
    CHECK(screen.receiver_side.source_height == 20.0);
    CHECK(screen.receiver_side.receiver_height == 1.2);
    CHECK(screen.receiver_side.ground_distance == 50.0);
    // A flush microphone: R' = R, and its side's height is the edge's above it.
    REQUIRE(rubber_band(heights, 9, 10.0, 100.0, 0.0, 30.0, points, hull, edges, screen));
    screen_sides(screen, edges, GroundLine {0.0, 0.0}, GroundLine {-0.5, 0.0}, true);
    CHECK(screen.delta_receiver_image == screen.delta);
    CHECK(screen.receiver_side.source_height == 20.0);
    CHECK(screen.receiver_side.receiver_height == 0.0);
    // A tilted source-side line: S' is S mirrored across it, so the distance
    // from S' to the line equals S's.
    const GroundLine tilted {-20.0, 0.3};
    REQUIRE(rubber_band(heights, 9, 10.0, 100.0, 1.2, 30.0, points, hull, edges, screen));
    screen_sides(screen, edges, tilted, GroundLine {0.0, 0.0}, false);
    const double norm = std::sqrt(1.0 + 0.09);
    const double s_above = (30.0 - (-20.0 + 0.3 * 100.0)) / norm;
    CHECK(screen.source_side.source_height == Catch::Approx(s_above).epsilon(1e-14));
}

TEST_CASE("a screened pair never exceeds the bound's screened term", "[screening][bounds]")
{
    // Random profiles, sources, side lines and side ground gains from silent to
    // a surface wave's: the screened power ratio 10^(-A/10) must be within
    // screen_gain_bound(G_S, G_R), because F is held to [0, 1] and
    // Delta_dif(S,R) >= 10 lg 3.  Rigid ground (G = 4) at grazing comes within
    // a few percent of it, so a bound without the 1/3, or without G_S, fails.
    std::mt19937 random(17);
    std::uniform_real_distribution<double> height(0.0, 30.0), gain(0.0, 6.0), slope(-0.6, 0.6), unit(0.0, 1.0);
    const std::vector<double> bands {20.0, 63.0, 250.0, 1000.0, 4000.0, 10000.0};
    std::vector<double> loss(bands.size()), gs(bands.size()), gr(bands.size());
    std::vector<ScreenPoint> points, edges;
    std::vector<std::size_t> hull;
    double closest = 0.0;
    int screened = 0, clamped = 0;
    for (int trial = 0; trial < 20000; trial++)
    {
        const std::size_t interior = 1 + trial % 60;
        std::vector<double> heights(interior + 2);
        for (auto& h : heights)
            h = height(random) * (trial % 5 == 0 ? 0.05 : 1.0);
        const double step = 10.0;
        const double distance = step * (static_cast<double>(interior) + 0.5 + 0.5 * unit(random));
        Screen screen;
        if (!rubber_band(heights, interior, step, distance, 1.2 * unit(random), 40.0 * unit(random), points, hull,
                    edges, screen))
            continue;
        screened++;
        screen_sides(screen, edges, GroundLine {10.0 * (unit(random) - 0.5), slope(random)},
                GroundLine {10.0 * (unit(random) - 0.5), slope(random)}, trial % 4 == 0);
        const bool rigid = trial % 3 == 0;
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            gs[b] = rigid ? 4.0 : gain(random);
            gr[b] = rigid ? 4.0 : gain(random);
        }
        ScreenTally tally;
        cnossos_screen_loss(screen, bands, 340.0, gs.data(), gr.data(), 0, bands.size(), loss.data(), &tally);
        clamped += tally.source_clamped || tally.receiver_clamped;
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            const double ratio = std::pow(10.0, -loss[b] / 10.0);
            const double bound = screen_gain_bound(gs[b], gr[b]);
            REQUIRE(ratio <= bound * (1.0 + 1e-12));
            closest = std::max(closest, ratio / bound);
        }
    }
    CHECK(screened > 10000);
    CHECK(clamped > 100);
    CHECK(closest > 0.97);
}

TEST_CASE("Delta_dif against the exact half-plane", "[screening]")
{
    // Sommerfeld's exact solution for a rigid half-plane lit by a plane wave,
    // at kr = 2000 from the edge, against the CNOSSOS term at the same path
    // difference, delta = r (1 + cos(theta - theta0)).  On the shadow
    // boundary the exact field is half the incident one, 6.02 dB down, where
    // CNOSSOS gives 10 lg 3 = 4.77 dB: 1.25 dB low at grazing.  Just inside the
    // shadow CNOSSOS reads up to 2.7 dB low (a degree in, where its slope in
    // delta is the shallower), and by 60 degrees, where the exact solution's
    // second term counts, it reads 0.7 dB high.
    const double kr = 2000.0, theta0 = pi / 4.0;
    const double lambda = 1.0, r = kr * lambda / (2.0 * pi);
    const double boundary = pi + theta0;
    const double grazing_exact = half_plane_loss(kr, boundary, theta0);
    CHECK(grazing_exact == Catch::Approx(20.0 * std::log10(2.0)).margin(0.1));
    CHECK(grazing_exact - cnossos_dif(1e-12, 0.0, 1, lambda) == Catch::Approx(1.25).margin(0.1));
    for (const double degrees : {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0})
    {
        const double theta = boundary + degrees * pi / 180.0;
        const double delta = r * (1.0 + std::cos(theta - theta0));
        const double exact = half_plane_loss(kr, theta, theta0);
        const double cnossos = cnossos_dif(delta, 0.0, 1, lambda);
        INFO(degrees << " degrees into the shadow: exact " << exact << " dB, CNOSSOS " << cnossos << " dB");
        CHECK(exact - cnossos > -1.0);
        CHECK(exact - cnossos < 2.8);
    }
}
