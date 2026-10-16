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

// The CNOSSOS-EU conformance cases of ISO/TR 17534-4:2020 that a vertical-plane
// terrain model can take: no buildings, barriers or reflections.  The expected
// values are the TR's, as NoiseModelling's AttenuationComputeOutputCnossosTest
// transcribes them, rounded to 0.01 dB or 0.01 m; homogeneous conditions only,
// since NICEOPS's screening is straight-ray.  Octave bands 63 Hz-8 kHz at the
// TR's c = 340 m/s.  Each profile is the TR's polyline, sampled every 0.01 m
// (which puts every vertex on a sample), seen from the receiver as NICEOPS
// sees a path: s from the microphone toward the source.
namespace {

const std::vector<double> tr_bands {63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0};
constexpr double tr_sound_speed = 340.0;
constexpr double tr_step = 0.01;

struct TrCase
{
    std::vector<double> heights; // heights[k] at s = k tr_step, on the microphone's ground
    TerrainRaster raster;        // the same profile along x, for the line fits
    double distance = 0.0;       // s of the source
};

// `vertices` are the TR's (x from the source, z); the receiver stands at the
// last x, its ground the profile's there.
TrCase tr_case(const std::vector<ScreenPoint>& vertices)
{
    TrCase c;
    const double length = vertices.back().s;
    const double base = vertices.back().z;
    const std::size_t n = static_cast<std::size_t>(std::llround(length / tr_step));
    c.distance = length;
    const auto z_at = [&](double x) {
        for (std::size_t v = 1; v < vertices.size(); v++)
        {
            if (x <= vertices[v].s + 1e-9)
            {
                const double t = (x - vertices[v - 1].s) / (vertices[v].s - vertices[v - 1].s);
                return vertices[v - 1].z + std::clamp(t, 0.0, 1.0) * (vertices[v].z - vertices[v - 1].z) - base;
            }
        }
        return vertices.back().z - base;
    };
    c.heights.resize(n + 1);
    for (std::size_t k = 0; k <= n; k++)
        c.heights[k] = z_at(length - static_cast<double>(k) * tr_step);
    // Three rows, so the bilinear reads along y = 0 are the profile itself.
    std::vector<double> grid(3 * (n + 1));
    for (std::size_t j = 0; j < 3; j++)
        std::copy(c.heights.begin(), c.heights.end(), grid.begin() + static_cast<std::ptrdiff_t>(j * (n + 1)));
    c.raster = TerrainRaster(0.0, -tr_step, tr_step, n + 1, 3, std::move(grid));
    return c;
}

// The least-squares line over samples k0..k1 (fit_ground_line, as fit_range
// approximates it), with its intercept moved to the microphone, s = 0.
GroundLine tr_line(const TrCase& c, std::size_t k0, std::size_t k1)
{
    const double s0 = static_cast<double>(k0) * tr_step, s1 = static_cast<double>(k1) * tr_step;
    GroundLine line = fit_ground_line(c.raster, s0, 0.0, s1, 0.0, k1 - k0);
    line.intercept -= line.slope * s0;
    return line;
}

// The TR's per-band values against ours, to `tolerance` dB.
void check_bands(const char* what, const std::vector<double>& expected, const std::vector<double>& actual,
        double tolerance)
{
    for (std::size_t b = 0; b < expected.size(); b++)
    {
        INFO(what << " at " << tr_bands[b] << " Hz");
        CHECK(actual[b] == Catch::Approx(expected[b]).margin(tolerance));
    }
}

} // namespace

TEST_CASE("ISO/TR 17534-4 TC23: diffraction over an earth berm", "[screening][iso17534]")
{
    // TR tables 264 and 268-; the source 1 m and the receiver 4 m above level
    // ground, a berm 5 m high between, ground G 0, 1, 0 along the path.
    const TrCase c = tr_case({{0.0, 0.0}, {14.21, 0.0}, {22.64, 5.0}, {23.98, 5.0}, {32.30, 0.0}, {70.03, 0.0}});
    std::vector<ScreenPoint> points, edges;
    std::vector<std::size_t> hull;
    Screen screen;
    const std::size_t interior = c.heights.size() - 2;
    REQUIRE(rubber_band(c.heights, interior, tr_step, c.distance, 4.0, 1.0, points, hull, edges, screen));
    // The berm's two crest corners, 1.34 m apart.
    CHECK(screen.edges == 2);
    CHECK(screen.e == Catch::Approx(1.34).margin(1e-6));

    std::vector<double> dif(tr_bands.size());
    for (std::size_t b = 0; b < tr_bands.size(); b++)
        dif[b] = cnossos_dif(screen.delta, screen.e, screen.edges, tr_sound_speed / tr_bands[b]);
    // Uncapped: 29.05 dB at 8 kHz, which the cap below takes to 25.
    check_bands("Delta_dif(S,R)", {7.17, 8.69, 10.78, 13.46, 16.98, 21.34, 25.53, 29.05}, dif, 0.02);

    // Each side's mean line, as the TR's segment planes (table 268): the
    // source's from its edge to the footprint, the receiver's from the
    // microphone to its edge.
    screen_sides(screen, edges, tr_line(c, screen.source_edge, interior + 1), tr_line(c, 0, screen.receiver_edge),
            false);
    // z_s, z_r and d_p of the source side (S over it, its edge O over it) ...
    CHECK(screen.source_side.source_height == Catch::Approx(2.13).margin(0.01));
    CHECK(screen.source_side.receiver_height == Catch::Approx(1.94).margin(0.01));
    CHECK(screen.source_side.ground_distance == Catch::Approx(22.99).margin(0.01));
    // ... and of the receiver side (its edge O', then R).
    CHECK(screen.receiver_side.source_height == Catch::Approx(3.35).margin(0.01));
    CHECK(screen.receiver_side.receiver_height == Catch::Approx(4.73).margin(0.01));
    CHECK(screen.receiver_side.ground_distance == Catch::Approx(46.04).margin(0.01));
    std::vector<double> source_image(tr_bands.size()), receiver_image(tr_bands.size());
    for (std::size_t b = 0; b < tr_bands.size(); b++)
    {
        const double lambda = tr_sound_speed / tr_bands[b];
        source_image[b] = cnossos_dif(screen.delta_source_image, screen.e, screen.edges, lambda);
        receiver_image[b] = cnossos_dif(screen.delta_receiver_image, screen.e, screen.edges, lambda);
    }
    check_bands("Delta_dif(S',R)", {10.54, 12.93, 15.67, 18.78, 22.56, 27.04, 31.28, 34.81}, source_image, 0.03);
    check_bands("Delta_dif(S,R')", {10.86, 13.28, 16.05, 19.18, 22.97, 27.46, 31.70, 35.22}, receiver_image, 0.03);

    // With the TR's own A_ground on each side (its G-based ground, not
    // NICEOPS's model, injected as G = 10^(-A/10)), the combination gives the
    // TR's Delta_ground terms and A_dif, the cap included.
    const std::vector<double> a_ground_so {-1.02, -0.08, -0.73, -2.79, -2.79, -2.79, -2.79, -2.79};
    const std::vector<double> a_ground_or {-0.28, -1.00, -2.46, -2.46, -2.46, -2.46, -2.46, -2.46};
    std::vector<double> gain_so(tr_bands.size()), gain_or(tr_bands.size()), unity(tr_bands.size(), 1.0);
    for (std::size_t b = 0; b < tr_bands.size(); b++)
    {
        gain_so[b] = std::pow(10.0, -a_ground_so[b] / 10.0);
        gain_or[b] = std::pow(10.0, -a_ground_or[b] / 10.0);
    }
    std::vector<double> loss(tr_bands.size()), source_only(tr_bands.size()), receiver_only(tr_bands.size());
    ScreenTally tally;
    cnossos_screen_loss(screen, tr_bands, tr_sound_speed, gain_so.data(), gain_or.data(), 0, tr_bands.size(),
            loss.data(), &tally);
    CHECK(tally.capped);
    CHECK(!tally.source_clamped);
    CHECK(!tally.receiver_clamped);
    cnossos_screen_loss(screen, tr_bands, tr_sound_speed, gain_so.data(), unity.data(), 0, tr_bands.size(),
            source_only.data());
    cnossos_screen_loss(screen, tr_bands, tr_sound_speed, unity.data(), gain_or.data(), 0, tr_bands.size(),
            receiver_only.data());
    std::vector<double> delta_ground_so(tr_bands.size()), delta_ground_or(tr_bands.size());
    for (std::size_t b = 0; b < tr_bands.size(); b++)
    {
        const double capped = std::min(dif[b], cnossos_cap);
        delta_ground_so[b] = source_only[b] - capped;
        delta_ground_or[b] = receiver_only[b] - capped;
    }
    check_bands("Delta_ground(S,O)", {-0.71, -0.05, -0.42, -1.62, -1.58, -1.56, -1.55, -1.55}, delta_ground_so, 0.02);
    check_bands("Delta_ground(O,R)", {-0.18, -0.60, -1.42, -1.36, -1.32, -1.30, -1.30, -1.29}, delta_ground_or, 0.02);
    check_bands("A_dif", {6.28, 8.04, 8.93, 10.48, 14.09, 18.48, 22.15, 22.16}, loss, 0.03);
}

TEST_CASE("ISO/TR 17534-4 TC06: a path just over a crest is not screened", "[screening][iso17534]")
{
    // TR table 23: the source 1 m up, the receiver 1.5 m above a plateau 10 m
    // up, whose edge the line of sight clears by 0.67 m.  The TR screens the
    // 500 Hz and 1 kHz bands all the same (Delta_dif(S,R) 3.16 and 0.56 dB,
    // A_dif 4.31 and -0.83 dB), by the Directive's -lambda/20 rule; NICEOPS
    // screens only a profile above the line of sight (ground.md, Screening),
    // so this path is unscreened, and the difference is the documented one.
    const TrCase c = tr_case({{0.0, 0.0}, {112.41, 0.0}, {178.84, 10.0}, {194.16, 10.0}});
    std::vector<ScreenPoint> points, edges;
    std::vector<std::size_t> hull;
    Screen screen;
    CHECK(!rubber_band(c.heights, c.heights.size() - 2, tr_step, c.distance, 1.5, 1.0 - 10.0, points, hull, edges,
            screen));
    // What the TR does there: a path clearing its edge (delta < 0) is
    // diffracted in a band only if delta > -lambda/20 and delta > lambda/4 -
    // delta', delta' the path difference by the ground images.  The second
    // clause leaves 63-250 Hz to the ground alone; the first, 2-8 kHz.  In the
    // two bands it diffracts, the term is the Directive's formula at the
    // negative path difference.
    const ScreenPoint receiver {0.0, 1.5}, source {c.distance, -9.0};
    const ScreenPoint crest {194.16 - 178.84, 0.0};
    const double delta = path_difference(receiver, &crest, 1, source);
    REQUIRE(delta < 0.0);
    const std::vector<double> tr_dif {0.0, 0.0, 0.0, 3.16, 0.56, 0.0, 0.0, 0.0};
    for (std::size_t b = 0; b < tr_bands.size(); b++)
    {
        const double lambda = tr_sound_speed / tr_bands[b];
        INFO(tr_bands[b] << " Hz");
        if (tr_dif[b] > 0.0)
            CHECK(cnossos_dif(delta, 0.0, 1, lambda) == Catch::Approx(tr_dif[b]).margin(0.02));
        else if (tr_bands[b] > 1000.0)
            CHECK(delta <= -lambda / 20.0);
    }
}

TEST_CASE("ISO/TR 17534-4 TC05 and TC20: the mean plane over uneven ground", "[screening][iso17534]")
{
    // Unobstructed paths over a rise: the whole path's least-squares line
    // gives the heights and distance a table over terrain reads (z_s, z_r,
    // d_p).  The source 1 m and the receiver 4 m up.  Checked against the
    // polyline's exact (continuous) least-squares plane to 2 mm, and against
    // the TR's rounded values: TC05's (table 18) to 0.01 m; TC20's, which
    // NoiseModelling takes from the table it labels TC21's and checks only to
    // 0.1 m, sit 0.014 and 0.038 m from the exact plane, and are held to 0.05.
    struct Expect
    {
        std::vector<ScreenPoint> profile;
        double zs, zr, dp;          // the exact plane's
        double tr_zs, tr_zr, tr_dp; // the TR's
        double tr_margin;
    };
    for (const Expect& e : {Expect {{{0.0, 0.0}, {112.41, 0.0}, {178.84, 10.0}, {194.16, 10.0}}, 3.826, 6.159, 194.581,
                                    3.83, 6.16, 194.59, 0.01},
                 Expect {{{0.0, 0.0}, {110.34, 0.0}, {175.54, 10.0}, {190.59, 10.0}}, 3.827, 6.158, 191.019, 3.84, 6.12,
                         191.02, 0.05}})
    {
        const TrCase c = tr_case(e.profile);
        const GroundLine line = tr_line(c, 0, c.heights.size() - 1);
        const PlaneGeometry g = plane_geometry(c.distance, 1.0 - 10.0, 4.0, line);
        INFO("profile ending at x = " << c.distance);
        CHECK(g.source_height == Catch::Approx(e.zs).margin(0.002));
        CHECK(g.receiver_height == Catch::Approx(e.zr).margin(0.002));
        CHECK(g.ground_distance == Catch::Approx(e.dp).margin(0.002));
        CHECK(g.source_height == Catch::Approx(e.tr_zs).margin(e.tr_margin));
        CHECK(g.receiver_height == Catch::Approx(e.tr_zr).margin(e.tr_margin));
        CHECK(g.ground_distance == Catch::Approx(e.tr_dp).margin(e.tr_margin));
    }
}

// --terrain_screening kirchhoff: the knife edge.
TEST_CASE("the Fresnel integrals against reference values", "[screening][kirchhoff]")
{
    // C(x) and S(x) from SciPy's scipy.special.fresnel (Cephes), a separate
    // implementation: both branches, either side of the switch at 1.5, and far
    // out where the continued fraction carries them.
    struct Row
    {
        double x, c, s;
    };
    const Row rows[] = {
        {-3.0, -0.60572078929768569, -0.49631299896737502},
        {-0.5, -0.49234422587144638, -0.064732432859999287},
        {0.0, 0.0, 0.0},
        {0.001, 0.00099999999999975326, 5.2359877559820663e-10},
        {0.1, 0.099997532627085065, 0.00052358954761221076},
        {0.5, 0.49234422587144638, 0.064732432859999287},
        {1.0, 0.77989340037682298, 0.43825914739035471},
        {1.4999, 0.44535357300621869, 0.69754320665599856},
        {1.5, 0.44526117603982146, 0.69750496008209295},
        {1.5001, 0.445168797106956, 0.69746666997139084},
        {2.0, 0.48825340607534073, 0.34341567836369824},
        {3.0, 0.60572078929768569, 0.49631299896737502},
        {5.0, 0.56363118870401219, 0.49919138191711687},
        {10.0, 0.49989869420551575, 0.46816997858488224},
        {50.0, 0.49999918943072796, 0.49363380258593875},
        {200.0, 0.49999998733485207, 0.49840845056938343},
    };
    for (const Row& row : rows)
    {
        double c = 0.0, s = 0.0;
        fresnel_integrals(row.x, c, s);
        INFO("x = " << row.x);
        CHECK(c == Catch::Approx(row.c).margin(1e-13));
        CHECK(s == Catch::Approx(row.s).margin(1e-13));
    }
}

TEST_CASE("the knife edge's power: grazing, the shadow and its asymptote", "[screening][kirchhoff]")
{
    // |F(v)|^2 from the same SciPy values.
    const double reference[][2] = {{0.0, 0.25}, {0.5, 0.094758232940741233}, {1.0, 0.041076124227733234},
            {2.0, 0.012328316125565821}, {5.0, 0.002024791019544787}, {20.0, 0.00012665107854758868},
            {100.0, 5.0660591564519346e-06}};
    for (const auto& row : reference)
    {
        INFO("v = " << row[0]);
        CHECK(knife_edge_power(row[0]) == Catch::Approx(row[1]).epsilon(1e-11));
    }
    // It falls throughout the shadow, toward 1/(2 pi^2 v^2), and in the lit
    // region it nears the free field, oscillating about it.
    double previous = knife_edge_power(0.0);
    for (double v = 0.05; v < 50.0; v += 0.05)
    {
        const double p = knife_edge_power(v);
        CHECK(p < previous);
        previous = p;
    }
    CHECK(knife_edge_power(200.0) * 2.0 * pi * pi * 200.0 * 200.0 == Catch::Approx(1.0).epsilon(1e-4));
    CHECK(knife_edge_power(-20.0) == Catch::Approx(1.0).margin(0.02));
}

TEST_CASE("the knife edge against the exact half-plane", "[screening][kirchhoff]")
{
    // Sommerfeld's rigid half-plane at kr = 2000, as the CNOSSOS case above,
    // against the Fresnel knife edge at the same path difference, v = 2
    // sqrt(delta / lambda).  On the boundary both are 6.02 dB down.  The
    // Fresnel-Kirchhoff edge is a black screen, where the rigid half-plane
    // also sends its faces' reflection into the shadow, so it is always the
    // deeper: by 0.11 dB half a degree in, 0.50 dB at 5 degrees, 0.94 dB at 10
    // and 3.6 dB at 60.
    const double kr = 2000.0, theta0 = pi / 4.0;
    const double lambda = 1.0, r = kr * lambda / (2.0 * pi);
    const double boundary = pi + theta0;
    CHECK(-10.0 * std::log10(knife_edge_power(0.0)) == Catch::Approx(20.0 * std::log10(2.0)).margin(1e-9));
    for (const double degrees : {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0})
    {
        const double theta = boundary + degrees * pi / 180.0;
        const double delta = r * (1.0 + std::cos(theta - theta0));
        const double exact = half_plane_loss(kr, theta, theta0);
        const double knife = -10.0 * std::log10(knife_edge_power(2.0 * std::sqrt(delta / lambda)));
        INFO(degrees << " degrees into the shadow: exact " << exact << " dB, knife edge " << knife << " dB");
        CHECK(knife - exact > 0.0);
        CHECK(knife - exact < (degrees <= 5.0 ? 0.6 : 4.0));
    }
}

TEST_CASE("the knife edge's band loss is at least 10 lg 4 and rises with delta", "[screening][kirchhoff]")
{
    const std::vector<double> bands {50.0, 63.0, 80.0, 100.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0};
    std::vector<double> loss(bands.size()), previous(bands.size(), 0.0);
    for (const double delta : {1e-9, 1e-6, 1e-3, 0.01, 0.1, 1.0, 10.0, 100.0})
    {
        KnifeEdge edge;
        edge.delta = delta;
        knife_edge_loss(edge, bands, 343.0, 0, bands.size(), loss.data());
        for (std::size_t band = 0; band < bands.size(); band++)
        {
            INFO("delta " << delta << " m, " << bands[band] << " Hz: " << loss[band] << " dB");
            CHECK(loss[band] > 20.0 * std::log10(2.0));
            CHECK(loss[band] > previous[band]);
            if (band > 0)
                CHECK(loss[band] > loss[band - 1]);
            // The band's mean lies between its edges' values.
            const double lo = knife_edge_power(2.0 * std::sqrt(delta * bands[band] * std::pow(2.0, -1.0 / 6.0) / 343.0));
            const double hi = knife_edge_power(2.0 * std::sqrt(delta * bands[band] * std::pow(2.0, 1.0 / 6.0) / 343.0));
            CHECK(loss[band] >= -10.0 * std::log10(lo) - 1e-12);
            CHECK(loss[band] <= -10.0 * std::log10(hi) + 1e-12);
            previous[band] = loss[band];
        }
    }
}

TEST_CASE("the knife edge's band loss table against direct quadrature", "[screening][kirchhoff]")
{
    // The quadrature is the band's mean: against composite Simpson in
    // frequency over the band, 4000 panels.
    const double lo = std::pow(2.0, -1.0 / 6.0);
    const double hi = std::pow(2.0, 1.0 / 6.0);
    for (const double x : {1e-7, 1e-4, 0.01, 0.1, 0.5, 1.0, 3.0, 30.0, 1e3})
    {
        const int panels = 4000;
        const double h = (hi - lo) / panels;
        double sum = 0.0;
        for (int i = 0; i <= panels; i++)
        {
            const double weight = i == 0 || i == panels ? 1.0 : i % 2 == 1 ? 4.0 : 2.0;
            sum += weight * knife_edge_power(2.0 * std::sqrt(x * (lo + h * i)));
        }
        const double simpson = sum * h / 3.0 / (hi - lo);
        INFO("x " << x);
        CHECK(knife_edge_band_power(x) == Catch::Approx(simpson).epsilon(1e-10));
    }
    // The table is the quadrature: within 1e-8 dB over the range sound,
    // terrain and distance give (x = delta f_c / c from 1e-6 to 1e3), within
    // 1e-7 dB to its ends, and the quadrature itself past them.
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> typical(std::log(1e-6), std::log(1e3));
    std::uniform_real_distribution<double> whole(std::log(1e-9), std::log(1e7));
    double worst_typical = 0.0, worst = 0.0;
    for (int i = 0; i < 20000; i++)
    {
        const double a = std::exp(typical(rng));
        worst_typical = std::max(worst_typical,
                std::fabs(knife_edge_band_loss(a) + 10.0 * std::log10(knife_edge_band_power(a))));
        const double b = std::exp(whole(rng));
        worst = std::max(worst, std::fabs(knife_edge_band_loss(b) + 10.0 * std::log10(knife_edge_band_power(b))));
    }
    CHECK(worst_typical < 1e-8);
    CHECK(worst < 1e-7);
    for (const double x : {1e-9, 2e6})
        CHECK(knife_edge_band_loss(x) == -10.0 * std::log10(knife_edge_band_power(x)));
    // knife_edge_loss is the table at x = delta f_c / c.
    const std::vector<double> bands {63.0, 1000.0, 8000.0};
    std::vector<double> loss(bands.size());
    KnifeEdge edge;
    edge.delta = 0.37;
    knife_edge_loss(edge, bands, 340.0, 0, bands.size(), loss.data());
    for (std::size_t band = 0; band < bands.size(); band++)
        CHECK(loss[band] == knife_edge_band_loss(0.37 * bands[band] / 340.0));
}

namespace {

// The lowest point of the aperture plane through `deepest` that both ends see
// over every sample, by bisection along the plane's normal: the brute-force
// reading of knife_edge's construction.
double brute_force_edge(const std::vector<ScreenPoint>& p, std::size_t deepest)
{
    const ScreenPoint& r = p.front();
    const ScreenPoint& s = p.back();
    const double distance = std::hypot(s.s - r.s, s.z - r.z);
    const double ts = (s.s - r.s) / distance, tz = (s.z - r.z) / distance;
    const double a = (p[deepest].s - r.s) * ts + (p[deepest].z - r.z) * tz;
    const auto point = [&](double b) { return ScreenPoint {r.s + a * ts - b * tz, r.z + a * tz + b * ts}; };
    // The segment from u to w clears every interior sample between them.
    const auto clear = [&](const ScreenPoint& u, const ScreenPoint& w) {
        for (std::size_t k = 1; k + 1 < p.size(); k++)
        {
            const double lo = std::min(u.s, w.s), hi = std::max(u.s, w.s);
            if (p[k].s <= lo || p[k].s >= hi)
                continue;
            const double z = u.z + (w.z - u.z) * (p[k].s - u.s) / (w.s - u.s);
            if (p[k].z > z)
                return false;
        }
        return true;
    };
    double lo = 0.0, hi = 1.0;
    while (!(clear(r, point(hi)) && clear(point(hi), s)))
        hi *= 2.0;
    for (int i = 0; i < 200; i++)
    {
        const double mid = 0.5 * (lo + hi);
        if (clear(r, point(mid)) && clear(point(mid), s))
            hi = mid;
        else
            lo = mid;
    }
    return hi;
}

} // namespace

TEST_CASE("the equivalent edge is the lowest point both ends see", "[screening][kirchhoff]")
{
    // A microphone 1.2 m up, a source 300 m up 2000 m away, over profiles with
    // one sharp peak, a thick hill and random bumps: the edge from the hull's
    // horizon rays equals the brute-force bisection, and its path difference
    // is the two legs less the line of sight.
    const auto run = [](std::vector<double> heights, double step) {
        const std::size_t interior = heights.size() - 1;
        const double distance = step * static_cast<double>(interior + 1);
        std::vector<ScreenPoint> points, edges;
        std::vector<std::size_t> hull;
        Screen screen;
        if (!rubber_band(heights, interior, step, distance, 1.2, 300.0, points, hull, edges, screen))
            return false;
        const KnifeEdge edge = knife_edge(screen, points);
        std::size_t deepest = 1;
        const double ts = (points.back().s - points.front().s), tz = (points.back().z - points.front().z);
        const auto above = [&](const ScreenPoint& q) { return -(q.s - points.front().s) * tz + (q.z - points.front().z) * ts; };
        for (std::size_t k = 2; k + 1 < points.size(); k++)
            if (above(points[k]) > above(points[deepest]))
                deepest = k;
        CHECK(edge.height == Catch::Approx(brute_force_edge(points, deepest)).epsilon(1e-9).margin(1e-9));
        CHECK(edge.height > 0.0);
        CHECK(edge.d1 + edge.d2 == Catch::Approx(std::hypot(distance, 300.0 - 1.2)).epsilon(1e-12));
        const double legs = std::hypot(edge.d1, edge.height) + std::hypot(edge.d2, edge.height) - edge.d1 - edge.d2;
        CHECK(edge.delta == Catch::Approx(legs).epsilon(1e-9));
        return true;
    };
    const double step = 10.0;
    std::vector<double> peak(200, 0.0), hill(200, 0.0);
    peak[30] = 80.0;                             // one sharp edge: the edge is the peak
    for (std::size_t k = 20; k < 60; k++)        // a broad hill: two horizon rays
        hill[k] = 100.0 - 0.1 * (double(k) - 40.0) * (double(k) - 40.0);
    CHECK(run(peak, step));
    CHECK(run(hill, step));
    std::mt19937 rng(20261015);
    std::uniform_real_distribution<double> bump(0.0, 45.0);
    int blocked = 0;
    for (int trial = 0; trial < 200; trial++)
    {
        std::vector<double> heights(200, 0.0);
        for (std::size_t k = 1; k < heights.size(); k++)
            heights[k] = (k > 5 && k < 120) ? bump(rng) * std::exp(-0.02 * double(k)) : 0.0;
        blocked += run(heights, step);
    }
    CHECK(blocked > 100);
}
