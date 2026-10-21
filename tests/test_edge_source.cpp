#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <random>
#include <vector>

#include "buildings.hpp"
#include "edge_source.hpp"
#include "point_paths.hpp"
#include "urban.hpp"
#include "util/constants.hpp"
#include "wedge.hpp"
#include "wedge_series.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace)

namespace {

using complex = std::complex<double>;

// A point about the edge, in its own frame.
struct Place
{
    double r, theta, z;
};

EdgeFrameEnds ends_of(const Place& s, const Place& r)
{
    return {s.r, s.theta, s.z, r.r, r.theta, r.z};
}

// The geometrical-optics field about the wedge (wedge_series.hpp).
complex geometrical(double n, const Place& s, const Place& r, double k)
{
    return test::wedge_geometrical(n, ends_of(s, r), k);
}

complex total(double n, const Place& s, const Place& r, double k, double z_lo = -1e15, double z_hi = 1e15)
{
    return geometrical(n, s, r, k) + edge_source_field(n, ends_of(s, r), k, z_lo, z_hi);
}

double level(const complex& p)
{
    return 20.0 * std::log10(std::abs(p));
}

// The edge-source integral by brute force: Gauss-Legendre along the edge in
// z, panels a fifth of a radian of phase or less, over [z_lo, z_hi].
complex direct(double n, const Place& s, const Place& r, double k, double z_lo, double z_hi)
{
    static const double x[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
            -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
    static const double w[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
            0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
    const double nu = 1.0 / n;
    const double phis[4] = {pi + s.theta + r.theta, pi + s.theta - r.theta, pi - s.theta + r.theta,
            pi - s.theta - r.theta};
    const auto integrand = [&](double z) {
        const double m = std::hypot(s.r, z - s.z), l = std::hypot(r.r, z - r.z);
        const double ch = (m * l + (z - s.z) * (z - r.z)) / (s.r * r.r);
        const double eta = std::acosh(std::max(ch, 1.0));
        double beta = 0.0;
        for (double phi : phis)
            beta += std::sin(nu * phi) / (std::cosh(nu * eta) - std::cos(nu * phi));
        return beta * std::polar(1.0 / (m * l), k * (m + l));
    };
    const int panels = static_cast<int>(std::ceil(k * (z_hi - z_lo) / 0.2)) + 200;
    const double h = (z_hi - z_lo) / panels;
    complex sum = 0.0;
    for (int p = 0; p < panels; p++)
        for (int q = 0; q < 8; q++)
            sum += w[q] * integrand(z_lo + h * (p + 0.5 * (x[q] + 1.0)));
    return -nu / (4.0 * pi) * 0.5 * h * sum;
}

} // namespace

// The field of a unit plane wave from theta0 at k r = kr from the edge, as
// the limit of a point source's field over its own at the edge as the source
// recedes along theta0: taken at 1e8 and 2e8 and extrapolated, which leaves
// (k r^2 / 2 R)^2 of the curvature and none of a far source's round-off.
complex plane_wave_limit(double n, double kr, double theta0, double theta)
{
    const auto at = [&](double far) {
        const Place s {far, theta0 * pi / 180.0, 0.0}, r {kr, theta * pi / 180.0, 0.0};
        return total(n, s, r, 1.0) * far * std::polar(1.0, -far);
    };
    return 2.0 * at(2e8) - at(1e8);
}

// The plane-wave field of a rigid half-plane, Sommerfeld's.
complex half_plane(double kr, double theta, double theta0)
{
    const auto term = [&](double angle) {
        return std::polar(1.0, -kr * std::cos(angle)) * test::sommerfeld_f(std::sqrt(2.0 * kr) * std::cos(0.5 * angle));
    };
    return term(theta - theta0) + term(theta + theta0);
}

TEST_CASE("far from its source the edge source is Sommerfeld's half-plane", "[edge_source]")
{
    // A point source far out along theta0 is a plane wave at the edge; the
    // field over the source's own is the plane wave's.  The grid lands on
    // shadow and reflection boundaries.
    double worst = 0.0;
    for (double kr : {0.5, 2.0, 10.0, 40.0, 200.0})
        for (double theta0 : {50.0, 120.0, 200.0, 300.0})
            for (double theta = 0.0; theta <= 360.0; theta += 7.5)
            {
                const complex u = plane_wave_limit(2.0, kr, theta0, theta);
                const complex exact = half_plane(kr, theta * pi / 180.0, theta0 * pi / 180.0);
                const double error = std::abs(u - exact) / std::abs(exact);
                worst = std::max(worst, error);
            }
    INFO("worst relative error " << worst);
    // 8e-7, at k r = 40.
    REQUIRE(worst < 2e-6);
}

TEST_CASE("far from its source the edge source is the exact right-angle wedge", "[edge_source]")
{
    // A point source far out along theta0 is a plane wave at the edge; the
    // field over the source's own is the plane wave's.
    std::size_t i = 0;
    double worst = 0.0;
    for (double kr : {10.0, 40.0})
        for (double theta0 : {50.0, 120.0})
            for (double theta : test::directions_degrees)
            {
                const complex u = plane_wave_limit(1.5, kr, theta0, theta);
                const double error = std::abs(level(u) - test::wedge_exact[i++]);
                worst = std::max(worst, error);
            }
    INFO("worst " << worst << " dB");
    // 1.5e-6 dB, the series' own nine figures being 1e-8.
    REQUIRE(worst < 5e-6);
}

TEST_CASE("the edge source's descent is the direct integral along the edge", "[edge_source]")
{
    std::mt19937 random(11);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    double worst = 0.0;
    for (int trial = 0; trial < 200; trial++)
    {
        const double n = trial % 2 == 0 ? 1.5 : 1.0 + unit(random);
        const Place s {0.5 + 10.0 * unit(random), n * pi * unit(random), 10.0 * unit(random)};
        const Place r {0.5 + 10.0 * unit(random), n * pi * unit(random), 10.0 * unit(random)};
        const double z_lo = 10.0 * (unit(random) - 0.3), z_hi = z_lo + 0.5 + 10.0 * unit(random);
        for (double k : {0.5, 5.0, 40.0})
        {
            const complex a = edge_source_field(n, ends_of(s, r), k, z_lo, z_hi);
            const complex b = direct(n, s, r, k, z_lo, z_hi);
            const double error = std::abs(a - b) / std::max(std::abs(b), 1e-3 * std::abs(geometrical(n, s, r, k)));
            worst = std::max(worst, error);
        }
    }
    INFO("worst relative error " << worst);
    // 1.4e-8.
    REQUIRE(worst < 1e-6);
}


// The ray the codebase's UTD takes for the same source and receiver
// (point_paths.cpp): diffracted at the apex, the coefficient at k sin(beta)
// and distance L sin(beta).
complex utd(double n, const Place& s, const Place& r, double k)
{
    const double apex = s.z + (r.z - s.z) * s.r / (s.r + r.r);
    const double s_in = std::hypot(s.r, apex - s.z), s_out = std::hypot(r.r, r.z - apex);
    const double sin_beta = s.r / s_in;
    return std::polar(1.0 / std::sqrt(s_in * s_out * (s_in + s_out)), k * (s_in + s_out))
            * hard_wedge_coefficient(n, r.theta, s.theta, k * sin_beta, s_in * s_out * sin_beta / (s_in + s_out))
            / std::sqrt(sin_beta);
}

TEST_CASE("at high frequency the edge source is UTD", "[edge_source]")
{
    // Away from the boundaries UTD is the integral's stationary point, to
    // O(1 / kL); near them its transition function holds the same uniform
    // limit.  So the total fields agree, closer as k grows.
    std::mt19937 random(5);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    double worst[3] = {0.0, 0.0, 0.0};
    const double ks[3] = {5.0, 50.0, 500.0};
    for (int trial = 0; trial < 300; trial++)
    {
        const double n = trial % 2 == 0 ? 1.5 : 1.0 + unit(random);
        const Place s {1.0 + 10.0 * unit(random), n * pi * unit(random), 10.0 * unit(random)};
        Place r {1.0 + 10.0 * unit(random), n * pi * unit(random), 10.0 * unit(random)};
        if (trial % 3 == 0)
        {
            const double boundary = trial % 2 == 0 ? s.theta + pi : pi - s.theta;
            if (boundary > 0.0 && boundary < n * pi)
                r.theta = boundary + 0.05 * (unit(random) - 0.5);
        }
        for (int i = 0; i < 3; i++)
        {
            const complex go = geometrical(n, s, r, ks[i]);
            const double btm = level(go + edge_source_field(n, ends_of(s, r), ks[i], -1e15, 1e15));
            worst[i] = std::max(worst[i], std::abs(btm - level(go + utd(n, s, r, ks[i]))));
        }
    }
    INFO("worst at k 5, 50, 500: " << worst[0] << ", " << worst[1] << ", " << worst[2] << " dB");
    REQUIRE(worst[1] < worst[0]);
    REQUIRE(worst[2] < worst[1]);
    REQUIRE(worst[2] < 0.01);
}

TEST_CASE("the edge source is reciprocal", "[edge_source]")
{
    std::mt19937 random(3);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    double worst = 0.0;
    for (int trial = 0; trial < 100; trial++)
    {
        const double n = 1.0 + unit(random);
        const Place s {0.5 + 10.0 * unit(random), n * pi * unit(random), 10.0 * unit(random)};
        const Place r {0.5 + 10.0 * unit(random), n * pi * unit(random), 10.0 * unit(random)};
        const double z_lo = 10.0 * (unit(random) - 0.3), z_hi = z_lo + 0.5 + 10.0 * unit(random);
        for (double k : {0.5, 5.0, 40.0})
        {
            const complex a = edge_source_field(n, ends_of(s, r), k, z_lo, z_hi);
            const complex b = edge_source_field(n, ends_of(r, s), k, z_lo, z_hi);
            worst = std::max(worst, std::abs(a - b) / std::abs(a));
        }
    }
    INFO("worst relative difference " << worst);
    REQUIRE(worst < 1e-8);
}

TEST_CASE("the edge source is continuous across a boundary and an end", "[edge_source]")
{
    // Across the shadow and either reflection boundary the geometrical-optics
    // wave jumps by its whole amplitude, and the diffracted wave makes it up;
    // as the apex crosses an end of the edge the integral loses half its
    // stationary point smoothly.
    const double n = 1.5, k = 20.0, step = 1e-9;
    const Place s {4.0, 0.6, 1.0};
    double worst = 0.0;
    for (double boundary : {s.theta + pi, pi - s.theta, (2.0 * n - 1.0) * pi - s.theta})
        for (double side : {-1.0, 1.0})
        {
            const Place a {6.0, boundary - side * step, 3.0}, b {6.0, boundary + side * step, 3.0};
            const double jump = std::abs(total(n, s, a, k, -5.0, 8.0) - total(n, s, b, k, -5.0, 8.0))
                    / std::abs(total(n, s, a, k, -5.0, 8.0));
            worst = std::max(worst, jump);
        }
    INFO("worst jump across a boundary " << worst);
    // 8e-8 and 4e-8 below: the field's own slope over the 2e-9 step.
    REQUIRE(worst < 1e-6);

    // The receiver moved along the edge so that the apex crosses z_hi, the
    // geometrical-optics field held still.
    const double z_hi = 8.0;
    const double z_r = s.z + (z_hi - s.z) * (s.r + 6.0) / s.r;
    double end_jump = 0.0;
    for (double theta : {2.0, 3.5, 4.4})
    {
        const Place a {6.0, theta, z_r - step}, b {6.0, theta, z_r + step};
        const complex fa = edge_source_field(n, ends_of(s, a), k, -5.0, z_hi);
        const complex fb = edge_source_field(n, ends_of(s, b), k, -5.0, z_hi);
        end_jump = std::max(end_jump, std::abs(fa - fb) / std::abs(fa));
    }
    INFO("worst jump at the end " << end_jump);
    REQUIRE(end_jump < 1e-6);
}


TEST_CASE("the edge-source spectrum is the field at each wavenumber", "[edge_source]")
{
    // 72 wavenumbers, third-octave sub-bands from 50 Hz to 10 kHz.
    std::vector<double> ks;
    for (int i = 0; i < 72; i++)
        ks.push_back(2.0 * pi * 50.0 * std::pow(2.0, i / 9.0) / 343.0);
    std::mt19937 random(13);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    double worst = 0.0;
    for (int trial = 0; trial < 300; trial++)
    {
        const double n = trial % 2 == 0 ? 1.5 : 1.0 + unit(random);
        const Place s {0.3 + 20.0 * unit(random), n * pi * unit(random), 20.0 * unit(random)};
        Place r {0.3 + 20.0 * unit(random), n * pi * unit(random), 20.0 * unit(random)};
        if (trial % 3 == 0)
        {
            const double boundary = trial % 2 == 0 ? s.theta + pi : pi - s.theta;
            if (boundary > 0.0 && boundary < n * pi)
                r.theta = boundary + 0.01 * (unit(random) - 0.5);
        }
        const double z_lo = 20.0 * (unit(random) - 0.3), z_hi = z_lo + 0.2 + 30.0 * unit(random);
        const std::vector<complex> spectrum = edge_source_spectrum(n, ends_of(s, r), z_lo, z_hi, ks);
        for (std::size_t i = 0; i < ks.size(); i++)
        {
            const complex exact = edge_source_field(n, ends_of(s, r), ks[i], z_lo, z_hi);
            const double scale = std::abs(exact) + 1e-3 * std::abs(geometrical(n, s, r, ks[i]));
            worst = std::max(worst, std::abs(spectrum[i] - exact) / scale);
        }
    }
    INFO("worst relative error " << worst);
    // 4.4e-5, 0.0004 dB.
    REQUIRE(worst < 1e-4);
}

TEST_CASE("away from transitions the twice-diffracted cascade is the edge-source integral", "[edge_source]")
{
    // Two screens: the first's top along x at z = 4; the second's edge through
    // (5, 5, 3) in the plane y = 5, at alpha from the vertical (90 degrees:
    // parallel to the first), the screen on its +x / below side.  The source
    // and the receiver low on either side, both edges well out of their
    // transitions.  The cascade's spreading for edges that are not parallel
    // was once the parallel one, 1.45 dB high on the perpendicular pair.
    Edge e1;
    e1.start = {-50.0, 0.0, 4.0};
    e1.direction = {1.0, 0.0, 0.0};
    e1.length = 100.0;
    e1.face0 = {0.0, 0.0, -1.0};
    e1.air = {0.0, -1.0, 0.0};
    e1.n = 2.0;
    const Vec3 s {9.0, -4.0, 1.0}, r {9.0, 11.0, 1.0};
    double worst = 0.0;
    for (double alpha_degrees : {0.0, 45.0, 90.0})
    {
        const double alpha = alpha_degrees * pi / 180.0;
        Edge e2;
        e2.direction = {std::sin(alpha), 0.0, std::cos(alpha)};
        e2.start = Vec3(5.0, 5.0, 3.0) - 50.0 * e2.direction;
        e2.length = 100.0;
        e2.face0 = {std::cos(alpha), 0.0, -std::sin(alpha)};
        e2.air = {0.0, -1.0, 0.0};
        e2.n = 2.0;
        // The stationary point, by alternating golden-section searches.
        double z1 = 50.0, z2 = 50.0;
        const auto path = [&](double a, double b) {
            const Vec3 p = e1.start + a * e1.direction, q = e2.start + b * e2.direction;
            return (s - p).norm() + (q - p).norm() + (r - q).norm();
        };
        for (int iteration = 0; iteration < 100; iteration++)
            for (int which = 0; which < 2; which++)
            {
                double lo = 0.0, hi = 100.0;
                for (int i = 0; i < 100; i++)
                {
                    const double c = lo + 0.381966 * (hi - lo), d = hi - 0.381966 * (hi - lo);
                    const bool left = which == 0 ? path(c, z2) < path(d, z2) : path(z1, c) < path(z1, d);
                    (left ? hi : lo) = left ? d : c;
                }
                (which == 0 ? z1 : z2) = 0.5 * (lo + hi);
            }
        const PointRay ray = twice_diffracted_ray(e1, e2, s, r, e1.start + z1 * e1.direction,
                e2.start + z2 * e2.direction);
        const std::vector<double> ks {50.0, 100.0};
        const std::vector<complex> spectrum = localized_twice_diffracted_spectrum(e1, e2, s, r, 1.0, ks,
                PartnerModel::keller);
        for (std::size_t j = 0; j < ks.size(); j++)
        {
            const double k = ks[j];
            const complex cascade = point_ray_amplitude(ray, ray.spreading, k) * std::polar(1.0, k * ray.length);
            const complex exact = localized_twice_diffracted_field(e1, e2, s, r, k, 1.0, 80.0, 1e-6);
            // The spectrum takes the edges' whole lengths, 50 m each way, where
            // the window takes a few Fresnel zones: their far ends' waves
            // differ, by a few thousandths.
            INFO("spectrum " << level(spectrum[j]) << " dB");
            REQUIRE(std::abs(spectrum[j] - exact) < 0.01 * std::abs(exact));
            const double error = std::abs(level(cascade) - level(exact));
            INFO("alpha " << alpha_degrees << " k " << k << ": cascade " << level(cascade) << " dB, edge source "
                          << level(exact) << " dB, phase " << std::arg(cascade / exact));
            worst = std::max(worst, error);
            REQUIRE(std::abs(std::arg(cascade / exact)) < 0.01);
        }
    }
    INFO("worst " << worst << " dB");
    // 0.007 dB.
    REQUIRE(worst < 0.03);
}

namespace {

// An edge along +x at height h and across-position y, its faces running
// down (face0) and, n pi on, toward +y (n = 1.5, a block's roof edge) or
// down again (n = 2, a screen).
Edge edge_along_x(double y, double h, double n, double length, bool face_toward_minus_y = false)
{
    Edge e;
    e.start = {-0.5 * length, y, h};
    e.direction = {1.0, 0.0, 0.0};
    e.length = length;
    e.face0 = {0.0, 0.0, -1.0};
    e.air = {0.0, face_toward_minus_y ? 1.0 : -1.0, 0.0};
    e.n = n;
    return e;
}

} // namespace

TEST_CASE("the level-curve spectrum is the nested integral", "[edge_source]")
{
    // Two screens (parallel edges), a block's roof (a shared face, the plane
    // factor a half) and a screen's top with a vertical screen's edge
    // (perpendicular, whose boundaries cross the other edge: the level curves'
    // steps), on edges a few metres long, so their ends count.
    struct Pair
    {
        const char* name;
        Edge e1, e2;
        Vec3 s, r;
        double plane;
    };
    std::vector<Pair> pairs;
    pairs.push_back({"screens", edge_along_x(0.0, 4.0, 2.0, 10.0), edge_along_x(5.0, 4.0, 2.0, 10.0),
            {-2.0, -6.0, 1.5}, {3.0, 13.0, 1.0}, 1.0});
    pairs.push_back({"roof", edge_along_x(0.0, 5.18, 1.5, 6.0), edge_along_x(2.44, 5.18, 1.5, 6.0, true),
            {0.5, -3.0, 3.18}, {-1.0, 5.44, 3.18}, 0.5});
    Edge top;
    top.start = {-5.0, 0.0, 4.0};
    top.direction = {1.0, 0.0, 0.0};
    top.length = 20.0;
    top.face0 = {0.0, 0.0, -1.0};
    top.air = {0.0, -1.0, 0.0};
    top.n = 2.0;
    Edge side;
    side.start = {5.0, 5.0, -5.0};
    side.direction = {0.0, 0.0, 1.0};
    side.length = 15.0;
    side.face0 = {1.0, 0.0, 0.0};
    side.air = {0.0, -1.0, 0.0};
    side.n = 2.0;
    pairs.push_back({"perpendicular", top, side, {9.0, -4.0, 1.0}, {9.0, 11.0, 3.0}, 1.0});
    const std::vector<double> ks {1.0, 10.0};
    for (const Pair& p : pairs)
    {
        const std::vector<complex> spectrum = twice_diffracted_spectrum(p.e1, p.e2, p.s, p.r, {0.0, p.e1.length},
                {0.0, p.e2.length}, p.plane, ks);
        for (std::size_t i = 0; i < ks.size(); i++)
        {
            const complex nested = twice_diffracted_field(p.e1, p.e2, p.s, p.r, ks[i], {0.0, p.e1.length},
                    {0.0, p.e2.length}, p.plane, 1e-7);
            INFO(p.name << " at k " << ks[i] << ": " << level(spectrum[i]) << " dB against " << level(nested));
            // Within 1.3e-4 at the integral's own tolerances (1e-3 and 1e-4,
            // edge_source.cpp), 0.001 dB; 4e-6 a hundred times tighter.
            REQUIRE(std::abs(spectrum[i] - nested) < 5e-4 * std::abs(nested));
        }
    }
}



namespace {

Footprint turned_box(double x0, double y0, double x1, double y1, double top, double degrees)
{
    const double c = std::cos(degrees * pi / 180.0), s = std::sin(degrees * pi / 180.0);
    Footprint f;
    for (const auto& [x, y] : {std::make_pair(x0, y0), std::make_pair(x1, y0), std::make_pair(x1, y1),
                 std::make_pair(x0, y1)})
    {
        f.x.push_back(c * x - s * y);
        f.y.push_back(s * x + c * y);
    }
    f.top = top;
    return f;
}

} // namespace

// A district turned about the vertical, its source and microphone with it, is
// the same field under the edge-source integral: a point in an edge's face's
// plane is on the face whatever the rounding, where the partner search once
// put it at 0 or just under 2 pi as it went, and moved bands by up to 1 dB.
TEST_CASE("the edge source is the same with the district turned", "[edge_source]")
{
    std::vector<double> bands;
    for (int i = 0; i <= 19; i++)
        bands.push_back(125.0 * std::pow(2.0, i / 3.0));
    std::vector<double> first;
    for (double degrees : {0.0, 17.0, 61.3})
    {
        const double c = std::cos(degrees * pi / 180.0), s = std::sin(degrees * pi / 180.0);
        const auto turn = [&](const Vec3& p) { return Vec3(c * p.x - s * p.y, s * p.x + c * p.y, p.z); };
        const Buildings city({turned_box(0.0, 0.0, 0.96, 1.0, 0.492, degrees)});
        PointPathOptions options;
        options.edge_source = true;
        const PointSourceField field(city, turn(Vec3(-2.0, 0.5, 1.5)), turn(Vec3(1.46, 0.5, 0.0032)), options);
        std::vector<double> reflection(city.facades().size() * bands.size(), 1.0), out(bands.size());
        field.band_levels(reflection, 1.0, bands, 343.2, out.data());
        if (first.empty())
            first = out;
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            INFO(degrees << " degrees, " << bands[b] << " Hz");
            REQUIRE(std::abs(out[b] - first[b]) < 0.005);
        }
    }
}

// A partner's part along an edge that runs to a vertex the two edges share
// reaches the end exactly: the search skips the shared vertex and bisected to
// a hair short of it, so a wave over the whole edge was taken as a partner and
// added in whole, one ray 10-16 dB high.
TEST_CASE("a partner's part reaches a shared vertex", "[edge_source]")
{
    Footprint a, b;
    a.x = {0.0, 10.0, 10.0, 0.0};
    a.y = {0.0, 0.0, 10.0, 10.0};
    a.top = 12.0;
    b.x = {20.0, 30.0, 30.0, 20.0};
    b.y = {3.0, 3.0, 13.0, 13.0};
    b.top = 8.0;
    const Buildings city({a, b});
    PointPathOptions options;
    options.edge_source = true;
    const PointSourceField field(city, Vec3(-40.0, 5.0, 25.0), Vec3(19.0, 1.0, 7.0), options);
    std::size_t twice = 0, parts = 0;
    for (const PointRay& ray : field.rays())
    {
        if (ray.wedges != 2)
            continue;
        twice++;
        const Edge g1 = continued(city.edges()[ray.edge_source.edges[0]]);
        const Edge g2 = continued(city.edges()[ray.edge_source.edges[1]]);
        for (const PartnerWave& p : twice_diffracted_partners(g1, g2, ray.edge_source.source, ray.edge_source.microphone))
        {
            const double length = (p.over_second ? g2 : g1).length;
            parts++;
            INFO("part " << p.part.lo << " to " << p.part.hi << " of " << length);
            REQUIRE((p.part.lo == 0.0 || p.part.lo > 1e-6 * length));
            REQUIRE((p.part.hi == length || p.part.hi < length * (1.0 - 1e-6)));
        }
    }
    REQUIRE(twice > 0);
    REQUIRE(parts > 0);
}

// Image theory: a low block whose two roof edges run into a tall wall that
// fills the half-space y < 0.  The wall is a rigid plane, so the field is the
// free field of the block doubled through it, radiated by the source and by
// its mirror together.  The model's rays -- over the edges, off the wall before
// the first edge, between the edges (the second mirrored), after the second --
// must add up to the doubled block's.  The edges' ends at the wall are no
// ends: each edge goes on in its image there (continued), so a ray takes the
// doubled edge's far ends, and its mirror past the wall is the same ray and is
// not traced.  By UTD, ends and all, and by the edge source the rays add up to
// round-off (0.02 dB measured); when the ends at the wall were ends they
// missed by 1.3 dB (UTD) and 1.4 dB (the edge source).  Without the rays that
// reflect between the edges the sum misses by over 3 dB.
TEST_CASE("a wall the edges run into is a mirror, by UTD and by the edge source", "[edge_source]")
{
    const std::vector<double> bands {125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0};
    const double c = 343.2;
    Footprint wall;
    wall.x = {-400.0, 400.0, 400.0, -400.0};
    wall.y = {-400.0, -400.0, 0.0, 0.0};
    wall.top = 300.0;
    Footprint low;
    low.x = {0.0, 10.0, 10.0, 0.0};
    low.y = {0.0, 0.0, 30.0, 30.0};
    low.top = 5.0;
    Footprint doubled = low;
    doubled.y = {-30.0, -30.0, 30.0, 30.0};
    const Buildings with_wall({wall, low}), open({doubled});
    const Vec3 microphone(20.0, 3.0, 1.2);
    // The worst difference, dB, between the wall's field and the doubled
    // block's, the source `ys` from the wall.
    auto worst_with = [&](bool edge_source, std::size_t between, double ys) {
        PointPathOptions options;
        options.edge_source = edge_source;
        options.max_reflections = 1;
        options.max_edge_reflections = 1;
        options.max_twice_reflections = 1;
        options.max_between_reflections = between;
        options.triple_diffraction = PointPathOptions::Triple::none;
        using Problems = std::vector<std::vector<std::pair<std::complex<double>, double>>>;
        auto gather = [&](const Buildings& city, const Vec3& s, Problems& into) {
            const PointSourceField field(city, s, microphone, options);
            const std::vector<double> reflection(city.facades().size() * bands.size(), 1.0);
            into.resize(3 * bands.size());
            field.sub_band_problems(reflection, 1.0, bands, c,
                    [&](std::size_t band, int part, const std::vector<std::complex<double>>& amplitude,
                            const std::vector<double>& length, const SubBand&) {
                        for (std::size_t i = 0; i < amplitude.size(); i++)
                            into[3 * band + part].push_back({amplitude[i], length[i]});
                    });
        };
        Problems mirrored, unfolded;
        gather(with_wall, Vec3(-10.0, ys, 1.5), mirrored);
        gather(open, Vec3(-10.0, ys, 1.5), unfolded);
        gather(open, Vec3(-10.0, -ys, 1.5), unfolded);
        double worst = 0.0;
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            const auto parts = sub_bands(bands[b], c);
            double energy[2][3];
            for (int which = 0; which < 2; which++)
                for (int part = 0; part < 3; part++)
                {
                    auto rays = (which == 0 ? mirrored : unfolded)[3 * b + part];
                    std::sort(rays.begin(), rays.end(),
                            [](const auto& x, const auto& y) { return x.second < y.second; });
                    std::vector<std::complex<double>> amplitude;
                    std::vector<double> length;
                    for (const auto& [a, l] : rays)
                    {
                        amplitude.push_back(a);
                        length.push_back(l);
                    }
                    energy[which][part] = band_averaged_energy(amplitude, length, parts[part].k1, parts[part].k2);
                }
            worst = std::max(worst,
                    std::abs(band_level(bands[b], parts, energy[0]) - band_level(bands[b], parts, energy[1])));
        }
        return worst;
    };
    for (double ys : {2.0, 8.0})
    {
        INFO("the source " << ys << " m from the wall");
        REQUIRE(worst_with(false, 1, ys) < 0.05);
        REQUIRE(worst_with(true, 1, ys) < 0.05);
    }
    REQUIRE(worst_with(true, 0, 2.0) > 3.0);
}
