#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <limits>
#include <random>
#include <vector>

#include "faddeeva.hpp"
#include "util/constants.hpp"
#include "wedge.hpp"
#include "wedge_series.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace)

namespace {

double half_plane_level(double kr, double theta, double theta0)
{
    const auto term = [&](double angle) {
        return std::exp(std::complex<double>(0.0, -kr * std::cos(angle)))
                * test::sommerfeld_f(std::sqrt(2.0 * kr) * std::cos(0.5 * angle));
    };
    return 20.0 * std::log10(std::abs(term(theta - theta0) + term(theta + theta0)));
}

} // namespace

TEST_CASE("a flat plane diffracts nothing", "[wedge]")
{
    for (double theta0 : {20.0, 70.0, 140.0})
        for (double theta = 5.0; theta < 180.0; theta += 17.0)
        {
            const double t = theta * pi / 180.0, t0 = theta0 * pi / 180.0;
            REQUIRE(std::abs(hard_wedge_coefficient(1.0, t, t0, 2.0, 10.0)) < 1e-12);
        }
}

TEST_CASE("the transition function's limits", "[wedge]")
{
    // F -> 1 far from a boundary, and sqrt(pi X) exp(-i pi/4) close to it.
    REQUIRE(std::abs(utd_transition(400.0) - std::complex<double>(1.0, 0.0)) < 2e-3);
    const double x = 1e-6;
    const auto small = utd_transition(x);
    REQUIRE(std::abs(small / (std::sqrt(pi * x) * std::polar(1.0, -pi / 4.0)) - 1.0) < 1e-2);
    // The asymptotic series and the table meet at X = 50, within the
    // series' sixth term (1e-8); F's own slope there is 2e-4.
    REQUIRE(std::abs(utd_transition(50.0 - 1e-9) - utd_transition(50.0 + 1e-9)) < 3e-8);
}

// Below X = 50 F is tabulated, a polynomial in sqrt(X) per interval: the
// Faddeeva form it was fitted to, from X = 1e-12 up, on a log and a uniform
// sweep and either side of every interval's ends (32 across sqrt(50)).
TEST_CASE("the tabulated transition function is the Faddeeva form", "[wedge]")
{
    const auto faddeeva_form = [](double x) {
        const double root = std::sqrt(x);
        return std::complex<double>(0.0, -1.0) * std::sqrt(pi) * root * std::polar(1.0, pi / 4.0)
                * faddeeva(std::polar(root, pi / 4.0));
    };
    std::vector<double> xs;
    for (double e = -12.0; e < std::log10(50.0); e += 0.001)
        xs.push_back(std::pow(10.0, e));
    for (int i = 1; i <= 50000; i++)
        xs.push_back(50.0 * i / 50000.0);
    for (int i = 1; i <= 32; i++)
    {
        const double u = std::sqrt(50.0) * i / 32.0;
        for (double d : {-1e-12, 0.0, 1e-12})
            if ((u + d) * (u + d) <= 50.0)
                xs.push_back((u + d) * (u + d));
    }
    double worst = 0.0;
    for (double x : xs)
        worst = std::max(worst, std::abs(utd_transition(x) - faddeeva_form(x)) / std::abs(faddeeva_form(x)));
    INFO("worst relative " << worst);
    REQUIRE(worst < 2e-11);
    // NaN passes through rather than indexing the table.
    const auto nan = utd_transition(std::numeric_limits<double>::quiet_NaN());
    REQUIRE(std::isnan(nan.real()));
    REQUIRE(std::isnan(nan.imag()));
}

// UTD against Sommerfeld's half-plane, across the shadow and reflection
// boundaries (the test's grid lands on some of them exactly).
TEST_CASE("the hard wedge field matches the exact half-plane", "[wedge]")
{
    double worst[3] = {0.0, 0.0, 0.0};
    const double krs[3] = {10.0, 40.0, 200.0};
    for (int k = 0; k < 3; k++)
        for (double theta0 : {50.0, 120.0, 200.0, 300.0})
            for (double theta = 0.0; theta <= 360.0; theta += 7.5)
            {
                const double t = theta * pi / 180.0, t0 = theta0 * pi / 180.0;
                const double utd = 20.0 * std::log10(std::abs(hard_wedge_field(2.0, krs[k], t, t0)));
                worst[k] = std::max(worst[k], std::abs(utd - half_plane_level(krs[k], t, t0)));
            }
    INFO("worst at kr 10, 40, 200: " << worst[0] << ", " << worst[1] << ", " << worst[2] << " dB");
    // For a half-plane under plane-wave incidence UTD's Fresnel form is
    // Sommerfeld's solution itself: 4e-7, 2e-6 and 6e-6 dB.
    REQUIRE(worst[0] < 1e-4);
    REQUIRE(worst[1] < 1e-4);
    REQUIRE(worst[2] < 1e-4);
}

// UTD against the exact series for a building's corner (n = 3/2).
TEST_CASE("the hard wedge field matches the exact right-angle wedge", "[wedge]")
{
    std::size_t i = 0;
    double worst[2] = {0.0, 0.0};
    const double krs[2] = {10.0, 40.0};
    for (int k = 0; k < 2; k++)
        for (double theta0 : {50.0, 120.0})
            for (double theta : test::directions_degrees)
            {
                const double utd = 20.0 * std::log10(std::abs(
                        hard_wedge_field(1.5, krs[k], theta * pi / 180.0, theta0 * pi / 180.0)));
                worst[k] = std::max(worst[k], std::abs(utd - test::wedge_exact[i++]));
            }
    INFO("worst at kr 10, 40: " << worst[0] << ", " << worst[1] << " dB");
    // 0.016 and 0.0045 dB.
    REQUIRE(worst[0] < 0.03);
    REQUIRE(worst[1] < 0.01);
}

// The sweep's angle-addition form against the coefficient computed directly,
// over random wedges and directions: to round-off away from the boundaries,
// and finite and close on them.
TEST_CASE("the wedge sweep is the coefficient to round-off", "[wedge]")
{
    std::mt19937 generator(7);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double worst = 0.0, worst_near = 0.0;
    for (int i = 0; i < 200000; i++)
    {
        const double n = i % 3 == 0 ? 1.5 : 1.0 + uniform(generator);
        const double theta = uniform(generator) * n * pi;
        double theta0 = uniform(generator) * n * pi;
        const bool near = i % 5 == 0;
        if (near) // within a micro-radian of a shadow or reflection boundary
            theta0 = std::abs(theta - pi) + 1e-6 * (uniform(generator) - 0.5);
        const double k = 1.0 + 200.0 * uniform(generator), r = 0.1 + 50.0 * uniform(generator);
        const std::complex<double> direct = hard_wedge_coefficient(n, theta, theta0, k, r);
        const std::complex<double> swept = hard_wedge_coefficient(WedgeSweep(n, theta).at(theta0), k, r);
        REQUIRE(std::isfinite(std::abs(swept)));
        const double difference = std::abs(swept - direct) / std::abs(direct);
        (near ? worst_near : worst) = std::max(near ? worst_near : worst, difference);
    }
    INFO("worst relative difference " << worst << ", near a boundary " << worst_near);
    // 1.3e-9; near a boundary the coefficient is a difference of large terms.
    REQUIRE(worst < 1e-8);
    REQUIRE(worst_near < 1e-5);
}

// fresnel_tail against its integral, taken down the steepest-descent path:
// from t = a = sqrt(x), t = a + e^(i pi/4) u turns exp(i t^2) into
// exp(i x) exp(-u^2 + sqrt(2) a (i - 1) u), which Simpson's rule takes on
// [0, 12].  Either side of the table's switch at x = 50, and 1/2 at 0.
TEST_CASE("the Fresnel tail is the stationary-phase integral beyond a point", "[wedge]")
{
    const auto integral = [](double x) {
        const double a = std::sqrt(2.0 * x);
        const int n = 24000;
        const double h = 12.0 / n;
        std::complex<double> sum(0.0, 0.0);
        for (int i = 0; i <= n; i++)
        {
            const double u = i * h;
            sum += (i == 0 || i == n ? 1.0 : i % 2 ? 4.0 : 2.0) * std::exp(std::complex<double>(-u * u - a * u, a * u));
        }
        return std::polar(1.0, x) * sum * (h / 3.0) / std::sqrt(pi);
    };
    REQUIRE(std::abs(fresnel_tail(0.0) - 0.5) < 1e-15);
    double worst = 0.0;
    for (double x : {1e-8, 1e-3, 0.1, 0.7, 2.0, 6.0, 20.0, 49.999, 50.001, 150.0, 800.0})
        worst = std::max(worst, std::abs(fresnel_tail(x) - integral(x)) / std::abs(integral(x)));
    INFO("worst relative difference " << worst);
    // 1.0e-8: the transition function's table.
    REQUIRE(worst < 3e-8);
}

// An edge's ends take each term's share beyond them (ended_wedge_coefficient):
// with none, D itself; far off, D; with the stationary point on one, half of
// D, which the ray past it takes there too; and at a shadow boundary the term
// that jumps keeps D's whole jump, the geometrical-optics wave's, which the
// ends do not touch.
TEST_CASE("an edge's ends take their share of each term", "[wedge]")
{
    const double n = 1.5, k = 12.0, r = 4.0;
    const WedgeAngles angles = wedge_angles(n, 2.1, 0.6);
    const std::complex<double> d = hard_wedge_coefficient(angles, k, r);
    REQUIRE(ended_wedge_coefficient(angles, k, r, EdgeEnds {}, k) == d);
    EdgeEnds far;
    far.lo = far.hi = 1e7;
    REQUIRE(std::abs(ended_wedge_coefficient(angles, k, r, far, k) / d - 1.0) < 1e-4);
    EdgeEnds on;
    on.lo = 0.0;
    REQUIRE(std::abs(ended_wedge_coefficient(angles, k, r, on, k) - 0.5 * d) < 1e-14 * std::abs(d));
    EdgeEnds past = on;
    past.past = -1;
    REQUIRE(std::abs(ended_wedge_coefficient(angles, k, r, past, k) - 0.5 * d) < 1e-14 * std::abs(d));
    // Across the source's shadow boundary, theta = theta0 + pi.
    EdgeEnds ends;
    ends.lo = 0.3;
    ends.hi = 0.8;
    const double theta0 = 0.6, epsilon = 1e-9;
    const WedgeAngles before = wedge_angles(n, theta0 + pi - epsilon, theta0);
    const WedgeAngles after = wedge_angles(n, theta0 + pi + epsilon, theta0);
    const std::complex<double> jump = hard_wedge_coefficient(after, k, r) - hard_wedge_coefficient(before, k, r);
    const std::complex<double> ended_jump
            = ended_wedge_coefficient(after, k, r, ends, k) - ended_wedge_coefficient(before, k, r, ends, k);
    REQUIRE(std::abs(jump) > 0.1 * std::abs(d));
    REQUIRE(std::abs(ended_jump - jump) < 1e-6 * std::abs(jump));
}
