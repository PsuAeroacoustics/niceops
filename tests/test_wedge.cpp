#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <vector>

#include "faddeeva.hpp"
#include "util/constants.hpp"
#include "wedge.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace)

namespace {

// Sommerfeld's exact field of a rigid half-plane lit by a plane wave, as
// tests/test_screening.cpp has it: the screen along theta = 0 (both faces),
// the wave from theta0.
std::complex<double> sommerfeld_f(double a)
{
    const std::complex<double> tail = 0.5 * std::exp(std::complex<double>(0.0, a * a))
            * faddeeva(std::polar(std::fabs(a), pi / 4.0));
    return a < 0.0 ? tail : 1.0 - tail;
}

double half_plane_level(double kr, double theta, double theta0)
{
    const auto term = [&](double angle) {
        return std::exp(std::complex<double>(0.0, -kr * std::cos(angle)))
                * sommerfeld_f(std::sqrt(2.0 * kr) * std::cos(0.5 * angle));
    };
    return 20.0 * std::log10(std::abs(term(theta - theta0) + term(theta + theta0)));
}

// The exact field of a rigid wedge with faces at 0 and 270 degrees (n = 3/2),
// 20 lg |u| for a unit plane wave, from the eigenfunction series
//   u = (2/n) sum_m eps_m exp(-i m pi / 2n) J_{m/n}(kr) cos(m theta/n) cos(m theta0/n)
// summed in Wolfram Language to 40 digits (m to 1.5 (kr + 60)), which
// reproduces the plane's incident and reflected waves exactly for n = 1.
// Order: kr in {10, 40}, theta0 in {50, 120} degrees, theta in the list below.
const std::vector<double> directions_degrees {0, 30, 60, 90, 135, 170, 180, 190, 200, 215, 230, 245, 260, 270};
const std::vector<double> wedge_exact {6.07975967, 3.23217207, 5.36821019, -6.10199301, 0.492007079, 0.660698739,
        -0.169349897, -0.0997909593, -0.121105345, -1.53458835, -3.91251894, -6.12857803, -7.46280710, -7.72382281,
        5.24715305, 0.104094917, -2.73126094, 1.48015474, 0.0884629544, -0.0852715861, -0.970003969, 0.131067558,
        1.11104573, -1.87171833, 2.16739812, -2.48987719, 3.47917479, 6.63013133, 6.22504617, 5.61474535,
        -7.92289553, 3.73607172, 2.53768510, -0.128890316, -0.190481530, 0.0100532322, -0.268604562, 0.296269637,
        -4.96236734, -10.3178652, -12.8820876, -13.3161881, 6.18507161, -14.8988505, -5.20667036, 0.900786106,
        -0.283589287, -0.384043060, -0.219059255, -0.272030345, -0.565647275, 0.847842077, 0.236661492, -2.01113085,
        6.55169491, 5.04557179};

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
    // The asymptotic series and the Faddeeva form meet at X = 50, within the
    // series' sixth term (1e-8); F's own slope there is 2e-4.
    REQUIRE(std::abs(utd_transition(50.0 - 1e-9) - utd_transition(50.0 + 1e-9)) < 3e-8);
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
            for (double theta : directions_degrees)
            {
                const double utd = 20.0 * std::log10(std::abs(
                        hard_wedge_field(1.5, krs[k], theta * pi / 180.0, theta0 * pi / 180.0)));
                worst[k] = std::max(worst[k], std::abs(utd - wedge_exact[i++]));
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
