#pragma once

// Reference values shared by the wedge tests.

#include <cmath>
#include <complex>
#include <vector>

#include "edge_source.hpp"
#include "faddeeva.hpp"
#include "util/constants.hpp"

namespace niceops::test {

// Sommerfeld's exact field of a rigid half-plane lit by a plane wave, as
// tests/test_screening.cpp has it: the screen along theta = 0 (both faces),
// the wave from theta0.
inline std::complex<double> sommerfeld_f(double a)
{
    const std::complex<double> tail = 0.5 * std::exp(std::complex<double>(0.0, a * a))
            * faddeeva(std::polar(std::fabs(a), pi / 4.0));
    return a < 0.0 ? tail : 1.0 - tail;
}

// The exact field of a rigid wedge with faces at 0 and 270 degrees (n = 3/2),
// 20 lg |u| for a unit plane wave, from the eigenfunction series
//   u = (2/n) sum_m eps_m exp(-i m pi / 2n) J_{m/n}(kr) cos(m theta/n) cos(m theta0/n)
// summed in Wolfram Language to 40 digits (m to 1.5 (kr + 60)), which
// reproduces the plane's incident and reflected waves exactly for n = 1.
// Order: kr in {10, 40}, theta0 in {50, 120} degrees, theta in the list below.
inline const std::vector<double> directions_degrees {0, 30, 60, 90, 135, 170, 180, 190, 200, 215, 230, 245, 260, 270};
inline const std::vector<double> wedge_exact {6.07975967, 3.23217207, 5.36821019, -6.10199301, 0.492007079, 0.660698739,
        -0.169349897, -0.0997909593, -0.121105345, -1.53458835, -3.91251894, -6.12857803, -7.46280710, -7.72382281,
        5.24715305, 0.104094917, -2.73126094, 1.48015474, 0.0884629544, -0.0852715861, -0.970003969, 0.131067558,
        1.11104573, -1.87171833, 2.16739812, -2.48987719, 3.47917479, 6.63013133, 6.22504617, 5.61474535,
        -7.92289553, 3.73607172, 2.53768510, -0.128890316, -0.190481530, 0.0100532322, -0.268604562, 0.296269637,
        -4.96236734, -10.3178652, -12.8820876, -13.3161881, 6.18507161, -14.8988505, -5.20667036, 0.900786106,
        -0.283589287, -0.384043060, -0.219059255, -0.272030345, -0.565647275, 0.847842077, 0.236661492, -2.01113085,
        6.55169491, 5.04557179};

// The geometrical-optics field of a unit point source about a rigid wedge of
// open angle n pi, the two ends as edge_frame_ends gives them: the incident
// wave where the wedge leaves it, and either face's reflection where the face
// sends one.  On a boundary, to round-off, the wave is there:
// edge_source_field takes its limit from that side.
inline std::complex<double> wedge_geometrical(double n, const EdgeFrameEnds& e, double k)
{
    std::complex<double> sum = 0.0;
    const auto wave = [&](double theta) {
        const double d = std::sqrt(e.r_s * e.r_s + e.r_r * e.r_r - 2.0 * e.r_s * e.r_r * std::cos(theta - e.theta_r)
                + (e.z_s - e.z_r) * (e.z_s - e.z_r));
        return std::polar(1.0 / d, k * d);
    };
    const double slack = 1e-12;
    if (std::abs(e.theta_r - e.theta_s) < pi + slack)
        sum += wave(e.theta_s);
    if (e.theta_r + e.theta_s < pi + slack)
        sum += wave(-e.theta_s);
    if (e.theta_r + e.theta_s > (2.0 * n - 1.0) * pi - slack)
        sum += wave(2.0 * n * pi - e.theta_s);
    return sum;
}

} // namespace niceops::test
