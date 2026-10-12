#include "wedge.hpp"

#include <cmath>

#include "faddeeva.hpp"
#include "util/constants.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

// The two terms of one beta (theta - theta0 or theta + theta0):
//   cot((pi + beta) / 2n) F(k L a+(beta)) + cot((pi - beta) / 2n) F(k L a-(beta)),
// a+-(beta) = 2 cos^2((2 n pi N+- - beta) / 2), N+- the integer nearest
// (beta +- pi) / (2 n pi).  On a boundary a cotangent's pole meets F's zero;
// the product's limit is taken a hair off it.
complex pair(double n, double beta, double kl)
{
    auto one = [&](int s) {
        double b = beta;
        if (std::abs(std::sin((pi + s * b) / (2.0 * n))) < 1e-12)
            b += 1e-9;
        const double nearest = std::round((b + s * pi) / (2.0 * n * pi));
        const double c = std::cos((2.0 * n * pi * nearest - b) / 2.0);
        return (1.0 / std::tan((pi + s * b) / (2.0 * n))) * utd_transition(kl * 2.0 * c * c);
    };
    return one(+1) + one(-1);
}

} // namespace

complex utd_transition(double x)
{
    if (x <= 0.0)
        return complex(0.0, 0.0);
    // Past X = 50 the asymptotic series, sum over m of (-i)^m (2m - 1)!! / (2X)^m
    // (Kouyoumjian and Pathak's, conjugated), to its sixth term, which leaves
    // 1e-8; it keeps w's argument within its domain.
    if (x > 50.0)
    {
        const double i1 = 1.0 / x, i2 = i1 * i1;
        return complex(1.0 - 0.75 * i2 + 105.0 / 16.0 * i2 * i2,
                -0.5 * i1 + 15.0 / 8.0 * i2 * i1 - 945.0 / 32.0 * i2 * i2 * i1);
    }
    const double root = std::sqrt(x);
    const complex z = std::polar(root, pi / 4.0);
    return complex(0.0, -1.0) * std::sqrt(pi) * root * std::polar(1.0, pi / 4.0) * faddeeva(z);
}

complex hard_wedge_coefficient(double n, double theta, double theta0, double k, double r)
{
    const double kl = k * r;
    const double minus = theta - theta0;
    const double plus = theta + theta0;
    // A hard wedge adds the incident (theta - theta0) and the reflected
    // (theta + theta0) pairs.
    const complex sum = pair(n, minus, kl) + pair(n, plus, kl);
    return -std::polar(1.0, pi / 4.0) / (2.0 * n * std::sqrt(2.0 * pi * k)) * sum;
}

complex hard_wedge_field(double n, double kr, double theta, double theta0)
{
    // On a shadow or reflection boundary the geometrical-optics test and the
    // coefficient's limit must take the same side; a direction within a hair
    // of one is moved off it once, for both.
    for (double boundary : {theta0 + pi, theta0 - pi, pi - theta0, (2.0 * n - 1.0) * pi - theta0})
        if (std::abs(theta - boundary) < 1e-9)
            theta = boundary + 1e-8;
    complex go(0.0, 0.0);
    auto wave = [&](double from) { return std::polar(1.0, -kr * std::cos(theta - from)); };
    // The incident wave, wherever the wedge does not shadow it.
    if (std::abs(theta - theta0) < pi)
        go += wave(theta0);
    // The reflection off face 0, for a source that sees it.
    if (theta0 < pi && theta < pi - theta0)
        go += wave(-theta0);
    // The reflection off face n.
    if (theta0 > (n - 1.0) * pi && theta > (2.0 * n - 1.0) * pi - theta0)
        go += wave(2.0 * n * pi - theta0);
    // k = 1, r = kr: the coefficient scales as 1/sqrt(k) and the spreading
    // as 1/sqrt(r), so only k r matters.
    return go + hard_wedge_coefficient(n, theta, theta0, 1.0, kr) * std::polar(1.0, kr) / std::sqrt(kr);
}

} // namespace niceops
