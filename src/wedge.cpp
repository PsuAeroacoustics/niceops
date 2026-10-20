#include "wedge.hpp"

#include <algorithm>
#include <cmath>

#include "faddeeva.hpp"
#include "util/constants.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

// One term of D's sum, for one beta (theta - theta0 or theta + theta0) and
// sign s:
//   cot((pi + s beta) / 2n) F(k L a_s(beta)),
// a_s(beta) = 2 cos^2((2 n pi N_s - beta) / 2), N_s the integer nearest
// (beta + s pi) / (2 n pi): its cotangent and cosine.  On a boundary a
// cotangent's pole meets F's zero; the product's limit is taken a hair off it,
// on the side where that boundary's geometrical-optics wave exists, because
// that is the side the path tests put a ray grazing it on (a path through an
// edge is open): the incident wave where |beta| < pi (`reflected` false),
// the reflection off face 0 where beta < pi and off face n where
// beta > (2n - 1) pi.  On the other side the ray would be kept and its
// edge's term would add to it, not stand in for it: 1 + 1/2 where the
// field is 1/2.
void term_angles(double n, double beta, int s, bool reflected, double& cot, double& c)
{
    double b = beta;
    if (std::abs(std::sin((pi + s * b) / (2.0 * n))) < 1e-12)
    {
        const bool down = reflected ? std::abs(b - pi) <= std::abs(b - (2.0 * n - 1.0) * pi) : b > 0.0;
        b += down ? -1e-9 : 1e-9;
    }
    const double nearest = std::round((b + s * pi) / (2.0 * n * pi));
    c = std::cos((2.0 * n * pi * nearest - b) / 2.0);
    cot = 1.0 / std::tan((pi + s * b) / (2.0 * n));
}

// Below X = 50, F(X) = -i sqrt(pi) e^(i pi/4) u w(u e^(i pi/4)), u = sqrt(X).
// The argument is real, so F is a function of one variable, and F / u is w
// along a line in the upper half-plane, where w is smooth and slowly varying:
// a polynomial in u on each of 32 intervals of [0, sqrt(50)], of degree 10,
// fitted to the Faddeeva form at the interval's Chebyshev nodes once, the
// first time it is wanted.  It agrees with the Faddeeva form to 7e-12
// relative, which is that form's own error along the line (against SciPy's
// wofz, 5e-12 at u below 1 and near 7); the fit's own is under 1e-14.  It
// costs 3 ns a call where the Faddeeva form costs 15, and that form was a
// third of a --buildings run.  Built on first use rather than at load,
// because the Faddeeva form reads a constant faddeeva.cpp initializes at
// load, in an order across translation units the language leaves open.
struct TransitionTable
{
    static constexpr int intervals = 32;
    static constexpr int degree = 10;
    double width = 0.0, inverse_width = 0.0;
    // Each interval's polynomial in t in [-1, 1] across it, lowest power first.
    double re[intervals][degree + 1] = {}, im[intervals][degree + 1] = {};

    TransitionTable()
    {
        width = std::sqrt(50.0) / intervals;
        inverse_width = 1.0 / width;
        // T_k's coefficients, to turn the Chebyshev series into one in t.
        double chebyshev[degree + 1][degree + 1] = {};
        chebyshev[0][0] = 1.0;
        chebyshev[1][1] = 1.0;
        for (int k = 2; k <= degree; k++)
            for (int p = 0; p <= k; p++)
                chebyshev[k][p] = (p > 0 ? 2.0 * chebyshev[k - 1][p - 1] : 0.0) - chebyshev[k - 2][p];
        const complex scale = complex(0.0, -1.0) * std::sqrt(pi) * std::polar(1.0, pi / 4.0);
        for (int i = 0; i < intervals; i++)
        {
            complex value[degree + 1];
            for (int j = 0; j <= degree; j++)
            {
                const double u = (i + 0.5 + 0.5 * std::cos(pi * (j + 0.5) / (degree + 1))) * width;
                value[j] = scale * faddeeva(std::polar(u, pi / 4.0));
            }
            for (int k = 0; k <= degree; k++)
            {
                complex c(0.0, 0.0);
                for (int j = 0; j <= degree; j++)
                    c += value[j] * std::cos(pi * k * (j + 0.5) / (degree + 1));
                c *= (k == 0 ? 1.0 : 2.0) / (degree + 1);
                for (int p = 0; p <= k; p++)
                {
                    re[i][p] += c.real() * chebyshev[k][p];
                    im[i][p] += c.imag() * chebyshev[k][p];
                }
            }
        }
    }

    complex operator()(double x) const
    {
        const double u = std::sqrt(x);
        const int i = std::min(static_cast<int>(u * inverse_width), intervals - 1);
        const double t = (u - (i + 0.5) * width) * (2.0 * inverse_width);
        const double* a = re[i];
        const double* b = im[i];
        double r = a[degree], m = b[degree];
        for (int p = degree - 1; p >= 0; p--)
        {
            r = std::fma(r, t, a[p]);
            m = std::fma(m, t, b[p]);
        }
        return complex(r * u, m * u);
    }
};

} // namespace

complex utd_transition(double x)
{
    // NaN goes through as NaN, as the Faddeeva form took it; it must not
    // index the table.
    if (std::isnan(x))
        return complex(x, x);
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
    static const TransitionTable table;
    return table(x);
}

WedgeAngles wedge_angles(double n, double theta, double theta0)
{
    // A hard wedge adds the incident (theta - theta0) and the reflected
    // (theta + theta0) pairs.
    WedgeAngles a;
    a.n = n;
    term_angles(n, theta - theta0, +1, false, a.cot[0], a.c[0]);
    term_angles(n, theta - theta0, -1, false, a.cot[1], a.c[1]);
    term_angles(n, theta + theta0, +1, true, a.cot[2], a.c[2]);
    term_angles(n, theta + theta0, -1, true, a.cot[3], a.c[3]);
    return a;
}

complex hard_wedge_coefficient(const WedgeAngles& a, double k, double r)
{
    const double kl = k * r;
    auto term = [&](int i) { return a.cot[i] * utd_transition(kl * 2.0 * a.c[i] * a.c[i]); };
    const complex sum = (term(0) + term(1)) + (term(2) + term(3));
    return -std::polar(1.0, pi / 4.0) / (2.0 * a.n * std::sqrt(2.0 * pi * k)) * sum;
}

WedgeSweep::WedgeSweep(double n, double theta) : n_(n), theta_(theta)
{
    for (int i = 0; i < 2; i++)
    {
        const double s = i == 0 ? 1.0 : -1.0;
        sin_x_[i] = std::sin((pi + s * theta) / (2.0 * n));
        cos_x_[i] = std::cos((pi + s * theta) / (2.0 * n));
    }
    sin_half_ = std::sin(theta / 2.0);
    cos_half_ = std::cos(theta / 2.0);
    for (int N = -3; N <= 4; N++)
    {
        sin_npn_[N + 3] = std::sin(n * pi * N);
        cos_npn_[N + 3] = std::cos(n * pi * N);
    }
}

WedgeAngles WedgeSweep::at(double theta0) const
{
    const double phi = theta0 / (2.0 * n_);
    const double sin_phi = std::sin(phi), cos_phi = std::cos(phi);
    const double sin_h0 = std::sin(theta0 / 2.0), cos_h0 = std::cos(theta0 / 2.0);
    WedgeAngles a;
    a.n = n_;
    // Terms 0 and 1: beta = theta - theta0, s = +1 and -1; 2 and 3: beta =
    // theta + theta0.  Each term's cotangent argument is (pi + s theta) / 2n
    // plus e phi, e = -s for the first beta and +s for the second.
    for (int t = 0; t < 4; t++)
    {
        const bool plus = t >= 2;
        const int i = t % 2;
        const double s = i == 0 ? 1.0 : -1.0;
        const double beta = plus ? theta_ + theta0 : theta_ - theta0;
        const double e = plus ? s : -s;
        const double sin_x = sin_x_[i] * cos_phi + e * cos_x_[i] * sin_phi;
        const double cos_x = cos_x_[i] * cos_phi - e * sin_x_[i] * sin_phi;
        const double nearest = std::round((beta + s * pi) / (2.0 * n_ * pi));
        const int index = static_cast<int>(nearest) + 3;
        if (std::abs(sin_x) < 1e-8 || index < 0 || index >= 8)
        {
            term_angles(n_, beta, i == 0 ? +1 : -1, plus, a.cot[t], a.c[t]);
            continue;
        }
        // c = cos(n pi N - beta / 2), beta / 2 = theta / 2 -+ theta0 / 2.
        const double sin_hb = plus ? sin_half_ * cos_h0 + cos_half_ * sin_h0 : sin_half_ * cos_h0 - cos_half_ * sin_h0;
        const double cos_hb = plus ? cos_half_ * cos_h0 - sin_half_ * sin_h0 : cos_half_ * cos_h0 + sin_half_ * sin_h0;
        a.c[t] = cos_npn_[index] * cos_hb + sin_npn_[index] * sin_hb;
        a.cot[t] = cos_x / sin_x;
    }
    return a;
}

complex hard_wedge_coefficient(double n, double theta, double theta0, double k, double r)
{
    return hard_wedge_coefficient(wedge_angles(n, theta, theta0), k, r);
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
