#include "edge_source.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "util/constants.hpp"
#include "wedge.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

// Gauss-Legendre's eight nodes and weights on [-1, 1].
constexpr std::array<double, 8> legendre_nodes {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
        -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
constexpr std::array<double, 8> legendre_weights {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
        0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};

// Gauss-Laguerre's nodes and weights, for the integral over [0, inf) of
// exp(-t) g(t): Newton's method on the Laguerre polynomial from the usual
// first guesses (Press et al., gaulag).
constexpr int laguerre_count = 20;
struct Laguerre
{
    std::array<double, laguerre_count> t {}, w {};
};

Laguerre laguerre_rule()
{
    Laguerre rule;
    const int n = laguerre_count;
    double z = 0.0;
    for (int i = 0; i < n; i++)
    {
        if (i == 0)
            z = 3.0 / (1.0 + 2.4 * n);
        else if (i == 1)
            z += 15.0 / (1.0 + 2.5 * n);
        else
        {
            const double a = i - 1;
            z += (1.0 + 2.55 * a) / (1.9 * a) * (z - rule.t[i - 2]);
        }
        double p1 = 1.0, p2 = 0.0, derivative = 1.0;
        for (int iteration = 0; iteration < 100; iteration++)
        {
            p1 = 1.0;
            p2 = 0.0;
            for (int j = 1; j <= n; j++)
            {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1 - z) * p2 - (j - 1) * p3) / j;
            }
            derivative = n * (p1 - p2) / z;
            const double step = p1 / derivative;
            z -= step;
            if (std::abs(step) <= 1e-15 * std::max(1.0, z))
                break;
        }
        rule.t[i] = z;
        rule.w[i] = -1.0 / (derivative * n * p2);
    }
    return rule;
}

const Laguerre laguerre = laguerre_rule();

// The integral's geometry for one source and receiver about an edge.
struct Geometry
{
    double r_s = 0.0, r_r = 0.0, z_s = 0.0, z_r = 0.0;
    double length0 = 0.0;  // R at the apex: the least path by the edge's line
    double rr = 0.0;       // r_s r_r

    double apex = 0.0;     // the apex's parameter z

    // eta at the edge's point z, from cosh(eta) - 1 in a form that cancels
    // neither near the apex nor far out along the edge:
    // (r_S b + r_R a)^2 / (r_S r_R (m l - a b + r_S r_R)), a = z - z_S,
    // b = z - z_R, where r_S b + r_R a = (r_S + r_R)(z - apex), and where a b
    // > 0, m l - a b = (r_S^2 r_R^2 + r_S^2 b^2 + r_R^2 a^2) / (m l + a b).
    double eta(double z) const
    {
        const double a = z - z_s, b = z - z_r;
        const double m = std::sqrt(r_s * r_s + a * a), l = std::sqrt(r_r * r_r + b * b);
        const double t = (r_s + r_r) * (z - apex);
        const double ab = a * b;
        const double excess = ab > 0.0 ? (rr * rr + r_s * r_s * b * b + r_r * r_r * a * a) / (m * l + ab) : m * l - ab;
        const double q = t * t / (rr * (excess + rr));
        return 2.0 * std::asinh(std::sqrt(0.5 * q));
    }
    // R(eta) = m + l = sqrt(R0^2 + 4 r_s r_r sinh^2(eta / 2)).
    double path(double eta) const
    {
        const double s = std::sinh(0.5 * eta);
        return std::sqrt(length0 * length0 + 4.0 * rr * s * s);
    }
    // The edge's point at eta, above the apex (`up`) or below it: where the
    // line meets the ellipsoid m + l = d = R(eta), the roots of a quadratic in
    // z - z_S whose discriminant is 16 d^2 (d^2 - R0^2)(d^2 - R0^2 + 4 r_S r_R),
    // which is (8 d r_S r_R sinh eta)^2: so they are
    // c (d^2 - r_R^2 + r_S^2 - c^2) / 2 (d^2 - c^2) +- d r_S r_R sinh(eta) / (d^2 - c^2),
    // c = z_R - z_S, with nothing to cancel near the apex.
    double point(double eta, bool up) const
    {
        const double d = path(eta), c = z_r - z_s;
        const double denominator = d * d - c * c;
        const double centre = c * (d * d - r_r * r_r + r_s * r_s - c * c) / (2.0 * denominator);
        const double half = d * rr * std::sinh(eta) / denominator;
        return z_s + centre + (up ? half : -half);
    }
    // dR / d eta.
    double rate(double eta) const
    {
        return rr * std::sinh(eta) / path(eta);
    }
    // The eta at which R has grown by `growth` past R(from).
    double grown(double from, double growth) const
    {
        const double d = path(from) + growth;
        const double s2 = (d - length0) * (d + length0) / (4.0 * rr);
        return 2.0 * std::asinh(std::sqrt(s2));
    }
};

Geometry geometry(const EdgeFrameEnds& ends)
{
    Geometry g;
    g.r_s = ends.r_s;
    g.r_r = ends.r_r;
    g.z_s = ends.z_s;
    g.z_r = ends.z_r;
    const double w = ends.r_s + ends.r_r, dz = ends.z_r - ends.z_s;
    g.length0 = std::sqrt(w * w + dz * dz);
    g.rr = ends.r_s * ends.r_r;
    g.apex = ends.z_s + (ends.z_r - ends.z_s) * ends.r_s / (ends.r_s + ends.r_r);
    return g;
}

// One of the directivity's four terms, sin(a) / (cosh(nu eta) - cos(a)) with
// a = nu phi_i, its denominator written 2 sinh^2(nu eta / 2) + 2 sin^2(a / 2)
// so that it does not cancel at the apex near a boundary, where cos(a) rounds
// to 1 within 1.5e-8 of it.
struct Term
{
    double s = 0.0;               // sin a
    double sh = 0.0, ch = 0.0;    // sin a/2, cos a/2
    double q = 0.0;               // sin^2 a/2

    // The integral of beta from eta_a to eta_b, both >= 0: the difference of
    // 2/nu atan(tanh(nu eta / 2) / tan(a / 2)) at the two, by atan's
    // subtraction formula, which holds the difference to within a quarter
    // turn without branch cuts.
    double integral(double nu, double eta_a, double eta_b) const
    {
        const double ta = std::tanh(0.5 * nu * eta_a), tb = std::tanh(0.5 * nu * eta_b);
        return 2.0 / nu * std::atan2((tb - ta) * sh * ch, sh * sh + ta * tb * ch * ch);
    }
};

struct Directivity
{
    double nu = 1.0;
    std::array<Term, 4> terms;
    // The narrowest peak's width in eta, the distance of a term's a from the
    // nearest pole 2 pi m, over nu.
    double width = 1.0;

    double operator()(double eta) const
    {
        const double h = std::sinh(0.5 * nu * eta), h2 = h * h;
        double sum = 0.0;
        for (const Term& t : terms)
            sum += t.s / (2.0 * (h2 + t.q));
        return sum;
    }
    complex operator()(const complex& eta) const
    {
        const complex h = std::sinh(0.5 * nu * eta), h2 = h * h;
        complex sum = 0.0;
        for (const Term& t : terms)
            sum += t.s / (2.0 * (h2 + t.q));
        return sum;
    }
    double integral(double eta_a, double eta_b) const
    {
        double sum = 0.0;
        for (const Term& t : terms)
            sum += t.integral(nu, eta_a, eta_b);
        return sum;
    }
};

Directivity directivity(double n, double theta_s, double theta_r)
{
    Directivity d;
    d.nu = 1.0 / n;
    const double phis[4] = {pi + theta_s + theta_r, pi + theta_s - theta_r, pi - theta_s + theta_r,
            pi - theta_s - theta_r};
    d.width = 1.0;
    for (int i = 0; i < 4; i++)
    {
        double a = d.nu * phis[i];
        // On a boundary, a at 2 pi m, the term is the limit from above, which
        // for every term is the side its geometrical-optics wave exists on
        // (wedge.cpp's term_angles has why that is the side wanted).  The
        // nearest 2 pi m by rounding: a stays within a few turns.
        double off = a - 2.0 * pi * std::floor(a / (2.0 * pi) + 0.5);
        if (std::abs(off) < 1e-12)
        {
            a += 1e-12 - off;
            off = 1e-12;
        }
        // sin a = 2 sin(a/2) cos(a/2): one sincos.
        const double sh = std::sin(0.5 * a), ch = std::cos(0.5 * a);
        d.terms[i] = {2.0 * sh * ch, sh, ch, sh * sh};
        d.width = std::min(d.width, std::abs(off) / d.nu);
    }
    return d;
}

// The directivity from its half-angles' parts, for a run of many points that
// share an edge: a_i / 2 = nu (pi +- theta_S +- theta_R) / 2 is P +- S +- R,
// so the four terms take their sines and cosines from those of P, S and R by
// the addition formulas, one sincos per angle that varies.  A term within
// 1e-12 of its pole takes directivity()'s own path.
struct HalfAngle
{
    double angle = 0.0;   // nu theta / 2
    double s = 0.0, c = 1.0;
};

HalfAngle half_angle(double n, double theta)
{
    HalfAngle h;
    h.angle = 0.5 * theta / n;
    h.s = std::sin(h.angle);
    h.c = std::cos(h.angle);
    return h;
}

Directivity directivity(double n, const HalfAngle& p, const HalfAngle& source, const HalfAngle& receiver,
        double theta_s, double theta_r)
{
    Directivity d;
    d.nu = 1.0 / n;
    // P + s S + r R for (s, r) = (+, +), (+, -), (-, +), (-, -).
    const double signs[4][2] = {{1.0, 1.0}, {1.0, -1.0}, {-1.0, 1.0}, {-1.0, -1.0}};
    for (int i = 0; i < 4; i++)
    {
        const double half = p.angle + signs[i][0] * source.angle + signs[i][1] * receiver.angle;
        const double off = 2.0 * half - 2.0 * pi * std::floor(half / pi + 0.5);
        if (std::abs(off) < 1e-12)
            return directivity(n, theta_s, theta_r);
        const double su = p.s * source.c + signs[i][0] * p.c * source.s;
        const double cu = p.c * source.c - signs[i][0] * p.s * source.s;
        const double sh = su * receiver.c + signs[i][1] * cu * receiver.s;
        const double ch = cu * receiver.c - signs[i][1] * su * receiver.s;
        d.terms[i] = {2.0 * sh * ch, sh, ch, sh * sh};
    }
    return d;
}

// Past the window the integrand is taken along its steepest-descent path,
// where it decays as exp(-k p): the integral from eta to infinity of
// beta exp(i k R) / R along it, by Gauss-Laguerre.  On the path
// R(h) = R(eta) + i p, so sinh^2(h/2) = sinh^2(eta/2) + (2 i R p - p^2) / (4 r_S r_R),
// which lies in the first quadrant: h stays right of the imaginary axis,
// where beta's poles are, and below Im h = pi.
complex descent(const Geometry& g, const Directivity& beta, double k, double eta)
{
    const double r = g.path(eta);
    const double s0 = std::sinh(0.5 * eta);
    complex sum = 0.0;
    for (int j = 0; j < laguerre_count; j++)
    {
        const double p = laguerre.t[j] / k;
        const complex rc(r, p);
        const complex s2 = s0 * s0 + complex(-p * p, 2.0 * r * p) / (4.0 * g.rr);
        const complex s = std::sqrt(s2);
        const complex h = 2.0 * std::asinh(s);
        // dh/dp = i R / (r_S r_R sinh h), sinh h = 2 sinh(h/2) cosh(h/2).
        const complex dh = complex(0.0, 1.0) * rc / (g.rr * 2.0 * s * std::sqrt(1.0 + s2));
        sum += laguerre.w[j] * beta(h) * dh / rc;
    }
    return std::polar(1.0 / k, k * r) * sum;
}

// The window: the integrand on the real line, out to where its phase has
// grown by this many radians past the branch's start (or the branch's end,
// if sooner).  Every singularity of the descent's integrand -- beta's poles
// at h = i (a - 2 pi m) / nu, the branch points at h = 0 and i pi -- has
// R(h) real and at most R0, so lies at least this far from the Laguerre
// rule's real t = k p axis, however near a boundary the terms are.
constexpr double window_phase = 12.0;

// The integral of beta exp(i k R) / R over [from, to] of eta.
complex branch(const Geometry& g, const Directivity& beta, double k, double from, double to)
{
    if (!(to > from))
        return 0.0;
    const double window = std::min(to, g.grown(from, window_phase / k));
    const double r0 = g.path(from);
    const complex f0 = std::polar(1.0 / r0, k * r0);
    // The peaks' closed part, against the integrand's value at the start (the
    // apex, when it lies on the stretch); the smooth rest by Gauss-Legendre,
    // panels graded out from the start at the narrowest peak's width and no
    // wider than three radians of phase.
    complex sum = f0 * beta.integral(from, window);
    double eta = from;
    double h = 0.25 * std::max(beta.width, 1e-4);
    while (eta < window)
    {
        h = std::min({h, 0.5, window - eta});
        const double limit = 3.0 / (k * g.rate(eta + h));
        if (h > limit)
            h = std::max(limit, 1e-12);
        complex panel = 0.0;
        for (int q = 0; q < 8; q++)
        {
            const double x = eta + 0.5 * h * (legendre_nodes[q] + 1.0);
            const double r = g.path(x);
            panel += legendre_weights[q] * beta(x) * (std::polar(1.0 / r, k * r) - f0);
        }
        sum += 0.5 * h * panel;
        eta += h;
        h *= 2.0;
    }
    if (window < to)
        sum += descent(g, beta, k, window) - descent(g, beta, k, to);
    return sum;
}

// Filon's rule of degree four: on a panel of path length [a, a + h], g at
// the five Chebyshev-Lobatto points x h fixes a quartic in x, whose integral
// against exp(i k R) is taken exactly however many wavelengths the panel
// spans.  For the same nodes at many wavenumbers: g does not depend on k.
constexpr int filon_count = 5;
struct Filon
{
    std::array<double, filon_count> x {};
    // The quartic's coefficients from the values: c_j = sum_i inverse[j][i] g_i.
    std::array<std::array<double, filon_count>, filon_count> inverse {};
};

Filon filon_rule()
{
    Filon rule;
    for (int i = 0; i < filon_count; i++)
        rule.x[i] = 0.5 * (1.0 - std::cos(pi * i / (filon_count - 1)));
    // Invert the Vandermonde matrix x_i^j by Gauss-Jordan elimination.
    double a[filon_count][2 * filon_count] = {};
    for (int i = 0; i < filon_count; i++)
    {
        double power = 1.0;
        for (int j = 0; j < filon_count; j++, power *= rule.x[i])
            a[i][j] = power;
        a[i][filon_count + i] = 1.0;
    }
    for (int c = 0; c < filon_count; c++)
    {
        int pivot = c;
        for (int r = c + 1; r < filon_count; r++)
            if (std::abs(a[r][c]) > std::abs(a[pivot][c]))
                pivot = r;
        for (int j = 0; j < 2 * filon_count; j++)
            std::swap(a[c][j], a[pivot][j]);
        const double d = a[c][c];
        for (int j = 0; j < 2 * filon_count; j++)
            a[c][j] /= d;
        for (int r = 0; r < filon_count; r++)
            if (r != c)
            {
                const double f = a[r][c];
                for (int j = 0; j < 2 * filon_count; j++)
                    a[r][j] -= f * a[c][j];
            }
    }
    // a's right half is V^-1, with V's rows the nodes: c = V^-1 g.
    for (int j = 0; j < filon_count; j++)
        for (int i = 0; i < filon_count; i++)
            rule.inverse[j][i] = a[j][filon_count + i];
    return rule;
}

const Filon filon = filon_rule();

// mu_j = the integral over [0, 1] of x^j exp(i phi x), j < filon_count: by
// the recurrence mu_j = (exp(i phi) - j mu_(j-1)) / (i phi) from a radian
// up, where it loses at most j!/phi^j (24) of its precision, and below by the
// power series sum (i phi)^n / (n! (n + j + 1)), whose even terms are the
// real part and odd ones the imaginary: polynomials in phi^2 with fixed
// coefficients, eighteen terms leaving under 1e-17.
constexpr int moment_terms = 18;
struct MomentSeries
{
    double c[filon_count][moment_terms];
};

constexpr MomentSeries moment_series()
{
    MomentSeries series {};
    double factorial = 1.0;
    for (int n = 0; n < moment_terms; n++)
    {
        if (n > 0)
            factorial *= n;
        for (int j = 0; j < filon_count; j++)
            series.c[j][n] = ((n / 2) % 2 == 0 ? 1.0 : -1.0) / (factorial * (n + j + 1));
    }
    return series;
}

constexpr MomentSeries moment_coefficients = moment_series();

std::array<complex, filon_count> filon_moments(double phi)
{
    std::array<complex, filon_count> mu;
    if (std::abs(phi) >= 1.0)
    {
        const complex e = std::polar(1.0, phi), by_i_phi(0.0, -1.0 / phi);
        mu[0] = (e - 1.0) * by_i_phi;
        for (int j = 1; j < filon_count; j++)
            mu[j] = (e - static_cast<double>(j) * mu[j - 1]) * by_i_phi;
        return mu;
    }
    const double u = phi * phi;
    for (int j = 0; j < filon_count; j++)
    {
        const double* c = moment_coefficients.c[j];
        double even = c[moment_terms - 2], odd = c[moment_terms - 1];
        for (int n = moment_terms - 4; n >= 0; n -= 2)
        {
            even = even * u + c[n];
            odd = odd * u + c[n + 1];
        }
        mu[j] = {even, phi * odd};
    }
    return mu;
}

// A panel of the transform: its start and width in path length and the
// quartic's coefficients in x = (R - a) / h.
struct FilonPanel
{
    double a = 0.0, h = 0.0;
    std::array<double, filon_count> c {};
};

FilonPanel filon_panel(double a, double h, const std::array<double, filon_count>& g)
{
    FilonPanel p;
    p.a = a;
    p.h = h;
    for (int j = 0; j < filon_count; j++)
    {
        double c = 0.0;
        for (int i = 0; i < filon_count; i++)
            c += filon.inverse[j][i] * g[i];
        p.c[j] = c;
    }
    return p;
}

complex filon_integral(const FilonPanel& p, double k)
{
    const std::array<complex, filon_count> mu = filon_moments(k * p.h);
    complex sum = 0.0;
    for (int j = 0; j < filon_count; j++)
        sum += p.c[j] * mu[j];
    return std::polar(p.h, k * p.a) * sum;
}

// A branch of the edge-source integral for many wavenumbers: the window
// near its start, Gauss-Legendre nodes resolving the largest wavenumber,
// whose peak part is closed against exp(i k R0) / R0 at each; past it Filon
// panels in R, geometric in R - R(start).
struct BranchSpectrum
{
    double r0 = 0.0;       // R at the branch's start
    double closed = 0.0;   // the integral of beta over the window
    std::vector<double> weight, length;  // the window's nodes: weight beta, R
    std::vector<FilonPanel> panels;

    complex at(double k) const
    {
        const complex f0 = std::polar(1.0 / r0, k * r0);
        complex sum = f0 * closed;
        for (std::size_t j = 0; j < weight.size(); j++)
            sum += weight[j] * (std::polar(1.0 / length[j], k * length[j]) - f0);
        for (const FilonPanel& p : panels)
            sum += filon_integral(p, k);
        return sum;
    }
};

// The window's phase at the largest wavenumber, and the Filon panels' ratio.
constexpr double spectrum_window_phase = 6.0;
constexpr double spectrum_ratio = 1.25;

BranchSpectrum branch_spectrum(const Geometry& g, const Directivity& beta, double k_max, double from, double to)
{
    BranchSpectrum b;
    if (!(to > from))
        return b;
    b.r0 = g.path(from);
    const double window = std::min(to, g.grown(from, spectrum_window_phase / k_max));
    b.closed = beta.integral(from, window);
    double eta = from;
    double h = 0.25 * std::max(beta.width, 1e-4);
    while (eta < window)
    {
        h = std::min({h, 0.5, window - eta});
        const double limit = 3.0 / (k_max * g.rate(eta + h));
        if (h > limit)
            h = std::max(limit, 1e-12);
        for (int q = 0; q < 8; q++)
        {
            const double x = eta + 0.5 * h * (legendre_nodes[q] + 1.0);
            const double r = g.path(x);
            b.weight.push_back(0.5 * h * legendre_weights[q] * beta(x));
            b.length.push_back(r);
        }
        eta += h;
        h *= 2.0;
    }
    if (!(window < to))
        return b;
    // In R the integrand is beta / (r_S r_R sinh eta): d eta / R = dR / (R dR/d eta).
    const double end = g.path(to), start = g.path(window);
    const double r_apex = g.length0;
    const auto integrand = [&](double r) {
        // eta from R: sinh^2(eta/2) = (R - R0)(R + R0) / 4 r_S r_R.
        const double d = r - r_apex;
        const double s2 = d * (2.0 * r_apex + d) / (4.0 * g.rr);
        const double e = 2.0 * std::asinh(std::sqrt(s2));
        return beta(e) / (g.rr * std::sinh(e));
    };
    // Each panel spans a fixed fraction of its distance from the apex, where
    // the integrand's singularity is.
    double a = start, width = std::max((start - r_apex) * (spectrum_ratio - 1.0), 1e-9 * end);
    while (a < end)
    {
        const double w = std::min(width, end - a);
        std::array<double, filon_count> values;
        for (int i = 0; i < filon_count; i++)
            values[i] = integrand(a + filon.x[i] * w);
        b.panels.push_back(filon_panel(a, w, values));
        a += w;
        width *= spectrum_ratio;
    }
    return b;
}

// Gauss-Kronrod's 15-point rule and its embedded 7-point Gauss rule
// (Piessens et al., QUADPACK): nodes on [0, 1] of the half-interval, the
// last the middle; the Gauss rule takes the odd ones and the middle.
constexpr double kronrod_nodes[8] = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
        0.864864423359769072789712788640926, 0.741531185599394439863864773280788, 0.586087235467691130294144845693013,
        0.405845151377397166906606412076961, 0.207784955007898467600689403773245, 0.0};
constexpr double kronrod_weights[8] = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
        0.104790010322250183839876322541518, 0.140653259715525918745189590510238, 0.169004726639267902826583426598550,
        0.190350578064785409913256402421014, 0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
constexpr double gauss7_weights[4] = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
        0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

struct Panel
{
    double a = 0.0, b = 0.0;
    complex value;
    double error = 0.0;
};

template <class F>
Panel kronrod(const F& f, double a, double b)
{
    const double c = 0.5 * (a + b), h = 0.5 * (b - a);
    const complex middle = f(c);
    complex k = kronrod_weights[7] * middle, g = gauss7_weights[3] * middle;
    for (int j = 0; j < 7; j++)
    {
        const double x = h * kronrod_nodes[j];
        const complex pair = f(c - x) + f(c + x);
        k += kronrod_weights[j] * pair;
        if (j % 2 == 1)
            g += gauss7_weights[(j - 1) / 2] * pair;
    }
    return {a, b, h * k, std::abs(h * (k - g))};
}

// The integral of f over the cuts' span, globally adaptive: the panel with
// the largest error estimate is halved until the estimates sum to under
// `tolerance` of the total, or `floor`, or `limit` panels.
template <class F>
complex adaptive(const F& f, const std::vector<double>& cuts, double tolerance, double floor, std::size_t limit)
{
    std::vector<Panel> panels;
    for (std::size_t i = 0; i + 1 < cuts.size(); i++)
        if (cuts[i + 1] > cuts[i])
            panels.push_back(kronrod(f, cuts[i], cuts[i + 1]));
    while (true)
    {
        complex total = 0.0;
        double error = 0.0;
        std::size_t worst = 0;
        for (std::size_t i = 0; i < panels.size(); i++)
        {
            total += panels[i].value;
            error += panels[i].error;
            if (panels[i].error > panels[worst].error)
                worst = i;
        }
        if (panels.empty() || error <= std::max(tolerance * std::abs(total), floor) || panels.size() >= limit)
            return total;
        const Panel p = panels[worst];
        const double middle = 0.5 * (p.a + p.b);
        if (!(middle > p.a && middle < p.b))
            return total;
        panels[worst] = kronrod(f, p.a, middle);
        panels.push_back(kronrod(f, middle, p.b));
    }
}

// A point about an edge: its distance from the line, its angle (put on the
// nearer face if round-off leaves it just inside the wedge), and its
// parameter along the line.
struct About
{
    double r = 0.0, theta = 0.0, z = 0.0;
};

// An angle about the edge in its air: one past n pi, which only round-off
// puts a point in a face's plane at, goes onto the nearer face, 0 or n pi.
// Without it a point in face 0's plane came out at 0 or just under 2 pi as
// the rounding went, and a rotated district's partners changed.
double onto_air(const Edge& edge, double theta)
{
    const double open = edge.n * pi;
    if (theta > open)
        theta = theta - open < 2.0 * pi - theta ? open : 0.0;
    return theta;
}

About about(const Edge& edge, const Vec3& p)
{
    const Edge::Projection q = edge.project(p);
    return {q.normal.norm(), onto_air(edge, edge.angle(q.normal)), q.along};
}

} // namespace

EdgeFrameEnds edge_frame_ends(const Edge& edge, const Vec3& source, const Vec3& receiver)
{
    const Edge::Projection s = edge.project(source), r = edge.project(receiver);
    EdgeFrameEnds e;
    e.r_s = s.normal.norm();
    e.theta_s = onto_air(edge, edge.angle(s.normal));
    e.z_s = s.along;
    e.r_r = r.normal.norm();
    e.theta_r = onto_air(edge, edge.angle(r.normal));
    e.z_r = r.along;
    return e;
}

complex edge_source_field(double n, const EdgeFrameEnds& ends, double k, double z_lo, double z_hi)
{
    if (!(z_hi > z_lo) || !(ends.r_s > 0.0) || !(ends.r_r > 0.0))
        return 0.0;
    const Geometry g = geometry(ends);
    const Directivity beta = directivity(n, ends.theta_s, ends.theta_r);

    // Out from the apex each side where it lies on the stretch, otherwise
    // from the nearer end to the farther.  The two sides share eta's
    // integrand.
    const double apex = g.apex;
    complex sum;
    if (apex <= z_lo)
        sum = branch(g, beta, k, g.eta(z_lo), g.eta(z_hi));
    else if (apex >= z_hi)
        sum = branch(g, beta, k, g.eta(z_hi), g.eta(z_lo));
    else
        sum = branch(g, beta, k, 0.0, g.eta(z_hi)) + branch(g, beta, k, 0.0, g.eta(z_lo));
    return -beta.nu / (4.0 * pi) * sum;
}

std::vector<complex> edge_source_spectrum(double n, const EdgeFrameEnds& ends, double z_lo, double z_hi,
        const std::vector<double>& wavenumbers)
{
    std::vector<complex> out(wavenumbers.size(), 0.0);
    if (wavenumbers.empty() || !(z_hi - z_lo > 1e-12 * std::max({1.0, std::abs(z_lo), std::abs(z_hi)}))
            || !(ends.r_s > 0.0) || !(ends.r_r > 0.0))
        return out;
    const double k_max = *std::max_element(wavenumbers.begin(), wavenumbers.end());
    const Geometry g = geometry(ends);
    const Directivity beta = directivity(n, ends.theta_s, ends.theta_r);
    std::vector<BranchSpectrum> branches;
    if (g.apex <= z_lo)
        branches.push_back(branch_spectrum(g, beta, k_max, g.eta(z_lo), g.eta(z_hi)));
    else if (g.apex >= z_hi)
        branches.push_back(branch_spectrum(g, beta, k_max, g.eta(z_hi), g.eta(z_lo)));
    else
    {
        branches.push_back(branch_spectrum(g, beta, k_max, 0.0, g.eta(z_hi)));
        branches.push_back(branch_spectrum(g, beta, k_max, 0.0, g.eta(z_lo)));
    }
    for (std::size_t i = 0; i < wavenumbers.size(); i++)
    {
        complex sum = 0.0;
        for (const BranchSpectrum& b : branches)
            sum += b.at(wavenumbers[i]);
        out[i] = -beta.nu / (4.0 * pi) * sum;
    }
    return out;
}

complex edge_source_field(const Edge& edge, const Vec3& source, const Vec3& receiver, double k)
{
    return edge_source_field(edge.n, edge_frame_ends(edge, source, receiver), k, 0.0, edge.length);
}

double window_taper(double path, const PathWindow& window)
{
    if (!(window.reach > 0.0))
        return 1.0;
    const double x = (path - window.reference) / window.reach * 2.0 - 1.0;
    if (x <= 0.0)
        return 1.0;
    if (x >= 1.0)
        return 0.0;
    return 1.0 / (1.0 + std::exp(1.0 / (1.0 - x) - 1.0 / x));
}

complex edge_source_field(double n, const EdgeFrameEnds& ends, double k, double z_lo, double z_hi,
        const PathWindow& window, double tolerance)
{
    if (!(z_hi > z_lo) || !(ends.r_s > 0.0) || !(ends.r_r > 0.0))
        return 0.0;
    const Geometry g = geometry(ends);
    const Directivity beta = directivity(n, ends.theta_s, ends.theta_r);
    const double most = window.reference + window.reach;
    if (!(most > g.length0))
        return 0.0;
    const double limit = g.grown(0.0, most - g.length0);
    const auto branch_integral = [&](double from, double to) -> complex {
        to = std::min(to, limit);
        if (!(to > from))
            return 0.0;
        const double r0 = g.path(from);
        const complex start = std::polar(window_taper(r0, window) / r0, k * r0);
        std::vector<double> cuts {from};
        double h = 0.25 * std::max(beta.width, 1e-4);
        while (cuts.back() < to)
        {
            const double t = cuts.back();
            h = std::min({h, 0.5, 2.0 / (k * g.rate(t + h))});
            cuts.push_back(std::min(to, t + h));
            h *= 2.0;
        }
        const auto rest = [&](double eta) {
            const double d = g.path(eta);
            return beta(eta) * (std::polar(window_taper(d, window) / d, k * d) - start);
        };
        const double closed = beta.integral(from, to);
        const double scale = std::abs(start) * (std::abs(closed) + 1.0);
        return start * closed + adaptive(rest, cuts, tolerance, 1e-3 * tolerance * scale, 4000);
    };
    complex sum;
    if (g.apex <= z_lo)
        sum = branch_integral(g.eta(z_lo), g.eta(z_hi));
    else if (g.apex >= z_hi)
        sum = branch_integral(g.eta(z_hi), g.eta(z_lo));
    else
        sum = branch_integral(0.0, g.eta(z_hi)) + branch_integral(0.0, g.eta(z_lo));
    return -beta.nu / (4.0 * pi) * sum;
}

double least_twice_path(const Edge& first, const Edge& second, const Vec3& source, const Vec3& receiver,
        const Stretch& along_first)
{
    const About r = about(second, receiver);
    const auto least = [&](double z1) {
        const Vec3 p1 = first.start + z1 * first.direction;
        const About q = about(second, p1);
        return (source - p1).norm() + std::hypot(q.r + r.r, r.z - q.z);
    };
    double a = along_first.lo, b = along_first.hi;
    for (int i = 0; i < 200 && b - a > 1e-12 * (1.0 + std::abs(a)); i++)
    {
        const double c = a + 0.381966 * (b - a), d = b - 0.381966 * (b - a);
        if (least(c) < least(d))
            b = d;
        else
            a = c;
    }
    return least(0.5 * (a + b));
}

complex twice_diffracted_field(const Edge& first, const Edge& second, const Vec3& source, const Vec3& receiver,
        double k, const Stretch& along_first, const Stretch& along_second, double plane_factor, double tolerance,
        double window_phase)
{
    const About s = about(first, source), r = about(second, receiver);
    if (!(s.r > 0.0) || !(r.r > 0.0) || !(along_first.hi > along_first.lo) || !(along_second.hi > along_second.lo))
        return 0.0;
    const double inner_tolerance = 0.1 * tolerance;

    // m + the least n + l along the second edge's line, at the first edge's
    // point z1: convex, as the least over z2 of a jointly convex path.
    const auto least = [&](double z1) {
        const Vec3 p1 = first.start + z1 * first.direction;
        const About q = about(second, p1);
        return (source - p1).norm() + std::hypot(q.r + r.r, r.z - q.z);
    };
    // The window: where the path is within window_phase / k of its least
    // over the stretches, its outer half tapered (Planck's taper, smooth to
    // every order, so that what it leaves out is an integral of an
    // oscillating function with no stationary point).
    double z1_lo = along_first.lo, z1_hi = along_first.hi, shortest = 0.0, reach = 0.0;
    if (window_phase > 0.0)
    {
        double a = along_first.lo, b = along_first.hi;
        for (int i = 0; i < 200 && b - a > 1e-12 * (1.0 + std::abs(a)); i++)
        {
            const double c = a + 0.381966 * (b - a), d = b - 0.381966 * (b - a);
            if (least(c) < least(d))
                b = d;
            else
                a = c;
        }
        const double middle = 0.5 * (a + b);
        shortest = least(middle);
        reach = window_phase / k;
        const auto edge_of = [&](double inside, double outside) {
            if (least(outside) <= shortest + reach)
                return outside;
            for (int i = 0; i < 200; i++)
            {
                const double c = 0.5 * (inside + outside);
                (least(c) <= shortest + reach ? inside : outside) = c;
            }
            return 0.5 * (inside + outside);
        };
        z1_lo = edge_of(middle, along_first.lo);
        z1_hi = edge_of(middle, along_first.hi);
    }
    const auto taper = [&](double path) {
        if (!(window_phase > 0.0))
            return 1.0;
        const double x = (path - shortest) / reach * 2.0 - 1.0;
        if (x <= 0.0)
            return 1.0;
        if (x >= 1.0)
            return 0.0;
        return 1.0 / (1.0 + std::exp(1.0 / (1.0 - x) - 1.0 / x));
    };

    // At the first edge's point z1: exp(i k m) / m times the integral along
    // the second edge of beta1 beta2 exp(i k (n + l)) / (n l), in the second
    // edge's eta for the point and the receiver.
    const auto inner = [&](double z1) -> complex {
        const Vec3 p1 = first.start + z1 * first.direction;
        const double m = (source - p1).norm();
        const About q = about(second, p1);
        if (!(q.r > 0.0))
            return 0.0;
        const EdgeFrameEnds ends {q.r, q.theta, q.z, r.r, r.theta, r.z};
        const Geometry g = geometry(ends);
        const Directivity beta2 = directivity(second.n, q.theta, r.theta);
        // The first edge's directivity toward the second's point z2.
        const auto beta1 = [&](double z2) {
            const About t = about(first, second.start + z2 * second.direction);
            if (!(t.r > 0.0))
                return 0.0;
            const Geometry g1 = geometry({s.r, s.theta, s.z, t.r, t.theta, t.z});
            return directivity(first.n, s.theta, t.theta)(g1.eta(z1));
        };
        // Within the window, eta no further than where m + n + l reaches its
        // edge.
        double limit = std::numeric_limits<double>::infinity();
        if (window_phase > 0.0)
        {
            const double most = shortest + reach - m;
            if (!(most > g.length0))
                return 0.0;
            limit = g.grown(0.0, most - g.length0);
        }
        const auto branch_integral = [&](double from, double to, bool up) -> complex {
            to = std::min(to, limit);
            if (!(to > from))
                return 0.0;
            const double r0 = g.path(from);
            const complex start = beta1(g.point(from, up)) * std::polar(taper(m + r0) / r0, k * r0);
            std::vector<double> cuts {from};
            double h = 0.25 * std::max(beta2.width, 1e-4);
            while (cuts.back() < to)
            {
                const double t = cuts.back();
                h = std::min({h, 0.5, 2.0 / (k * g.rate(t + h))});
                cuts.push_back(std::min(to, t + h));
                h *= 2.0;
            }
            const auto rest = [&](double eta) {
                const double d = g.path(eta);
                return beta2(eta) * (beta1(g.point(eta, up)) * std::polar(taper(m + d) / d, k * d) - start);
            };
            const double closed = beta2.integral(from, to);
            const double scale = std::abs(start) * (std::abs(closed) + 1.0);
            return start * closed + adaptive(rest, cuts, inner_tolerance, 1e-3 * inner_tolerance * scale, 4000);
        };
        complex sum;
        if (g.apex <= along_second.lo)
            sum = branch_integral(g.eta(along_second.lo), g.eta(along_second.hi), true);
        else if (g.apex >= along_second.hi)
            sum = branch_integral(g.eta(along_second.hi), g.eta(along_second.lo), false);
        else
            sum = branch_integral(0.0, g.eta(along_second.hi), true)
                    + branch_integral(0.0, g.eta(along_second.lo), false);
        return std::polar(1.0 / m, k * m) * sum;
    };

    // The outer cuts: no wider than two radians of the least path's phase.
    std::vector<double> cuts {z1_lo};
    const double span = z1_hi - z1_lo;
    if (!(span > 0.0))
        return 0.0;
    while (cuts.back() < z1_hi)
    {
        const double z = cuts.back(), dz = 1e-6 * std::max(span, 1.0);
        const double rate = std::abs(least(z + dz) - least(z)) / dz;
        const double h = std::min(0.125 * span, 2.0 / (k * std::max(rate, 1e-3)));
        cuts.push_back(std::min(z1_hi, z + std::max(h, 1e-6 * span)));
    }
    const complex sum = adaptive(inner, cuts, tolerance, 0.0, 4000);
    return plane_factor / (first.n * second.n * 16.0 * pi * pi) * sum;
}

namespace {

// The mirror of p in the plane through the edge with unit normal `normal`.
Vec3 mirror(const Edge& e, const Vec3& normal, const Vec3& p)
{
    return p - 2.0 * (p - e.start).dot(normal) * normal;
}

// The second face's unit normal toward the air.
Vec3 face_n_normal(const Edge& e)
{
    return std::sin(e.n * pi) * e.face0 - std::cos(e.n * pi) * e.air;
}

double angle_about(const Edge& e, const Vec3& p)
{
    const Vec3 d = p - e.start;
    return onto_air(e, e.angle(d - d.dot(e.direction) * e.direction));
}

// Where along `other` a condition on the angle about `edge` of other's point
// holds: an interval, the angle being monotonic along a line, sampled and
// its ends found by bisection; lo > hi where it holds nowhere.
template <class C>
Stretch where_along(const Edge& edge, const Edge& other, const C& condition)
{
    const int samples = 400;
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    // A point on the edge's own line -- where two edges meet at a corner --
    // has no angle about it.
    const double hair = 1e-9 * std::max(edge.length, other.length);
    const auto off_line = [&](double z) { return edge.project(other.start + z * other.direction).normal.norm() > hair; };
    for (int i = 0; i <= samples; i++)
    {
        const double z = other.length * i / samples;
        if (off_line(z) && condition(angle_about(edge, other.start + z * other.direction)))
        {
            lo = std::min(lo, z);
            hi = std::max(hi, z);
        }
    }
    if (lo > hi)
        return {1.0, 0.0};
    const auto refine = [&](double inside, double outside) {
        if (outside < 0.0 || outside > other.length)
            return inside;
        for (int i = 0; i < 60; i++)
        {
            const double c = 0.5 * (inside + outside);
            (condition(angle_about(edge, other.start + c * other.direction)) ? inside : outside) = c;
        }
        return 0.5 * (inside + outside);
    };
    const double step = other.length / samples;
    return {refine(lo, lo - step), refine(hi, hi + step)};
}

} // namespace

complex keller_once_field(const Edge& edge, const Vec3& source, const Vec3& receiver, double k, double lo, double hi,
        const PathWindow& window)
{
    const EdgeFrameEnds f = edge_frame_ends(edge, source, receiver);
    if (!(f.r_s > 0.0) || !(f.r_r > 0.0))
        return 0.0;
    const double apex = f.z_s + (f.z_r - f.z_s) * f.r_s / (f.r_s + f.r_r);
    if (apex < lo || apex > hi)
        return 0.0;
    const double s_in = std::hypot(f.r_s, apex - f.z_s), s_out = std::hypot(f.r_r, f.z_r - apex);
    const double sin_beta = f.r_s / s_in;
    return std::polar(window_taper(s_in + s_out, window) / std::sqrt(s_in * s_out * (s_in + s_out)), k * (s_in + s_out))
            * hard_wedge_coefficient(edge.n, f.theta_r, f.theta_s, k * sin_beta, s_in * s_out * sin_beta / (s_in + s_out))
            / std::sqrt(sin_beta);
}

std::vector<PartnerWave> twice_diffracted_partners(const Edge& first, const Edge& second, const Vec3& source,
        const Vec3& receiver)
{
    std::vector<PartnerWave> partners;
    // A part reaching within a hair of an end reaches the end: the search
    // skips a vertex the two edges share, and its bisection toward it stops a
    // hair short, which made a wave over the whole edge look partial.
    const auto add = [&](bool over_second, const Vec3& image, Stretch part, const Edge& along) {
        const double hair = 1e-9 * along.length;
        if (part.lo <= hair)
            part.lo = 0.0;
        if (part.hi >= along.length - hair)
            part.hi = along.length;
        if (part.hi - part.lo > hair && (part.lo > 0.0 || part.hi < along.length))
            partners.push_back({over_second, image, part});
    };
    // The first edge's waves at the second edge's points: the incident wave,
    // and the reflections off its faces, each from the source or its image.
    const double ts = angle_about(first, source), n1 = first.n;
    add(true, source, where_along(first, second, [&](double t) { return std::abs(t - ts) < pi; }), second);
    add(true, mirror(first, first.air, source), where_along(first, second, [&](double t) { return t + ts < pi; }),
            second);
    add(true, mirror(first, face_n_normal(first), source),
            where_along(first, second, [&](double t) { return t + ts > (2.0 * n1 - 1.0) * pi; }), second);
    // The second edge's waves toward the receiver, from the first edge's
    // points: their receivers the receiver or its images in its faces.
    const double tr = angle_about(second, receiver), n2 = second.n;
    add(false, receiver, where_along(second, first, [&](double t) { return std::abs(t - tr) < pi; }), first);
    add(false, mirror(second, second.air, receiver),
            where_along(second, first, [&](double t) { return t + tr < pi; }), first);
    add(false, mirror(second, face_n_normal(second), receiver),
            where_along(second, first, [&](double t) { return t + tr > (2.0 * n2 - 1.0) * pi; }), first);
    return partners;
}

complex localized_twice_diffracted_field(const Edge& first, const Edge& second, const Vec3& source,
        const Vec3& receiver, double k, double plane_factor, double window_phase, double tolerance)
{
    const Stretch whole1 {0.0, first.length}, whole2 {0.0, second.length};
    const PathWindow window {least_twice_path(first, second, source, receiver, whole1), window_phase / k};
    complex sum = twice_diffracted_field(first, second, source, receiver, k, whole1, whole2, plane_factor, tolerance,
            window_phase);
    for (const PartnerWave& p : twice_diffracted_partners(first, second, source, receiver))
    {
        const Edge& edge = p.over_second ? second : first;
        const Vec3& s = p.over_second ? p.image : source;
        const Vec3& r = p.over_second ? receiver : p.image;
        sum += edge_source_field(edge.n, edge_frame_ends(edge, s, r), k, p.part.lo, p.part.hi, window, 0.1 * tolerance)
                - keller_once_field(edge, s, r, k, p.part.lo, p.part.hi, window);
    }
    return sum;
}

namespace {

// Golden-section search for the least of a convex function on [a, b].
template <class F>
double convex_least(const F& f, double a, double b)
{
    for (int i = 0; i < 200 && b - a > 1e-13 * (1.0 + std::abs(a) + std::abs(b)); i++)
    {
        const double c = a + 0.381966011250105 * (b - a), d = b - 0.381966011250105 * (b - a);
        if (f(c) < f(d))
            b = d;
        else
            a = c;
    }
    return 0.5 * (a + b);
}

// Where a continuous function crosses `level` on [inside, outside], given
// f(inside) <= level < f(outside): by bisection.
// By the Illinois method: regula falsi that halves the stale end's weight.
template <class F>
double crossing(const F& f, double level, double inside, double outside)
{
    double fa = f(inside) - level, fb = f(outside) - level;
    double a = inside, b = outside;
    int side = 0;
    for (int i = 0; i < 100; i++)
    {
        const double c = fb - fa != 0.0 ? b - fb * (b - a) / (fb - fa) : 0.5 * (a + b);
        if (!(c > std::min(a, b) && c < std::max(a, b)) || std::abs(b - a) < 1e-14 * (1.0 + std::abs(a)))
            break;
        const double fc = f(c) - level;
        if (fc == 0.0)
            return c;
        if ((fc > 0.0) == (fb > 0.0))
        {
            b = c;
            fb = fc;
            if (side == -1)
                fa *= 0.5;
            side = -1;
        }
        else
        {
            a = c;
            fa = fc;
            if (side == 1)
                fb *= 0.5;
            side = 1;
        }
    }
    return std::abs(fa) < std::abs(fb) ? a : b;
}

// The second-order integral's level curves, G(L) = the integral of
// beta1 beta2 / (m n l) over the curve m + n + l = L in the stretches'
// rectangle, divided by |grad (m + n + l)|; the field is
// (f nu1 nu2 / 16 pi^2) times the integral of G(L) exp(i k L) dL, so G does
// not depend on k.  The curve is taken as z1 runs over the interval where it
// exists, and for each z1 the two points of the second edge -- above and
// below its apex for the first edge's point and the receiver -- where
// n + l = L - m, which the second edge's eta gives in closed form.
class LevelCurves
{
public:
    LevelCurves(const Edge& first, const Edge& second, const Vec3& source, const Vec3& receiver, const Stretch& one,
            const Stretch& two)
        : first_(first), second_(second), source_(source), receiver_(receiver), one_(one), two_(two),
          s_(about(first, source)), r_(about(second, receiver))
    {
        // For parallel edges each directivity's angles are the same along the
        // other edge: build them once.
        p1_ = half_angle(first.n, pi);
        p2_ = half_angle(second.n, pi);
        source_half_ = half_angle(first.n, s_.theta);
        receiver_half_ = half_angle(second.n, r_.theta);
        parallel_ = std::abs(first.direction.dot(second.direction)) > 1.0 - 1e-12;
        if (parallel_)
        {
            const About t = about(first_, second_.start);
            const About q = about(second_, first_.start);
            beta1_ = directivity(first_.n, s_.theta, t.theta);
            beta2_ = directivity(second_.n, q.theta, r_.theta);
            narrow1_ = beta1_.width;
            narrow2_ = beta2_.width;
        }
        else
            // The narrowest each directivity's peak gets along the other edge.
            for (int i = 0; i <= 32; i++)
            {
                const double z2 = two_.lo + (two_.hi - two_.lo) * i / 32.0;
                const double z1 = one_.lo + (one_.hi - one_.lo) * i / 32.0;
                const double t = about(first_, second_.start + z2 * second_.direction).theta;
                const double q = about(second_, first_.start + z1 * first_.direction).theta;
                narrow1_ = std::min(narrow1_, directivity(first_.n, s_.theta, t).width);
                narrow2_ = std::min(narrow2_, directivity(second_.n, q, r_.theta).width);
            }
        const auto line = [&](double z1) { return psi(z1); };
        psi_least_ = convex_least(line, one_.lo, one_.hi);
        for (int side = 0; side < 2; side++)
        {
            const double z2 = side == 0 ? two_.lo : two_.hi;
            side_least_[side] = convex_least([&](double z1) { return path(z1, z2); }, one_.lo, one_.hi);
        }
    }

    // The narrowest either directivity's peak gets in eta, and the least
    // distance between the edges' points.
    double least() const
    {
        return clamped(convex_least([&](double z1) { return clamped(z1); }, one_.lo, one_.hi));
    }
    double narrowest() const
    {
        return std::min(narrow1_, narrow2_);
    }
    double gap() const
    {
        double least = std::numeric_limits<double>::infinity();
        for (int i = 0; i <= 8; i++)
        {
            const double z1 = one_.lo + (one_.hi - one_.lo) * i / 8.0;
            least = std::min(least, about(second_, first_.start + z1 * first_.direction).r);
        }
        return least;
    }

    // m + n + l.
    double path(double z1, double z2) const
    {
        const Vec3 p1 = first_.start + z1 * first_.direction, p2 = second_.start + z2 * second_.direction;
        return (source_ - p1).norm() + (p2 - p1).norm() + (receiver_ - p2).norm();
    }
    // m + the least n + l along the second edge's line, and along its stretch.
    double psi(double z1) const
    {
        const Vec3 p1 = first_.start + z1 * first_.direction;
        const About q = about(second_, p1);
        const double w = q.r + r_.r, dz = r_.z - q.z;
        return (source_ - p1).norm() + std::sqrt(w * w + dz * dz);
    }
    double clamped(double z1) const
    {
        const Vec3 p1 = first_.start + z1 * first_.direction;
        const About q = about(second_, p1);
        const double apex = q.r + r_.r > 0.0 ? q.z + (r_.z - q.z) * q.r / (q.r + r_.r) : q.z;
        const double z2 = std::clamp(apex, two_.lo, two_.hi);
        const Vec3 p2 = second_.start + z2 * second_.direction;
        return (source_ - p1).norm() + (p2 - p1).norm() + (receiver_ - p2).norm();
    }

    // The least path over the rectangle, and the values where G is not
    // smooth: the least along each side, where the curve first meets it, and
    // the corners, where it leaves one.  Sorted, from the least.
    std::vector<double> breaks(double resolution) const
    {
        std::vector<double> b;
        const double least_z1 = convex_least([&](double z1) { return clamped(z1); }, one_.lo, one_.hi);
        const double least = clamped(least_z1);
        b.push_back(least);
        for (int side = 0; side < 2; side++)
            b.push_back(path(side_least_[side], side == 0 ? two_.lo : two_.hi));
        for (double z1 : {one_.lo, one_.hi})
        {
            const double z2 = convex_least([&](double z) { return path(z1, z); }, two_.lo, two_.hi);
            b.push_back(path(z1, z2));
        }
        for (double z1 : {one_.lo, one_.hi})
            for (double z2 : {two_.lo, two_.hi})
                b.push_back(path(z1, z2));
        // Where a boundary of the first edge's (seen from the source) crosses
        // the second edge, or of the second's (seen from the receiver) the
        // first: a directivity term's pole meets the real line there, on its
        // edge's ridge, and G steps.
        for (double z2 : boundary_points(true))
        {
            const About t = about(first_, second_.start + z2 * second_.direction);
            const double z1 = s_.z + (t.z - s_.z) * s_.r / (s_.r + t.r);
            if (z1 >= one_.lo && z1 <= one_.hi)
                b.push_back(path(z1, z2));
        }
        for (double z1 : boundary_points(false))
        {
            const About q = about(second_, first_.start + z1 * first_.direction);
            const double z2 = q.z + (r_.z - q.z) * q.r / (q.r + r_.r);
            if (z2 >= two_.lo && z2 <= two_.hi)
                b.push_back(path(z1, z2));
        }
        std::sort(b.begin(), b.end());
        // Breaks closer than `resolution` are one: the segment from each holds
        // G still that close to its start anyway.
        const double close = std::min(1e-6 * (b.back() - least), resolution);
        std::vector<double> out;
        for (double x : b)
            if (x >= least && (out.empty() || x > out.back() + close))
                out.push_back(x);
        if (out.size() > 1 && b.back() > out.back())
            out.back() = b.back();
        return out;
    }

    double operator()(double level) const
    {
        // The first edge's interval where the second edge's line has points
        // at the level: psi is convex.
        if (!(psi(psi_least_) < level))
            return 0.0;
        double a = one_.lo, b = one_.hi;
        bool singular_a = false, singular_b = false;
        const auto line = [&](double z1) { return psi(z1); };
        if (psi(a) > level)
        {
            a = crossing(line, level, psi_least_, a);
            singular_a = true;
        }
        if (psi(b) > level)
        {
            b = crossing(line, level, psi_least_, b);
            singular_b = true;
        }
        if (!(b > a))
            return 0.0;
        // Cuts: where either branch meets a side of the rectangle, and where
        // it crosses the first edge's ridge (its apex for the source and the
        // second edge's point, beta1's peak).
        std::vector<double> cuts {a, b}, ridges;
        for (int side = 0; side < 2; side++)
        {
            const double z2 = side == 0 ? two_.lo : two_.hi;
            const auto along = [&](double z1) { return path(z1, z2); };
            const double m = std::clamp(side_least_[side], a, b);
            if (along(m) < level)
            {
                if (along(a) > level)
                    cuts.push_back(crossing(along, level, m, a));
                if (along(b) > level)
                    cuts.push_back(crossing(along, level, m, b));
            }
        }
        for (int branch = 0; branch < 2 && narrow1_ < 0.3; branch++)
        {
            const bool up = branch == 0;
            const int samples = 16;
            double previous = std::numeric_limits<double>::quiet_NaN(), previous_z = a;
            for (int i = 0; i <= samples; i++)
            {
                const double z1 = a + (b - a) * (0.5 - 0.5 * std::cos(pi * i / samples));
                const double f = ridge(z1, level, up);
                if (std::isfinite(f) && std::isfinite(previous) && (f > 0.0) != (previous > 0.0))
                {
                    double lo = previous_z, hi = z1, flo = previous;
                    for (int j = 0; j < 100; j++)
                    {
                        const double c = 0.5 * (lo + hi);
                        const double fc = ridge(c, level, up);
                        if (!std::isfinite(fc))
                            break;
                        if ((fc > 0.0) == (flo > 0.0))
                        {
                            lo = c;
                            flo = fc;
                        }
                        else
                            hi = c;
                    }
                    cuts.push_back(0.5 * (lo + hi));
                    ridges.push_back(cuts.back());
                }
                previous = f;
                previous_z = z1;
            }
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
        const auto is_ridge = [&](double z) { return std::find(ridges.begin(), ridges.end(), z) != ridges.end(); };
        double sum = 0.0;
        for (std::size_t i = 0; i + 1 < cuts.size(); i++)
        {
            const double c = cuts[i], d = cuts[i + 1];
            if (!(d > c))
                continue;
            const bool sc = singular_a && c == a, sd = singular_b && d == b;
            for (int branch = 0; branch < 2; branch++)
            {
                const bool up = branch == 0;
                // Whether this branch lies in the rectangle on the piece.
                double z2;
                if (!point(0.5 * (c + d), level, up, z2) || z2 < two_.lo || z2 > two_.hi)
                    continue;
                sum += piece(c, d, sc, sd, is_ridge(c), is_ridge(d), level, up);
            }
        }
        return sum;
    }

private:
    // The points of one edge's stretch where a term of the other edge's
    // directivity has its pole on the real line: nu phi at 2 pi m, phi = pi +-
    // theta_S +- theta of the point (over_second: the second edge's points
    // about the first, with the source; else the first's about the second,
    // with the receiver).
    std::vector<double> boundary_points(bool over_second) const
    {
        const Edge& edge = over_second ? first_ : second_;
        const Edge& other = over_second ? second_ : first_;
        const Stretch& along = over_second ? two_ : one_;
        const double fixed = over_second ? s_.theta : r_.theta;
        const double nu = 1.0 / edge.n;
        const auto offsets = [&](double z, double* out) {
            const double t = about(edge, other.start + z * other.direction).theta;
            const double phis[4] = {pi + fixed + t, pi + fixed - t, pi - fixed + t, pi - fixed - t};
            for (int i = 0; i < 4; i++)
                out[i] = std::remainder(nu * phis[i], 2.0 * pi);
        };
        std::vector<double> points;
        const int samples = 64;
        double previous[4], current[4];
        offsets(along.lo, previous);
        double z_previous = along.lo;
        for (int j = 1; j <= samples; j++)
        {
            const double z = along.lo + (along.hi - along.lo) * j / samples;
            offsets(z, current);
            for (int i = 0; i < 4; i++)
                // A sign change through zero, not a wrap through +-pi.
                if ((previous[i] > 0.0) != (current[i] > 0.0) && std::abs(previous[i] - current[i]) < pi)
                {
                    double lo = z_previous, hi = z, flo = previous[i];
                    for (int n = 0; n < 100; n++)
                    {
                        const double c = 0.5 * (lo + hi);
                        double f[4];
                        offsets(c, f);
                        if ((f[i] > 0.0) == (flo > 0.0))
                        {
                            lo = c;
                            flo = f[i];
                        }
                        else
                            hi = c;
                    }
                    points.push_back(0.5 * (lo + hi));
                }
            std::copy(current, current + 4, previous);
            z_previous = z;
        }
        return points;
    }

    // The second edge's point at the level for the first edge's point z1, on
    // the branch, and the integrand there.
    bool point(double z1, double level, bool up, double& z2) const
    {
        const Vec3 p1 = first_.start + z1 * first_.direction;
        const About q = about(second_, p1);
        if (!(q.r > 0.0))
            return false;
        const Geometry g = geometry({q.r, q.theta, q.z, r_.r, r_.theta, r_.z});
        const double rho = level - (source_ - p1).norm();
        const double s2 = (rho - g.length0) * (rho + g.length0) / (4.0 * g.rr);
        if (!(s2 > 0.0))
            return false;
        z2 = g.point(2.0 * std::asinh(std::sqrt(s2)), up);
        return true;
    }
    double integrand(double z1, double level, bool up) const
    {
        const Vec3 p1 = first_.start + z1 * first_.direction;
        const double m = (source_ - p1).norm();
        const About q = about(second_, p1);
        if (!(q.r > 0.0))
            return 0.0;
        const Geometry g = geometry({q.r, q.theta, q.z, r_.r, r_.theta, r_.z});
        const double rho = level - m;
        const double s2 = (rho - g.length0) * (rho + g.length0) / (4.0 * g.rr);
        if (!(s2 > 0.0))
            return 0.0;
        const double eta2 = 2.0 * std::asinh(std::sqrt(s2));
        const double z2 = g.point(eta2, up);
        if (z2 < two_.lo || z2 > two_.hi)
            return 0.0;
        // sinh(eta2) = 2 sinh(eta2/2) cosh(eta2/2).
        const double sinh2 = 2.0 * std::sqrt(s2) * std::sqrt(1.0 + s2);
        const About t = about(first_, second_.start + z2 * second_.direction);
        if (!(t.r > 0.0))
            return 0.0;
        const Geometry g1 = geometry({s_.r, s_.theta, s_.z, t.r, t.theta, t.z});
        const double b1 = parallel_ ? beta1_(g1.eta(z1))
                                    : directivity(first_.n, p1_, source_half_, half_angle(first_.n, t.theta), s_.theta,
                                              t.theta)(g1.eta(z1));
        const double b2 = parallel_ ? beta2_(eta2)
                                    : directivity(second_.n, p2_, half_angle(second_.n, q.theta), receiver_half_,
                                              q.theta, r_.theta)(eta2);
        return b1 * b2 / (m * g.rr * sinh2);
    }
    // z1 less the first edge's apex for the source and the branch's point of
    // the second edge: zero on beta1's ridge.  NaN off the curve.
    double ridge(double z1, double level, bool up) const
    {
        double z2;
        if (!point(z1, level, up, z2))
            return std::numeric_limits<double>::quiet_NaN();
        const About t = about(first_, second_.start + z2 * second_.direction);
        return z1 - (s_.z + (t.z - s_.z) * s_.r / (s_.r + t.r));
    }
    // The integral over [c, d] of z1, a singular end (where the curve turns,
    // as 1/sqrt) taken out by z1 = c + (d - c) w(t).
    double piece(double c, double d, bool sc, bool sd, bool rc, bool rd, double level, bool up) const
    {
        const double span = d - c;
        const auto mapped = [&](double t) -> complex {
            double w, dw;
            if (sc && sd)
            {
                w = 0.5 - 0.5 * std::cos(pi * t);
                dw = 0.5 * pi * std::sin(pi * t);
            }
            else if (sc)
            {
                w = t * t;
                dw = 2.0 * t;
            }
            else if (sd)
            {
                w = 1.0 - (1.0 - t) * (1.0 - t);
                dw = 2.0 * (1.0 - t);
            }
            else
            {
                w = t;
                dw = 1.0;
            }
            return integrand(c + span * w, level, up) * span * dw;
        };
        // Five panels, graded toward the ends; toward a ridge end down to
        // 1e-10, where beta1's peak narrows without limit as the level nears a
        // point where its pole meets the ridge.  Gauss-Kronrod's |K - G|
        // overstates K's own error by orders: 1e-5 of it leaves K near 1e-8.
        std::vector<double> cuts {0.0};
        if (rc)
            for (double x = 1e-8; x < 0.04; x *= 30.0)
                cuts.push_back(x);
        if (rc || rd || ((sc || sd) && narrow2_ < 0.3))
            cuts.insert(cuts.end(), {0.08, 0.5, 0.92});
        else
            cuts.push_back(0.5);
        if (rd)
        {
            std::vector<double> top;
            for (double x = 1e-8; x < 0.04; x *= 30.0)
                top.push_back(1.0 - x);
            cuts.insert(cuts.end(), top.rbegin(), top.rend());
        }
        cuts.push_back(1.0);
        return adaptive(mapped, cuts, 1e-5, 0.0, 400).real();
    }

    const Edge& first_;
    const Edge& second_;
    Vec3 source_, receiver_;
    Stretch one_, two_;
    About s_, r_;
    bool parallel_ = false;
    Directivity beta1_, beta2_;
    HalfAngle p1_, p2_, source_half_, receiver_half_;
    double narrow1_ = 1.0, narrow2_ = 1.0;
    double psi_least_ = 0.0;
    double side_least_[2] = {0.0, 0.0};
};

// Chebyshev interpolation of u G on [B + u0^2, B + u1^2] in u = sqrt(L - B):
// there u G is smooth, the square root G has where a curve first meets a side
// or shrinks to a point, and the 1/sqrt it has past a boundary point, both
// taken out by the variable and the factor u.
struct Interpolant
{
    double base = 0.0, u0 = 0.0, u1 = 0.0;
    std::vector<double> u, value, weight;  // value: u G

    double operator()(double level) const
    {
        const double x = std::sqrt(std::max(level - base, 0.0));
        double num = 0.0, den = 0.0;
        for (std::size_t i = 0; i < u.size(); i++)
        {
            const double d = x - u[i];
            if (d == 0.0)
                return value[i] / x;
            const double w = weight[i] / d;
            num += w * value[i];
            den += w;
        }
        return num / den / x;
    }
};

// Chebyshev points of the first kind, which avoid the ends: at u = 0 the curve
// is a point, where G is only a limit.
Interpolant chebyshev(const LevelCurves& g, double base, double u0, double u1, int count)
{
    Interpolant f;
    f.base = base;
    f.u0 = u0;
    f.u1 = u1;
    for (int i = 0; i < count; i++)
    {
        const double angle = pi * (2 * i + 1) / (2 * count);
        const double x = 0.5 * (u0 + u1) - 0.5 * (u1 - u0) * std::cos(angle);
        f.u.push_back(x);
        f.value.push_back(x * g(base + x * x));
        f.weight.push_back((i % 2 == 0 ? 1.0 : -1.0) * std::sin(angle));
    }
    return f;
}

} // namespace

std::vector<complex> twice_diffracted_spectrum(const Edge& first, const Edge& second, const Vec3& source,
        const Vec3& receiver, const Stretch& along_first, const Stretch& along_second, double plane_factor,
        const std::vector<double>& wavenumbers)
{
    std::vector<complex> out(wavenumbers.size(), 0.0);
    if (wavenumbers.empty() || !(along_first.hi > along_first.lo) || !(along_second.hi > along_second.lo))
        return out;
    const About s = about(first, source), r = about(second, receiver);
    if (!(s.r > 0.0) || !(r.r > 0.0))
        return out;
    const double k_max = *std::max_element(wavenumbers.begin(), wavenumbers.end());
    const LevelCurves curves(first, second, source, receiver, along_first, along_second);
    // What is held still near each break, and how near two breaks are one: a
    // hundredth of a radian at the largest wavenumber, and a thousandth of the
    // scale a directivity's peak has in the path near the least, where its
    // width in eta is psi, which is about psi^2 r / 4 for a leg r off its edge.
    const double peak = 0.25 * std::pow(curves.narrowest(), 2) * std::min({s.r, r.r, curves.gap()});
    // ... and no finer than 1e-10 of the least path, round-off's.
    const double hold = std::max(std::min(1e-2 / k_max, 1e-3 * peak), 1e-10 * curves.least());
    const std::vector<double> breaks = curves.breaks(hold);
    const double scale = plane_factor / (first.n * second.n * 16.0 * pi * pi);
    // The integral of |G| over all the segments, roughly, from three points
    // in each: what a piece's error is measured against, the transform seeing
    // all of them.  (Raised as the pieces find more.)
    double measure = 0.0;
    for (std::size_t i = 0; i + 1 < breaks.size(); i++)
    {
        const double span = breaks[i + 1] - breaks[i];
        double largest = 0.0;
        for (double t : {0.25, 0.5, 0.75})
            largest = std::max(largest, std::abs(curves(breaks[i] + t * t * span)));
        measure += largest * span;
    }
    // Per segment between breaks: G in u = sqrt(L - B), interpolated on
    // Chebyshev pieces, halved until the interpolant at the pieces' middles
    // agrees with G; then the transform, the first stretch (under a tenth of
    // a radian at the largest wavenumber) by Gauss-Legendre in u, the rest by
    // Filon panels geometric in L - B.
    for (std::size_t i = 0; i + 1 < breaks.size(); i++)
    {
        const double base = breaks[i], top = breaks[i + 1];
        const double span = top - base;
        // The interpolants.
        // A piece is good when the interpolant's error, over its width in L,
        // is small against the segment's integral of |G|: what the transform
        // sees.  So a sharp layer at a break -- where a boundary point lies, G
        // steps through a layer microns wide -- is taken coarsely.
        // Below u_min = 1e-3 of the segment's u, a millionth of its span in L,
        // G is held at its value there: its weight below is a millionth of
        // the segment's, and at a boundary point the layer G steps through
        // lies there.  From the least path, where a grazing transition's
        // structure lies on the directivity peak's own scale, u_min is no more
        // than the hold's root, and the pieces start geometric in u, a factor
        // of eight apiece: what G does near the least is on the scale of how
        // near.
        std::vector<Interpolant> pieces;
        const bool first = i == 0;
        const double u_min = first ? std::min(1e-3 * std::sqrt(span), std::sqrt(hold)) : 1e-3 * std::sqrt(span);
        std::vector<std::pair<double, double>> todo;
        if (first)
            for (double u = std::sqrt(span); u > u_min;)
            {
                const double next = std::max(0.125 * u, u_min);
                todo.push_back({next, u});
                u = next;
            }
        else
            todo.push_back({u_min, std::sqrt(span)});
        while (!todo.empty())
        {
            const auto [u0, u1] = todo.back();
            todo.pop_back();
            // Five points, checked at one; nine, checked at two; seventeen,
            // checked at four, where nine fall short.
            const double width = u1 * u1 - u0 * u0;
            Interpolant f;
            double error = 0.0;
            for (int count : {9, 17})
            {
                f = chebyshev(curves, base, u0, u1, count);
                double largest = 0.0;
                for (std::size_t j = 0; j < f.value.size(); j++)
                    largest = std::max(largest, std::abs(f.value[j] / f.u[j]));
                measure = std::max(measure, largest * width);
                error = 0.0;
                const std::vector<double> checks = count == 5 ? std::vector<double> {0.62}
                        : count == 9                            ? std::vector<double> {0.31, 0.77}
                                                                : std::vector<double> {0.137, 0.42, 0.71, 0.93};
                for (double t : checks)
                {
                    const double x = u0 + t * (u1 - u0);
                    error = std::max(error, std::abs(f(base + x * x) - curves(base + x * x)));
                }
                if (error * width <= 1e-6 * measure)
                    break;
            }
            if (error * width <= 1e-6 * measure || u1 - u0 < 1e-4 * u1 || pieces.size() > 400)
                pieces.push_back(std::move(f));
            else
            {
                todo.push_back({0.5 * (u0 + u1), u1});
                todo.push_back({u0, 0.5 * (u0 + u1)});
            }
        }
        std::sort(pieces.begin(), pieces.end(), [](const Interpolant& a, const Interpolant& b) { return a.u0 < b.u0; });
        const auto g_at = [&](double level) {
            const double x = std::max(std::sqrt(std::max(level - base, 0.0)), u_min);
            for (const Interpolant& f : pieces)
                if (x <= f.u1)
                    return f(base + x * x);
            return pieces.back()(base + x * x);
        };
        // The first stretch: Gauss-Legendre in u over [0, u_first].
        const double u_first = std::min(std::sqrt(0.1 / k_max), std::sqrt(span));
        std::vector<double> g_nodes, l_nodes;
        for (int q = 0; q < 8; q++)
        {
            const double x = 0.5 * u_first * (legendre_nodes[q] + 1.0);
            g_nodes.push_back(0.5 * u_first * legendre_weights[q] * g_at(base + x * x) * 2.0 * x);
            l_nodes.push_back(base + x * x);
        }
        // Filon panels from there.
        std::vector<FilonPanel> panels;
        double a = base + u_first * u_first;
        double width = (a - base) * (spectrum_ratio - 1.0);
        while (a < top * (1.0 - 1e-15))
        {
            const double w = std::min(width, top - a);
            std::array<double, filon_count> values;
            for (int j = 0; j < filon_count; j++)
                values[j] = g_at(a + filon.x[j] * w);
            panels.push_back(filon_panel(a, w, values));
            a += w;
            width *= spectrum_ratio;
        }
        for (std::size_t n = 0; n < wavenumbers.size(); n++)
        {
            const double k = wavenumbers[n];
            complex sum = 0.0;
            for (std::size_t q = 0; q < g_nodes.size(); q++)
                sum += g_nodes[q] * std::polar(1.0, k * l_nodes[q]);
            for (const FilonPanel& p : panels)
                sum += filon_integral(p, k);
            out[n] += scale * sum;
        }
    }
    return out;
}

std::vector<complex> localized_twice_diffracted_spectrum(const Edge& first, const Edge& second, const Vec3& source,
        const Vec3& receiver, double plane_factor, const std::vector<double>& wavenumbers, PartnerModel model)
{
    const Stretch whole1 {0.0, first.length}, whole2 {0.0, second.length};
    std::vector<complex> sum = twice_diffracted_spectrum(first, second, source, receiver, whole1, whole2, plane_factor,
            wavenumbers);
    for (const PartnerWave& p : twice_diffracted_partners(first, second, source, receiver))
    {
        const Edge& edge = p.over_second ? second : first;
        const Vec3& s = p.over_second ? p.image : source;
        const Vec3& r = p.over_second ? receiver : p.image;
        const EdgeFrameEnds ends = edge_frame_ends(edge, s, r);
        const std::vector<complex> part = edge_source_spectrum(edge.n, ends, p.part.lo, p.part.hi, wavenumbers);
        // The partner's own ray, where its Keller point lies on the part: the
        // model's.
        const double apex = ends.z_s + (ends.z_r - ends.z_s) * ends.r_s / (ends.r_s + ends.r_r);
        const bool lit = apex >= p.part.lo && apex <= p.part.hi && apex >= 0.0 && apex <= edge.length;
        std::vector<complex> own;
        if (lit && model == PartnerModel::edge_source)
            own = edge_source_spectrum(edge.n, ends, 0.0, edge.length, wavenumbers);
        for (std::size_t i = 0; i < wavenumbers.size(); i++)
        {
            sum[i] += part[i];
            if (lit)
                sum[i] -= model == PartnerModel::edge_source
                        ? own[i]
                        : keller_once_field(edge, s, r, wavenumbers[i], 0.0, edge.length);
        }
    }
    return sum;
}

} // namespace niceops
