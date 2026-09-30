#include "pade_propagator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "util/constants.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

// Solves A x = b in place by Gaussian elimination with partial pivoting; false
// if A is singular to working precision.
bool solve_linear(std::vector<complex>& a, std::vector<complex>& b)
{
    const std::size_t n = b.size();
    for (std::size_t col = 0; col < n; col++)
    {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < n; row++)
            if (std::norm(a[row * n + col]) > std::norm(a[pivot * n + col]))
                pivot = row;
        if (!(std::norm(a[pivot * n + col]) > 0.0))
            return false;
        if (pivot != col)
        {
            for (std::size_t k = 0; k < n; k++)
                std::swap(a[col * n + k], a[pivot * n + k]);
            std::swap(b[col], b[pivot]);
        }
        const complex inverse = 1.0 / a[col * n + col];
        for (std::size_t row = col + 1; row < n; row++)
        {
            const complex factor = a[row * n + col] * inverse;
            for (std::size_t k = col; k < n; k++)
                a[row * n + k] -= factor * a[col * n + k];
            b[row] -= factor * b[col];
        }
    }
    for (std::size_t row = n; row-- > 0;)
    {
        complex sum = b[row];
        for (std::size_t k = row + 1; k < n; k++)
            sum -= a[row * n + k] * b[k];
        b[row] = sum / a[row * n + row];
    }
    return true;
}

// A Chebyshev series and its derivative at t.
void chebyshev(const std::vector<complex>& coefficients, complex t, complex& value, complex& derivative)
{
    complex t_previous = 1.0, t_current = t; // T_0, T_1
    complex d_previous = 0.0, d_current = 1.0; // T_0', T_1'
    value = coefficients[0];
    derivative = 0.0;
    for (std::size_t k = 1; k < coefficients.size(); k++)
    {
        value += coefficients[k] * t_current;
        derivative += coefficients[k] * d_current;
        const complex t_next = 2.0 * t * t_current - t_previous;
        const complex d_next = 2.0 * t_current + 2.0 * t * d_current - d_previous;
        t_previous = t_current;
        t_current = t_next;
        d_previous = d_current;
        d_current = d_next;
    }
}

// All roots of a Chebyshev series at once, by the Aberth-Ehrlich iteration:
// Newton's step for each, corrected for the pull of the others.  False if it
// does not settle; the caller's check of the finished propagator is what
// vouches for the roots either way.
bool chebyshev_roots(const std::vector<complex>& coefficients, std::vector<complex>& roots)
{
    const std::size_t degree = coefficients.size() - 1;
    roots.resize(degree);
    // Starting points on a circle below the interval, where the propagator's
    // poles lie, at angles that are not symmetric about either axis.
    for (std::size_t k = 0; k < degree; k++)
        roots[k] = complex(0.0, -1.5)
                + 1.2 * std::polar(1.0, 2.0 * pi * (static_cast<double>(k) + 0.3) / static_cast<double>(degree));
    for (int iteration = 0; iteration < 500; iteration++)
    {
        double largest_step = 0.0;
        for (std::size_t k = 0; k < degree; k++)
        {
            complex value, derivative;
            chebyshev(coefficients, roots[k], value, derivative);
            if (value == 0.0)
                continue;
            const complex newton = value / derivative;
            complex repulsion = 0.0;
            for (std::size_t other = 0; other < degree; other++)
                if (other != k)
                    repulsion += 1.0 / (roots[k] - roots[other]);
            const complex step = newton / (1.0 - newton * repulsion);
            roots[k] -= step;
            largest_step = std::max(largest_step, std::abs(step) / (1.0 + std::abs(roots[k])));
        }
        if (!std::isfinite(largest_step))
            return false;
        // Rounding in a high-degree series leaves steps of 1e-13 or so that
        // never shrink further.
        if (largest_step < 1e-11)
            return true;
    }
    return false;
}

// Two nodes among the evanescent modes, where the exact propagator has decayed
// to nothing.  Without them the interpolant of a function of modulus one comes
// out all-pass -- |r| = 1 along the whole real axis -- with a phase that is
// right only on the accurate interval.  The modes steeper than max_angle are
// then neither absorbed nor carried off: their vertical group velocity, the
// phase's derivative, is a twentieth of the true one, so what the starter
// puts into them lingers near the source's height and arrives downrange as
// noise several dB strong.  Pinned to zero out here, r falls away below the
// interval instead and those modes die within a few steps.  This is Collins'
// stability constraint.  Two is what it takes: a third and fourth overshoot,
// |r| > 1, between the interval and the first of them.
const std::vector<double> evanescent_nodes = {-1.5, -6.0};

// The type (terms, terms) interpolant, or false if it does not meet the
// request.
bool fit(const PadeRequest& request, int terms, PadePropagator& out)
{
    const std::size_t n = static_cast<std::size_t>(terms);
    const std::size_t unknowns = 2 * n + 1; // p_0 .. p_n, q_1 .. q_n; q_0 = 1
    const double mid = 0.5 * (request.accurate_high + request.accurate_low);
    const double half = 0.5 * (request.accurate_high - request.accurate_low);

    // P(x_k) - f_k Q(x_k) = 0 at each node, in the Chebyshev basis of
    // t = (x - mid) / half: Chebyshev nodes on the accurate interval, and the
    // evanescent ones.  Out there the basis is enormous, so each row is scaled
    // to its largest entry.
    const std::size_t inside = unknowns > evanescent_nodes.size() ? unknowns - evanescent_nodes.size() : 1;
    if (inside + evanescent_nodes.size() != unknowns)
        return false; // too few terms to carry the evanescent nodes
    std::vector<complex> matrix(unknowns * unknowns), rhs(unknowns);
    std::vector<double> basis(n + 1);
    for (std::size_t k = 0; k < unknowns; k++)
    {
        const double x = k < inside
                ? mid + half * std::cos(pi * (static_cast<double>(k) + 0.5) / static_cast<double>(inside))
                : evanescent_nodes[k - inside];
        const double t = (x - mid) / half;
        const complex f = exact_propagator(x, request.sigma);
        basis[0] = 1.0;
        basis[1] = t;
        for (std::size_t j = 2; j <= n; j++)
            basis[j] = 2.0 * t * basis[j - 1] - basis[j - 2];
        const double scale = 1.0 / std::max(1.0, std::abs(basis[n]));
        for (std::size_t j = 0; j <= n; j++)
            matrix[k * unknowns + j] = scale * basis[j];
        for (std::size_t j = 1; j <= n; j++)
            matrix[k * unknowns + n + j] = -f * (scale * basis[j]);
        rhs[k] = f * scale;
    }
    if (!solve_linear(matrix, rhs))
        return false;
    std::vector<complex> numerator(rhs.begin(), rhs.begin() + static_cast<std::ptrdiff_t>(n) + 1);
    std::vector<complex> denominator(n + 1);
    denominator[0] = 1.0;
    std::copy(rhs.begin() + static_cast<std::ptrdiff_t>(n) + 1, rhs.end(), denominator.begin() + 1);
    if (!(std::norm(denominator[n]) > 0.0))
        return false;

    // Partial fractions: the denominator's roots, and P / Q' at each.
    std::vector<complex> roots;
    if (!chebyshev_roots(denominator, roots))
        return false;
    out.constant = numerator[n] / denominator[n];
    out.poles.resize(n);
    out.residues.resize(n);
    for (std::size_t j = 0; j < n; j++)
    {
        complex p, p_derivative, q, q_derivative;
        chebyshev(numerator, roots[j], p, p_derivative);
        chebyshev(denominator, roots[j], q, q_derivative);
        out.poles[j] = mid + half * roots[j];
        out.residues[j] = half * p / q_derivative;
        // Every solve is against Q - p_j, whose spectrum lies on and above the
        // real axis: a pole there, or near it, is a singular solve.
        if (!std::isfinite(std::norm(out.residues[j])) || !(out.poles[j].imag() < -1e-3))
            return false;
    }

    // The measured error, of the partial fractions actually used.
    out.error = 0.0;
    const int accuracy_samples = 64 * terms;
    for (int k = 0; k <= accuracy_samples; k++)
    {
        const double x = request.accurate_low
                + (request.accurate_high - request.accurate_low) * k / accuracy_samples;
        out.error = std::max(out.error, std::abs(out(x) - exact_propagator(x, request.sigma)));
    }
    if (!(out.error <= request.tolerance))
        return false;

    // And that it amplifies nothing anywhere in the spectrum.  The poles sit a
    // few hundredths and more below the axis, so r's features are that wide:
    // fine samples near the accurate interval, coarser ones far from it.
    out.largest_modulus = 0.0;
    const double near_low = std::min(request.accurate_low - 2.0, -3.0);
    for (const double fraction : {0.0, 0.05, 0.25, 0.5, 1.0})
    {
        const double imag = fraction * request.spectrum_imag;
        for (double x = request.accurate_high; x >= near_low; x -= 0.002)
            out.largest_modulus = std::max(out.largest_modulus, std::abs(out(complex(x, imag))));
        for (double x = near_low; x >= request.spectrum_low; x -= 0.02)
            out.largest_modulus = std::max(out.largest_modulus, std::abs(out(complex(x, imag))));
        if (request.spectrum_imag <= 0.0)
            break;
    }
    return out.largest_modulus <= 1.0 + 10.0 * request.tolerance;
}

} // namespace

std::complex<double> PadePropagator::operator()(std::complex<double> x) const
{
    complex sum = constant;
    for (std::size_t j = 0; j < poles.size(); j++)
        sum += residues[j] / (x - poles[j]);
    return sum;
}

std::complex<double> exact_propagator(std::complex<double> x, double sigma)
{
    // sqrt(1 + x) with its cut turned onto the negative imaginary axis:
    // exp(i pi/4) sqrt(-i (1 + x)).
    const complex root = std::polar(1.0, 0.25 * pi) * std::sqrt(complex(0.0, -1.0) * (1.0 + x));
    return std::exp(complex(0.0, sigma) * (root - 1.0));
}

PadePropagator pade_propagator(const PadeRequest& request)
{
    if (!(request.sigma > 0.0) || !(request.accurate_high > request.accurate_low)
            || !(request.accurate_low > -1.0))
        throw std::invalid_argument("pade_propagator: the step or the accurate interval is not usable");
    PadePropagator propagator;
    for (int terms = 1; terms <= request.max_terms; terms++)
        if (fit(request, terms, propagator))
            return propagator;
    throw std::runtime_error("pade_propagator: no stable approximation within the tolerance");
}

} // namespace niceops
