#include "smoothing_spline.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

// Reinsch's construction.  With knot spacings h_i = x_{i+1} - x_i, the natural
// cubic smoothing spline's interior second derivatives gamma solve
//
//     (R + lambda * Q^T W^-1 Q) gamma = Q^T y
//     f = y - lambda * W^-1 Q gamma
//
// where Q is the (n x n-2) second-difference operator and R the (n-2 x n-2)
// tridiagonal Gram matrix of the hat functions.  Column j of Q (for interior
// knot j+1) has three nonzeros:
//
//     Q[j,   j] =  1/h_j
//     Q[j+1, j] = -1/h_j - 1/h_{j+1}
//     Q[j+2, j] =  1/h_{j+1}
//
// so R + lambda Q^T W^-1 Q is symmetric with bandwidth 2, and an LDL^T sweep
// solves it in O(n).

struct Bands
{
    std::vector<double> diagonal;   // size m
    std::vector<double> first_off;  // size m-1
    std::vector<double> second_off; // size m-2
};

// L has unit diagonal; l1/l2 are its two subdiagonals, d the pivots.  Passed in
// so a caller that solves repeatedly is not reallocating them every time.
struct Factorization
{
    std::vector<double> d, l1, l2;
};

// Solve a symmetric positive-definite pentadiagonal system in place via LDL^T.
// Returns false if the factorisation breaks down (non-positive pivot).
bool solve_pentadiagonal(const Bands& bands, Factorization& lu, std::vector<double>& rhs)
{
    const std::size_t m = bands.diagonal.size();
    if (m == 0)
        return true;

    lu.d.assign(m, 0.0);
    lu.l1.assign(m, 0.0);
    lu.l2.assign(m, 0.0);
    auto& d = lu.d;
    auto& l1 = lu.l1;
    auto& l2 = lu.l2;
    for (std::size_t i = 0; i < m; i++)
    {
        double value = bands.diagonal[i];
        if (i >= 1)
            value -= l1[i - 1] * l1[i - 1] * d[i - 1];
        if (i >= 2)
            value -= l2[i - 2] * l2[i - 2] * d[i - 2];
        if (!(value > 0.0) || !std::isfinite(value))
            return false;
        d[i] = value;

        if (i + 1 < m)
        {
            double off = bands.first_off[i];
            if (i >= 1)
                off -= l2[i - 1] * l1[i - 1] * d[i - 1];
            l1[i] = off / d[i];
        }
        if (i + 2 < m)
        {
            l2[i] = bands.second_off[i] / d[i];
        }
    }

    // Forward substitution (L z = rhs), then D, then back substitution (L^T).
    for (std::size_t i = 0; i < m; i++)
    {
        double value = rhs[i];
        if (i >= 1)
            value -= l1[i - 1] * rhs[i - 1];
        if (i >= 2)
            value -= l2[i - 2] * rhs[i - 2];
        rhs[i] = value;
    }
    for (std::size_t i = 0; i < m; i++)
        rhs[i] /= d[i];
    for (std::size_t i = m; i-- > 0;)
    {
        double value = rhs[i];
        if (i + 1 < m)
            value -= l1[i] * rhs[i + 1];
        if (i + 2 < m)
            value -= l2[i] * rhs[i + 2];
        rhs[i] = value;
    }
    return true;
}

void validate(const std::vector<double>& x, const std::vector<double>& y,
        const std::vector<double>& weights)
{
    if (x.size() != y.size())
        input_error("SmoothingSpline: x and y differ in length (", x.size(), " vs ", y.size(), ").");
    if (x.size() < 2)
        input_error("SmoothingSpline: need at least two points, got ", x.size(), ".");
    for (std::size_t i = 1; i < x.size(); i++)
    {
        if (!(x[i] > x[i - 1]))
            input_error("SmoothingSpline: x must be strictly increasing (x[", i - 1, "] = ", x[i - 1],
                    ", x[", i, "] = ", x[i], ").");
    }
    if (!weights.empty())
    {
        if (weights.size() != x.size())
            input_error("SmoothingSpline: weights length ", weights.size(), " does not match ", x.size(), ".");
        for (std::size_t i = 0; i < weights.size(); i++)
        {
            if (!(weights[i] > 0.0))
                input_error("SmoothingSpline: weights must be positive (weight[", i, "] = ", weights[i], ").");
        }
    }
}

// Weighted squared residual of the fit.
double weighted_residual(const std::vector<double>& y, const std::vector<double>& fitted,
        const std::vector<double>& w)
{
    double total = 0.0;
    for (std::size_t i = 0; i < y.size(); i++)
    {
        const double r = y[i] - fitted[i];
        total += w[i] * r * r;
    }
    return total;
}

// Everything in Reinsch's construction that does not depend on lambda: the knot
// spacings, the three nonzeros of each Q column, the Gram bands R, the penalty
// bands Q^T W^-1 Q and the right-hand side Q^T y.  within_tolerance refits some
// sixty-odd times hunting for its lambda, and each refit used to rebuild all of
// this -- six divisions per interior knot among the rest -- when only the
// scalar lambda had changed.  The scratch vectors are members too, so the
// refits do not reallocate them either.
//
// Holds references to y and w; both factory functions keep them alive for the
// lifetime of the system.
class PenalizedSystem
{
public:
    PenalizedSystem(const std::vector<double>& x, const std::vector<double>& y,
            const std::vector<double>& w)
        : y_(y)
        , w_(w)
        , n_(x.size())
    {
        if (n_ < 3)
            return;
        const std::size_t m = n_ - 2; // interior knots

        std::vector<double> h(n_ - 1);
        for (std::size_t i = 0; i + 1 < n_; i++)
            h[i] = x[i + 1] - x[i];

        // Column j of Q, rows j, j+1, j+2.
        a_.resize(m);
        b_.resize(m);
        c_.resize(m);
        for (std::size_t j = 0; j < m; j++)
        {
            a_[j] = 1.0 / h[j];
            b_[j] = -1.0 / h[j] - 1.0 / h[j + 1];
            c_[j] = 1.0 / h[j + 1];
        }

        r_diag_.resize(m);
        r_off_.assign(m > 0 ? m - 1 : 0, 0.0);
        p_diag_.resize(m);
        p_off1_.assign(m > 0 ? m - 1 : 0, 0.0);
        p_off2_.assign(m > 1 ? m - 2 : 0, 0.0);
        qty_.resize(m);
        for (std::size_t j = 0; j < m; j++)
        {
            const double a = a_[j], b = b_[j], c = c_[j];
            r_diag_[j] = (h[j] + h[j + 1]) / 3.0;
            p_diag_[j] = a * a / w[j] + b * b / w[j + 1] + c * c / w[j + 2];
            if (j + 1 < m)
            {
                r_off_[j] = h[j + 1] / 6.0;
                p_off1_[j] = b * a_[j + 1] / w[j + 1] + c * b_[j + 1] / w[j + 2];
            }
            if (j + 2 < m)
            {
                p_off2_[j] = c * a_[j + 2] / w[j + 2];
            }
            qty_[j] = a * y[j] + b * y[j + 1] + c * y[j + 2]; // (Q^T y)_j
        }

        bands_.diagonal.resize(m);
        bands_.first_off.resize(m > 0 ? m - 1 : 0);
        bands_.second_off.resize(m > 1 ? m - 2 : 0);
        gamma_.resize(m);
    }

    // Fills `fitted` and `curvature` for this lambda.
    void fit(double lambda, std::vector<double>& fitted, std::vector<double>& curvature)
    {
        fitted = y_;
        curvature.assign(n_, 0.0);

        // Two points are already a straight line, with no interior knot to
        // carry curvature.
        //
        // lambda == 0 is NOT a special case and must not short-circuit here.
        // It interpolates, so `fitted` is y either way -- but the curvature
        // still has to be solved for, because the lambda = 0 limit of this
        // system is R gamma = Q^T y, the natural cubic *interpolating* spline.
        // Returning zero curvature instead made the "interpolating" spline
        // piecewise linear: derivative(x, 1) came back as the cell's chord
        // slope and derivative(x, 2) as identically zero.  That reached the
        // shipped geodetic reader, whose --path_smoothing default is 0 and
        // whose load factor is the second derivative, so every WGS84 track
        // flew at exactly 1 g and its CT never moved.
        if (n_ < 3)
            return;

        const std::size_t m = n_ - 2;
        for (std::size_t j = 0; j < m; j++)
        {
            // R + lambda * (Q^T W^-1 Q)
            bands_.diagonal[j] = r_diag_[j] + lambda * p_diag_[j];
            if (j + 1 < m)
                bands_.first_off[j] = r_off_[j] + lambda * p_off1_[j];
            if (j + 2 < m)
                bands_.second_off[j] = lambda * p_off2_[j];
        }
        gamma_ = qty_;

        if (!solve_pentadiagonal(bands_, lu_, gamma_))
        {
            input_error("SmoothingSpline: penalised system is not positive definite; check that x is "
                  "strictly increasing and weights are positive.");
        }

        // f = y - lambda * W^-1 Q gamma, and f'' = gamma at the interior knots.
        for (std::size_t i = 0; i < n_; i++)
        {
            double qg = 0.0;
            if (i >= 2)
                qg += c_[i - 2] * gamma_[i - 2];
            if (i >= 1 && i - 1 < m)
                qg += b_[i - 1] * gamma_[i - 1];
            if (i < m)
                qg += a_[i] * gamma_[i];
            fitted[i] = y_[i] - lambda * qg / w_[i];
        }
        for (std::size_t j = 0; j < m; j++)
            curvature[j + 1] = gamma_[j];
    }

private:
    const std::vector<double>& y_;
    const std::vector<double>& w_;
    std::size_t n_;
    std::vector<double> a_, b_, c_;                 // the nonzeros of Q's columns
    std::vector<double> r_diag_, r_off_;            // R
    std::vector<double> p_diag_, p_off1_, p_off2_;  // Q^T W^-1 Q
    std::vector<double> qty_;
    Bands bands_;                                   // scratch, rebuilt per lambda
    Factorization lu_;
    std::vector<double> gamma_;
};

} // namespace

namespace niceops {

SmoothingSpline SmoothingSpline::fit(const std::vector<double>& x, const std::vector<double>& y,
        const std::vector<double>& w, double lambda)
{
    SmoothingSpline spline;
    spline.knots_ = x;
    spline.lambda_ = lambda;
    PenalizedSystem(x, y, w).fit(lambda, spline.fitted_, spline.curvature_);
    spline.residual_ = weighted_residual(y, spline.fitted_, w);
    return spline;
}

SmoothingSpline SmoothingSpline::with_smoothing(const std::vector<double>& x,
        const std::vector<double>& y, double lambda, const std::vector<double>& weights)
{
    validate(x, y, weights);
    if (!(lambda >= 0.0) || !std::isfinite(lambda))
        input_error("SmoothingSpline: smoothing parameter must be finite and non-negative, got ", lambda, ".");

    // Only materialize the unit weights when the caller gave none; binding a
    // reference otherwise avoids copying the caller's vector.
    std::vector<double> unit_weights;
    if (weights.empty())
        unit_weights.assign(x.size(), 1.0);
    const std::vector<double>& w = weights.empty() ? unit_weights : weights;

    return fit(x, y, w, lambda);
}

SmoothingSpline SmoothingSpline::within_tolerance(const std::vector<double>& x,
        const std::vector<double>& y, double tolerance, const std::vector<double>& weights)
{
    validate(x, y, weights);
    if (!(tolerance >= 0.0) || !std::isfinite(tolerance))
        input_error("SmoothingSpline: tolerance must be finite and non-negative, got ", tolerance, ".");

    std::vector<double> unit_weights;
    if (weights.empty())
        unit_weights.assign(x.size(), 1.0);
    const std::vector<double>& w = weights.empty() ? unit_weights : weights;

    if (tolerance == 0.0 || x.size() < 3)
        return fit(x, y, w, 0.0);

    // The residual grows monotonically with lambda, from 0 (interpolation) to the
    // residual of the weighted least-squares line.  Bracket, then bisect on
    // log(lambda), which behaves far better than lambda itself: the useful range
    // spans many orders of magnitude and the response is roughly logarithmic.
    std::vector<double> fitted, curvature;
    PenalizedSystem system(x, y, w);

    auto residual_at = [&](double lambda) {
        system.fit(lambda, fitted, curvature);
        return weighted_residual(y, fitted, w);
    };

    // Expand upward until the residual exceeds the tolerance.
    double upper = 1.0;
    double upper_residual = residual_at(upper);
    int expansions = 0;
    while (upper_residual < tolerance && expansions < 200)
    {
        upper *= 10.0;
        upper_residual = residual_at(upper);
        expansions++;
    }
    if (upper_residual < tolerance)
    {
        // Even a nearly rigid line stays inside the tolerance; that line is the
        // smoothest answer available, so take it.
        trace("SmoothingSpline: tolerance ", tolerance, " exceeds the least-squares residual; "
              "returning the straight-line fit.");
        return fit(x, y, w, upper);
    }

    // Walk down if the first probe already overshot.  Keeping the residual in a
    // variable matters here: evaluating it in the loop condition refits the very
    // lambda that was just fitted, an entire wasted solve every pass.
    double lower = upper / 10.0;
    double lower_residual = residual_at(lower);
    int reductions = 0;
    while (lower_residual > tolerance && reductions < 200)
    {
        upper = lower;
        lower /= 10.0;
        lower_residual = residual_at(lower);
        reductions++;
    }

    // Bisect in log space.  Sixty halvings of a one-decade bracket is far below
    // double precision, so this is a tight bound, not an approximation.
    double log_lower = std::log(lower), log_upper = std::log(upper);
    for (int i = 0; i < 60; i++)
    {
        const double log_mid = 0.5 * (log_lower + log_upper);
        if (residual_at(std::exp(log_mid)) > tolerance)
            log_upper = log_mid;
        else
            log_lower = log_mid;
    }

    // Take the lower end, so the constraint is satisfied rather than just missed.
    return fit(x, y, w, std::exp(log_lower));
}

std::size_t SmoothingSpline::interval_for(double xi) const
{
    // Index of the cubic piece containing xi, clamped so that evaluation outside
    // the data extends the end pieces.
    const std::size_t n = knots_.size();
    if (xi <= knots_.front())
        return 0;
    if (xi >= knots_.back())
        return n - 2;
    const auto it = std::upper_bound(knots_.begin(), knots_.end(), xi);
    return static_cast<std::size_t>(it - knots_.begin()) - 1;
}

double SmoothingSpline::derivative(double xi, int order) const
{
    if (order < 0 || order > 2)
        input_error("SmoothingSpline: derivative order must be 0, 1 or 2, got ", order, ".");

    // Outside the data, continue linearly.  A natural spline has zero curvature
    // at its ends, so the line tangent there IS the spline's continuation -- and
    // extending the end cubic instead would let the curve diverge sharply just
    // past the last sample, which for a trajectory is worse than useless.
    if (xi < knots_.front() || xi > knots_.back())
    {
        const bool past_end = xi > knots_.back();
        const double edge = past_end ? knots_.back() : knots_.front();
        if (order == 2)
            return 0.0;
        const double slope = derivative(edge, 1);
        if (order == 1)
            return slope;
        return derivative(edge, 0) + (xi - edge) * slope;
    }

    const std::size_t i = interval_for(xi);
    const double h = knots_[i + 1] - knots_[i];
    const double t = xi - knots_[i];
    const double f0 = fitted_[i], f1 = fitted_[i + 1];
    const double c0 = curvature_[i], c1 = curvature_[i + 1];

    // Cubic through (f0, f1) with prescribed end curvatures.
    const double slope = (f1 - f0) / h - h * (2.0 * c0 + c1) / 6.0;
    const double cubic = (c1 - c0) / (6.0 * h);

    switch (order)
    {
    case 0:
        return f0 + t * (slope + t * (0.5 * c0 + t * cubic));
    case 1:
        return slope + t * (c0 + 3.0 * t * cubic);
    default:
        return c0 + 6.0 * t * cubic;
    }
}

double SmoothingSpline::operator()(double xi) const
{
    return derivative(xi, 0);
}

std::vector<double> SmoothingSpline::operator()(const std::vector<double>& xi) const
{
    std::vector<double> out(xi.size());
    for (std::size_t i = 0; i < xi.size(); i++)
        out[i] = derivative(xi[i], 0);
    return out;
}

} // namespace niceops
