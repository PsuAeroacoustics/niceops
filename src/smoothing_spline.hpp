#pragma once

// Natural cubic smoothing spline, replacing the vendored SPLINTER library.
//
// Fits the function f minimising
//
//     sum_i w_i (y_i - f(x_i))^2  +  lambda * integral f''(t)^2 dt
//
// over the data abscissae.  lambda = 0 interpolates; lambda -> infinity gives
// the weighted least-squares straight line, since only a line has zero
// curvature everywhere.  Between those, larger lambda buys smoothness at the
// cost of fidelity.
//
// Two ways to choose it:
//
//   with_smoothing(x, y, lambda)      you pick lambda directly.
//   within_tolerance(x, y, tolerance) you state how far the curve may sit from
//                                     the data, and lambda is solved for.  This
//                                     is MATLAB's spaps(): the smoothest curve
//                                     whose weighted squared residual does not
//                                     exceed `tolerance`.
//
// The tolerance form is usually the one you want for measured data, because the
// bound is in the units of the measurement -- "stay within the noise" -- rather
// than an abstract knob whose meaning shifts with sampling rate and units.
//
// Implementation is Reinsch's algorithm: the interior second derivatives solve a
// symmetric pentadiagonal system, factored in O(n).

#include <cstddef>
#include <vector>

namespace niceops {

class SmoothingSpline
{
public:
    // Fit with an explicit smoothing parameter.  `weights` may be empty, meaning
    // all ones; otherwise it must match x and y in length and be positive.
    // Requires x strictly increasing and at least two points.
    static SmoothingSpline with_smoothing(const std::vector<double>& x,
            const std::vector<double>& y, double lambda,
            const std::vector<double>& weights = {});

    // Fit the smoothest curve whose weighted squared residual
    //     sum_i w_i (y_i - f(x_i))^2
    // is at most `tolerance` (MATLAB spaps).  A tolerance of zero interpolates.
    // A tolerance at or above the residual of the least-squares straight line
    // returns that line, which is the smoothest curve available.
    static SmoothingSpline within_tolerance(const std::vector<double>& x,
            const std::vector<double>& y, double tolerance,
            const std::vector<double>& weights = {});

    // Value at xi.  Outside [x.front(), x.back()] the end cubics are extended,
    // which for a natural spline means linear extrapolation.
    double operator()(double xi) const;

    // order 0, 1 or 2: value, slope or curvature.
    double derivative(double xi, int order) const;

    std::vector<double> operator()(const std::vector<double>& xi) const;

    // The smoothing parameter used -- the one supplied, or the one solved for.
    double smoothing_parameter() const
    {
        return lambda_;
    }

    // Achieved weighted squared residual, sum_i w_i (y_i - f(x_i))^2.
    double residual() const
    {
        return residual_;
    }

    // Fitted values at the data abscissae.
    const std::vector<double>& fitted() const
    {
        return fitted_;
    }

private:
    SmoothingSpline() = default;

    std::vector<double> knots_;        // x
    std::vector<double> fitted_;       // f(x_i), NOT the raw data
    std::vector<double> curvature_;    // f''(x_i); zero at both ends
    double lambda_ = 0.0;
    double residual_ = 0.0;

    std::size_t interval_for(double xi) const;
};

} // namespace niceops
