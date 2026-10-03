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

    // within_tolerance with the ends freed from the natural condition.
    //
    // A natural spline has zero curvature at both ends, and a smoothing spline's
    // boundary bias is order one for about a smoothing length from each: with the
    // stiff fits the trajectory code asks for (a bandwidth of ten seconds) the last
    // several seconds of a track have their second derivative -- the acceleration
    // -- pulled towards zero.  On a track that ends in a deceleration it cut the
    // last four seconds' deceleration by 43 percent, left the last state at 5.4 m/s
    // where the aircraft had stopped, and tilted the thrust vector 7 to 11 degrees
    // wrong, which moved the footprint under the landing end by up to 1.8 dB.
    //
    // The cure is to give the fit something to bend into: each end is extended by
    // `extension` (same units as x) of virtual points on the quadratic that
    // least-squares fits the last `window` of the data, at the data's mean spacing;
    // the spline is fitted through all of them, and only the real range is
    // evaluated.  `rms` is a root-mean-square residual per point (so the budget
    // keeps the meaning it has on the plain fit): the tolerance is rms^2 times the
    // number of points, extension included.  An end with too few points in its
    // window (under 8) or a window that is not spanned is left natural.
    //
    // Returns a spline over the extended range; evaluate it inside [x.front(),
    // x.back()].  residual() and fitted() refer to the extended points.
    static SmoothingSpline within_tolerance_extended(const std::vector<double>& x,
            const std::vector<double>& y, double rms, double window, double extension);

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

    // Build a spline for one lambda from already-validated inputs and an
    // explicit weight vector.  Both factory functions go through this, so
    // within_tolerance does not re-validate and re-copy the weights for the
    // final fit it has already been holding all along.
    static SmoothingSpline fit(const std::vector<double>& x, const std::vector<double>& y,
            const std::vector<double>& w, double lambda);
};

} // namespace niceops
