#pragma once

// Trajectory smoothing, built on the natural cubic smoothing spline in
// smoothing_spline.hpp.
//
// This used to wrap the vendored SPLINTER library through an extern "C" shim
// (a leftover of the original D FFI boundary).  SPLINTER brought its own
// bundled Eigen 3.2.8, which predates C++17 and forced build workarounds, all
// to supply one penalised-B-spline fit.  The replacement is self-contained.
//
// Two ways to choose how hard to smooth:
//
//   Spline(x, f, lambda)                  an explicit smoothing parameter.
//   Spline::within_tolerance(x, f, tol)   state how far the curve may sit from
//                                         the data and let lambda be solved
//                                         for, as MATLAB's spaps() does.
//
// Prefer the tolerance form for measured data: the bound is in the units of the
// measurement, so "stay within the noise" is expressible directly, whereas a
// bare lambda changes meaning with sampling rate and units.

#include <vector>

#include "smoothing_spline.hpp"

namespace niceops {

class Spline
{
public:
    // Smooth with an explicit penalty weight.  lambda = 0 interpolates.
    Spline(const std::vector<double>& x, const std::vector<double>& f, double lambda = 0.3)
        : spline_(SmoothingSpline::with_smoothing(x, f, lambda))
    {
    }

    // Smooth just enough to stay within `tolerance` of the data, measured as
    // the weighted sum of squared residuals.
    static Spline within_tolerance(const std::vector<double>& x, const std::vector<double>& f,
            double tolerance, const std::vector<double>& weights = {})
    {
        return Spline(SmoothingSpline::within_tolerance(x, f, tolerance, weights));
    }

    // order 0, 1 or 2: value, slope or curvature.
    double eval(double xi, int order = 0) const
    {
        return spline_.derivative(xi, order);
    }

    // Convenience function to apply eval to entire array
    std::vector<double> eval(const std::vector<double>& xi, int order = 0) const
    {
        std::vector<double> out(xi.size());
        for (std::size_t i = 0; i < xi.size(); i++)
            out[i] = spline_.derivative(xi[i], order);
        return out;
    }

    // The penalty weight used -- supplied, or solved for by within_tolerance.
    double smoothing_parameter() const
    {
        return spline_.smoothing_parameter();
    }

    // Weighted sum of squared residuals actually achieved.
    double residual() const
    {
        return spline_.residual();
    }

private:
    explicit Spline(SmoothingSpline spline) : spline_(std::move(spline)) {}

    SmoothingSpline spline_;
};

} // namespace niceops
