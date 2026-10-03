// Natural cubic smoothing spline.
//
// The checks here are the defining properties of the object rather than
// hard-coded numbers, so they pin the behavior without encoding one particular
// implementation:
//
//   lambda = 0            interpolates the data exactly
//   lambda -> infinity    collapses to the weighted least-squares straight line
//   linear data           reproduced exactly at every lambda (a line is already
//                         curvature-free, so the penalty term costs it nothing)
//   residual(lambda)      monotonically non-decreasing
//   within_tolerance(t)   achieves a residual at or just under t
//
// Continuity of value, slope and curvature across knots is checked directly,
// since a discontinuity there is the classic symptom of a botched cubic
// assembly and would not show up in the residual.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numeric>
#include <vector>

#include "smoothing_spline.hpp"

using namespace niceops;

namespace {

std::vector<double> linspace(double a, double b, std::size_t n)
{
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; i++)
        v[i] = a + (b - a) * static_cast<double>(i) / static_cast<double>(n - 1);
    return v;
}

// A deterministic wiggle, so "noisy" data is reproducible.
std::vector<double> noisy_curve(const std::vector<double>& x)
{
    std::vector<double> y(x.size());
    for (std::size_t i = 0; i < x.size(); i++)
        y[i] = std::sin(x[i]) + 0.05 * std::sin(37.0 * x[i] + 1.0);
    return y;
}

double max_abs_difference(const std::vector<double>& a, const std::vector<double>& b)
{
    double worst = 0.0;
    for (std::size_t i = 0; i < a.size(); i++)
        worst = std::max(worst, std::fabs(a[i] - b[i]));
    return worst;
}

} // namespace

TEST_CASE("zero smoothing interpolates the data", "[spline]")
{
    const auto x = linspace(0.0, 6.0, 25);
    const auto y = noisy_curve(x);
    const auto spline = SmoothingSpline::with_smoothing(x, y, 0.0);

    REQUIRE(max_abs_difference(spline.fitted(), y) < 1e-12);
    REQUIRE(spline.residual() < 1e-20);
    for (std::size_t i = 0; i < x.size(); i++)
        REQUIRE(std::fabs(spline(x[i]) - y[i]) < 1e-10);
}

TEST_CASE("heavy smoothing collapses to the least-squares line", "[spline]")
{
    const auto x = linspace(0.0, 6.0, 25);
    const auto y = noisy_curve(x);
    const auto spline = SmoothingSpline::with_smoothing(x, y, 1e12);

    // Compare against the ordinary least-squares line computed independently.
    const double n = static_cast<double>(x.size());
    const double sx = std::accumulate(x.begin(), x.end(), 0.0);
    const double sy = std::accumulate(y.begin(), y.end(), 0.0);
    double sxx = 0.0, sxy = 0.0;
    for (std::size_t i = 0; i < x.size(); i++)
    {
        sxx += x[i] * x[i];
        sxy += x[i] * y[i];
    }
    const double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    const double intercept = (sy - slope * sx) / n;

    for (std::size_t i = 0; i < x.size(); i++)
        REQUIRE(std::fabs(spline(x[i]) - (intercept + slope * x[i])) < 1e-6);

    // A line has no curvature anywhere.
    REQUIRE(std::fabs(spline.derivative(3.0, 2)) < 1e-6);
}

TEST_CASE("linear data is reproduced at every smoothing level", "[spline]")
{
    const auto x = linspace(-2.0, 5.0, 17);
    std::vector<double> y(x.size());
    for (std::size_t i = 0; i < x.size(); i++)
        y[i] = 3.5 * x[i] - 1.25;

    for (double lambda : {0.0, 1e-6, 1.0, 1e3, 1e9})
    {
        const auto spline = SmoothingSpline::with_smoothing(x, y, lambda);
        REQUIRE(max_abs_difference(spline.fitted(), y) < 1e-9);
        REQUIRE(std::fabs(spline.derivative(1.0, 1) - 3.5) < 1e-8);
    }
}

TEST_CASE("residual grows monotonically with the smoothing parameter", "[spline]")
{
    const auto x = linspace(0.0, 6.0, 40);
    const auto y = noisy_curve(x);

    double previous = -1.0;
    for (double lambda : {0.0, 1e-6, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6})
    {
        const auto spline = SmoothingSpline::with_smoothing(x, y, lambda);
        REQUIRE(spline.residual() >= previous - 1e-12);
        previous = spline.residual();
    }
}

TEST_CASE("value, slope and curvature are continuous across knots", "[spline]")
{
    const auto x = linspace(0.0, 6.0, 20);
    const auto y = noisy_curve(x);
    const auto spline = SmoothingSpline::with_smoothing(x, y, 1e-3);

    const double eps = 1e-7;
    for (std::size_t i = 1; i + 1 < x.size(); i++)
    {
        for (int order = 0; order <= 2; order++)
        {
            const double left = spline.derivative(x[i] - eps, order);
            const double right = spline.derivative(x[i] + eps, order);
            REQUIRE(std::fabs(left - right) < 1e-4);
        }
    }
}

TEST_CASE("natural boundary: curvature vanishes at the ends", "[spline]")
{
    const auto x = linspace(0.0, 6.0, 20);
    const auto y = noisy_curve(x);
    const auto spline = SmoothingSpline::with_smoothing(x, y, 1e-2);

    REQUIRE(std::fabs(spline.derivative(x.front(), 2)) < 1e-10);
    REQUIRE(std::fabs(spline.derivative(x.back(), 2)) < 1e-10);
}

TEST_CASE("analytic derivatives agree with finite differences", "[spline]")
{
    const auto x = linspace(0.0, 6.0, 30);
    const auto y = noisy_curve(x);
    const auto spline = SmoothingSpline::with_smoothing(x, y, 1e-3);

    const double h = 1e-5;
    for (double xi : {0.7, 2.2, 3.9, 5.1})
    {
        const double numeric_slope = (spline(xi + h) - spline(xi - h)) / (2.0 * h);
        REQUIRE(std::fabs(spline.derivative(xi, 1) - numeric_slope) < 1e-4);

        const double numeric_curvature = (spline(xi + h) - 2.0 * spline(xi) + spline(xi - h)) / (h * h);
        REQUIRE(std::fabs(spline.derivative(xi, 2) - numeric_curvature) < 1e-2);
    }
}

TEST_CASE("within_tolerance honors the residual bound", "[spline][spaps]")
{
    const auto x = linspace(0.0, 6.0, 60);
    const auto y = noisy_curve(x);

    for (double tolerance : {1e-4, 1e-3, 1e-2, 1e-1})
    {
        const auto spline = SmoothingSpline::within_tolerance(x, y, tolerance);
        // Satisfied, and not wastefully tight -- it should sit close to the bound,
        // because any slack means a smoother curve was available.
        REQUIRE(spline.residual() <= tolerance * (1.0 + 1e-6));
        REQUIRE(spline.residual() > tolerance * 0.5);
        REQUIRE(spline.smoothing_parameter() > 0.0);
    }
}

TEST_CASE("a looser tolerance yields a smoother curve", "[spline][spaps]")
{
    const auto x = linspace(0.0, 6.0, 60);
    const auto y = noisy_curve(x);

    const auto tight = SmoothingSpline::within_tolerance(x, y, 1e-4);
    const auto loose = SmoothingSpline::within_tolerance(x, y, 1e-1);

    REQUIRE(loose.smoothing_parameter() > tight.smoothing_parameter());
    REQUIRE(loose.residual() > tight.residual());

    // Total absolute curvature at the knots is a serviceable roughness measure.
    auto roughness = [&](const SmoothingSpline& s) {
        double total = 0.0;
        for (double xi : x)
            total += std::fabs(s.derivative(xi, 2));
        return total;
    };
    REQUIRE(roughness(loose) < roughness(tight));
}

TEST_CASE("zero tolerance interpolates", "[spline][spaps]")
{
    const auto x = linspace(0.0, 6.0, 20);
    const auto y = noisy_curve(x);
    const auto spline = SmoothingSpline::within_tolerance(x, y, 0.0);
    REQUIRE(max_abs_difference(spline.fitted(), y) < 1e-12);
}

TEST_CASE("a tolerance beyond the line's residual returns the line", "[spline][spaps]")
{
    const auto x = linspace(0.0, 6.0, 20);
    const auto y = noisy_curve(x);
    const auto spline = SmoothingSpline::within_tolerance(x, y, 1e9);
    // Nothing is smoother than a straight line, so curvature should vanish.
    REQUIRE(std::fabs(spline.derivative(3.0, 2)) < 1e-6);
}

TEST_CASE("weights steer the fit toward trusted points", "[spline]")
{
    auto x = linspace(0.0, 6.0, 21);
    auto y = noisy_curve(x);
    // Corrupt one sample badly, then tell the fit not to believe it.
    const std::size_t bad = 10;
    y[bad] += 5.0;

    std::vector<double> weights(x.size(), 1.0);
    weights[bad] = 1e-6;

    const auto trusting = SmoothingSpline::with_smoothing(x, y, 1e-3);
    const auto skeptical = SmoothingSpline::with_smoothing(x, y, 1e-3, weights);

    const double clean = std::sin(x[bad]);
    REQUIRE(std::fabs(skeptical(x[bad]) - clean) < std::fabs(trusting(x[bad]) - clean));
}

TEST_CASE("two points give the straight line through them", "[spline]")
{
    const std::vector<double> x = {1.0, 3.0};
    const std::vector<double> y = {2.0, 8.0};
    const auto spline = SmoothingSpline::with_smoothing(x, y, 1.0);
    REQUIRE(std::fabs(spline(2.0) - 5.0) < 1e-12);
    REQUIRE(std::fabs(spline.derivative(2.0, 1) - 3.0) < 1e-12);
}

TEST_CASE("evaluation outside the data extends linearly", "[spline]")
{
    const auto x = linspace(0.0, 6.0, 20);
    const auto y = noisy_curve(x);
    const auto spline = SmoothingSpline::with_smoothing(x, y, 1e-3);

    // A natural spline has zero end curvature, so extrapolation is a straight
    // line continuing the end slope.
    const double slope = spline.derivative(x.back(), 1);
    const double predicted = spline(x.back()) + 2.0 * slope;
    REQUIRE(std::fabs(spline(x.back() + 2.0) - predicted) < 1e-6);
}

TEST_CASE("zero smoothing gives the interpolating cubic, not a polyline", "[spline]")
{
    // lambda = 0 interpolates, which is true of the values and was once taken
    // to settle the curvature too: fit() short-circuited and left it zero, so
    // the "interpolating spline" was piecewise linear.  Its first derivative
    // came back as the cell's chord slope and its second as identically zero,
    // which is what the geodetic reader's default --path_smoothing 0 was
    // handing the load factor.
    //
    // A parabola is the case that says so: an interpolating cubic carries a
    // constant second derivative exactly, away from the ends where the natural
    // boundary condition pins it to zero.
    const auto x = linspace(0.0, 4.0, 81);
    std::vector<double> y(x.size());
    for (std::size_t i = 0; i < x.size(); i++)
        y[i] = x[i] * x[i];

    for (const auto& spline : {SmoothingSpline::with_smoothing(x, y, 0.0),
                 SmoothingSpline::within_tolerance(x, y, 0.0)})
    {
        REQUIRE(spline.smoothing_parameter() == 0.0);
        // Interpolation is still exact.
        REQUIRE(max_abs_difference(spline.fitted(), y) < 1e-12);
        for (double xi : {1.0, 2.0, 3.0})
        {
            REQUIRE(std::fabs(spline.derivative(xi, 1) - 2.0 * xi) < 1e-6);
            REQUIRE(std::fabs(spline.derivative(xi, 2) - 2.0) < 1e-6);
        }
    }
}
