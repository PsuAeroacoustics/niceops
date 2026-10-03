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

#include <catch2/catch_approx.hpp>
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

namespace {

// A deterministic stand-in for measurement noise: a few hundred microns, white
// enough for a smoother, and repeatable.
double jitter(std::size_t i)
{
    return 5e-4 * std::sin(12.9898 * static_cast<double>(i) + 78.233 * std::sin(0.37 * static_cast<double>(i)));
}

} // namespace

TEST_CASE("freed ends recover the acceleration the natural condition pins to zero",
        "[spline][ends]")
{
    // y = v0 t + a t^2 / 2 at 50 Hz for 131 s, smoothed to a 3.048 m RMS budget as
    // the geodetic reader does.  The natural spline's second derivative is exactly
    // 0 at the ends and recovers over seconds; the extended fit has a at the end.
    const double a = 1.5, v0 = 20.0;
    std::vector<double> t, y;
    for (std::size_t i = 0; i <= 6550; i++)
    {
        const double ti = 0.02 * static_cast<double>(i);
        t.push_back(ti);
        y.push_back(v0 * ti + 0.5 * a * ti * ti + jitter(i));
    }
    const double budget = 3.048;
    const auto natural = SmoothingSpline::within_tolerance(t, y, budget * budget * static_cast<double>(t.size()));
    const auto freed = SmoothingSpline::within_tolerance_extended(t, y, budget, 3.0, 10.0);

    REQUIRE(std::fabs(natural.derivative(t.back(), 2)) < 1e-9);
    for (double from_end : {0.0, 1.0, 2.0, 4.0, 8.0})
    {
        const double at = t.back() - from_end;
        INFO("seconds from the end: " << from_end);
        // Curvature at 0, 0.5, 1, 2, 4, 8 s from the end, for a = 1.5: natural 0, 0.01,
        // 0.05, 0.16, 0.50, 1.14; freed 1.33, 1.37, 1.40, 1.46, 1.53, 1.56.  The end
        // point itself keeps about a tenth of the bias, where the natural fit has all.
        const double tolerance = from_end == 0.0 ? 0.15 : from_end <= 1.0 ? 0.10 : 0.06;
        REQUIRE(freed.derivative(at, 2) == Catch::Approx(a).epsilon(tolerance));
        REQUIRE(freed.derivative(t.front() + from_end, 2) == Catch::Approx(a).epsilon(tolerance));
        // And the speed: the natural fit's end speed is off where it matters.
        REQUIRE(freed.derivative(at, 1) == Catch::Approx(v0 + a * at).epsilon(0.01));
    }
    // The natural fit really was biased here, so the test is not vacuous.
    REQUIRE(std::fabs(natural.derivative(t.back() - 2.0, 2) - a) > 0.15 * a);
}

TEST_CASE("freed ends carry a track that ends in a stop", "[spline][ends]")
{
    // 30 m/s decelerating at 1.5 m/s^2 to rest, then the track ends: the case the
    // reviewers measured.  The natural fit's last state was still doing 5.4 m/s and
    // its last four seconds had 43 percent too little deceleration.
    std::vector<double> t, y;
    const double a = -1.5, v0 = 30.0, t_stop = 20.0;
    for (std::size_t i = 0; i <= static_cast<std::size_t>(t_stop / 0.02); i++)
    {
        const double ti = 0.02 * static_cast<double>(i);
        t.push_back(ti);
        y.push_back(v0 * ti + 0.5 * a * ti * ti + jitter(i));
    }
    const double budget = 3.048;
    const auto natural = SmoothingSpline::within_tolerance(t, y, budget * budget * static_cast<double>(t.size()));
    const auto freed = SmoothingSpline::within_tolerance_extended(t, y, budget, 3.0, 10.0);

    const double end = t.back();
    REQUIRE(std::fabs(natural.derivative(end, 1)) > 3.0); // the bias being removed
    REQUIRE(std::fabs(freed.derivative(end, 1)) < 1.0);
    for (double from_end : {1.0, 2.0, 4.0})
    {
        INFO("seconds from the end: " << from_end);
        REQUIRE(freed.derivative(end - from_end, 2) == Catch::Approx(a).epsilon(0.15));
    }
}

TEST_CASE("freed ends keep the fit's budget and leave a line a line", "[spline][ends]")
{
    // The RMS budget keeps its meaning: the extended fit stays within it over the
    // real points too.  A straight line has no end bias to remove and comes back a
    // straight line (zero curvature everywhere, the end included).
    std::vector<double> t, y;
    for (std::size_t i = 0; i <= 3000; i++)
    {
        t.push_back(0.02 * static_cast<double>(i));
        y.push_back(50.0 + 12.0 * t.back() + jitter(i));
    }
    const auto line = SmoothingSpline::within_tolerance_extended(t, y, 3.048, 3.0, 10.0);
    for (double at : {0.0, 1.0, 30.0, 59.0, 60.0})
    {
        REQUIRE(std::fabs(line.derivative(at, 2)) < 2e-3);
        REQUIRE(line.derivative(at, 1) == Catch::Approx(12.0).margin(1e-2));
    }
    double squares = 0.0;
    for (std::size_t i = 0; i < t.size(); i++)
        squares += (line(t[i]) - y[i]) * (line(t[i]) - y[i]);
    REQUIRE(std::sqrt(squares / static_cast<double>(t.size())) <= 3.048);
}

TEST_CASE("an end with too little data in its window stays natural", "[spline][ends]")
{
    // Sparse samples (under 8 in the 3 s window) cannot support a quadratic, so the
    // fit is the plain one, bit for bit.
    std::vector<double> t, y;
    for (std::size_t i = 0; i < 40; i++)
    {
        t.push_back(2.0 * static_cast<double>(i));
        y.push_back(0.5 * t.back() * t.back() * 0.1 + jitter(i));
    }
    const double budget = 1.0;
    const auto plain = SmoothingSpline::within_tolerance(t, y, budget * budget * static_cast<double>(t.size()));
    const auto freed = SmoothingSpline::within_tolerance_extended(t, y, budget, 3.0, 10.0);
    for (double at : {0.0, 10.0, 40.0, 78.0})
    {
        REQUIRE(freed(at) == plain(at));
        REQUIRE(freed.derivative(at, 2) == plain.derivative(at, 2));
    }
}

TEST_CASE("a tolerance beyond the line's residual returns the exact weighted line", "[spline][spaps]")
{
    // Closed form, not a limit: the line the penalty tends to is the weighted
    // least-squares fit, and that is what comes back, curvature exactly zero.
    const auto x = linspace(0.0, 6.0, 40);
    const auto y = noisy_curve(x);
    std::vector<double> w(x.size(), 1.0);
    for (std::size_t i = 0; i < w.size(); i += 3)
        w[i] = 4.0;
    const auto spline = SmoothingSpline::within_tolerance(x, y, 1e9, w);
    // The normal equations, written out.
    double sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < x.size(); i++)
    {
        sw += w[i];
        sx += w[i] * x[i];
        sy += w[i] * y[i];
        sxx += w[i] * x[i] * x[i];
        sxy += w[i] * x[i] * y[i];
    }
    const double slope = (sw * sxy - sx * sy) / (sw * sxx - sx * sx);
    const double intercept = (sy - slope * sx) / sw;
    for (std::size_t i = 0; i < x.size(); i++)
    {
        REQUIRE(spline.fitted()[i] == Catch::Approx(intercept + slope * x[i]).margin(1e-12));
        REQUIRE(spline.derivative(x[i], 2) == 0.0);
    }
    // And it is what a very stiff penalty approaches.
    const auto stiff = SmoothingSpline::with_smoothing(x, y, 1e8, w);
    REQUIRE(max_abs_difference(stiff.fitted(), spline.fitted()) < 1e-6);
}

TEST_CASE("long, straight and level tracks are fitted, not refused", "[spline][spaps]")
{
    // Fifty-hertz level flight with a millimeter of noise and a budget of
    // meters.  The search used to climb lambda to 1e200 and, from about 30,000
    // samples, lose positive definiteness on the way and reject the input.
    for (const std::size_t n : {std::size_t(6551), std::size_t(30000), std::size_t(60000)})
    {
        std::vector<double> t(n), y(n);
        for (std::size_t i = 0; i < n; i++)
        {
            t[i] = 0.02 * static_cast<double>(i);
            y[i] = 65.0 + 1e-3 * jitter(i);
        }
        const double budget = 3.048 * 3.048 * static_cast<double>(n);
        const auto spline = SmoothingSpline::within_tolerance(t, y, budget);
        INFO("n = " << n);
        REQUIRE(spline.residual() <= budget);
        REQUIRE(std::fabs(spline.derivative(0.5 * t.back(), 2)) < 1e-12);
        REQUIRE(std::fabs(spline.derivative(0.5 * t.back(), 1)) < 1e-5);
    }
    // A single tiny time gap used to do the same at n = 6551.
    std::vector<double> t(6551), y(6551);
    for (std::size_t i = 0; i < t.size(); i++)
    {
        t[i] = 0.02 * static_cast<double>(i) + (i > 3000 ? 1e-6 : 0.0);
        y[i] = 65.0 + 1e-3 * jitter(i);
    }
    REQUIRE_NOTHROW(SmoothingSpline::within_tolerance(t, y, 3.048 * 3.048 * 6551.0));
}

TEST_CASE("a budget the factorization cannot reach on a long track is kept, not refused", "[spline][spaps]")
{
    // The residual approaches the line's slowly, so a tolerance just under it
    // asks for a lambda past what the n^4-conditioned system can solve at this
    // length (about 40,000 points).  The answer is then the last lambda that
    // could be solved: inside the tolerance, which is the constraint, if less
    // smooth than asked for.
    const std::size_t n = 50000;
    std::vector<double> t(n), y(n);
    for (std::size_t i = 0; i < n; i++)
    {
        t[i] = 0.02 * static_cast<double>(i);
        y[i] = 0.3 * std::sin(0.01 * t[i]) + 1e-3 * jitter(i);
    }
    const auto line = SmoothingSpline::within_tolerance(t, y, 1e30);
    const double hair = line.residual() * (1.0 - 1e-9);
    const auto spline = SmoothingSpline::within_tolerance(t, y, hair);
    REQUIRE(spline.residual() <= hair);
    REQUIRE(spline.residual() > 0.0);
    REQUIRE(std::isfinite(spline.derivative(500.0, 2)));
}

TEST_CASE("unevenly spaced knots give a continuous, natural spline", "[spline][nonuniform]")
{
    // Every other test uses linspace, where the spacing either side of a knot is the
    // same, so a transposition of the two in the assembly of R, Q or the penalty is
    // invisible.  Irregular timestamps are real input (fit_path fits the raw time
    // column), and there a wrong assembly gives a spline that is not C2 across knots.
    const std::vector<double> x = {0.0, 0.3, 1.0, 1.1, 2.5, 2.6, 4.0, 7.0, 7.2, 9.0};
    std::vector<double> y(x.size());
    for (std::size_t i = 0; i < x.size(); i++)
        y[i] = std::sin(1.3 * x[i]) + 0.2 * std::cos(5.0 * x[i]);
    for (const double lambda : {0.0, 1e-3, 1e-1})
    {
        const auto spline = SmoothingSpline::with_smoothing(x, y, lambda);
        INFO("lambda " << lambda);
        // Natural: no curvature at either end.
        CHECK(std::fabs(spline.derivative(x.front(), 2)) < 1e-9);
        CHECK(std::fabs(spline.derivative(x.back(), 2)) < 1e-9);
        // Value, slope and curvature continue across every interior knot.
        const double eps = 1e-7;
        for (std::size_t k = 1; k + 1 < x.size(); k++)
        {
            for (int order = 0; order <= 2; order++)
            {
                const double left = order == 0 ? spline(x[k] - eps) : spline.derivative(x[k] - eps, order);
                const double right = order == 0 ? spline(x[k] + eps) : spline.derivative(x[k] + eps, order);
                INFO("knot " << x[k] << " order " << order);
                CHECK(left == Catch::Approx(right).margin(1e-4));
            }
        }
        if (lambda == 0.0)
        {
            for (std::size_t k = 0; k < x.size(); k++)
                CHECK(spline(x[k]) == Catch::Approx(y[k]).margin(1e-9));
        }
    }
    // A line is reproduced exactly at any lambda, however the knots fall.
    std::vector<double> line(x.size());
    for (std::size_t i = 0; i < x.size(); i++)
        line[i] = 2.0 * x[i] + 1.0;
    for (const double lambda : {0.0, 1e-3, 1.0, 1e3})
    {
        const auto spline = SmoothingSpline::with_smoothing(x, line, lambda);
        for (double at : {0.0, 0.5, 1.05, 3.3, 7.1, 9.0})
        {
            INFO("lambda " << lambda << " at " << at);
            CHECK(spline(at) == Catch::Approx(2.0 * at + 1.0).margin(1e-8));
            CHECK(spline.derivative(at, 1) == Catch::Approx(2.0).margin(1e-8));
        }
    }
}
