// The split-step Pade march (pade_propagator.hpp, pe.hpp).
//
// The exact solutions in test_pe.cpp, test_pe_shadow.cpp and test_benchmark.cpp
// are what say the march is right, and they run through it by default.  These
// check the two things it rests on -- that the rational propagator is the
// exact one over its interval and amplifies nothing anywhere else, and that
// the result therefore does not depend on the range step -- and hold it
// against the Crank-Nicolson march it replaced, which shares the grid and the
// starter but none of the range arithmetic.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "exact_solutions.hpp"
#include "ground_effect.hpp"
#include "pade_propagator.hpp"
#include "pe.hpp"
#include "util/constants.hpp"

using namespace niceops;

namespace {

using complex = std::complex<double>;

// An upward-refracting column over grass: a shadow past about 500 m for a
// 10 m source.
PeSetup refracting(double f)
{
    PeSetup setup;
    setup.frequency = f;
    setup.source_height = 10.0;
    setup.receiver_heights = {1.2, 20.0};
    setup.max_range = 1200.0;
    setup.output_step = 20.0;
    setup.z = {0.0, 1000.0};
    setup.c = {340.0, 340.0 - 0.05 * 1000.0};
    GroundModel ground;
    ground.flow_resistance = 200.0;
    setup.admittance = ground_impedance(ground, {f}, 340.0)[0];
    return setup;
}

} // namespace

TEST_CASE("the Pade propagator is the exact one on its interval and amplifies nothing", "[pade]")
{
    for (const double wavelengths : {0.25, 1.0, 4.0, 8.0, 16.0})
    {
        PadeRequest request;
        request.sigma = 2.0 * pi * wavelengths;
        // 35 degrees, and a little above zero for a surface wave.
        request.accurate_low = -0.33;
        request.accurate_high = 0.05;
        request.spectrum_low = -42.0; // 20 points a wavelength
        request.spectrum_imag = 0.2; // the absorbing layer
        const PadePropagator r = pade_propagator(request);
        INFO("step of " << wavelengths << " wavelengths, " << r.poles.size() << " terms");
        CHECK(r.error <= request.tolerance);
        // Points the fit did not see.
        for (double x = request.accurate_low; x <= request.accurate_high; x += 0.00137)
            CHECK(std::abs(r(x) - exact_propagator(x, request.sigma)) < 2.0 * request.tolerance);
        // Nothing grows: along the real axis to the grid's last mode, and
        // above it as far as the absorbing layer's eigenvalues go.
        double largest = 0.0;
        for (const double imag : {0.0, 0.003, 0.03, 0.2})
            for (double x = request.spectrum_low; x <= request.accurate_high; x += 0.0031)
                largest = std::max(largest, std::abs(r(complex(x, imag))));
        CHECK(largest <= 1.0 + 1e-6);
        // Every solve is against Q - pole, and Q's spectrum is on and above
        // the axis.
        for (const complex pole : r.poles)
            CHECK(pole.imag() < -0.01);
        // About a term a wavelength, and a few to start with.
        CHECK(static_cast<double>(r.poles.size()) <= 5.0 + 0.6 * wavelengths);
    }
}

TEST_CASE("the Pade propagator absorbs what it cannot carry", "[pade]")
{
    // Interpolating a function of modulus one gives an all-pass r unless it is
    // pinned among the evanescent modes (pade_propagator.cpp): the modes
    // steeper than the interval then linger instead of leaving.  At an
    // eight-wavelength step the horizontal mode, x = -1, is down by 30 dB a
    // step and the evanescent ones by more.
    PadeRequest request;
    request.sigma = 2.0 * pi * 8.0;
    request.accurate_low = -0.33;
    request.accurate_high = 0.05;
    request.spectrum_low = -42.0;
    request.spectrum_imag = 0.2;
    const PadePropagator r = pade_propagator(request);
    CHECK(std::abs(r(-1.0)) < 0.05);
    for (const double x : {-1.5, -3.0, -10.0, -40.0})
        CHECK(std::abs(r(x)) < 0.02);
    // And it follows the exact propagator's branch: decaying, not growing,
    // for evanescent modes.
    CHECK(std::abs(exact_propagator(-2.0, 1.0)) == Catch::Approx(std::exp(-1.0)).epsilon(1e-12));
    CHECK(std::abs(exact_propagator(complex(-0.2, 0.1), 5.0)) < 1.0);
}

TEST_CASE("a step too long for its interval is refused, not approximated", "[pade]")
{
    PadeRequest request;
    request.sigma = 2.0 * pi * 300.0;
    request.accurate_low = -0.33;
    request.accurate_high = 0.05;
    request.spectrum_low = -42.0;
    request.spectrum_imag = 0.2;
    CHECK_THROWS_AS(pade_propagator(request), std::runtime_error);
}

TEST_CASE("the split-step march does not depend on its range step", "[pade]")
{
    // Through the shadow boundary and 50 dB into the shadow: a two-wavelength
    // step and a sixteen-wavelength one, which share no propagator.
    PeSetup setup = refracting(500.0);
    setup.range_step_wavelengths = 2.0;
    const auto fine = solve_pe(setup);
    setup.range_step_wavelengths = 16.0;
    const auto coarse = solve_pe(setup);
    CHECK(coarse.range_step > 4.0 * fine.range_step);
    CHECK(coarse.pade_terms > fine.pade_terms);
    REQUIRE(fine.range.size() == coarse.range.size());
    int compared = 0;
    double deepest = 0.0;
    for (std::size_t k = 0; k < setup.receiver_heights.size(); k++)
    {
        for (std::size_t q = 0; q < fine.range.size(); q++)
        {
            // Past the first few steps.  Within them the modes steeper than
            // max_angle are still on their way out, and how fast they go is
            // the one thing that does depend on the step: 0.1-0.4 dB at 20
            // and 40 m here.
            if (std::isnan(fine.relative_level[k][q]) || fine.range[q] < 100.0)
                continue;
            INFO("r = " << fine.range[q] << ", z_r = " << setup.receiver_heights[k]);
            CHECK(coarse.relative_level[k][q] == Catch::Approx(fine.relative_level[k][q]).margin(0.005));
            deepest = std::min(deepest, fine.relative_level[k][q]);
            compared++;
        }
    }
    CHECK(compared > 100);
    CHECK(deepest < -40.0);
}

TEST_CASE("the split-step march is Crank-Nicolson's with a fine enough step", "[pade]")
{
    // Low angles, where Crank-Nicolson's Pade(1,1) equation is good and only
    // its range step stands between the two: forty steps a wavelength.
    // Like for like: Crank-Nicolson's own second-order grid and Gaussian
    // starter, so that the range arithmetic is all that differs.
    PeSetup setup = refracting(250.0);
    setup.vertical_order = 2;
    setup.starter = PeStarter::gaussian;
    const auto pade = solve_pe(setup);
    CHECK(!pade.self_started);
    setup.scheme = PeScheme::crank_nicolson;
    setup.range_points_per_wavelength = 40.0;
    const auto reference = solve_pe(setup);
    REQUIRE(pade.range.size() == reference.range.size());
    int compared = 0;
    for (std::size_t k = 0; k < setup.receiver_heights.size(); k++)
    {
        for (std::size_t q = 0; q < pade.range.size(); q++)
        {
            if (pade.range[q] < 100.0 || reference.relative_level[k][q] < -40.0)
                continue;
            INFO("r = " << pade.range[q] << ", z_r = " << setup.receiver_heights[k]);
            CHECK(pade.relative_level[k][q] == Catch::Approx(reference.relative_level[k][q]).margin(0.05));
            compared++;
        }
    }
    CHECK(compared > 80);
}

TEST_CASE("the split-step march reports no level steeper than it reaches", "[pade]")
{
    // A 100 m source over a 1.2 m receiver: 50 degrees is 83 m out.
    PeSetup setup;
    setup.frequency = 250.0;
    setup.source_height = 100.0;
    setup.receiver_heights = {1.2};
    setup.max_range = 300.0;
    setup.output_step = 20.0;
    setup.z = {0.0};
    setup.c = {340.0};
    const auto pade = solve_pe(setup);
    for (std::size_t q = 0; q < pade.range.size(); q++)
    {
        INFO("r = " << pade.range[q]);
        if (pade.range[q] < 82.9)
            CHECK(std::isnan(pade.relative_level[0][q]));
        else
            CHECK(std::abs(pade.relative_level[0][q]) < 10.0);
    }
    // Crank-Nicolson, which carries those modes, reports everything.
    setup.scheme = PeScheme::crank_nicolson;
    const auto cn = solve_pe(setup);
    for (const double level : cn.relative_level[0])
        CHECK(std::isfinite(level));
}

TEST_CASE("the self-started, fourth-order march is the exact field out to its angle", "[pade]")
{
    // A 100 m source over soft ground in still air at 1 kHz, from directly
    // below the angle limit outward: the direct and reflected rays with the
    // spherical-wave reflection coefficient.  A Gaussian starter on a
    // second-order grid is 0.5-1.2 dB off between 30 and 50 degrees here; the
    // two changes together bring that under 0.1.
    const double c = 340.0, f = 1000.0, z_s = 100.0, z_r = 1.0;
    GroundModel ground;
    ground.flow_resistance = 200.0;
    const auto admittance = ground_impedance(ground, {f}, c)[0];
    PeSetup setup;
    setup.frequency = f;
    setup.source_height = z_s;
    setup.receiver_heights = {z_r};
    setup.max_range = 600.0;
    setup.output_step = 5.0;
    setup.z = {0.0};
    setup.c = {c};
    setup.admittance = admittance;
    auto exact = [&](double r) {
        const double r1 = std::hypot(r, z_s - z_r), r2 = std::hypot(r, z_s + z_r);
        const auto q = spherical_reflection_coefficient((z_s + z_r) / r2, pi * f * r2 / c, admittance, true);
        const complex field = 1.0 + complex(q[0], q[1]) * (r1 / r2) * std::exp(complex(0.0, 2.0 * pi * f / c * (r2 - r1)));
        return 20.0 * std::log10(std::abs(field));
    };
    // rms over the ranges whose direct path lies between the two elevations,
    // leaving out the interference minima.
    auto rms = [&](const PeResult& pe, double steepest, double shallowest) {
        double sum = 0.0;
        int count = 0;
        for (std::size_t q = 0; q < pe.range.size(); q++)
        {
            const double angle = std::atan2(z_s - z_r, pe.range[q]) * 180.0 / pi;
            const double level = exact(pe.range[q]);
            if (angle > steepest || angle < shallowest || level < -10.0)
                continue;
            const double d = pe.relative_level[0][q] - level;
            sum += d * d;
            count++;
        }
        REQUIRE(count > 8);
        return std::sqrt(sum / count);
    };
    const auto wide = solve_pe(setup);
    CHECK(wide.self_started);
    const double wide_steep = rms(wide, 50.0, 30.0), wide_shallow = rms(wide, 30.0, 9.0);
    // What it replaces, for the size of the difference.
    setup.starter = PeStarter::gaussian;
    setup.vertical_order = 2;
    const auto narrow = solve_pe(setup);
    const double narrow_steep = rms(narrow, 50.0, 30.0);
    INFO("30-50 degrees: " << wide_steep << " dB rms, Gaussian on second order " << narrow_steep
                           << "; 9-30 degrees: " << wide_shallow);
    CHECK(wide_steep < 0.15);
    CHECK(wide_shallow < 0.05);
    CHECK(narrow_steep > 3.0 * wide_steep);
}

TEST_CASE("fourth order in height converges from fewer points", "[pade]")
{
    // The refracting case at 10 and at 20 points a wavelength: fourth order
    // moves less between them than second order does.
    auto spread = [](int order) {
        PeSetup setup = refracting(500.0);
        setup.vertical_order = order;
        setup.points_per_wavelength = 10.0;
        const auto coarse = solve_pe(setup);
        setup.points_per_wavelength = 20.0;
        const auto fine = solve_pe(setup);
        double sum = 0.0;
        int count = 0;
        for (std::size_t q = 0; q < fine.range.size(); q++)
        {
            if (fine.range[q] < 100.0 || std::isnan(fine.relative_level[1][q]))
                continue;
            const double d = coarse.relative_level[1][q] - fine.relative_level[1][q];
            sum += d * d;
            count++;
        }
        return std::sqrt(sum / count);
    };
    const double second = spread(2), fourth = spread(4);
    INFO("second order " << second << " dB rms, fourth order " << fourth);
    CHECK(fourth < 0.5 * second);
    CHECK(fourth < 0.1);
}

TEST_CASE("the self-starter gives an elevated source its own strength", "[pade]")
{
    // The benchmark's upward-refracting column (test_benchmark.cpp) with the
    // source 100 m up, where k is 3% above its value at the ground, against
    // the residue series.  The Gaussian starter is normalised with the
    // ground's wavenumber and comes out 0.11-0.15 dB high at every range; the
    // self-starter launches each mode with its own k_n^(-1/2).
    const double c0 = 343.23, radius = c0 / 0.1, f = 100.0;
    ExactProblem problem;
    problem.frequency = f;
    problem.sound_speed = c0;
    problem.admittance = 1.0 / complex(12.81, 11.62);
    problem.source_height = 100.0;
    problem.receiver_height = 1.0;
    PeSetup setup;
    setup.frequency = f;
    setup.source_height = 100.0;
    setup.receiver_heights = {1.0};
    setup.max_range = 2000.0;
    setup.output_step = 50.0;
    setup.admittance = problem.admittance;
    for (double z = 0.0; z <= 1600.0; z += 0.1)
    {
        setup.z.push_back(z);
        setup.c.push_back(c0 / std::sqrt(1.0 + 2.0 * z / radius));
    }
    const auto self = solve_pe(setup);
    setup.starter = PeStarter::gaussian;
    const auto gaussian = solve_pe(setup);
    const auto series = upward_pressure(problem, radius, self.range);
    int compared = 0;
    for (std::size_t q = 0; q < self.range.size(); q++)
    {
        const double r = self.range[q];
        if (r < 400.0 || std::isnan(series[q].real()))
            continue;
        const double exact = 20.0 * std::log10(std::abs(series[q]) * std::hypot(r, 99.0));
        INFO("r = " << r << ", exact " << exact);
        CHECK(self.relative_level[0][q] == Catch::Approx(exact).margin(0.05));
        CHECK(gaussian.relative_level[0][q] - exact > 0.08);
        compared++;
    }
    CHECK(compared > 20);
}

TEST_CASE("with the ground and the interpolation at fourth order the grid converges at fourth order", "[pade]")
{
    // The refracting case over grass, receivers 1.2 and 20 m up, at 6 and at
    // 12 points a wavelength against 24.  A fourth-order scheme gains a
    // factor of 16 from the doubling, a second-order one 4 -- and with the
    // operator at fourth order but the ground condition and the
    // interpolation at second, that is all it gained.
    auto error = [](double points, int boundary) {
        PeSetup setup = refracting(500.0);
        setup.boundary_order = 4;
        setup.points_per_wavelength = 24.0;
        const auto reference = solve_pe(setup);
        setup.boundary_order = boundary;
        setup.points_per_wavelength = points;
        const auto pe = solve_pe(setup);
        double sum = 0.0;
        int count = 0;
        for (std::size_t k = 0; k < setup.receiver_heights.size(); k++)
        {
            for (std::size_t q = 0; q < pe.range.size(); q++)
            {
                if (pe.range[q] < 100.0 || std::isnan(reference.relative_level[k][q]))
                    continue;
                const double d = pe.relative_level[k][q] - reference.relative_level[k][q];
                sum += d * d;
                count++;
            }
        }
        return std::sqrt(sum / count);
    };
    const double fourth_6 = error(6.0, 4), fourth_12 = error(12.0, 4);
    const double second_6 = error(6.0, 2), second_12 = error(12.0, 2);
    INFO("fourth-order boundaries: " << fourth_6 << " dB rms at 6 points, " << fourth_12 << " at 12; second: "
                                     << second_6 << ", " << second_12);
    CHECK(fourth_6 / fourth_12 > 10.0);
    CHECK(second_6 / second_12 < 6.0);
    CHECK(fourth_12 < 0.005);
    CHECK(second_12 > 10.0 * fourth_12);
}

TEST_CASE("the grid defaults to ten points a wavelength only where it is fourth order throughout", "[pade]")
{
    PeSetup setup = refracting(250.0);
    setup.max_range = 400.0;
    const auto chosen = solve_pe(setup);
    setup.points_per_wavelength = 10.0;
    const auto ten = solve_pe(setup);
    CHECK(chosen.relative_level[0] == ten.relative_level[0]);
    setup.points_per_wavelength = 0.0;
    setup.boundary_order = 2;
    const auto second = solve_pe(setup);
    setup.points_per_wavelength = 20.0;
    const auto twenty = solve_pe(setup);
    CHECK(second.relative_level[0] == twenty.relative_level[0]);
}

TEST_CASE("a long first step is still self-started", "[pade]")
{
    // The starter's fit is refused at sixteen wavelengths on this grid; the
    // first step is then taken in parts rather than handed to the Gaussian,
    // which would show as a step-dependent offset of 0.007 dB here.
    PeSetup setup = refracting(500.0);
    setup.range_step_wavelengths = 16.0;
    CHECK(solve_pe(setup).self_started);
}
