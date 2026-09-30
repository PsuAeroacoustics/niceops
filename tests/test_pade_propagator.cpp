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
    PeSetup setup = refracting(250.0);
    const auto pade = solve_pe(setup);
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
    // A 100 m source over a 1.2 m receiver: 35 degrees is 141 m out.
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
        if (pade.range[q] < 141.0)
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
