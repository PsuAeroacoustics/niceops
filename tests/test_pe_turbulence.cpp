// The PE's turbulence (pe.hpp): its phase screens, against the von Karman
// spectrum they are drawn from, and the ensemble-mean field, against the
// exponential decay the Markov approximation gives it.
//
// Each screen is the phase a slab `dx` deep accumulates, so its variance is
// 2 pi k^2 dx times the 3-D spectrum integrated over the plane normal to the
// path: 0.782 k^2 dx C^2 L0^(5/3) with C^2 = gamma_t / 4.  After r metres of
// screens the phase along any path has variance 0.782 k^2 r C^2 L0^(5/3), and
// with it Gaussian the mean field is the still-air field times exp(-sigma^2/2),
// both rays of it alike.  The mean field's decay therefore checks the screens'
// amplitude and their application in the march together.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include "pe.hpp"
#include "refraction.hpp"
#include "util/constants.hpp"

using namespace niceops;

namespace {

double screen_variance(double k, double gamma_t, double outer_scale, double dx)
{
    return 0.7817 * k * k * dx * 0.25 * gamma_t * std::pow(outer_scale, 5.0 / 3.0);
}

} // namespace

TEST_CASE("pe_phase_screen_variance", "[pe_turbulence]")
{
    const double k = 2.0 * pi * 500.0 / 340.0;
    const double gamma_t = 1e-5, outer_scale = 10.0, dx = 10.0, dz = 0.068;
    const std::size_t count = 8000; // a 544 m column
    double sum = 0.0, sum2 = 0.0;
    std::size_t n = 0;
    for (unsigned seed = 1; seed <= 200; seed++)
    {
        const auto phi = phase_screen(k, gamma_t, outer_scale, dx, dz, count, seed);
        for (const double v : phi)
        {
            sum += v;
            sum2 += v * v;
            n++;
        }
    }
    const double mean = sum / n;
    const double variance = sum2 / n - mean * mean;
    // The finite column leaves out kappa below 2 pi / (2 column): about 4% of
    // the variance at this outer scale.
    CHECK(variance == Catch::Approx(0.96 * screen_variance(k, gamma_t, outer_scale, dx)).epsilon(0.08));
    CHECK(std::fabs(mean) < 0.05 * std::sqrt(variance));
}

TEST_CASE("pe_turbulent_mean_field_decays", "[pe_turbulence]")
{
    // 100 Hz over rigid ground, source and receivers 40 m up; gamma_t chosen
    // so that sigma^2 = 1 at 300 m.
    PeSetup setup;
    setup.frequency = 100.0;
    setup.source_height = 40.0;
    setup.receiver_heights = {40.0};
    setup.max_range = 300.0;
    setup.output_step = 50.0;
    setup.z = {0.0};
    setup.c = {340.0};
    const double k = 2.0 * pi * setup.frequency / 340.0;
    setup.outer_scale = 10.0;
    setup.gamma_t = 1.0 / (0.7817 * k * k * 300.0 * 0.25 * std::pow(setup.outer_scale, 5.0 / 3.0));
    const auto still = [&] {
        PeSetup s = setup;
        s.gamma_t = 0.0;
        return solve_pe(s);
    }();
    constexpr int realisations = 300;
    std::vector<std::complex<double>> mean(still.range.size(), 0.0);
    for (int i = 0; i < realisations; i++)
    {
        setup.seed = 7 + static_cast<unsigned>(i);
        const auto one = solve_pe(setup);
        for (std::size_t q = 0; q < mean.size(); q++)
            mean[q] += one.psi[0][q] / static_cast<double>(realisations);
    }
    for (std::size_t q = 0; q < mean.size(); q++)
    {
        const double r = still.range[q];
        const double expected = std::exp(-0.5 * screen_variance(k, setup.gamma_t, setup.outer_scale, r));
        const double ratio = std::abs(mean[q] / still.psi[0][q]);
        INFO("range " << r << " m: |<psi>| / |psi_0| = " << ratio << ", expected " << expected);
        // Sampling error of the mean of 300 unit phasors is about 0.04.
        CHECK(std::fabs(ratio - expected) < 0.08);
    }
}

TEST_CASE("deep in a shadow the turbulent PE is at the scattering model's level", "[pe_turbulence]")
{
    // A measured day-284 upwind profile (the one tests/test_shadow_transition.cpp
    // uses), source 50 m up, microphone 3 mm up over the site's ground, 500 Hz,
    // gamma_t 3e-5: past a kilometre the still field is 20+ dB below the
    // turbulent one, which ShadowScattering was fitted to on twelve such
    // profiles.  Twelve realisations, spread across the band.
    PeSetup setup;
    setup.source_height = 50.0;
    setup.receiver_heights = {0.003};
    setup.max_range = 1600.0;
    setup.output_step = 100.0;
    setup.z = {0.0, 0.3, 0.91, 3.05, 9.14, 30.5, 60.0, 122.0, 300.0};
    setup.c = {339.19, 336.44, 335.10, 333.56, 331.97, 330.77, 330.30, 330.53, 331.20};
    // Variable porosity, sigma_e 200 kPa s/m^2, at 500 Hz: beta = 1/Z,
    // Z = (1 + i) sqrt(R_se / (pi gamma rho0 f)).
    const double z_mag = std::sqrt(200e3 / (pi * 1.4 * 1.2 * 500.0));
    setup.admittance = 1.0 / std::complex<double>(z_mag, z_mag);
    setup.gamma_t = 3e-5;
    setup.outer_scale = 10.0;
    const auto band = solve_pe_band(setup, 500.0, 12);
    const auto model = ShadowScattering::for_receiver(setup.gamma_t, 0.003, false);
    std::vector<double> difference;
    for (std::size_t q = 0; q < band.range.size(); q++)
    {
        const double r = band.range[q];
        if (r < 1000.0)
            continue;
        const double slant = std::hypot(r, setup.source_height);
        difference.push_back(band.relative_level[0][q] - model.level(500.0, slant, setup.source_height));
        INFO(r << " m: PE " << band.relative_level[0][q] << " dB, model " << model.level(500.0, slant, 50.0));
    }
    REQUIRE(difference.size() >= 5);
    std::sort(difference.begin(), difference.end());
    const double median = difference[difference.size() / 2];
    INFO("median PE - model " << median << " dB");
    CHECK(std::fabs(median) < 3.0);
}
