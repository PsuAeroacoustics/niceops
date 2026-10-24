// PE propagation tables (pe_table.hpp): reciprocity, and turbulence in them.
//
// A table can be filled by one march per source height (the receivers read
// from the field) or by one march per receiver height with the source put
// there (the sources read from the field).  In the effective-sound-speed model
// these are the same field, so the two must agree; on the grid that holds only
// if a source near the ground is the adjoint of the receivers' read, which is
// what pe.cpp's place_ground_source is for (tests/test_pe.cpp holds it against
// exact fields).

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "atmosphere_profile.hpp"
#include "ground_effect.hpp"
#include "pe.hpp"
#include "pe_table.hpp"
#include "stratified.hpp"
#include "util/constants.hpp"

using namespace niceops;

namespace {

// A small table through the windy example: downwind and upwind, a flush and a
// 1.2 m microphone, two low sources, 125 Hz and 500 Hz over grass.
PeTableRequest windy_request()
{
    auto profile = read_atmosphere_profile("atmospheres/example_windy.csv", feet, 0.0);
    if (!profile.has_pressure())
        profile = profile.with_hydrostatic_pressure(reference_pressure);
    const LayeredSoundSpeed layered(profile, 90.0);
    PeTableRequest request;
    request.bearings = {90.0, 270.0};
    request.z = layered.heights();
    for (const double bearing : request.bearings)
    {
        std::vector<double> c(layered.levels());
        layered.along(std::sin(bearing * pi / 180.0), std::cos(bearing * pi / 180.0), c.data());
        request.c_by_bearing.push_back(c);
    }
    request.source_heights = {3.0, 15.0};
    request.receiver_heights = {0.0, 1.2};
    request.max_range = 600.0;
    request.step = 50.0;
    request.bands = {125.0, 500.0};
    GroundModel ground;
    ground.flow_resistance = 200.0;
    request.admittance = ground_impedance(ground, request.bands, profile.at(0.0).sound_speed());
    request.samples = 2;
    return request;
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b)
{
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

} // namespace

TEST_CASE("a table marched from the receivers is the table marched from the sources", "[pe_table]")
{
    auto request = windy_request();
    request.march = PeTableMarch::source;
    const auto by_source = build_pe_table(request);
    request.march = PeTableMarch::receiver;
    const auto by_receiver = build_pe_table(request);
    REQUIRE(by_source.size() == by_receiver.size());
    double worst = 0.0;
    std::size_t compared = 0;
    for (std::size_t k = 0; k < by_source.size(); k++)
    {
        // NaN (steeper than max_angle) in the same cells either way.
        REQUIRE(std::isnan(by_source[k]) == std::isnan(by_receiver[k]));
        if (std::isnan(by_source[k]))
            continue;
        worst = std::max(worst, std::fabs(by_source[k] - by_receiver[k]));
        compared++;
    }
    INFO("largest difference " << worst << " dB over " << compared << " cells");
    CHECK(compared > 150);
    // Measured 0.003 dB; the remainder is the domain's top, which follows the
    // highest height in each march.
    CHECK(worst < 0.02);
}

TEST_CASE("a still table does not depend on the turbulence options", "[pe_table]")
{
    auto request = windy_request();
    request.bands = {250.0};
    request.admittance.resize(1);
    request.samples = 1;
    const auto still = build_pe_table(request);
    request.realizations = 3;
    request.seed = 17;
    request.outer_scale = 3.0;
    CHECK(same_bits(still, build_pe_table(request)));
    // A profile with no reference strength is refused rather than ignored.
    request.turbulence_z = {0.0, 100.0};
    request.turbulence_gamma = {1e-5, 1e-6};
    CHECK_THROWS(build_pe_table(request));
}

TEST_CASE("a screen scaled by a strength profile has each height's variance", "[pe_turbulence]")
{
    // Screens drawn at the reference gamma_t, each row scaled by
    // sqrt(gamma(z) / gamma_t): the variance follows gamma(z), here falling
    // tenfold from the ground to 50 m and held above.
    const double k = 2.0 * pi * 500.0 / 340.0;
    const double gamma_ref = 1e-5, outer_scale = 10.0, dx = 10.0, dz = 0.068;
    const std::size_t count = 8000;
    PeSetup setup;
    setup.gamma_t = gamma_ref;
    setup.turbulence_z = {0.0, 50.0};
    setup.turbulence_gamma = {1e-5, 1e-6};
    const auto scale = turbulence_row_scale(setup, dz, count);
    CHECK(scale[0] == Catch::Approx(1.0));
    CHECK(scale[count - 1] == Catch::Approx(std::sqrt(0.1)));
    const auto expected = [&](double gamma) {
        return 0.96 * 0.7817 * k * k * dx * 0.25 * gamma * std::pow(outer_scale, 5.0 / 3.0);
    };
    for (const double z : {5.0, 25.0, 45.0, 200.0, 400.0})
    {
        const std::size_t j = static_cast<std::size_t>(z / dz);
        double sum2 = 0.0;
        int n = 0;
        for (unsigned seed = 1; seed <= 400; seed++)
        {
            const auto phi = phase_screen(k, gamma_ref, outer_scale, dx, dz, count, seed);
            const double v = scale[j] * phi[j];
            sum2 += v * v;
            n++;
        }
        const double gamma = z >= 50.0 ? 1e-6 : 1e-5 + (1e-6 - 1e-5) * (static_cast<double>(j) * dz) / 50.0;
        INFO("z = " << z << " m");
        // 400 draws: a variance's sampling error is sqrt(2/400) = 7%.
        CHECK(sum2 / n == Catch::Approx(expected(gamma)).epsilon(0.25));
    }
}

TEST_CASE("a strength profile scales the march's screens by its square root", "[extended][pe_turbulence]")
{
    // The mean-field decay of pe_turbulent_mean_field_decays with a uniform
    // profile at a quarter of the reference: the field must decay as gamma_t / 4
    // would.  Scaling the phase by the strength ratio instead of its square
    // root decays as gamma_t / 16, scaling it not at all as gamma_t.
    PeSetup setup;
    setup.frequency = 100.0;
    setup.source_height = 40.0;
    setup.receiver_heights = {40.0};
    setup.max_range = 300.0;
    setup.output_step = 100.0;
    setup.z = {0.0};
    setup.c = {340.0};
    setup.points_per_wavelength = 10.0;
    const double k = 2.0 * pi * setup.frequency / 340.0;
    setup.outer_scale = 10.0;
    // sigma^2 = 1 at 300 m for gamma = gamma_t / 4.
    const double gamma = 1.0 / (0.7817 * k * k * 300.0 * 0.25 * std::pow(setup.outer_scale, 5.0 / 3.0));
    setup.gamma_t = 4.0 * gamma;
    setup.turbulence_z = {0.0};
    setup.turbulence_gamma = {gamma};
    const auto still = [&] {
        PeSetup s = setup;
        s.turbulence_gamma = {1e-30};
        return solve_pe(s);
    }();
    constexpr int realizations = 600;
    std::vector<std::complex<double>> mean(still.range.size(), 0.0);
    for (int i = 0; i < realizations; i++)
    {
        setup.seed = 11 + static_cast<unsigned>(i);
        const auto one = solve_pe(setup);
        for (std::size_t q = 0; q < mean.size(); q++)
            mean[q] += one.psi[0][q] / static_cast<double>(realizations);
    }
    for (std::size_t q = 0; q < mean.size(); q++)
    {
        const double r = still.range[q];
        const double expected = std::exp(-0.5 * r / 300.0);
        const double ratio = std::abs(mean[q] / still.psi[0][q]);
        INFO("range " << r << " m: |<psi>| / |psi_0| = " << ratio << ", expected " << expected);
        // sqrt((1 - |mean|^2) / 2N) = 0.025 at N = 600: 0.09 is 3.5 sigma.
        CHECK(std::fabs(ratio - expected) < 0.09);
    }
}

TEST_CASE("deep in a shadow the scattered level is linear in gamma_t", "[extended][pe_turbulence]")
{
    // Single scattering puts energy into a shadow in proportion to gamma_t:
    // three times the strength is 4.77 dB more, where the still field is far
    // below the scattered one.  The day-284 upwind profile of the scattering
    // test in test_pe_turbulence.cpp, 500 Hz, four frequencies times three
    // realizations each, through solve_pe_band's realizations.
    PeSetup setup;
    setup.source_height = 50.0;
    setup.receiver_heights = {0.003};
    setup.max_range = 1600.0;
    setup.output_step = 100.0;
    setup.z = {0.0, 0.3, 0.91, 3.05, 9.14, 30.5, 60.0, 122.0, 300.0};
    setup.c = {339.19, 336.44, 335.10, 333.56, 331.97, 330.77, 330.30, 330.53, 331.20};
    const double z_mag = std::sqrt(200e3 / (pi * 1.4 * 1.2 * 500.0));
    setup.admittance = 1.0 / std::complex<double>(z_mag, z_mag);
    setup.outer_scale = 10.0;
    setup.gamma_t = 1e-5;
    const auto weak = solve_pe_band(setup, 500.0, 4, 3);
    setup.gamma_t = 3e-5;
    const auto strong = solve_pe_band(setup, 500.0, 4, 3);
    std::vector<double> difference;
    for (std::size_t q = 0; q < weak.range.size(); q++)
        if (weak.range[q] >= 1000.0)
        {
            difference.push_back(strong.relative_level[0][q] - weak.relative_level[0][q]);
            INFO(weak.range[q] << " m: " << difference.back() << " dB");
        }
    REQUIRE(difference.size() >= 5);
    std::sort(difference.begin(), difference.end());
    const double median = difference[difference.size() / 2];
    INFO("median strong - weak " << median << " dB, theory 4.77");
    CHECK(std::fabs(median - 4.77) < 1.0);
}

TEST_CASE("realizations average in energy and leave one realization as it was", "[pe_turbulence]")
{
    PeSetup setup;
    setup.source_height = 20.0;
    setup.receiver_heights = {1.2};
    setup.max_range = 400.0;
    setup.output_step = 100.0;
    setup.z = {0.0};
    setup.c = {340.0};
    setup.gamma_t = 1e-5;
    setup.seed = 5;
    const auto one = solve_pe_band(setup, 250.0, 2, 1);
    const auto two = solve_pe_band(setup, 250.0, 2, 2);
    // The second realization's seeds are seed + samples + sample: the band of
    // two realizations is the energy mean of seed 5's band and seed 7's.
    PeSetup other = setup;
    other.seed = 7;
    const auto second = solve_pe_band(other, 250.0, 2, 1);
    for (std::size_t q = 0; q < one.range.size(); q++)
    {
        const double mean = 10.0 * std::log10(0.5 * (std::pow(10.0, 0.1 * one.relative_level[0][q])
                                                        + std::pow(10.0, 0.1 * second.relative_level[0][q])));
        CHECK(two.relative_level[0][q] == Catch::Approx(mean).margin(1e-9));
    }
}
