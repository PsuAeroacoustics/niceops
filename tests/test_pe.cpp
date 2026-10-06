// The wide-angle parabolic equation (pe.hpp), against exact solutions.
//
// Over rigid ground in still air the field is two rays, direct and image.
// Over impedance ground it is the Weyl-van der Pol solution, whose spherical-
// wave reflection coefficient ground_effect.cpp already computes (and pins
// against panam).  The PE shares neither's derivation -- it marches a wave
// equation on a grid, with the impedance condition in the grid -- so agreement
// checks its normalization, its ground condition's sign convention and its
// absorbing top together.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <complex>
#include <vector>
#include <stdexcept>

#include "ground_effect.hpp"
#include "pe.hpp"
#include "util/constants.hpp"

using namespace niceops;

namespace {

// 20 log10 |1 + Q (R1/R2) exp(i k (R2 - R1))|: the level relative to the
// direct ray's free field, for reflection coefficient Q.
double two_ray(double r, double z_s, double z_r, double f, double c, std::complex<double> admittance,
        bool rigid)
{
    const double r1 = std::hypot(r, z_r - z_s);
    const double r2 = std::hypot(r, z_r + z_s);
    const double k = 2.0 * pi * f / c;
    std::complex<double> q(1.0, 0.0);
    if (!rigid)
    {
        const auto coefficient = spherical_reflection_coefficient((z_s + z_r) / r2, pi * f * r2 / c,
                admittance, true);
        q = std::complex<double>(coefficient[0], coefficient[1]);
    }
    const auto field = 1.0 + q * (r1 / r2) * std::exp(std::complex<double>(0.0, k * (r2 - r1)));
    return 20.0 * std::log10(std::abs(field));
}

} // namespace

TEST_CASE("the PE is two rays over rigid ground in still air", "[pe]")
{
    PeSetup setup;
    setup.frequency = 250.0;
    setup.source_height = 5.0;
    setup.receiver_heights = {1.2, 10.0};
    setup.max_range = 800.0;
    setup.output_step = 10.0;
    setup.z = {0.0};
    setup.c = {340.0};
    const auto result = solve_pe(setup);
    int compared = 0;
    for (std::size_t q = 0; q < result.range.size(); q++)
    {
        const double r = result.range[q];
        // Past 100 m: nearer in, the image ray's angle heads toward the
        // wide-angle equation's limit.
        if (r < 100.0)
            continue;
        for (std::size_t k = 0; k < setup.receiver_heights.size(); k++)
        {
            const double exact = two_ray(r, 5.0, setup.receiver_heights[k], 250.0, 340.0, 0.0, true);
            if (exact < -10.0)
                continue; // the deep nulls are where tiny phase errors show; skip them
            INFO("r = " << r << ", z_r = " << setup.receiver_heights[k]);
            CHECK(result.relative_level[k][q] == Catch::Approx(exact).margin(0.3));
            compared++;
        }
    }
    CHECK(compared > 100);
}

TEST_CASE("the PE is Weyl-van der Pol over impedance ground in still air", "[pe]")
{
    for (const double f : {125.0, 500.0})
    {
        GroundModel ground;
        ground.flow_resistance = 200.0;
        const auto admittance = ground_impedance(ground, {f}, 340.0)[0];
        PeSetup setup;
        setup.frequency = f;
        setup.source_height = 10.0;
        setup.receiver_heights = {1.2};
        setup.max_range = 1000.0;
        setup.z = {0.0};
        setup.c = {340.0};
        setup.admittance = admittance;
        const auto result = solve_pe(setup);
        int compared = 0;
        for (std::size_t q = 0; q < result.range.size(); q++)
        {
            const double r = result.range[q];
            if (r < 100.0)
                continue;
            const double exact = two_ray(r, 10.0, 1.2, f, 340.0, admittance, false);
            if (exact < -12.0)
                continue;
            INFO("f = " << f << ", r = " << r);
            CHECK(result.relative_level[0][q] == Catch::Approx(exact).margin(0.5));
            compared++;
        }
        CHECK(compared > 30);
    }
}

TEST_CASE("a source near or on the ground is where it is, as a receiver there is", "[pe]")
{
    // Within three grid points of the ground the self-starter's delta function
    // is the receivers' read turned around through the ground rows' own
    // quadrature (pe.cpp, place_ground_source).  Taking those rows' weights as
    // one put a source at node 1 1.34 dB high and one at node 2 0.37 dB low,
    // and a source at z = 0 -- a reciprocal table's flush microphone -- had no
    // self-starter at all.  Against the exact field, at grazing (200-1000 m)
    // and at 17-45 degrees (60-200 m) from a 60 m receiver, over rigid ground
    // and grass, the source agrees with the receiver read at the same height
    // (reciprocity) within 0.003 dB at grazing and 0.1 at the steep angles,
    // and with the exact field (Weyl-van der Pol over grass) within 0.02 and
    // 0.1.  The steep worst, 0.06-0.08 dB, is a source at 2-2.5 dz seen at
    // 45 degrees: node 2's weight, a second-order series in a bounded
    // variable (pe.cpp).  At node 0, the flush microphone, it is 0.03.
    for (const bool rigid : {true, false})
        for (const double f : {125.0, 500.0})
        {
            GroundModel ground;
            ground.flow_resistance = 200.0;
            const auto admittance = rigid ? std::complex<double>(0.0) : ground_impedance(ground, {f}, 340.0)[0];
            const double dz = 340.0 / f / 10.0; // the default grid
            for (const double u : {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.5})
                for (const bool steep : {false, true})
                {
                    const double low = u * dz, high = steep ? 60.0 : 10.0;
                    PeSetup setup;
                    setup.frequency = f;
                    setup.source_height = low;
                    setup.receiver_heights = {high};
                    setup.max_range = steep ? 200.0 : 1000.0;
                    setup.output_step = 10.0;
                    setup.z = {0.0};
                    setup.c = {340.0};
                    setup.admittance = admittance;
                    const auto result = solve_pe(setup);
                    CHECK(result.self_started);
                    PeSetup swapped = setup;
                    swapped.source_height = high;
                    swapped.receiver_heights = {low};
                    const auto reciprocal = solve_pe(swapped);
                    int compared = 0;
                    for (std::size_t q = 0; q < result.range.size(); q++)
                    {
                        const double r = result.range[q];
                        if (r < (steep ? 60.0 : 200.0) || std::isnan(result.relative_level[0][q]))
                            continue;
                        const double exact = two_ray(r, low, high, f, 340.0, admittance, rigid);
                        if (exact < -10.0)
                            continue;
                        INFO((rigid ? "rigid" : "grass") << ", " << f << " Hz, source " << u << " dz, r " << r);
                        CHECK(result.relative_level[0][q] == Catch::Approx(reciprocal.relative_level[0][q])
                                        .margin(steep ? 0.1 : 0.003));
                        CHECK(result.relative_level[0][q] == Catch::Approx(exact).margin(steep ? 0.1 : 0.02));
                        compared++;
                    }
                    CHECK(compared > 5);
                }
        }
}

TEST_CASE("the PE fills an upward-refraction shadow only by diffraction", "[pe]")
{
    // c falling 0.05 /s with height: the geometric shadow for a 10 m source
    // and a 1.2 m receiver starts at sqrt(2R)(sqrt(10) + sqrt(1.2)), R = c/|a|,
    // about 0.5 km.  Well past it, 1 kHz is down by tens of dB on free field;
    // well short of it, it is not.
    PeSetup setup;
    setup.frequency = 1000.0;
    setup.source_height = 10.0;
    setup.receiver_heights = {1.2};
    setup.max_range = 1500.0;
    setup.output_step = 50.0;
    setup.z = {0.0, 1000.0};
    setup.c = {340.0, 340.0 - 0.05 * 1000.0};
    GroundModel ground;
    setup.admittance = ground_impedance(ground, {1000.0}, 340.0)[0];
    const auto band = solve_pe_band(setup, 1000.0, 3);
    auto level_at = [&](double r) {
        for (std::size_t q = 0; q < band.range.size(); q++)
            if (band.range[q] >= r)
                return band.relative_level[0][q];
        return band.relative_level[0].back();
    };
    CHECK(level_at(200.0) > -10.0);
    CHECK(level_at(1400.0) < level_at(200.0) - 20.0);
}

TEST_CASE("a band average is the mean over the band, not over its edges", "[pe][band]")
{
    // solve_pe_band averages `samples` frequencies in energy across the third-octave
    // band.  They once sat at k/(N-1) of its width, which put the first and last on
    // the band edges with a full 1/N weight each: 0.4 to 0.6 dB rms off the band's
    // true mean in a case like this one (bias up to +0.3 dB), where the midpoints
    // are within 0.06.  The exact band mean over rigid ground in still air is the
    // two-ray energy integrated over log frequency.
    const double center = 1000.0, c = 340.0, z_s = 30.0, z_r = 1.2;
    PeSetup setup;
    setup.source_height = z_s;
    setup.receiver_heights = {z_r};
    setup.max_range = 500.0;
    setup.output_step = 10.0;
    setup.z = {0.0};
    setup.c = {c};
    const auto band = solve_pe_band(setup, center, 5);

    double sum_squares = 0.0, bias = 0.0;
    int compared = 0;
    for (std::size_t q = 0; q < band.range.size(); q++)
    {
        const double r = band.range[q];
        if (r < 100.0)
            continue; // the image ray's angle is toward the wide-angle limit nearer in
        double energy = 0.0;
        const int n = 2001;
        for (int i = 0; i < n; i++)
        {
            const double f = center * std::pow(2.0, ((i + 0.5) / n - 0.5) / 3.0);
            energy += std::pow(10.0, 0.1 * two_ray(r, z_s, z_r, f, c, 0.0, true));
        }
        const double exact = 10.0 * std::log10(energy / n);
        const double error = band.relative_level[0][q] - exact;
        sum_squares += error * error;
        bias += error;
        compared++;
    }
    REQUIRE(compared > 30);
    const double rms = std::sqrt(sum_squares / compared);
    INFO("rms error " << rms << " dB, bias " << bias / compared << " dB");
    REQUIRE(rms < 0.15);
    REQUIRE(std::fabs(bias / compared) < 0.1);
}

TEST_CASE("one band sample is the center frequency", "[pe][band]")
{
    PeSetup setup;
    setup.source_height = 10.0;
    setup.receiver_heights = {1.2};
    setup.max_range = 300.0;
    setup.output_step = 10.0;
    setup.z = {0.0};
    setup.c = {340.0};
    setup.frequency = 500.0;
    const auto single = solve_pe(setup);
    const auto band = solve_pe_band(setup, 500.0, 1);
    REQUIRE(band.range.size() == single.range.size());
    for (std::size_t q = 0; q < single.range.size(); q++)
        REQUIRE(band.relative_level[0][q] == Catch::Approx(single.relative_level[0][q]).margin(1e-9));
}

TEST_CASE("PE band averaging rejects nonpositive sample counts", "[pe][band]")
{
    PeSetup setup;
    REQUIRE_THROWS_AS(solve_pe_band(setup, 500.0, 0), std::invalid_argument);
    REQUIRE_THROWS_AS(solve_pe_band(setup, 500.0, -1), std::invalid_argument);
    REQUIRE_THROWS_AS(solve_pe_bands({}, {}, 0), std::invalid_argument);
    REQUIRE_THROWS_AS(solve_pe_bands({}, {}, -1), std::invalid_argument);
}

TEST_CASE("band reduction preserves single-frequency levels and invalid samples", "[pe][band]")
{
    for (const auto scheme : {PeScheme::split_step_pade, PeScheme::crank_nicolson})
    {
        PeSetup setup;
        setup.scheme = scheme;
        setup.source_height = 5.0;
        setup.receiver_heights = {0.0, 1.2, 30.0};
        setup.max_range = 40.0;
        setup.output_step = 10.0;
        setup.z = {0.0};
        setup.c = {340.0};
        const double center = 250.0;
        constexpr int samples = 3;
        std::vector<PeResult> tones;
        for (int k = 0; k < samples; k++)
        {
            setup.frequency = center * std::pow(2.0, ((k + 0.5) / samples - 0.5) / 3.0);
            tones.push_back(solve_pe(setup));
            REQUIRE(tones.back().psi.size() == setup.receiver_heights.size());
        }
        for (unsigned threads : {1u, 2u})
        {
            setup.threads = threads;
            const auto band = solve_pe_band(setup, center, samples);
            REQUIRE(band.psi.empty());
            REQUIRE(band.range == tones.front().range);
            for (std::size_t rcv = 0; rcv < setup.receiver_heights.size(); rcv++)
                for (std::size_t q = 0; q < band.range.size(); q++)
                {
                    double energy = 0.0;
                    for (const auto& tone : tones)
                        energy += std::pow(10.0, 0.1 * tone.relative_level[rcv][q]) / samples;
                    const double expected = 10.0 * std::log10(energy + 1e-300);
                    if (std::isnan(expected))
                        CHECK(std::isnan(band.relative_level[rcv][q]));
                    else
                        CHECK(band.relative_level[rcv][q] == expected);
                }
        }
    }
}
