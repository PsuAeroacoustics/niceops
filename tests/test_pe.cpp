// The wide-angle parabolic equation (pe.hpp), against exact solutions.
//
// Over rigid ground in still air the field is two rays, direct and image.
// Over impedance ground it is the Weyl-van der Pol solution, whose spherical-
// wave reflection coefficient ground_effect.cpp already computes (and pins
// against panam).  The PE shares neither's derivation -- it marches a wave
// equation on a grid, with the impedance condition in the grid -- so agreement
// checks its normalisation, its ground condition's sign convention and its
// absorbing top together.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <complex>
#include <vector>

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
