// The parabolic equation (pe.hpp) against the exact solutions of the
// Attenborough et al. benchmark (JASA 97, 173 (1995); exact_solutions.hpp,
// tools/benchmark.cpp), on the benchmark's own ground: source 5 m, receiver
// 1 m, the impedances their four-parameter model gives at 10, 100 and 1000 Hz.
//
//   case 1  uniform air: the direct and image rays with the spherical-wave
//           reflection coefficient;
//   case 2  downward refraction: the normal-mode series;
//   case 3  upward refraction: the residue series, which is where their
//           printed eq. (11) and the derivation in exact_solutions.cpp part
//           company.  With the ground condition as printed the series misses
//           the PE by 1.5 dB at 100 m and 12 dB at 1 km at 100 Hz; as derived,
//           by hundredths.
//
// The refracting cases hand the PE the exact solutions' profiles,
// c0 / sqrt(1 -+ 2 z / R), rather than the benchmark's linear ones, so that
// both sides solve one problem.  The mode sums themselves are checked first:
// their roots and Hankel function against values from Mathematica.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <complex>
#include <vector>

#include "airy.hpp"
#include "exact_solutions.hpp"
#include "pe.hpp"
#include "util/constants.hpp"

using namespace niceops;

namespace {

using complex = std::complex<double>;

constexpr double c0 = 343.23;
constexpr double radius = c0 / 0.1; // R for |dc/dz| = 0.1 /s
constexpr double z_s = 5.0, z_r = 1.0;

complex benchmark_admittance(double f)
{
    if (f == 10.0)
        return 1.0 / complex(38.79, 38.41);
    if (f == 100.0)
        return 1.0 / complex(12.81, 11.62);
    return 1.0 / complex(5.96, 2.46);
}

ExactProblem problem_at(double f)
{
    ExactProblem p;
    p.frequency = f;
    p.sound_speed = c0;
    p.admittance = benchmark_admittance(f);
    p.source_height = z_s;
    p.receiver_height = z_r;
    return p;
}

PeSetup pe_at(double f, double max_range, double step)
{
    PeSetup s;
    s.frequency = f;
    s.source_height = z_s;
    s.receiver_heights = {z_r};
    s.max_range = max_range;
    s.output_step = step;
    s.admittance = benchmark_admittance(f);
    s.points_per_wavelength = std::max(10.0, c0 / f / 0.5);
    return s;
}

// c0 / sqrt(1 + sign 2 z / R), every 10 cm to `top`.
void airy_profile(PeSetup& s, double sign, double top)
{
    for (double z = 0.0; z <= top; z += 0.1)
    {
        s.z.push_back(z);
        s.c.push_back(c0 / std::sqrt(1.0 + sign * 2.0 * z / radius));
    }
}

double relative(complex p, double r)
{
    return 20.0 * std::log10(std::abs(p) * std::hypot(r, z_r - z_s));
}

} // namespace

TEST_CASE("H0 matches Mathematica on both sides of its switch", "[benchmark]")
{
    struct Value
    {
        complex z, h;
    };
    for (const Value v : {Value {1.0, {0.76519768655796655, 0.08825696421567696}},
                 Value {10.0, {-0.24593576445134834, 0.05567116728359939}},
                 Value {{0.0, 1.0}, {0.0, -0.26803248203398855}},
                 Value {{9.5, 2.0}, {-0.02341011560860242, 0.02540327080644915}},
                 Value {{10.5, 2.0}, {-0.03238163874809728, -0.00604197949931342}},
                 Value {{3.0, 6.0}, {-0.00006218617372164, 0.00074900048139330}},
                 Value {{25.0, 0.1}, {0.08687335541787914, -0.11531119073987443}},
                 Value {{0.2, 0.05}, {0.82761108075528167, -1.06725436608890377}}})
    {
        INFO("z = " << v.z);
        CHECK(std::abs(hankel1_0(v.z) - v.h) <= 1e-9 * std::abs(v.h));
    }
}

TEST_CASE("the zeros of Ai' + q Ai are found, and none is missed", "[benchmark]")
{
    // q = 0: the zeros of Ai' themselves.
    const auto rigid = impedance_airy_zeros(0.0, 5);
    REQUIRE(rigid.size() == 5);
    for (int n = 1; n <= 5; n++)
        CHECK(std::abs(rigid[n - 1] - airy_prime_zero(n)) < 1e-12);

    // Case 2's ground at 100 Hz, q = i k0 l beta: the lowest zero has
    // traveled half-way to the imaginary axis.  Mathematica's, by FindRoot.
    const complex q1(0.57, 0.628);
    const auto zeros = impedance_airy_zeros(q1, 3);
    for (const complex expected : {complex(-0.55900779524144920, 1.00743162226423679),
                 complex(-3.05769022524916337, 0.19095592088349130),
                 complex(-4.69512495599364548, 0.12774278824241452)})
    {
        double nearest = 1e300;
        for (const complex z : zeros)
            nearest = std::min(nearest, std::abs(z - expected));
        INFO("expected " << expected);
        CHECK(nearest < 1e-11);
    }

    // At 1 kHz a zero appears that continues from no zero of Ai': the
    // surface wave, near q^2.
    const complex q2(1.86, 4.51);
    const auto soft = impedance_airy_zeros(q2, 3);
    double nearest = 1e300;
    for (const complex z : soft)
        nearest = std::min(nearest, std::abs(z - complex(-16.91957380202532, 16.87217198826913)));
    CHECK(nearest < 1e-10);

    // Every zero returned is one, over the reach the finder is built for
    // (|q| <= 6; the benchmark's grounds come to 4.9).
    for (const complex q : {q1, q2, complex(-2.97, 3.86), complex(4.0, -4.4), complex(-3.0, -5.0)})
    {
        for (const complex t : impedance_airy_zeros(q, 40))
        {
            complex ai, aip;
            airy(t, ai, aip);
            INFO("q = " << q << ", t = " << t);
            CHECK(std::abs(aip + q * ai) <= 1e-9 * (std::abs(aip) + std::abs(q * ai)));
        }
    }
}

TEST_CASE("the PE's absorption is the medium's", "[benchmark]")
{
    // Rigid ground, uniform air, 5.4 dB/km: the level falls by the absorption
    // along the paths, which near grazing is the range's.
    PeSetup s = pe_at(1000.0, 1000.0, 50.0);
    s.admittance = 0.0;
    s.z = {0.0};
    s.c = {c0};
    const auto lossless = solve_pe(s);
    s.absorption = 5.4e-3 / 8.685889638;
    const auto absorbing = solve_pe(s);
    for (std::size_t q = 0; q < lossless.range.size(); q++)
    {
        const double r = lossless.range[q];
        if (r < 100.0)
            continue;
        INFO("r = " << r);
        CHECK(absorbing.relative_level[0][q] - lossless.relative_level[0][q]
                == Catch::Approx(-5.4e-3 * r).margin(0.02));
    }
}

TEST_CASE("case 1: the PE is the direct and image rays over the benchmark ground", "[benchmark]")
{
    for (const double f : {10.0, 100.0, 1000.0})
    {
        PeSetup s = pe_at(f, 1000.0, 5.0);
        s.z = {0.0};
        s.c = {c0};
        const auto pe = solve_pe(s);
        const auto p = problem_at(f);
        int compared = 0;
        for (std::size_t q = 0; q < pe.range.size(); q++)
        {
            const double r = pe.range[q];
            // Nearer than 20 m the image ray is steeper than the wide-angle
            // limit.  At 10 Hz that is still inside a wavelength and a half
            // of the source, where the field is not yet the far field the
            // starter stands for: the march, which no longer adds an error of
            // its own, is 0.29 dB over at 25 m and inside 0.2 from 45.
            if (r < (f == 10.0 ? 50.0 : 20.0))
                continue;
            const double exact = relative(homogeneous_pressure(p, r), r);
            INFO("f = " << f << ", r = " << r);
            CHECK(pe.relative_level[0][q] == Catch::Approx(exact).margin(0.2));
            compared++;
        }
        CHECK(compared > 150);
    }
}

TEST_CASE("case 2: the PE is the normal-mode series in downward refraction", "[benchmark]")
{
    // Their fig. 4 is 100 Hz to 200 m; this goes to 1 km.  The series is cut
    // at 75 degrees (exact_solutions.hpp).  Uncut, the near-vertical rays that
    // turn at z = R/2 fall back among these ranges -- 5 dB about the PE by
    // 1 km at 100 Hz -- and the PE rightly has none of them.  What reaches
    // 1 km in one hop leaves below 8.4 degrees, what falls back there above
    // 81.6, so the cut takes the one and leaves the other.
    //
    // At 10 Hz there are only ~60 modes below the cut, too few for it to be
    // clean within 6 wavelengths of the source: there the full series, whose
    // returns are still weak, and the cut one only from 500 m; and from 50 m,
    // past the source's near field (case 1).  At 100 Hz the grid is ten points a
    // wavelength, 0.34 m under a 1 m receiver, and the PE closes on the series
    // as 1/r: 0.16 dB under at 15 m, 0.10 at 20, 0.07 at 25, where the
    // comparison starts.  At 1 kHz the
    // PE's own grid is the limit past 500 m -- 0.8 dB at ten points a
    // wavelength by 1 km, 0.3 at twenty -- so the comparison stops there.
    struct Leg
    {
        double f, from, to, cut;
    };
    for (const Leg leg : {Leg {10.0, 50.0, 200.0, 90.0}, Leg {10.0, 500.0, 1000.0, 75.0},
                 Leg {100.0, 25.0, 1000.0, 75.0}, Leg {1000.0, 20.0, 500.0, 75.0}})
    {
        PeSetup s = pe_at(leg.f, leg.to, 5.0);
        airy_profile(s, -1.0, 1600.0);
        const auto pe = solve_pe(s);
        const auto modes = downward_pressure(problem_at(leg.f), radius, pe.range, leg.cut);
        int compared = 0;
        for (std::size_t q = 0; q < pe.range.size(); q++)
        {
            const double r = pe.range[q];
            if (r < leg.from || std::isnan(modes[q].real()))
                continue;
            const double exact = relative(modes[q], r);
            INFO("f = " << leg.f << ", r = " << r << ", cut at " << leg.cut);
            CHECK(pe.relative_level[0][q] == Catch::Approx(exact).margin(leg.f == 10.0 ? 0.2 : 0.1));
            compared++;
        }
        CHECK(compared > 30);
    }
}

TEST_CASE("case 3: the PE is the residue series over impedance ground", "[benchmark]")
{
    // Past the geometric shadow boundary, sqrt(2R)(sqrt(z_s) + sqrt(z_r)) =
    // 270 m, to 60 dB down (1 kHz, on this soft a ground, gets there in
    // 600 m).
    for (const double f : {100.0, 1000.0})
    {
        PeSetup s = pe_at(f, 1000.0, 25.0);
        airy_profile(s, 1.0, 1200.0);
        const auto pe = solve_pe(s);
        const auto series = upward_pressure(problem_at(f), radius, pe.range);
        int compared = 0;
        for (std::size_t q = 0; q < pe.range.size(); q++)
        {
            const double r = pe.range[q];
            if (r < 100.0 || std::isnan(series[q].real()))
                continue;
            const double exact = relative(series[q], r);
            if (exact < -60.0)
                continue;
            INFO("f = " << f << ", r = " << r << ", exact = " << exact);
            CHECK(pe.relative_level[0][q] == Catch::Approx(exact).margin(0.1));
            compared++;
        }
        CHECK(compared > 15);
    }
}
