// The parabolic equation (pe.hpp) in an upward-refraction shadow, against the
// exact residue series.
//
// For k^2(z) = k0^2 (1 + 2 z / R) over rigid ground the field is exactly a sum of
// creeping-wave modes (Pierce, Acoustics, sec. 9-5; exact_solutions.cpp):
//   p(r, z) = pi i sum_n g_n(z) g_n(z_s) H0(kappa_n r) / N_n,
//   g_n(z) = Ai(a'_n - w z / l),  l = (R / 2 k0^2)^(1/3),  w = exp(2 pi i / 3),
//   kappa_n^2 = k0^2 + a'_n / (w l^2),  N_n = l a'_n Ai(a'_n)^2 / w,
// a'_n the zeros of Ai'.  The rigid ground puts g_n'(0) = 0, which is why the
// zeros of Ai' and not of Ai; kappa_n's imaginary part is each mode's decay
// with range, so deep in the shadow the first few modes are all of it.  With
// p normalized to a source whose free field is exp(ikR)/R, 20 log10(|p| R) is
// the level relative to free field, which is what solve_pe() reports.
// tests/test_benchmark.cpp takes the same series over impedance ground.
//
// This checks the PE's refraction and its diffraction into the shadow -- the
// part of it that arbitrates between the ray models -- against a solution
// that shares none of its approximations.  The PE is handed the same profile,
// c(z) = c0 / sqrt(1 + 2 z / R), sampled finely enough that its piecewise-linear
// interpolation is exact for the purpose.
//
// The residue series and the complex Airy function it needs live in
// exact_solutions.hpp and airy.hpp; the Airy function is checked here against
// values computed independently (Bessel-integral quadrature), and its two
// expansions against each other where they meet.

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

// Level relative to free field by the residue series.
double residue_level(double f, double c0, double R, double z_s, double z_r, double r)
{
    ExactProblem problem;
    problem.frequency = f;
    problem.sound_speed = c0;
    problem.source_height = z_s;
    problem.receiver_height = z_r;
    const complex p = upward_pressure(problem, R, {r})[0];
    return 20.0 * std::log10(std::abs(p) * std::hypot(r, z_r - z_s));
}

} // namespace

TEST_CASE("the complex Airy function matches tabulated values and itself", "[pe_shadow]")
{
    struct Value
    {
        double x, ai, aip;
    };
    // By quadrature of Ai = sqrt(x/3) K_{1/3}(zeta) / pi (x > 0) and
    // Ai(-y) = sqrt(y) (J_{1/3} + J_{-1/3})(zeta) / 3, and their derivatives,
    // zeta = 2/3 |x|^(3/2) -- independent of either expansion used here.
    for (const Value v : {Value {0.0, 0.355028053887817239, -0.258819403792806798},
                 Value {1.0, 1.352924163128726e-01, -1.591474412967858e-01},
                 Value {-1.0, 5.355608832923321e-01, -1.016056711663267e-02},
                 Value {5.0, 1.083444281360635e-04, -2.474138908684387e-04},
                 Value {-5.0, 3.507610090241173e-01, 3.271928185544626e-01},
                 Value {10.0, 1.104753255289774e-10, -3.520633676738718e-10},
                 Value {-10.0, 4.024123848644436e-02, 9.962650441327958e-01}})
    {
        complex ai, aip;
        airy(v.x, ai, aip);
        INFO("x = " << v.x);
        CHECK(ai.real() == Catch::Approx(v.ai).epsilon(1e-9));
        CHECK(aip.real() == Catch::Approx(v.aip).epsilon(1e-9));
        CHECK(std::fabs(ai.imag()) < 1e-12 * (std::fabs(v.ai) + 1e-30) + 1e-18);
    }
    // Series against asymptotics just past where they hand over, all the way
    // round (so through the connection formula's sector too).
    for (int k = 0; k < 24; k++)
    {
        const complex z = std::polar(5.2, -pi + 2.0 * pi * (k + 0.5) / 24.0);
        complex a1, d1, a2, d2;
        airy_series(z, a1, d1);
        airy(z, a2, d2); // |z| > 5: asymptotic (and connection)
        INFO("arg = " << std::arg(z));
        CHECK(std::abs(a1 - a2) <= 2e-6 * std::abs(a1));
        CHECK(std::abs(d1 - d2) <= 2e-6 * std::abs(d1));
    }
    // The zeros of Ai' (DLMF table 9.9.1).
    CHECK(airy_prime_zero(1) == Catch::Approx(-1.018792971647471).epsilon(1e-12));
    CHECK(airy_prime_zero(2) == Catch::Approx(-3.248197582179837).epsilon(1e-12));
    CHECK(airy_prime_zero(3) == Catch::Approx(-4.820099211178736).epsilon(1e-12));
    CHECK(airy_prime_zero(10) == Catch::Approx(-12.384788371845747).epsilon(1e-10));
}

TEST_CASE("the PE matches the residue series in an upward-refraction shadow", "[pe_shadow]")
{
    const double c0 = 340.0;
    const double R = 6800.0; // |dc/dz| ~ 0.05 /s near the ground
    const double z_s = 5.0, z_r = 1.5;
    // The geometric shadow begins at sqrt(2R)(sqrt(z_s) + sqrt(z_r)), ~0.4 km.
    const double boundary = std::sqrt(2.0 * R) * (std::sqrt(z_s) + std::sqrt(z_r));

    for (const double f : {250.0, 500.0, 1000.0})
    {
        PeSetup setup;
        setup.frequency = f;
        setup.source_height = z_s;
        setup.receiver_heights = {z_r};
        // Out to four times the boundary: 25 to 47 dB down, deep enough that
        // a spurious floor from the absorbing layer would show.
        setup.max_range = 4.0 * boundary;
        setup.output_step = 25.0;
        setup.admittance = 0.0; // rigid
        // c = c0 / sqrt(1 + 2z/R): k^2 exactly linear in z, sampled every 10 cm
        // up past the top of the PE's domain.
        for (double z = 0.0; z <= 800.0; z += 0.1)
        {
            setup.z.push_back(z);
            setup.c.push_back(c0 / std::sqrt(1.0 + 2.0 * z / R));
        }
        const auto pe = solve_pe(setup);
        int compared = 0;
        for (std::size_t q = 0; q < pe.range.size(); q++)
        {
            const double r = pe.range[q];
            // From just past the boundary, where the series converges quickly.
            if (r < 1.25 * boundary)
                continue;
            const double exact = residue_level(f, c0, R, z_s, z_r, r);
            INFO("f = " << f << ", r = " << r << ", exact = " << exact);
            CHECK(pe.relative_level[0][q] == Catch::Approx(exact).margin(0.1));
            compared++;
        }
        CHECK(compared > 10);
    }
}
