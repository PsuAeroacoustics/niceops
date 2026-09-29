// The shadow zone's creeping-wave loss over impedance ground (ShadowAttenuation):
// the zero it is built on, its rigid limit, and its rate against the exact
// residue series it is the leading term of.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <complex>
#include <vector>

#include "airy.hpp"
#include "exact_solutions.hpp"
#include "ground_effect.hpp"
#include "refraction.hpp"
#include "util/constants.hpp"

using namespace niceops; // NOLINT

namespace {

using complex = std::complex<double>;

const complex w = std::polar(1.0, 2.0 * pi / 3.0);

// q = s d upward (exact_solutions.cpp): d = -i w^2 beta.
complex upward_direction(complex beta)
{
    return complex(0.0, -1.0) * w * w * beta;
}

// The benchmark's grounds (Attenborough et al. 1995), as tests/test_benchmark.cpp has them.
const complex ground_10 = 1.0 / complex(38.79, 38.41);
const complex ground_100 = 1.0 / complex(12.81, 11.62);
const complex ground_1k = 1.0 / complex(5.96, 2.46);

} // namespace

TEST_CASE("first_impedance_airy_zero follows the first zero as far as asked", "[shadow]")
{
    // q = 0 is the first zero of Ai'.
    CHECK(first_impedance_airy_zero(ground_1k, {0.0})[0] == complex(airy_prime_zero(1)));

    // Within the full finder's reach it is the zero that finder continues from
    // a'_1, which it returns first.
    for (const complex beta : {ground_10, ground_100, ground_1k})
    {
        const complex d = upward_direction(beta);
        const std::vector<double> s = {0.5, 2.0, 5.0 / std::abs(d)};
        const auto followed = first_impedance_airy_zero(d, s);
        for (std::size_t i = 0; i < s.size(); i++)
        {
            INFO("beta " << beta << ", |q| " << s[i] * std::abs(d));
            CHECK(std::abs(followed[i] - impedance_airy_zeros(s[i] * d, 1)[0]) < 1e-10);
        }
    }

    // Beyond it, every value is a zero, and they close on the first zero of
    // Ai as a_1 - 1/q - a_1/(3 q^3) + O(q^-4) -- down to the 1e-14 relative
    // Newton converges to, which is where 10 / |q|^4 goes by |q| ~ 1e4.
    constexpr double a_1 = -2.338107410459767;
    for (const complex beta : {ground_100, ground_1k, complex(0.3, -0.1), complex(1.0, 0.0)})
    {
        const complex d = upward_direction(beta);
        std::vector<double> s;
        for (double x = 1.0; x <= 1e5; x *= 3.0)
            s.push_back(x);
        const auto zeros = first_impedance_airy_zero(d, s);
        for (std::size_t i = 0; i < s.size(); i++)
        {
            const complex q = s[i] * d;
            complex ai, aip;
            airy(zeros[i], ai, aip);
            INFO("beta " << beta << ", q " << q << ", t " << zeros[i]);
            CHECK(std::abs(aip + q * ai) <= 1e-9 * (std::abs(aip) + std::abs(q * ai)));
            if (std::abs(q) > 100.0)
                CHECK(std::abs(zeros[i] - (a_1 - 1.0 / q - a_1 / (3.0 * q * q * q)))
                        < 10.0 / std::pow(std::abs(q), 4.0) + 1e-13);
        }
    }
}

TEST_CASE("over rigid ground ShadowAttenuation is the old rate", "[shadow]")
{
    const std::vector<double> bands = {50.0, 1000.0, 8000.0};
    const ShadowAttenuation rigid(bands, 0, bands.size());
    REQUIRE(rigid.rigid());
    // The rigid constructor is the old code: rates from shadow_attenuation_rates,
    // and the loss formed as it always was.
    const auto rates = shadow_attenuation_rates(bands, 0, bands.size());
    RefractedRay ray;
    ray.profile = {340.0, -0.037};
    ray.shadow_depth = 250.0;
    std::vector<double> loss(bands.size());
    shadow_loss(ray, rigid, 1e9, 0, bands.size(), loss.data());
    const double scale = std::cbrt(0.037 * 0.037) / 340.0 * 250.0;
    for (std::size_t b = 0; b < bands.size(); b++)
        CHECK(loss[b] == rates[b] * scale);

    // A table over zero admittance is the same rate, to the 9 digits the rigid
    // formula gives tau_1 (1.01879297; a'_1 is -1.0187929716...).
    const ShadowAttenuation zero(bands, std::vector<complex>(bands.size(), 0.0), 0, bands.size());
    REQUIRE_FALSE(zero.rigid());
    for (const double a : {1e-4, 0.01, 0.3, 5.0})
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            const LinearProfile profile {343.0, -a};
            CHECK(zero.rate(b, profile) == Catch::Approx(rigid.rate(b, profile)).epsilon(1e-8));
        }
}

TEST_CASE("the creeping wave decays as the exact residue series does", "[shadow]")
{
    // Benchmark case 3's geometry (source 5 m, receiver 1 m, |dc/dz| = 0.1 /s),
    // deep enough in the shadow that the first mode is all that is left: the
    // level re cylindrical spreading, 20 log|p| + 10 log r, falls by the rate.
    constexpr double c0 = 343.23, gradient = 0.1;
    struct Case
    {
        double frequency;
        complex beta;
        double r1, r2;
    };
    for (const Case c : {Case {100.0, 0.0, 1000.0, 1500.0}, Case {100.0, ground_100, 1000.0, 1500.0},
                 Case {250.0, ground_100, 800.0, 1200.0}, Case {1000.0, 0.0, 600.0, 900.0},
                 Case {1000.0, ground_1k, 600.0, 900.0}})
    {
        ExactProblem problem;
        problem.frequency = c.frequency;
        problem.sound_speed = c0;
        problem.admittance = c.beta;
        problem.source_height = 5.0;
        problem.receiver_height = 1.0;
        std::vector<double> ranges;
        for (int i = 0; i <= 10; i++)
            ranges.push_back(c.r1 + (c.r2 - c.r1) * i / 10.0);
        const auto field = upward_pressure(problem, c0 / gradient, ranges);
        double mean_r = 0.0, mean_level = 0.0;
        std::vector<double> level;
        for (std::size_t i = 0; i < ranges.size(); i++)
        {
            REQUIRE(std::isfinite(field[i].real()));
            level.push_back(20.0 * std::log10(std::abs(field[i])) + 10.0 * std::log10(ranges[i]));
            mean_r += ranges[i] / ranges.size();
            mean_level += level.back() / ranges.size();
        }
        double sxy = 0.0, sxx = 0.0;
        for (std::size_t i = 0; i < ranges.size(); i++)
        {
            sxy += (ranges[i] - mean_r) * (level[i] - mean_level);
            sxx += (ranges[i] - mean_r) * (ranges[i] - mean_r);
        }
        const double exact = -sxy / sxx;
        const ShadowAttenuation table({c.frequency}, {c.beta}, 0, 1);
        const double rate = table.rate(0, {c0, -gradient});
        INFO(c.frequency << " Hz, beta " << c.beta << ": exact " << exact << " dB/m, table " << rate);
        // Measured 0.02-0.3%: Pierce's first-order form, and the next mode's remainder.
        CHECK(rate == Catch::Approx(exact).epsilon(0.005));
        if (c.beta != 0.0)
            CHECK(rate > 1.1 * ShadowAttenuation({c.frequency}, 0, 1).rate(0, {c0, -gradient}));
    }
}

TEST_CASE("the decay factor leaves rigid to first order and ends at pressure release", "[shadow]")
{
    // Im(b / w): rigid's tau_1 sin 60 at q = 0, shifted to first order by
    // -q / a'_1, and tending to the first zero of Ai's 2.338 sin 60 as |q|
    // grows.  Not monotonic between: an ordinary ground first attenuates less
    // than rigid (the shift's sign), then more.  Measured over every model in
    // the next test, it stays within [0.774, 2.025].
    const double a1p = airy_prime_zero(1);
    const double rigid = -a1p * std::sqrt(3.0) / 2.0;
    const double release = 2.338107410459767 * std::sqrt(3.0) / 2.0;
    const std::vector<double> bands = {1000.0};
    for (const complex beta : {ground_10, ground_100, ground_1k})
    {
        const ShadowAttenuation table(bands, {beta}, 0, 1);
        const complex d = upward_direction(beta);
        // First order in q, at the table's first node: O(|q|^2) off.
        const double s0 = ShadowAttenuation::s_min;
        const double first_order = ((a1p - s0 * d / a1p) / w).imag();
        CHECK(std::abs(table.decay_factor(0, s0) - first_order) < 2.0 * std::norm(s0 * d));
        for (double s = s0; s <= ShadowAttenuation::s_max; s *= 2.0)
        {
            const double g = table.decay_factor(0, s);
            INFO("beta " << beta << ", s " << s << ": " << g);
            CHECK(g > 0.87 * rigid);
            CHECK(g < release + 1e-9);
        }
        // |q| = 10^5 |beta| >= 1800: within 1/|q| of the limit.
        CHECK(table.decay_factor(0, ShadowAttenuation::s_max)
                == Catch::Approx(release).epsilon(1.0 / (1e5 * std::abs(beta))));
    }
}

TEST_CASE("the usual impedance models leave no surface wave in the shadow", "[shadow]")
{
    // What ShadowAttenuation leaves out needs arg beta < -60 degrees.
    std::vector<double> bands;
    for (double f = 10.0; f <= 10000.0 * 1.001; f *= std::pow(10.0, 0.1))
        bands.push_back(f);
    for (const auto model : {ImpedanceModel::delany_bazley, ImpedanceModel::miki,
                 ImpedanceModel::variable_porosity})
        for (const double sigma : {10.0, 30.0, 100.0, 200.0, 1000.0, 10000.0})
        {
            GroundModel ground;
            ground.impedance_model = model;
            ground.flow_resistance = sigma;
            const auto beta = ground_impedance(ground, bands, 343.0);
            const ShadowAttenuation table(bands, beta, 0, bands.size());
            for (std::size_t b = 0; b < bands.size(); b++)
            {
                INFO("model " << static_cast<int>(model) << ", sigma " << sigma << ", " << bands[b]
                              << " Hz: arg beta " << std::arg(beta[b]) * 180.0 / pi);
                CHECK(std::arg(beta[b]) > -60.0 * pi / 180.0);
                CHECK_FALSE(table.admits_surface_wave(b));
            }
        }

    // A strong porosity rate does admit one, in the low bands, and says so.
    // There the first zero turns into the surface wave and cannot be followed,
    // so the band takes the rigid rate.
    GroundModel graded;
    graded.impedance_model = ImpedanceModel::variable_porosity;
    graded.flow_resistance = 10.0;
    graded.porosity_rate = 100.0;
    const ShadowAttenuation table(bands, ground_impedance(graded, bands, 343.0), 0, bands.size());
    const ShadowAttenuation rigid(bands, 0, bands.size());
    CHECK(table.admits_surface_wave(0));
    CHECK_FALSE(table.admits_surface_wave(bands.size() - 1));
    const LinearProfile profile {343.0, -0.05};
    CHECK(table.rate(0, profile) == rigid.rate(0, profile));
    CHECK(table.rate(bands.size() - 1, profile) != rigid.rate(bands.size() - 1, profile));
}
