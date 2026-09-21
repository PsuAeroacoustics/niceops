#include "faddeeva.hpp"

#include <cmath>
#include <iterator>

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

// Weideman's L for N = 32: sqrt(N / sqrt(2)).
constexpr double weideman_L = 4.7568284600108841;

// The 32 series coefficients, highest power first, as Horner wants them.
// Generated once from Weideman's construction -- an FFT of exp(-t^2)(L^2 + t^2)
// sampled at t = L tan(theta/2) -- rather than computed at startup, which would
// mean carrying an FFT for 32 numbers that never change.
constexpr double weideman_a[32] = {
    -1.30317978630500875e-12, +3.74107755472269789e-12,
    +8.03033538228436012e-12, -2.15435373143760700e-11,
    -5.54423561143486540e-11, +1.16582546820233235e-10,
    +4.15374271614691897e-10, -5.23102143724924582e-10,
    -3.20801514078875819e-09, +8.12488896568210460e-10,
    +2.37975567593842971e-08, +2.29304391194392872e-08,
    -1.48130789124396946e-07, -4.18407637053084035e-07,
    +4.25583313792292626e-07, +4.40153173159092505e-06,
    +6.82103194402497841e-06, -2.14096192018698432e-05,
    -1.30754492546239731e-04, -2.45329802700202065e-04,
    +3.92591360700705167e-04, +4.51954110534923473e-03,
    +1.90061557848454077e-02, +5.73044035298371918e-02,
    +1.40607162268937713e-01, +2.95444510715087261e-01,
    +5.46013972063934205e-01, +9.01925489364799993e-01,
    +1.34554416923454490e+00, +1.82566962963248147e+00,
    +2.26353729990026764e+00, +2.57225340812456960e+00,
};

const double one_over_sqrt_pi = 1.0 / std::sqrt(M_PI);

// The series proper, valid for Im(z) >= 0.
//
// Written in real arithmetic rather than std::complex.  Every std::complex
// multiply carries a branch to recover the NaN/infinity cases the standard
// requires, and the Horner loop below does 31 of them per call; the recovery
// path cannot trigger here, because the argument is bounded by construction.
// The two divides at the end are also folded into one:
//   2p/d^2 + c/d  ==  (2p + c d)/d^2.
std::complex<double> weideman(std::complex<double> z)
{
    // denominator d = L - i z, numerator n = L + i z.
    const double dr = weideman_L + z.imag();
    const double di = -z.real();
    const double nr = weideman_L - z.imag();
    const double ni = z.real();

    // t = n / d, the Moebius transform the series is a polynomial in.
    const double d_norm = dr * dr + di * di;
    const double tr = (nr * dr + ni * di) / d_norm;
    const double ti = (ni * dr - nr * di) / d_norm;

    // Four interleaved Horner chains over v = t^4 rather than one over t.
    //
    // The series is 31 complex multiply-adds deep and every one of them needs
    // the previous result, so the loop is bound by the latency of that chain
    // rather than by how many operations the core can issue -- it was running
    // at about a quarter of the machine's floating-point throughput.  Splitting
    // the polynomial by the power of t modulo four gives four chains a quarter
    // as deep that all run at once:
    //
    //     P(t) = C0(v) + t C1(v) + t^2 C2(v) + t^3 C3(v),    v = t^4
    //
    // where Cp collects the coefficients of the powers congruent to p.  The
    // table is left exactly as Weideman's construction produced it and indexed
    // by stride instead, so it stays checkable against the source: the chain
    // for t^p starts at a[3 - p] and steps by four from there.
    //
    // Same operations, regrouped, so the result moves by rounding only -- the
    // reference values in tests/test_ground_effect.cpp still hold to 1e-11.
    const double ur = tr * tr - ti * ti; // u = t^2
    const double ui = 2.0 * tr * ti;
    const double vr = ur * ur - ui * ui; // v = t^4
    const double vi = 2.0 * ur * ui;

    // The four chains only tile the table if its length is a multiple of four.
    static_assert(std::size(weideman_a) % 4 == 0,
            "the interleaved Horner chains tile weideman_a four at a time");
    double chain_r[4] = {weideman_a[3], weideman_a[2], weideman_a[1], weideman_a[0]};
    double chain_i[4] = {0.0, 0.0, 0.0, 0.0};
    for (std::size_t term = 4; term < std::size(weideman_a); term += 4)
    {
        for (std::size_t power = 0; power < 4; power++)
        {
            const double next
                    = chain_r[power] * vr - chain_i[power] * vi + weideman_a[term + 3 - power];
            chain_i[power] = chain_r[power] * vi + chain_i[power] * vr;
            chain_r[power] = next;
        }
    }

    // Recombine: (C0 + t C1) + u (C2 + t C3).
    const double lo_r = chain_r[0] + (tr * chain_r[1] - ti * chain_i[1]);
    const double lo_i = chain_i[0] + (tr * chain_i[1] + ti * chain_r[1]);
    const double hi_r = chain_r[2] + (tr * chain_r[3] - ti * chain_i[3]);
    const double hi_i = chain_i[2] + (tr * chain_i[3] + ti * chain_r[3]);
    const double pr = lo_r + (ur * hi_r - ui * hi_i);
    const double pi_ = lo_i + (ur * hi_i + ui * hi_r);

    // numerator = 2p + (1/sqrt(pi)) d
    const double num_r = 2.0 * pr + one_over_sqrt_pi * dr;
    const double num_i = 2.0 * pi_ + one_over_sqrt_pi * di;

    // d^2
    const double d2r = dr * dr - di * di;
    const double d2i = 2.0 * dr * di;

    const double d2_norm = d2r * d2r + d2i * d2i;
    return {(num_r * d2r + num_i * d2i) / d2_norm, (num_i * d2r - num_r * d2i) / d2_norm};
}

} // namespace

namespace niceops {

std::complex<double> faddeeva(std::complex<double> z)
{
    if (z.imag() >= 0.0)
        return weideman(z);

    // Reflect: w(z) = 2 exp(-z^2) - w(-z).  See the header on why this is safe
    // for the only caller.
    return 2.0 * std::exp(-(z * z)) - weideman(-z);
}

} // namespace niceops
