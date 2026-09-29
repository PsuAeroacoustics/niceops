#include "faddeeva.hpp"

#include "util/constants.hpp"

#include <cmath>
#include <iterator>
#include <vector>

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

const double one_over_sqrt_pi = 1.0 / std::sqrt(pi);

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

std::complex<double> weideman_everywhere(std::complex<double> z)
{
    if (z.imag() >= 0.0)
        return weideman(z);
    // Reflect: w(z) = 2 exp(-z^2) - w(-z).  See the header on why this is safe
    // for the only caller.
    return 2.0 * std::exp(-(z * z)) - weideman(-z);
}

// ---- Fast paths over the domain the ground model reaches -------------------
//
// Sampled from real runs (the Be407 demo, point and plate receivers, 480,000
// arguments), z lies in the first quadrant, Re(z) in [0.3, 22.4] and
// Im(z) >= -0.07, with |z| >= 6 for 44% of calls.  Two methods cover it at a
// fraction of the Weideman series' cost, each measured against SciPy's wofz
// on those arguments:
//
//  * |z| >= 6: the asymptotic series
//        w(z) ~ i/(sqrt(pi) z) sum_n (2n-1)!! / (2 z^2)^n,
//    12 terms, 1e-11 relative.  Valid off the negative imaginary axis; the
//    exp(-z^2) it omits is below 1e-15 wherever it is used (x^2 - y^2 >= 36
//    or y >= 0).
//  * |z| < 6 inside the table: Taylor series of order 6 about the nearest node
//    of a 0.05 grid over Re in [0, 6], Im in [-0.25, 6], 3.6e-12 relative
//    (further below the real axis, where w grows, it loses digits: 1e-10 at
//    Im = -1.1).
//    Derivatives by w' = -2 z w + 2i/sqrt(pi) and
//    w^(n+1) = -2 z w^(n) - 2 n w^(n-1), from the node's w alone, which the
//    Weideman series fills once at start-up (15,246 nodes, 240 KB).
//
// Anything else falls back to the Weideman series.
constexpr int asymptotic_terms = 12;
constexpr double asymptotic_radius_squared = 36.0;
constexpr double table_step = 0.05;
constexpr double table_x_max = 6.0;
constexpr double table_y_min = -0.25;
constexpr double table_y_max = 6.0;
constexpr int table_nx = static_cast<int>(table_x_max / table_step + 0.5) + 1;
constexpr int table_ny = static_cast<int>((table_y_max - table_y_min) / table_step + 0.5) + 1;
constexpr int taylor_order = 6;

struct NodeTable
{
    std::vector<double> re;
    std::vector<double> im;
    NodeTable() : re(table_nx * table_ny), im(table_nx * table_ny)
    {
        for (int j = 0; j < table_ny; j++)
        {
            for (int i = 0; i < table_nx; i++)
            {
                const auto w = weideman_everywhere({i * table_step, table_y_min + j * table_step});
                re[j * table_nx + i] = w.real();
                im[j * table_nx + i] = w.imag();
            }
        }
    }
};

const NodeTable& node_table()
{
    static const NodeTable table;
    return table;
}

inline std::complex<double> asymptotic(double x, double y)
{
    // u = 1 / (2 z^2); sum by Horner in u, then times i / (sqrt(pi) z).
    const double zr2 = x * x - y * y;
    const double zi2 = 2.0 * x * y;
    const double norm2 = zr2 * zr2 + zi2 * zi2;
    const double ur = 0.5 * zr2 / norm2;
    const double ui = -0.5 * zi2 / norm2;
    // (2n-1)!! for n = 0..asymptotic_terms
    static constexpr double coefficient[asymptotic_terms + 1] = {1.0, 1.0, 3.0, 15.0, 105.0, 945.0,
            10395.0, 135135.0, 2027025.0, 34459425.0, 654729075.0, 13749310575.0,
            316234143225.0};
    double sr = coefficient[asymptotic_terms];
    double si = 0.0;
    for (int n = asymptotic_terms - 1; n >= 0; n--)
    {
        const double tr = sr * ur - si * ui + coefficient[n];
        si = sr * ui + si * ur;
        sr = tr;
    }
    // i s / (sqrt(pi) z) = i s conj(z) / (sqrt(pi) |z|^2)
    const double scale = one_over_sqrt_pi / (x * x + y * y);
    const double pr = sr * x + si * y;          // s * conj(z)
    const double pi_ = si * x - sr * y;
    return {-pi_ * scale, pr * scale};
}

inline std::complex<double> taylor(double x, double y)
{
    const auto& table = node_table();
    const int i = static_cast<int>(x / table_step + 0.5);
    const int j = static_cast<int>((y - table_y_min) / table_step + 0.5);
    const double x0 = i * table_step;
    const double y0 = table_y_min + j * table_step;
    const double dr = x - x0;
    const double di = y - y0;
    // Derivatives w^(n) / n! by the recurrence, divided as they go.
    double cr[taylor_order + 1];
    double ci[taylor_order + 1];
    cr[0] = table.re[j * table_nx + i];
    ci[0] = table.im[j * table_nx + i];
    // w' = -2 z0 w + 2i/sqrt(pi)
    cr[1] = -2.0 * (x0 * cr[0] - y0 * ci[0]);
    ci[1] = -2.0 * (x0 * ci[0] + y0 * cr[0]) + 2.0 * one_over_sqrt_pi;
    // c_{n+1} = (-2 z0 c_n - 2 c_{n-1}) / (n + 1), with c_n = w^(n)/n!
    for (int n = 1; n < taylor_order; n++)
    {
        const double inv = 1.0 / (n + 1);
        cr[n + 1] = (-2.0 * (x0 * cr[n] - y0 * ci[n]) - 2.0 * cr[n - 1]) * inv;
        ci[n + 1] = (-2.0 * (x0 * ci[n] + y0 * cr[n]) - 2.0 * ci[n - 1]) * inv;
    }
    double sr = cr[taylor_order];
    double si = ci[taylor_order];
    for (int n = taylor_order - 1; n >= 0; n--)
    {
        const double tr = sr * dr - si * di + cr[n];
        si = sr * di + si * dr + ci[n];
        sr = tr;
    }
    return {sr, si};
}

} // namespace

namespace niceops {

std::complex<double> faddeeva_reference(std::complex<double> z)
{
    return weideman_everywhere(z);
}

std::complex<double> faddeeva(std::complex<double> z)
{
    const double x = z.real();
    const double y = z.imag();
    const double r2 = x * x + y * y;
    if (r2 >= asymptotic_radius_squared)
    {
        if (y >= 0.0 || x * x - y * y >= asymptotic_radius_squared)
            return asymptotic(x, y);
    }
    else if (x >= 0.0 && y >= table_y_min && y <= table_y_max)
    {
        return taylor(x, y);
    }
    return weideman_everywhere(z);
}

} // namespace niceops
