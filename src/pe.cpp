#include "pe.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>

#include "pade_propagator.hpp"
#include "util/constants.hpp"
#include "util/parallel.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

double speed_at(const PeSetup& s, double z)
{
    if (s.z.size() == 1 || z <= s.z.front())
        return s.c.front();
    if (z >= s.z.back())
        return s.c.back();
    const auto upper = std::upper_bound(s.z.begin(), s.z.end(), z);
    const std::size_t i = static_cast<std::size_t>(upper - s.z.begin()) - 1;
    return s.c[i] + (z - s.z[i]) / (s.z[i + 1] - s.z[i]) * (s.c[i + 1] - s.c[i]);
}

// In-place radix-2 transform, x_j <- sum_n x_n exp(+2 pi i n j / N); N a power of 2.
void inverse_fft(std::vector<complex>& x)
{
    const std::size_t n = x.size();
    for (std::size_t i = 1, j = 0; i < n; i++)
    {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            std::swap(x[i], x[j]);
    }
    for (std::size_t len = 2; len <= n; len <<= 1)
    {
        const double angle = 2.0 * pi / static_cast<double>(len);
        const complex w_len(std::cos(angle), std::sin(angle));
        for (std::size_t i = 0; i < n; i += len)
        {
            complex w(1.0, 0.0);
            for (std::size_t k = 0; k < len / 2; k++)
            {
                const complex u = x[i + k];
                const complex v = x[i + k + len / 2] * w;
                x[i + k] = u + v;
                x[i + k + len / 2] = u - v;
                w *= w_len;
            }
        }
    }
}

// Screens for one march: the spectrum's amplitudes once, a fresh draw of
// phases per screen.
class ScreenSource
{
public:
    ScreenSource(double k, double gamma_t, double outer_scale, double thickness, double dz,
            std::size_t count, unsigned long long seed)
        : count_(count), rng_(seed)
    {
        // Twice the column, so the screen's period does not tie its top to its
        // bottom.
        std::size_t size = 1;
        while (size < 2 * count)
            size <<= 1;
        const double dkappa = 2.0 * pi / (static_cast<double>(size) * dz);
        const double c2 = 0.25 * gamma_t;
        const double kappa0 = 1.0 / outer_scale;
        // S_phi(kappa_z) = 2 pi k^2 dx 0.033 C^2 sqrt(pi) Gamma(4/3)/Gamma(11/6)
        //                  (kappa_z^2 + kappa0^2)^(-4/3), two-sided.
        const double scale = 2.0 * pi * k * k * thickness * 0.033 * c2 * std::sqrt(pi)
                * std::tgamma(4.0 / 3.0) / std::tgamma(11.0 / 6.0);
        amplitude_.assign(size / 2, 0.0);
        for (std::size_t n = 1; n < size / 2; n++)
        {
            const double kappa = static_cast<double>(n) * dkappa;
            const double spectrum = scale * std::pow(kappa * kappa + kappa0 * kappa0, -4.0 / 3.0);
            // A cosine of amplitude a has variance a^2 / 2; each carries both
            // signs of kappa_z.
            amplitude_[n] = 2.0 * std::sqrt(spectrum * dkappa);
        }
        buffer_.resize(size);
    }

    const std::vector<double>& next()
    {
        std::uniform_real_distribution<double> phase(0.0, 2.0 * pi);
        std::fill(buffer_.begin(), buffer_.end(), complex(0.0));
        for (std::size_t n = 1; n < amplitude_.size(); n++)
            buffer_[n] = std::polar(amplitude_[n], phase(rng_));
        inverse_fft(buffer_);
        screen_.resize(count_);
        for (std::size_t j = 0; j < count_; j++)
            screen_[j] = buffer_[j].real();
        return screen_;
    }

private:
    std::size_t count_;
    std::mt19937_64 rng_;
    std::vector<double> amplitude_;
    std::vector<complex> buffer_;
    std::vector<double> screen_;
};

} // namespace

std::vector<double> phase_screen(double k, double gamma_t, double outer_scale, double thickness,
        double dz, std::size_t count, unsigned long long seed)
{
    ScreenSource source(k, gamma_t, outer_scale, thickness, dz, count, seed);
    return source.next();
}

namespace {

// One range step, psi(r) -> psi(r + dr), on rows psi_1 .. psi_{n-1}.
using March = std::function<void(std::vector<complex>&)>;

// The vertical operator both schemes share: Q psi_j = a (psi_{j+1} - 2 psi_j +
// psi_{j-1}) + epsilon_j psi_j on rows j = 1 .. n-1, with psi_0 = g1 psi_1 +
// g2 psi_2 from the ground condition and psi_n = 0.
struct VerticalOperator
{
    std::vector<complex> epsilon; // (k_j^2 - k_a^2) / k_a^2 at j = 0 .. n
    double a = 0.0; // 1 / (k_a dz)^2
    complex g1, g2;
    std::size_t m = 0; // unknowns psi_1 .. psi_{n-1}
    // The fourth-order ground condition (pade_march): psi_0 is then carried
    // as psi[m], and the condition is psi' = (ground_slope / dz) psi.
    bool fourth_ground = false;
    complex ground_slope;
};

// Crank-Nicolson on the Pade(1,1) wide-angle equation:
//   [1 + (1/4 - i k_a dr/4) Q] psi(r + dr) = [1 + (1/4 + i k_a dr/4) Q] psi(r).
March crank_nicolson_march(const VerticalOperator& op, double k_a, double dr)
{
    const std::vector<complex>& epsilon = op.epsilon;
    const double a = op.a;
    const complex g1 = op.g1, g2 = op.g2;
    const std::size_t m = op.m;
    const complex i(0.0, 1.0);
    const complex alpha_lhs = 0.25 - i * k_a * dr / 4.0;
    const complex alpha_rhs = 0.25 + i * k_a * dr / 4.0;
    auto diagonal = [&](complex alpha, std::size_t j) { return 1.0 + alpha * (-2.0 * a + epsilon[j]); };

    // The left side does not change with range: factor it once, from both ends
    // at once (a twisted, or burn-at-both-ends, factorization).  Rows 0 .. p-1
    // are eliminated downward and p+1 .. m-1 upward, meeting at row p; the two
    // elimination chains, and then the two back-substitution chains, are
    // independent, so a step's two sweeps each carry two recurrences the
    // processor can overlap.  Each recurrence is one complex multiply-add per
    // row -- the pivot is folded into the stored coefficients -- where the
    // plain Thomas sweep chained two multiplies through every row and ran at
    // the latency of both.
    const complex off_lhs = alpha_lhs * a; // the sub- and super-diagonal
    std::vector<complex> main(m);
    for (std::size_t row = 0; row < m; row++)
        main[row] = diagonal(alpha_lhs, row + 1);
    main[0] += off_lhs * g1;
    const complex upper_0 = off_lhs + off_lhs * g2;
    auto lower_of = [&](std::size_t row) { return row > 0 ? off_lhs : complex(0.0); };
    auto upper_of = [&](std::size_t row) { return row + 1 < m ? (row == 0 ? upper_0 : off_lhs) : complex(0.0); };
    if (m < 3)
        throw std::invalid_argument("solve_pe: the grid is too small");
    const std::size_t p = (m - 1) / 2;
    // Downward, rows r < p:  w_r = b_r inv_r - lower_r inv_r w_{r-1},
    // x_r = w_r - (upper_r inv_r) x_{r+1}.  Upward, rows r > p: the mirror.
    std::vector<complex> inv(m), chain(m), back(m);
    {
        complex previous = 0.0; // upper_{r-1} inv_{r-1}
        for (std::size_t row = 0; row < p; row++)
        {
            const complex d = main[row] - lower_of(row) * previous;
            inv[row] = 1.0 / d;
            chain[row] = lower_of(row) * inv[row];
            back[row] = upper_of(row) * inv[row];
            previous = back[row];
        }
    }
    {
        complex previous = 0.0; // lower_{r+1} inv_{r+1}
        for (std::size_t row = m; row-- > p + 1;)
        {
            const complex d = main[row] - upper_of(row) * previous;
            inv[row] = 1.0 / d;
            chain[row] = upper_of(row) * inv[row];
            back[row] = lower_of(row) * inv[row];
            previous = back[row];
        }
    }
    // The meeting row takes both.
    inv[p] = 1.0 / (main[p] - lower_of(p) * (p > 0 ? back[p - 1] : complex(0.0))
                    - upper_of(p) * (p + 1 < m ? back[p + 1] : complex(0.0)));
    const complex lower_p = lower_of(p), upper_p = upper_of(p);
    const complex rhs_off = alpha_rhs * a;
    std::vector<complex> rhs_main(m);
    for (std::size_t row = 0; row < m; row++)
        rhs_main[row] = diagonal(alpha_rhs, row + 1);
    const complex rhs_main_0 = rhs_main[0] + alpha_rhs * a * g1;
    const complex rhs_upper_0 = alpha_rhs * a + alpha_rhs * a * g2;

    std::vector<complex> w(m);
    return [=](std::vector<complex>& psi) mutable {
        // Complex arithmetic written out: std::complex's multiply carries a branch
        // for NaN and infinity that cannot arise here, and it sits in the chain.
        auto mul = [](complex x, complex y) {
            return complex(x.real() * y.real() - x.imag() * y.imag(), x.real() * y.imag() + x.imag() * y.real());
        };
        // e - c x, as two dependent fused multiply-adds per part: the recurrences
        // run at their latency, and the compiler otherwise leaves a subtract at the
        // end of each.
        auto minus_mul = [](complex e, complex c, complex x) {
            return complex(std::fma(c.imag(), x.imag(), std::fma(-c.real(), x.real(), e.real())),
                    std::fma(-c.imag(), x.real(), std::fma(-c.real(), x.imag(), e.imag())));
        };
        // Right side [1 + alpha_rhs Q] psi at a row, psi_0 from the ground
        // condition and psi_n = 0.
        auto rhs = [&](std::size_t row) {
            if (row == 0)
                return mul(rhs_main_0, psi[0]) + mul(rhs_upper_0, psi[1]);
            if (row + 1 == m)
                return mul(rhs_off, psi[m - 2]) + mul(rhs_main[m - 1], psi[m - 1]);
            return mul(rhs_off, psi[row - 1] + psi[row + 1]) + mul(rhs_main[row], psi[row]);
        };
        const std::size_t down = p;         // rows 0 .. p-1
        const std::size_t up = m - 1 - p;   // rows m-1 .. p+1
        const std::size_t both = std::min(down, up);
        // Eliminate: both chains a row at a time, from the ends inward.  The
        // end rows carry the boundary conditions and go first; inside, the
        // right side is the plain three-point stencil.
        complex top = mul(rhs(0), inv[0]);
        complex bottom = mul(rhs(m - 1), inv[m - 1]);
        w[0] = top;
        w[m - 1] = bottom;
        for (std::size_t k = 1; k < both; k++)
        {
            const std::size_t rt = k, rb = m - 1 - k;
            const complex bt = mul(rhs_off, psi[rt - 1] + psi[rt + 1]) + mul(rhs_main[rt], psi[rt]);
            const complex bb = mul(rhs_off, psi[rb - 1] + psi[rb + 1]) + mul(rhs_main[rb], psi[rb]);
            top = minus_mul(mul(bt, inv[rt]), chain[rt], top);
            bottom = minus_mul(mul(bb, inv[rb]), chain[rb], bottom);
            w[rt] = top;
            w[rb] = bottom;
        }
        for (std::size_t k = both; k < down; k++)
            w[k] = top = minus_mul(mul(rhs(k), inv[k]), chain[k], top);
        for (std::size_t k = both; k < up; k++)
        {
            const std::size_t rb = m - 1 - k;
            w[rb] = bottom = minus_mul(mul(rhs(rb), inv[rb]), chain[rb], bottom);
        }
        // The meeting row, then substitute back outward, again two chains.
        const complex x_p = mul(rhs(p) - mul(lower_p, p > 0 ? w[p - 1] : complex(0.0))
                        - mul(upper_p, p + 1 < m ? w[p + 1] : complex(0.0)), inv[p]);
        psi[p] = x_p;
        top = x_p;
        bottom = x_p;
        for (std::size_t k = 0; k < both; k++)
        {
            const std::size_t rt = p - 1 - k, rb = p + 1 + k;
            top = minus_mul(w[rt], back[rt], top);
            bottom = minus_mul(w[rb], back[rb], bottom);
            psi[rt] = top;
            psi[rb] = bottom;
        }
        for (std::size_t k = both; k < down; k++)
        {
            const std::size_t rt = p - 1 - k;
            psi[rt] = top = minus_mul(w[rt], back[rt], top);
        }
        for (std::size_t k = both; k < up; k++)
        {
            const std::size_t rb = p + 1 + k;
            psi[rb] = bottom = minus_mul(w[rb], back[rb], bottom);
        }
    };
}

// Split-step Pade: psi(r + dr) = [c_0 + sum_j c_j (Q - p_j)^-1] psi(r), the
// propagator of pade_propagator.hpp, so a step is one tridiagonal solve per
// term and may be several wavelengths long.
//
// Each solve is the Thomas algorithm against a matrix factored once.  Written
// for a x' = x, the scaled unknown, the factorization leaves one stored number
// per row, e_r = a / d_r with d_r the pivot: forward w_r = e_r (psi_r -
// w_{r-1}), back x'_r = w_r - e_r x'_{r+1}.  The terms are independent, so the
// loops run term-innermost and the sweeps are as many interleaved recurrences
// as there are terms -- throughput-bound, where Crank-Nicolson's single solve
// waits on its own latency.  Storage is term-innermost to match, two doubles
// of factor and two of w per term per row.
// Two terms' worth of one real quantity.  The compiler would not vectorize the
// march across terms from scalar code -- it came out one term at a time, a
// third as fast -- so the pairing is explicit: the GCC/Clang vector type where
// there is one, and a plain struct that says the same thing where there is not.
#if defined(__GNUC__)

typedef double Lanes __attribute__((vector_size(2 * sizeof(double))));

#else

struct Lanes
{
    double v[2];

    double operator[](std::size_t lane) const
    {
        return v[lane];
    }
    double& operator[](std::size_t lane)
    {
        return v[lane];
    }
};

inline Lanes operator+(Lanes a, Lanes b)
{
    return Lanes {{a.v[0] + b.v[0], a.v[1] + b.v[1]}};
}

inline Lanes operator-(Lanes a, Lanes b)
{
    return Lanes {{a.v[0] - b.v[0], a.v[1] - b.v[1]}};
}

inline Lanes operator*(Lanes a, Lanes b)
{
    return Lanes {{a.v[0] * b.v[0], a.v[1] * b.v[1]}};
}

#endif

inline Lanes both(double value)
{
    return Lanes {value, value};
}

constexpr std::size_t max_pade_pairs = 12;

// The factored solves, two terms to a Lanes; an odd term count is padded with
// a term of zero weight.
struct PadeFactors
{
    std::size_t m = 0;
    std::size_t terms = 0;
    // [(row * pairs + pair) * 2 + part], part 0 real and 1 imaginary.
    std::vector<Lanes> e, w;
    std::vector<complex> first_back; // row 0's back coefficient, upper_0 / d_0
    // Fourth order in height only: each term's weight row by row, laid out as
    // e is, and the right side (1 + delta^2 / 12) psi, which the terms share.
    std::vector<Lanes> row_weight;
    std::vector<complex> smoothed;
    complex g1, g2;
    // The fourth-order ground condition only, per term: the ground row's
    // right side is kappa R short of the others', and psi_0 is rebuilt from
    // the two rows above it and R (pade_march).
    bool fourth_ground = false;
    std::vector<complex> kappa, from_first, from_second;
    complex from_right;
    std::vector<complex> weight; // c_j / a
    complex constant; // c_0
};

// One step, with the number of pairs a compile-time constant so that the
// inner loops have a known, short trip count and unroll.
//
// Both recurrences are written as e - c x, which contracts to two fused
// multiply-adds per part; e (psi - w) would put a subtract on the chain as
// well.  With the terms interleaved the sweeps are no longer waiting on that
// chain, though -- eight terms is sixteen cycles of arithmetic a row against
// sixteen of latency -- which is why the single Crank-Nicolson solve's twisted
// factorization is not repeated here.
template <std::size_t pairs, bool fourth_order>
void pade_step(PadeFactors& factors, std::vector<complex>& psi)
{
    const std::size_t m = factors.m;
    if (fourth_order)
    {
        // b_j = (psi_{j-1} + 10 psi_j + psi_{j+1}) / 12, psi_0 from the ground
        // condition and psi_n = 0.
        // In parts: std::complex's arithmetic would not vectorize here.
        double* __restrict const b = reinterpret_cast<double*>(factors.smoothed.data());
        const double* __restrict const f = reinterpret_cast<const double*>(psi.data());
        const complex ground = factors.fourth_ground ? psi[m] : factors.g1 * psi[0] + factors.g2 * psi[1];
        const double tenth = 10.0 / 12.0, twelfth = 1.0 / 12.0;
        b[0] = tenth * f[0] + twelfth * (ground.real() + f[2]);
        b[1] = tenth * f[1] + twelfth * (ground.imag() + f[3]);
        for (std::size_t k = 2; k + 2 < 2 * m; k++)
            b[k] = tenth * f[k] + twelfth * (f[k - 2] + f[k + 2]);
        b[2 * m - 2] = tenth * f[2 * m - 2] + twelfth * f[2 * m - 4];
        b[2 * m - 1] = tenth * f[2 * m - 1] + twelfth * f[2 * m - 3];
    }
    const complex* __restrict const right = fourth_order ? factors.smoothed.data() : psi.data();
    const Lanes* __restrict const row_weight = factors.row_weight.data();
    const Lanes* __restrict const e = factors.e.data();
    Lanes* __restrict const w = factors.w.data();
    complex* __restrict const field = psi.data();
    // The row in hand for each term, and the terms' weights, on the stack
    // where nothing else can alias them.
    Lanes x_re[pairs], x_im[pairs], weight_re[pairs], weight_im[pairs];
    for (std::size_t pair = 0; pair < pairs; pair++)
    {
        x_re[pair] = both(0.0);
        x_im[pair] = both(0.0);
        for (std::size_t lane = 0; lane < 2; lane++)
        {
            const std::size_t term = 2 * pair + lane;
            weight_re[pair][lane] = term < factors.terms ? factors.weight[term].real() : 0.0;
            weight_im[pair][lane] = term < factors.terms ? factors.weight[term].imag() : 0.0;
        }
    }
    // Forward: w_r = e_r psi_r - e_r w_{r-1}.
    // Row 0 a term at a time: under the fourth-order ground condition each
    // term's right side there is its own.
    const complex ground_before = psi[m];
    const complex ground_right = factors.fourth_ground ? 7.0 * psi[m] + 6.0 * psi[0] - psi[1] : complex(0.0);
    for (std::size_t term = 0; term < factors.terms; term++)
    {
        const std::size_t pair = term / 2, lane = term % 2;
        complex first = right[0];
        if (factors.fourth_ground)
            first -= factors.kappa[term] * ground_right;
        const complex value = complex(e[2 * pair][lane], e[2 * pair + 1][lane]) * first;
        x_re[pair][lane] = value.real();
        x_im[pair][lane] = value.imag();
        w[2 * pair][lane] = value.real();
        w[2 * pair + 1][lane] = value.imag();
    }
    for (std::size_t row = 1; row < m; row++)
    {
        const Lanes p_re = both(right[row].real()), p_im = both(right[row].imag());
        const Lanes* const e_row = e + row * pairs * 2;
        Lanes* const w_row = w + row * pairs * 2;
        for (std::size_t pair = 0; pair < pairs; pair++)
        {
            const Lanes e_re = e_row[2 * pair], e_im = e_row[2 * pair + 1];
            const Lanes t_re = e_re * p_re - e_im * p_im;
            const Lanes t_im = e_re * p_im + e_im * p_re;
            const Lanes re = t_re - e_re * x_re[pair] + e_im * x_im[pair];
            const Lanes im = t_im - e_re * x_im[pair] - e_im * x_re[pair];
            x_re[pair] = re;
            x_im[pair] = im;
            w_row[2 * pair] = re;
            w_row[2 * pair + 1] = im;
        }
    }
    // Back: x_r = w_r - e_r x_{r+1}, summed into psi as it is formed.  The top
    // row's x is its w, which the forward sweep left in x.
    const complex constant = factors.constant;
    auto combine = [&](std::size_t row) {
        Lanes sum_re = both(0.0), sum_im = both(0.0);
        for (std::size_t pair = 0; pair < pairs; pair++)
        {
            const Lanes g_re = fourth_order ? row_weight[(row * pairs + pair) * 2] : weight_re[pair];
            const Lanes g_im = fourth_order ? row_weight[(row * pairs + pair) * 2 + 1] : weight_im[pair];
            sum_re = sum_re + g_re * x_re[pair] - g_im * x_im[pair];
            sum_im = sum_im + g_re * x_im[pair] + g_im * x_re[pair];
        }
        const double p_re = field[row].real(), p_im = field[row].imag();
        field[row] = complex(constant.real() * p_re - constant.imag() * p_im + (sum_re[0] + sum_re[1]),
                constant.real() * p_im + constant.imag() * p_re + (sum_im[0] + sum_im[1]));
    };
    combine(m - 1);
    for (std::size_t row = m - 1; row-- > 1;)
    {
        const Lanes* const e_row = e + row * pairs * 2;
        const Lanes* const w_row = w + row * pairs * 2;
        for (std::size_t pair = 0; pair < pairs; pair++)
        {
            const Lanes e_re = e_row[2 * pair], e_im = e_row[2 * pair + 1];
            const Lanes re = w_row[2 * pair] - e_re * x_re[pair] + e_im * x_im[pair];
            const Lanes im = w_row[2 * pair + 1] - e_re * x_im[pair] - e_im * x_re[pair];
            x_re[pair] = re;
            x_im[pair] = im;
        }
        combine(row);
    }
    // Row 0, whose back coefficient carries the ground; and psi_0 itself,
    // where it is carried.
    complex ground_after = constant * ground_before + factors.from_right * ground_right;
    for (std::size_t term = 0; term < factors.terms; term++)
    {
        const std::size_t pair = term / 2, lane = term % 2;
        const complex x_above(x_re[pair][lane], x_im[pair][lane]);
        const complex x = complex(w[2 * pair][lane], w[2 * pair + 1][lane]) - factors.first_back[term] * x_above;
        x_re[pair][lane] = x.real();
        x_im[pair][lane] = x.imag();
        if (factors.fourth_ground)
            ground_after += factors.from_first[term] * x + factors.from_second[term] * x_above;
    }
    combine(0);
    if (factors.fourth_ground)
        psi[m] = ground_after;
}

using PadeStep = void (*)(PadeFactors&, std::vector<complex>&);

template <bool fourth_order, std::size_t... counts>
constexpr std::array<PadeStep, sizeof...(counts)> pade_steps(std::index_sequence<counts...>)
{
    return {&pade_step<counts + 1, fourth_order>...};
}

// `fourth_order` takes d^2/dz^2 as Numerov does, delta^2 / (1 + delta^2 / 12)
// over dz^2, so that (Q - p) y = psi becomes
//   [a delta^2 + (1 + delta^2 / 12) D] y = (1 + delta^2 / 12) psi,  D = epsilon - p.
// Its off-diagonal entry in column s is c_s = a + D_s / 12, so in y' = c y the
// off-diagonals are one again and the same one-number-a-row factorization
// serves; what it costs is a weight per row, c_0's share over c, where second
// order has one per term.
//
// The ground at fourth order (op.fourth_ground).  Eliminating psi_0 with a
// three-point one-sided derivative is second order, and it was what the
// fourth-order operator converged at.  A one-sided derivative of higher order
// would reach past the tridiagonal, so the equation supplies the accuracy
// instead: with F = y'' = k_a^2 (psi - D y) and y' = lambda y at the ground,
//   y_1 - (1 + dz lambda) y_0 = dz^2 (7 F_0 + 6 F_1 - F_2) / 24 + O(dz^5),
// Taylor's series with F' and F'' taken by differences of F.  It is a row in
// y_0, y_1, y_2 that is linear in the pole, like Numerov's, so the solves are
// still resolvents of one operator -- on the grid with psi_0 in it, which is
// therefore carried (as psi[m]) instead of being derived from psi_1 and psi_2.
// Solving that row for y_0,
//   y_0 = G1 y_1 + G2 y_2 + H R,  R = 7 psi_0 + 6 psi_1 - psi_2,
// puts G1 and G2 where g1 and g2 were and takes c_0 H R off row 1's right
// side, a different amount for each term.
March pade_march(const VerticalOperator& op, const PadePropagator& propagator, bool fourth_order)
{
    const std::size_t m = op.m;
    const std::size_t terms = propagator.poles.size();
    const std::size_t pairs = (terms + 1) / 2;
    if (terms == 0 || pairs > max_pade_pairs)
        throw std::invalid_argument("solve_pe: the propagator has no terms, or too many");
    const double a = op.a;
    auto factors = std::make_shared<PadeFactors>();
    factors->m = m;
    factors->terms = terms;
    factors->e.assign(m * pairs * 2, both(0.0));
    factors->w.assign(m * pairs * 2, both(0.0));
    factors->first_back.resize(terms);
    factors->weight.resize(terms);
    factors->constant = propagator.constant;
    factors->g1 = op.g1;
    factors->g2 = op.g2;
    factors->fourth_ground = fourth_order && op.fourth_ground;
    factors->kappa.assign(terms, 0.0);
    factors->from_first.assign(terms, 0.0);
    factors->from_second.assign(terms, 0.0);
    if (fourth_order)
    {
        factors->row_weight.assign(m * pairs * 2, both(0.0));
        factors->smoothed.resize(m);
    }
    for (std::size_t j = 0; j < terms; j++)
    {
        const complex pole = propagator.poles[j];
        factors->weight[j] = propagator.residues[j] / a;
        // In units where the off-diagonals are one: second order divides the
        // row by a, fourth order by the column's c.  Row 0 carries the ground:
        // its diagonal gains psi_0's g1 share, its upper neighbor's
        // coefficient psi_0's g2 share.
        auto column = [&](std::size_t grid) { return a + (op.epsilon[grid] - pole) / 12.0; };
        complex upper_0, e = 0.0;
        for (std::size_t row = 0; row < m; row++)
        {
            complex d;
            if (fourth_order)
            {
                const complex c = column(row + 1);
                d = (-2.0 * a + (10.0 / 12.0) * (op.epsilon[row + 1] - pole)) / c;
                if (row == 0)
                {
                    complex share_1 = op.g1, share_2 = op.g2; // psi_0's shares of psi_1, psi_2
                    if (factors->fourth_ground)
                    {
                        const complex d0 = op.epsilon[0] - pole, d1 = op.epsilon[1] - pole, d2 = op.epsilon[2] - pole;
                        const complex alpha = -(1.0 + op.ground_slope) + 7.0 * d0 / (24.0 * a);
                        share_1 = -(1.0 + d1 / (4.0 * a)) / alpha;
                        share_2 = d2 / (24.0 * a) / alpha;
                        const complex h = 1.0 / (24.0 * a * alpha);
                        factors->kappa[j] = column(0) * h;
                        factors->from_first[j] = propagator.residues[j] * share_1 / c;
                        factors->from_second[j] = propagator.residues[j] * share_2 / column(2);
                        factors->from_right += propagator.residues[j] * h;
                    }
                    d += column(0) * share_1 / c;
                    upper_0 = 1.0 + column(0) * share_2 / column(2);
                }
                const complex weight = propagator.residues[j] / c;
                factors->row_weight[(row * pairs + j / 2) * 2][j % 2] = weight.real();
                factors->row_weight[(row * pairs + j / 2) * 2 + 1][j % 2] = weight.imag();
            }
            else
            {
                d = -2.0 + (op.epsilon[row + 1] - pole) / a;
                if (row == 0)
                {
                    d += op.g1;
                    upper_0 = 1.0 + op.g2;
                }
            }
            if (row > 0)
                d -= (row == 1 ? upper_0 : complex(1.0)) * e;
            e = 1.0 / d;
            if (!std::isfinite(std::norm(e)))
                throw std::runtime_error("solve_pe: a Pade term's factorization broke down");
            factors->e[(row * pairs + j / 2) * 2][j % 2] = e.real();
            factors->e[(row * pairs + j / 2) * 2 + 1][j % 2] = e.imag();
            if (row == 0)
                factors->first_back[j] = upper_0 * e;
        }
    }
    static constexpr auto second = pade_steps<false>(std::make_index_sequence<max_pade_pairs>());
    static constexpr auto fourth = pade_steps<true>(std::make_index_sequence<max_pade_pairs>());
    const PadeStep step = (fourth_order ? fourth : second)[pairs - 1];
    return [factors, step](std::vector<complex>& psi) { step(*factors, psi); };
}

// Band averages need levels only; single-frequency callers also keep the field.
PeResult solve_pe_impl(const PeSetup& s, bool keep_field)
{
    if (s.z.empty() || s.z.size() != s.c.size())
        throw std::invalid_argument("solve_pe: the sound-speed profile is empty or uneven");
    if (!(s.frequency > 0.0) || !(s.max_range > 0.0) || s.receiver_heights.empty())
        throw std::invalid_argument("solve_pe: frequency, range and receivers are required");

    const double omega = 2.0 * pi * s.frequency;
    const double c_ref = speed_at(s, 0.0);
    const double k_a = omega / c_ref;
    const double wavelength = c_ref / s.frequency;
    // Ten points a wavelength where everything on the grid is fourth order,
    // twenty where anything is second (pe.hpp).
    const bool all_fourth = s.scheme == PeScheme::split_step_pade && s.vertical_order == 4
            && s.boundary_order == 4;
    const double points = s.points_per_wavelength > 0.0 ? s.points_per_wavelength : (all_fourth ? 10.0 : 20.0);
    const double dz = wavelength / points;

    // The domain: everything the source and receivers need, a clear margin,
    // then an absorbing layer that takes what leaves the top.  The margin is
    // generous on purpose.  In a shadow the field near the ground is a weak
    // remainder, while energy leaks upward past it; what the layer sends back
    // down lands on that remainder as a floor.  With 60 wavelengths of margin,
    // 25 m at 1 kHz, the exact residue series (tests/test_pe_shadow.cpp) found
    // the PE 3.8 dB high at -44 dB; with 250 m it agrees to hundredths.
    const double highest = std::max(s.source_height,
            *std::max_element(s.receiver_heights.begin(), s.receiver_heights.end()));
    const double physical_top = s.domain_height > 0.0
            ? s.domain_height
            : highest + std::max({60.0 * wavelength, 0.5 * highest, 250.0});
    const double layer = s.layer_wavelengths * wavelength;
    const double top = physical_top + layer;
    const std::size_t n = static_cast<std::size_t>(std::ceil(top / dz)); // psi_n = 0

    // epsilon_j = (k_j^2 - k_a^2) / k_a^2, complex with the absorption, and
    // with the absorbing layer: a quadratic ramp to Im k ~ k_a / 10 at the top.
    VerticalOperator op;
    std::vector<complex>& epsilon = op.epsilon;
    epsilon.resize(n + 1);
    for (std::size_t j = 0; j <= n; j++)
    {
        const double z = static_cast<double>(j) * dz;
        const complex k(omega / speed_at(s, z), s.absorption);
        double ramp = 0.0;
        if (z > physical_top)
        {
            const double t = (z - physical_top) / layer;
            ramp = 0.2 * t * t;
        }
        epsilon[j] = (k * k - k_a * k_a) / (k_a * k_a) + complex(0.0, ramp);
    }

    // Rows j = 1 .. n-1 (psi_0 from the ground condition, psi_n = 0).
    // Q psi_j = a (psi_{j+1} - 2 psi_j + psi_{j-1}) + epsilon_j psi_j.
    const double a = 1.0 / (k_a * k_a * dz * dz);
    const complex i(0.0, 1.0);
    // Ground: (-3 psi_0 + 4 psi_1 - psi_2) / (2 dz) = -i k beta psi_0, so
    // psi_0 = (4 psi_1 - psi_2) / (3 - 2 i k beta dz) = g1 psi_1 + g2 psi_2.
    const complex denominator = 3.0 - 2.0 * i * k_a * s.admittance * dz;
    const complex g1 = 4.0 / denominator;
    const complex g2 = -1.0 / denominator;
    const std::size_t m = n - 1; // unknowns psi_1 .. psi_{n-1}
    if (m < 3)
        throw std::invalid_argument("solve_pe: the grid is too small");
    op.a = a;
    op.g1 = g1;
    op.g2 = g2;
    op.m = m;
    if (s.boundary_order != 2 && s.boundary_order != 4)
        throw std::invalid_argument("solve_pe: boundary_order is 2 or 4");
    const bool split_step = s.scheme == PeScheme::split_step_pade;
    // Fourth order at the boundaries: the ground condition, which needs the
    // fourth-order operator, and the interpolation onto and off the grid.
    const bool cubic = split_step && s.boundary_order == 4;
    op.fourth_ground = cubic && s.vertical_order == 4;
    op.ground_slope = -i * k_a * s.admittance * dz;
    // Cubic weights for a point u grid spacings up, on the four nodes from
    // `first`: the two either side and one beyond each, or the lowest four.
    auto cubic_weights = [](double u, std::size_t& first, double weight[4]) {
        const std::size_t below = static_cast<std::size_t>(u);
        first = below > 0 ? below - 1 : 0;
        const double t = u - static_cast<double>(first);
        for (int k = 0; k < 4; k++)
        {
            weight[k] = 1.0;
            for (int other = 0; other < 4; other++)
                if (other != k)
                    weight[k] *= (t - other) / static_cast<double>(k - other);
        }
    };

    const bool turbulent = s.gamma_t > 0.0;
    const double spacing = s.screen_spacing > 0.0 ? s.screen_spacing : s.outer_scale;

    const bool fourth_order = s.vertical_order == 4;
    if (s.vertical_order != 2 && s.vertical_order != 4)
        throw std::invalid_argument("solve_pe: vertical_order is 2 or 4");
    PeResult result;
    double dr = 0.0;
    double steepest = 0.0; // tan of the split-step's max_angle; 0 for Crank-Nicolson
    March march;
    March first_step; // the self-starter, when it is used
    if (s.scheme == PeScheme::crank_nicolson)
    {
        dr = wavelength / s.range_points_per_wavelength;
        march = crank_nicolson_march(op, k_a, dr);
    }
    else
    {
        // A whole number of steps per output interval, so results are taken
        // where the march lands and never interpolated: over steps this long
        // psi turns through radians between them.  No step spans more than one
        // phase screen.
        double longest = s.range_step_wavelengths * wavelength;
        if (turbulent)
            longest = std::min(longest, spacing);
        std::size_t per_output = static_cast<std::size_t>(std::ceil(s.output_step / longest - 1e-9));
        per_output = std::max<std::size_t>(per_output, 1);

        // Where the propagator must be accurate: every mode traveling within
        // max_angle of horizontal, anywhere below the absorbing layer.  A wave
        // at elevation theta where the index is k has x = epsilon - (1 +
        // epsilon) sin^2 theta.  Above zero, room for the ground's surface
        // wave, whose horizontal wavenumber exceeds k: x = -beta^2.
        double epsilon_low = epsilon[0].real(), epsilon_high = epsilon[0].real();
        double spectrum_low = epsilon[0].real(), spectrum_imag = 0.0;
        for (std::size_t j = 0; j <= n; j++)
        {
            if (static_cast<double>(j) * dz <= physical_top)
            {
                epsilon_low = std::min(epsilon_low, epsilon[j].real());
                epsilon_high = std::max(epsilon_high, epsilon[j].real());
            }
            spectrum_low = std::min(spectrum_low, epsilon[j].real());
            spectrum_imag = std::max(spectrum_imag, epsilon[j].imag());
        }
        const double cosine = std::cos(s.max_angle * pi / 180.0);
        steepest = std::tan(s.max_angle * pi / 180.0);
        PadeRequest request;
        // Short of x = -1, the propagator's branch point, which no rational
        // function follows: a column whose index falls that far keeps less
        // than max_angle at its top.
        request.accurate_low = std::max((1.0 + epsilon_low) * cosine * cosine - 1.0, -0.75);
        request.accurate_high = epsilon_high + 0.02 + std::max(0.0, -(s.admittance * s.admittance).real());
        // The grid's last mode: -4a at second order, -6a at fourth.
        request.spectrum_low = spectrum_low - (fourth_order ? 6.0 : 4.0) * a - 0.5;
        request.spectrum_imag = spectrum_imag;
        request.tolerance = s.propagator_tolerance;
        PadePropagator propagator;
        for (int attempt = 0;; attempt++)
        {
            dr = s.output_step / static_cast<double>(per_output);
            request.sigma = k_a * dr;
            try
            {
                propagator = pade_propagator(request);
                break;
            }
            catch (const std::runtime_error&)
            {
                // Too long a step for this interval: shorten it.
                if (attempt == 8)
                    throw;
                per_output += (per_output + 1) / 2;
            }
        }
        result.pade_terms = static_cast<int>(propagator.poles.size());
        march = pade_march(op, propagator, fourth_order);
        // The self-starter: the same step with the source's modal amplitudes
        // folded in.  A source within a grid point of the ground has no node
        // below it to share the delta with, and takes the Gaussian.
        if (s.starter == PeStarter::self && s.source_height >= dz)
        {
            // Its fit is harder than the propagator's -- the amplitude factor
            // rises toward the branch point -- and at a long step it can be
            // refused where the propagator's was not.  Then the first step is
            // taken in halves, or quarters: the starter over the first part
            // and that length's own propagator over the rest.
            PadeRequest first = request;
            first.starter = true;
            first.max_terms = static_cast<int>(2 * max_pade_pairs);
            for (const int parts : {1, 2, 4})
            {
                first.sigma = request.sigma / parts;
                try
                {
                    const March start = pade_march(op, pade_propagator(first), fourth_order);
                    March rest;
                    if (parts > 1)
                    {
                        PadeRequest part = request;
                        part.sigma = first.sigma;
                        rest = pade_march(op, pade_propagator(part), fourth_order);
                    }
                    first_step = [start, rest, parts](std::vector<complex>& field) {
                        start(field);
                        for (int k = 1; k < parts; k++)
                            rest(field);
                    };
                    result.self_started = true;
                    break;
                }
                catch (const std::runtime_error&)
                {
                    first_step = nullptr;
                }
            }
        }
    }
    result.range_step = dr;

    // Salomons' second-order Gaussian starter, with its image in the ground
    // at the normal-incidence reflection coefficient.  Normalized so that
    // p = psi exp(i k_a r) / sqrt(r) is a unit point source, exp(ikR)/R.
    // psi_1 .. psi_m, and psi_0 after them where the ground condition carries it.
    std::vector<complex> psi(m + 1);
    const complex reflection = (1.0 - s.admittance) / (1.0 + s.admittance);
    const complex scale = std::sqrt(i * k_a);
    auto starter = [&](double z) {
        const double u = k_a * k_a * z * z;
        return (1.3717 - 0.3701 * u) * std::exp(-u / 3.0);
    };
    for (std::size_t row = 0; row < m && !first_step; row++)
    {
        const double z = static_cast<double>(row + 1) * dz;
        psi[row] = scale * (starter(z - s.source_height) + reflection * starter(z + s.source_height));
    }
    if (!first_step)
        psi[m] = scale * (1.0 + reflection) * starter(s.source_height);
    if (first_step)
    {
        // The self-starter's delta function, shared between the two grid
        // points either side of the source.  psi one step out is
        // sqrt(2 pi i / k_a) (1 + Q)^(-1/4) exp(i k_a dr (sqrt(1 + Q) - 1))
        // delta(z - z_s): the modal sum of i pi H0(k_n r), each Hankel
        // function in its far field.  The Gaussian's area is the same
        // constant.
        const double u = s.source_height / dz;
        const std::size_t below = std::min(static_cast<std::size_t>(u), m - 1); // grid index j >= 1
        const double t = u - static_cast<double>(below);
        const complex weight = std::sqrt(2.0 * pi * i / k_a) / dz;
        std::size_t first;
        double share[4];
        cubic_weights(u, first, share);
        if (cubic && (first > 0 || op.fourth_ground) && first + 3 <= m)
        {
            // Four nodes, so that the source sits where it is to fourth order
            // and not second; node 0 is psi[m].
            for (std::size_t k = 0; k < 4; k++)
                psi[first + k > 0 ? first + k - 1 : m] = weight * share[k];
        }
        else
        {
            psi[below - 1] = weight * (1.0 - t);
            psi[below] = weight * t;
        }
    }

    // Turbulence: a screen every `spacing` meters, below the absorbing layer.
    const std::size_t screened = std::min(m, static_cast<std::size_t>(physical_top / dz));
    std::unique_ptr<ScreenSource> screens;
    if (turbulent)
        screens = std::make_unique<ScreenSource>(k_a, s.gamma_t, s.outer_scale, spacing, dz, screened + 1,
                s.seed);
    double next_screen = 0.5 * spacing;

    result.relative_level.assign(s.receiver_heights.size(), {});
    if (keep_field)
        result.psi.assign(s.receiver_heights.size(), {});
    // Crank-Nicolson's step follows the wavelength; the split-step's divides
    // the output step, up to rounding.
    const std::size_t steps = static_cast<std::size_t>(std::ceil(s.max_range / dr - 1e-9));
    // psi at each receiver, cubic or linear between grid points; psi_0 from
    // the ground condition.
    struct ReceiverStencil
    {
        std::size_t first;
        std::array<double, 4> share {};
    };
    std::vector<ReceiverStencil> receiver_stencils;
    receiver_stencils.reserve(s.receiver_heights.size());
    for (double zr : s.receiver_heights)
    {
        const double u = zr / dz;
        ReceiverStencil stencil;
        if (cubic)
            cubic_weights(u, stencil.first, stencil.share.data());
        else
        {
            stencil.first = std::min(static_cast<std::size_t>(u), m - 1);
            const double t = u - static_cast<double>(stencil.first);
            stencil.share[0] = 1.0 - t;
            stencil.share[1] = t;
        }
        receiver_stencils.push_back(stencil);
    }
    auto at_receiver = [&](std::size_t receiver) {
        const auto& stencil = receiver_stencils[receiver];
        auto at = [&](std::size_t jj) {
            if (jj == 0)
                return op.fourth_ground ? psi[m] : g1 * psi[0] + g2 * psi[1];
            return jj <= m ? psi[jj - 1] : complex(0.0);
        };
        if (cubic)
        {
            complex value = 0.0;
            for (std::size_t k = 0; k < 4; k++)
                value += stencil.share[k] * at(stencil.first + k);
            return value;
        }
        return stencil.share[0] * at(stencil.first) + stencil.share[1] * at(stencil.first + 1);
    };
    // Results at exact multiples of the output step, whatever the range step,
    // so every frequency lands on one grid: psi interpolated between the two
    // range steps either side.  The split-step march lands on them, and the
    // interpolation then takes the later step whole.
    std::vector<complex> previous_value(s.receiver_heights.size());
    for (std::size_t k = 0; k < s.receiver_heights.size(); k++)
        previous_value[k] = at_receiver(k);
    std::size_t next_output = 1;
    std::vector<complex> value(s.receiver_heights.size());
    for (std::size_t step = 1; step <= steps; step++)
    {
        if (step == 1 && first_step)
        {
            first_step(psi);
            first_step = nullptr; // and its factors with it
        }
        else
            march(psi);

        const double r = static_cast<double>(step) * dr;
        // Each screen stands for the slab of turbulence centered on it.
        while (turbulent && r >= next_screen)
        {
            const auto& phi = screens->next();
            for (std::size_t row = 0; row < screened; row++)
                psi[row] *= std::polar(1.0, phi[row + 1]);
            psi[m] *= std::polar(1.0, phi[0]);
            next_screen += spacing;
        }
        for (std::size_t k = 0; k < s.receiver_heights.size(); k++)
            value[k] = at_receiver(k);
        while (next_output * s.output_step <= r + 1e-9 && next_output * s.output_step <= s.max_range + 1e-9)
        {
            const double target = next_output * s.output_step;
            const double t = (target - (r - dr)) / dr;
            result.range.push_back(target);
            for (std::size_t k = 0; k < s.receiver_heights.size(); k++)
            {
                const complex v = (1.0 - t) * previous_value[k] + t * value[k];
                const double zr = s.receiver_heights[k];
                const double slant = std::hypot(target, zr - s.source_height);
                // Steeper than the propagator reaches there is no field to
                // report: it absorbs those modes (pade_propagator.cpp), so the
                // number would be hundreds of dB down and mean nothing.
                const bool too_steep = steepest > 0.0 && std::abs(zr - s.source_height) > steepest * target;
                result.relative_level[k].push_back(too_steep
                        ? std::numeric_limits<double>::quiet_NaN()
                        : 20.0 * std::log10(std::abs(v) * slant / std::sqrt(target) + 1e-300));
                if (keep_field)
                    result.psi[k].push_back(v);
            }
            next_output++;
        }
        previous_value = value;
    }
    return result;
}

} // namespace

PeResult solve_pe(const PeSetup& setup)
{
    return solve_pe_impl(setup, true);
}

PeResult solve_pe_band(PeSetup setup, double center, int samples)
{
    return solve_pe_bands({setup}, {center}, samples)[0];
}

std::vector<PeResult> solve_pe_bands(const std::vector<PeSetup>& setups, const std::vector<double>& centers,
        int samples)
{
    if (setups.size() != centers.size())
        throw std::invalid_argument("solve_pe_bands: one setup per band");
    if (samples < 1)
        throw std::invalid_argument("solve_pe_bands: samples must be positive");
    if (setups.empty())
        return {};
    // Evenly in log frequency across each band's edges, 2^(+-1/6).  Every
    // frequency reports on the same range grid (solve_pe), so the energies add
    // point by point.
    const std::size_t per_band = static_cast<std::size_t>(samples);
    struct Task
    {
        std::size_t band, sample;
        double frequency;
    };
    std::vector<Task> tasks;
    for (std::size_t b = 0; b < centers.size(); b++)
        for (std::size_t k = 0; k < per_band; k++)
        {
            // The midpoint of each of `samples` equal slices of the band's width in
            // log frequency, a third of an octave.  The samples used to sit at
            // k/(N-1), which puts the first and last on the band edges with a full
            // 1/N weight each: 0.16 to 0.47 dB rms off the true band mean downwind
            // (bias +0.1 to +0.35 dB) and 0.4 to 1.5 upwind, against 0.015 to 0.04
            // and 0.1 to 1 with the midpoints.  One sample is the center.
            const double offset = ((static_cast<double>(k) + 0.5) / static_cast<double>(per_band) - 0.5) / 3.0;
            tasks.push_back({b, k, centers[b] * std::pow(2.0, offset)});
        }
    std::stable_sort(tasks.begin(), tasks.end(),
            [](const Task& x, const Task& y) { return x.frequency > y.frequency; });
    std::vector<std::vector<PeResult>> each(centers.size(), std::vector<PeResult>(per_band));
    const unsigned threads = setups[0].threads > 0 ? setups[0].threads : default_thread_count();
    parallel_chunks(tasks.size(), threads, [&](std::size_t lo, std::size_t hi) {
        for (std::size_t t = lo; t < hi; t++)
        {
            PeSetup one = setups[tasks[t].band];
            one.frequency = tasks[t].frequency;
            // With turbulence each frequency is also its own realization.
            one.seed = setups[tasks[t].band].seed + static_cast<unsigned>(tasks[t].sample);
            each[tasks[t].band][tasks[t].sample] = solve_pe_impl(one, false);
        }
    });
    // Summed in frequency order whatever order they ran in, so the result is
    // the same on any thread count.
    std::vector<PeResult> out(centers.size());
    for (std::size_t b = 0; b < centers.size(); b++)
    {
        PeResult& sum = out[b];
        sum.range = each[b][0].range;
        sum.relative_level.assign(each[b][0].relative_level.size(), std::vector<double>(sum.range.size(), 0.0));
        for (const PeResult& one : each[b])
        {
            for (std::size_t rcv = 0; rcv < sum.relative_level.size(); rcv++)
            {
                for (std::size_t q = 0; q < std::min(sum.range.size(), one.range.size()); q++)
                    sum.relative_level[rcv][q] += std::pow(10.0, 0.1 * one.relative_level[rcv][q]) / samples;
            }
        }
        // Release this band before reducing the next, keeping sample order intact.
        each[b].clear();
        for (auto& row : sum.relative_level)
            for (double& v : row)
                v = 10.0 * std::log10(v + 1e-300);
    }
    return out;
}

} // namespace niceops
