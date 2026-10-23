#include "profile_pe.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "exact_solutions.hpp"
#include "pade_propagator.hpp"
#include "util/constants.hpp"
#include "wedge.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

// The grid: twenty points a wavelength in height, the spacing shortened by up
// to a tenth where that puts the roofs on nodes (aligned_spacing).
constexpr double points_per_wavelength = 20.0;
constexpr double most_shortened = 1.1;
// Above the highest of the roofs and the ends, the physical domain reaches
// this far, and an absorbing layer as thick lies on it: each the larger of
// 15 wavelengths and sqrt(8 lambda L), L the march's length.  A layer of a
// fixed number of wavelengths reflects the waves meeting it near grazing,
// whose vertical wavelength lambda / sin(theta) is far longer: with 15, a
// street's parts at 8 kHz were off their converged values by 3 dB, and by
// 0.3 dB at 1 kHz.  What has to be absorbed is what would come back down to
// the roofs within the march, which leaves at sin(theta) ~ reach / L, so the
// layer goes as lambda L / reach and the two together as the root.  Within
// 0.09 dB of a domain 250 wavelengths high on streets, two streets, near and
// far sources and a march from the microphone, 500 Hz to 8 kHz (twice the
// scale would take it to 0.02, for 14% more march at 8 kHz).
constexpr double least_margin_wavelengths = 15.0;
constexpr double margin_scale = 8.0;
// The Huygens line: this far before the last edge, at most half its roof,
// and up to this far above the highest roof or end (within the physical
// domain).  15 was 0.24 dB short marched from a microphone at 2 kHz.
constexpr double huygens_wavelengths = 4.0;
constexpr double line_wavelengths = 30.0;
// The split-step Pade march (pade_propagator.hpp): every mode within this
// angle of horizontal carried to the propagator's tolerance, the steps from
// this set, longest first, each a whole number of times: a roof or a
// street is then covered to an eighth of a wavelength, the rest by one
// Crank-Nicolson step, and the propagators are fitted once for every run.
// Against the Crank-Nicolson march at a tenth of a wavelength it moves the
// parts by 0.04 dB at most; 50 degrees against 60 by 0.03 at most, sources
// to 40 degrees over the first roof, for a sixth fewer terms.
constexpr double max_angle_degrees = 50.0;
constexpr std::array<double, 7> step_wavelengths {8.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.125};
constexpr std::size_t max_poles = 16;

// The 2-D Green's function, (i/4) H0(k r).
complex green(double k, double r)
{
    return complex(0.0, 0.25) * hankel1_0(complex(k * r, 0.0));
}

// A rigid right-angled corner's Green's function by UTD: the corner at c,
// its face 0 straight down (the back wall) and its face n = 3/2 along -u (the
// roof).  The direct wave and each face's reflection where they reach R, and
// the corner's diffraction.
struct CornerGreen
{
    double k, cu, cz, ru, rz;
    static constexpr double n = 1.5;

    // Angles from face 0, counterclockwise: down 0, +u pi/2, up pi, -u 3pi/2.
    double angle(double u, double z) const
    {
        const double a = std::atan2(u - cu, -(z - cz));
        return a < 0.0 ? a + 2.0 * pi : a;
    }
    // R's angle and Q's.  On a shadow or reflection boundary the waves' test
    // and the coefficient's limit must take the same side; R within a hair of
    // one is moved off it, for both, as hard_wedge_field does.  A Huygens
    // node lands on one exactly at some wavenumber (its height is a whole
    // number of steps of a twentieth of a wavelength), and was 30 dB loud.
    std::array<double, 2> angles(double qu, double qz) const
    {
        const double th0 = angle(qu, qz);
        double th = angle(ru, rz);
        for (double boundary : {th0 + pi, th0 - pi, pi - th0, (2.0 * n - 1.0) * pi - th0})
            if (std::abs(th - boundary) < 1e-9)
                th = boundary + 1e-8;
        return {th, th0};
    }
    // Which geometrical-optics waves reach R from Q: direct, off face 0, off
    // face n.
    std::array<bool, 3> lit(double qu, double qz) const
    {
        const auto [th, th0] = angles(qu, qz);
        return {std::abs(th - th0) < pi, th + th0 < pi, th + th0 > (2.0 * n - 1.0) * pi};
    }
    // G(Q, R) with the waves `lit` says reach R from Q.
    complex operator()(double qu, double qz, const std::array<bool, 3>& waves) const
    {
        const auto [th, th0] = angles(qu, qz);
        const double rho0 = std::hypot(qu - cu, qz - cz), rho = std::hypot(ru - cu, rz - cz);
        complex out(0.0, 0.0);
        if (waves[0])
            out += green(k, std::hypot(ru - qu, rz - qz));
        if (waves[1])
            out += green(k, std::sqrt(rho * rho + rho0 * rho0 - 2.0 * rho * rho0 * std::cos(th + th0)));
        if (waves[2])
            out += green(k, std::sqrt(rho * rho + rho0 * rho0 - 2.0 * rho * rho0 * std::cos(2.0 * n * pi - th0 - th)));
        const double L = rho * rho0 / (rho + rho0);
        out += green(k, rho0) * hard_wedge_coefficient(n, th, th0, k, L) * std::polar(1.0, k * rho) / std::sqrt(rho);
        return out;
    }
};

// Complex arithmetic written out: std::complex's multiply carries a branch for
// NaN and infinity that cannot arise here, and it sits in the recurrences.
inline complex mul(complex x, complex y)
{
    return complex(x.real() * y.real() - x.imag() * y.imag(), x.real() * y.imag() + x.imag() * y.real());
}

// e - c x, as two dependent fused multiply-adds per part (pe.cpp's
// minus_mul): the recurrences run at their latency.  Without a fused
// multiply-add in the target std::fma is a library call, so there it is the
// plain products.
inline complex minus_mul(complex e, complex c, complex x)
{
#if defined(__FMA__) || defined(__aarch64__) || defined(__ARM_FEATURE_FMA)
    return complex(std::fma(c.imag(), x.imag(), std::fma(-c.real(), x.real(), e.real())),
            std::fma(-c.imag(), x.real(), std::fma(-c.real(), x.imag(), e.imag())));
#else
    return complex(e.real() - c.real() * x.real() + c.imag() * x.imag(),
            e.imag() - c.real() * x.imag() - c.imag() * x.real());
#endif
}

// The propagators for step_wavelengths, fitted once.  c, the second
// difference's weight 1 / (k dz)^2, sets the grid's most negative eigenvalue,
// -4c, which the fit must not amplify: the largest c the spacing allows.
const std::vector<PadePropagator>& propagators()
{
    static const std::vector<PadePropagator> fitted = [] {
        const double c = std::pow(points_per_wavelength * most_shortened / (2.0 * pi), 2);
        const double cosine = std::cos(max_angle_degrees * pi / 180.0);
        std::vector<PadePropagator> out;
        for (const double steps : step_wavelengths)
        {
            PadeRequest request;
            request.sigma = 2.0 * pi * steps;
            request.accurate_low = std::max(cosine * cosine - 1.0, -0.75);
            request.accurate_high = 0.02;
            request.spectrum_low = -4.0 * c - 0.5;
            request.spectrum_imag = 1.0;
            out.push_back(pade_propagator(request));
            if (out.back().poles.size() > max_poles)
                throw std::logic_error("profile_pe: a propagator has more terms than the march takes");
        }
        return out;
    }();
    return fitted;
}

// Every pole of one propagator, (q - p) y = b on the nodes from a ground
// node g up (rigid there, its ghost the node above) and zero above the top,
// factored from the top down: y_i = w_i - e_i y_{i-1}, w_i = b_i inv_i -
// e_i w_{i+1}.  Rows above g do not depend on it, so one factorization serves
// every roof and street, g's own pivot formed per solve.  Below the layer the
// coefficients settle to the recurrence's fixed point and are that from there
// down, so only the layer and the rows to settle are factored.  Row-major
// over the poles, so that a sweep takes them all a row at a time.
struct Poles
{
    std::size_t count = 0, rows = 0; // rows from the top; the last holds from there down
    std::vector<double> inv_re, inv_im, e_re, e_im; // [row from the top * count + pole]
};

Poles factor_from_top(const std::vector<complex>& poles, std::size_t layer_nodes, std::size_t n, double c)
{
    const std::size_t count = poles.size();
    std::vector<std::vector<complex>> inv(count), e(count);
    std::size_t rows = 1;
    for (std::size_t j = 0; j < count; j++)
    {
        complex e_above = 0.0;
        for (std::size_t k = 0; k < n; k++)
        {
            complex a = 0.0;
            if (k < layer_nodes)
            {
                const double s = static_cast<double>(layer_nodes - k) / static_cast<double>(layer_nodes);
                a = complex(0.0, s * s);
            }
            const complex iv = 1.0 / (-2.0 * c + a - poles[j] - c * e_above);
            const complex ev = c * iv;
            inv[j].push_back(iv);
            e[j].push_back(ev);
            const bool settled = k > layer_nodes && std::abs(ev - e_above) <= 1e-15 * std::abs(ev);
            e_above = ev;
            if (settled)
                break;
        }
        rows = std::max(rows, inv[j].size());
    }
    Poles f;
    f.count = count;
    f.rows = rows;
    f.inv_re.resize(rows * count);
    f.inv_im.resize(rows * count);
    f.e_re.resize(rows * count);
    f.e_im.resize(rows * count);
    for (std::size_t j = 0; j < count; j++)
        for (std::size_t k = 0; k < rows; k++)
        {
            const std::size_t from = std::min(k, inv[j].size() - 1);
            f.inv_re[k * count + j] = inv[j][from].real();
            f.inv_im[k * count + j] = inv[j][from].imag();
            f.e_re[k * count + j] = e[j][from].real();
            f.e_im[k * count + j] = e[j][from].imag();
        }
    return f;
}

// One split step, psi <- c0 psi + sum_j c_j (q - p_j)^-1 psi, for both fields
// from the ground node g up; `work` holds the downward sweep,
// [row][field][real, imaginary][pole].
template <std::size_t P>
void split_step(const PadePropagator& propagator, const Poles& f, const std::vector<complex>& a, double c,
        std::size_t g, std::vector<complex>& f0, std::vector<complex>& f1, std::vector<double>& work)
{
    const std::size_t n = a.size();
    work.resize(n * 4 * P);
    double t[4][P] = {};
    for (std::size_t i = n; i-- > g + 1;)
    {
        const std::size_t r = std::min(n - 1 - i, f.rows - 1) * P;
        const double* __restrict ir = &f.inv_re[r];
        const double* __restrict ii = &f.inv_im[r];
        const double* __restrict er = &f.e_re[r];
        const double* __restrict ei = &f.e_im[r];
        double* __restrict out = &work[i * 4 * P];
        const double ar = f0[i].real(), ai = f0[i].imag(), br = f1[i].real(), bi = f1[i].imag();
        for (std::size_t j = 0; j < P; j++)
        {
            const double x0r = ar * ir[j] - ai * ii[j] - (er[j] * t[0][j] - ei[j] * t[1][j]);
            const double x0i = ar * ii[j] + ai * ir[j] - (er[j] * t[1][j] + ei[j] * t[0][j]);
            const double x1r = br * ir[j] - bi * ii[j] - (er[j] * t[2][j] - ei[j] * t[3][j]);
            const double x1i = br * ii[j] + bi * ir[j] - (er[j] * t[3][j] + ei[j] * t[2][j]);
            t[0][j] = x0r;
            t[1][j] = x0i;
            t[2][j] = x1r;
            t[3][j] = x1i;
            out[j] = x0r;
            out[P + j] = x0i;
            out[2 * P + j] = x1r;
            out[3 * P + j] = x1i;
        }
    }
    // The ground's row, its ghost the row above, and the sums there.
    double y[4][P];
    complex sum0 = mul(propagator.constant, f0[g]), sum1 = mul(propagator.constant, f1[g]);
    const std::size_t above = std::min(n - 2 - g, f.rows - 1) * P;
    for (std::size_t j = 0; j < P; j++)
    {
        const complex e1(f.e_re[above + j], f.e_im[above + j]);
        const complex inv = 1.0 / (-2.0 * c + a[g] - propagator.poles[j] - 2.0 * c * e1);
        const complex y0 = mul(f0[g] - 2.0 * c * complex(t[0][j], t[1][j]), inv);
        const complex y1 = mul(f1[g] - 2.0 * c * complex(t[2][j], t[3][j]), inv);
        y[0][j] = y0.real();
        y[1][j] = y0.imag();
        y[2][j] = y1.real();
        y[3][j] = y1.imag();
        sum0 += mul(propagator.residues[j], y0);
        sum1 += mul(propagator.residues[j], y1);
    }
    f0[g] = sum0;
    f1[g] = sum1;
    double cr[P], ci[P];
    for (std::size_t j = 0; j < P; j++)
    {
        cr[j] = propagator.residues[j].real();
        ci[j] = propagator.residues[j].imag();
    }
    const double c0r = propagator.constant.real(), c0i = propagator.constant.imag();
    // Back up, each pole's solution added into the step's result.
    for (std::size_t i = g + 1; i < n; i++)
    {
        const std::size_t r = std::min(n - 1 - i, f.rows - 1) * P;
        const double* __restrict er = &f.e_re[r];
        const double* __restrict ei = &f.e_im[r];
        const double* __restrict in = &work[i * 4 * P];
        double s0r = c0r * f0[i].real() - c0i * f0[i].imag(), s0i = c0r * f0[i].imag() + c0i * f0[i].real();
        double s1r = c0r * f1[i].real() - c0i * f1[i].imag(), s1i = c0r * f1[i].imag() + c0i * f1[i].real();
        for (std::size_t j = 0; j < P; j++)
        {
            const double n0r = in[j] - (er[j] * y[0][j] - ei[j] * y[1][j]);
            const double n0i = in[P + j] - (er[j] * y[1][j] + ei[j] * y[0][j]);
            const double n1r = in[2 * P + j] - (er[j] * y[2][j] - ei[j] * y[3][j]);
            const double n1i = in[3 * P + j] - (er[j] * y[3][j] + ei[j] * y[2][j]);
            y[0][j] = n0r;
            y[1][j] = n0i;
            y[2][j] = n1r;
            y[3][j] = n1i;
            s0r += cr[j] * n0r - ci[j] * n0i;
            s0i += cr[j] * n0i + ci[j] * n0r;
            s1r += cr[j] * n1r - ci[j] * n1i;
            s1i += cr[j] * n1i + ci[j] * n1r;
        }
        f0[i] = complex(s0r, s0i);
        f1[i] = complex(s1r, s1i);
    }
}

// split_step for each number of poles, the count a template parameter so the
// loops over the poles unroll and vectorize.
using SplitStep = void (*)(const PadePropagator&, const Poles&, const std::vector<complex>&, double, std::size_t,
        std::vector<complex>&, std::vector<complex>&, std::vector<double>&);

template <std::size_t... counts>
constexpr std::array<SplitStep, sizeof...(counts)> split_steps(std::index_sequence<counts...>)
{
    return {&split_step<counts + 1>...};
}

// One Crank-Nicolson step, (1 + alpha_l q) psi' = (1 + alpha_r q) psi, on the
// nodes from g up: the rest of a roof or a street under an eighth of a
// wavelength.  The pivot is folded into the stored coefficients, so each
// recurrence is one complex multiply-add a row, and both fields go together,
// two independent recurrences in each sweep.
class CrankNicolson
{
public:
    CrankNicolson(const std::vector<complex>& a, std::size_t g, double c, complex alpha_l, complex alpha_r)
        : g_(g), off_(alpha_r * c)
    {
        const std::size_t m = a.size() - g;
        inverse_.resize(m);
        chain_.resize(m);
        back_.resize(m);
        main_.resize(m);
        const complex off = alpha_l * c;
        complex previous = 0.0; // upper_{j-1} / pivot_{j-1}
        for (std::size_t j = 0; j < m; j++)
        {
            const complex lower = j > 0 ? off : complex(0.0);
            const complex upper = j + 1 < m ? (j == 0 ? 2.0 * off : off) : complex(0.0);
            const complex inverse = 1.0 / (1.0 + alpha_l * (-2.0 * c + a[g + j]) - lower * previous);
            inverse_[j] = inverse;
            chain_[j] = lower * inverse;
            back_[j] = upper * inverse;
            previous = back_[j];
            main_[j] = 1.0 + alpha_r * (-2.0 * c + a[g + j]);
        }
    }
    void operator()(std::vector<complex>& f0, std::vector<complex>& f1)
    {
        const std::size_t m = inverse_.size();
        w0_.resize(m);
        w1_.resize(m);
        const complex* p0 = f0.data() + g_;
        const complex* p1 = f1.data() + g_;
        // The right side's row j: main psi_j + off (psi_{j-1} + psi_{j+1}),
        // the ghost below row 0 being row 1 and the node above the top zero.
        complex t0 = mul(mul(main_[0], p0[0]) + mul(2.0 * off_, p0[1]), inverse_[0]);
        complex t1 = mul(mul(main_[0], p1[0]) + mul(2.0 * off_, p1[1]), inverse_[0]);
        w0_[0] = t0;
        w1_[0] = t1;
        for (std::size_t j = 1; j < m; j++)
        {
            const complex next0 = j + 1 < m ? p0[j + 1] : complex(0.0), next1 = j + 1 < m ? p1[j + 1] : complex(0.0);
            const complex b0 = mul(main_[j], p0[j]) + mul(off_, p0[j - 1] + next0);
            const complex b1 = mul(main_[j], p1[j]) + mul(off_, p1[j - 1] + next1);
            w0_[j] = t0 = minus_mul(mul(b0, inverse_[j]), chain_[j], t0);
            w1_[j] = t1 = minus_mul(mul(b1, inverse_[j]), chain_[j], t1);
        }
        complex* q0 = f0.data() + g_;
        complex* q1 = f1.data() + g_;
        q0[m - 1] = t0;
        q1[m - 1] = t1;
        for (std::size_t j = m - 1; j-- > 0;)
        {
            q0[j] = t0 = minus_mul(w0_[j], back_[j], t0);
            q1[j] = t1 = minus_mul(w1_[j], back_[j], t1);
        }
    }

private:
    std::size_t g_;
    complex off_;
    std::vector<complex> inverse_, chain_, back_, main_, w0_, w1_;
};

// The system (1 + alpha q) x = rhs on the nodes from g up, rigid at g: the
// Huygens line's derivative.
class Factored
{
public:
    Factored(const std::vector<complex>& a, std::size_t g, double c, complex alpha)
    {
        const std::size_t m = a.size() - g;
        upper_.assign(m, alpha * c);
        upper_[0] = alpha * 2.0 * c;
        upper_[m - 1] = 0.0;
        multiplier_.resize(m);
        inverse_.resize(m);
        complex pivot = 1.0 + alpha * (-2.0 * c + a[g]);
        inverse_[0] = 1.0 / pivot;
        multiplier_[0] = 0.0;
        for (std::size_t j = 1; j < m; j++)
        {
            multiplier_[j] = alpha * c / pivot;
            pivot = 1.0 + alpha * (-2.0 * c + a[g + j]) - multiplier_[j] * upper_[j - 1];
            inverse_[j] = 1.0 / pivot;
        }
    }
    // x in place of rhs.
    void solve(std::vector<complex>& rhs) const
    {
        const std::size_t m = inverse_.size();
        for (std::size_t j = 1; j < m; j++)
            rhs[j] -= multiplier_[j] * rhs[j - 1];
        rhs[m - 1] *= inverse_[m - 1];
        for (std::size_t j = m - 1; j-- > 0;)
            rhs[j] = (rhs[j] - upper_[j] * rhs[j + 1]) * inverse_[j];
    }

private:
    std::vector<complex> upper_, multiplier_, inverse_;
};

// The vertical spacing: lambda / 20, or up to a tenth less where that puts
// the roofs on nodes.  A roof between nodes is taken at the nearer, and the
// level stepped up and down by up to 0.2 dB as the frequency moved it from
// one side to the other.  Every spacing that puts some roof exactly on a node
// is tried, and the one leaving the others nearest theirs kept, the coarser
// of equals.
double aligned_spacing(const std::vector<ProfileBlock>& blocks, double dz)
{
    std::vector<double> heights;
    for (const ProfileBlock& b : blocks)
        if (std::none_of(heights.begin(), heights.end(), [&](double h) { return h == b.h; }))
            heights.push_back(b.h);
    auto off = [&](double spacing) {
        double worst = 0.0;
        for (double h : heights)
            worst = std::max(worst, std::abs(h / spacing - std::round(h / spacing)));
        return worst;
    };
    double best = dz, best_off = off(dz);
    for (double h : heights)
        for (double count = std::ceil(h / dz); count <= std::floor(h * most_shortened / dz); count++)
        {
            const double spacing = h / count;
            const double o = off(spacing);
            if (o < best_off - 1e-9 || (o < best_off + 1e-9 && spacing > best))
            {
                best = spacing;
                best_off = o;
            }
        }
    return best;
}

} // namespace

bool profile_applies(const ProfileProblem& p)
{
    if (p.blocks.empty())
        return false;
    for (std::size_t i = 0; i < p.blocks.size(); i++)
    {
        const ProfileBlock& b = p.blocks[i];
        if (!(b.u1 > b.u0) || !(b.h > 0.0) || (i > 0 && b.u0 < p.blocks[i - 1].u1))
            return false;
    }
    return p.su < p.blocks.front().u0 && p.ru > p.blocks.back().u1 && p.sz >= 0.0 && p.rz >= 0.0;
}

ProfileParts profile_parts(const ProfileProblem& p, double k)
{
    if (!profile_applies(p))
        throw std::invalid_argument("profile_parts: the blocks must be ascending and apart, the source before them "
                                    "and the receiver behind");
    const double lambda = 2.0 * pi / k;
    const double dz = aligned_spacing(p.blocks, lambda / points_per_wavelength);
    const double c = 1.0 / ((k * dz) * (k * dz));
    const ProfileBlock& last = p.blocks.back();
    const double u_start = p.blocks.front().u0;
    const double huygens = std::min(huygens_wavelengths * lambda, 0.5 * (last.u1 - last.u0));
    const double u_stop = last.u1 - huygens;
    // The column holds what the Huygens line needs -- the roofs, the receiver
    // and where the line from the source to it crosses -- and, at the start,
    // the source's wave that comes down to that within the march: below its
    // ray to the top of it, where that crosses the first front face.  Not the
    // source's own height, which the start does not need: a source 122 m up
    // took a march over a terrace 160 m across to 11.8 s one way.
    double end = std::max(p.rz, p.sz + (p.rz - p.sz) * (u_stop - p.su) / (p.ru - p.su));
    for (const ProfileBlock& b : p.blocks)
        end = std::max(end, b.h);
    const double top = std::max(end, p.sz + (end - p.sz) * (u_start - p.su) / (u_stop - p.su));
    // The physical domain's top on a node and the layer a whole number of
    // nodes over it (factor_from_top's rows).
    const double margin
            = std::max(least_margin_wavelengths * lambda, std::sqrt(margin_scale * lambda * (u_stop - u_start)));
    const std::size_t physical = static_cast<std::size_t>(std::ceil((top + margin) / dz));
    const std::size_t layer_nodes = std::max<std::size_t>(1, static_cast<std::size_t>(std::lround(margin / dz)));
    const std::size_t n = physical + layer_nodes + 1;
    const double z_phys = static_cast<double>(physical) * dz;
    std::vector<double> z(n);
    std::vector<complex> a(n, complex(0.0, 0.0));
    for (std::size_t i = 0; i < n; i++)
    {
        z[i] = static_cast<double>(i) * dz;
        if (i > physical)
        {
            const double s = static_cast<double>(i - physical) / static_cast<double>(layer_nodes);
            a[i] = complex(0.0, s * s);
        }
    }
    auto node = [&](double h) { return static_cast<std::size_t>(std::lround(h / dz)); };

    // The starts: the source's wave on the first block's front face, above
    // its roof, and its floor image's.
    std::vector<complex> psi[2] = {std::vector<complex>(n, 0.0), std::vector<complex>(n, 0.0)};
    const std::size_t g0 = node(p.blocks.front().h);
    const complex back = std::polar(1.0, -k * u_start);
    for (std::size_t i = g0; i < n; i++)
    {
        psi[0][i] = green(k, std::hypot(u_start - p.su, z[i] - p.sz)) * back;
        psi[1][i] = green(k, std::hypot(u_start - p.su, z[i] + p.sz)) * back;
    }

    // The march: over each roof, and each street between, to u_stop.
    struct Segment
    {
        double length, h;
    };
    std::vector<Segment> segments;
    for (std::size_t b = 0; b < p.blocks.size(); b++)
    {
        const ProfileBlock& block = p.blocks[b];
        if (b + 1 == p.blocks.size())
            segments.push_back({u_stop - block.u0, block.h});
        else
        {
            segments.push_back({block.u1 - block.u0, block.h});
            if (p.blocks[b + 1].u0 > block.u1)
                segments.push_back({p.blocks[b + 1].u0 - block.u1, 0.0});
        }
    }
    static constexpr auto steps = split_steps(std::make_index_sequence<max_poles>());
    const std::vector<PadePropagator>& propagator = propagators();
    std::vector<Poles> factored(step_wavelengths.size());
    std::vector<double> work;
    for (const Segment& segment : segments)
    {
        // Below the ground the field is gone (a wall stepped up into) or not
        // yet there (a street's floor stepped down to).
        const std::size_t g = node(segment.h);
        for (std::vector<complex>& f : psi)
            std::fill(f.begin(), f.begin() + static_cast<std::ptrdiff_t>(g), complex(0.0, 0.0));
        double left = segment.length / lambda;
        for (std::size_t level = 0; level < step_wavelengths.size(); level++)
        {
            const double length = step_wavelengths[level];
            // The longest as many times as it fits, each shorter once at most.
            std::size_t times = level == 0 ? static_cast<std::size_t>(std::floor(left / length + 1e-12))
                                           : (left >= length - 1e-12 ? 1 : 0);
            if (times == 0)
                continue;
            if (factored[level].count == 0)
                factored[level] = factor_from_top(propagator[level].poles, layer_nodes, n, c);
            for (; times > 0; times--)
            {
                steps[factored[level].count - 1](propagator[level], factored[level], a, c, g, psi[0], psi[1], work);
                left -= length;
            }
        }
        if (left > 1e-9)
        {
            const double du = left * lambda;
            CrankNicolson rest(a, g, c, (1.0 - complex(0.0, k * du)) / 4.0, (1.0 + complex(0.0, k * du)) / 4.0);
            rest(psi[0], psi[1]);
        }
    }
    const double u = u_stop;

    // The Huygens line: p and dp/du = i k Q psi exp(i k u), Q = (1 + 3q/4) / (1 + q/4),
    // to the receiver and to its image.
    const std::size_t g = node(last.h);
    const std::size_t m = n - g;
    const Factored quarter(a, g, c, complex(0.25, 0.0));
    std::vector<complex> dpsi[2];
    for (int f = 0; f < 2; f++)
    {
        dpsi[f].resize(m);
        for (std::size_t j = 0; j < m; j++)
        {
            const std::size_t i = g + j;
            const complex below = j == 0 ? psi[f][i + 1] : psi[f][i - 1];
            const complex above = j + 1 < m ? psi[f][i + 1] : complex(0.0, 0.0);
            dpsi[f][j] = psi[f][i] + 0.75 * (c * (below - 2.0 * psi[f][i] + above) + a[i] * psi[f][i]);
        }
        quarter.solve(dpsi[f]);
    }
    const complex phase = std::polar(1.0, k * u);
    const double dl = 1e-3 * lambda;
    ProfileParts out {};
    const double extent = std::min(z_phys, top + line_wavelengths * lambda) - last.h, t0 = 0.5 * extent;
    for (int side = 0; side < 2; side++)
    {
        const CornerGreen green_w {k, last.u1, last.h, p.ru, side == 0 ? p.rz : -p.rz};
        for (std::size_t j = 0; j < m; j++)
        {
            const double zj = z[g + j], t = zj - last.h;
            if (t > extent)
                break;
            double w = j == 0 ? 0.5 * dz : dz;
            if (t > t0)
                w *= 0.5 * (1.0 + std::cos(pi * (t - t0) / (extent - t0)));
            // dG/du by a central difference, or a one-sided one where a
            // shadow boundary passes between the neighbours: there the
            // diffraction cancels the wave's step only to O(1 / k rho), and
            // the remainder over 2 dl was a spike in the receiver's level.
            const std::array<bool, 3> waves = green_w.lit(u, zj), plus = green_w.lit(u + dl, zj),
                                      minus = green_w.lit(u - dl, zj);
            const complex gq = green_w(u, zj, waves);
            complex dgdu;
            if (plus == waves && minus == waves)
                dgdu = (green_w(u + dl, zj, plus) - green_w(u - dl, zj, minus)) / (2.0 * dl);
            else if (plus == waves)
                dgdu = (green_w(u + dl, zj, plus) - gq) / dl;
            else
                dgdu = (gq - green_w(u - dl, zj, minus)) / dl;
            for (int f = 0; f < 2; f++)
            {
                const complex pu = psi[f][g + j] * phase;
                const complex dpu = complex(0.0, k) * dpsi[f][j] * phase;
                out.part[f][side] += w * (-gq * dpu + pu * dgdu);
            }
        }
    }
    // From (i/4) H0(k r) ~ exp(i k r) / sqrt(r) times sqrt(1 / (8 pi k)) exp(i pi / 4).
    const complex scale = std::sqrt(8.0 * pi * k) * std::polar(1.0, -0.25 * pi);
    for (auto& row : out.part)
        for (complex& v : row)
            v *= scale;
    return out;
}

complex profile_field(const ProfileProblem& p, double k)
{
    const ProfileParts parts = profile_parts(p, k);
    return parts.part[0][0] + p.floor * (parts.part[1][0] + parts.part[0][1]) + p.floor * p.floor * parts.part[1][1];
}

std::vector<ProfileParts> profile_spectrum(const ProfileProblem& problem, const std::vector<double>& wavenumbers)
{
    std::vector<ProfileParts> out;
    out.reserve(wavenumbers.size());
    for (double k : wavenumbers)
        out.push_back(profile_parts(problem, k));
    return out;
}

double profile_path(const ProfileProblem& p, bool source_image, bool receiver_image,
        std::vector<std::array<double, 2>>* points)
{
    // The upper hull of the ends and the roofs' corners.
    std::vector<std::array<double, 2>> hull {{p.su, source_image ? -p.sz : p.sz}};
    auto add = [&](std::array<double, 2> q) {
        while (hull.size() >= 2)
        {
            const auto& o = hull[hull.size() - 2];
            const auto& b = hull.back();
            if ((b[0] - o[0]) * (q[1] - o[1]) - (b[1] - o[1]) * (q[0] - o[0]) < 0.0)
                break;
            hull.pop_back();
        }
        hull.push_back(q);
    };
    for (const ProfileBlock& b : p.blocks)
    {
        add({b.u0, b.h});
        add({b.u1, b.h});
    }
    add({p.ru, receiver_image ? -p.rz : p.rz});
    double length = 0.0;
    for (std::size_t i = 1; i < hull.size(); i++)
        length += std::hypot(hull[i][0] - hull[i - 1][0], hull[i][1] - hull[i - 1][1]);
    if (points)
        *points = std::move(hull);
    return length;
}

} // namespace niceops
