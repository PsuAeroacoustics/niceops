#include "exact_solutions.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

#include "airy.hpp"
#include "ground_effect.hpp"
#include "util/constants.hpp"

// The mode sums, derived here rather than copied, because one of the
// published forms is wrong (below).
//
// In cylindrical coordinates the field of a unit point source is
//   p(r, z) = int_{-inf}^{inf} G(kappa; z, z_s) H0(kappa r) kappa dkappa,
//   G'' + (k^2(z) - kappa^2) G = -delta(z - z_s),
//   G' = -i k0 beta G at z = 0 (the PE's condition, and the paper's eq. 40),
// and closing the contour in the upper half plane leaves one residue per
// eigenvalue kappa_n^2 of the height problem.  With eigenfunctions u_n and
// N_n = int_0^inf u_n^2 dz (not |u_n|^2: the problem is complex symmetric),
//   p = pi i sum_n u_n(z) u_n(z_s) H0(kappa_n r) / N_n.
// Take l = (R / 2 k0^2)^(1/3) and w = exp(2 pi i / 3).
//
// Downward, k^2 = k0^2 (1 - 2z/R): the height equation is Airy's in t + z/l,
// so u = Ai(t + z/l), decaying upward, kappa^2 = k0^2 + t / l^2.  The ground
// condition is Ai'(t) + i k0 l beta Ai(t) = 0, and
// N = l (Ai'(t)^2 - t Ai(t)^2).  This is their eq. (10) up to its overall sign.
//
// Upward, k^2 = k0^2 (1 + 2z/R): Airy's in (kappa^2 - k0^2) l^2 - z/l, and the
// solution that radiates upward is u = Ai(b - w z / l), b = w (kappa^2 - k0^2) l^2.
// Then u'(0) = -(w / l) Ai'(b), so the ground condition is
//   Ai'(b) - i k0 l beta w^2 Ai(b) = 0,
// and N = (l / w)(b Ai(b)^2 - Ai'(b)^2).  Their eq. (11) matches this except in
// the ground condition, which it writes Ai'(b) + i k0 l beta w Ai(b) = 0: a
// factor exp(i pi / 3) in the admittance term.  The PE sides with the form
// here (tests/test_benchmark.cpp), and with beta = 0 the two coincide, which is
// why the rigid case never showed it.

namespace niceops {

namespace {

using complex = std::complex<double>;

const complex I(0.0, 1.0);
const complex w_rot = std::polar(1.0, 2.0 * pi / 3.0);

// Ai'(t) + q Ai(t), and its derivative, from Ai'' = t Ai.
complex mode_equation(complex t, complex q, complex& derivative)
{
    complex ai, aip;
    airy(t, ai, aip);
    derivative = t * ai + q * aip;
    return aip + q * ai;
}

bool newton(complex& t, complex q)
{
    for (int i = 0; i < 60; i++)
    {
        complex derivative;
        const complex step = mode_equation(t, q, derivative) / derivative;
        t -= step;
        if (!std::isfinite(t.real()) || !std::isfinite(t.imag()) || std::abs(t) > 1e7)
            return false;
        if (std::abs(step) <= 1e-14 * std::max(1.0, std::abs(t)))
            return true;
    }
    return false;
}

// Past this |a'_n| a zero of Ai' has moved by less than a quarter of the gap
// to its neighbours, pi / sqrt|a'|, so first-order perturbation finds it:
// the shift is -q / a', and |q| / |a'| < pi / (4 sqrt|a'|) once |a'| > 2 |q|^2.
double settled(complex q)
{
    return 4.0 * std::norm(q) + 10.0;
}

// Follows a zero as q moves from `from` to `to`, halving the step whenever
// Newton lands farther than a fraction of the local spacing away -- the sign
// of having jumped to a neighbour.
void follow(complex& t, complex from, complex to, int depth = 0)
{
    complex trial = t;
    const double spacing = pi / std::sqrt(std::max(1.0, std::abs(t)));
    if (newton(trial, to) && std::abs(trial - t) < 0.3 * spacing)
    {
        t = trial;
        return;
    }
    if (depth > 30)
        throw std::runtime_error("impedance_airy_zeros: lost a zero while following it");
    const complex middle = 0.5 * (from + to);
    follow(t, from, middle, depth + 1);
    follow(t, middle, to, depth + 1);
}

complex ordinary_zero(complex q, int n)
{
    const double a = airy_prime_zero(n);
    if (-a > settled(q))
    {
        complex t = a - q / a;
        if (newton(t, q))
            return t;
    }
    complex t = a;
    const int steps = 32;
    for (int j = 1; j <= steps; j++)
        follow(t, q * (static_cast<double>(j - 1) / steps), q * (static_cast<double>(j) / steps));
    return t;
}

// Zeros of Ai' + q Ai inside the box [x0, x1] x [-h, h], by the argument
// principle: the equation is entire, so its argument turns once round the
// edge for every zero inside.
int zeros_inside(complex q, double x0, double x1, double h)
{
    const complex corners[] = {{x0, -h}, {x1, -h}, {x1, h}, {x0, h}, {x0, -h}};
    double turned = 0.0;
    for (int side = 0; side < 4; side++)
    {
        const complex a = corners[side];
        const complex b = corners[side + 1];
        const double length = std::abs(b - a);
        complex unused;
        complex previous = mode_equation(a, q, unused);
        double s = 0.0, step = 0.05;
        while (s < length)
        {
            const double next = std::min(length, s + step);
            const complex value = mode_equation(a + (b - a) * (next / length), q, unused);
            const double change = std::arg(value / previous);
            if (std::fabs(change) > 0.3 && step > 1e-9)
            {
                step *= 0.5;
                continue;
            }
            turned += change;
            previous = value;
            s = next;
            step = std::min(0.5, step * 1.5);
        }
    }
    return static_cast<int>(std::lround(turned / (2.0 * pi)));
}

// The zeros up to the first settled one, and any extra, checked by the
// argument principle.  `ordinary` receives how many continue zeros of Ai'.
// The benchmark's grounds come to |q| = 4.9.  Much past 6 the zeros go over to
// those of Ai rather than Ai' across a band thousands wide, where following
// them from q = 0 loses its way, and the counting box's edge reaches values
// of Ai past the range of a double -- both measured, at |q| = 8.5 and 10.8.
void check_reach(complex q)
{
    if (std::abs(q) > 6.0)
        throw std::invalid_argument("impedance_airy_zeros: |q| = " + std::to_string(std::abs(q))
                + " is past the 6 the root finder is built for");
}

std::vector<complex> checked_zeros(complex q, int& ordinary)
{
    check_reach(q);
    int m = 1;
    while (-airy_prime_zero(m) <= settled(q))
        m++;
    std::vector<complex> zeros;
    for (int n = 1; n <= m; n++)
        zeros.push_back(ordinary_zero(q, n));
    ordinary = m;
    for (std::size_t i = 1; i < zeros.size(); i++)
        for (std::size_t j = 0; j < i; j++)
            if (std::abs(zeros[i] - zeros[j]) < 1e-6 * std::max(1.0, std::abs(zeros[i])))
                throw std::runtime_error("impedance_airy_zeros: two zeros of Ai' led to one zero");

    // A zero with no counterpart among the zeros of Ai' sits where
    // Ai'/Ai ~ -sqrt(t) meets -q, near t = q^2, when Re q > 0.
    if (q.real() > 0.0)
    {
        complex t = q * q;
        bool distinct = newton(t, q);
        for (const complex z : zeros)
            distinct = distinct && std::abs(t - z) > 1e-6 * std::max(1.0, std::abs(t));
        if (distinct)
            zeros.push_back(t);
    }

    // Count them.  The left edge falls midway between the m-th and next zero
    // of Ai', which the settled zeros cannot cross; the rest of the box is
    // wide enough to take in everything found.
    const double x0 = 0.5 * (airy_prime_zero(m) + airy_prime_zero(m + 1));
    double x1 = std::max(6.0, (q * q).real() + 6.0);
    double h = std::max(6.0, std::fabs((q * q).imag()) + 6.0);
    int found = 0;
    for (const complex z : zeros)
    {
        if (z.real() <= x0)
            continue;
        x1 = std::max(x1, z.real() + 4.0);
        h = std::max(h, std::fabs(z.imag()) + 4.0);
        found++;
    }
    const int counted = zeros_inside(q, x0, x1, h);
    if (counted != found)
        throw std::runtime_error("impedance_airy_zeros: the argument principle counts "
                + std::to_string(counted) + " zeros where " + std::to_string(found) + " were found");
    return zeros;
}

struct Mode
{
    complex kappa; // horizontal wavenumber, Im >= 0
    complex weight; // p += weight H0(kappa r)
};

complex upper(complex kappa_squared)
{
    const complex kappa = std::sqrt(kappa_squared);
    return kappa.imag() < 0.0 ? -kappa : kappa;
}

// Sums modes made from the zeros of Ai' + q Ai, taking as many as the
// nearest range needs: until the terms there have fallen 17 orders below the
// largest, sixteen in a row.
std::vector<complex> modal_sum(complex q, const std::function<Mode(complex)>& mode,
        const std::vector<double>& ranges)
{
    if (ranges.empty())
        return {};
    const double nearest = *std::min_element(ranges.begin(), ranges.end());
    if (!(nearest > 0.0))
        throw std::invalid_argument("modal_sum: ranges must be positive");

    constexpr std::size_t max_modes = 12000;
    int ordinary = 0;
    std::vector<Mode> modes;
    double largest = 0.0;
    int quiet = 0;
    auto add = [&](complex t) {
        const Mode m = mode(t);
        modes.push_back(m);
        const double size = std::abs(m.weight * hankel1_0(m.kappa * nearest));
        largest = std::max(largest, size);
        quiet = size < 1e-17 * largest ? quiet + 1 : 0;
    };
    for (const complex t : checked_zeros(q, ordinary))
        add(t);
    for (int n = ordinary + 1; (quiet < 16 || modes.size() < 20) && modes.size() < max_modes; n++)
        add(ordinary_zero(q, n));

    // Each range on its own: NaN where the terms cancel by more than nine
    // orders, or where the last few have not died away.
    std::vector<complex> field(ranges.size());
    const std::size_t tail_from = modes.size() > 16 ? modes.size() - 16 : 0;
    for (std::size_t k = 0; k < ranges.size(); k++)
    {
        complex sum = 0.0;
        double biggest = 0.0, tail = 0.0;
        for (std::size_t n = 0; n < modes.size(); n++)
        {
            const complex term = modes[n].weight * hankel1_0(modes[n].kappa * ranges[k]);
            sum += term;
            biggest = std::max(biggest, std::abs(term));
            if (n >= tail_from)
                tail = std::max(tail, std::abs(term));
        }
        const double size = std::abs(sum);
        const bool trusted = biggest <= 1e9 * size && tail <= 1e-9 * size;
        field[k] = trusted ? sum : complex(std::numeric_limits<double>::quiet_NaN(), 0.0);
    }
    return field;
}

complex wavenumber(const ExactProblem& p)
{
    return {2.0 * pi * p.frequency / p.sound_speed, p.absorption};
}

} // namespace

std::vector<complex> impedance_airy_zeros(complex q, int count)
{
    check_reach(q);
    int ordinary = 0;
    auto checked = checked_zeros(q, ordinary);
    std::vector<complex> zeros(checked.begin(), checked.begin() + std::min(ordinary, count));
    for (int n = ordinary + 1; n <= count; n++)
        zeros.push_back(ordinary_zero(q, n));
    zeros.insert(zeros.end(), checked.begin() + ordinary, checked.end());
    return zeros;
}

complex hankel1_0(complex z)
{
    // The series loses digits as exp(|z| + Im z) -- J0 grows where H0 decays
    // -- and the expansion's smallest term is exp(-2|z|); they cross near
    // |z| = 12 on the real axis and 9 on the imaginary one.
    const double size = std::abs(z);
    if (size < 12.0 - 3.0 * z.imag() / std::max(size, 1e-300))
    {
        // J0 and Y0 by their power series (DLMF 10.2.2, 10.8.2).
        constexpr double euler_gamma = 0.577215664901532861;
        const complex quarter = -z * z / 4.0;
        complex term = 1.0, j0 = 1.0, series = 0.0;
        double harmonic = 0.0;
        for (int k = 1; k < 200; k++)
        {
            term *= quarter / static_cast<double>(k * k);
            harmonic += 1.0 / k;
            j0 += term;
            series += harmonic * term;
            if (std::abs(term) * harmonic < 1e-18)
                break;
        }
        const complex y0 = (2.0 / pi) * ((std::log(z / 2.0) + euler_gamma) * j0 - series);
        return j0 + I * y0;
    }
    // Its large-argument expansion (DLMF 10.17.5), to the smallest term.
    complex sum = 1.0, term = 1.0;
    double last = 1.0;
    for (int k = 1; k < 60; k++)
    {
        const double odd = 2.0 * k - 1.0;
        term *= I * (-odd * odd) / (8.0 * k * z);
        if (std::abs(term) > last)
            break;
        last = std::abs(term);
        sum += term;
        if (last < 1e-17)
            break;
    }
    return std::sqrt(2.0 / (pi * z)) * std::exp(I * (z - pi / 4.0)) * sum;
}

complex homogeneous_pressure(const ExactProblem& p, double range)
{
    const double z_s = p.source_height, z_r = p.receiver_height;
    const double r1 = std::hypot(range, z_r - z_s);
    const double r2 = std::hypot(range, z_r + z_s);
    complex reflection = 1.0;
    if (p.admittance != complex(0.0))
    {
        const auto q = spherical_reflection_coefficient((z_s + z_r) / r2,
                pi * p.frequency * r2 / p.sound_speed, p.admittance, true);
        reflection = complex(q[0], q[1]);
    }
    const complex k = wavenumber(p);
    return std::exp(I * k * r1) / r1 + reflection * std::exp(I * k * r2) / r2;
}

std::vector<complex> upward_pressure(const ExactProblem& p, double radius,
        const std::vector<double>& ranges)
{
    const complex k0 = wavenumber(p);
    const complex ell = std::pow(radius / (2.0 * k0 * k0), 1.0 / 3.0);
    const complex q = -I * k0 * ell * p.admittance * w_rot * w_rot;
    auto mode = [&](complex b) {
        complex a0, d0, us, ur, unused;
        airy(b, a0, d0);
        airy(b - w_rot * p.source_height / ell, us, unused);
        airy(b - w_rot * p.receiver_height / ell, ur, unused);
        const complex norm = ell / w_rot * (b * a0 * a0 - d0 * d0);
        return Mode {upper(k0 * k0 + b / (w_rot * ell * ell)), pi * I * us * ur / norm};
    };
    return modal_sum(q, mode, ranges);
}

std::vector<complex> downward_pressure(const ExactProblem& p, double radius,
        const std::vector<double>& ranges, double max_elevation)
{
    const complex k0 = wavenumber(p);
    const complex ell = std::pow(radius / (2.0 * k0 * k0), 1.0 / 3.0);
    const complex gamma = I * k0 * ell * p.admittance;
    // A mode is the family of rays leaving the ground at cos(elevation) =
    // kappa / k0; evanescent ones leave straight up.  A raised cosine over the
    // last five degrees keeps the cut from ringing.
    constexpr double taper = 5.0;
    auto kept = [&](complex kappa) {
        if (max_elevation >= 90.0)
            return 1.0;
        const double cosine = kappa.real() / k0.real();
        const double elevation = cosine >= 1.0 ? 0.0 : std::acos(std::max(cosine, 0.0)) * 180.0 / pi;
        if (elevation <= max_elevation - taper)
            return 1.0;
        if (elevation >= max_elevation)
            return 0.0;
        return 0.5 * (1.0 + std::cos(pi * (elevation - max_elevation + taper) / taper));
    };
    auto mode = [&](complex t) {
        complex a0, d0, us, ur, unused;
        airy(t, a0, d0);
        airy(t + p.source_height / ell, us, unused);
        airy(t + p.receiver_height / ell, ur, unused);
        const complex norm = ell * (d0 * d0 - t * a0 * a0);
        const complex kappa = upper(k0 * k0 + t / (ell * ell));
        return Mode {kappa, kept(kappa) * pi * I * us * ur / norm};
    };
    return modal_sum(gamma, mode, ranges);
}

} // namespace niceops
