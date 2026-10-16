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

// |c| as sqrt(re^2 + im^2) where the squares can neither overflow nor
// underflow, and std::abs (a hypot) only outside that: the library's is a
// call on every term of the mode sums.
double magnitude(complex c)
{
    const double n2 = c.real() * c.real() + c.imag() * c.imag();
    return n2 > 1e-290 && n2 < 1e290 ? std::sqrt(n2) : std::abs(c);
}

// a / b in real arithmetic: the library's divide is a call (__divdc3) that
// rescales to recover infinities and NaN, which cannot arise here.
complex quotient(complex a, complex b)
{
    const double n2 = b.real() * b.real() + b.imag() * b.imag();
    return {(a.real() * b.real() + a.imag() * b.imag()) / n2, (a.imag() * b.real() - a.real() * b.imag()) / n2};
}

// The principal square root in real arithmetic (std::sqrt on a complex goes
// through polar form, an atan2 and a sincos).
complex principal_root(complex w)
{
    const double a = w.real(), b = w.imag();
    const double t = std::sqrt(0.5 * (magnitude(w) + std::fabs(a)));
    if (t == 0.0)
        return 0.0;
    return a >= 0.0 ? complex(t, b / (2.0 * t)) : complex(std::fabs(b) / (2.0 * t), std::copysign(t, b));
}

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
        const complex step = quotient(mode_equation(t, q, derivative), derivative);
        t -= step;
        if (!std::isfinite(t.real()) || !std::isfinite(t.imag()) || magnitude(t) > 1e7)
            return false;
        if (magnitude(step) <= 1e-14 * std::max(1.0, magnitude(t)))
            return true;
    }
    return false;
}

// Past this |a'_n| a zero of Ai' has moved by less than a quarter of the gap
// to its neighbors, pi / sqrt|a'|, so first-order perturbation finds it:
// the shift is -q / a', and |q| / |a'| < pi / (4 sqrt|a'|) once |a'| > 2 |q|^2.
double settled(complex q)
{
    return 4.0 * std::norm(q) + 10.0;
}

// Follows a zero as q moves from `from` to `to`, halving the step whenever
// Newton lands farther than a fraction of the local spacing away -- the sign
// of having jumped to a neighbor.
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

// Past |q| = 6, where following zeros out from q = 0 loses its way (below),
// but only where there is no surface wave (Re q <= 0): each zero seeded from
// the mode equation's large-argument form and polished by Newton.  With
// b = -x and zeta = (2/3) x^(3/2), Ai'(b) + q Ai(b) = 0 becomes
// tan(zeta + pi/4) = sqrt(x) / q, so zeta = (n - 3/4) pi - atan(q / sqrt(x)):
// the n-th zero of Ai' as q -> 0 and, Re q < 0 taking the principal atan to
// -pi/2, the n-th zero of Ai as |q| grows -- the one continuation keeps.
// Within |q| <= 6 these agree with the followed zeros to 4e-14 over the first
// forty (tests/test_shadow_attenuation.cpp).
bool seeds_reach(complex q)
{
    return q.real() <= 0.0;
}

complex seeded_zero(complex q, int n)
{
    complex x = std::pow(1.5 * (n - 0.75) * pi, 2.0 / 3.0);
    for (int i = 0; i < 100; i++)
    {
        const complex next = std::pow(1.5 * ((n - 0.75) * pi - std::atan(q / std::sqrt(x))), 2.0 / 3.0);
        const bool settled_here = std::abs(next - x) <= 1e-13 * std::abs(x);
        x = next;
        if (settled_here)
            break;
    }
    complex t = -x;
    if (!newton(t, q))
        throw std::runtime_error("impedance_airy_zeros: Newton did not settle on zero "
                + std::to_string(n) + " at |q| = " + std::to_string(std::abs(q)));
    return t;
}

// Consecutive zeros step by pi in zeta; a step far from it means one was
// skipped or found twice.
void check_step(complex previous, complex t, int n)
{
    const double step = (2.0 / 3.0 * (std::pow(-t, 1.5) - std::pow(-previous, 1.5))).real() / pi;
    if (!(step > 0.5 && step < 1.5))
        throw std::runtime_error("impedance_airy_zeros: zeros " + std::to_string(n - 1) + " and "
                + std::to_string(n) + " are " + std::to_string(step) + " pi apart in zeta, not one");
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
    if (std::abs(q) > 6.0 && !seeds_reach(q))
        throw std::invalid_argument("impedance_airy_zeros: |q| = " + std::to_string(std::abs(q))
                + " with Re q > 0 is past the 6 the root finder is built for");
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
    const complex kappa = principal_root(kappa_squared);
    return kappa.imag() < 0.0 ? -kappa : kappa;
}

// One range's sum of `count` modes, term(n) the n-th: NaN where the terms
// cancel by more than nine orders, or where the last sixteen have not died
// away.  The mode count serves the nearest range; a farther one's terms die
// sooner, since H0(kappa r) falls as exp(-Im kappa r).  Once a range's terms
// have stayed under 1e-17 of its largest sixteen in a row -- the rule that
// set the count -- the rest up to the last sixteen are skipped: each is
// smaller still, so neither the largest term nor the sum moves beyond
// rounding.  The last sixteen are always taken, for the tail test.
template <class Term>
complex range_sum(std::size_t count, const Term& term)
{
    const std::size_t tail_from = count > 16 ? count - 16 : 0;
    complex sum = 0.0;
    double biggest = 0.0, tail = 0.0;
    int quiet = 0;
    for (std::size_t n = 0; n < count; n++)
    {
        const complex t = term(n);
        sum += t;
        const double size = magnitude(t);
        biggest = std::max(biggest, size);
        if (n >= tail_from)
            tail = std::max(tail, size);
        quiet = size < 1e-17 * biggest ? quiet + 1 : 0;
        if (quiet >= 16 && n + 1 < tail_from)
            n = tail_from - 1;
    }
    const double size = std::abs(sum);
    const bool trusted = biggest <= 1e9 * size && tail <= 1e-9 * size;
    return trusted ? sum : complex(std::numeric_limits<double>::quiet_NaN(), 0.0);
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
        const double size = magnitude(m.weight * hankel1_0(m.kappa * nearest));
        largest = std::max(largest, size);
        quiet = size < 1e-17 * largest ? quiet + 1 : 0;
    };
    if (std::abs(q) <= 6.0)
    {
        for (const complex t : checked_zeros(q, ordinary))
            add(t);
        for (int n = ordinary + 1; (quiet < 16 || modes.size() < 20) && modes.size() < max_modes; n++)
            add(ordinary_zero(q, n));
    }
    else
    {
        check_reach(q);
        complex previous = 0.0;
        for (int n = 1; (quiet < 16 || modes.size() < 20) && modes.size() < max_modes; n++)
        {
            const complex t = seeded_zero(q, n);
            if (n > 1)
                check_step(previous, t, n);
            add(t);
            previous = t;
        }
    }

    std::vector<complex> field(ranges.size());
    for (std::size_t k = 0; k < ranges.size(); k++)
        field[k] = range_sum(modes.size(), [&](std::size_t n) {
            return modes[n].weight * hankel1_0(modes[n].kappa * ranges[k]);
        });
    return field;
}

// modal_sum for several weight sets sharing the modes' wavenumbers -- one
// source height each -- so the zeros are formed once.  Modes are taken until
// every set's terms at the nearest range have died away.  Each range is
// summed for its own set only, set_of[k], since a range belongs to one height
// and H0(kappa r) is all the sets would share there.
std::vector<complex> modal_sum_heights(complex q,
        const std::function<complex(complex, std::vector<complex>&)>& mode, std::size_t sets,
        const std::vector<double>& ranges, const std::vector<std::size_t>& set_of)
{
    std::vector<complex> field(ranges.size());
    if (ranges.empty() || sets == 0)
        return field;
    if (set_of.size() != ranges.size())
        throw std::invalid_argument("modal_sum: one set per range");
    for (const std::size_t m : set_of)
        if (m >= sets)
            throw std::invalid_argument("modal_sum: a range's set is out of range");
    const double nearest = *std::min_element(ranges.begin(), ranges.end());
    if (!(nearest > 0.0))
        throw std::invalid_argument("modal_sum: ranges must be positive");
    constexpr std::size_t max_modes = 12000;
    int ordinary = 0;
    std::vector<complex> kappas;
    std::vector<std::vector<complex>> weights;
    std::vector<double> largest(sets, 0.0);
    std::vector<int> quiet(sets, 0);
    std::vector<complex> w(sets);
    auto all_quiet = [&]() {
        for (const int qq : quiet)
            if (qq < 16)
                return false;
        return true;
    };
    auto add = [&](complex t) {
        const complex kappa = mode(t, w);
        kappas.push_back(kappa);
        weights.push_back(w);
        const double h = magnitude(hankel1_0(kappa * nearest));
        for (std::size_t m = 0; m < sets; m++)
        {
            const double size = magnitude(w[m]) * h;
            largest[m] = std::max(largest[m], size);
            quiet[m] = size < 1e-17 * largest[m] ? quiet[m] + 1 : 0;
        }
    };
    if (std::abs(q) <= 6.0)
    {
        for (const complex t : checked_zeros(q, ordinary))
            add(t);
        for (int n = ordinary + 1; (!all_quiet() || kappas.size() < 20) && kappas.size() < max_modes; n++)
            add(ordinary_zero(q, n));
    }
    else
    {
        check_reach(q);
        complex previous = 0.0;
        for (int n = 1; (!all_quiet() || kappas.size() < 20) && kappas.size() < max_modes; n++)
        {
            const complex t = seeded_zero(q, n);
            if (n > 1)
                check_step(previous, t, n);
            add(t);
            previous = t;
        }
    }
    for (std::size_t k = 0; k < ranges.size(); k++)
    {
        const std::size_t m = set_of[k];
        field[k] = range_sum(kappas.size(), [&](std::size_t n) {
            return weights[n][m] * hankel1_0(kappas[n] * ranges[k]);
        });
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
    if (std::abs(q) > 6.0)
    {
        std::vector<complex> zeros;
        for (int n = 1; n <= count; n++)
        {
            zeros.push_back(seeded_zero(q, n));
            if (n > 1)
                check_step(zeros[n - 2], zeros[n - 1], n);
        }
        return zeros;
    }
    int ordinary = 0;
    auto checked = checked_zeros(q, ordinary);
    std::vector<complex> zeros(checked.begin(), checked.begin() + std::min(ordinary, count));
    for (int n = ordinary + 1; n <= count; n++)
        zeros.push_back(ordinary_zero(q, n));
    zeros.insert(zeros.end(), checked.begin() + ordinary, checked.end());
    return zeros;
}

std::vector<complex> first_impedance_airy_zero(complex direction, const std::vector<double>& s)
{
    std::vector<complex> zeros;
    zeros.reserve(s.size());
    complex t = airy_prime_zero(1);
    double previous = 0.0;
    for (const double next : s)
    {
        if (!(next >= previous))
            throw std::invalid_argument("first_impedance_airy_zero: s must be non-negative and increasing");
        // Steps of at most a tenth of s itself (0.1 near q = 0), so that Newton
        // starts near the zero it follows; follow() halves them if not.
        const double reach = std::max(0.1, 0.1 * previous);
        const int steps = std::max(1, static_cast<int>(std::ceil((next - previous) / reach)));
        for (int j = 1; j <= steps; j++)
        {
            const double from = previous + (next - previous) * (j - 1) / steps;
            const double to = previous + (next - previous) * j / steps;
            follow(t, from * direction, to * direction);
        }
        zeros.push_back(t);
        previous = next;
    }
    return zeros;
}

complex hankel1_0(complex z)
{
    // The series loses digits as exp(|z| + Im z) -- J0 grows where H0 decays
    // -- and the expansion's smallest term is exp(-2|z|); they cross near
    // |z| = 12 on the real axis and 9 on the imaginary one.
    const double size = magnitude(z);
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
    // Its large-argument expansion (DLMF 10.17.5), to the smallest term.  In
    // real arithmetic: every std::complex multiply and divide carries a branch
    // to recover NaN and infinity that cannot arise here, and the divide is a
    // call (__divdc3) that this loop made up to 60 times a call -- 56% of the
    // ShadowTransition table build.  1/z is formed once, and the step from one
    // term to the next is i f_k / z with f_k = -(2k-1)^2 / 8k real.
    const double zr = z.real(), zi = z.imag();
    const double n2 = zr * zr + zi * zi;
    const double inv_r = zr / n2, inv_i = -zi / n2; // 1/z
    double sum_r = 1.0, sum_i = 0.0, term_r = 1.0, term_i = 0.0;
    double last2 = 1.0;
    for (int k = 1; k < 60; k++)
    {
        const double odd = 2.0 * k - 1.0;
        const double f = -odd * odd / (8.0 * k);
        // i f / z = f (zi, zr) / |z|^2, as (-inv_i, inv_r) f.
        const double c_r = -f * inv_i, c_i = f * inv_r;
        const double next_r = term_r * c_r - term_i * c_i;
        const double next_i = term_r * c_i + term_i * c_r;
        term_r = next_r;
        term_i = next_i;
        const double size2 = term_r * term_r + term_i * term_i;
        if (size2 > last2)
            break;
        last2 = size2;
        sum_r += term_r;
        sum_i += term_i;
        if (last2 < 1e-34)
            break;
    }
    // sqrt(2 / (pi z)), the principal root, then exp(i (z - pi/4)).
    const complex root = principal_root(complex((2.0 / pi) * inv_r, (2.0 / pi) * inv_i));
    const double root_r = root.real(), root_i = root.imag();
    const double decay = std::exp(-zi);
    const double phase = zr - pi / 4.0;
    const double e_r = decay * std::cos(phase), e_i = decay * std::sin(phase);
    const double p_r = root_r * e_r - root_i * e_i;
    const double p_i = root_r * e_i + root_i * e_r;
    return {p_r * sum_r - p_i * sum_i, p_r * sum_i + p_i * sum_r};
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

namespace {

// An upward-refracting problem's terms: the impedance's q, and what a mode
// needs, formed once so that per mode there are only products: w / l (a
// height's shift in the Airy argument), l / w (the norm's scale), 1 / (w l^2)
// (the wavenumber's slope in the zero) and k0^2.
struct UpwardTerms
{
    complex q, shift, scale, curvature, k0_squared;

    UpwardTerms(const ExactProblem& p, double radius)
    {
        const complex k0 = wavenumber(p);
        const complex ell = std::pow(radius / (2.0 * k0 * k0), 1.0 / 3.0);
        q = -I * k0 * ell * p.admittance * w_rot * w_rot;
        shift = quotient(w_rot, ell);
        scale = quotient(ell, w_rot);
        curvature = quotient(1.0, w_rot * ell * ell);
        k0_squared = k0 * k0;
    }
};

} // namespace

std::vector<complex> upward_pressure(const ExactProblem& p, double radius,
        const std::vector<double>& ranges)
{
    const UpwardTerms u(p, radius);
    auto mode = [&](complex b) {
        complex a0, d0, us, ur, unused;
        airy(b, a0, d0);
        airy(b - u.shift * p.source_height, us, unused);
        airy(b - u.shift * p.receiver_height, ur, unused);
        const complex norm = u.scale * (b * a0 * a0 - d0 * d0);
        return Mode {upper(u.k0_squared + b * u.curvature), quotient(pi * I * us * ur, norm)};
    };
    return modal_sum(u.q, mode, ranges);
}

std::vector<complex> upward_pressure_heights(const ExactProblem& p, double radius,
        const std::vector<double>& source_heights, const std::vector<double>& ranges,
        const std::vector<std::size_t>& height_of)
{
    const UpwardTerms u(p, radius);
    // A source on the ground takes Ai(b) itself.
    auto mode = [&](complex b, std::vector<complex>& weights) {
        complex a0, d0, us, ur, unused;
        airy(b, a0, d0);
        airy(b - u.shift * p.receiver_height, ur, unused);
        const complex norm = u.scale * (b * a0 * a0 - d0 * d0);
        const complex common = quotient(pi * I * ur, norm);
        for (std::size_t m = 0; m < source_heights.size(); m++)
        {
            if (source_heights[m] == 0.0)
                us = a0;
            else
                airy(b - u.shift * source_heights[m], us, unused);
            weights[m] = common * us;
        }
        return upper(u.k0_squared + b * u.curvature);
    };
    return modal_sum_heights(u.q, mode, source_heights.size(), ranges, height_of);
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
