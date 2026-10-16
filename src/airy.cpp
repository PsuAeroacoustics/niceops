#include "airy.hpp"

#include <cmath>

#include "util/constants.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

const complex w_rot = std::polar(1.0, 2.0 * pi / 3.0);

} // namespace

void airy_series(complex z, complex& ai, complex& aip)
{
    constexpr double c1 = 0.355028053887817239; // Ai(0)
    constexpr double c2 = 0.258819403792806798; // -Ai'(0)
    const complex z3 = z * z * z;
    const complex z2 = z * z;
    complex f = 1.0, g = z, fp = 0.0, gp = 1.0;
    complex t = 1.0, u = z;
    for (int k = 0; k < 200; k++)
    {
        const double kk = static_cast<double>(k);
        const complex t_next = t * z3 / ((3 * kk + 2) * (3 * kk + 3));
        const complex u_next = u * z3 / ((3 * kk + 3) * (3 * kk + 4));
        f += t_next;
        g += u_next;
        // d/dz of t_next (degree 3k+3) and u_next (degree 3k+4), written
        // without dividing by z.
        fp += t * z2 / (3 * kk + 2);
        gp += u * z2 / (3 * kk + 3);
        t = t_next;
        u = u_next;
        if (std::abs(t) + std::abs(u) < 1e-18 * (std::abs(f) + std::abs(g)))
            break;
    }
    ai = c1 * f - c2 * g;
    aip = c1 * fp - c2 * gp;
}

namespace {

// The principal square root, in real arithmetic (std::sqrt on a complex goes
// through polar form, an atan2 and a sincos).
complex principal_sqrt(complex w)
{
    const double a = w.real(), b = w.imag();
    const double t = std::sqrt(0.5 * (std::hypot(a, b) + std::fabs(a)));
    if (t == 0.0)
        return 0.0;
    return a >= 0.0 ? complex(t, b / (2.0 * t)) : complex(std::fabs(b) / (2.0 * t), std::copysign(t, b));
}

// Complex product without the NaN-recovery branch of the library's.
complex product(complex a, complex b)
{
    return {a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real()};
}

// The asymptotic expansions' terms (DLMF 9.7.5, 9.7.6) in `inverse` = -1/zeta:
// add(k, u_k, v_k, inverse^k) for k from 1, stopping before the terms grow
// (the series is asymptotic: it stops at its smallest term) or once they are
// negligible.
template <class Add>
void asymptotic_terms(complex inverse, const Add& add)
{
    double u_k = 1.0;
    double last2 = 1e300;
    complex power = 1.0;
    for (int k = 1; k < 60; k++)
    {
        const double kk = static_cast<double>(k);
        u_k *= (6 * kk - 5) * (6 * kk - 3) * (6 * kk - 1) / ((2 * kk - 1) * 216.0 * kk);
        const double v_k = -(6 * kk + 1) / (6 * kk - 1) * u_k;
        power = product(power, inverse);
        const double size2 = u_k * u_k * (power.real() * power.real() + power.imag() * power.imag());
        if (size2 > last2)
            break;
        last2 = size2;
        add(k, u_k, v_k, power);
        if (last2 < 1e-34)
            break;
    }
}

} // namespace

void airy_asymptotic(complex z, complex& ai, complex& aip)
{
    // In real arithmetic, as hankel1_0 is: the powers z^(3/2) and z^(1/4) on the
    // principal branch are z sqrt(z) and sqrt(sqrt(z)), the divides are by a
    // reciprocal formed once, and the series stops on squared magnitudes.  The
    // library's complex pow, exp, divide and abs were half the ShadowTransition
    // table build once hankel1_0 was fixed.
    const complex root = principal_sqrt(z);
    const complex zeta = (2.0 / 3.0) * product(z, root);
    const double zeta2 = std::norm(zeta);
    const complex inverse(-zeta.real() / zeta2, zeta.imag() / zeta2); // -1/zeta
    double sum_u_r = 1.0, sum_u_i = 0.0, sum_v_r = 1.0, sum_v_i = 0.0;
    asymptotic_terms(inverse, [&](int, double u_k, double v_k, complex power) {
        sum_u_r += u_k * power.real();
        sum_u_i += u_k * power.imag();
        sum_v_r += v_k * power.real();
        sum_v_i += v_k * power.imag();
    });
    // exp(-zeta) / (2 sqrt(pi)), and z^(1/4) with its reciprocal conj / |.|^2.
    const double magnitude = std::exp(-zeta.real()) / (2.0 * std::sqrt(pi));
    const complex e(magnitude * std::cos(zeta.imag()), -magnitude * std::sin(zeta.imag()));
    const complex quarter = principal_sqrt(root);
    const complex quarter_inverse = std::conj(quarter) / std::norm(quarter);
    ai = product(product(e, quarter_inverse), complex(sum_u_r, sum_u_i));
    aip = -product(product(e, quarter), complex(sum_v_r, sum_v_i));
}

namespace {

// The connection formula's two asymptotic expansions, at w z and w^2 z, for
// z in the sector |arg z| > 2 pi / 3 (where arg(w z) lies in (-2 pi/3, 0)
// and arg(w^2 z) = arg(w z) + 2 pi / 3, both clear of the cut).  They share
// all but signs: zeta(w^2 z) = -zeta(w z), so the series in -1/zeta are the
// same even and odd sums with the odd one's sign flipped, stopping at the
// same term; exp(-zeta) changes the sign of its exponent; and (w^2 z)^(1/4)
// is (w z)^(1/4) turned by pi / 6.  One series, one sincos and two square
// roots, where airy_asymptotic twice takes two of each.
void airy_connected(complex z, complex& ai, complex& aip)
{
    const complex x = product(w_rot, z);
    const complex root = principal_sqrt(x);
    const complex zeta = (2.0 / 3.0) * product(x, root);
    const double zeta2 = std::norm(zeta);
    const complex inverse(-zeta.real() / zeta2, zeta.imag() / zeta2); // -1/zeta
    double even_u_r = 1.0, even_u_i = 0.0, odd_u_r = 0.0, odd_u_i = 0.0;
    double even_v_r = 1.0, even_v_i = 0.0, odd_v_r = 0.0, odd_v_i = 0.0;
    asymptotic_terms(inverse, [&](int k, double u_k, double v_k, complex power) {
        if (k % 2 == 0)
        {
            even_u_r += u_k * power.real();
            even_u_i += u_k * power.imag();
            even_v_r += v_k * power.real();
            even_v_i += v_k * power.imag();
        }
        else
        {
            odd_u_r += u_k * power.real();
            odd_u_i += u_k * power.imag();
            odd_v_r += v_k * power.real();
            odd_v_i += v_k * power.imag();
        }
    });
    const complex u1(even_u_r + odd_u_r, even_u_i + odd_u_i), u2(even_u_r - odd_u_r, even_u_i - odd_u_i);
    const complex v1(even_v_r + odd_v_r, even_v_i + odd_v_i), v2(even_v_r - odd_v_r, even_v_i - odd_v_i);
    const double scale = 1.0 / (2.0 * std::sqrt(pi));
    const double m1 = std::exp(-zeta.real()) * scale, m2 = std::exp(zeta.real()) * scale;
    const double c = std::cos(zeta.imag()), sn = std::sin(zeta.imag());
    const complex e1(m1 * c, -m1 * sn), e2(m2 * c, m2 * sn);
    const complex quarter1 = principal_sqrt(root);
    const complex quarter2 = product(quarter1, complex(0.5 * std::sqrt(3.0), 0.5));
    const double q2 = std::norm(quarter1);
    const complex inverse1 = std::conj(quarter1) / q2, inverse2 = std::conj(quarter2) / q2;
    const complex a1 = product(product(e1, inverse1), u1), d1 = -product(product(e1, quarter1), v1);
    const complex a2 = product(product(e2, inverse2), u2), d2 = -product(product(e2, quarter2), v2);
    const complex w2 = product(w_rot, w_rot);
    ai = -product(w_rot, a1) - product(w2, a2);
    aip = -product(w2, d1) - product(w_rot, d2);
}

} // namespace

void airy(complex z, complex& ai, complex& aip)
{
    // |z| <= 5, and |arg z| <= 2 pi / 3 as |Im z| >= -sqrt(3) Re z: the
    // library's abs and arg would be a hypot and an atan2 on every call.
    if (std::norm(z) <= 25.0)
    {
        airy_series(z, ai, aip);
        return;
    }
    if (std::fabs(z.imag()) >= -std::sqrt(3.0) * z.real())
    {
        airy_asymptotic(z, ai, aip);
        return;
    }
    // Near the negative real axis: Ai(z) = -w Ai(w z) - w^2 Ai(w^2 z), and
    // Ai'(z) = -w^2 Ai'(w z) - w Ai'(w^2 z), both arguments inside the sector.
    airy_connected(z, ai, aip);
}

double airy_prime_zero(int n)
{
    // DLMF 9.9.19, then Newton on Ai' with Ai'' = x Ai.
    const double t = 3.0 * pi * (4.0 * n - 3.0) / 8.0;
    double x = -std::pow(t, 2.0 / 3.0) * (1.0 - 7.0 / (48.0 * t * t));
    for (int i = 0; i < 50; i++)
    {
        complex ai, aip;
        airy(x, ai, aip);
        const double step = aip.real() / (x * ai.real());
        x -= step;
        if (std::fabs(step) < 1e-15 * std::fabs(x))
            break;
    }
    return x;
}

} // namespace niceops
