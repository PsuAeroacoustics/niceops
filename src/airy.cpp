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
            break; // the series is asymptotic: stop at its smallest term
        last2 = size2;
        sum_u_r += u_k * power.real();
        sum_u_i += u_k * power.imag();
        sum_v_r += v_k * power.real();
        sum_v_i += v_k * power.imag();
        if (last2 < 1e-34)
            break;
    }
    // exp(-zeta) / (2 sqrt(pi)), and z^(1/4) with its reciprocal conj / |.|^2.
    const double magnitude = std::exp(-zeta.real()) / (2.0 * std::sqrt(pi));
    const complex e(magnitude * std::cos(zeta.imag()), -magnitude * std::sin(zeta.imag()));
    const complex quarter = principal_sqrt(root);
    const complex quarter_inverse = std::conj(quarter) / std::norm(quarter);
    ai = product(product(e, quarter_inverse), complex(sum_u_r, sum_u_i));
    aip = -product(product(e, quarter), complex(sum_v_r, sum_v_i));
}

void airy(complex z, complex& ai, complex& aip)
{
    if (std::abs(z) <= 5.0)
    {
        airy_series(z, ai, aip);
        return;
    }
    if (std::abs(std::arg(z)) <= 2.0 * pi / 3.0)
    {
        airy_asymptotic(z, ai, aip);
        return;
    }
    // Near the negative real axis: Ai(z) = -w Ai(w z) - w^2 Ai(w^2 z), and
    // Ai'(z) = -w^2 Ai'(w z) - w Ai'(w^2 z), both arguments inside the sector.
    complex a1, d1, a2, d2;
    airy_asymptotic(w_rot * z, a1, d1);
    airy_asymptotic(w_rot * w_rot * z, a2, d2);
    ai = -w_rot * a1 - w_rot * w_rot * a2;
    aip = -w_rot * w_rot * d1 - w_rot * d2;
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
