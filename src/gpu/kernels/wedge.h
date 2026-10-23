// The hard wedge's UTD coefficient (wedge.hpp) and Filon's panel
// (double_wedge.cpp) in float, for the kernels.  Their angles' parts -- each
// term's cotangent and the cosine in its transition function's argument --
// are worked out in double on the host and arrive as floats, which keeps them
// accurate however near a shadow boundary: a cotangent near its pole and a
// cosine near its zero each keep their relative precision.

// Weideman's 32-term rational approximation of the Faddeeva function
// w(z) = exp(-z^2) erfc(-i z), for Im z >= 0 (src/faddeeva.cpp): its
// coefficients, highest power first.
NK_CONSTANT float nk_weideman_a[32] = {
    -1.30317978630500875e-12f, +3.74107755472269789e-12f, +8.03033538228436012e-12f, -2.15435373143760700e-11f,
    -5.54423561143486540e-11f, +1.16582546820233235e-10f, +4.15374271614691897e-10f, -5.23102143724924582e-10f,
    -3.20801514078875819e-09f, +8.12488896568210460e-10f, +2.37975567593842971e-08f, +2.29304391194392872e-08f,
    -1.48130789124396946e-07f, -4.18407637053084035e-07f, +4.25583313792292626e-07f, +4.40153173159092505e-06f,
    +6.82103194402497841e-06f, -2.14096192018698432e-05f, -1.30754492546239731e-04f, -2.45329802700202065e-04f,
    +3.92591360700705167e-04f, +4.51954110534923473e-03f, +1.90061557848454077e-02f, +5.73044035298371918e-02f,
    +1.40607162268937713e-01f, +2.95444510715087261e-01f, +5.46013972063934205e-01f, +9.01925489364799993e-01f,
    +1.34554416923454490e+00f, +1.82566962963248147e+00f, +2.26353729990026764e+00f, +2.57225340812456960e+00f,
};

// w(z) for Im z >= 0: the series in t = (L + i z) / (L - i z), then
// w = (2 p(t) + (L - i z) / sqrt(pi)) / (L - i z)^2.
NK_FN nk_complex nk_faddeeva(nk_complex z)
{
    const float L = 4.75682846001088410f;
    const float dr = L + z.im, di = -z.re;
    const float nr = L - z.im, ni = z.re;
    const float dn = dr * dr + di * di;
    const float tr = (nr * dr + ni * di) / dn, ti = (ni * dr - nr * di) / dn;
    float pr = nk_weideman_a[0], pi = 0.0f;
    for (int j = 1; j < 32; j++)
    {
        const float r = pr * tr - pi * ti + nk_weideman_a[j];
        pi = pr * ti + pi * tr;
        pr = r;
    }
    const float one_over_sqrt_pi = 0.564189583547756287f;
    const float numr = 2.0f * pr + one_over_sqrt_pi * dr, numi = 2.0f * pi + one_over_sqrt_pi * di;
    const float d2r = dr * dr - di * di, d2i = 2.0f * dr * di;
    const float d2n = d2r * d2r + d2i * d2i;
    return nk_make((numr * d2r + numi * d2i) / d2n, (numi * d2r - numr * d2i) / d2n);
}

// The UTD transition function F(x), x >= 0 (utd_transition in wedge.cpp):
// sqrt(pi x) exp(-i pi/4) w(sqrt(x) exp(i pi/4)), and its asymptotic series
// past x = 50.
NK_FN nk_complex nk_utd_transition(float x)
{
    if (x <= 0.0f)
        return nk_make(0.0f, 0.0f);
    if (x > 50.0f)
    {
        const float i1 = 1.0f / x, i2 = i1 * i1;
        return nk_make(1.0f - 0.75f * i2 + 105.0f / 16.0f * i2 * i2,
                -0.5f * i1 + 15.0f / 8.0f * i2 * i1 - 945.0f / 32.0f * i2 * i2 * i1);
    }
    const float root = sqrt(x), half_root2 = 0.707106781186547524f;
    const nk_complex w = nk_faddeeva(nk_make(root * half_root2, root * half_root2));
    // sqrt(pi) root exp(-i pi/4) w
    const float s = 1.772453850905516027f * root * half_root2;
    return nk_make(s * (w.re + w.im), s * (w.im - w.re));
}

// fresnel_tail (wedge.hpp): 1/2 e^(ix) w(sqrt(x) e^(i pi/4)), regular at
// x = 0, and past x = 50 F(x) e^(i (x + pi/4)) / (2 sqrt(pi x)).
NK_FN nk_complex nk_fresnel_tail(float x)
{
    if (!(x > 0.0f))
        return nk_make(0.5f, 0.0f);
    if (x > 50.0f)
        return nk_scale(nk_mul(nk_utd_transition(x), nk_cis(x + 0.785398163397448310f)),
                1.0f / (2.0f * sqrt(3.14159265358979323846f * x)));
    const float root = sqrt(x) * 0.707106781186547524f;
    return nk_scale(nk_mul(nk_faddeeva(nk_make(root, root)), nk_cis(x)), 0.5f);
}

// trust_taper (point_paths.cpp): 1 to lo, 0 from hi, cosine-squared in log x.
NK_FN float nk_trust_taper(float x, float lo, float hi)
{
    if (!(x > lo))
        return 1.0f;
    if (!(x < hi))
        return 0.0f;
    const float c = cos(1.57079632679489662f * log(x / lo) / log(hi / lo));
    return c * c;
}

// The hard wedge's coefficient D (hard_wedge_coefficient(WedgeAngles, k, r)):
// -exp(i pi/4) / (2 n sqrt(2 pi k)) times the sum over its four terms of
// cot_i F(2 k r c_i^2).
NK_FN nk_complex nk_wedge_coefficient(NK_GLOBAL_IN float* cot, NK_GLOBAL_IN float* c, float n, float k, float r)
{
    const float kl = k * r;
    nk_complex sum = nk_make(0.0f, 0.0f);
    for (int i = 0; i < 4; i++)
        sum = nk_add(sum, nk_scale(nk_utd_transition(2.0f * kl * c[i] * c[i]), cot[i]));
    const float scale = -1.0f / (2.0f * n * sqrt(6.283185307179586477f * k));
    // scale exp(i pi/4) sum
    const float h = 0.707106781186547524f * scale;
    return nk_make(h * (sum.re - sum.im), h * (sum.re + sum.im));
}

// ended_wedge_coefficient (point_paths.cpp): D at k_plane and r with its
// edge's ends -- excesses lo and hi, m, -1 for none, past (-1, 0, +1), and
// past an end its scale -- each term losing each end's fresnel_tail weighted
// by 1 / sqrt(1 + k Delta / X); past an end faded 15-30 radians out and taken
// relative to the ray through the near end.
NK_FN nk_complex nk_ended_wedge_coefficient(NK_GLOBAL_IN float* cot, NK_GLOBAL_IN float* c, float n, float k_plane,
        float r, float lo, float hi, float past, float scale, float k)
{
    const bool has_lo = lo >= 0.0f, has_hi = hi >= 0.0f;
    if (!has_lo && !has_hi)
        return nk_wedge_coefficient(cot, c, n, k_plane, r);
    const float klo = has_lo ? k * lo : 0.0f, khi = has_hi ? k * hi : 0.0f;
    const nk_complex tl = has_lo ? nk_fresnel_tail(klo) : nk_make(0.0f, 0.0f);
    const nk_complex th = has_hi ? nk_fresnel_tail(khi) : nk_make(0.0f, 0.0f);
    const bool near_lo = past < 0.0f;
    const float near = near_lo ? klo : khi;
    const nk_complex shift = past == 0.0f ? nk_make(1.0f, 0.0f)
                                          : nk_scale(nk_cis(-near), scale * nk_trust_taper(near, 15.0f, 30.0f));
    const float kl = k_plane * r;
    nk_complex sum = nk_make(0.0f, 0.0f);
    for (int i = 0; i < 4; i++)
    {
        const float x = 2.0f * kl * c[i] * c[i];
        const float wl = klo > 0.0f ? 1.0f / sqrt(1.0f + klo / x) : 1.0f;
        const float wh = khi > 0.0f ? 1.0f / sqrt(1.0f + khi / x) : 1.0f;
        const nk_complex a = nk_scale(tl, wl), b = nk_scale(th, wh);
        nk_complex share;
        if (past == 0.0f)
            share = nk_make(1.0f - a.re - b.re, -a.im - b.im);
        else
            share = nk_mul(near_lo ? nk_make(a.re - b.re, a.im - b.im) : nk_make(b.re - a.re, b.im - a.im), shift);
        sum = nk_add(sum, nk_scale(nk_mul(nk_utd_transition(x), share), cot[i]));
    }
    const float scale_d = -1.0f / (2.0f * n * sqrt(6.283185307179586477f * k_plane));
    const float h = 0.707106781186547524f * scale_d;
    return nk_make(h * (sum.re - sum.im), h * (sum.re + sum.im));
}

// m_j / s^(j+1), j = 0, 1, 2, for |phase| < 1 (filon_moment_series): the sum
// over n of (i phase)^n / (n! (n + j + 1)), its even terms the real part and
// its odd ones the imaginary, Horner in phase^2; fourteen terms leave under
// 1e-12.
NK_CONSTANT float nk_filon_c[3][14] = {
    {1.0f, 0.5f, -0.166666667f, -0.0416666667f, 0.00833333333f, 0.00138888889f, -0.000198412698f,
            -2.48015873e-05f, 2.75573192e-06f, 2.75573192e-07f, -2.50521084e-08f, -2.08767570e-09f,
            1.60590438e-10f, 1.14707456e-11f},
    {0.5f, 0.333333333f, -0.125f, -0.0333333333f, 0.00694444444f, 0.00119047619f, -0.000173611111f,
            -2.20458554e-05f, 2.48015873e-06f, 2.50521084e-07f, -2.29644327e-08f, -1.92708526e-09f,
            1.49119692e-10f, 1.07060292e-11f},
    {0.333333333f, 0.25f, -0.1f, -0.0277777778f, 0.00595238095f, 0.00104166667f, -0.000154320988f,
            -1.98412698e-05f, 2.25468976e-06f, 2.29644327e-07f, -2.11979379e-08f, -1.78943632e-09f,
            1.39178379e-10f, 1.00369024e-11f},
};

NK_FN nk_complex nk_filon_moment(int j, float phase)
{
    const float u = phase * phase;
    float even = nk_filon_c[j][12], odd = nk_filon_c[j][13];
    for (int n = 10; n >= 0; n -= 2)
    {
        even = even * u + nk_filon_c[j][n];
        odd = odd * u + nk_filon_c[j][n + 1];
    }
    return nk_make(even, phase * odd);
}

// Filon's panel (filon_panel in double_wedge.cpp): the integral over the
// panel of g(u) exp(i k u), g the quadratic through its three nodes, taken as
// the panel's start u0 and its other nodes' distances d1 and d2 from it --
// those formed in double before they were rounded, so that a panel of a few
// micrometres beside a cut keeps them.  k and u0 come as float-float pairs,
// for the panel's phase k u0 (nk_phase): tens of metres at 10 kHz would lose
// 6e-4 rad in float.
NK_FN nk_complex nk_filon_panel(float k, float k_lo, float u0, float u0_lo, float d1, float d2, nk_complex g0,
        nk_complex g1, nk_complex g2)
{
    const float span = d2;
    const nk_complex c1 = nk_scale(nk_make(g1.re - g0.re, g1.im - g0.im), 1.0f / d1);
    const nk_complex c2 = nk_scale(
            nk_make((g2.re - g0.re) / span - c1.re, (g2.im - g0.im) / span - c1.im), 1.0f / (d2 - d1));
    const float phase = k * span;
    nk_complex m0, m1, m2;
    if (fabs(phase) < 1.0f)
    {
        m0 = nk_scale(nk_filon_moment(0, phase), span);
        m1 = nk_scale(nk_filon_moment(1, phase), span * span);
        m2 = nk_scale(nk_filon_moment(2, phase), span * span * span);
    }
    else
    {
        // Dividing by i k is multiplying by -i / k.
        const nk_complex e = nk_cis(phase);
        const float by_k = 1.0f / k;
        m0 = nk_make(e.im * by_k, -(e.re - 1.0f) * by_k);
        m1 = nk_make((span * e.im - m0.im) * by_k, -(span * e.re - m0.re) * by_k);
        m2 = nk_make((span * span * e.im - 2.0f * m1.im) * by_k, -(span * span * e.re - 2.0f * m1.re) * by_k);
    }
    // g0 m0 + c1 m1 + c2 (m2 - d1 m1)
    const nk_complex sum = nk_add(nk_add(nk_mul(g0, m0), nk_mul(c1, m1)),
            nk_mul(c2, nk_make(m2.re - d1 * m1.re, m2.im - d1 * m1.im)));
    return nk_mul(nk_cis(nk_phase(k, k_lo, u0, u0_lo)), sum);
}
