// What every kernel program shares, in the dialect (prelude_*.h): complex
// floats and a group-wide sum.

// A complex float as two floats, since the backends' vector types differ.
struct nk_complex
{
    float re;
    float im;
};

NK_FN nk_complex nk_make(float re, float im)
{
    nk_complex z;
    z.re = re;
    z.im = im;
    return z;
}

NK_FN nk_complex nk_add(nk_complex a, nk_complex b)
{
    return nk_make(a.re + b.re, a.im + b.im);
}

NK_FN nk_complex nk_mul(nk_complex a, nk_complex b)
{
    return nk_make(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

NK_FN nk_complex nk_scale(nk_complex a, float s)
{
    return nk_make(a.re * s, a.im * s);
}

NK_FN nk_complex nk_div(nk_complex a, nk_complex b)
{
    const float d = b.re * b.re + b.im * b.im;
    return nk_make((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}

// exp(i x).
NK_FN nk_complex nk_cis(float x)
{
    float c;
    const float s = nk_sincos(x, c);
    return nk_make(c, s);
}

// The sum over a group of each thread's `value`, in every thread; `scratch`
// holds a float per thread, and the group's size is a power of two.
NK_FN float nk_group_sum(NK_SHARED float* scratch, float value, uint thread_index, uint threads)
{
    scratch[thread_index] = value;
    NK_BARRIER();
    for (uint s = threads / 2; s > 0; s >>= 1)
    {
        if (thread_index < s)
            scratch[thread_index] += scratch[thread_index + s];
        NK_BARRIER();
    }
    const float total = scratch[0];
    NK_BARRIER();
    return total;
}

// A phase k L reduced to [-pi, pi], from k and L each as a float-float pair
// (high + low): the product in float-float by an exact two-product, then less
// the nearest multiple of 2 pi (itself a float-float pair).  About 3e-7 rad of
// error for any phase a float's range allows, where k L in float would lose
// 0.03 rad at 4e5.
NK_FN float nk_phase(float k_hi, float k_lo, float l_hi, float l_lo)
{
    const float product = k_hi * l_hi;
    const float error = fma(k_hi, l_hi, -product) + (k_hi * l_lo + k_lo * l_hi);
    const float two_pi_hi = 6.28318548202514648f, two_pi_lo = -1.74845553e-07f;
    const float turns = rint(product * 0.159154943091895336f);
    const float reduced = fma(-turns, two_pi_hi, product);
    return (reduced - turns * two_pi_lo) + error;
}
