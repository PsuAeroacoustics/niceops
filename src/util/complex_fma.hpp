#pragma once

// Complex arithmetic written out, for the PE's recurrences (pe.cpp's
// Crank-Nicolson sweep, profile_pe.cpp's).  std::complex's multiply carries a
// branch for NaN and infinity that cannot arise there, and its divide is a call
// to __divdc3; both sit on latency-bound chains.
//
// exact_solutions.cpp's quotient() is deliberately left out.  It is the
// textbook form with no fused multiply-add, not matched to __divdc3, and the
// mode sums have carried its bits since it replaced the library's divide;
// divide() below is __divdc3's to the bit on Apple arm64.  Merging the two
// would move one set of results or the other.

#include <cmath>
#include <complex>

// Where the target has a fused multiply-add instruction.  Without one (x86-64
// without -mfma, CI's GCC build) std::fma is a library call, so the helpers
// fall back to the plain products there.
#if defined(__FMA__) || defined(__aarch64__) || defined(__ARM_FEATURE_FMA)
#define NICEOPS_HAS_FMA 1
#else
#define NICEOPS_HAS_FMA 0
#endif

namespace niceops::complex_fma {

using complex = std::complex<double>;

inline complex mul(complex x, complex y)
{
    return complex(x.real() * y.real() - x.imag() * y.imag(), x.real() * y.imag() + x.imag() * y.real());
}

// e - c x, as two dependent fused multiply-adds per part: the recurrences run
// at their latency, and the compiler otherwise leaves a subtract at the end of
// each.
inline complex minus_mul(complex e, complex c, complex x)
{
#if NICEOPS_HAS_FMA
    return complex(std::fma(c.imag(), x.imag(), std::fma(-c.real(), x.real(), e.real())),
            std::fma(-c.imag(), x.real(), std::fma(-c.real(), x.imag(), e.imag())));
#else
    return complex(e.real() - c.real() * x.real() + c.imag() * x.imag(),
            e.imag() - c.real() * x.imag() - c.imag() * x.real());
#endif
}

// w / z.  __divdc3 scales the operands by a power of two (exact in range) and
// recovers NaN and infinity, neither of which can arise here: the denominators
// are O(1) to O(1e6).  This is compiler-rt's arithmetic with the contraction
// its build applies, so on Apple arm64 the result is __divdc3's to the bit
// (checked on 4e7 samples; the operand order matters, fma(d, d, c c) differing
// in a sixth of them).  libgcc's __divdc3 is Smith's method, a ULP or so away.
inline complex divide(complex w, complex z)
{
#if NICEOPS_HAS_FMA
    const double a = w.real(), b = w.imag(), c = z.real(), d = z.imag();
    const double den = std::fma(c, c, d * d);
    return complex(std::fma(a, c, b * d) / den, std::fma(b, c, -(a * d)) / den);
#else
    return w / z;
#endif
}

inline complex reciprocal(complex z)
{
    return divide(complex(1.0, 0.0), z);
}

} // namespace niceops::complex_fma
