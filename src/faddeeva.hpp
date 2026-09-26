#pragma once

// The Faddeeva function w(z) = exp(-z^2) erfc(-i z), for complex z.
//
// This exists for one caller: the boundary-loss term of the excess ground
// attenuation model in ground_effect.cpp.  panam writes that term as
//
//     exp(-w) * (1 - erf(-i sqrt(w)))
//
// using SciPy's complex erf.  With z = sqrt(w) that is exactly w(z), which is
// the form with a stable algorithm behind it -- evaluating it as written
// overflows, since exp(-w) and erf of a large imaginary argument blow up in
// opposite directions.
//
// The reference implementation is Weideman's rational approximation (J.A.C. Weideman,
// "Computation of the Complex Error Function", SIAM J. Numer. Anal. 31 (1994)
// 1497-1518), which is a polynomial in a Moebius transform of z.  32 terms give
// a maximum relative error of 3e-13 over the domain the ground model actually
// reaches -- Re(z) in [0.09, 22.3], Im(z) in [-1.2, 15.6] -- measured against
// SciPy's wofz over 60,000 samples.
//
// Weideman's series converges in the upper half plane.  The lower half is
// reached by w(z) = 2 exp(-z^2) - w(-z), which is safe here only because the
// model never sends Im(z) far below zero: the exp(-z^2) in that identity grows
// like exp(Im(z)^2 - Re(z)^2), so a large negative imaginary part with a small
// real part would overflow.  The function does not guard against that, because
// no caller can produce it.

#include <complex>

namespace niceops {

std::complex<double> faddeeva(std::complex<double> z);

// Weideman's series everywhere: faddeeva() uses it outside its fast paths
// (an asymptotic series for |z| >= 6, a Taylor table below; see
// faddeeva.cpp).  Exposed so the tests can hold the fast paths to it.
std::complex<double> faddeeva_reference(std::complex<double> z);

} // namespace niceops
