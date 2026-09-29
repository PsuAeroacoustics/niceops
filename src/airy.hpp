#pragma once

// The Airy function Ai and its derivative for complex argument.
//
// This exists for the exact solutions (exact_solutions.hpp) the parabolic
// equation is validated against: in an atmosphere whose k^2 is linear in
// height, every normal mode is an Airy function of a complex argument, and no
// standard library provides one.
//
// Two expansions meet at |z| = 5: the Maclaurin series inside, the asymptotic
// expansion outside, carried into the sector around the negative real axis
// by the connection formula Ai(z) = -w Ai(w z) - w^2 Ai(w^2 z), w = exp(2 pi i/3).
// Five is where their errors cross.  Near the positive real axis, where Ai is
// exponentially small, the series cancels exp(4/3 |z|^1.5) of its digits away;
// the asymptotic series, stopped at its smallest term, cannot do better than
// exp(-4/3 |z|^1.5).  At |z| = 5 both are ~1e-7 relative, and much better away
// from that ring.  tests/test_pe_shadow.cpp holds the function to values
// computed independently, by Bessel-integral quadrature, and the two halves to
// each other where they meet.

#include <complex>

namespace niceops {

// Ai(z) and Ai'(z) together: every caller wants both.
void airy(std::complex<double> z, std::complex<double>& ai, std::complex<double>& ai_prime);

// The two expansions on their own, for the tests.
void airy_series(std::complex<double> z, std::complex<double>& ai, std::complex<double>& ai_prime);
void airy_asymptotic(std::complex<double> z, std::complex<double>& ai,
        std::complex<double>& ai_prime); // |arg z| <= 2 pi / 3

// The n-th zero of Ai' (n = 1, 2, ...), all on the negative real axis: the
// asymptotic estimate refined by Newton.
double airy_prime_zero(int n);

} // namespace niceops
