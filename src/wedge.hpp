#pragma once

// Diffraction by a rigid wedge, by the uniform theory of diffraction
// (Kouyoumjian and Pathak 1974): the urban tables' edge diffraction, at the
// roof edges of a street canyon (urban.hpp).
//
// The wedge's faces lie at theta = 0 and theta = n pi, the air between them
// (n = 3/2 for a building's corner, 2 for a thin screen, 1 for a flat plane).
// A plane wave of unit amplitude arrives from direction theta0, and the field
// at distance r from the edge toward theta is the geometrical-optics field --
// the incident wave where the wedge does not shadow it, and the reflection off
// either face where that face sends one -- plus the diffracted wave
//
//   u_d = D(theta, theta0) exp(i k r) / sqrt(r),
//
// for the time convention exp(-i w t) the ground model uses.  The coefficient
// for a hard (Neumann) wedge, conjugated from Kouyoumjian and Pathak's
// exp(j w t) form, is
//
//   D = -exp(i pi/4) / (2 n sqrt(2 pi k)) sum of four terms
//       cot((pi -+ beta) / 2n) F(k L a+-(beta)),  beta = theta -+ theta0,
//
// with L = r for a plane wave, a+-(beta) = 2 cos^2((2 n pi N+- - beta) / 2) and
// N+- the integers nearest (beta +- pi) / (2 n pi).  The transition function,
//
//   F(X) = -2i sqrt(X) exp(-iX) integral from sqrt(X) to infinity of exp(i t^2) dt
//        = -i sqrt(pi X) exp(i pi/4) w(exp(i pi/4) sqrt(X)),
//
// through the Faddeeva function w, tends to 1 away from the shadow and
// reflection boundaries, where the geometrical-optics field jumps, and keeps
// the total field continuous across them.  For a half-plane under plane-wave
// incidence this is Sommerfeld's exact solution; for a right-angle wedge it is
// within 0.016 dB of the exact eigenfunction series at k r = 10 and 0.005 dB at
// 40 (tests/test_wedge.cpp).

#include <complex>

namespace niceops {

// The diffraction coefficient D above, in m^(1/2).
std::complex<double> hard_wedge_coefficient(double n, double theta, double theta0, double k, double r);

// The transition function F(X) for X >= 0.
std::complex<double> utd_transition(double x);

// The total field, geometrical optics and diffraction, of a unit plane wave
// at (r, theta), the edge at the origin; k r is the only scale.
std::complex<double> hard_wedge_field(double n, double kr, double theta, double theta0);

} // namespace niceops
