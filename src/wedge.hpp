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

// The part of D that depends on the angles alone -- each of the four terms'
// cotangent and the cosine whose square is its a+- / 2 -- for a ray whose
// coefficient is wanted at many wavenumbers: three per band, every band.
struct WedgeAngles
{
    double n = 1.5;
    double cot[4] = {0.0, 0.0, 0.0, 0.0};
    double c[4] = {0.0, 0.0, 0.0, 0.0};
};

WedgeAngles wedge_angles(double n, double theta, double theta0);

// D from its angles' part: hard_wedge_coefficient(n, theta, theta0, k, r) to
// the bit, without its twelve trigonometric calls.
std::complex<double> hard_wedge_coefficient(const WedgeAngles& angles, double k, double r);

// The angles' part for one wedge and direction theta and many theta0, as
// along a Kirchhoff aperture (double_wedge.hpp): each term's angle is theta's
// part, worked out once, plus or minus theta0's, so by the angle-addition
// formulas each theta0 costs two sincos calls rather than twelve calls.  The
// same as wedge_angles to round-off (1e-8 relative in D); a term within 1e-8
// rad of its boundary, where the sum of products cancels, takes
// wedge_angles' own path.
class WedgeSweep
{
public:
    WedgeSweep(double n, double theta);
    WedgeAngles at(double theta0) const;

private:
    double n_ = 1.5;
    double theta_ = 0.0;
    // (pi + s theta) / 2n for s = +1, -1; theta / 2.
    double sin_x_[2] = {0.0, 0.0}, cos_x_[2] = {0.0, 0.0};
    double sin_half_ = 0.0, cos_half_ = 0.0;
    // n pi N for N from -3 to 4.
    double sin_npn_[8] = {}, cos_npn_[8] = {};
};

// The transition function F(X) for X >= 0.
std::complex<double> utd_transition(double x);

// The total field, geometrical optics and diffraction, of a unit plane wave
// at (r, theta), the edge at the origin; k r is the only scale.
std::complex<double> hard_wedge_field(double n, double kr, double theta, double theta0);

} // namespace niceops
