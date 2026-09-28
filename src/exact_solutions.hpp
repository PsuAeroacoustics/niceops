#pragma once

// Exact solutions for a point source over an impedance plane: the references
// the parabolic equation (pe.hpp) is validated against, and those of the
// Attenborough et al. benchmark (JASA 97, 173 (1995); tools/benchmark.cpp).
//
//   homogeneous_pressure  still, uniform air: the direct ray and the image ray
//                         with the spherical-wave reflection coefficient
//                         (their eqs. 12-16 with B = 1; ground_effect.hpp's Q).
//   upward_pressure       k^2 = k0^2 (1 + 2 z / R): upward refraction, the
//                         residue series of creeping waves (their eq. 11).
//   downward_pressure     k^2 = k0^2 (1 - 2 z / R): downward refraction, the
//                         normal-mode series (their eq. 10).
//
// The two refracting profiles are c0 / sqrt(1 +- 2 z / R), which is
// c0 (1 -+ z / R) near the ground, R = c0 / |dc/dz|: the forms whose modes are
// Airy functions.  The benchmark's linear profiles depart from them only far
// above the ground.
//
// Everything is exp(-i omega t), and a source whose free field is
// exp(i k R) / R, so 20 log10 |p| is the level re 1 m from the source and
// 20 log10 (|p| R) the level re free field.  Absorption enters as the
// imaginary part of the wavenumber, the admittance as in ground_effect.hpp.

#include <complex>
#include <vector>

namespace niceops {

struct ExactProblem
{
    double frequency = 0.0; // Hz
    double sound_speed = 343.0; // c0, at the ground, m/s
    double absorption = 0.0; // Np/m, of the pressure amplitude
    std::complex<double> admittance {0.0, 0.0}; // beta = 1/Z; 0 is rigid
    double source_height = 0.0; // m
    double receiver_height = 0.0; // m
};

std::complex<double> homogeneous_pressure(const ExactProblem& problem, double range);

// The mode sums at each of `ranges` (m).  `radius` is R, m.  A range where the
// series cannot be trusted to about 1e-6 is NaN: the residue series near the
// source, whose terms grow before they decay and cancel, or any range the
// mode limit leaves unconverged.
std::vector<std::complex<double>> upward_pressure(const ExactProblem& problem, double radius,
        const std::vector<double>& ranges);
std::vector<std::complex<double>> downward_pressure(const ExactProblem& problem, double radius,
        const std::vector<double>& ranges, double max_elevation = 90.0);

// On `max_elevation` (degrees): the downward profile's k^2 reaches zero at
// z = R/2, so every ray turns below it -- the near-vertical ones as well, which
// fall back to the ground within R of the source like mortar rounds, landing
// at R sin(2 elevation).  The full series carries them; the parabolic equation,
// good to 30-45 degrees and absorbing at its top, cannot, and nor does a
// profile that stays linear in c.  Modes leaving the ground steeper than
// max_elevation are dropped (tapered over its last 5 degrees), which at a
// range r removes those returns and nothing else so long as the cut lies
// between the steepest path that reaches r in one hop, asin(r / R) / 2, and
// its complement.  90 keeps everything, as their eq. (10) does.

// The zeros of Ai'(t) + q Ai(t), the mode equation of both: the first `count`
// in the order of the zeros of Ai' they continue from, then any that no zero
// of Ai' continues to (the surface wave of a soft enough ground).  Checked
// against the argument principle, and throws if that finds one missing.
// Built for |q| <= 6 and refuses more: q = i k0 l beta (downward), and the
// benchmark's grounds reach 4.9.
std::vector<std::complex<double>> impedance_airy_zeros(std::complex<double> q, int count);

// The one zero of Ai'(t) + q Ai(t) that continues the first zero of Ai' from
// q = 0, followed out along q = s * direction, for each of the increasing,
// non-negative `s`.  impedance_airy_zeros stops at |q| = 6 because it accounts
// for every zero, and the high ones crowd; the first stays well apart from its
// neighbours, so it can be followed as far as asked, to the first zero of Ai,
// a_1 - 1/q, as |q| grows.  It is the least attenuated creeping wave of an
// upward-refracting atmosphere over ground of any admittance (ShadowAttenuation).
std::vector<std::complex<double>> first_impedance_airy_zero(std::complex<double> direction,
        const std::vector<double>& s);

// H0^(1)(z), for Im z >= 0.
std::complex<double> hankel1_0(std::complex<double> z);

} // namespace niceops
