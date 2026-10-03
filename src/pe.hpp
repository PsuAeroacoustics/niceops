#pragma once

// A wide-angle parabolic equation: the full-wave reference the ray models are
// checked against (roadmap step 7; tools/pe_compare.cpp).
//
// Sound from a point source over a locally reacting ground, in one vertical
// plane through a horizontally stratified atmosphere, taken in the
// effective-sound-speed approximation like the ray models.  The one-way wave
// equation for psi = p sqrt(r) exp(-i k_a r),
//   d psi / dr = i k_a (sqrt(1 + Q) - 1) psi,
//   Q = (d^2/dz^2 + k^2(z) - k_a^2) / k_a^2,
// on a grid in height, with the ground's impedance condition
// dpsi/dz = -i k beta psi at z = 0 written into the grid, so the ground wave
// comes out of it with no special handling.  That is why this rather than the
// Green's-function PE the roadmap named: the GFPE's surface-wave term is where
// a sign convention goes wrong unseen, and it splits refraction from
// diffraction, which limits its step where the gradient is strong -- near the
// ground, where the shadows are.  The time convention is exp(-i omega t), as in
// ground_effect.cpp, whose admittance it takes.
//
// Two ways to march it in range (PeSetup::scheme):
//
//   split_step_pade (the default).  psi(r + dr) = exp(i k_a dr (sqrt(1 + Q) -
//   1)) psi(r) exactly, the exponential replaced by a sum of resolvents of Q
//   (Collins 1993; pade_propagator.hpp), so a step is a handful of tridiagonal
//   solves and is eight wavelengths long: about ten solves where
//   Crank-Nicolson takes eighty steps.  It is accurate to 1e-7 a step for
//   every mode within max_angle of horizontal, 50 degrees, whatever the step:
//   from 2 to 32 wavelengths a step the result does not move by 0.001 dB,
//   60 dB into a shadow.  Steeper modes it absorbs, where Crank-Nicolson
//   carries them at the wrong speed; a level on a path steeper than max_angle
//   from the source is reported as NaN.
//
//   A propagator good to 50 degrees is not by itself a field good to 50
//   degrees.  Two other things had to change, and with max_angle alone raised
//   the error at 30-50 degrees stayed where it was, 0.6-1.7 dB rms against the
//   exact two-ray field at 1 kHz:
//     - the starter.  The Gaussian is a far-field, narrow-angle source, and it
//       is normalised with the ground's wavenumber.  The self-starter
//       (PeStarter::self) is the first range step itself with the point
//       source's modal amplitudes in it, applied to a delta function on the
//       grid: exact to max_angle like every other step.  It also gives a
//       source aloft its own strength.  Where k at the source is 3% from k at
//       the ground the Gaussian is 0.11-0.15 dB off at every range, against
//       the residue series; the self-starter is within 0.04.
//     - the grid.  Second-order differences misplace a wave at elevation
//       theta by (k dz sin theta)^2 / 12 of its vertical wavenumber squared,
//       a radian of phase by 35 degrees at 1 kHz from a 100 m source.
//       vertical_order 4 takes d^2/dz^2 as Numerov does, which costs a weight
//       per term per row and a third more time.
//     - the grid's edges (boundary_order 4).  With the operator alone at
//       fourth order the grid still converged at second, because the ground
//       condition was a three-point one-sided derivative and the source and
//       receivers met the grid by linear interpolation.  The ground row is
//       now fourth order too (pe.cpp, pade_march), with psi_0 carried rather
//       than derived, and both interpolations are cubic.
//   Together, at ten points a wavelength: within 0.07 dB rms of the exact
//   field at every elevation to 60 degrees, and within 0.01 to 50 at twenty.
//
//   crank_nicolson.  The Pade(1,1) wide-angle equation marched by
//   Crank-Nicolson (Salomons, Computational Atmospheric Acoustics, 2001,
//   appendix F), ten steps a wavelength:
//     [1 + (1/4 - i k_a dr/4) Q] psi(r + dr) = [1 + (1/4 + i k_a dr/4) Q] psi(r).
//   What this was before the split-step, kept as an independent check on it.
//   Good to about 30 degrees, and its phase error builds with range: five
//   steps a wavelength fails the Attenborough benchmark by 0.2-0.4 dB at long
//   range.  The two agree to 0.03-0.12 dB rms on a windy profile out to
//   1.5 km, all of it Crank-Nicolson's.
//
// Either covers the long, grazing paths where refraction matters -- the shadow
// zones -- and not the steep ones near the source, where refraction does not
// matter and the straight ray is exact anyway.
//
// It never runs in a footprint.  One march is 39 times faster than
// Crank-Nicolson at the defaults.  A band's frequencies, and several bands at
// once (solve_pe_bands), run in parallel.
//
// The grid is 10 points per wavelength in height where it is fourth order
// throughout, which is the default, and 20 otherwise.  Fourth order is what
// allows it: on a windy profile out to 1.5 km, against 80 points,
//                                          20 points      10 points
//   second order throughout                0.02-0.04 dB   0.06-0.19 dB rms
//   operator fourth, edges second          0.01-0.02      0.05-0.09
//   fourth order throughout                < 0.0005       0.001-0.002
// and 6 points are still within 0.01.  (Crank-Nicolson at 10, for a
// microphone 1.2 m up at 1 kHz, was 0.47 dB rms and 2.2 at worst from a
// 30 x 30 reference; rerunning the stored calibration cases at 20 moved the PE
// by 0.16-0.39 dB rms and left the shadow-edge calibration against it
// unchanged.)

#include <complex>
#include <vector>

namespace niceops {

// How the equation is marched in range (the header comment has both).
enum class PeScheme
{
    split_step_pade,
    crank_nicolson,
};

// The field the march starts from.
enum class PeStarter
{
    // Collins' self-starter: the first step is the propagator with the point
    // source's modal amplitudes, applied to a delta function on the grid.  It
    // is exact to max_angle like every later step, and the ground enters
    // through the grid's own condition.
    self,
    // Salomons' second-order Gaussian at r = 0, with its image in the ground
    // at the normal-incidence reflection coefficient: a far-field, narrow-
    // angle source.
    gaussian,
};

struct PeSetup
{
    double frequency = 0.0; // Hz
    double source_height = 0.0; // m
    std::vector<double> receiver_heights; // m
    double max_range = 0.0; // m
    double output_step = 10.0; // m, range spacing of the results
    // Normalised ground admittance 1/Z, exp(-i omega t); 0 is rigid.
    std::complex<double> admittance {0.0, 0.0};
    // The effective sound speed along the plane: linear between (z, c)
    // nodes, constant outside them.
    std::vector<double> z;
    std::vector<double> c;
    // Atmospheric absorption, Np/m of the pressure amplitude (dB/m / 8.686):
    // the imaginary part of the wavenumber, so it acts along the path the
    // sound actually takes.  0, as in the comparisons with the ray models,
    // leaves it out.
    double absorption = 0.0;
    PeScheme scheme = PeScheme::split_step_pade;
    // Grid points per wavelength in height.  0 chooses: 10 for the split-step
    // march with vertical_order and boundary_order both 4, 20 otherwise.
    double points_per_wavelength = 0.0;
    // Split-step Pade: the longest range step, in wavelengths -- the step is
    // the largest that divides output_step and does not exceed it -- the
    // elevation angle the propagator is accurate to, in degrees, and its
    // tolerance per step, as |r - exact| for every mode within that angle.
    double range_step_wavelengths = 8.0;
    double max_angle = 50.0;
    double propagator_tolerance = 1e-7;
    // Split-step Pade: how the field is started.  Crank-Nicolson always takes
    // the Gaussian.
    PeStarter starter = PeStarter::self;
    // Split-step Pade: the order of d^2/dz^2 on the grid, 2 or 4.
    int vertical_order = 4;
    // Split-step Pade: the order of what meets the grid at its edges, 2 or 4
    // -- the ground condition (4 only with vertical_order 4), and the
    // interpolation that puts the source on the grid and takes the receivers
    // off it.
    int boundary_order = 4;
    // Crank-Nicolson: range steps per wavelength.
    double range_points_per_wavelength = 10.0;
    // Height of the physical domain, below the absorbing layer; 0 chooses it.
    double domain_height = 0.0;
    // Thickness of the absorbing layer, in wavelengths.
    double layer_wavelengths = 100.0;
    // Turbulence, by random phase screens (Gilbert, Raspet and Di 1990;
    // Salomons 2001, appendix I): off at 0.  gamma_t is the HARMONOISE
    // strength, C_T^2/T0^2 + 22 C_v^2 / (3 c0^2) (ground_effect.hpp), so the
    // effective index's structure parameter is gamma_t / 4 -- the value whose
    // spherical-wave coherence Ostashev's form (Attenborough & Taherzadeh 2026,
    // eq. 20) and HARMONOISE's (eq. 18) both reduce to.  The spectrum is von
    // Karman's with outer scale `outer_scale`, 0.033 C^2 (kappa^2 + 1/L0^2)^(-11/6),
    // uniform in height.  One screen every `screen_spacing` metres of range (0:
    // the outer scale) carries the phase a slab that thick accumulates, which
    // is the Markov approximation; its vertical spectrum is the 3-D one
    // integrated across the plane, so the vertical coherence is the 3-D
    // medium's.  Each realisation draws its screens from `seed`.
    double gamma_t = 0.0;
    double outer_scale = 10.0; // m
    double screen_spacing = 0.0; // m
    unsigned seed = 1;
    // Threads solve_pe_band spreads its frequencies over (0: one per hardware
    // thread).  The result does not depend on it.
    unsigned threads = 0;
};

struct PeResult
{
    std::vector<double> range; // m
    // 20 log10 |p / p_free| per receiver height, [receiver][range]: the level
    // relative to a point source's field in free space at the straight-line
    // distance.  Rigid ground in still air gives two rays' interference.
    std::vector<std::vector<double>> relative_level;
    // psi at each receiver and output range, [receiver][range]: for ensemble
    // means of the complex field (tests).  solve_pe_band leaves it empty.
    std::vector<std::vector<std::complex<double>>> psi;
    // What the march took: its range step, m, and the split-step propagator's
    // terms (0 for Crank-Nicolson).  solve_pe_band leaves them unset.
    double range_step = 0.0;
    int pade_terms = 0;
    // Whether the self-starter was used.  It is not for Crank-Nicolson, for a
    // source within a grid point of the ground, or if its fit was refused.
    bool self_started = false;
};

PeResult solve_pe(const PeSetup& setup);

// The same averaged in energy over `samples` frequencies spread evenly in
// log frequency across the one-third-octave band centred on `centre`: what a
// band level sees, rather than one tone's interference pattern.  The frequencies
// are the midpoints of `samples` equal slices of the band (one sample is the
// centre), so each stands for an equal share of it; they used to include both
// band edges at full weight, which read 0.4 to 0.6 dB off the band's mean.
// With turbulence each frequency is also its own realisation (seed + k), so the
// average is over the band and the ensemble at once.  The frequencies are
// independent marches and run in parallel (PeSetup::threads); they are summed
// in order afterwards, so the result is the same on any thread count.
PeResult solve_pe_band(PeSetup setup, double centre, int samples = 5);

// solve_pe_band for several bands at once, band b from setups[b] at centres[b]:
// every (band, frequency) march is one task in a single pool, the highest
// frequencies -- the longest marches, their cost going as f^2 -- dealt out
// first.  The same result as a band at a time, sooner.  Threads from
// setups[0].threads.
std::vector<PeResult> solve_pe_bands(const std::vector<PeSetup>& setups, const std::vector<double>& centres,
        int samples = 5);

// The PE's phase screen: phi(z_j) at z_j = j dz, j = 0 .. count - 1, for a slab
// `thickness` metres deep at wavenumber k.  Its variance is
// 0.782 k^2 thickness (gamma_t / 4) outer_scale^(5/3), less what the finite
// column leaves out (tests).
std::vector<double> phase_screen(double k, double gamma_t, double outer_scale, double thickness,
        double dz, std::size_t count, unsigned long long seed);

} // namespace niceops
