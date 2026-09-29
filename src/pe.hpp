#pragma once

// A wide-angle parabolic equation: the full-wave reference the ray models are
// checked against (roadmap step 7; tools/pe_compare.cpp).
//
// Sound from a point source over a locally reacting ground, in one vertical
// plane through a horizontally stratified atmosphere, taken in the
// effective-sound-speed approximation like the ray models.  The Pade(1,1)
// wide-angle equation marched in range by Crank-Nicolson (Salomons,
// Computational Atmospheric Acoustics, 2001, appendix F):
//   [1 + (1/4 - i k_a dr/4) Q] psi(r + dr) = [1 + (1/4 + i k_a dr/4) Q] psi(r),
//   Q = (d^2/dz^2 + k^2(z) - k_a^2) / k_a^2,
// with the ground's impedance condition dpsi/dz = -i k beta psi at z = 0 written
// into the grid, so the ground wave comes out of it with no special handling.
// That is why this rather than the Green's-function PE the roadmap named:
// the GFPE's surface-wave term is where a sign convention goes wrong unseen,
// and speed does not matter for a reference.  The time convention is
// exp(-i omega t), as in ground_effect.cpp, whose admittance it takes.
//
// It is valid for elevation angles up to about 30 degrees from the source.
// That covers the long, grazing paths where refraction matters -- the shadow
// zones -- and not the steep ones near the source, where refraction does not
// matter and the straight ray is exact anyway.
//
// It is slow by the observer loop's standards and never runs in a footprint:
// 2.7 ns per grid point per range step, so a kilometre at 1 kHz from a 100 m
// source is about a second.  A band's frequencies, and several bands at once
// (solve_pe_bands), run in parallel.
//
// The grid is 20 points per wavelength in height and 10 per wavelength in
// range.  Against a 30 x 30 reference, 10 in height was within about 0.1 dB rms
// for a microphone near the ground but 0.47 dB rms (2.2 dB at worst) for one
// 1.2 m up at 1 kHz; 20 brings that to 0.09, for twice the cost.  Rerunning
// the stored calibration cases at 20 moved the PE by 0.16-0.39 dB rms (1.1 at
// worst) and left the shadow-edge calibration against it unchanged.  Five per
// wavelength in range would halve the cost and moves most cases by only
// 0.02-0.03 dB rms, but it fails the Attenborough benchmark by 0.2-0.4 dB at
// long range, where Crank-Nicolson's phase error builds up over hundreds of
// wavelengths, so it is not the default.

#include <complex>
#include <vector>

namespace niceops {

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
    // Grid points per wavelength, vertically and in range.
    double points_per_wavelength = 20.0;
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
};

PeResult solve_pe(const PeSetup& setup);

// The same averaged in energy over `samples` frequencies spread evenly in
// log frequency across the one-third-octave band centred on `centre`: what a
// band level sees, rather than one tone's interference pattern.
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
