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
// It is slow by the observer loop's standards, a second or so per frequency,
// and never runs in a footprint.

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
    double points_per_wavelength = 10.0;
    double range_points_per_wavelength = 10.0;
    // Height of the physical domain, below the absorbing layer; 0 chooses it.
    double domain_height = 0.0;
    // Thickness of the absorbing layer, in wavelengths.
    double layer_wavelengths = 100.0;
};

struct PeResult
{
    std::vector<double> range; // m
    // 20 log10 |p / p_free| per receiver height, [receiver][range]: the level
    // relative to a point source's field in free space at the straight-line
    // distance.  Rigid ground in still air gives two rays' interference.
    std::vector<std::vector<double>> relative_level;
};

PeResult solve_pe(const PeSetup& setup);

// The same averaged in energy over `samples` frequencies spread evenly in
// log frequency across the one-third-octave band centred on `centre`: what a
// band level sees, rather than one tone's interference pattern.
PeResult solve_pe_band(PeSetup setup, double centre, int samples = 5);

} // namespace niceops
