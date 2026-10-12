#pragma once

// A street canyon's response by the wave equation: the low band of the urban
// tables (urban.hpp), where specular rays stop describing a street a few
// wavelengths wide.
//
// The canyon is two-dimensional: a slot of width W and depth H in a rigid
// plane, the roofs flush at z = H on both sides and the floor at z = 0, rigid,
// and the walls rigid or of a real impedance.
// The run is reciprocal.  A line source at the microphone radiates a pulse,
// solved by a staggered pressure-velocity FDTD with a split-field PML; the
// vertical velocity across the opening is recorded; and the far field in every
// direction of the upper half-plane follows from the opening alone, because the
// roofs are rigid.  In the region z > H the Neumann Green's function of the
// half-space represents any outgoing field by the normal derivative on the
// plane z = H, which vanishes on the roofs:
//
//   p(r) = (i/2) integral over the opening of H0(k|r - r'|) (-dp/dz') dy',
//
// so the far field relative to the line source's own free field (H0/4) is
//
//   G(f, phi) = |F / F_free|^2 = 4 | h sum_y' (V/Q)(y') exp(-i k y' cos phi) |^2,
//
// V the opening's vertical velocity and Q the source's volume velocity, both
// transformed in time.  Open ground gives G_open = 4 cos^2(k z_m sin phi) for a
// microphone z_m up, and the table's entry is the band average of G over the
// band average of G_open.  phi is the direction toward the source in the
// cross-section, from the +y axis (wall A toward wall B).
//
// Oblique incidence is the same problem at a lower frequency.  A plane wave
// with a component along the street varies as exp(i k u_a x) there, so in the
// cross-section it obeys the two-dimensional Helmholtz equation at
// k |u_perp|, |u_perp| = sqrt(u_c^2 + u_z^2), arriving from
// phi = atan2(u_z, u_c); rigid walls do not mix the two.  One run therefore
// covers every direction of the three-dimensional table (canyon_wave_ratio).

#include <complex>
#include <limits>
#include <cstddef>
#include <string>
#include <vector>

namespace niceops {

struct CanyonWaveSetup
{
    double width = 20.0;              // wall to wall, m
    double depth = 20.0;              // floor to the roofs, m
    double offset = 10.0;             // microphone from wall A, m
    double microphone_height = 1.5;   // above the floor, m
    // The walls' absorption at normal incidence, as a real impedance
    // zeta = (1 + sqrt(1 - alpha)) / (1 - sqrt(1 - alpha)) times rho c, locally
    // reacting, the same in every band.  The roofs stay rigid (the far field
    // needs them so) and so does the floor, against which the open-ground
    // reference is taken.  Zero is rigid; a rigid canyon rings for seconds in
    // its cross modes (n c / 2W), which real facades damp.
    double wall_absorption = 0.05;
    // Wall B's, when it differs from wall A's (wall_absorption); NaN takes A's.
    double wall_absorption_b = std::numeric_limits<double>::quiet_NaN();
    double max_frequency = 500.0;     // the highest frequency the spectrum must resolve, Hz
    double points_per_wavelength = 15.0;
    double roof_margin = 8.0;         // roof inside the grid each side before the PML, m
    double sky_margin = 8.0;          // air above the roofs before the PML, m
    int pml_cells = 40;
    double sound_speed = 343.2;       // m/s
    double density = 1.2;             // kg/m^3
    double decay_db = 30.0;           // run until the field has fallen this far
    double max_duration = 4.0;        // s, whatever the decay
    double frequency_step = 0.25;     // Hz, at most, after zero padding
    unsigned threads = 1;
};

// What a run leaves: the opening's velocity per unit source volume velocity,
// at each frequency bin and recording point, and what the far field needs.
struct ApertureSpectrum
{
    double cell = 0.0;             // grid step, m
    double opening_height = 0.0;   // z of the recording line (H as gridded), m
    double source_height = 0.0;    // the microphone's z as gridded, m
    std::vector<double> y;         // recording points, from wall A (or the source's column), m
    std::vector<double> weight;    // each point's quadrature weight (a taper for the open test), m
    double sound_speed = 343.2;
    double frequency_step = 0.0;   // Hz
    std::size_t bins = 0;          // 0 .. bins-1, at bin * frequency_step
    std::vector<std::complex<double>> spectrum; // [bin * y.size() + point], 1/m
    double duration = 0.0;         // s the run took in simulated time
    double residual_db = 0.0;      // energy left at the end, dB re its peak
    std::size_t cells = 0;         // grid cells updated per step
    std::size_t steps = 0;

    // G at bin `bin` toward `phi` (radians, from +y, in (0, pi)), relative to
    // the line source's free field.
    double gain(std::size_t bin, double phi) const;
    // G at any frequency up to the last bin, linear in frequency between bins.
    double gain_at(double frequency, double phi) const;
    // Open ground's G for the same microphone: 4 cos^2(k z_m sin phi).
    double open_gain(double frequency, double phi) const;
    // The band average of G over the band average of G_open, frequencies in
    // [f_lo, f_hi]; the bin nearest the band's center when the band holds none.
    double band_ratio(double f_lo, double f_hi, double phi) const;
};

ApertureSpectrum solve_canyon(const CanyonWaveSetup& setup);

// The same solver over flat rigid ground, recording across a horizontal line
// `recording_height` up and `line_length` long, centered on the source, with
// a Tukey taper over `taper` of it at each end.  For checking the solver and
// the far-field transform against the image solution; `width` and `offset` of
// the setup are not used.
ApertureSpectrum solve_open_ground(const CanyonWaveSetup& setup, double recording_height, double line_length,
        double taper = 0.2);

// The table entry for a three-dimensional arrival direction whose components
// are u_c across the canyon (toward wall B), u_z upward, in the band
// [f_lo, f_hi]: the two-dimensional band ratio at the reduced frequencies
// f |u_perp|, toward phi = atan2(u_z, u_c).  Below `floor_frequency` (the band's
// reduced upper edge) the canyon is taken as transparent, 1.
double canyon_wave_ratio(const ApertureSpectrum& spectrum, double u_c, double u_z, double f_lo, double f_hi,
        double floor_frequency);

// The canyon's table gains for the receiver of `spectrum`, a canyon whose axis
// runs along compass bearing `axis_bearing`: [azimuth][elevation][band] in dB,
// each the band ratio over the third-octave band around bands[b] (edges at
// 2^(+-1/6)) at the reduced frequencies of that direction.  Only the first
// `band_count` bands are filled; the rest are left at 0 dB for the caller.
std::vector<double> canyon_wave_gains(const ApertureSpectrum& spectrum, double axis_bearing,
        const std::vector<double>& azimuths, const std::vector<double>& elevations, const std::vector<double>& bands,
        std::size_t band_count, double floor_frequency);

} // namespace niceops
