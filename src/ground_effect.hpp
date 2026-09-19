#pragma once

// Excess ground attenuation for a non-directional point source over level,
// finite-impedance ground.
//
// Ported from panam's ega() in flight_acoustics.py, which follows:
//   Chien and Soroka, "Sound Propagation Along an Impedance Plane",
//     J. Sound Vib. 43(1) 9-20, 1975 (corrected 1980)
//   Delany and Bazley, "Acoustical Properties of Fibrous Absorbent Materials",
//     Appl. Acoust. 3, 105-116, 1970
//   Chessell, "Propagation of noise along a finite impedance boundary",
//     JASA 62(4), 825-834, 1977
//   Daigle, Embleton and Piercy, JASA 66(3), 918-919, 1979
//   Stusnick, Plotkin and Sutherland, Wyle Research Report WR 85-19, 1985
//
// Assumes level terrain, no wind, and straight rays.
//
// **The result is a gain, not a loss.**  It is added to the direct-wave level:
// over hard ground it tends to +6 dB, the pressure doubling of a rigid
// reflector, and it goes sharply negative in the ground-interference dip, where
// the reflected path arrives out of phase.  Calling it an "attenuation" is the
// literature's convention, not a sign.
//
// This replaces nothing by default.  The existing broadband path multiplies the
// database's EAA by a constant 4.0 under a "TODO Add ground effects" -- inherited
// verbatim from the D original (see c74bab3) and encoded in the goldens -- which
// is not a ground effect at all but a fourfold overstatement of atmospheric
// absorption.  This model is what that TODO was asking for, and it is wired into
// the spectral path, which is opt-in.

#include <complex>
#include <cstddef>
#include <vector>

namespace niceops {

// Specific flow resistance of common surfaces, kPa s/m^2, from the AAM
// technical reference.  The parameter spans five orders of magnitude and the
// answer is sensitive to it, so name the surface rather than guessing a number.
namespace flow_resistance {
inline constexpr double snow_cover = 30.0;
inline constexpr double grassy_field = 225.0;
inline constexpr double roadside_soil = 650.0;
inline constexpr double packed_sand = 1650.0;
inline constexpr double hard_packed_dirt = 3000.0;
inline constexpr double exposed_rock = 6000.0;
inline constexpr double concrete = 10000.0;
inline constexpr double asphalt = 50000.0;
inline constexpr double water = 1.0e6;
} // namespace flow_resistance

struct GroundModel
{
    // Specific flow resistance of the surface, kPa s/m^2.
    double flow_resistance = flow_resistance::grassy_field;
    // Receiver height above the ground, m.  Zero puts the microphone on the
    // surface, where the direct and reflected paths coincide and the model
    // tends to its +6 dB limit.
    double receiver_height = 1.2;
    // Turbulence parameter, rad s / sqrt(m).  Typically 0 to 16e-4.  It
    // decorrelates the two paths, filling in the interference nulls; zero means
    // a perfectly coherent atmosphere, which gives unrealistically deep dips.
    double turbulence = 0.0;
    // Include Chessell's boundary-loss correction to the plane-wave reflection
    // coefficient.  It matters at grazing incidence, which is most of a noise
    // footprint.
    bool boundary_loss_correction = true;
};

// Per-band ground impedance ratio, Delany-Bazley.
//
// It depends only on the band centres and the surface, both fixed for a run,
// but it costs two pow() and a complex divide per band -- so it is computed
// once here rather than inside the observer loop, which reaches this 94 million
// times on the demo case.
std::vector<std::complex<double>> ground_impedance(const GroundModel& ground,
        const std::vector<double>& frequencies);

// Excess ground attenuation into a caller-provided buffer, with the impedance
// precomputed. This is the form the observer loop uses: no allocation, no
// repeated pow(). `attenuation` must have room for frequencies.size() entries.
void excess_ground_attenuation(const GroundModel& ground,
        const std::vector<std::complex<double>>& impedance, double source_height,
        double ground_distance, double sound_speed, const std::vector<double>& frequencies,
        double* attenuation);

// Excess ground attenuation in dB for each band centre, for a source at
// `source_height` metres and a receiver `ground_distance` metres away
// horizontally.  `sound_speed` is m/s.
//
// Broadband form (Chessell's equations 20 and 21): the third-octave average
// rather than the pure-tone interference pattern, which is what a band level
// wants.  Pure-tone mode is deliberately not ported -- the database carries
// bands, and the pure-tone expression at a band centre would put a full
// interference null where the band average has a shallow dip.
std::vector<double> excess_ground_attenuation(const GroundModel& ground, double source_height,
        double ground_distance, double sound_speed, const std::vector<double>& frequencies);

} // namespace niceops
