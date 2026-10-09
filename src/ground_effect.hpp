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
// This replaces nothing by default.  The broadband path used to multiply the
// database's EAA by a constant 4.0 under a "TODO Add ground effects" --
// inherited verbatim from the D original (see c74bab3) -- which was not a
// ground effect at all but a fourfold overstatement of atmospheric absorption.
// That factor is gone and golden/197421_LAE.csv was regenerated without it.
// This model is what that TODO was asking for, and it is wired into the
// spectral path (the default since the 2017 approach validation).

#include <array>
#include <complex>
#include <cstddef>
#include <string>
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

// Ground impedance models (Attenborough & Taherzadeh, "Prediction of outdoor
// ground effect", Applied Acoustics 242 (2026) 111067, sec. 1.3), the same
// formulas as panam's ground_plane.surface_admittance:
//   delany_bazley        one parameter, sigma: Z = 1 + 9.08 X^-0.75 + i 11.9 X^-0.73, X = f/sigma
//   miki                 one parameter, sigma: Z = 1 + 5.50 X^-0.632 + i 8.43 X^-0.632
//   variable_porosity    sigma_e and alpha_e (1/m): their eq. (10),
//                        Z = (1+i) sqrt(1000 sigma_e/(pi gamma rho0 f)) + i c0 alpha_e/(8 pi gamma f)
//                        (alpha_e = alpha/Omega0, the WKB form; the often-quoted 19.48 i alpha_e/f
//                        has twice the coefficient, so halve an alpha_e fitted in that form)
//   delany_bazley_layer  sigma and a depth d (m): a hard-backed layer, eq. (5),
//                        Z = Z_c coth(-i k d), Delany-Bazley's Z_c and k (eqs. 4a, 4b)
// Delany-Bazley is poor for grassland (their Table 2); on the 2017 Amedee data
// variable porosity with sigma_e ~ 200 fitted best (panam
// docs/ground_plane_corrections.md).  A thin, low-resistivity layer is not
// passive at low frequency -- ground_impedance() refuses it.
enum class ImpedanceModel
{
    delany_bazley,
    miki,
    variable_porosity,
    delany_bazley_layer,
};

// Turbulence: Chessell's factor exp(-(C/2 f sqrt(R))^2), with C the
// `turbulence` parameter, or the HARMONOISE coherence
// exp(-3/8 B k^2 gamma_p^(5/3) R gamma_T), gamma_p = h_s h_r/(h_s + h_r), B = 0.364
// (Attenborough & Taherzadeh eq. 18), with `gamma_t` ~1e-6 (moderate) to ~1e-5.
enum class TurbulenceModel
{
    chessell,
    harmonoise,
};

struct GroundModel
{
    // Specific flow resistance of the surface, kPa s/m^2 (the effective one,
    // sigma_e, for variable_porosity).
    double flow_resistance = flow_resistance::grassy_field;
    ImpedanceModel impedance_model = ImpedanceModel::delany_bazley;
    // variable_porosity: rate of change of porosity with depth, alpha_e, 1/m.
    double porosity_rate = 0.0;
    // delany_bazley_layer: depth of the hard-backed layer, m.
    double layer_depth = 0.05;
    // Receiver height above the ground, m.  Zero puts the microphone on the
    // surface, where the direct and reflected paths coincide and the model
    // tends to its +6 dB limit.
    double receiver_height = 1.2;
    // Turbulence parameter, rad s / sqrt(m).  Typically 0 to 16e-4.  It
    // decorrelates the two paths, filling in the interference nulls; zero means
    // a perfectly coherent atmosphere, which gives unrealistically deep dips.
    double turbulence = 0.0;
    TurbulenceModel turbulence_model = TurbulenceModel::chessell;
    // HARMONOISE gamma_T = C_T^2/T0^2 + 22 C_v^2/(3 c0^2).
    double gamma_t = 0.0;
    // RMS roughness height, m.  The coherent reflection falls by the Kirchhoff
    // (Ament) factor exp(-2 (k sigma cos theta)^2); the scattered part is taken
    // as lost.  It acts near normal incidence and vanishes toward grazing.
    double roughness = 0.0;
    // Include Chessell's boundary-loss correction to the plane-wave reflection
    // coefficient.  It matters at grazing incidence, which is most of a noise
    // footprint.
    bool boundary_loss_correction = true;
};

// Per-band normalized surface admittance beta = 1/Z of the ground, from its
// impedance model.  (Named "impedance" throughout for history; it is the
// admittance the reflection coefficient takes.)
//
// It depends only on the band centers and the surface, both fixed for a run,
// but it costs two pow() and a complex divide per band -- so it is computed
// once here rather than inside the observer loop, which reaches this 94 million
// times on the demo case.  `sound_speed` (m/s) enters the variable-porosity and
// layer models.  Refuses a non-passive result (Re beta <= 0).
std::vector<std::complex<double>> ground_impedance(const GroundModel& ground,
        const std::vector<double>& frequencies, double sound_speed = 343.2);

// Parse an impedance or turbulence model name (as spelled in the enums).
ImpedanceModel parse_impedance_model(const std::string& name);
TurbulenceModel parse_turbulence_model(const std::string& name);

// Excess ground attenuation into a caller-provided buffer, with the impedance
// precomputed. This is the form the observer loop uses: no allocation, no
// repeated pow(). `attenuation` must have room for frequencies.size() entries.
// With `as_power_ratio`, it receives 10^(attenuation/10) instead -- what the
// band sum multiplies by -- which saves a log per band.
void excess_ground_attenuation(const GroundModel& ground,
        const std::vector<std::complex<double>>& impedance, double source_height,
        double ground_distance, double sound_speed, const std::vector<double>& frequencies,
        double* attenuation, std::size_t first_band = 0,
        std::size_t last_band = static_cast<std::size_t>(-1), bool as_power_ratio = false);

// Per band, an upper bound on 10^(attenuation/10) -- the energy the ground can
// add -- over every geometry, for this ground and impedance.  Bands outside
// [first_band, last_band) are left at 1.  Used to extend the spectral level
// bounds (Database::spectral_direction_bounds) to runs with the ground effect:
// multiplied into the band gains, it bounds the level with the ground included.
std::vector<double> ground_gain_bounds(const GroundModel& ground,
        const std::vector<std::complex<double>>& impedance, std::size_t first_band = 0,
        std::size_t last_band = static_cast<std::size_t>(-1));

// Upper bound on |Q|, the spherical-wave reflection coefficient's size, per
// entry of `impedance`, over every geometry (the sweep ground_gain_bounds()
// describes, with its 1% margin).  Entries outside [first_band, last_band) are
// left at 1.  The ground-plane receiver (plate_receiver.hpp) bounds its level
// with it.
std::vector<double> reflection_coefficient_bounds(const GroundModel& ground,
        const std::vector<std::complex<double>>& impedance, std::size_t first_band = 0,
        std::size_t last_band = static_cast<std::size_t>(-1));

// The spherical-wave reflection coefficient Q = R_p + F(w)(1 - R_p) the model
// uses, as (real, imaginary), at grazing cosine `cos_grazing` (from the surface
// normal), kappa = pi f R_image / c, and normalized surface admittance
// `admittance`.  No roughness or turbulence: those enter the model outside Q.
std::array<double, 2> spherical_reflection_coefficient(double cos_grazing, double kappa,
        std::complex<double> admittance, bool boundary_loss = true);

// Per band, what `bounces` earlier reflections at grazing angle `grazing`
// (radians) take off a ray that arrives after them, dB: -10 lg |R_p|^2 each,
// the plane-wave coefficient.  Added into `loss` for bands in
// [first_band, last_band).  A passive ground has |R_p| <= 1, so it is never a
// gain.  For --refraction stratified's bounced arrivals (refraction.hpp,
// Arrival); the last bounce, beside the receiver, is the ground model's.
void add_bounce_loss(double grazing, int bounces,
        const std::vector<std::complex<double>>& admittance, std::size_t first_band,
        std::size_t last_band, double* loss);

// The point-receiver ground effect tabulated over the geometry, for the
// observer loop.  The model depends on the geometry only through the source
// height h_s and horizontal distance d (the receiver height is fixed for a
// run), so it is computed exactly at the nodes of a grid in
// s = sqrt(h_s / r) and ln r, r = hypot(h_s, d) -- s resolves grazing
// incidence, where the model changes fastest -- and interpolated bilinearly,
// as a power ratio, in between.  The interpolant is a convex combination of
// exact values, so ground_gain_bounds() bounds it as it bounds the model.
// Ranges outside [r_min, r_max] are clamped to the table's edge.
//
// Why: the model costs ~20 ns per band (a complex reflection coefficient and
// a Faddeeva evaluation) and the observer loop reaches it for every observer
// at every state -- 94 million times on the demo case, where it was 5.5 of
// the 7.1 s.  A lookup is four reads per band.  With the defaults (641
// angles, 40 per decade), level histories agree with the exact model to
// 0.001 dB; on the demo grid 99.97% of SELs agree to 0.01 dB, the rest being
// points where the 10 dB-down window's first-to-last rule flips on a
// secondary hump at the threshold (a perturbation of any size can do that).
// --exact_ground bypasses the table.
class GroundTable
{
public:
    GroundTable() = default;
    GroundTable(const GroundModel& ground, const std::vector<std::complex<double>>& impedance,
            double sound_speed, const std::vector<double>& frequencies, std::size_t first_band,
            std::size_t last_band, double r_min, double r_max, int angle_nodes = 641,
            int nodes_per_decade = 40);

    // 10^(attenuation/10) per band in [first_band, last_band) of the table.
    void power_ratio(double source_height, double ground_distance, double* ratio) const;

    bool empty() const
    {
        return values_.empty();
    }
    std::size_t first_band() const
    {
        return first_band_;
    }
    std::size_t last_band() const
    {
        return last_band_;
    }

private:
    std::size_t first_band_ = 0;
    std::size_t last_band_ = 0;
    std::size_t bands_ = 0;
    int n_s_ = 0;
    int n_r_ = 0;
    double log_r_min_ = 0.0;
    double log_r_step_ = 1.0;
    // [r][s][band], float: 1e-7 relative is 4e-7 dB, and half the memory
    // keeps more of a large table in cache.
    std::vector<float> values_;
};

// Excess ground attenuation in dB for each band center, for a source at
// `source_height` meters and a receiver `ground_distance` meters away
// horizontally.  `sound_speed` is m/s.
//
// Broadband form (Chessell's equations 20 and 21): the third-octave average
// rather than the pure-tone interference pattern, which is what a band level
// wants.  Pure-tone mode is deliberately not ported -- the database carries
// bands, and the pure-tone expression at a band center would put a full
// interference null where the band average has a shallow dip.
std::vector<double> excess_ground_attenuation(const GroundModel& ground, double source_height,
        double ground_distance, double sound_speed, const std::vector<double>& frequencies);

} // namespace niceops
