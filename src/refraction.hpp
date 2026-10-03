#pragma once

// Refraction by a linear effective sound-speed profile: step 4 of the
// propagation roadmap (--refraction linear).
//
// Where the sound speed along the path is c(z) = c0 + a z, every ray is an arc
// of a circle centred on the line where c would vanish, z = -c0/a.  So the
// eigenrays between a source and a receiver have closed forms:
//   - the direct ray is the one circle through both points with its centre on
//     that line;
//   - the ground-reflected ray is found by a one-dimensional solve on the
//     launch angle, along which its range is monotonic;
//   - under upward refraction (a < 0), past the ray that grazes the ground
//     there is no ray at all: the receiver is in the shadow zone.
// This is how Harmonoise and Nord2000 treat refraction.  A real profile
// becomes an equivalent linear one per path, by a least-squares fit over the
// heights the path spans (EffectiveSoundSpeed::fit).  Fitted over a straight
// chord, uniform weighting in height is uniform weighting along the path.
//
// The effective sound speed is c(z) + w(z).h, with w the wind and h the
// horizontal direction of propagation.  This is the effective-sound-speed
// approximation.  It overstates the wind's effect by 1/cos(elevation), which is
// worst on the steep paths where refraction matters least.  The general
// stratified model, planned for later, uses the moving-medium invariants
// instead, and this one is the oracle it will be tested against.
//
// Everything here is metres, seconds and radians.

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <vector>

#include "atmosphere_profile.hpp"
#include "propagation.hpp"
#include "vec3.hpp"

namespace niceops {

// c(z) = sound_speed + gradient * z along one path: m/s and 1/s.
struct LinearProfile
{
    double sound_speed = 0.0;
    double gradient = 0.0;
};

// The direct eigenray from a source at height z_s to a receiver at z_r, a
// horizontal distance r > 0 away.  Past the limiting ray (a < 0 only) there is
// no direct ray; the fields then describe the limiting ray to the shadow
// boundary, continued `shadow_depth` metres along the ground, which is where
// the diffracted and scattered energy that does arrive comes from.
struct DirectEigenray
{
    bool lit = true;
    // Elevation of the ray at the source and at the receiver, radians above
    // the horizontal (negative descending).
    double launch = 0.0;
    double arrival = 0.0;
    // The distance spherical spreading is equivalent to: r / cos(launch).  In
    // a linear profile the ray tube's area works out to exactly that; it is
    // the slant range when a = 0.
    double spreading_range = 0.0;
    // Arc length, and travel time along it.
    double path_length = 0.0;
    double travel_time = 0.0;
    // Metres past the shadow boundary, 0 when lit.
    double shadow_depth = 0.0;
    // The shadow boundary itself: the horizontal range from the source at which
    // the limiting ray, grazing the ground, comes back up to the receiver's
    // height.  Lit rays carry it too, for the transition into the shadow
    // (ShadowTransition).  Infinite where the profile makes no such shadow.
    double boundary = std::numeric_limits<double>::infinity();
    // The path's linear fit's own boundary.  For --refraction linear it is
    // `boundary`; through layers it lies beyond, since a log law's steep lowest
    // metres bring the layered boundary in.  ShadowTransition puts each band's
    // effective boundary between the two.
    double linear_boundary = std::numeric_limits<double>::infinity();
    // The gradient (1/s, negative) of the linear profile whose grazing ray makes
    // `boundary`, r_b = sqrt(2 R z_s) + sqrt(2 R z_r): for --refraction linear the
    // profile's own.  0 where there is no boundary.
    double boundary_gradient = 0.0;
    // --refraction stratified: the boundary was set although air above the
    // source outruns the ground, so rays turned back from aloft (or bounced
    // along a duct) may arrive beyond it.  It stands only where none does;
    // layered_direct_ray takes it back where one does.
    bool boundary_under_aloft = false;
};

DirectEigenray direct_eigenray(const LinearProfile& profile, double z_s, double z_r, double r);

// The ground-reflected eigenray from z_s to a receiver z_r > 0 above the
// ground, a horizontal distance r away.  `exists` is false past the limiting
// ray, where nothing reflects to the receiver either.
struct ReflectedEigenray
{
    bool exists = false;
    double launch = 0.0;
    // Angle between the ray and the ground where it reflects, radians.
    double grazing = 0.0;
    // Total arc length, source to ground to receiver.
    double path_length = 0.0;
};

ReflectedEigenray reflected_eigenray(const LinearProfile& profile, double z_s, double z_r,
        double r);

// Horizontal distance from the source to the shadow boundary, for a receiver
// at z_r: infinite unless the gradient is negative.
double shadow_boundary(const LinearProfile& profile, double z_s, double z_r);

// The effective sound speed of a stratified atmosphere, tabulated so that the
// equivalent linear profile of any path costs a handful of operations: the
// running integrals of c(z), and of the two local wind components, and of z
// times each.  A least-squares line over [z1, z2] needs nothing else.
class EffectiveSoundSpeed
{
public:
    // Running integrals from the table's bottom: of c, z c, w_x, z w_x, w_y, z w_y.
    using Moments = std::array<double, 6>;

    EffectiveSoundSpeed() = default;
    // Tabulate over [z_low, z_high] metres.  The wind is turned into the local
    // frame, whose +x axis points along compass bearing `frame_bearing`, in
    // degrees.  A profile without wind is still air.  The profile must carry
    // pressure (with_hydrostatic_pressure).
    EffectiveSoundSpeed(const AtmosphereProfile& profile, double frame_bearing, double z_low,
            double z_high);

    Moments moments(double z) const;

    // The equivalent linear profile along horizontal direction (hx, hy), a
    // unit vector, for a path between heights z_a and z_b whose moments are
    // given.  Spans under min_fit_span are widened upward, so a level path fits
    // the air it actually crosses rather than a single point.
    LinearProfile fit(double z_a, const Moments& at_a, double z_b, const Moments& at_b, double hx,
            double hy) const;

    // The local-frame wind at height z, m/s.
    std::array<double, 2> wind(double z) const;

    static constexpr double min_fit_span = 10.0; // m

private:
    void locate(double z, std::size_t& node, double& fraction) const;

    double z_low_ = 0.0;
    double step_ = 1.0;
    double inv_step_ = 1.0;
    std::size_t nodes_ = 0;
    // [node][quantity]: c, w_x, w_y at the node.
    std::vector<std::array<double, 3>> values_;
    // [node]: the six running integrals up to the node.
    std::vector<Moments> moments_;
};

// Which ground receiver the level is for, and so what ground geometry a ray
// hands on: a point receiver's ground model needs the reflected ray, a flush
// plate microphone the direct ray's arrival, and no ground model neither.
enum class GroundGeometry
{
    none,
    point,
    plate,
};

// What ray theory leaves out of an arrival that --tunnelling or --skip_fill
// adds (stratified.hpp, layered_wave_arrivals), as a per-band loss on top of its
// ray tube.  Both are Airy tails, and both grow with frequency:
//   - evanescent: the ray turns at z_t above the receiver, which it reaches as
//     the tail of the field below its turning point, Ai^2(x) with
//     x = (3/2 omega I)^(2/3), I = int_{z_r}^{z_t} sqrt(p^2 - 1/c(z)^2) dz in
//     seconds;
//   - caustic: the receiver lies on the dark side of a caustic where a family's
//     range X(p) is stationary, a distance d short of (or past) it, and takes
//     Ai^2(x) with x = omega^(2/3) d (2 / X'')^(1/3) -- caustic_depth holds the
//     d (2 / X'')^(1/3), in s^(2/3).
// The tube the arrival carries is the one at wave_reference_frequency (see
// add_wave_loss), so the loss is never negative and the peak search's bound on
// the tube still bounds the band.
struct WaveCorrection
{
    double evanescent = 0.0;
    bool caustic = false;
    double caustic_depth = 0.0;

    bool active() const { return evanescent > 0.0 || caustic; }
};

inline constexpr double wave_reference_frequency = 10000.0; // Hz

// Add the correction's loss, dB, to loss[band] for bands in [first_band,
// last_band), and hold the band's gain over spherical spreading -- its tube,
// `spreading_range`, against the straight-line `chord` -- to
// max_focusing_gain_db, the caustic cap every ray tube takes.
void add_wave_loss(const WaveCorrection& wave, double chord, double spreading_range,
        const std::vector<double>& frequencies, std::size_t first_band, std::size_t last_band,
        double* loss);

// One eigenray besides a pair's primary, under --refraction stratified: in a
// duct, rays that have bounced off the ground `bounces` times before arriving
// from above, and the direct rays of families other than the primary's.  Each
// is a level of its own -- its own launch direction on the sphere, its own ray
// tube -- summed with the primary's in energy, since a band averages away
// their phases.  The ground term at the receiver is the same model as the
// primary's, through the same equivalent straight geometry; each earlier bounce
// costs the plane-wave |R|^2 at `grazing`.
struct Arrival
{
    int bounces = 0;
    double launch = 0.0;  // rad above the horizontal at the source
    double arrival = 0.0; // rad, negative descending
    // Angle to the ground where the ray meets it, rad; 0 for one that turns
    // above the ground and never reaches it, which takes no ground term.
    double grazing = 0.0;
    double spreading_range = 0.0;
    double path_length = 0.0;
    double travel_time = 0.0;
    // --tunnelling, --skip_fill: an arrival ray theory does not have, and its
    // per-band loss.  One that tunnels takes no ground term, plate or not.
    WaveCorrection wave;
    // Formed by layered_direct_ray; its ground geometry by the ground's finish.
    RayPath path;
};

// One source-observer pair under --refraction linear: the RayPath for the
// level, and what else the run needs from the ray.
struct RefractedRay
{
    RayPath path;
    double travel_time = 0.0; // s
    double shadow_depth = 0.0; // m
    LinearProfile profile; // the path's equivalent linear profile
    // What finish_ground() works from: whether the ray bends at all (a
    // straight one's ground geometry is already the real one), the source's
    // height and horizontal distance, and the direct eigenray.
    bool curved = false;
    double source_z = 0.0;
    double receiver_z = 0.0;
    double horizontal = 0.0;
    DirectEigenray direct;
    // Set by layered_direct_ray when it hands back the linear fit's ray past
    // the layered shadow boundary (see there); finished as a linear ray.
    bool from_linear_fit = false;
    // --refraction stratified in a surface duct, past the reach of every ray
    // that arrives without touching the ground: the primary is then the
    // strongest bounced arrival, having met the ground `bounces` times at
    // `ground_grazing` before its last descent.
    int bounces = 0;
    double ground_grazing = 0.0;
    // The pair's other eigenrays (Arrival).  Empty but for --refraction
    // stratified, and there only where more than one ray arrives.
    std::vector<Arrival> extra;
};

// The refracted path from `source` to the microphone at `receiver` --
// the observer's point raised by the microphone's height above the ground,
// which is where the level is evaluated whichever the path model.
// `at_source` and `at_receiver` are the table's moments at their heights.
// Paths too near vertical, or through air too close to uniform to bend them
// (|gradient| under 1e-9 /s), are straight, to the same microphone, so uniform
// air gives back straight_ray() exactly.
//
// The ground is the plane z = 0, as in the straight model.  With a point
// receiver, the reflected ray runs from the source to a microphone
// `microphone_height` above that plane, whatever the observer's own z, just as
// the ground models place it.
RefractedRay refracted_ray(const EffectiveSoundSpeed& air, const Vec3& source,
        const EffectiveSoundSpeed::Moments& at_source, const Vec3& receiver,
        const EffectiveSoundSpeed::Moments& at_receiver, GroundGeometry ground,
        double microphone_height);

// refracted_ray() in two halves.  The first gives everything the level's bound
// needs -- the launch direction and the ray tube's range, so a state can be
// ruled out by DirectionBounds exactly as a straight one is (see
// windowed_sound_exposure_level) -- plus the travel time.  Its ground geometry
// is the straight one until the second half, which solves for the reflected
// ray when a point receiver needs it, and so is left until a state survives
// the bound.
RefractedRay direct_ray(const EffectiveSoundSpeed& air, const Vec3& source,
        const EffectiveSoundSpeed::Moments& at_source, const Vec3& receiver,
        const EffectiveSoundSpeed::Moments& at_receiver);
void finish_ground(RefractedRay& ray, GroundGeometry ground, double microphone_height);

// Per band, the shadow zone's attenuation rate for a unit gradient and sound
// speed.  This is the first creeping-wave mode over a rigid ground (Pierce,
// Acoustics, sec. 9-5):
//   alpha = 8.686 tau_1 sin(pi/3) (pi f)^(1/3) |a|^(2/3) / c0   dB/m,
// with tau_1 = 1.0188, the first zero of Ai'.  This returns everything but the
// |a|^(2/3)/c0.  Rigid ground attenuates least, so this errs toward higher
// levels.  Bands outside [first_band, last_band) are 0.
std::vector<double> shadow_attenuation_rates(const std::vector<double>& frequencies,
        std::size_t first_band, std::size_t last_band);

// The shadow zone's attenuation per band: the least attenuated creeping wave of
// the path's equivalent linear profile, over rigid ground or over ground of a
// given admittance.
//
// Near the ground c0 + a z is k^2 = k0^2 (1 + 2 z / R), R = c0 / |a|, whose
// waves are Ai(b - w z / l), w = exp(2 pi i / 3), l = (R / 2 k0^2)^(1/3), with
// Ai'(b) + q Ai(b) = 0 and q = -i w^2 beta k0 l (exact_solutions.cpp derives
// it).  The first decays at
//   alpha = Im(b / w) k0 / (2 s^2)  Np/m,   s = k0 l = (pi f / |a|)^(1/3).
// Over rigid ground b = a'_1 and this is shadow_attenuation_rates().  Over any
// other, b is first_impedance_airy_zero() along q = s (-i w^2 beta): beta is
// fixed per band for a run and s carries the path, so each band gets one table
// in log s.  Im(b / w) is 0.882 over rigid ground and tends to 2.025, the first
// zero of Ai (a pressure-release ground), as |q| grows.  It is not monotonic:
// the first-order shift -q / a'_1 lowers it wherever arg beta lies in (-90,
// -30) degrees -- always for variable porosity (-45), and for Delany-Bazley and
// Miki except on soft ground at high frequency, where they reach -4 -- so at
// small |q| real ground can attenuate less than rigid, by up to 12% (Miki,
// |q| = 0.47), before it attenuates up to 2.3 times as much.  Both keep Pierce's first-order form of
// alpha, so a ground tending to rigid tends to the rigid rate.  Against the
// exact residue series' slope deep in the shadow it is within 0.3%
// (tests/test_shadow_attenuation.cpp).
//
// Left out: the surface wave, a zero near q^2 that no zero of Ai' continues
// to.  Upward it needs Re q > 0, which is arg beta < -60 degrees.  Delany-Bazley,
// Miki and variable porosity without a porosity rate never go below -57, at
// 10 Hz - 10 kHz and 10 - 10000 kPa s/m^2 (tests/test_shadow_attenuation.cpp).
// A porosity rate of tens per metre, or a thin hard-backed layer, does, in the
// low bands.  Along such a q the first zero turns into the surface wave (it
// runs out toward q^2) and cannot be followed, so those bands take the rigid
// rate, and admits_surface_wave() names them for the run's log.
class ShadowAttenuation
{
public:
    ShadowAttenuation() = default;
    // Rigid ground: exactly shadow_attenuation_rates().
    ShadowAttenuation(const std::vector<double>& frequencies, std::size_t first_band,
            std::size_t last_band);
    // Ground whose normalised admittance in each band is admittance[band]
    // (ground_impedance()); 0 is rigid.
    ShadowAttenuation(const std::vector<double>& frequencies,
            const std::vector<std::complex<double>>& admittance, std::size_t first_band,
            std::size_t last_band);

    bool rigid() const { return table_.empty(); }
    // dB/m in `band` along a path whose equivalent linear profile is `profile`.
    double rate(std::size_t band, const LinearProfile& profile) const;
    // Im(b / w) in `band` at s: what the table holds (tests).
    double decay_factor(std::size_t band, double s) const;
    const std::vector<double>& rigid_rates() const { return rigid_rates_; }
    // Whether the ground admits the surface wave left out above, in `band`,
    // which then takes the rigid rate.
    bool admits_surface_wave(std::size_t band) const
    {
        return !surface_wave_.empty() && surface_wave_[band];
    }

    // s's table: log s from log(s_min) in steps of ln(10) / per_decade.
    static constexpr double s_min = 0.5;
    static constexpr double s_max = 1e5;
    static constexpr int per_decade = 32;

private:
    double factor_at(std::size_t band, double log_s) const;

    std::vector<double> frequencies_;
    std::vector<double> rigid_rates_; // shadow_attenuation_rates(), when rigid
    std::vector<double> log_pi_f_; // log(pi f) per band
    std::vector<std::vector<double>> table_; // [band][i]: Im(b / w) at s_min 10^(i / per_decade)
    std::vector<bool> surface_wave_;
};

// Fill `loss` (dB per band, bands in [first_band, last_band)) for a ray
// `ray.shadow_depth` metres into the shadow: the creeping-wave decay, capped at
// `floor` dB, where energy scattered in by turbulence takes over from it.
void shadow_loss(const RefractedRay& ray, const ShadowAttenuation& attenuation, double floor,
        std::size_t first_band, std::size_t last_band, double* loss);

// The field through the shadow boundary, for a point receiver over ground of a
// given admittance: what replaces the ray's ground term and creeping-wave loss
// there.
//
// Ray theory breaks down within a few diffraction lengths F = 2 k0 l^2 of the
// boundary, on both sides.  Its ground term, taken with the grazing ray's
// straight equivalent geometry, sinks toward the grazing ground wave, 10 dB
// below the exact field over grass at 1 kHz; over rigid ground it doubles
// where the exact field does not.  For k^2 = k0^2 (1 + 2 z / R) the exact field
// is the residue series (upward_pressure).  Relative to the spreading of the
// model's own ray it depends on the signed distance past the boundary in
// diffraction lengths, xi = (r - r_b) / F, and hardly on the source's height:
// within about 1 dB for sources 5 l to 100 l up, from xi = -1.5 on (Fock's
// result, for a source well above the ground).  So each band gets one table,
// V(xi, s), built from the series with the source reference_height l up, and s
// = (pi f / |a|)^(1/3) carries the path as in ShadowAttenuation.  The level is
// the ray's own, ground term included, up to xi = xi_ray; spreading plus V from
// xi_modal on; and a blend between, where the two agree within about a dB.
//
// Past the table, V continues at its last slope, the first mode's.  It is
// held at V(0) - floor at worst, turbulence's scattered floor, as
// --shadow_floor held the creeping-wave loss; and at 20 log10 2 at best, the
// most two coherent rays give, which keeps it under the ground's bound.  A
// band whose s falls below s_min keeps the old model, ground term and
// creeping-wave loss: there the source is a large fraction of R up and the
// profile near the ground is not what the table assumes.
//
// Through layers there are two boundaries: the layered model's r_s, set by a
// log law's steep lowest metres, and the path's linear fit's r_l, beyond it.
// Against the parabolic equation the shadow begins between them, near r_l at
// low frequency and moving toward r_s as the frequency rises, but reaching it
// by no means even at 2 kHz: the wave does not resolve the thin fast layer
// until its diffraction length is short against the gap.  Each band therefore
// takes r_b = r_s + w (r_l - r_s), w = x / (x + boundary_resolution), with x
// the band's diffraction length at the fit's curvature over r_l - r_s, and is
// scaled by the curvature that makes r_b.  Fitted on twelve log laws (u* 0.2 -
// 0.6 m/s, z0 0.01 - 0.1 m, sources 30 and 100 m, 125 Hz - 2 kHz) and held out
// on atmospheres/example_windy.csv: shadow-edge error against the PE 18.9 dB rms
// with r_b = r_s, 6.7 dB with the law.  With r_s = r_l, a linear profile, it is
// the one boundary.
//
// A flush plate microphone gets the field at the plate's top, receiver_height
// its thickness, plus the plate's level over bare ground there, K.  Near
// grazing K barely depends on elevation -- over the GR1425 on the Amedee
// ground, within 0.1 dB from 0.5 to 10 degrees, about 0 dB to 1 kHz, +0.9 at
// 2 kHz, +2.6 at 4 kHz, +7 at 8 kHz -- but it does on azimuth (+-1 dB at 1-4 kHz,
// +-2.5 at 8 kHz), so it is tabulated per band and azimuth at 2 degrees and
// applied at each pair's own.
class PlateReceiver;
struct GroundModel;

// Turbulent scattering into a refraction shadow (--scattering): the incoherent
// energy that eddies scatter down past the boundary, which the coherent field
// there -- ray theory, the residue series, a PE in still air -- has none of.
// Its level re free field at the slant range R from a source z_s up is
//   S = a_db + 10 lg gamma_t + per_decade_frequency lg(f / 1 kHz)
//       + per_decade_range lg(R / 100 m) + per_decade_height lg(z_s / 100 m),
// gamma_t the HARMONOISE strength (ground_effect.hpp), fitted to the wide-angle
// PE with von Karman phase screens (pe.hpp) deep in the shadows of measured
// upwind profiles, wherever the turbulent field stands 10 dB clear of the
// still one (tests/test_pe_turbulence.cpp holds it to the PE).  Single
// scattering is linear in the medium's variance, hence 10 lg gamma_t, and the
// PE agrees: a third of the strength moves it 4.82 dB.  R and z_s are held to the fit's span, and S to 0 dB.  The
// constants depend on the receiver: for_receiver() picks them.
struct ShadowScattering
{
    double gamma_t = 0.0; // 0: none, and the shadow keeps --shadow_floor
    double a_db = 0.0;
    double per_decade_frequency = 0.0;
    double per_decade_range = 0.0;
    double per_decade_height = 0.0;

    bool on() const { return gamma_t > 0.0; }
    double level(double frequency, double slant, double source_height) const;
    // The fitted constants for a flush microphone (under 5 cm up, or on a
    // plate) or one receiver_height metres up.
    static ShadowScattering for_receiver(double gamma_t, double receiver_height, bool plate);
};

class ShadowTransition
{
public:
    ShadowTransition() = default;
    // A point receiver `receiver_height` up; with `plate` (bound to these
    // bands, over `ground`), a flush microphone on that plate, receiver_height
    // its thickness.  The bands are built on `threads` threads (0: one per
    // hardware thread), as the observer loop's -j.
    ShadowTransition(const std::vector<double>& frequencies,
            const std::vector<std::complex<double>>& admittance, double receiver_height,
            double sound_speed, std::size_t first_band, std::size_t last_band,
            const PlateReceiver* plate = nullptr, const GroundModel* ground = nullptr,
            unsigned threads = 0);

    bool empty() const { return table_.empty(); }
    double receiver_height() const { return receiver_height_; }
    // Scattered energy to add in the shadow (off by default).  With it the
    // coherent field is not held at V(0) - floor: the scattered energy is what
    // the floor stood in for.
    void set_scattering(const ShadowScattering& scattering) { scattering_ = scattering; }
    const ShadowScattering& scattering() const { return scattering_; }
    // The largest 10^(V/10) in `band`'s table, at least 1: what the peak
    // search's bound on the ground's gain has to allow for.
    double max_power_ratio(std::size_t band) const
    {
        return band < max_power_ratio_.size() ? max_power_ratio_[band] : 1.0;
    }

    // For an upward-refracting ray that carries a boundary, per band in
    // [first_band, last_band): the weight on the ray's own ground term (1 on
    // the lit side, 0 from xi_modal on) and the loss to take off (dB).  Bands
    // too far from the boundary get weight 1 and, in the shadow, the
    // creeping-wave loss of `attenuation`, which is also what a band below
    // s_min keeps.  False, with nothing written, when the ray is nowhere near
    // a boundary in any band: then the caller's ordinary path stands.
    // `plate_azimuth` is the horizontal propagation direction in the plate's
    // frame, as PlateReceiver takes it; ignored for a point receiver.
    bool apply(const RefractedRay& ray, const ShadowAttenuation& attenuation, double floor,
            std::size_t first_band, std::size_t last_band, double* ground_weight,
            double* loss, double plate_azimuth = 0.0) const;

    // K, dB, in `band` at `azimuth` (radians, the plate's frame); 0 for a point
    // receiver (tests).
    double plate_offset(std::size_t band, double azimuth) const;

    // One band's effective boundary between the layered r_s and the linear
    // fit's r_l > r_s (see above), for a source and receiver z_s and z_r up;
    // `flush` for a plate or flush microphone.  Where the fit casts no shadow
    // (r_l infinite) there is nothing to interpolate toward, and the band takes
    // a fitted multiple of r_s instead.
    static double band_boundary(double r_s, double r_l, double z_s, double z_r, double frequency,
            double sound_speed, bool flush = false);

    // V, dB, in `band` at xi and s: the table (tests).  NaN where it has no
    // value, which is only on the lit side, below xi_modal.
    double level(std::size_t band, double xi, double s) const;

    static constexpr double xi_ray = -2.0;
    static constexpr double xi_modal = -1.0;
    static constexpr double xi_max = 6.0;
    static constexpr double xi_step = 0.1;
    static constexpr double s_min = 5.0;
    static constexpr double s_max = 3000.0;
    static constexpr int s_per_decade = 8;
    // The table's source height in l, or 0.1 s^2 where that is lower, which
    // keeps the source below R / 10.
    static constexpr double reference_height = 20.0;
    // Where a band's diffraction length at the fit's curvature equals this
    // fraction of r_l - r_s, its boundary is half-way between them.  Fitted for
    // a microphone 1.2 m up (0.387 refitted with day-284 profiles added); a
    // flush one resolves the thin fast layer later, and wants about twice it:
    // 0.47 on log laws, 0.94 on the measured day-284 profiles.
    static constexpr double boundary_resolution = 0.363;
    static constexpr double boundary_resolution_flush = 0.75;
    // Where the fit casts no shadow: r_b = m (f / 1 kHz)^fallback_per_decade r_s,
    // the per-band optimum against the PE pooled over the calibration cases.
    static constexpr double fallback_multiple_flush = 1.60;
    static constexpr double fallback_multiple_raised = 1.90;
    static constexpr double fallback_per_decade = -0.08;
    // Receivers below this are flush for the purposes above, m.
    static constexpr double flush_height = 0.05;

private:
    double column(std::size_t band, std::size_t j, double xi) const;

    std::vector<double> frequencies_;
    double receiver_height_ = 0.0;
    std::size_t xi_nodes_ = 0;
    std::size_t s_nodes_ = 0;
    // [band][j * xi_nodes_ + i]: V at s_j, xi_i.
    std::vector<std::vector<double>> table_;
    // [band][j]: V at xi = 0, for the floor.
    std::vector<std::vector<double>> at_boundary_;
    std::vector<double> max_power_ratio_;
    // [band][i]: K at azimuth i * plate_azimuth_step; empty for a point receiver.
    std::vector<std::vector<double>> plate_offset_;
    static constexpr int plate_azimuths = 36;
    ShadowScattering scattering_;
};

} // namespace niceops
