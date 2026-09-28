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
    // The gradient (1/s, negative) the diffraction near the boundary is scaled
    // by: for --refraction linear, the profile's own.  0 where there is no
    // boundary, or where ShadowTransition should not be used -- which for now
    // is --refraction stratified (see layered_direct_eigenray).
    double boundary_gradient = 0.0;
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
    double horizontal = 0.0;
    DirectEigenray direct;
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
class ShadowTransition
{
public:
    ShadowTransition() = default;
    ShadowTransition(const std::vector<double>& frequencies,
            const std::vector<std::complex<double>>& admittance, double receiver_height,
            double sound_speed, std::size_t first_band, std::size_t last_band);

    bool empty() const { return table_.empty(); }
    double receiver_height() const { return receiver_height_; }
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
    bool apply(const RefractedRay& ray, const ShadowAttenuation& attenuation, double floor,
            std::size_t first_band, std::size_t last_band, double* ground_weight,
            double* loss) const;

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
};

} // namespace niceops
