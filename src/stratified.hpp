#pragma once

// Refraction through the profile as given: step 6 of the propagation roadmap
// (--refraction stratified).
//
// --refraction linear fits one straight line to the effective sound speed of
// each path.  That averages away exactly what shapes near-ground propagation:
// a log-law wind's steep lowest meters, an inversion.  This takes the effective
// sound speed c(z) + w(z).h as piecewise linear between the profile's own
// levels, held constant below the lowest and above the highest.  Within each
// layer a ray is still a circular arc, so crossing a layer has a closed form,
// and a ray through the whole profile is a sum over its layers.
//
// Rays are parametrized by launch elevation theta.  In a horizontally
// stratified medium the Snell invariant p = cos(theta) / c(z_s) fixes the
// whole ray.  Two branches reach a receiver below the source from above:
//   A  launched downward, reaching the receiver before turning;
//   B  launched upward, turning where c(z) = 1/p, and coming back down.
// When the sound speed below the source never exceeds its value at the source,
// the two join at horizontal launch into one continuous family.  When it does
// -- upward refraction somewhere below -- A ends at the critical angle, whose
// ray turns without reaching the receiver.  B, if the air aloft is fast
// enough, starts beyond, and between them lies a shadow.  The range is taken to
// be monotonic along each branch.  It is for monotonic profiles; where it is
// not (caustics), the solve finds one of the arrivals, not all of them.
//
// The ray tube -- the spreading -- comes from finite differences of the range
// in launch angle:
//   r_eq^2 = r |dr/dtheta| sin|arrival| / cos(launch),
// which for a linear profile is exactly (r / cos(launch))^2.  Where neighboring
// rays cross, at a caustic, it goes to zero.  The focusing it implies is capped
// at max_focusing_gain_db, a stand-in for the Airy function's finite peak.
//
// As with --refraction linear, this is the effective-sound-speed approximation:
// rays stay in the vertical plane of the bearing.

#include <cstddef>
#include <limits>
#include <vector>

#include "atmosphere_profile.hpp"
#include "refraction.hpp"
#include "util/constants.hpp"

namespace niceops {

// The focusing cap per ray over spherical spreading, dB. Distinct arrivals
// add in energy, so the summed field can exceed this per-ray cap.
inline constexpr double max_focusing_gain_db = 10.0;

// The highest frequency the layered solve must answer for, Hz, unless told
// otherwise: the top third-octave band of every shipped database.  A run
// passes its own highest band (LayeredSoundSpeed's constructor).
//
// It sets the finest vertical structure the rays are allowed to see at the
// fastest air.  Ray theory sends the limiting ray of a ground-grazing shadow
// along a flat or nearly flat stretch of that air without limit -- a calm,
// isothermal floor 3 cm thick put the shadow boundary 20 km out, and the pairs
// short of it got near-grazing rays with tubes tens of times their range,
// 22-28 dB low -- where a wave confined to a stretch h tall diffracts out of
// it within about k h^2.  So a downward ray is not followed past the invariant
// at which its horizontal travel through the limiting ray's turning region
// (the heights within one radian of vertical phase of the fastest air) would
// exceed 1.5 k h^2, k at this frequency.  In any linear gradient ray theory
// itself travels 8/9 of that there, so a resolved profile never reaches it;
// the highest band gives the largest k, so the limit removes only what no
// band of the run resolves (resolved_pass in stratified.cpp).
inline constexpr double default_resolution_frequency = 10000.0; // Hz

// The profile's levels and, per level, the sound speed and local-frame wind:
// the effective sound speed along any bearing is c + hx w_x + hy w_y there.
class LayeredSoundSpeed
{
public:
    LayeredSoundSpeed() = default;
    // Wind turned into the local frame whose +x axis points along compass
    // bearing `frame_bearing`, degrees, as EffectiveSoundSpeed does.
    // `resolution_frequency`: the highest frequency the run evaluates, Hz
    // (default_resolution_frequency).
    LayeredSoundSpeed(const AtmosphereProfile& profile, double frame_bearing,
            double resolution_frequency = default_resolution_frequency);

    std::size_t levels() const
    {
        return heights_.size();
    }
    const std::vector<double>& heights() const
    {
        return heights_;
    }
    // The effective sound speed at each level along horizontal direction
    // (hx, hy), written into `out` (room for levels()).
    void along(double hx, double hy, double* out) const;
    double resolution_frequency() const
    {
        return resolution_frequency_;
    }

private:
    std::vector<double> heights_;
    std::vector<double> sound_speed_;
    std::vector<double> wind_x_;
    std::vector<double> wind_y_;
    double resolution_frequency_ = default_resolution_frequency;
};

// One path's effective sound speed: piecewise linear between `levels` heights
// `z` (ascending), constant outside them.
struct Layers
{
    const double* z = nullptr;
    const double* c = nullptr;
    std::size_t levels = 0;
    // The highest frequency the rays answer for, Hz: it limits the vertical
    // structure they see at the fastest air (default_resolution_frequency).
    // 0 is plain ray theory, followed along flat air without limit.
    double resolution_frequency = default_resolution_frequency;

    // `upper`, if known, is the first level above `height` (std::upper_bound's
    // answer), so that a walk that has found it does not search again.
    double at(double height, const double* upper = nullptr) const;
    // Whether the effective sound speed is the same everywhere on the path's
    // heights, to rounding: then every ray is straight.
    bool uniform() const;
};

// The direct and ground-reflected eigenrays in a layered profile, in the
// linear model's terms (refraction.hpp).  `guess` is a launch angle to start
// the solve from -- the linear model's, usually within a fraction of a degree
// -- and matters only for speed.
DirectEigenray layered_direct_eigenray(const Layers& layers, double z_s, double z_r, double r,
        double guess);
ReflectedEigenray layered_reflected_eigenray(const Layers& layers, double z_s, double z_r,
        double r, double guess);

// Every ray that reaches the receiver from above, but the primary: the direct
// rays of the families after the one layered_direct_eigenray takes (`primary`
// is its result, which counts only if lit), and, in a surface duct, the rays
// that bounce off the ground first.  With `bounced_only`, the dipping family's
// direct rays are dropped (the linear fit stands for them); the bounced rays
// and those turned back from aloft are kept.
//
// A duct is air above the source faster than all of it below, so that a ray
// that reaches the ground turns aloft and comes back down.  A ray of
// invariant p that bounces m times before its last descent travels
//   X_m(p) = X_0(p) + m L(p),
// X_0 its direct range and L the hop from the ground up to its turning height
// and back.  Only rays arriving from above are listed: the one rising to the
// receiver from the next bounce is the partner the ground model pairs each
// with, as the primary's reflected ray is the primary's.  X_0 rises with p
// where L falls, so X_m is not monotonic and has caustics: the arms are
// sampled once, together, in launch angle (the hops reuse the samples) and
// every crossing of r is solved for.  A ray whose hops alone already outrun r
// is skipped without a trace, which leaves nothing to do for a source above
// the duct or a path shorter than its first hop.
//
// Past the reach of every direct ray in a duct -- the single-hop skip
// distance, ~2.4 km for a 30 m source under day 286's morning inversion --
// these are the only rays there are.  Leaving them out, the level sat 8-54 dB
// under the parabolic equation's there (which stays within a few dB of free
// field) and fell into the gaps between hops.
//
// The launch angles are sampled 32 times, plus 8 for each bounce the path can
// hold (at most 512 in all), bunched toward the shallow end.  `samples`, if
// positive, replaces that count, and exists so that a test can ask whether the
// answer has converged in it. All solved arrivals are retained.
std::vector<Arrival> layered_arrivals(const Layers& layers, double z_s, double z_r, double r,
        const DirectEigenray& primary, bool bounced_only = false, int samples = 0);

// Two wave effects the rays above leave out, each behind its own option
// (WaveCorrection, refraction.hpp, gives the loss):
//
// tunneling: a ray turned back from aloft that turns again just above the
//   receiver never reaches it, and ray theory leaves the receiver in silence.
//   The wave reaches it as the field's tail below the turning point.  Under day
//   286's last profile the turning height is 1-20 ft, inside the log-law wind
//   assumed below the lidar's first gate, and the parabolic equation hardly
//   notices the layer (<1 dB to 250 Hz): the rays that turn lowest carry the
//   duct's return to the ground.  Listed: rays launched up, and rays launched
//   down that have turned aloft at least once, that turn at z_t above z_r after
//   m = 0, 1, ... cycles between their two turning heights, where the range to
//   that turning point is r.  Rays launched down that have not turned aloft are
//   left to ShadowTransition, whose field past a ground-grazing boundary is
//   theirs.  Their tube is taken a reference Airy length above the turning
//   point (at wave_reference_frequency), where the ray is still a ray.
//
// skip_fill: short of where a family from aloft first lands -- a caustic, where
//   its range X(p) is least -- is a skip zone that no ray reaches, and the
//   level there fell to the creeping-wave floor or to nothing.  Each interior
//   extremum of a family's range with the receiver on its dark side adds the
//   caustic's Airy tail: I = C Ai^2(x), C = 4 pi omega^(1/3) /
//   (2^(1/3) G X''^(2/3)), G = r |dp/dtheta| sin|arrival| / cos(launch), which
//   matches the two rays' summed energy on the lit side.  The families are the
//   ones from aloft (direct, bounced along a duct, and with tunneling the
//   tunneled ones); a downward family's interior extrema are a log law's, and
//   the hand-over to the linear fit stands for those.
struct WaveOptions
{
    bool tunneling = false;
    bool skip_fill = false;

    bool any() const { return tunneling || skip_fill; }
};

std::vector<Arrival> layered_wave_arrivals(const Layers& layers, double z_s, double z_r, double r,
        const WaveOptions& options);

// --shadow_model resolved: what each band's shadow boundary is set by.
//
// Ray theory sets a ground-grazing shadow boundary by c(z_s) / c(0) alone, so
// the log law's lowest foot -- 0.25 m/s on the R44's light-wind mornings,
// extrapolated below the lidar's first gate -- casts it.  A wave does not
// resolve a layer that gathers less than about a radian of vertical phase, so
// the band at frequency f sees, below the height h_f where
//   k_f int_0^h sqrt(2 max(c(0) - c(z), 0) / c(0)) dz = phase,
// only the profile's tangent at h_f (resolved_height, tangent_below).  The
// radian is the one resolved_pass already takes for flat floors; on a linear
// profile h_f = (3 phase / (2 k sqrt(2 |a| / c0)))^(2/3), which at one radian
// is 1.31 l, l = (c0 / (2 |a| k^2))^(1/3) the creeping-wave layer thickness.
inline constexpr double default_resolved_phase = 1.0; // rad

// h_f for a band at `frequency` on these layers (above), meters.  Infinite
// where the integral never reaches `phase`: the band resolves no upward
// refraction near the ground.
double resolved_height(const Layers& layers, double frequency, double phase);

// The layers with the effective sound speed below `height` replaced by its
// tangent there -- the slope of the layer containing it -- written into z and
// c (levels 0, height, then the levels above).  `gradient` gets that slope.
Layers tangent_below(const Layers& layers, double height, std::vector<double>& z,
        std::vector<double>& c, double& gradient);

// The ground-grazing shadow boundary layered_direct_eigenray sets on these
// layers for a source z_s above a receiver z_r: infinite where it sets none.
double layered_shadow_boundary(const Layers& layers, double z_s, double z_r);

// One band's shadow under --shadow_model resolved: ResolvedBand (refraction.hpp).


// Every band in [first_band, last_band) for one path: `layers` its effective
// sound speed, `fit` its linear fit over [min(z_s, z_r), max(z_s, z_r)]
// widened to EffectiveSoundSpeed::min_fit_span (as EffectiveSoundSpeed::fit),
// whose boundary is corrected by the tangent's change to the profile for each
// band's r_l(f) (or, without `resolve_linear`, the fit's own boundary for
// every band: a sensitivity).  `flush` as for band_boundary.  Writes out[band].
void resolved_shadow_bands(const Layers& layers, const LinearProfile& fit, double z_s, double z_r,
        const std::vector<double>& frequencies, std::size_t first_band, std::size_t last_band,
        double phase, bool flush, ResolvedBand* out, bool resolve_linear = true, bool boundary_law = true);

// resolved_shadow_bands over a run: per receiver height, a table over bearing
// (every bearing_step degrees of the local frame's
// direction, counterclockwise from +x) and source height (nodes even in
// sqrt(z)), interpolated bilinearly in 1/r_b (so a band with no boundary is 0)
// and in the gradient.  Built once per trajectory; a pair's lookup is a few
// dozen multiplies per band.
class ResolvedShadowTable
{
public:
    ResolvedShadowTable() = default;
    // `air`/`linear` as layered_direct_ray takes them; source heights over
    // [z_low, z_high]; receivers at `receiver_heights` (flush below
    // ShadowTransition::flush_height, or every one with `plate`).
    ResolvedShadowTable(const LayeredSoundSpeed& air, const EffectiveSoundSpeed& linear,
            const std::vector<double>& frequencies, std::size_t first_band, std::size_t last_band,
            double phase, double z_low, double z_high, const std::vector<double>& receiver_heights,
            const std::vector<bool>& flush, unsigned threads = 0);

    bool empty() const { return bands_ == 0; }
    // The bands for a path along local-frame direction (hx, hy) from a source
    // z_s up to receiver `receiver` (an index into the constructor's heights).
    void lookup(double hx, double hy, double z_s, std::size_t receiver, ResolvedBand* out) const;

    static constexpr double bearing_step = 2.0; // degrees
    static constexpr int height_nodes_per_meter_root = 4; // nodes per unit of sqrt(z / 1 m)

private:
    std::size_t first_band_ = 0, bands_ = 0;
    std::size_t bearings_ = 0, heights_ = 0;
    double root_low_ = 0.0, root_step_ = 1.0;
    // [receiver][(bearing * heights_ + height) * bands_ + band]: 1 / r_b and the gradient.
    std::vector<std::vector<float>> inverse_boundary_;
    std::vector<std::vector<float>> gradient_;
};

// --refraction stratified's counterparts of direct_ray() and finish_ground().
// `linear` supplies the warm start; `scratch` holds the path's effective sound
// speed at each level between the two calls, so is per thread.
//
// With `linear_past_shadow`, a pair beyond handover_fraction of the layered
// boundary whose linear fit's boundary lies further out gets the fit's ray
// instead, lit or in the fit's own shadow, carrying both boundaries
// (from_linear_fit set).  That is for ShadowTransition, whose bands' boundaries
// lie between the two: past the layered one the lit level is the fit's.
// Formed here, not later, so that the peak search bounds the ray whose level
// it will take.
//
// handover_fraction: against the parabolic equation on the lit side, over the
// last half of the way to the layered boundary, the layered ray is 6-9 dB out
// at worst and the fit's 1-3; handing over at 0.9 of the boundary is 0.6 dB rms
// on measured day-284 profiles (1.4 at the boundary itself) and 1.7 on log laws
// (3.3).
inline constexpr double handover_fraction = 0.9;

// What a pair's layered solves give -- the direct eigenray, the hand-over,
// the arrivals (ray and wave), and for a point receiver the reflected ray's
// grazing angle and the length the ground term takes it over -- kept, so that
// a state between two solved ones can take it interpolated
// (interpolate_layered) and skip the solves (--ray_stride).
struct LayeredSolution
{
    bool solved = false;
    // The pair's straight distance, source to microphone: interpolation takes
    // the ray tubes, path lengths and travel times as fractions of it, so
    // that a tube's growth with range is exact and only its focusing is
    // interpolated.
    double chord = 0.0;
    DirectEigenray direct;
    bool handover = false;
    std::vector<Arrival> arrivals;
    std::vector<Arrival> waves;
    bool ground_set = false;
    double grazing = 0.0;
    double length = 0.0;
};

// `a` and `b` (two states' solutions) at t in (0, 1) between them, for a
// pair whose straight distance is `chord`, into `out`: angles and boundaries
// linear in t, tubes, lengths and times linear as fractions of the chord.
// False, and `out` left as it may be, unless both were solved and alike: the
// same lit or shadowed primary, boundary or none, hand-over, and arrivals of
// the same bounces and wave kinds, every launch and arrival angle within
// ray_interpolation_angle of its partner, and every tube, as a fraction of
// its chord, within ray_interpolation_spreading -- a caustic, a family's
// switch or a boundary passing makes them differ, and the states between are
// solved.  The ground's part only where both have it, the reflected rays
// alike.
bool interpolate_layered(const LayeredSolution& a, const LayeredSolution& b, double t, double chord,
        LayeredSolution& out);
bool interpolate_ground(const LayeredSolution& a, const LayeredSolution& b, double t, double chord,
        LayeredSolution& out);
inline constexpr double ray_interpolation_angle = 0.5 * deg_to_rad; // rad
inline constexpr double ray_interpolation_spreading = 0.02;          // relative

// A pair's straight distance as LayeredSolution::chord takes it: the
// microphone held on the plane, as the layered solve is.
double layered_chord(const Vec3& source, const Vec3& receiver);

// `memo`: a pair's solves' results.  Where it holds them (`solved`, and
// `ground_set` for the ground's part) they are taken instead of solving --
// the direct ray's only if this pair hands over as they did; where it does
// not, the solves fill it.
RefractedRay layered_direct_ray(const LayeredSoundSpeed& air, const EffectiveSoundSpeed& linear,
        const Vec3& source, const EffectiveSoundSpeed::Moments& at_source, const Vec3& receiver,
        const EffectiveSoundSpeed::Moments& at_receiver, std::vector<double>& scratch,
        bool linear_past_shadow = false, const WaveOptions& waves = {}, LayeredSolution* memo = nullptr);
void layered_finish_ground(const LayeredSoundSpeed& air, RefractedRay& ray, GroundGeometry ground,
        double microphone_height, std::vector<double>& scratch, LayeredSolution* memo = nullptr);

// --ray_stride's default.
inline constexpr std::size_t default_ray_stride = 4;

// --ray_stride: one observer's layered rays, solved at its anchor states --
// every stride-th, and the last -- and between two anchors interpolated from
// their solves where they agree (interpolate_layered, interpolate_ground).
// Each anchor keeps the unfinished ray its ground finish starts from.
// Anchors are fixed states, each solved alike whichever state asks first, so
// a ray does not depend on the order the states are asked for: a windowed
// run and a marched one take the same.  A stride of 1, or a trajectory of two
// states or fewer, solves every state.  One per thread.
class StridedLayeredRays
{
public:
    // Every state's source and its moments, which must outlive this.
    StridedLayeredRays(const LayeredSoundSpeed& air, const EffectiveSoundSpeed& linear,
            const std::vector<Vec3>& sources, const std::vector<EffectiveSoundSpeed::Moments>& moments,
            std::size_t stride, const WaveOptions& waves);
    // The observer the rays end at from now on; the last one's anchors go.
    // `linear_past_shadow` as layered_direct_ray takes it.
    void begin_observer(const Vec3& receiver, const EffectiveSoundSpeed::Moments& at_receiver,
            bool linear_past_shadow);
    // State it's ray, as layered_direct_ray gives it.
    RefractedRay ray(std::size_t it);
    // Its ground finish, as layered_finish_ground's: an anchor's from its own
    // solve, an interpolated ray's from its anchors' where they agree.
    void finish(std::size_t it, RefractedRay& ray, GroundGeometry ground, double microphone_height);

private:
    struct Anchor
    {
        bool formed = false;
        bool finished = false;
        LayeredSolution solution;
        RefractedRay ray;
    };
    // The anchors at or before `it` and after it (one or both `it` itself
    // when it is one), and its fraction of the way between them.
    void block_of(std::size_t it, std::size_t& a, std::size_t& b, double& t) const;
    Anchor& slot(std::size_t anchor) { return anchors_[(anchor + stride_ - 1) / stride_]; }
    Anchor& formed(std::size_t anchor);
    RefractedRay solve(std::size_t it, LayeredSolution* memo);

    const LayeredSoundSpeed& air_;
    const EffectiveSoundSpeed& linear_;
    const std::vector<Vec3>& sources_;
    const std::vector<EffectiveSoundSpeed::Moments>& moments_;
    std::size_t stride_;
    WaveOptions waves_;
    Vec3 receiver_;
    EffectiveSoundSpeed::Moments at_receiver_ {};
    bool linear_past_shadow_ = false;
    std::vector<Anchor> anchors_;
    std::vector<double> scratch_;
};

} // namespace niceops
