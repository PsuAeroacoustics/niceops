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
#include <vector>

#include "atmosphere_profile.hpp"
#include "refraction.hpp"

namespace niceops {

// The most a caustic may focus a ray over spherical spreading, dB.
inline constexpr double max_focusing_gain_db = 10.0;

// The profile's levels and, per level, the sound speed and local-frame wind:
// the effective sound speed along any bearing is c + hx w_x + hy w_y there.
class LayeredSoundSpeed
{
public:
    LayeredSoundSpeed() = default;
    // Wind turned into the local frame whose +x axis points along compass
    // bearing `frame_bearing`, degrees, as EffectiveSoundSpeed does.
    LayeredSoundSpeed(const AtmosphereProfile& profile, double frame_bearing);

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

private:
    std::vector<double> heights_;
    std::vector<double> sound_speed_;
    std::vector<double> wind_x_;
    std::vector<double> wind_y_;
};

// One path's effective sound speed: piecewise linear between `levels` heights
// `z` (ascending), constant outside them.
struct Layers
{
    const double* z = nullptr;
    const double* c = nullptr;
    std::size_t levels = 0;

    double at(double height) const;
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
// that bounce off the ground first.  With `bounced_only`, only those.
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
// answer has converged in it.
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

RefractedRay layered_direct_ray(const LayeredSoundSpeed& air, const EffectiveSoundSpeed& linear,
        const Vec3& source, const EffectiveSoundSpeed::Moments& at_source, const Vec3& receiver,
        const EffectiveSoundSpeed::Moments& at_receiver, std::vector<double>& scratch,
        bool linear_past_shadow = false, const WaveOptions& waves = {});
void layered_finish_ground(const LayeredSoundSpeed& air, RefractedRay& ray, GroundGeometry ground,
        double microphone_height, std::vector<double>& scratch);

} // namespace niceops
