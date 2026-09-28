#pragma once

// Refraction through the profile as given: step 6 of the propagation roadmap
// (--refraction stratified).
//
// --refraction linear fits one straight line to the effective sound speed of
// each path.  That averages away exactly what shapes near-ground propagation:
// a log-law wind's steep lowest metres, an inversion.  This takes the effective
// sound speed c(z) + w(z).h as piecewise linear between the profile's own
// levels, held constant below the lowest and above the highest.  Within each
// layer a ray is still a circular arc, so crossing a layer has a closed form,
// and a ray through the whole profile is a sum over its layers.
//
// Rays are parametrised by launch elevation theta.  In a horizontally
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
// which for a linear profile is exactly (r / cos(launch))^2.  Where neighbouring
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

// --refraction stratified's counterparts of direct_ray() and finish_ground().
// `linear` supplies the warm start; `scratch` holds the path's effective sound
// speed at each level between the two calls, so is per thread.
//
// With `linear_past_shadow`, a pair in the layered shadow whose linear fit is
// still lit gets the fit's ray instead, carrying both boundaries
// (from_linear_fit set).  That is for ShadowTransition, whose bands' boundaries
// lie between the two: past the layered one the lit level is the fit's.
// Formed here, not later, so that the peak search bounds the ray whose level
// it will take.
RefractedRay layered_direct_ray(const LayeredSoundSpeed& air, const EffectiveSoundSpeed& linear,
        const Vec3& source, const EffectiveSoundSpeed::Moments& at_source, const Vec3& receiver,
        const EffectiveSoundSpeed::Moments& at_receiver, std::vector<double>& scratch,
        bool linear_past_shadow = false);
void layered_finish_ground(const LayeredSoundSpeed& air, RefractedRay& ray, GroundGeometry ground,
        double microphone_height, std::vector<double>& scratch);

} // namespace niceops
