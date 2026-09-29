#pragma once

// The propagation seam: what the level at an observer takes from the path
// between source and observer, as distinct from what it takes from the source.
//
// straight_ray() below is the straight line; refracted_ray() (refraction.hpp)
// is the curved ray of a linear sound-speed profile (--refraction linear).
// Both fill the same record:
//
//   direction        The unit wave normal leaving the source, which is what the
//                    sphere is looked up by (lookup_angles).  For a straight
//                    ray it is the emission-time line to the observer, which is
//                    also how panam's hemigen defined the spheres' angles, so a
//                    refracting model keys the lookup on its launch direction.
//   range            The distance spherical spreading is equivalent to, from
//                    the source's centre: the level falls by
//                    20 log10(range / r_sphere).  For a curved ray it comes
//                    from the ray tube's area, not the arc.
//   path_length      The length absorbed along, from the source's centre
//                    (eval_spectrum absorbs over path_length - r_sphere): the
//                    slant range for a straight ray, the arc for a curved one.
//   source_height,   The ground model's geometry (GroundTable, the exact model,
//   ground_distance  the plate receiver).  A curved ray supplies the
//                    straight-ray geometry that has its grazing angle and path
//                    length, so the ground models stay as they are.
//   separation       The straight displacement, source to observer.  What
//                    depends on the bearing alone reads it -- the plate
//                    microphone's orientation -- since a ray in a stratified
//                    atmosphere keeps its bearing unless a crosswind deflects
//                    it, and the effective-sound-speed model has no crosswind.
//
// What the path absorbs per band is the caller's to supply: eval_spectrum
// multiplies a per-band coefficient by the path.  For a homogeneous atmosphere
// that coefficient is one array for the whole run; through a layered one
// (LayeredAbsorption) it is the mean over the heights the path spans.
//
// Two further things differ between the models and do not live here, because
// moving them in would cost the straight path time or bits:
//   - Reception time.  reception_factor (util/sel.hpp) differentiates the
//     straight-line range; a curved ray supplies its travel time instead, as
//     an equivalent range (sound speed times travel time).
//   - The bounds.  DirectionBounds and the spectral bounds hold the loudest
//     sphere in a direction bin, spread spherically.  A curved ray is bounded
//     by the same tables taken along its launch direction at its ray tube's
//     range, so the peak search and --db_cutoff hand them that geometry
//     instead of the straight displacement (windowed_exposure in niceops.cpp).
//     Absorption, the shadow's loss and the ground's gain are bounded apart.
// The broadband (EAA) path does not use this record at all.  Its windowed
// evaluation shares its geometry with the bound test (level_from_geometry), and
// EAA is itself a straight-line fit.

#include "vec3.hpp"

namespace niceops {

struct RayPath
{
    Vec3 direction;
    double range = 0.0;
    double path_length = 0.0;
    double source_height = 0.0;
    double ground_distance = 0.0;
    Vec3 separation;
};

// The straight line from `source` to `observer`.  Coincident points leave the
// direction as the zero vector rather than dividing by zero.
inline RayPath straight_ray(const Vec3& source, const Vec3& observer)
{
    RayPath path;
    const Vec3 separation = observer - source;
    path.range = separation.norm();
    path.path_length = path.range;
    path.direction = path.range > 0.0 ? separation / path.range : separation;
    path.separation = separation;
    // Height above the observer plane, not above this observer: the ground is
    // the plane z = 0, whatever height the observer stands at.
    path.source_height = source.z;
    path.ground_distance = separation.horizontal_norm();
    return path;
}

} // namespace niceops
