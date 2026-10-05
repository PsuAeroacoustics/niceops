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
//                    the source's center: the level falls by
//                    20 log10(range / r_sphere).  For a curved ray it comes
//                    from the ray tube's area, not the arc.
//   path_length      The length absorbed along, from the source's center
//                    (eval_spectrum absorbs over path_length - r_sphere): the
//                    slant range for a straight ray, the arc for a curved one.
//   source_height,   The ground model's geometry (GroundTable, the exact model,
//   ground_distance  the plate receiver).  A curved ray supplies the
//                    straight-ray geometry that has its grazing angle and path
//                    length, so the ground models stay as they are.  Over
//                    terrain (terrain_ray) they are taken over the path's mean
//                    ground line instead of the plane z = 0.
//   receiver_height  The point receiver's height above the ground, for the
//                    ground model; NaN, the default, is the model's own
//                    (--receiver_height, or the observer file's), and only that
//                    height can be read from its GroundTable.  terrain_ray sets
//                    it to the microphone's height above the mean ground line.
//   separation       The straight displacement, source to observer.  What
//                    depends on the bearing alone reads it -- the plate
//                    microphone's orientation -- since a ray in a stratified
//                    atmosphere keeps its bearing unless a crosswind deflects
//                    it, and the effective-sound-speed model has no crosswind.
//                    terrain_ray turns it round where the source's foot on the
//                    mean ground line lies behind the microphone's.
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

#include <algorithm>
#include <cmath>
#include <limits>

#include "terrain.hpp"
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
    double receiver_height = std::numeric_limits<double>::quiet_NaN();
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

// The straight line from `source` to `microphone` over terrain (--terrain): the
// direction, range and path length are the straight line's, through the real
// three-dimensional geometry, and only the ground model's geometry changes.
// It is taken over `line`, the least-squares mean ground line along the path's
// ground projection (GroundPlaneCache, from the observer toward the source's
// footprint): the source's height above the line and the distance between the
// two feet along it (plane_geometry): the mean ground plane of CNOSSOS-EU
// (Directive (EU) 2015/996, Annex II 2.5.6), which NORAH2 and ECAC Doc 32 use
// over terrain.
//
// A point receiver (`flush` false) also gets its own height above the line,
// which is no longer its pole's height wherever the ground at its foot sits
// off the mean line.  A flush microphone on a plate keeps its plate: the plate
// lies along the line at the microphone, so its source height is the source's
// height above the microphone measured perpendicular to the line.  A height
// that comes out below the line -- a receiver in a dip, or a source behind
// rising ground, which is screening and is not modeled -- is taken as zero, as
// CNOSSOS-EU takes a negative equivalent height.  On a level line through z = 0
// the record is straight_ray()'s to the bit, receiver height included.
// `below`, when given, is set to whether a height was raised to zero.
inline RayPath terrain_ray(const Vec3& source, const Vec3& microphone, const GroundLine& line, bool flush,
        bool* below = nullptr)
{
    RayPath path = straight_ray(source, microphone);
    const PlaneGeometry g = plane_geometry(path.ground_distance, source.z, microphone.z, line);
    bool raised = false;
    if (flush)
    {
        const double above_microphone = g.source_height - g.receiver_height;
        raised = above_microphone < 0.0;
        path.source_height = std::max(0.0, above_microphone);
    }
    else
    {
        raised = g.source_height < 0.0 || g.receiver_height < 0.0;
        path.source_height = std::max(0.0, g.source_height);
        path.receiver_height = std::max(0.0, g.receiver_height);
    }
    // The distance along the line is negative when the source's foot on it lies
    // behind the microphone's -- a nearly overhead source over ground falling
    // toward it.  The image construction depends on its square, so its size is
    // the distance (zero would shorten the ground model's direct path below the
    // real range), and the arrival crosses the plate from the other side, so the
    // bearing the plate reads is turned round.
    path.ground_distance = std::fabs(g.ground_distance);
    if (g.ground_distance < 0.0)
    {
        path.separation.x = -path.separation.x;
        path.separation.y = -path.separation.y;
    }
    if (below != nullptr)
        *below = raised;
    return path;
}

} // namespace niceops
