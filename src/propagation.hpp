#pragma once

// The propagation seam: what the level at an observer takes from the path
// between source and observer, as distinct from what it takes from the source.
//
// There is one path model today, the straight line through a homogeneous
// atmosphere, and straight_ray() below is it.  A refracting model returns the
// same record, once per eigenray:
//
//   direction        The unit wave normal leaving the source, which is what the
//                    sphere is looked up by (lookup_angles).  For a straight
//                    ray it is the emission-time line to the observer, which is
//                    also how panam's hemigen defined the spheres' angles, so a
//                    refracting model keys the lookup on its launch direction.
//   range            The distance that spherical spreading and absorption are
//                    both taken over (eval_spectrum), measured from the
//                    source's centre.  A curved ray needs the two apart --
//                    ray-tube spreading on one side, arc length on the other --
//                    and eval_spectrum split to match.
//   source_height,   The ground model's geometry (GroundTable, the exact model,
//   ground_distance  the plate receiver).  A curved ray supplies the
//                    straight-ray geometry that has its grazing angle and path
//                    difference, so the ground models stay as they are.
//
// What the path absorbs per band is the caller's to supply: eval_spectrum
// multiplies a per-band coefficient by the path.  For a homogeneous atmosphere
// that coefficient is one array for the whole run; through a layered one
// (LayeredAbsorption) it is the mean over the heights the path spans.
//
// A refracting model must change two further things, and neither lives here,
// because moving them in would cost the straight path time or bits:
//   - Reception time.  reception_factor (util/sel.hpp) differentiates the
//     straight-line range; a curved ray supplies its travel time instead.
//   - The bounds.  DirectionBounds and the spectral bounds subtract
//     straight-line spreading from the loudest sphere in a direction bin.
//     Refraction can add gain by focusing, and a bound has to cover that.
// The broadband (EAA) path does not use this record at all.  Its windowed
// evaluation shares its geometry with the bound test (level_from_geometry), and
// EAA is itself a straight-line fit.

#include "vec3.hpp"

namespace niceops {

struct RayPath
{
    Vec3 direction;
    double range = 0.0;
    double source_height = 0.0;
    double ground_distance = 0.0;
};

// The straight line from `source` to `observer`.  Coincident points leave the
// direction as the zero vector rather than dividing by zero.
inline RayPath straight_ray(const Vec3& source, const Vec3& observer)
{
    RayPath path;
    const Vec3 separation = observer - source;
    path.range = separation.norm();
    path.direction = path.range > 0.0 ? separation / path.range : separation;
    // Height above the observer plane, not above this observer: the ground is
    // the plane z = 0, whatever height the observer stands at.
    path.source_height = source.z;
    path.ground_distance = separation.horizontal_norm();
    return path;
}

} // namespace niceops
