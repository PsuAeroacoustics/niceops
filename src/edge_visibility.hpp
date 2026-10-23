#pragma once

// Which stretch of an edge can be seen from a point, along the folded path
// folded_path_open() (building_paths.hpp) checks one point of it at a time:
// for pruning the point-source search, where a point of an edge a source's or
// microphone's image cannot see diffracts nothing from it.
//
// The answer errs one way only.  Every point of the edge outside [from, to]
// is proved hidden, and false means every point is, so a pruned point is one
// folded_path_open would refuse.  The proof: unfolded, the
// path from the edge's point q to the end is the straight segment from q to
// the end's image (its mirror in the facade, and in the floor for a floor
// bounce); a building blocks the path's leg in front of the facade where that
// segment's front part passes through it, and the leg behind where its back
// part passes through the building's mirror.  For a convex building, the
// points q of the edge whose segment passes through it form one interval of
// the edge (a linear-fractional programme over the building's faces), so the
// hidden points are those the intervals of some buildings cover, and those
// outside the stretch whose path meets the facade at all; from and to are
// the first and last point left.  Buildings are shrunk by a margin
// first, which makes every such passage one blocked() counts; a footprint
// that is not convex goes in as its convex pieces (Buildings::convex_pieces),
// each under the building's own top and roof and so within its solid, and one
// with courtyards is left out, either of which can only leave a point
// visible.  For a chain of facades the path is not one segment, so the
// proof goes leg by leg in the real world instead: the leg to facade k lies on
// the line from q, mirrored in the facades after k, to the end's image through
// the facades up to k, and between the planes of the facades either side of
// it, in front of both.  Each facade's window is the one-facade window on its
// leg's line, and each leg's buildings are the real ones cut to that slab, so
// the intervals are found the same way.

#include <cstddef>
#include <vector>

#include "buildings.hpp"
#include "vec3.hpp"

namespace niceops {

bool visible_stretch(const Buildings& buildings, const Edge& edge, const Vec3& end,
        const std::vector<std::size_t>& facades, bool below_floor, double& from, double& to);

} // namespace niceops
