#pragma once

// LoD1 buildings for the urban tables' general builder (urban_paths.hpp): each
// a footprint polygon extruded from the floor, z = 0, to a flat roof.  What
// the builder needs of them is here: their walls as facades (urban.hpp), the
// edges they diffract at, and whether a straight segment passes through one.
//
// The floor is a mirror for the builder, as for the ideal canyon: a ray the
// floor reflects is traced straight on into the floor's image of the world,
// where every building stands from -top to top.  `below_floor` selects that
// world; above the floor the two are the same.
//
// The visibility test is exact for extruded footprints, which is why the
// builder needs no ray-tracing library at LoD1: a segment is blocked when some
// stretch of it longer than `contact` lies strictly inside a building.  A
// segment that only touches a building -- ends on a wall where it reflects,
// starts on an edge where it diffracts, or grazes along a wall -- is not.

#include <cstddef>
#include <vector>

#include "urban.hpp"
#include "vec3.hpp"

namespace niceops {

// A footprint: its vertices in the run's local frame (m), in either
// orientation, the first not repeated at the end, and its roof's height above
// the floor.
struct Footprint
{
    std::vector<double> x;
    std::vector<double> y;
    double top = 0.0;
};

// A straight edge where a building diffracts: the roof edge along each side of
// its footprint, and the vertical edge at each convex corner.  The edge runs
// from `start` for `length` along the unit `direction`.  In the plane normal
// to it the wedge's first face leaves the edge along the unit `face0`, the air
// lies toward the unit `air` (normal to both the edge and face0), and the
// second face is at the angle n pi from the first, measured through the air:
// n = 3/2 for a roof edge or a right-angled corner (wedge.hpp).
struct Edge
{
    Vec3 start;
    Vec3 direction;
    double length = 0.0;
    Vec3 face0;
    Vec3 air;
    double n = 1.5;
    std::size_t building = 0;
    bool vertical = false;

    // The angle of `v` about the edge, from face0 through the air, in
    // [0, 2 pi): the wedge's own angle (wedge.hpp).
    double angle(const Vec3& v) const;
};

class Buildings
{
public:
    // Distances under this are contact, not passage, in the visibility test
    // and in the tests for lying on a wall (m).
    static constexpr double contact = 1e-6;

    Buildings() = default;

    // Throws std::invalid_argument for a footprint with fewer than three
    // vertices, a zero or negative area or top, or sides that cross.
    explicit Buildings(std::vector<Footprint> footprints);

    std::size_t size() const
    {
        return footprints_.size();
    }

    // The footprints, each counterclockwise.
    const std::vector<Footprint>& footprints() const
    {
        return footprints_;
    }

    // Every side of every footprint as a facade standing from the floor to the
    // roof, its base direction clockwise around the building so that the
    // outside is on its left (facade_side > 0).
    const std::vector<Facade>& facades() const
    {
        return facades_;
    }

    // The building each facade belongs to.
    std::size_t facade_building(std::size_t facade) const
    {
        return facade_building_[facade];
    }

    const std::vector<Edge>& edges() const
    {
        return edges_;
    }

    // The tallest roof.
    double max_top() const
    {
        return max_top_;
    }

    // Whether the segment from a to b passes through a building: some stretch
    // longer than `contact` lies strictly inside one.
    bool blocked(const Vec3& a, const Vec3& b, bool below_floor = false) const;

    // Whether a ray from `a` along `direction` reaches infinity without
    // passing through a building.  A ray that does not rise escapes once past
    // every footprint.
    bool escapes(const Vec3& a, const Vec3& direction, bool below_floor = false) const;

    // Whether (x, y) lies strictly inside a footprint, `contact` or more from
    // its sides; the building's index, or size() for none.
    std::size_t inside(double x, double y) const;

private:
    bool blocked_by(std::size_t building, const Vec3& a, const Vec3& b, bool below_floor) const;
    bool strictly_inside(std::size_t building, double x, double y) const;

    std::vector<Footprint> footprints_;
    std::vector<Facade> facades_;
    std::vector<std::size_t> facade_building_;
    std::vector<Edge> edges_;
    // Per building: x min, x max, y min, y max.
    std::vector<double> bounds_;
    double max_top_ = 0.0;
    double x_min_ = 0.0, x_max_ = 0.0, y_min_ = 0.0, y_max_ = 0.0;
};

} // namespace niceops
