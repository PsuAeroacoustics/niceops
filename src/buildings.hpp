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

#include "plan_grid.hpp"
#include "urban.hpp"
#include "vec3.hpp"

namespace niceops {

// A ring of a footprint: its vertices, in either orientation, the first not
// repeated at the end.
struct Ring
{
    std::vector<double> x;
    std::vector<double> y;
};

// A footprint: its outer ring's vertices in the run's local frame (m), in
// either orientation, the first not repeated at the end; its roof's height
// above the floor; and the holes it surrounds (courtyards, open to the sky).
// Buildings orients every ring so that the solid is on its left: the outer
// one counterclockwise, the holes clockwise.
struct Footprint
{
    std::vector<double> x;
    std::vector<double> y;
    double top = 0.0;
    std::vector<Ring> holes;
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
    // The second face's unit normal toward the air, sin(n pi) face0 -
    // cos(n pi) air, as `air` is the first's: which side of each face a
    // direction lies on, without its angle.
    Vec3 air1;
    std::size_t building = 0;
    bool vertical = false;

    // The angle of `v` about the edge, from face0 through the air, in
    // [0, 2 pi): the wedge's own angle (wedge.hpp).
    double angle(const Vec3& v) const;

    // Where p stands about the edge's line: its parameter along the line from
    // `start`, and its part normal to the line.
    struct Projection
    {
        double along = 0.0;
        Vec3 normal;
    };
    Projection project(const Vec3& p) const
    {
        const Vec3 d = p - start;
        const double along = d.dot(direction);
        return {along, d - along * direction};
    }

    // Whether a direction `v` normal to the edge lies in the wedge's air, from
    // the signs of its components about the two faces (`air`, `air1`): 1 or 0
    // where they decide it, 4e-7 rad or more from a face -- wider than any
    // slack a caller lets the angle have (in_air's 1e-9, onto_face's 1e-7),
    // so that the angle would decide it the same way -- and -1 where they do
    // not, and the angle must.  On a building's edges (1 < n < 2) the air is
    // the half-turn above face0 and whatever lies on the air side of the
    // second face below it.
    int air_by_sign(const Vec3& v) const;
};

class Buildings
{
public:
    // Distances under this are contact, not passage, in the visibility test
    // and in the tests for lying on a wall (m).
    static constexpr double contact = 1e-6;

    Buildings() = default;

    // Throws std::invalid_argument for a ring with fewer than three vertices
    // or no area, a top not above the floor, sides that cross (within a ring
    // or between rings), or a hole outside its footprint or inside another.
    //
    // Buildings that share a side exactly -- the same two vertices, as
    // terraced houses in OpenStreetMap do -- stand against each other there:
    // the lower one's facade, roof edge and corners on that side are inside
    // the taller (or as tall) neighbour and are dropped, and the taller one's
    // keep only what stands above the lower roof.  A side shared only in part
    // keeps its facade and edges, which `buried` then keeps out of the rays.
    //
    // Facades whose normals agree within `orientation_tolerance` (rad) share
    // one orientation (`reflect`).
    explicit Buildings(std::vector<Footprint> footprints, double orientation_tolerance = 1e-5);

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
    // its sides (and outside its holes); the building's index, or size() for
    // none.
    std::size_t inside(double x, double y) const;

    // Whether `p` lies strictly inside a building's solid: inside a footprint
    // and between the floor and its roof, by `contact` or more.
    bool buried(const Vec3& p) const;

    // The sides found shared exactly with another building.
    std::size_t shared_sides() const
    {
        return shared_sides_;
    }

    // The first facade a ray from `pos` along the unit `dir` meets from
    // outside, and where; facades().size() for none.  A facade `pos` lies on
    // (within contact) is behind the ray: it is where the ray reflected or
    // diffracted.  Below the floor a facade's image stands from -top.  A
    // rising ray enters a building only through a wall, so for one this is
    // its first obstruction.
    std::size_t first_hit(const Vec3& pos, const Vec3& dir, bool below_floor, Vec3& hit) const;

    // `v` reflected in a facade's orientation: its horizontal component
    // mirrored.  Facades within the orientation tolerance of one another
    // reflect a direction alike.  Projection round-off leaves walls that
    // were parallel 1e-10 rad apart, and taken one by one they multiply the
    // specular traces' candidates; the default 1e-5 rad turns a ray by a
    // millimetre in 100 m.  A wider tolerance treats the walls of a street a
    // degree out of parallel as parallel, which is what lets a ray bounce up
    // it as it would up the canyon it nearly is.
    Vec3 reflect(std::size_t facade, const Vec3& v) const
    {
        const FacadeFrame& fr = frames_[orientation_of_[facade]];
        const double d = v.x * fr.nx + v.y * fr.ny;
        return {v.x - 2.0 * d * fr.nx, v.y - 2.0 * d * fr.ny, v.z};
    }

    // One facade per orientation.
    const std::vector<std::size_t>& orientations() const
    {
        return orientations_;
    }

    // The facade whose orientation a facade shares.
    std::size_t orientation(std::size_t facade) const
    {
        return orientation_of_[facade];
    }

private:
    bool blocked_by(std::size_t building, const Vec3& a, const Vec3& b, bool below_floor) const;
    bool strictly_inside(std::size_t building, double x, double y) const;

    // A facade's base point, unit base direction and length, and the unit
    // normal to its outside.
    struct FacadeFrame
    {
        double x0, y0, ax, ay, length, nx, ny;
    };

    std::vector<Footprint> footprints_;
    std::vector<Facade> facades_;
    std::vector<FacadeFrame> frames_;
    std::vector<std::size_t> orientations_;
    std::vector<std::size_t> orientation_of_; // per facade, its orientation's facade
    // The buildings' and the facades' boxes, for the queries to walk.
    PlanGrid building_grid_;
    PlanGrid facade_grid_;
    std::vector<std::size_t> facade_building_;
    std::vector<Edge> edges_;
    // Per building: x min, x max, y min, y max.
    std::vector<double> bounds_;
    double max_top_ = 0.0;
    double orientation_tolerance_ = 1e-5;
    std::size_t shared_sides_ = 0;
    double x_min_ = 0.0, x_max_ = 0.0, y_min_ = 0.0, y_max_ = 0.0;
};

} // namespace niceops
