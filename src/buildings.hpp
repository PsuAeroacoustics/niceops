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

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <tuple>
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

// A plane of a pitched roof, z = a x + b y + c (m).
struct RoofPlane
{
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double at(double x, double y) const
    {
        return a * x + b * y + c;
    }
};

// A footprint: its outer ring's vertices in the run's local frame (m), in
// either orientation, the first not repeated at the end; its roof's height
// above the floor; and the holes it surrounds (courtyards, open to the sky).
// Buildings orients every ring so that the solid is on its left: the outer
// one counterclockwise, the holes clockwise.  `roof`, when given, makes the
// roof pitched (LoD2): its height over (x, y) is the least of its planes there
// -- a gable is two, a hip four, a mansard more -- so the roof is convex,
// every ridge between two planes an edge.  A roof that is not (a valley, a
// step) is several buildings sharing sides.  Buildings then sets `top` to the
// highest point.
struct Footprint
{
    std::vector<double> x;
    std::vector<double> y;
    double top = 0.0;
    std::vector<Ring> holes;
    std::vector<RoofPlane> roof;
};

// The roof's height over (x, y): `top`, or the least of the planes.
inline double roof_height(const Footprint& f, double x, double y)
{
    if (f.roof.empty())
        return f.top;
    double h = f.roof[0].at(x, y);
    for (std::size_t i = 1; i < f.roof.size(); i++)
        h = std::min(h, f.roof[i].at(x, y));
    return h;
}

// A pitched roof's highest point over its footprint (`roof` not empty).
double roof_top(const Footprint& f);

// The stretches of the segment from (ax, ay) to (bx, by) over which one plane
// of a pitched roof (`roof` not empty) is the least, as (t0, t1, plane), t
// from 0 at a to 1 at b, neighbours of one plane merged: the roof's profile
// along it is straight over each.
std::vector<std::tuple<double, double, std::size_t>> roof_stretches(const Footprint& f, double ax, double ay,
        double bx, double by);

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
    // -- a T-junction, or walls staggered along one line -- is first split at
    // the other building's vertices on it (footprints() then has them), so
    // that its shared part is shared exactly; a building touching a side at
    // one vertex leaves it whole.  A vertex within `contact` of the line
    // between its neighbours is not a corner.
    //
    // Facades whose normals agree within `orientation_tolerance` (rad) share
    // one orientation (`reflect`).
    explicit Buildings(std::vector<Footprint> footprints, double orientation_tolerance = 1e-5);

    // Which district this is, for a cache that outlives one query: each
    // construction takes a new one, and a copy keeps its original's.
    std::uint64_t id() const
    {
        return id_;
    }

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

    // The facade whose coefficients a building's roof takes -- its walls and
    // roof absorb and scatter alike -- its first; no_facade for a building
    // with none left standing.
    static constexpr std::size_t no_facade = static_cast<std::size_t>(-1);
    std::size_t roof_facade(std::size_t building) const
    {
        return roof_facade_[building];
    }

    const std::vector<Edge>& edges() const
    {
        return edges_;
    }

    // A building's box in plan: x min, x max, y min, y max.
    const double* bounds(std::size_t building) const
    {
        return &bounds_[4 * building];
    }
    // A building's edges: [edge_begin(b), edge_begin(b + 1)) of edges().
    std::size_t edge_begin(std::size_t building) const
    {
        return edge_begin_[building];
    }

    // The tallest roof.
    double max_top() const
    {
        return max_top_;
    }

    // Whether the segment from a to b passes through a building: some stretch
    // longer than `contact` lies strictly inside one.
    bool blocked(const Vec3& a, const Vec3& b, bool below_floor = false) const;
    // The same for one building.
    bool blocked_by(std::size_t building, const Vec3& a, const Vec3& b, bool below_floor) const;
    // The buildings whose boxes may meet the box [x_min, x_max] by
    // [y_min, y_max], each once, in no particular order: the grid's
    // candidates, which the caller tests.
    void near_box(double x_min, double x_max, double y_min, double y_max, std::vector<std::uint32_t>& out) const;
    // The plan segment from a to b across the buildings: each stretch of it
    // inside a footprint (outside its holes), from d0 to d1 metres along it
    // from a, and the building; in no particular order.
    struct SectionSpan
    {
        double d0 = 0.0;
        double d1 = 0.0;
        std::uint32_t building = 0;
    };
    void section(const Vec3& a, const Vec3& b, std::vector<SectionSpan>& out) const;
    // The same for the facades, ascending.
    void facades_in_box(double x_min, double x_max, double y_min, double y_max, std::vector<std::uint32_t>& out) const;

    // Whether a ray from `a` along `direction` reaches infinity without
    // passing through a building.  A ray that does not rise escapes once past
    // every footprint.
    bool escapes(const Vec3& a, const Vec3& direction, bool below_floor = false) const;

    // Whether (x, y) lies strictly inside a footprint, `contact` or more from
    // its sides (and outside its holes); the building's index, or size() for
    // none.
    std::size_t inside(double x, double y) const;

    // Whether (x, y) lies on a building's level roof: within its footprint, its
    // sides included to round-off, and under a pitched roof where its least
    // plane is level (a mansard's deck).  A roof's reflection ends where its
    // edge's reflection boundary is, and must not stop short of it as inside()
    // does.
    bool on_roof(std::size_t building, double x, double y) const;

    // Whether `p` lies strictly inside a building's solid: inside a footprint
    // and between the floor and its roof, by `contact` or more.
    bool buried(const Vec3& p) const;

    // What a point lies on, within `contact`, besides the floor (whose
    // reflection every ray's floor partner carries): a wall (its facade), a
    // level roof (the building: a flat roof, or a pitched one's level plane
    // where that is least), a pitched roof's slope, or an edge -- a building's
    // edge, or two surfaces meeting in a concave corner -- where no one
    // surface's reflection stands for it.
    struct Support
    {
        enum class Kind
        {
            none,
            facade,
            roof,
            slope,
            edge
        };
        Kind kind = Kind::none;
        std::size_t index = 0;
    };
    Support support(const Vec3& p) const;

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
    std::vector<std::size_t> roof_facade_;
    std::vector<Edge> edges_;
    std::vector<std::size_t> edge_begin_;
    // Per building: x min, x max, y min, y max.
    std::vector<double> bounds_;
    double max_top_ = 0.0;
    double orientation_tolerance_ = 1e-5;
    std::uint64_t id_ = 0;
    std::size_t shared_sides_ = 0;
    // Per building, the sides (x0, y0, x1, y1) it shares with an earlier
    // building whose roof is level with its own, which on_roof leaves to that
    // one.
    std::vector<std::vector<std::array<double, 4>>> ceded_;
    double x_min_ = 0.0, x_max_ = 0.0, y_min_ = 0.0, y_max_ = 0.0;
};

// The face of the edge's wedge that p lies on, within Buildings::contact of
// its plane and on its side of the edge: 0 for face0, 1 for the other, -1
// for none.
int edge_face(const Edge& e, const Vec3& p);

} // namespace niceops
