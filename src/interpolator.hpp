#pragma once

// Barycentric interpolation over a database's operating conditions.
//
// The triangulation is Delaunay, built by delaunay.hpp on the exact predicates
// in predicates.hpp.  This used to be CGAL's Delaunay_triangulation_3, which
// brought Boost, GMP and MPFR with it for one 3D triangulation of about a
// thousand points.
//
// It is an object, owned by the Database whose conditions it indexes.  It was
// briefly a set of file-scope globals -- the shape of the D original's "C" FFI
// seam -- which meant only one Database could be live at a time: constructing a
// second silently re-armed the interpolator under the first.  That cost a
// reload dance in main, a comment explaining the ordering constraint, and a
// stale-cache bug when a batch run switched databases.
//
// The dimension is a constructor argument rather than a fixed three, so adding
// a fourth operating parameter is a change to the caller and not to the
// geometry.

#include <cstddef>
#include <memory>
#include <vector>

namespace niceops {

class Delaunay;

// Dimension of the operating-condition space the database is indexed by
// (advance ratio, flight path angle, thrust coefficient).  The triangulation
// and predicates underneath are dimension-generic; this is where the caller
// fixes it, and what Database::Blend is sized from.  Adding a fourth parameter
// means changing this and the packing in Database::Database, not the geometry.
inline constexpr std::size_t condition_dims = 3;

class ConditionInterpolator
{
public:
    ConditionInterpolator();
    ~ConditionInterpolator();
    ConditionInterpolator(ConditionInterpolator&&) noexcept;
    ConditionInterpolator& operator=(ConditionInterpolator&&) noexcept;

    // Build over `count` points of `dim` coordinates, row-major.
    //
    // `extend_domain_by` widens the triangulation past the unit cube the
    // coordinates are normalised into: any point already sitting on a face of
    // that cube is duplicated out at the widened bound, carrying the same data
    // index, so that queries at the very edge of the operating envelope still
    // land inside a finite simplex instead of falling off the hull.
    void build(std::size_t dim, std::size_t count, const double* points,
            double extend_domain_by = 0.1);

    // Barycentric weights at `query` (`dim` coordinates), and the indices of the
    // input points they apply to.  Both outputs must have room for dim+1
    // entries.  The query is clamped into the unit cube first.
    void barycentric(const double* query, double* weights, std::size_t* indices) const;

private:
    std::size_t nearest(const double* query) const;

    std::unique_ptr<Delaunay> mesh_;
    // Which input point each triangulation vertex stands for.  Not the identity:
    // the domain extension above adds duplicates carrying their original's index.
    std::vector<std::size_t> vertex_to_point_;
    // Vertex coordinates, kept so the nearest-condition fallback works even when
    // there are too few points to triangulate at all.
    std::vector<double> vertex_coordinates_;
    std::size_t dimension_ = 0;
};

} // namespace niceops
