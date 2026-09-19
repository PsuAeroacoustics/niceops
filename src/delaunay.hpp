#pragma once

// Delaunay triangulation in arbitrary dimension, replacing CGAL's
// Delaunay_triangulation_3 and the Boost, GMP and MPFR that came with it.
//
// Built by Bowyer-Watson incremental insertion on top of the exact predicates
// in predicates.hpp.  The dimension is a constructor argument rather than a
// template parameter, so a fourth operating parameter can be added to the
// database without touching this file.
//
// Degeneracy is the whole difficulty.  The operating conditions are a prism --
// the same flight conditions repeated at each thrust coefficient -- so large
// families of points are exactly cospherical, and the Delaunay triangulation is
// genuinely non-unique.  Any consistent choice among them is a correct answer;
// what would not be correct is answering the same in-sphere question two
// different ways, which is what floating point does and what leaves the cavity
// non-star-shaped.  A simplex whose circumsphere passes exactly through the new
// point is left standing, since a point *on* a circumsphere does not violate
// the Delaunay condition.

#include <cstddef>
#include <vector>

namespace niceops {

class Delaunay
{
public:
    // `points` is `count` rows of `dim` coordinates, row-major.
    Delaunay(std::size_t dim, std::size_t count, const double* points);

    std::size_t dimension() const
    {
        return dim_;
    }

    std::size_t point_count() const
    {
        return count_;
    }

    // Barycentric coordinates of `query` in the simplex containing it.
    // Writes dim+1 weights and the dim+1 point indices they apply to, and
    // returns true.  Returns false when the query lies outside the convex hull,
    // leaving the outputs untouched.
    bool barycentric(const double* query, double* weights, std::size_t* indices) const;

    // Index of the input point closest to `query`.  Used for queries the
    // triangulation cannot bracket.
    std::size_t nearest_vertex(const double* query) const;

    // Number of simplices in the finished triangulation, excluding those still
    // attached to the artificial enclosing vertices.  Diagnostics and tests.
    std::size_t simplex_count() const;

private:
    // A simplex holds dim+1 vertices and, opposite each, the neighbouring
    // simplex across the facet that omits it.
    struct Simplex
    {
        std::vector<std::size_t> vertices;
        std::vector<long> neighbours;
        bool alive = true;
    };

    const double* coordinates(std::size_t vertex) const
    {
        return storage_.data() + vertex * dim_;
    }

    bool is_artificial(std::size_t vertex) const
    {
        return vertex >= count_;
    }

    void build_enclosing_simplex();
    void insert(std::size_t point);
    // Simplex containing `query`, or -1 if the walk leaves the triangulation.
    long locate(const double* query) const;
    int simplex_orientation(const Simplex& simplex) const;
    // Orientation of `simplex` with vertex `replaced` swapped for `query`.
    int orientation_against(const Simplex& simplex, std::size_t replaced,
            const double* query) const;

    std::size_t dim_ = 0;
    std::size_t count_ = 0;
    // Input points followed by the dim+1 artificial vertices of the enclosing
    // simplex, so a vertex index can address either.
    std::vector<double> storage_;
    std::vector<Simplex> simplices_;
    std::vector<std::size_t> free_list_;
    // Last simplex located, as the starting point of the next walk.  Queries
    // arrive in trajectory order, so consecutive ones are usually neighbours.
    mutable std::size_t hint_ = 0;
};

} // namespace niceops
