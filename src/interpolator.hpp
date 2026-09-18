#pragma once

// Barycentric interpolation over the database's operating conditions.
//
// The triangulation is Delaunay, built by delaunay.hpp on the exact predicates
// in predicates.hpp.  This used to be CGAL's Delaunay_triangulation_3, which
// brought Boost, GMP and MPFR with it for one 3D triangulation of about a
// thousand points.
//
// The interface takes the dimension as an argument rather than fixing it at
// three, so adding a fourth operating parameter to the database is a change to
// its caller and not to the geometry.

#include <cstddef>

namespace niceops {

// Build the interpolator over `count` points of `dim` coordinates, row-major.
// Returns non-zero on failure.
//
// `extend_domain_by` widens the triangulation past the unit cube the
// coordinates are normalised into: any point already sitting on a face of that
// cube is duplicated out at the widened bound, carrying the same data index, so
// that queries at the very edge of the operating envelope still land inside a
// finite simplex instead of falling off the hull.
int set_interpolator(std::size_t dim, std::size_t count, const double* points,
        double extend_domain_by = 0.1);

// Barycentric weights at `query` (`dim` coordinates), and the indices of the
// input points they apply to.  Both outputs must have room for dim+1 entries.
// The query is clamped into the unit cube first.
void barycentric_weights(const double* query, double* weights, std::size_t* indices);

} // namespace niceops
