#pragma once

// Exact sign evaluation for the determinants a Delaunay triangulation is built
// out of, in arbitrary dimension.
//
// Floating point is not good enough here, and not marginally so.  The database's
// operating conditions are a prism -- 164 distinct (mu, gamma) flight conditions
// repeated at each of 8 thrust coefficients -- so whole families of points are
// *exactly* cospherical.  An in-sphere test on them has a true determinant of
// zero, and a rounded one of whatever sign the round-off happens to land on.
// Answer such a test inconsistently and Bowyer-Watson's cavity stops being
// star-shaped: the triangulation comes out with holes, overlaps, or an infinite
// loop, rather than merely a slightly different tessellation.
//
// So the determinant sign is computed exactly, using Shewchuk's floating-point
// expansions ("Adaptive Precision Floating-Point Arithmetic and Fast Robust
// Geometric Predicates", Discrete & Computational Geometry 18:305-363, 1997).
// An expansion represents a real number exactly as an unevaluated sum of
// non-overlapping doubles; the primitives below add and multiply them with no
// rounding at all, so the final sign is the true sign.
//
// A fast path runs first: evaluate in plain double with a running error bound,
// and only fall back to the exact path when the bound admits zero.  For points
// in general position the fast path settles it.
//
// NOTE: these routines depend on strict IEEE-754 semantics.  Do not build this
// file with -ffast-math or -funsafe-math-optimizations; the two_sum trick relies
// on the compiler not "simplifying" (x - a) + (a - x) to zero.

#include <cstddef>
#include <vector>

namespace niceops {

// Sign of the determinant of a `size` x `size` matrix, given row-major.
// Returns -1, 0 or +1, and the zero is exact rather than a tolerance.
int exact_determinant_sign(const double* matrix, std::size_t size);

// Sign of the volume of the simplex on `dim + 1` points, each `dim` long:
// positive when the points are positively oriented, zero when degenerate.
// Points are given as a contiguous (dim + 1) x dim array.
int orientation(const double* points, std::size_t dim);

// Whether `query` lies inside the circumsphere of the simplex on `dim + 1`
// points: +1 inside, -1 outside, 0 exactly on it.  The sign is normalised so
// that it does not depend on the orientation of the simplex.
//
// This is the lifted orientation test: points rise to the paraboloid
// x -> (x, |x|^2) in one more dimension, where "inside the circumsphere"
// becomes "below the hyperplane through the lifted simplex".
int in_sphere(const double* points, const double* query, std::size_t dim);

} // namespace niceops
