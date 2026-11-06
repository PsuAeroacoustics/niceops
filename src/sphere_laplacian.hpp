#pragma once

// The Laplace-Beltrami operator on a database's (phi, theta) sphere grid, for the
// first-order near-field correction of a tone (THEORY 4.4): the field at range R
// is A + (i / (2 k R)) Lap[D A] for a far-field directivity A.

#include <complex>
#include <vector>

namespace niceops {

// The unit vector of a grid direction (forward, starboard, up), phi and theta in
// degrees as the database's grid holds them: theta from the nose, phi about it,
// 0 straight down and +90 to starboard.
struct GridDirection
{
    double forward, starboard, up;
};
GridDirection grid_direction(double phi_deg, double theta_deg);

// Lap f at every node of the grid, f indexed [i_phi * thetas.size() + i_theta]
// over the sorted grid (phis and thetas in degrees). A full circle in phi (its
// steps close the seam) is periodic; a hemisphere from phi = -90 to +90 (the
// lower half) is extended evenly across the horizontal plane, as Database reads
// its other half; any other range's first and last rows are an open edge.
// theta = 0 and 180, where present, are the grid's poles.
//
// `order` n (|n| >= 2) differentiates, at each node where it is the smoother of
// the two (the smaller differences to its grid neighbors, relative to the largest
// of the five values), the
// unwound function B = f exp(-i n alpha), alpha the direction's azimuth about the
// up axis toward starboard (the shaft, as ToneSource::unwind_order), and restores
// the winding exactly:
//
//     Lap(e^{i n alpha} B) = e^{i n alpha} [Lap B + 2 i n grad(alpha).grad(B) - n^2 |grad(alpha)|^2 B],
//
// alpha being harmonic away from the axis, so a coarse grid differentiates a
// smooth function; elsewhere f itself.  Order +-1 is never unwound (its B, like
// sin(beta), has a cone point on the shaft), nor a node on the shaft (alpha is
// undefined there), nor a tone whose content near the shaft does not wind (B
// would wind there instead), which the per-node choice catches.  Interior nodes
// take three-point differences in theta and phi (second order); the grid's poles
// and an open edge (in phi, or in theta without a pole) take a quadratic
// least-squares fit over the nodes within two grid steps, in the tangent plane
// (second order where they surround the node, first order on an edge).
std::vector<std::complex<double>> sphere_laplacian(const std::vector<double>& phis_deg,
        const std::vector<double>& thetas_deg, const std::vector<std::complex<double>>& f, int order);

} // namespace niceops
