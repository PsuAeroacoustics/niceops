#pragma once

// Diffraction at a rigid wedge's finite edge by Svensson's edge-source
// integral: the Biot-Tolstoy-Medwin solution as a line integral of secondary
// sources along the edge (Svensson, Fred and Vanderkooy 1999), in the
// frequency domain (Svensson, Calamia and Nakanishi 2009).  For a unit point
// source at S, whose free field is exp(i k r) / r (the time convention
// exp(-i w t), as wedge.hpp's), the diffracted field at R is
//
//   p_d = -(nu / 4 pi) sum_i int beta_i(z) exp(i k (m + l)) / (m l) dz,
//
// m and l the distances from S and R to the edge's point z, nu = 1 / n the
// wedge index (the open angle n pi), and the secondary sources' directivities
//
//   beta_i = sin(nu phi_i) / (cosh(nu eta) - cos(nu phi_i)),
//   phi_i = pi +- theta_S +- theta_R,
//   cosh(eta) = (m l + (z - z_S)(z - z_R)) / (r_S r_R),
//
// theta and r each end's angle from face 0 and distance from the edge's
// line.  Over the whole line this is the exact field of the infinite wedge;
// over a finite edge it holds the edge's ends, which UTD leaves out.  Unlike
// UTD it needs no high-frequency limit: it holds for edges a wavelength
// apart, or nearer.
//
// It is integrated in eta, which runs from 0 at the apex (the stationary
// point, the edge's point on the least path) out each side to the stretch's
// ends; the two sides share eta's integrand, and dz/(m l) = d eta / R with
// R(eta) = m + l = sqrt(r_S^2 + r_R^2 + (z_R - z_S)^2 + 2 r_S r_R cosh eta).
// Near a shadow or reflection boundary a beta_i is a narrow peak at eta = 0,
// and its integral is closed, 2/nu atan(tanh(nu eta / 2) / tan(nu phi / 2)):
// that part is taken analytically against the integrand's value at the
// branch's start, and the smooth rest by Gauss-Legendre, out to a dozen
// radians of phase.  On the boundary the closed part is -+ half the
// geometrical-optics wave, which is what makes the total continuous there.
// Past that the integrand oscillates, and is taken along its steepest-descent
// paths from the window's end and from the stretch's, where it decays
// instead (Huybrechs and Vandewalle 2006): a fixed number of points however
// high the frequency or long the edge.

#include <complex>
#include <vector>

#include "buildings.hpp"
#include "vec3.hpp"

namespace niceops {

// How the source S and receiver R stand about an edge, in its own frame:
// distances from its line, angles from face 0 (Edge::angle), and parameters
// along it from its start.
struct EdgeFrameEnds
{
    double r_s = 0.0, theta_s = 0.0, z_s = 0.0;
    double r_r = 0.0, theta_r = 0.0, z_r = 0.0;
};
EdgeFrameEnds edge_frame_ends(const Edge& edge, const Vec3& source, const Vec3& receiver);

// The diffracted field above, over the stretch [z_lo, z_hi] of the edge (from
// its start), for a wedge of open angle n pi with rigid faces.
std::complex<double> edge_source_field(double n, const EdgeFrameEnds& ends, double k, double z_lo, double z_hi);

// edge_source_field at each of a set of wavenumbers, which share the
// integrand's nodes: Gauss-Legendre in eta near the branches' starts, fine
// enough for the largest wavenumber, and past that Filon's rule in R, which
// is exact in the oscillation at any wavenumber.  The cost is a fixed number
// of points plus a few dozen operations per wavenumber.
std::vector<std::complex<double>> edge_source_spectrum(double n, const EdgeFrameEnds& ends, double z_lo, double z_hi,
        const std::vector<double>& wavenumbers);

// The same from the points themselves, over the whole edge.
std::complex<double> edge_source_field(const Edge& edge, const Vec3& source, const Vec3& receiver, double k);

// A smooth window on path length: 1 to reference + reach / 2, Planck's taper
// from there to 0 at reference + reach.
struct PathWindow
{
    double reference = 0.0, reach = 0.0;
};
double window_taper(double path, const PathWindow& window);

// edge_source_field with each path m + l weighted by the window, by
// Gauss-Kronrod on the real line to about `tolerance`.
std::complex<double> edge_source_field(double n, const EdgeFrameEnds& ends, double k, double z_lo, double z_hi,
        const PathWindow& window, double tolerance = 1e-7);

// A stretch of an edge, [lo, hi] in its parameter from its start.
struct Stretch
{
    double lo = 0.0, hi = 0.0;
};

// Diffraction at two edges in turn, by the edge-source integral's second
// order (Svensson, Fred and Vanderkooy 1999):
//
//   p_dd = f (nu1 nu2 / 16 pi^2) sum_i sum_j int int beta1_i beta2_j exp(i k (m + n + l)) / (m n l) dz1 dz2,
//
// m from the source to the first edge's point z1, n from there to the second
// edge's point z2, l from there to the receiver; beta1 the first edge's
// directivity for the source and the point z2, beta2 the second's for the
// point z1 and the receiver.  f is `plane_factor`: 1, or 1/2 where the path
// between the edges runs along a face that both share -- where each edge's
// directivity counts that face's reflection, and the pair would count it
// twice.
//
// This is the reference: nested Gauss-Kronrod quadrature on the real line to
// relative accuracy about `tolerance`, the inner integral in the second edge's
// eta with its peak closed as edge_source_field closes it.  Its cost grows as
// k^2 times the stretches' lengths.  With `window_phase` it takes only the
// paths within that many radians of the least, the outer half of that
// tapered smoothly: the far ends' waves are lost, but long edges are cheap.
std::complex<double> twice_diffracted_field(const Edge& first, const Edge& second, const Vec3& source,
        const Vec3& receiver, double k, const Stretch& along_first, const Stretch& along_second, double plane_factor,
        double tolerance = 1e-6, double window_phase = 0.0);

// twice_diffracted_field without a window at each of a set of wavenumbers,
// by its level curves: the integrand summed over each curve of equal path
// m + n + l, which does not depend on the wavenumber, then transformed by
// Filon's rule in the path.  Exact over the stretches, ends and corners
// included.
std::vector<std::complex<double>> twice_diffracted_spectrum(const Edge& first, const Edge& second, const Vec3& source,
        const Vec3& receiver, const Stretch& along_first, const Stretch& along_second, double plane_factor,
        const std::vector<double>& wavenumbers);

// The least of m + n + l over the first edge's stretch and the second edge's
// line: the window's reference in twice_diffracted_field.
double least_twice_path(const Edge& first, const Edge& second, const Vec3& source, const Vec3& receiver,
        const Stretch& along_first);

// The once-diffracted ray of UTD from source to receiver by the edge, as the
// ray model takes it (point_paths.cpp): its Keller point's, where that lies
// in [lo, hi], weighted by the window at its length, and 0 otherwise.
std::complex<double> keller_once_field(const Edge& edge, const Vec3& source, const Vec3& receiver, double k,
        double lo, double hi, const PathWindow& window = {});

// Where a geometrical-optics wave of one edge of a twice-diffracted pair
// reaches part of the other edge, not all of it: the first edge's incident
// wave from the source and its faces' reflections (from the source's images
// in them), at the second edge's points; the second edge's waves toward the
// receiver (its images in the second edge's faces) from the first edge's.
// Each one's boundary crosses the other edge, and there the second-order
// integral carries half that wave's once-diffracted field, which in the
// exact sum is cancelled by the once-diffracted field cut off where the wave
// stops.
struct PartnerWave
{
    bool over_second = true;  // a wave of the first edge, diffracted by the second
    Vec3 image;               // the source's image (over_second), or the receiver's
    Stretch part;             // the stretch of the other edge the wave reaches
};
std::vector<PartnerWave> twice_diffracted_partners(const Edge& first, const Edge& second, const Vec3& source,
        const Vec3& receiver);

// The twice-diffracted field as a ray model wants it: the second-order
// integral, windowed about its least path to `window_phase` radians, plus for
// each partner wave its once-diffracted integral over the part it reaches
// (in the same window) less the ray keller_once_field gives there.  The
// boundary terms cancel, and what is left belongs to this pair alone: a
// partner ray the model traces adds back the part subtracted, and one it does
// not trace was never in the field.  Over the edges' whole lengths.
std::complex<double> localized_twice_diffracted_field(const Edge& first, const Edge& second, const Vec3& source,
        const Vec3& receiver, double k, double plane_factor, double window_phase = 40.0, double tolerance = 1e-5);

// The ray model a partner's once-diffracted ray takes, which the localized
// field subtracts where that ray exists: UTD's (keller_once_field), or the
// edge-source integral over the whole edge.
enum class PartnerModel
{
    keller,
    edge_source
};

// localized_twice_diffracted_field without a window, at each of a set of
// wavenumbers, over the edges' whole lengths: the level-curve spectrum and the
// partners' once-diffracted spectra over the parts their waves reach, less
// the partner's own ray in the given model where its Keller point lies on
// the part.
std::vector<std::complex<double>> localized_twice_diffracted_spectrum(const Edge& first, const Edge& second,
        const Vec3& source, const Vec3& receiver, double plane_factor, const std::vector<double>& wavenumbers,
        PartnerModel model);

} // namespace niceops
