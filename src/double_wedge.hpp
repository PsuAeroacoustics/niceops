#pragma once

// Diffraction at two parallel wedges that share a face -- over a flat roof
// from one edge to the other, around a wall from one corner to the next --
// uniform where the second edge lies in the first edge's transition region.
//
// The cascade of two UTD coefficients (wedge.hpp) carries the first edge's
// diffracted wave to the second as a ray.  Near the first edge's shadow
// boundary that wave is a transition field, not a ray, and the cascade loses
// it: against an exact solution (2-D boundary elements, tools/thick_barrier_bem.py)
// it is 1-4 dB low over a building 2.44 m thick with the source 18 degrees
// below the roof, and up to 11 dB low with the source 1.6 degrees below it.
//
// Here the field above the first edge is instead a Kirchhoff aperture: the
// half-line normal to the shared face at the first edge, every point of it a
// line source of the incident wave, each diffracted by the second wedge with
// its own UTD coefficient.  The integral holds the transition, and its end
// point at the edge is the cascade, whose first coefficient is then put right
// by the ratio of the wedge's Keller coefficient to the aperture's.  It is
// evaluated from the side whose source lies below the face's plane, and from
// both when both do, the two averaged, so that it is reciprocal.  Against the
// same solution it is within 0.7 dB for buildings 2.44 m thick or more at
// 125 Hz to 2 kHz, and still within 2.3 dB at 1 m (THEORY 6.5).  Thinner than
// a wavelength or so it fails too, less badly than the cascade.

#include <complex>
#include <vector>

#include "wedge.hpp"

namespace niceops {

// In the plane normal to the edges: the first edge at the origin, the second
// at (width, 0), the shared face between them, the air above it (z > 0).
// Each wedge's angles are measured from the shared face through the air, its
// other face at n pi.
struct SharedFace
{
    double width = 0.0;
    double n1 = 1.5;
    double n2 = 1.5;
};

// Whether shared_face_field applies: the source below the face's plane and
// before the first edge (x < 0), or the receiver below it and past the second
// (x > width).
bool shared_face_applies(const SharedFace& face, double sx, double sz, double rx, double rz);

// The field at the receiver, diffracted at both edges, of the 2-D wave
// exp(i k r) / sqrt(r) from the source.  Below the face's plane it is the
// whole field that passes the second edge; above it, the part the second edge
// diffracts, the first edge's own wave reaching the receiver being a ray of
// its own.
std::complex<double> shared_face_field(const SharedFace& face, double sx, double sz, double rx, double rz, double k);

// shared_face_field at each of the wavenumbers, to the bit, for a ray taken at
// every sub-band of a run.  Each wavenumber's aperture keeps its own nodes, but
// most of them fall on points another's already has -- they are graded from
// the same edge and cuts -- and a node's geometry (its path lengths, and the
// second wedge's angles there) is worked out once for all that share it.
std::vector<std::complex<double>> shared_face_spectrum(const SharedFace& face, double sx, double sz, double rx,
        double rz, const std::vector<double>& wavenumbers);

// shared_face_field split into what does not depend on the wavenumber, worked
// out once for a set of wavenumbers, and what does: for a band's sub-bands,
// or a GPU, that take the same face at many wavenumbers.  The nodes run to the
// aperture's extent at the smallest, graded as shared_face_field grades them
// from the edge and the cuts, with a panel end wherever any of the set's
// apertures starts its taper or ends (the taper's curvature jumps there;
// without them the Filon test's faces are 0.019 dB from the fine midpoint rule
// at worst, with them 0.006).  The end nodes sit a hair of the shortest
// aperture inside each segment (aperture_geometry has why).
//
// One node of an aperture: its place t up the aperture, and per wave (the
// second wedge's diffraction, then the incident wave and its reflections off
// the shared face and off the wedge's other face) its path length less the
// aperture's reference, that length's rate along the aperture, and its
// amplitude less the launch factor and the taper -- zero where the wave does
// not reach the receiver.  The diffracted wave's amplitude also leaves out the
// second wedge's coefficient, which is wanted at each k: its angles at the
// node's theta0, and its distance parameter.
struct ApertureNode
{
    double t = 0.0;
    double length[4] = {0.0, 0.0, 0.0, 0.0};
    double rate[4] = {0.0, 0.0, 0.0, 0.0};
    double weight[4] = {0.0, 0.0, 0.0, 0.0};
    WedgeAngles angles;
    double distance = 0.0;
};

// One end's aperture.  Its segments run between the cuts, each 2p + 1 nodes:
// p Filon panels' ends and middles, the end nodes a hair inside the segment.
struct ApertureGeometry
{
    // The path length every node's lengths are measured from.
    double reference = 0.0;
    // extent(k) = max(3 sqrt(lambda scale), 30 lambda).
    double extent_scale = 0.0;
    std::vector<ApertureNode> nodes;
    // Where each segment's nodes start, and one past the last.
    std::vector<std::size_t> segments;
    // The end point's correction, the wedge's Keller coefficient over the
    // aperture's: the source's angle below the face, and the first wedge's
    // angles on the face.  Applied where sin_g > 1e-6.
    double sin_g = 0.0, cos_g = 0.0;
    WedgeAngles keller;
};

struct SharedFaceGeometry
{
    // Each end's weight (its source's angle below the face's plane, zero where
    // it is above) and its aperture; the field is the weighted mean.
    double forward_weight = 0.0, reverse_weight = 0.0;
    ApertureGeometry forward, reverse;
};

// Throws std::invalid_argument for no wavenumbers, or one not positive.
SharedFaceGeometry shared_face_geometry(const SharedFace& face, double sx, double sz, double rx, double rz,
        const std::vector<double>& wavenumbers);
// The field shared_face_field gives at one of the geometry's wavenumbers
// (or, less closely, any other at least the smallest), times
// exp(-i k reference) of the end it comes from: the phase of the path the
// reference measures is the caller's.  Both ends share one reference, the
// forward one's, so their phases add up.
std::complex<double> shared_face_field_at(const SharedFaceGeometry& geometry, double k);
// The extent of an aperture at wavenumber k.
double aperture_extent(double extent_scale, double k);

} // namespace niceops
