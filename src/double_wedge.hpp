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

} // namespace niceops
