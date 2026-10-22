#pragma once

// A row of flat-roofed blocks in cross-section, crossed near grazing: the
// field behind it by a parabolic equation over the roofs and the streets
// between them, finished through the last block's back corner (THEORY 6.5,
// "Over a row of roofs").  Where three or more roof edges lie in each other's
// transitions -- a street between two rows at grazing -- the UTD cascades of
// two, three and four edges do not converge (+8 dB at 500 Hz over one street,
// up to +14 over two, against boundary elements); the march carries the whole
// field over the row instead.
//
// Two-dimensional, in the plane normal to the blocks' edges: u across them, z
// up, a line source.  The march is split-step Pade (Collins; the propagators
// of pade_propagator.hpp) on the wide-angle operator,
//
//   psi(u + du) = r(q) psi(u) ~ exp(i k du (sqrt(1 + q) - 1)) psi(u),
//   q = d^2 / (k^2 dz^2) + i a(z),  p = psi exp(i k u),
//
// in steps of up to eight wavelengths (the rest of each roof or street under
// an eighth of one by a Crank-Nicolson step), with rigid ground (a ghost node)
// at the roof under it or at the street's floor and an absorbing layer at the
// top, both sized in profile_pe.cpp.  It starts on the first block's front
// face with the source's wave above the roof (Kirchhoff: what would meet the
// face below is dropped, as it is wherever the ground steps up), and stops a
// few wavelengths before the last block's back edge.  There the field and its
// u-derivative on that vertical line, a Huygens line, reach the receiver
// through the back corner's own Green's function -- the rigid wedge's, by
// UTD: its geometrical-optics waves where they see the receiver, and its
// diffraction:
//
//   p(R) = int ( -G_W(Q, R) dp/du + p dG_W/du ) dz.
//
// The source's floor image starts a second march on the same matrices, and
// the receiver's takes the integral again: four parts, which the point model
// carries as four rays (PointPathOptions::profile_pe).  Taking the line's
// field into the shadow with the free-space (Neumann) kernel instead missed
// the back corner and was 1.3 dB low behind one block.  What is left out: the
// waves the march cannot carry back, the street's reverberation between its
// walls, and the roofs' and the streets' floors' absorption, rigid here.
// Against boundary elements (FIDELITY): one street 1.6 dB rms (UTD 3.4), two
// streets 2.1 (6.7), 1 kHz 0.3 (3.5).

#include <array>
#include <complex>
#include <vector>

namespace niceops {

struct ProfileBlock
{
    double u0 = 0.0, u1 = 0.0; // its front and back faces
    double h = 0.0;            // its roof's height
};

struct ProfileProblem
{
    // Ascending in u, apart (a street between each pair) or touching.
    std::vector<ProfileBlock> blocks;
    // The source before the first block, the receiver behind the last.
    double su = 0.0, sz = 0.0;
    double ru = 0.0, rz = 0.0;
    // The floor's amplitude reflection at the two ends, profile_field's
    // weight on its images; the streets' floors and the roofs are rigid.
    double floor = 1.0;
};

// Whether the problem is one profile_field takes: blocks of positive width
// and height, ascending, the source before the first and the receiver behind
// the last, neither below the floor.
bool profile_applies(const ProfileProblem& problem);

// The field at the receiver in four parts, part[a][b]: from the source
// (a = 0) or its image in the floor (a = 1), to the receiver (b = 0) or its
// image (b = 1), the floor rigid.  Each goes as exp(i k l) for its own path
// length l (profile_path), so that a sum over a band can tell them apart as
// it does the traced rays off the floor.  A line source whose free field is
// exp(i k r) / sqrt(r) (shared_face_field's normalisation), at wavenumber k.
struct ProfileParts
{
    std::complex<double> part[2][2];
};
ProfileParts profile_parts(const ProfileProblem& problem, double k);

// The parts summed, the images weighed by problem.floor.
std::complex<double> profile_field(const ProfileProblem& problem, double k);

// profile_parts at each wavenumber.
std::vector<ProfileParts> profile_spectrum(const ProfileProblem& problem, const std::vector<double>& wavenumbers);

// A part's path in the plane: the taut string over the blocks' corners from
// the source or its floor image (source_image) to the receiver or its image,
// as (u, z) points from end to end, the images below the floor; and its
// length.
double profile_path(const ProfileProblem& problem, bool source_image, bool receiver_image,
        std::vector<std::array<double, 2>>* points = nullptr);

} // namespace niceops
