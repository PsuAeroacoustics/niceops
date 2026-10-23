#pragma once

// Rays a facade's relief sends into its grating orders (facade_relief.hpp,
// THEORY 6.5), from a point source to a microphone among LoD1 buildings: the
// non-specular orders, which the specular rays (point_paths.hpp) leave out --
// they carry the order-0 share of each relief they reflect in.
//
// A ray of order m off facade F reflects there with its direction's
// component along F's base direction raised by m lambda / period, and its
// component up F kept: it leaves on the cone of the incident elevation.
// Before F and after it, it reflects specularly in the facades the two ends'
// image trees reach F through: the source's images and the microphone's that
// end in F, each pair with at most `max_reflections` other facades between
// them, and the floor at either end or not.  The point on F's plane is where
// |P - S'| + |P - M'| + m (lambda / period) x is least, S' and M' the ends'
// images through the other facades, x along F: the generalized law of
// reflection, solved by Newton's method from the specular point.  The ray
// must meet F within its rectangle and every leg must be open.
//
// What it carries is energy, added to the coherent rays' band energy: the
// order's share at the incidence and in-plane frequency it meets F at, each
// other reflection's energy coefficient (a relief's order-0 share among
// them), the floor's, the air's over its length, and the spreading of its ray
// tube, worked out by perturbing its launch: |p|^2 = dOmega / dA, which is
// 1 / L^2 for order 0.  One ray per band, at the band's centre, since the
// geometry follows the wavelength, its share the order's mean over the band's
// three sub-bands, as the specular rays take order 0's in each.  The floor's
// coefficient goes with its phase, and at a plate the ray takes P_d, or P_r
// off the floor (P_d + R P_r on a roof).  Paths with two non-specular orders
// are left out.

#include <complex>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "buildings.hpp"
#include "vec3.hpp"

namespace niceops {

class FacadeRelief;
class PathEnd;

struct GratingRay
{
    // The band it carries energy in, and that energy, |p|^2 for a unit source
    // (1/r^2 in free field), before the source's own level along `launch`;
    // and the phase the floor's and a plate's coefficients give it, radians,
    // besides its length's.
    std::size_t band = 0;
    double energy = 0.0;
    double phase = 0.0;
    // The direction it leaves the source in (unit; toward the floor's image
    // where its first leg reflects in the floor), its length, the facade it
    // takes its order off, and the order.
    Vec3 launch;
    double length = 0.0;
    std::size_t facade = 0;
    int order = 0;
    // Whether its leg from the source, or to the microphone, reflects in the
    // floor.
    bool source_floor = false;
    bool microphone_floor = false;
    // Rays of one group are one order off the same facades in one band, with
    // the floor or without: near the floor they are one length, and they add
    // with their phases (grating_band_energies).
    std::size_t group = 0;
    // Where it goes: the source, its reflection points, the point on the
    // relief, the microphone (a leg off the floor drawn unfolded below it).
    std::vector<Vec3> path;
};

struct GratingInput
{
    // Per facade of the buildings, its relief, or null for a plain one.
    const std::vector<const FacadeRelief*>* reliefs = nullptr;
    // Band centres (Hz), the sound speed, and each facade's energy reflection
    // coefficient [facade][band] (a relief's order shares then divide it).
    std::vector<double> bands;
    double sound_speed = 343.2;
    const std::vector<double>* reflection = nullptr;
    // The floor's reflection coefficient for a bounce in band b, at the
    // bounce's grazing sine and the length of its leg (complex: its phase goes
    // with the ray's, GratingRay::phase); null for a rigid floor.
    std::function<std::complex<double>(double sin_grazing, double leg, std::size_t band)> floor;
    // A plate the microphone is flush on: P_d and P_r in band b for a ray
    // arriving at `elevation` above it from horizontal `azimuth` (FloorModel::
    // plate), and on a wall's or a roof's surface that surface's amplitude R
    // in band b (empty on the floor).  Null for none.
    std::function<std::pair<std::complex<double>, std::complex<double>>(
            double elevation, double azimuth, std::size_t band)>
            plate;
    std::vector<double> plate_surface;
    // The air's attenuation of amplitude, Np/m per band; empty for none.
    std::vector<double> air;
    // Facade reflections on a grating ray besides its order, at most, and of
    // those between its relief and the microphone: none in the hybrid, whose
    // diffuse field carries what a relief scatters onto other walls.  And the
    // least share of an order taken.
    std::size_t max_reflections = 2;
    std::size_t max_microphone_reflections = static_cast<std::size_t>(-1);
    double least_share = 1e-4;
    // Whether the source's side reflects in the floor (PointPathOptions::
    // source_reflections).
    bool source_reflections = true;
};

// The grating rays between the two ends (traced as point_paths.hpp traces
// them, with the same options), appended to `out`.
void grating_rays(const Buildings& buildings, const PathEnd& source, const PathEnd& microphone,
        const GratingInput& input, std::vector<GratingRay>& out);

// The grating rays' energy per band, |p|^2 for a unit source, into `out`
// (one per band, added to): each group's rays summed with their phases over
// the band (band_averaged_energy, their amplitudes the roots of their
// energies), the groups in energy.  `power`, when given, multiplies each
// ray's energy first (a directional source along its launch), one per ray.
void grating_band_energies(const std::vector<GratingRay>& rays, const std::vector<double>& bands, double sound_speed,
        const double* power, double* out);

// The point on a facade's plane where |P - a| + |P - b| + kick x is least, x
// along the facade's base direction from its start, z up: (x, z), and false
// where there is none (|kick| too large for the two to turn the ray).  `a`
// and `b` on the facade's front.
bool grating_point(const Facade& facade, const Vec3& a, const Vec3& b, double kick, double& x, double& z);

// The energy |p|^2 a ray tube carries from a unit source at `a` by the point
// `p` on the facade's plane to `b`, its along component raised there by
// `kick`: dOmega / dA at b, by perturbing its launch.  1 / |b - a'|^2 for
// kick 0, a' the mirror image.
double grating_spreading(const Facade& facade, const Vec3& a, const Vec3& p, const Vec3& b, double kick);

} // namespace niceops
