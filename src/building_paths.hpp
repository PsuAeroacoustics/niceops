#pragma once

// The urban table of a microphone among LoD1 buildings (buildings.hpp): the
// ideal canyon's model (urban.hpp, canyon_ratio) for any arrangement of
// extruded footprints, under the same plane-wave incidence and by the same
// reciprocity.  For a plane wave arriving from the unit direction u (toward
// the source) it finds the rays that reach the microphone and sums them
// coherently, band-averaged, as a ratio to open ground's direct ray and floor
// reflection (THEORY 6.5):
//
//   - specular rays, traced back from the microphone and from its image in the
//     floor (the floor's mirror world, buildings.hpp) along every direction
//     that up to `max_turns` changes of orientation can turn into u.  Facades
//     are vertical, so a reflection keeps the vertical component and turns the
//     horizontal one; each candidate is traced deterministically, through up
//     to `max_order` reflections, and kept if it leaves along u.  The
//     candidates are the reflections of u in the facades' orientations, four
//     for a district on a grid;
//   - with `diffraction`, the rays each lit edge diffracts (wedge.hpp, a hard
//     wedge of the edge's own angle) to the microphone, to its floor image and
//     to their images in facades within `radius`, to `max_diffracted_order`
//     reflections: on Keller's cone, the diffraction point where the
//     diffracted ray leaves the edge at the incident ray's angle to it.  The
//     images are found by beam tracing in plan view: each carries the part of
//     its last facade through which the microphone can be reached, a further
//     facade counts only where it lies in that beam, and an edge only if it
//     crosses it.
//
// What the run evaluates itself is left out, as in the canyon: a facade among
// `run_facades` with a first-order reflection open from above, and that
// reflection's floor partner when it is open too.
//
// Left out: diffraction at two edges in turn, an edge lit only by a reflected
// wave or by the floor, diffraction at an edge's ends (a ray whose
// diffraction point falls off the edge is dropped), reflections off roofs
// (the microphone is on the floor), and everything canyon_ratio leaves out.

#include <complex>
#include <cstddef>
#include <vector>

#include "buildings.hpp"
#include "urban.hpp"
#include "vec3.hpp"

namespace niceops {

struct BuildingPathOptions
{
    // Facade reflections on a specular ray, at most.  Bouncing between the
    // two walls of a street costs nothing to find, so this is set by the
    // grazing directions, which climb out of a deep street slowly.
    std::size_t max_order = 64;
    // Changes of orientation on a specular ray, at most: the depth of the
    // candidate departures (a run of reflections in parallel facades turns the
    // ray once or not at all), and at most `max_departures` of them, breadth
    // first, for a district whose facades face many ways.
    std::size_t max_turns = 4;
    std::size_t max_departures = 4096;
    // Facade reflections after a diffraction, at most.
    std::size_t max_diffracted_order = 6;
    // Edges, and the facades a diffracted ray reflects in, within this
    // horizontal distance of the microphone (m).
    double radius = 150.0;
    bool diffraction = true;
    // The images the diffracted rays are tried through, at most; past it the
    // search stops (BuildingReceiver::truncated).
    std::size_t max_images = 20000;
    // Prune the images by beam tracing.  Without, every chain of facades that
    // can see each other is tried, with every edge: slower, and the reference
    // the beams are tested against (they must find the same rays).
    bool beams = true;
};

// One ray of the field at the microphone.  Its phase is k times `length` (the
// plane wave's phase at the origin being 0), and its amplitude the product of
// sqrt(reflection) over `facades` and, `floor` times, sqrt(floor reflection),
// times for a diffracted ray the wedge's coefficient over the square root of
// the distance from the edge, both in the plane normal to it:
//
//   D(n, theta, theta0, k sin(beta), s sin(beta)) / sqrt(s sin(beta)).
struct BuildingRay
{
    double length = 0.0;
    std::vector<std::size_t> facades;
    int floor = 0;
    bool diffracted = false;
    std::size_t edge = 0;
    double n = 1.5;
    double theta = 0.0;
    double theta0 = 0.0;
    double sin_beta = 1.0;
    double s = 0.0;
    // The run's first-order reflection or its floor partner.
    bool run_evaluates = false;
};

// One microphone among the buildings, with what does not depend on the wave's
// direction worked out once: the facades' frames and orientations, the edges
// within reach, and the images a diffracted ray can reach it through.
class BuildingReceiver
{
public:
    // `microphone` in the local frame, z above the floor; `run_facades` are
    // indices into buildings.facades().  The buildings must outlive this.
    BuildingReceiver(const Buildings& buildings, const Vec3& microphone, std::vector<std::size_t> run_facades,
            const BuildingPathOptions& options);

    // Every ray for the plane wave from `u` (unit, toward the source; an
    // elevation under 1e-3 rad is taken at 1e-3, as the canyon does).
    std::vector<BuildingRay> rays(const Vec3& u) const;

    // The table entry: the energy ratio to open ground per band, every ray but
    // the run's summed coherently and band-averaged over three sub-bands (as
    // canyon_ratio does), or in energy without `options.diffraction`.
    // `reflection` is [facade][band] for buildings.facades(); `out` gets one
    // value per band.
    void ratio(const Vec3& u, const std::vector<double>& reflection, double floor_reflection,
            const std::vector<double>& bands, double sound_speed, double* out) const;

    // The images of the microphone the diffracted rays are tried through, the
    // microphone itself first: each with the facades it was mirrored in, the
    // microphone's side first.
    struct Image
    {
        Vec3 point;
        std::vector<std::size_t> facades;
    };
    const std::vector<Image>& images() const
    {
        return images_;
    }

    // Whether the image search stopped at options.max_images.
    bool truncated() const
    {
        return truncated_;
    }

private:
    struct Frame
    {
        double x0, y0, ax, ay, length, nx, ny;
    };
    std::size_t first_hit(const Vec3& pos, const Vec3& dir, bool below_floor, Vec3& hit) const;
    bool diffracted_path_open(const Vec3& q, const Image& image, bool below_floor) const;

    const Buildings* buildings_;
    Vec3 microphone_;
    std::vector<std::size_t> run_facades_;
    BuildingPathOptions options_;
    std::vector<Frame> frames_;
    std::vector<std::size_t> planes_;   // one facade per orientation
    std::vector<std::size_t> edges_;    // within reach
    std::vector<Image> images_;
    std::vector<std::vector<std::size_t>> image_edges_; // per image, the edges whose rays can reach it
    bool truncated_ = false;
};

// A ray's complex amplitude at wavenumber k, `reflection` being the facades'
// energy reflection coefficients for one band, [facade].
std::complex<double> ray_amplitude(const BuildingRay& ray, const std::vector<double>& reflection,
        double floor_reflection, double k);

// A receiver's table over `axes`, [azimuth][elevation][band] in dB (floored
// at urban_gain_floor_db), its azimuths compass bearings and the local +x
// along `frame_bearing`; the directions are shared out over `threads`.
// `reflection` is [facade][band] for the buildings' facades.
std::vector<double> building_gains(const BuildingReceiver& receiver, const UrbanTable::Axes& axes,
        double frame_bearing, const std::vector<double>& reflection, double floor_reflection, double sound_speed,
        unsigned threads = 1);

// The facades the run evaluates for a microphone: those within `radius`
// (horizontal, m) whose outer face looks at it.
std::vector<std::size_t> facing_facades(const Buildings& buildings, const Vec3& microphone, double radius);

} // namespace niceops
