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
#include <functional>
#include <vector>

#include "buildings.hpp"
#include "diffuse_facades.hpp"
#include "urban.hpp"
#include "vec3.hpp"
#include "wedge.hpp"

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
    // A diffracted ray's departure and arrival in the edge's own angle, as
    // wedge_angles takes them: the coefficient is wanted at every sub-band's
    // wavenumber.
    WedgeAngles angles;
    double sin_beta = 1.0;
    double s = 0.0;
    // The run's first-order reflection or its floor partner.
    bool run_evaluates = false;
};

// The departures a specular ray may leave a point along, for the plane wave
// from `u`: u's reflections in `planes` (facades, one per orientation), to
// `max_turns` changes of orientation and at most `max_departures` of them,
// breadth first.
std::vector<Vec3> specular_departures(const Buildings& buildings, const std::vector<std::size_t>& planes,
        const Vec3& u, std::size_t max_turns, std::size_t max_departures);

// One facade per orientation among the facades within `radius` (horizontal,
// m) of `point`: the planes its specular departures are made in.
std::vector<std::size_t> orientations_near(const Buildings& buildings, const Vec3& point, double radius);

// The edges within `radius` (horizontal, m) of `point`, or of any of
// `points`, ascending: the edges a ray from or to there may diffract at.
std::vector<std::size_t> edges_near(const Buildings& buildings, const Vec3& point, double radius);
std::vector<std::size_t> edges_near(const Buildings& buildings, const std::vector<Vec3>& points, double radius);

// A specular ray traced back from a point: the facades it reflected in, the
// point's side first, how far it went to its last reflection, and where that
// is (the point itself for none).
struct SpecularTrace
{
    std::vector<std::size_t> facades;
    double travelled = 0.0;
    Vec3 exit;
};

// Trace back from `start` along `departure`, in the floor's mirror world for
// `below_floor`, through up to `max_order` reflections; true, with `trace`,
// when it leaves along u.  Below the floor a rising ray can also enter a
// mirrored building through its mirrored roof -- the incident wave coming down
// through a real one -- so there each leg is tested whole.
bool trace_specular(const Buildings& buildings, const Vec3& start, const Vec3& departure, const Vec3& u,
        bool below_floor, std::size_t max_order, SpecularTrace& trace);

// A point mirrored in facades: the image, and the facades it was mirrored in,
// the point's side first.
struct PointImage
{
    Vec3 point;
    std::vector<std::size_t> facades;
};

// The images of a point through which a ray can reach it, found by beam
// tracing in plan view: breadth first over the facades within `radius`, to
// `max_order` reflections and at most `max_images` images, each carrying the
// part of its last facade through which the point can be reached; a further
// facade counts only where it lies in that beam, and of `candidate_edges` an
// image keeps those whose plan crosses it -- an image of more than
// `edge_order` reflections none, for a caller that diffracts nothing so deep.
// The point itself is the first.
struct ImageTree
{
    std::vector<PointImage> images;
    std::vector<std::vector<std::size_t>> edges;
    bool truncated = false; // stopped at max_images
};
ImageTree image_tree(const Buildings& buildings, const Vec3& point, double radius, std::size_t max_order,
        std::size_t max_images, bool beams, const std::vector<std::size_t>& candidate_edges,
        std::size_t edge_order = static_cast<std::size_t>(-1));

// Whether a ray from `from` reaches `point` (or, below the floor, its floor
// image) by reflecting in `facades` (the point's side first): the unfolded
// line crosses each, the last first, within it and from in front, and no leg
// of the folded path passes through a building.  `hits`, when given, receives
// the reflection points in order from `from` (below the floor, in its mirror
// world).
bool folded_path_open(const Buildings& buildings, const Vec3& from, const Vec3& point,
        const std::vector<std::size_t>& facades, bool below_floor, std::vector<Vec3>* hits = nullptr);
// The same, each leg's openness decided by `leg_open(leg, a, b)` instead of
// blocked(): leg 0 leaves `from` and leg facades.size() reaches the point,
// a and b its ends in the folded world (below the floor, its mirror's).
using LegTest = std::function<bool(std::size_t leg, const Vec3& a, const Vec3& b)>;
bool folded_path_open(const Buildings& buildings, const Vec3& from, const Vec3& point,
        const std::vector<std::size_t>& facades, bool below_floor, std::vector<Vec3>* hits, const LegTest& leg_open);

// The same with the ray's one horizontal reflection in the roof of `building`
// instead of the floor: the line from `from` to `point` mirrored in the roof's
// plane and then in `facades` crosses that plane once, on the roof
// (Buildings::on_roof), and every leg of the real path is open.  Among vertical walls a
// ray that leaves a roof keeps rising, so it reflects in one roof at most and
// never in the floor too, and both ends stand above the roof.  `hits` receives
// the reflection points in order from `from`, those past the roof in its
// plane's mirror world (as the floor's are below it), and `bounce` the point
// on the roof.
bool folded_roof_path_open(const Buildings& buildings, const Vec3& from, const Vec3& point,
        const std::vector<std::size_t>& facades, std::size_t building, std::vector<Vec3>* hits = nullptr,
        Vec3* bounce = nullptr);

// One microphone among the buildings, with what does not depend on the wave's
// direction worked out once: the facades' frames and orientations, the edges
// within reach, and the images a diffracted ray can reach it through.
class BuildingReceiver
{
public:
    // `microphone` in the local frame, z above the floor; `run_facades` are
    // indices into buildings.facades().  The buildings must outlive this.
    // A microphone on a wall (Buildings::support) takes the wall's reflection
    // at the microphone as one on the floor takes the floor's: every ray
    // arrives with its partner off the wall, of one length, 1 + R its factor
    // (R the wall's amplitude reflection), halved for a ray diffracted at an
    // edge with that wall as a face, whose coefficient holds the reflection
    // already.  The run must then leave it no first-order reflections of its
    // own (run_facades empty), since they too take the factor.
    BuildingReceiver(const Buildings& buildings, const Vec3& microphone, std::vector<std::size_t> run_facades,
            const BuildingPathOptions& options);

    // Every ray for the plane wave from `u` (unit, toward the source; an
    // elevation under 1e-3 rad is taken at 1e-3, as the canyon does).
    std::vector<BuildingRay> rays(const Vec3& u) const;

    // The table entry: the energy ratio to open ground per band, every ray but
    // the run's summed coherently and band-averaged over three sub-bands (as
    // canyon_ratio does), or in energy without `options.diffraction`, and the
    // facades' diffuse energy added where `diffuse` is given.  `reflection` is
    // the facades' specular energy reflection, [facade][band] for
    // buildings.facades() -- rho (1 - s) where they also scatter; `out` gets
    // one value per band.
    void ratio(const Vec3& u, const std::vector<double>& reflection, double floor_reflection,
            const std::vector<double>& bands, double sound_speed, double* out,
            const DiffuseFacades* diffuse = nullptr) const;

    // The images of the microphone the diffracted rays are tried through, the
    // microphone itself first.
    using Image = PointImage;
    const std::vector<Image>& images() const
    {
        return images_;
    }

    // Whether the image search stopped at options.max_images.
    bool truncated() const
    {
        return truncated_;
    }

    // The wall the microphone lies on, or no_surface.
    static constexpr std::size_t no_surface = static_cast<std::size_t>(-1);
    std::size_t surface() const
    {
        return surface_;
    }

private:
    const Buildings* buildings_;
    Vec3 microphone_;
    std::vector<std::size_t> run_facades_;
    BuildingPathOptions options_;
    std::vector<std::size_t> edges_;    // within reach
    std::vector<std::size_t> planes_;   // the orientations within reach
    std::vector<Image> images_;
    std::vector<std::vector<std::size_t>> image_edges_; // per image, the edges whose rays can reach it
    // Per image, edge of its list and floor (2 j + floor): where the image, or
    // its floor image, stands about the edge -- along it from its start and
    // across from its line -- which every direction's Keller point needs and
    // none changes.  `low` and `high` bound cot(beta) for a Keller point on
    // the edge, and `air` says whether the image can lie in the wedge's air,
    // which does not depend on where along the edge it is seen from (false
    // for an image on the edge's line); both are a hair generous, prefilters
    // for rays()'s own tests.
    struct Foot
    {
        double along = 0.0;
        double across = 0.0;
        double low = 0.0;
        double high = 0.0;
        bool air = false;
    };
    std::vector<std::vector<Foot>> image_feet_;
    bool truncated_ = false;
    std::size_t surface_ = no_surface;
    // A ray's factor from the surface, R its amplitude reflection: 1 + R, or
    // (1 + R) / 2 where the ray's last wedge has the surface as a face; 1 off
    // every surface.
    double surface_factor(const BuildingRay& ray, double r) const;
};

// A ray's complex amplitude at wavenumber k, `reflection` being the facades'
// energy reflection coefficients for one band, [facade]: what its
// reflections leave of it (reflected_amplitude), which is the same at every
// wavenumber of the band, times its wedge's coefficient.
std::complex<double> ray_amplitude(const BuildingRay& ray, const std::vector<double>& reflection,
        double floor_reflection, double k);
double reflected_amplitude(const BuildingRay& ray, const std::vector<double>& reflection, double floor_reflection);
std::complex<double> ray_amplitude(const BuildingRay& ray, double reflected, double k);

// A receiver's table over `axes`, [azimuth][elevation][band] in dB (floored
// at urban_gain_floor_db), its azimuths compass bearings and the local +x
// along `frame_bearing`; the directions are shared out over `threads`.
// `reflection` is [facade][band] for the buildings' facades, specular, and
// `diffuse` the facades' diffuse reflection, if any.
//
// Each entry is evaluated `azimuth_offset` (rad) clockwise of its azimuth.  A
// direction exactly parallel to a wall, or exactly through a corner, is a
// knife edge of the model -- the wall's grazing reflection and its corner's
// diffraction are there a hair to one side and not on the line -- and a table
// of buildings on a grid samples those lines; 1e-6 rad is a tenth of a
// millimetre in 100 m, and far more than a projection's round-off, so the
// same district read from different files lands on the same side of every
// one.
std::vector<double> building_gains(const BuildingReceiver& receiver, const UrbanTable::Axes& axes,
        double frame_bearing, const std::vector<double>& reflection, double floor_reflection, double sound_speed,
        unsigned threads = 1, const DiffuseFacades* diffuse = nullptr, double azimuth_offset = 1e-6);

// The facades the run evaluates for a microphone: those within `radius`
// (horizontal, m) whose outer face looks at it.
std::vector<std::size_t> facing_facades(const Buildings& buildings, const Vec3& microphone, double radius);

} // namespace niceops
