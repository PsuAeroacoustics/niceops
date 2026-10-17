#pragma once

// The field of a point source among LoD1 buildings at a microphone: the
// urban tables' ray model (building_paths.hpp) for a source near the
// buildings, where a plane wave will not do -- a loudspeaker in a street, a
// drone among the roofs (THEORY 6.5).  The rays:
//
//   - specular: the microphone's images in the facades, to `max_reflections`,
//     each with the floor or without;
//   - diffracted once (wedge.hpp, UTD for a spherical wave): at every edge,
//     between an image of the source and an image of the microphone, each to
//     `max_edge_reflections` and each with the floor or without -- so an edge
//     the source lights only by a reflection diffracts too;
//   - diffracted twice: between two edges, the source's and the
//     microphone's images to one reflection each.  The two diffraction
//     points are where the path is stationary, found by meeting the
//     single-edge condition at each in turn; the ray between the edges
//     reflects in nothing.  Two parallel edges that share a face -- a flat
//     roof between its two edges, a wall between two corners -- take the
//     uniform form of double_wedge.hpp, in the plane normal to them, at
//     k sin(beta), and spread along them as the ray does.  Any other pair,
//     or a shared face with both ends above its plane, takes the cascade of
//     two coefficients: the first wedge's wave reaches the second curved
//     about the first edge across it and about the source along it, so the
//     second edge's distance parameter is s2 s3 / (s2 + s3) and its
//     spreading sqrt((s1 + s2) / (s3 (s1 + s2 + s3))).
//
// Every ray is tested exactly (folded_path_open); the images are beam traced
// (image_tree).  Left out: three edges or more, a reflection between two
// edges, diffuse reflection, and an edge's ends.

#include <complex>
#include <cstddef>
#include <vector>

#include "building_paths.hpp"
#include "buildings.hpp"
#include "double_wedge.hpp"
#include "vec3.hpp"

namespace niceops {

struct PointPathOptions
{
    std::size_t max_reflections = 4;      // on a ray that diffracts nowhere
    std::size_t max_edge_reflections = 2; // on each side of a single diffraction
    bool diffraction = true;
    bool double_diffraction = true;
    // The uniform form over a shared face (double_wedge.hpp); off, the
    // cascade, for comparison.
    bool shared_face = true;
    // Reflections and the floor between the source and an edge.  Off, an
    // edge diffracts only what reaches it straight from the source, as in the
    // plane-wave tables (building_paths.hpp), for comparison with them.
    bool source_reflections = true;
    // Facades and edges within this horizontal distance of the source or the
    // microphone (m).
    double radius = 150.0;
    std::size_t max_images = 20000;
};

// One wedge on a ray: its exterior angle (n pi), the departure and arrival
// angles in its own angle, sin(beta) of the ray to the edge, the distance
// parameter in the plane normal to the edge, L / sin(beta), and whether the
// incident wave grazes a face.
struct WedgeTerm
{
    double n = 1.5;
    double theta = 0.0;
    double theta0 = 0.0;
    double sin_beta = 1.0;
    double distance = 0.0;
    bool grazing = false;
};

// A ray diffracted over a shared face: the face and the two ends in the plane
// normal to the edges (double_wedge.hpp), and sin(beta).
struct SharedFaceTerm
{
    SharedFace face;
    double sx = 0.0, sz = 0.0, rx = 0.0, rz = 0.0;
    double sin_beta = 1.0;
};

// A ray of the field: its length (the phase is k times it), its geometric
// spreading (1/r for a specular ray, 1/sqrt(s s' (s + s')) once diffracted,
// 1/sqrt(s1 s2 s3 (s1 + s2 + s3)) twice through the cascade), its reflections
// and its wedges.  A ray over a shared face carries the face instead of its
// two wedge terms, and its spreading is 1/sqrt(the in-plane length).
struct PointRay
{
    double length = 0.0;
    double spreading = 0.0;
    std::vector<std::size_t> facades;
    int floor = 0;
    int wedges = 0;
    WedgeTerm wedge[2];
    bool shared = false;
    SharedFaceTerm shared_face;
    // Where it goes, for drawing: the source, its reflection and diffraction
    // points in order, the microphone.  A leg off the floor is unfolded below
    // it, so only x and y are the ray's own there.
    std::vector<Vec3> path;
};

class PointSourceField
{
public:
    // Source and microphone in the local frame, z above the floor.  The
    // buildings must outlive this.
    PointSourceField(const Buildings& buildings, const Vec3& source, const Vec3& microphone,
            const PointPathOptions& options = {});

    const std::vector<PointRay>& rays() const
    {
        return rays_;
    }

    // The band levels: 10 lg of |p|^2 band-averaged (every ray summed
    // coherently, pair by pair, over three sub-bands, as the urban tables
    // are), for a unit point source, |p| = 1/r in free field: dB re the
    // source's level at 1 m.  `reflection` is the facades' energy reflection
    // coefficients, [facade][band].
    void band_levels(const std::vector<double>& reflection, double floor_reflection, const std::vector<double>& bands,
            double sound_speed, double* out) const;

private:
    std::vector<PointRay> rays_;
    std::size_t facade_count_ = 0;
};

// A ray's complex amplitude at wavenumber k; `reflection` is the facades'
// energy reflection coefficients for one band, [facade].
std::complex<double> point_ray_amplitude(const PointRay& ray, const std::vector<double>& reflection,
        double floor_reflection, double k);

} // namespace niceops
