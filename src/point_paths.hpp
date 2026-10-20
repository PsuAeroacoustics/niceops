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
//     microphone's images to `max_twice_reflections` each.  The two diffraction
//     points are where the path is stationary: a fixed point of meeting
//     the single-edge condition at each in turn, found by Steffensen's
//     method (in closed form for parallel edges), and a pair whose point
//     does not converge is no ray; the ray between the edges reflects in
//     nothing, or with max_between_reflections in one facade or the floor,
//     the second edge mirrored in it.  Diffracted three times: the cascade of three
//     coefficients, the middle edge's point a fixed point of meeting the
//     single-edge condition there with the other two at theirs (see
//     thrice_diffracted_ray), the ends' images to max_thrice_reflections and
//     nothing between the edges; three rows of roofs, or the partner a ray
//     that reflects between two edges needs at its facade's end.
//     Two parallel edges that share a face -- a flat
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
// (image_tree).  Left out: four edges or more, more than one reflection
// between two edges, any reflection between three,
// a shared face's aperture on a ray of three, diffuse reflection, and an
// edge's ends.

#include <array>
#include <complex>
#include <functional>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "atmosphere.hpp"
#include "building_paths.hpp"
#include "buildings.hpp"
#include "edge_source.hpp"
#include "double_wedge.hpp"
#include "vec3.hpp"
#include "wedge.hpp"

namespace niceops {

struct PointPathOptions
{
    std::size_t max_reflections = 4;      // on a ray that diffracts nowhere
    std::size_t max_edge_reflections = 2; // on each side of a single diffraction
    std::size_t max_twice_reflections = 2; // on each side of a double diffraction
    // Reflections between the two edges of a double diffraction, in a facade
    // or the floor: 0 or 1.  The ray unfolds about the facade or the floor,
    // the second edge mirrored in it.
    std::size_t max_between_reflections = 1;
    // Rays diffracted at three edges in turn, by the cascade of their three
    // coefficients, the source's and the microphone's images to
    // max_thrice_reflections each and nothing between the edges: over three
    // roof edges (`roofs`, the rows a street's sound crosses), over any three
    // (`all`), or none.  Any three take every chain of corners round a
    // block, thousands in a district and weak: eight times a street run's cost
    // for no level that moved, where roofs add a few rays.
    enum class Triple
    {
        none,
        roofs,
        all
    };
    Triple triple_diffraction = Triple::roofs;
    std::size_t max_thrice_reflections = 0;
    bool diffraction = true;
    bool double_diffraction = true;
    // The uniform form over a shared face (double_wedge.hpp); off, the
    // cascade, for comparison.
    bool shared_face = true;
    // Reflections and the floor between the source and an edge.  Off, an
    // edge diffracts only what reaches it straight from the source, as in the
    // plane-wave tables (building_paths.hpp), for comparison with them.
    bool source_reflections = true;
    // Roofs reflect, as the floor does.  Among vertical walls a ray reflects
    // in one horizontal surface at most -- once it leaves a roof it keeps
    // rising -- so a roof below both ends of the ray's stretch takes the
    // floor's place: on a specular ray, and on either side of a once-
    // diffracted one, where it counts as one of max_edge_reflections and the
    // edge must stand above it (the roof's own edges carry its reflection in
    // their coefficients already).  Off, for comparison.
    bool roofs = true;
    // Facades and edges within this horizontal distance of the source or the
    // microphone (m).
    double radius = 150.0;
    std::size_t max_images = 20000;
    // Rule out diffracted paths before solving or checking them: where an
    // image cannot see the point of the edge the path would take (from the
    // proof in edge_visibility.hpp), and for twice-diffracted pairs where the
    // source's image cannot light the first edge, or the microphone's the
    // second, or no point of either edge can put the other's stationary point
    // on the stretch of it its end sees.  Off, every path is solved and
    // checked: slower, and the reference the tests hold these to (they must
    // find the same rays).
    bool prune = true;
    // Make the diffracted partner the orders above do not trace at a
    // boundary (THEORY 6.5), for edges within this detour of a leg, m: 0 for
    // none.  Specular and once-diffracted rays; a twice-diffracted ray's
    // partner would need a third edge.
    double fade_reach = 0.0;
    // The partner's coefficient with the incident wave's terms only.
    bool fade_incident = false;
    // Make partners even where they are traced, for holding them to those.
    bool fade_traced = false;
    // Diffraction by the edge-source integral (edge_source.hpp) over each
    // edge's whole length, where UTD takes it infinite: once-diffracted rays
    // always, twice-diffracted ones where the cascade cannot be trusted (a
    // transition, an edge's end or a partner's boundary a few Fresnel zones
    // off), blended into it as each wavenumber leaves that.  CPU only.
    bool edge_source = false;
    // With edge_source, once-diffracted rays whose Keller point lies past an
    // end of the edge are kept while the path through the end is no longer
    // than this past the Keller path (m): the end's wave, which the integral
    // over the edge carries, so that a ray leaving the edge fades rather than
    // stops.  point_field sets it to 100 radians at the lowest band.
    double edge_source_reach = 10.0;
    // Over a shared face, the aperture's integral only where the cascade of
    // the two wedges cannot stand in for it (THEORY 6.5): at wavenumbers two
    // thirds of an octave apart, below where the first wedge's transition
    // parameter at the second reaches 30, the ratio of the two taken between
    // them in log k and the cascade itself beyond.  Within 0.1 dB of the
    // integral at every wavenumber.  Off, the integral at every wavenumber:
    // the reference, and what the GPU takes.
    bool sample_shared_faces = true;
};

// The diffracted ray the model does not trace at a boundary, made from the
// ray it partners (THEORY 6.5): one leg of that ray diffracted at an edge it
// passes, the rest kept, by the cascade.  The edge's coefficient (its angles
// at the edge, sin(beta) and distance parameter, as WedgeTerm's) and the
// spreading's ratio, the partner's over the ray's; its detour, which its
// length carries, and the fade's reach, beyond which it is tapered off.
// `incident` keeps only the incident wave's pair of terms.
struct MadeEdge
{
    WedgeAngles angles;
    double sin_beta = 1.0;
    double distance = 0.0;
    double ratio = 1.0;
    double detour = 0.0;
    double reach = 0.0;
    bool incident = false;
};
// The partner's factor on the ray's amplitude at wavenumber k: the edge's
// coefficient, the spreading's ratio and the taper.
std::complex<double> made_edge_factor(const MadeEdge& screen, double k);
// The fade's reach for a run whose lowest wavenumber is k_min: where the
// Fresnel parameter sqrt(2 k d / pi) reaches 3 there, at most `cap` m.
double fade_reach(double k_min, double cap = 0.5);

// One wedge on a ray: its exterior angle (n pi) and the departure and arrival
// angles in its own angle, as wedge_angles(n, theta, theta0) -- the part of
// its coefficient that every wavenumber the ray is taken at shares --
// sin(beta) of the ray to the edge, the distance parameter in the plane
// normal to the edge, L / sin(beta), and whether the incident wave grazes a
// face.
struct WedgeTerm
{
    WedgeAngles angles;
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

// One of the three sub-bands, equal in log frequency, that a band's level
// averages over: its frequencies, Hz, and wavenumbers, 1/m.
struct SubBand
{
    double g1 = 0.0, g2 = 0.0, k1 = 0.0, k2 = 0.0;
};
std::array<SubBand, 3> sub_bands(double band, double sound_speed);
// A band's level, dB, from its sub-bands' band-averaged energies, floored at
// an energy of band_level_floor.
constexpr double band_level_floor = 1e-30;
double band_level(double band, const std::array<SubBand, 3>& parts, const double* energy);
// The air's attenuation of a ray's amplitude per metre of its length, Np/m,
// at each sub-band's middle frequency (ISO 9613-1): [3 band + part].  A ray's
// amplitude in a sub-band is multiplied by exp(-a L).
std::vector<double> sub_band_absorption(const Atmosphere& air, const std::vector<double>& bands, double sound_speed);

// A ray of the field: its length (the phase is k times it), its geometric
// spreading (1/r for a specular ray, 1/sqrt(s s' (s + s')) once diffracted,
// 1/sqrt(s1 s2 s3 (s1 + s2 + s3)) twice through the cascade), its reflections
// and its wedges.  A ray over a shared face carries the face instead of its
// two wedge terms, and its spreading is 1/sqrt(the in-plane length).
// A diffracted ray as the edge-source integral (edge_source.hpp) takes it:
// its edges, in order, the images of the source and the microphone it is
// unfolded between (reflected in the facades and the floor it meets before
// the first edge and after the last), and the stretch of each edge outside
// which its end is proved blind to it.
struct EdgeSourceTerm
{
    std::size_t edges[2] = {0, 0};
    Vec3 source, microphone;
    Stretch stretches[2];
};

struct PointRay
{
    double length = 0.0;
    double spreading = 0.0;
    std::vector<std::size_t> facades;
    int floor = 0;
    int wedges = 0;
    WedgeTerm wedge[3];
    bool shared = false;
    SharedFaceTerm shared_face;
    // Over a shared face, the cascade's spreading, which `spreading` then
    // replaces: the cascade of its two wedge terms (`wedge`, kept) is what
    // the face's field is taken relative to where it is sampled.
    double cascade_spreading = 0.0;
    // For a diffracted ray (wedges > 0).
    EdgeSourceTerm edge_source;
    // The facade a twice-diffracted ray reflects in between its edges, the
    // last of its `facades` before the microphone's; no_facade for none; or
    // whether it reflects in the floor there (one of its `floor`).  Its wedge
    // terms are taken with the second edge mirrored, and it takes UTD under
    // the edge source, whose integral runs over the real edge.
    std::size_t between = no_facade;
    bool between_floor = false;
    // The height of the horizontal plane the ray reflects in on its source
    // side -- anywhere along it for a specular ray, before its first edge for
    // a diffracted one -- 0 for the floor, NaN for none (point_ray_launch).
    double source_plane = std::numeric_limits<double>::quiet_NaN();
    // For each roof it reflects in, a facade of the roof's building, whose
    // reflection coefficient the roof takes (the building's walls and roof
    // absorb alike); no_facade for a building with none, a rigid roof.
    static constexpr std::size_t no_facade = Buildings::no_facade;
    std::vector<std::size_t> roof_facades;
    // A partner made by the fade: the edges it diffracts at besides the
    // ray's own (its length includes their detours).
    std::vector<MadeEdge> made_edges;
    // Where it goes, for drawing: the source, its reflection and diffraction
    // points in order, the microphone.  A leg off the floor is unfolded below
    // it, so only x and y are the ray's own there.
    std::vector<Vec3> path;
};

// One end of the paths, a source or a microphone, with what does not depend
// on the other end, for a run that pairs it with many (a moving source): its
// image tree.  `partners` are the edges near every end it will be paired
// with -- edges_near over their points, at options.radius -- and its images
// keep those and its own near edges (`candidates`) that cross their beams.  A
// pair then keeps the ones near it, which are the lists the pair's own tree
// would have; a pair whose edges an end was not traced with is refused.
// Built and used with the same options.
class PathEnd
{
public:
    // Throws std::invalid_argument where `partners` is not ascending without
    // repeats.
    PathEnd(const Buildings& buildings, const Vec3& point, bool source, const PointPathOptions& options,
            const std::vector<std::size_t>& partners);

    const Vec3& point() const
    {
        return point_;
    }
    // The edges within options.radius of the point, and those and the
    // partners', ascending.
    const std::vector<std::size_t>& near() const
    {
        return near_;
    }
    const std::vector<std::size_t>& candidates() const
    {
        return candidates_;
    }
    const ImageTree& tree() const
    {
        return tree_;
    }
    // How an image, or its mirror in the floor, stands about one of its edges,
    // which does not depend on the partner: its parameter along the edge's
    // line and its distance off it, which place keller's point; the stretch
    // [from, to] of the edge outside which the end is proved not to see it
    // along the image's path (edge_visibility.hpp; the whole edge where
    // nothing is proved, and always without options.prune), and whether it
    // sees any of it; and whether it also lies in the wedge's air, strictly,
    // a hair wider than in_air for the angle taken from the edge's start.
    struct View
    {
        double along = 0.0;
        double off = 0.0;
        double from = 0.0;
        double to = 0.0;
        bool visible = true;
        bool lit = false;
    };
    // Per image, [2 k + floor] for the image's k-th edge, floor 1 the mirror.
    const std::vector<std::vector<View>>& views() const
    {
        return views_;
    }
    // Which end it is, and the options its tree was traced with.
    bool source() const
    {
        return source_;
    }
    const PointPathOptions& options() const
    {
        return options_;
    }

private:
    Vec3 point_;
    bool source_ = false;
    PointPathOptions options_;
    std::vector<std::size_t> near_;
    std::vector<std::size_t> candidates_;
    ImageTree tree_;
    std::vector<std::vector<View>> views_;
};

class PointSourceField
{
public:
    // Source and microphone in the local frame, z above the floor.  The
    // buildings must outlive this.
    PointSourceField(const Buildings& buildings, const Vec3& source, const Vec3& microphone,
            const PointPathOptions& options = {});

    // The same field from two ends traced beforehand, the same to the bit;
    // throws std::invalid_argument where the ends are not a source and a
    // microphone traced with these options, or one was not traced with the
    // edges near the other.
    PointSourceField(const Buildings& buildings, const PathEnd& source, const PathEnd& microphone,
            const PointPathOptions& options = {});

    const std::vector<PointRay>& rays() const
    {
        return rays_;
    }
    // Whether its shared faces are sampled (PointPathOptions::sample_shared_faces).
    bool samples_shared_faces() const
    {
        return sample_shared_faces_;
    }

    // The band levels: 10 lg of |p|^2 band-averaged (every ray summed
    // coherently, pair by pair, over three sub-bands, as the urban tables
    // are), for a unit point source, |p| = 1/r in free field: dB re the
    // source's level at 1 m.  `reflection` is the facades' energy reflection
    // coefficients, [facade][band].  `absorption` is the air's, per metre of
    // path (sub_band_absorption), or empty for none.  `weights`, when given,
    // multiplies each ray's amplitude in each band, [ray][band] in rays()'s
    // order: a directional source's amplitude along the ray's launch
    // (point_ray_launch), the same at every wavenumber of the band.
    void band_levels(const std::vector<double>& reflection, double floor_reflection, const std::vector<double>& bands,
            double sound_speed, double* out, const std::vector<double>& absorption = {},
            const double* weights = nullptr) const;

    // What band_levels sums, for a caller that sums it elsewhere (a GPU):
    // per band and sub-band, in order, the rays' amplitudes at the sub-band's
    // middle wavenumber and their lengths, ascending, with the sub-band.
    // band_level() turns the three sub-bands' energies into the band's level.
    // `order`, when given, gets each amplitude's index in rays() before the
    // first problem.
    using SubBandProblem = std::function<void(std::size_t band, int part, const std::vector<std::complex<double>>& amplitude,
            const std::vector<double>& length, const SubBand& sub_band)>;
    void sub_band_problems(const std::vector<double>& reflection, double floor_reflection, const std::vector<double>& bands,
            double sound_speed, const SubBandProblem& problem, const std::vector<double>& absorption = {},
            const double* weights = nullptr, std::vector<std::size_t>* order = nullptr) const;

private:
    // An end's images and, per image, the edges the pair keeps, where each
    // stands in the image's own list (`at`), and the end's views of them.
    struct Ends
    {
        const std::vector<PointImage>& images;
        const std::vector<std::vector<std::size_t>>& edges;
        const std::vector<std::vector<std::uint32_t>>& at;
        const std::vector<std::vector<PathEnd::View>>& views;
        // The view of the image's kept edge `kept`, from the floor or not.
        const PathEnd::View& view(std::size_t image, std::size_t kept, int floor) const
        {
            return views[image][2 * at[image][kept] + floor];
        }
    };
    void trace(const Buildings& buildings, const PathEnd& source, const PathEnd& microphone,
            const PointPathOptions& options);

    std::vector<PointRay> rays_;
    std::size_t facade_count_ = 0;
    const Buildings* buildings_ = nullptr;
    bool edge_source_ = false;
    bool sample_shared_faces_ = true;
};

// An edge as the edge-source integral takes it: a vertical edge standing on
// the floor (or ending on it) goes on below it as its mirror image, so the
// kept ray of each floor pair takes both halves and the floor is no end.
Edge through_floor(const Edge& g);

// The diffraction factor a ray's edge-source integral gives it at each
// wavenumber, what point_ray_amplitude's wedges give it at reflected = 1:
// its field over its spreading, the phase of its length taken out (a sub-
// band problem's lengths carry that).  Empty for a twice-diffracted ray the
// cascade is trusted on at every wavenumber, and for a ray that diffracts
// nowhere.  (A ray's share of the energy is no guide to skipping it: in
// DigiCity hundreds of twice-diffracted rays each 45 dB down move a
// microphone by half a dB together.)
std::vector<std::complex<double>> edge_source_factors(const Buildings& buildings, const PointRay& ray,
        const std::vector<double>& wavenumbers);

// A twice-diffracted ray's cascade from the source's image pa by g1 at q1 and
// g2 at q2 to the microphone's image pb: its length, spreading and wedge
// terms; the rest of the record is the caller's.
PointRay twice_diffracted_ray(const Edge& g1, const Edge& g2, const Vec3& pa, const Vec3& pb, const Vec3& q1,
        const Vec3& q2);
// The same through three edges.  Each edge's diffracted wave spreads from its
// own caustic at the edge and from a second along it, at the incident
// wavefront's radius of curvature along the edge (Kouyoumjian and Pathak's
// rho_e), which Euler's formula gives from the incident wave's two: s1 for the
// first edge, then for each next one rho' + s and s, the previous edge's
// rho_e and the leg between them, weighted by the cosine between the two
// edges as the ray between them sees them.  Parallel edges give
// 1 / sqrt(s1 s2 s3 s4 (s1 + s2 + s3 + s4)), and two edges the twice-diffracted
// spreading.  The distance parameters are the point-source ones, s s' / (s + s')
// for the legs either side of each edge.
PointRay thrice_diffracted_ray(const Edge& g1, const Edge& g2, const Edge& g3, const Vec3& pa, const Vec3& pb,
        const Vec3& q1, const Vec3& q2, const Vec3& q3);

// The unit direction a ray leaves the source in, which a directional source
// is read along.  A specular ray leaves toward its first hit, or toward the
// microphone or its floor image where it hits nothing; a diffracted one
// toward the first point on its source side, mirrored in the floor where that
// side reflects in it (its source image lies below the floor): the path holds
// a leg off the floor unfolded below it.  The ray's length is stationary in
// every point between, so this is minus its gradient in the source's
// position, which is how it is tested.
Vec3 point_ray_launch(const PointRay& ray);

// A ray's complex amplitude at wavenumber k; `reflection` is the facades'
// energy reflection coefficients for one band, [facade].  It is what the
// ray's spreading and reflections leave (reflected_amplitude), the same at
// every wavenumber of the band, times what its wedges do at k.
std::complex<double> point_ray_amplitude(const PointRay& ray, const std::vector<double>& reflection,
        double floor_reflection, double k);
double reflected_amplitude(const PointRay& ray, const std::vector<double>& reflection, double floor_reflection);
std::complex<double> point_ray_amplitude(const PointRay& ray, double reflected, double k);
// Where a shared face's integral is taken when it is sampled
// (PointPathOptions::sample_shared_faces) over the wavenumbers k: `probes`,
// the indices of those it is taken at, in ascending k, and `free`, the index
// of the lowest at which the cascade stands alone (k.size() for none).  The
// field at the rest is the cascade's times the integral's ratio to it,
// linear in log k in its log magnitude and unwrapped phase between the
// probes and 1 at `free` and above (THEORY 6.5).
struct SharedFaceSamples
{
    std::vector<std::size_t> probes;
    std::size_t free = 0;
};
SharedFaceSamples shared_face_samples(const SharedFaceTerm& face, const std::vector<double>& k);

// point_ray_amplitude for a ray over a shared face, given its face's field at
// k sin(beta) (shared_face_field, or shared_face_spectrum for many k at once).
std::complex<double> shared_face_amplitude(const PointRay& ray, double reflected, double k, std::complex<double> field);

} // namespace niceops
