#pragma once

#include <cstddef>
#include <limits>
#include <vector>

#include "terrain.hpp"

namespace niceops {

// Which terrain screening a run takes (--terrain_screening).
enum class TerrainScreening
{
    none,
    cnossos,
    kirchhoff,
};

// Terrain screening (--terrain_screening cnossos): the CNOSSOS-EU diffraction
// of a path over a terrain profile, in the vertical plane through the
// microphone and the source's footprint (Directive (EU) 2015/996, Annex II
// 2.5.6, as amended by Delegated Directive (EU) 2021/1226; the JRC reference
// report, Kephalopoulos, Paviotti & Anfosso-Ledee 2012, VI.4.4).  The path is
// the rubber band over the profile -- the upper convex hull of the
// microphone, the samples between, and the source -- and the loss is
//
//   A = Delta_dif(S,R) + Delta_ground(S,O) + Delta_ground(O,R),
//
// with the ground on each side of the diffracting edges from NICE-OPS's own
// ground model.  THEORY 6.4 has the equations as run.  Everything here is pure
// and allocation-free once its scratch has grown; process_trajectory holds the
// scratch per thread.

// Delta_dif = 10 lg(3 + 40 C'' delta / lambda) is at least 10 lg 3 whenever
// delta > 0, which is what bounds a screened pair's level for the windowed
// peak search (screen_gain_bound): the one place the 3 is written.
inline constexpr double cnossos_floor = 3.0;
// The cap on Delta_dif(S,R), dB.  It applies to that term only, not to the
// Delta_dif terms inside the ground terms (JRC 2012, VI.4.4.b).
inline constexpr double cnossos_cap = 25.0;

// A point in the vertical plane of a path: s the horizontal distance from the
// microphone's foot, z the height relative to the raster's ground there.
struct ScreenPoint
{
    double s = 0.0;
    double z = 0.0;
};

// The upper convex hull of points[0 .. count), whose s strictly increase, by
// one monotone-chain pass: the indices of its vertices from the first point to
// the last, into `hull` (room for count), and how many.  A point on a segment
// between two others (collinear) is not a vertex.
std::size_t upper_hull(const ScreenPoint* points, std::size_t count, std::size_t* hull);

// The path difference from `from` over `edges` [0, m) in order to `to`: the
// lengths of its segments, each a hypot, summed before the straight length is
// taken off.  Negative, by the same size, when no edge lies strictly above the
// straight line between the ends (the line of sight passes over them), which
// only an image path can give.
double path_difference(const ScreenPoint& from, const ScreenPoint* edges, std::size_t m, const ScreenPoint& to);

// Delta_dif for path difference `delta`, `edges` diffracting edges spanning a
// path length `e` between the first and the last, at wavelength `wavelength`
// (all meters), uncapped: 10 lg(3 + 40 C'' delta / lambda) where
// 40 C'' delta / lambda >= -2, and 0 otherwise.  C'' is 1 for one edge, or for
// several within 0.3 m of each other, and (1 + (5 lambda/e)^2) /
// (1/3 + (5 lambda/e)^2) otherwise.  C_h is 1.
double cnossos_dif(double delta, double e, std::size_t edges, double wavelength);

// A blocked pair's rubber band, and its two sides over their mean lines.
struct Screen
{
    ScreenPoint receiver;            // R, the microphone
    ScreenPoint source;              // S
    std::size_t edges = 0;           // diffracting edges, profile samples
    std::size_t receiver_edge = 0;   // the sample index of the edge nearest R
    std::size_t source_edge = 0;     // and of the edge nearest S
    double delta = 0.0;              // delta(S,R) > 0
    double e = 0.0;                  // path length from the first edge to the last
    // Filled by screen_sides(): delta(S',R) and delta(S,R'), and what the
    // ground models take on each side -- the source and its nearest edge over
    // the source side's line, and the edge nearest the microphone and the
    // microphone over the receiver side's -- with heights below a line taken
    // onto it and distances by their size.
    double delta_source_image = std::numeric_limits<double>::quiet_NaN();
    double delta_receiver_image = std::numeric_limits<double>::quiet_NaN();
    PlaneGeometry source_side;
    PlaneGeometry receiver_side;
};

// The rubber band over a blocked profile: the microphone `microphone_height`
// above the ground at s = 0, samples heights[1 .. interior] every `step`, and
// the source at (distance, source_height) (GroundPlaneCache::Profile).
// `points` and `hull` are scratch.  False if no sample lies above the line of
// sight after all, when `screen` is left alone.
bool rubber_band(const std::vector<double>& heights, std::size_t interior, double step, double distance,
        double microphone_height, double source_height, std::vector<ScreenPoint>& points,
        std::vector<std::size_t>& hull, std::vector<ScreenPoint>& edges, Screen& screen);

// The two sides of `screen` over their mean lines (heights on the profile's
// base, s from the microphone): the source side's fitted from the edge nearest
// the source to the footprint, the receiver side's from the microphone to the
// edge nearest it.  S' and R' are S and R mirrored in them, over the same
// edges.  A flush microphone (`flush`) is on its plate, which already holds
// the reflection under it: its side's source height is the edge's height
// above the microphone, measured perpendicular to the line, and R' is R.
void screen_sides(Screen& screen, const std::vector<ScreenPoint>& edges, const GroundLine& source_line,
        const GroundLine& receiver_line, bool flush);

// What a screened pair's evaluation came to, for the run's log.
struct ScreenTally
{
    bool capped = false;            // Delta_dif(S,R) reached the 25 dB cap in some band
    bool source_clamped = false;    // F_S was held at 1 in some band
    bool receiver_clamped = false;  // and F_R
};

// A per band, dB (a loss), into loss[band] for bands [first, last), at the
// wavelengths sound_speed / band_centers.  `source_gain` and `receiver_gain`
// are each side's ground power ratio G per band (10^(-A_ground/10)), or both
// null without a ground model, when A is Delta_dif(S,R) alone.  Each side's
// weight F = 10^(-(Delta_dif(image) - Delta_dif(S,R))/20) is held to [0, 1];
// past 1 (an image path shorter than the real one, from steeply tilted side
// lines) the Directive's formula has no bound, and this is the one departure
// from it.
void cnossos_screen_loss(const Screen& screen, const std::vector<double>& band_centers, double sound_speed,
        const double* source_gain, const double* receiver_gain, std::size_t first, std::size_t last,
        double* loss, ScreenTally* tally = nullptr);

// --terrain_screening kirchhoff: a blocked pair keeps the level it would have
// unscreened -- its own ground model, plate or table over the mean plane --
// less the Fresnel-Kirchhoff loss of one equivalent knife edge.  The edge
// stands in the aperture plane perpendicular to the line of sight above the
// profile's highest protrusion, at the lowest point of that plane both ends
// see over the terrain: the higher of where the source's horizon ray (through
// the hull vertex nearest it) and the microphone's (through the vertex
// nearest it) cross the plane.  THEORY 6.4 has the equations as run.
struct KnifeEdge
{
    double d1 = 0.0;      // along the line of sight, the source to the aperture plane
    double d2 = 0.0;      // and the microphone to it
    double height = 0.0;  // the edge above the line of sight, perpendicular to it (> 0)
    double delta = 0.0;   // the path difference through the edge, d1 and d2 to it less d1 + d2 (> 0)
};

// The equivalent edge of a blocked pair from rubber_band's `points` (the
// microphone, the profile samples and the source) and `screen`.
KnifeEdge knife_edge(const Screen& screen, const std::vector<ScreenPoint>& points);

// The Fresnel integrals C(x) and S(x), each the integral from 0 to x of the
// cosine and sine of pi t^2 / 2.
void fresnel_integrals(double x, double& c, double& s);

// |F(v)|^2, the power past a straight knife edge relative to the free field
// (the Fresnel-Kirchhoff half-plane): ((1/2 - C(v))^2 + (1/2 - S(v))^2) / 2,
// 1/4 at grazing (v = 0) and falling as 1/(2 pi^2 v^2) deep in the shadow.
double knife_edge_power(double v);

// The knife edge's loss, dB, into loss[band] for bands [first, last): -10 lg
// of |F(v)|^2 averaged uniformly in frequency over each band, from
// band_centers 2^(-1/6) to 2^(1/6), by eight-point Gauss-Legendre, with
// v = 2 sqrt(delta f / sound_speed).  At least 10 lg 4 for an edge above the
// line of sight, so a pair is never louder screened than unscreened.
void knife_edge_loss(const KnifeEdge& edge, const std::vector<double>& band_centers, double sound_speed,
        std::size_t first, std::size_t last, double* loss);

// The most a screened pair's level can exceed free field, as a power ratio,
// given bounds on its two sides' ground gains: max(1, G_S) max(1, G_R) / 3.
// Each side's factor |1 + (sqrt(G) - 1) F|^2 is at most max(1, G) for F in
// [0, 1], and Delta_dif(S,R) >= 10 lg 3.
inline double screen_gain_bound(double source_gain, double receiver_gain)
{
    const double s = source_gain > 1.0 ? source_gain : 1.0;
    const double r = receiver_gain > 1.0 ? receiver_gain : 1.0;
    return s * r / cnossos_floor;
}

} // namespace niceops
