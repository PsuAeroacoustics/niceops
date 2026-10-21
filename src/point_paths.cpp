#include "point_paths.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "edge_visibility.hpp"
#include "facade_relief.hpp"
#include "faddeeva.hpp"
#include "ground_effect.hpp"
#include "urban.hpp"
#include "util/constants.hpp"
#include "util/parallel.hpp"
#include "wedge.hpp"

namespace niceops {

namespace {

Vec3 below(const Vec3& p)
{
    return {p.x, p.y, -p.z};
}

// A bounce on the floor, counted on the ray with its geometry, on the leg
// from `from` to `to`, one of them the floor's image of its real point: its
// grazing angle's sine and the leg's length.  So no bounce is counted
// without the geometry floor_factors reads.
void add_floor_bounce(PointRay& ray, double sin_grazing, double leg, bool microphone)
{
    ray.floor_bounce.at(static_cast<std::size_t>(ray.floor++)) = {sin_grazing, leg, microphone};
}

void add_floor_bounce(PointRay& ray, const Vec3& from, const Vec3& to, bool microphone)
{
    const double leg = (to - from).norm();
    add_floor_bounce(ray, leg > 0.0 ? std::abs(to.z - from.z) / leg : 1.0, leg, microphone);
}

// How a ray meets a chain of facades (an image's, the point's side first)
// that it crosses toward the point, arriving at the chain's last facade along
// `d`: in the chain's order, appended to `out`.
void meeting_toward(const Buildings& buildings, const std::vector<std::size_t>& chain, Vec3 d,
        std::vector<FacadeHit>& out)
{
    const std::size_t at = out.size();
    out.resize(at + chain.size());
    for (std::size_t k = chain.size(); k-- > 0;)
    {
        const Facade& f = buildings.facades()[chain[k]];
        out[at + k] = facade_hit(f, d);
        d = mirror_direction(f, d);
    }
}

// Fold a facade chain and its floor partner into real coordinates. Unlike
// the older rays' unfolded drawings, a slope joins two independent floor
// images, so all points at that join must be in the real half-space. Return
// the interior points only, with a floor crossing inserted when it is not
// already an endpoint. A floor-contact source keeps its launch separately.
bool physical_chain(const Buildings& buildings, const Vec3& from, const Vec3& to,
        const std::vector<std::size_t>& chain, bool floor, std::vector<Vec3>& points)
{
    thread_local std::vector<Vec3> hits;
    if (!folded_path_open(buildings, from, to, chain, floor, &hits))
        return false;
    hits.push_back(floor ? below(to) : to);
    points.clear();
    Vec3 previous = from;
    for (std::size_t i = 0; i < hits.size(); i++)
    {
        const Vec3& next = hits[i];
        if (previous.z * next.z < 0.0)
        {
            const Vec3 bounce = previous + (previous.z / (previous.z - next.z)) * (next - previous);
            if ((bounce - from).norm() > Buildings::contact && (bounce - to).norm() > Buildings::contact)
                points.push_back(bounce);
        }
        if (i + 1 < hits.size())
            points.emplace_back(next.x, next.y, std::abs(next.z));
        previous = next;
    }
    return true;
}

// Where the path through a point on a line is stationary, the point leaving
// toward b at the angle a arrives at: from a's and b's places about the line
// (along it, and their distances off it), the parameter along it.  False
// where either lies on the line.
bool keller_point(double ta, double ra, double tb, double rb, double& t)
{
    if (ra < 1e-9 || rb < 1e-9)
        return false;
    t = ta + (tb - ta) * ra / (ra + rb);
    return true;
}

// An edge mirrored in a facade's plane, its faces' directions with it, so
// that its own angles (Edge::angle reads face0 and air) are unchanged.
Edge mirrored_edge(const Facade& facade, const Edge& g)
{
    Edge m = g;
    m.start = mirror_in(facade, g.start);
    m.direction = mirror_direction(facade, g.direction);
    m.face0 = mirror_direction(facade, g.face0);
    m.air = mirror_direction(facade, g.air);
    m.air1 = mirror_direction(facade, g.air1);
    return m;
}

// The same in the floor's plane.
Edge floor_mirrored_edge(const Edge& g)
{
    Edge m = g;
    m.start = below(g.start);
    m.direction = below(g.direction);
    m.face0 = below(g.face0);
    m.air = below(g.air);
    m.air1 = below(g.air1);
    return m;
}

// The component of v normal to the edge.
Vec3 normal_part(const Edge& e, const Vec3& v)
{
    return v - v.dot(e.direction) * e.direction;
}

// 1 at or below lo, 0 at or above hi, a cosine-squared taper in log x between.
double trust_taper(double x, double lo, double hi)
{
    if (!(x > lo))
        return 1.0;
    if (!(x < hi))
        return 0.0;
    const double c = std::cos(0.5 * pi * std::log(x / lo) / std::log(hi / lo));
    return c * c;
}

// A once-diffracted ray's Keller point, past an end of its edge, is kept
// while the path through the end, q, exceeds the Keller path by `reach` at
// most and the point is above the floor (below it the point is its mirror
// ray's, whose edge goes on there: continued).
bool end_within_reach(const Vec3& pa, const Vec3& pb, const Vec3& q, const Vec3& keller, double reach)
{
    return keller.z >= 0.0 && (pa - q).norm() + (pb - q).norm() - (pa - keller).norm() - (pb - keller).norm() <= reach;
}

// The edge's ends (EdgeEnds) as the ray from pa to pb diffracted at its
// Keller point t sees them, the edge taken through the floor.  Past an end the
// ray is the end's -- its length, path and launch through the end, as the
// edge-source integral takes it -- so its factor is put relative to that
// (EdgeEnds::scale).
EdgeEnds edge_ends(const Edge& edge, const Vec3& pa, const Vec3& pb, double t)
{
    const Edge g = continued(edge);
    const double tg = t + (edge.start - g.start).dot(g.direction);
    const Vec3 lo = g.start, hi = g.start + g.length * g.direction;
    const Vec3 keller = edge.start + t * edge.direction;
    const double length = (pa - keller).norm() + (pb - keller).norm();
    EdgeEnds e;
    if (std::isfinite(g.beyond[0]))
        e.lo = (pa - lo).norm() + (pb - lo).norm() - length;
    if (std::isfinite(g.beyond[1]))
        e.hi = (pa - hi).norm() + (pb - hi).norm() - length;
    e.past = tg < 0.0 ? -1 : tg > g.length ? 1 : 0;
    if (e.past != 0)
    {
        // The ray takes the Keller point's coefficient (its sin(beta) and
        // distance), and its spreading over the end's here, so that past the
        // end it is the Keller ray times the factor, phase and all, and
        // reciprocal.
        const Vec3 end = edge.start + std::clamp(t, 0.0, edge.length) * edge.direction;
        const auto spreading = [&](const Vec3& x) {
            const double s_in = (pa - x).norm(), s_out = (pb - x).norm();
            return 1.0 / std::sqrt(s_in * s_out * (s_in + s_out));
        };
        e.scale = spreading(keller) / spreading(end);
    }
    return e;
}

// Past an end the end's wave reaches wherever the end itself sees, which the
// wedge's own faces, infinite, do not decide: the end is a vertex, and a
// microphone beyond a face's plane may see it over the face's far edge.  Its
// coefficient there is the face's, the angle taken onto the nearer face.
double onto_air(double theta, double n)
{
    if (theta <= n * pi)
        return theta;
    return theta - n * pi < 2.0 * pi - theta ? n * pi : 0.0;
}

// How deep an angle lies in the wedge's solid side: 0 in its air.
double into_solid(double theta, double n)
{
    return theta <= n * pi ? 0.0 : std::min(theta - n * pi, 2.0 * pi - theta);
}

// A ray past an end with its source or microphone in the wedge's solid side
// (onto_air): the end wave there is the end's, which has to meet two limits.
// At the plane of the face it is the air side's, and at the end -- where the
// ray inside the edge, its other side, is no ray in the solid -- nothing.  The
// weight 1 / sqrt(1 + L 2 sin^2(depth / 2) / Delta) meets both, the depth's
// transition distance against the end's excess (as the coefficient's terms
// take their ends, ended_wedge_coefficient), and folds into the ends' scale.
void weigh_vertex(EdgeEnds& e, double theta0, double theta, double n, double s_in, double s_out)
{
    const double depth = std::max(into_solid(theta0, n), into_solid(theta, n));
    if (e.past == 0 || !(depth > 0.0))
        return;
    const double near = std::max(e.past < 0 ? e.lo : e.hi, 1e-300);
    const double sine = std::sin(0.5 * depth);
    e.scale /= std::sqrt(1.0 + s_in * s_out / (s_in + s_out) * 2.0 * sine * sine / near);
}

// The least path from p to q by a point of the edge's line: unfolded about it.
double via_line(const Edge& e, const Vec3& p, const Vec3& q)
{
    const Vec3 dp = p - e.start, dq = q - e.start;
    const double tp = dp.dot(e.direction), tq = dq.dot(e.direction);
    return std::hypot((dp - tp * e.direction).norm() + (dq - tq * e.direction).norm(), tq - tp);
}

// A twice-diffracted ray's edges' ends (EdgeEnds), from pa by g1 and g2 to pb,
// `length` long: each end's excess with the other edge's point free on its
// line, so that the coupling of the two (A - C^2 / B, end_ratios) is in it,
// and the first edge's at the second edge's shadow boundary is the once-
// diffracted ray's there, whose jump the pair carries.
void twice_ends(const Edge& g1, const Edge& g2, const Vec3& pa, const Vec3& pb, double length, EdgeEnds& e1,
        EdgeEnds& e2)
{
    const Edge f1 = continued(g1), f2 = continued(g2);
    const Vec3 lo1 = f1.start, hi1 = f1.start + f1.length * f1.direction;
    const Vec3 lo2 = f2.start, hi2 = f2.start + f2.length * f2.direction;
    const double none = std::numeric_limits<double>::infinity();
    e1.lo = std::isfinite(f1.beyond[0]) ? (pa - lo1).norm() + via_line(g2, lo1, pb) - length : none;
    e1.hi = std::isfinite(f1.beyond[1]) ? (pa - hi1).norm() + via_line(g2, hi1, pb) - length : none;
    e2.lo = std::isfinite(f2.beyond[0]) ? via_line(g1, pa, lo2) + (lo2 - pb).norm() - length : none;
    e2.hi = std::isfinite(f2.beyond[1]) ? via_line(g1, pa, hi2) + (hi2 - pb).norm() - length : none;
}

// A chain-diffracted ray's edges' ends (EdgeEnds), from pa by the edges g at
// points q to pb, `length` long: each end's excess with every other edge's
// point free on its line, as twice_ends takes a pair's, so that a chain's
// ends meet its own pair's where an edge of it leaves the path.
void chain_ends(const std::vector<const Edge*>& g, const Vec3& pa, const Vec3& pb, const std::vector<Vec3>& q,
        double length, PointRay& ray)
{
    const std::size_t n = g.size();
    // The least path from a over the lines g[from, to) to b, from the ray's
    // own points; by the points themselves where that does not settle.
    const auto least = [&](std::size_t from, std::size_t to, const Vec3& a, const Vec3& b) {
        if (from == to)
            return (b - a).norm();
        if (to == from + 1)
            return via_line(*g[from], a, b);
        const std::vector<const Edge*> lines(g.begin() + static_cast<std::ptrdiff_t>(from),
                g.begin() + static_cast<std::ptrdiff_t>(to));
        std::vector<double> t(lines.size());
        for (std::size_t i = 0; i < lines.size(); i++)
            t[i] = (q[from + i] - lines[i]->start).dot(lines[i]->direction);
        chain_points(lines, a, b, t, false);
        double l = 0.0;
        Vec3 p = a;
        for (std::size_t i = 0; i < lines.size(); i++)
        {
            const Vec3 r = lines[i]->start + t[i] * lines[i]->direction;
            l += (r - p).norm();
            p = r;
        }
        return l + (b - p).norm();
    };
    const double none = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < n; i++)
    {
        const Edge f = continued(*g[i]);
        EdgeEnds e;
        const Vec3 lo = f.start, hi = f.start + f.length * f.direction;
        e.lo = std::isfinite(f.beyond[0]) ? least(0, i, pa, lo) + least(i + 1, n, lo, pb) - length : none;
        e.hi = std::isfinite(f.beyond[1]) ? least(0, i, pa, hi) + least(i + 1, n, hi, pb) - length : none;
        (i < 3 ? ray.wedge[i] : ray.more_wedges[i - 3]).ends = e;
    }
}

// Whether an angle lies in the wedge's air, strictly or with its faces.
bool in_air(double theta, double n, bool faces)
{
    const double slack = 1e-9;
    return faces ? (theta > -slack && theta < n * pi + slack) || theta > 2.0 * pi - slack
                 : theta > slack && theta < n * pi - slack;
}

// From the direction `toward` in the plane normal to an edge, how far it is
// to the wedge's face reached turning toward -z, and to the one reached
// turning toward +z (radians).  The wedge's angle grows from face0 toward air.
void faces_either_side(const Edge& e, const Vec3& toward, const Vec3& z, double& below, double& above)
{
    const double a = e.angle(toward);
    const bool grows_up = toward.dot(e.face0) * z.dot(e.air) - toward.dot(e.air) * z.dot(e.face0) > 0.0;
    below = grows_up ? a : e.n * pi - a;
    above = grows_up ? e.n * pi - a : a;
}

// Two parallel edges a ray diffracts at in turn, sharing no face, as
// shared_face_field takes them (double_wedge.hpp): x along the line from the
// first edge (at q1) to the second in the plane normal to them, z off it on
// the side away from the first wedge's nearer face, and each wedge's face on
// the other side of the line tilted below it.  The line must touch both
// wedges from above, as a shared face does -- each wedge's other face half a
// turn or more round from it -- or the aperture's field past the second edge
// takes in waves the rays never would; no face may lie along it; and the two
// faces below it must meet, a step closed beneath the line, as the aperture's
// second wedge takes its face to run back to the first's.  Across a street
// the lower roof stops short and the street opens below.
bool tilted_pair(const Buildings& buildings, const Edge& g1, const Edge& g2, const Vec3& q1, double width,
        const Vec3& x, Vec3& z, SharedFace& face)
{
    double below1, above1, below2, above2;
    faces_either_side(g1, x, z, below1, above1);
    if (below1 > above1)
    {
        z = -1.0 * z;
        std::swap(below1, above1);
    }
    faces_either_side(g2, -1.0 * x, z, below2, above2);
    const double slack = 1e-9;
    if (!(below1 > slack && below2 > slack && above1 >= pi && above2 >= pi && below1 + below2 < pi))
        return false;
    // Where the faces' lines meet, C, and the notch between them and the line:
    // its middle in air, and the solid just behind each face, a quarter, half
    // and three quarters of the way from its edge to C.
    const Vec3 down1 = std::cos(below1) * x - std::sin(below1) * z;
    const Vec3 down2 = -std::cos(below2) * x - std::sin(below2) * z;
    const double reach1 = width * std::sin(below2) / std::sin(below1 + below2);
    const double reach2 = width * std::sin(below1) / std::sin(below1 + below2);
    const Vec3 q2 = q1 + width * x, corner = q1 + reach1 * down1;
    const Vec3 middle = (q1 + q2 + corner) / 3.0;
    if (buildings.buried(middle))
        return false;
    const double behind = 1e-3 * std::min(reach1, reach2);
    for (const auto& [edge, down, reach] : {std::tuple(q1, down1, reach1), std::tuple(q2, down2, reach2)})
        for (double part : {0.25, 0.5, 0.75})
        {
            const Vec3 p = edge + part * reach * down;
            Vec3 out = (p - middle) - (p - middle).dot(down) * down;
            out = out / out.norm();
            if (!buildings.buried(p + behind * out))
                return false;
        }
    face = {width, g1.n, g2.n, below1, below2};
    return true;
}

// An angle on a face (within round-off) taken exactly on it.
double onto_face(double theta, double n, bool& on_face)
{
    const double slack = 1e-7;
    on_face = true;
    if (theta < slack || theta > 2.0 * pi - slack)
        return 0.0;
    if (std::abs(theta - n * pi) < slack)
        return n * pi;
    on_face = false;
    return theta;
}

// The angle of v about the edge, or the face's own where the end lies on face
// `face` (edge_face; -1 none), which it then records in `on_face`.
double end_angle(const Edge& e, const Vec3& v, int face, bool& on_face)
{
    on_face = face >= 0;
    return face == 0 ? 0.0 : face == 1 ? e.n * pi : e.angle(normal_part(e, v));
}

// in_air(e.angle(v), e.n, faces), by the signs where they decide it.
bool air(const Edge& e, const Vec3& v, bool faces)
{
    const int sign = e.air_by_sign(v);
    return sign >= 0 ? sign == 1 : in_air(e.angle(v), e.n, faces);
}

// in_air(onto_face(e.angle(v), e.n, on_face), e.n, true): a direction on a
// face, within onto_face's 1e-7, counts as in the air.
bool air_or_face(const Edge& e, const Vec3& v)
{
    const int sign = e.air_by_sign(v);
    if (sign >= 0)
        return sign == 1;
    bool on_face;
    return in_air(onto_face(e.angle(v), e.n, on_face), e.n, true);
}

// The roof edges a path from a to b over the buildings bends at, in order,
// each with where along it: the corners of the taut string over the
// buildings' cross-section in the vertical plane through the two -- the
// upper hull, between them, of the corners of every roof profile the plan
// segment crosses -- each the roof edge it lies on.  Empty where the string
// bends nowhere, and where a corner lies on no edge in `near` (sorted), or at
// an edge's end: a footprint's vertex, or a building beyond the pair's reach.
std::vector<std::pair<std::size_t, double>> roof_chain(const Buildings& buildings, const Vec3& a, const Vec3& b,
        const std::vector<std::size_t>& near)
{
    std::vector<std::pair<std::size_t, double>> chain;
    const double dx = b.x - a.x, dy = b.y - a.y, span = std::hypot(dx, dy);
    if (!(span > 0.0))
        return chain;
    const double ux = dx / span, uy = dy / span;
    struct Corner
    {
        double d, z;
        std::uint32_t building;
    };
    thread_local std::vector<Buildings::SectionSpan> spans;
    thread_local std::vector<Corner> corners;
    buildings.section(a, b, spans);
    corners.clear();
    for (const Buildings::SectionSpan& s : spans)
    {
        const Footprint& f = buildings.footprints()[s.building];
        if (f.roof.empty())
        {
            corners.push_back({s.d0, f.top, s.building});
            corners.push_back({s.d1, f.top, s.building});
            continue;
        }
        // The roof over the span is straight over each stretch of one least
        // plane, and bends where the next one's slope along it differs.
        auto corner = [&](double d) {
            corners.push_back({d, roof_height(f, a.x + d * ux, a.y + d * uy), s.building});
        };
        auto slope = [&](const std::tuple<double, double, std::size_t>& stretch) {
            const RoofPlane& p = f.roof[std::get<2>(stretch)];
            return p.a * ux + p.b * uy;
        };
        const auto stretches = roof_stretches(f, a.x + s.d0 * ux, a.y + s.d0 * uy, a.x + s.d1 * ux, a.y + s.d1 * uy);
        corner(s.d0);
        for (std::size_t k = 1; k < stretches.size(); k++)
            if (std::abs(slope(stretches[k - 1]) - slope(stretches[k])) > 1e-12)
                corner(s.d0 + std::get<0>(stretches[k]) * (s.d1 - s.d0));
        corner(s.d1);
    }
    if (corners.empty())
        return chain;
    std::sort(corners.begin(), corners.end(),
            [](const Corner& p, const Corner& q) { return p.d < q.d || (p.d == q.d && p.z < q.z); });
    // The upper hull from a to b over them, dropping corners on a straight
    // stretch of it.
    thread_local std::vector<Corner> hull;
    hull.assign({{0.0, a.z, 0}});
    const double flat = 1e-12 * (span + std::abs(a.z) + std::abs(b.z)) * (span + std::abs(a.z) + std::abs(b.z));
    auto push = [&](const Corner& p) {
        while (hull.size() >= 2)
        {
            const Corner& o = hull[hull.size() - 2];
            const Corner& m = hull.back();
            if ((m.d - o.d) * (p.z - o.z) - (m.z - o.z) * (p.d - o.d) > -flat)
                hull.pop_back();
            else
                break;
        }
        hull.push_back(p);
    };
    for (const Corner& c : corners)
        push(c);
    push({span, b.z, 0});
    const auto& edges = buildings.edges();
    for (std::size_t h = 1; h + 1 < hull.size(); h++)
    {
        const Corner& c = hull[h];
        const Vec3 p(a.x + c.d * ux, a.y + c.d * uy, c.z);
        std::size_t best = edges.size();
        double nearest = 1e-6 + 1e-9 * span, along = 0.0;
        for (std::size_t e = buildings.edge_begin(c.building); e < buildings.edge_begin(c.building + 1); e++)
        {
            const Edge& g = edges[e];
            if (g.vertical || !std::binary_search(near.begin(), near.end(), e))
                continue;
            const double t = std::clamp((p - g.start).dot(g.direction), 0.0, g.length);
            const double distance = (g.start + t * g.direction - p).norm();
            if (distance < nearest)
            {
                nearest = distance;
                best = e;
                along = t;
            }
        }
        const double end = 1e-6 + 1e-9 * span;
        if (best == edges.size() || along < end || along > edges[best].length - end
                || (!chain.empty() && chain.back().first == best))
            return {};
        chain.emplace_back(best, along);
    }
    return chain;
}

} // namespace

namespace {

// Whether two edges share an end, within contact: a corner of a building.
bool meet(const Edge& a, const Edge& b)
{
    const Vec3 a1 = a.start + a.length * a.direction, b1 = b.start + b.length * b.direction;
    for (const Vec3& p : {a.start, a1})
        for (const Vec3& q : {b.start, b1})
            if ((p - q).norm() < Buildings::contact)
                return true;
    return false;
}

// Whether any point of edge e lies in edge f's air or on its faces: the solid
// between f's faces is convex (1 < n < 2), so e lies inside it whole only
// where both its ends do.
bool sees_air(const Edge& f, const Edge& e)
{
    return air_or_face(f, normal_part(f, e.start - f.start))
            || air_or_face(f, normal_part(f, e.start + e.length * e.direction - f.start));
}

// How a point stands about an edge (PathEnd::View), seeing [from, to] of it
// or nothing.
// `faces`: a point on a surface, whose leg to the edge may run along a face.
PathEnd::View view_from(const Edge& e, const Vec3& p, bool visible, double from, double to, bool faces = false)
{
    const Edge::Projection place = e.project(p);
    PathEnd::View s;
    s.along = place.along;
    s.off = place.normal.norm();
    s.from = from;
    s.to = to;
    s.visible = visible;
    const int sign = e.air_by_sign(place.normal);
    if (s.off < 1e-4 || sign == 1)
        s.lit = true;
    else if (sign == -1)
    {
        const double theta = e.angle(place.normal);
        const double slack = 1e-9 - 1e-10;
        bool on_face;
        s.lit = faces ? in_air(onto_face(theta, e.n, on_face), e.n, true)
                      : theta > slack && theta < e.n * pi - slack;
    }
    s.lit = s.lit && visible;
    return s;
}

// How one edge's points stand about another edge's line.  The point at
// parameter t of the segment lies along the line at c0 + c1 t and off it by
// |n0 + t nd|; over the segment those range over [along_min, along_max] and
// [off_min, off_max].
struct Side
{
    double c0 = 0.0, c1 = 0.0;
    Vec3 n0, nd;
    double along_min = 0.0, along_max = 0.0;
    double off_min = 0.0, off_max = 0.0;

    // keller_point's places for the point at parameter x.
    void at(double x, double& along, double& off) const
    {
        along = c0 + c1 * x;
        off = (n0 + x * nd).norm();
    }
};

Side side_of(const Edge& line, const Edge& segment)
{
    const Edge::Projection start = line.project(segment.start);
    Side s;
    s.c0 = start.along;
    s.c1 = segment.direction.dot(line.direction);
    s.n0 = start.normal;
    s.nd = normal_part(line, segment.direction);
    const double end = s.c0 + segment.length * s.c1;
    s.along_min = std::min(s.c0, end);
    s.along_max = std::max(s.c0, end);
    const double nn = s.nd.dot(s.nd);
    const double closest = nn > 0.0 ? std::clamp(-s.n0.dot(s.nd) / nn, 0.0, segment.length) : 0.0;
    s.off_min = (s.n0 + closest * s.nd).norm();
    s.off_max = std::max(s.n0.norm(), (s.n0 + segment.length * s.nd).norm());
    return s;
}

// Two edges about each other, for the twice-diffracted search: the second's
// points about the first's line and the first's about the second's.
struct EdgePair
{
    Side on1, on2;
    bool parallel = false;
};

EdgePair edge_pair(const Edge& g1, const Edge& g2)
{
    EdgePair p;
    p.on1 = side_of(g1, g2);
    p.on2 = side_of(g2, g1);
    p.parallel = std::abs(p.on2.c1) > 1.0 - 1e-12;
    return p;
}

// The edge pairs a field has met, by their two edges: worked out once per
// field however many images meet them.  Open addressing, and a generation per
// slot so that the next field starts empty without clearing; its memory
// follows the pairs met, not the square of the edges near.  A reference it
// returns lasts until the next get().
class PairTable
{
public:
    void start()
    {
        pairs_.clear();
        if (++generation_ == 0)
        {
            std::fill(stamp_.begin(), stamp_.end(), 0u);
            generation_ = 1;
        }
    }

    const EdgePair& get(const std::vector<Edge>& edges, std::size_t e1, std::size_t e2)
    {
        if (2 * (pairs_.size() + 1) > key_.size())
            grow();
        const std::uint64_t key = (static_cast<std::uint64_t>(e1) << 32) | static_cast<std::uint64_t>(e2);
        std::size_t i = slot(key);
        while (stamp_[i] == generation_)
        {
            if (key_[i] == key)
                return pairs_[value_[i]];
            i = (i + 1) & (key_.size() - 1);
        }
        stamp_[i] = generation_;
        key_[i] = key;
        value_[i] = static_cast<std::uint32_t>(pairs_.size());
        pairs_.push_back(edge_pair(edges[e1], edges[e2]));
        return pairs_.back();
    }

private:
    std::size_t slot(std::uint64_t key) const
    {
        return static_cast<std::size_t>((key * 0x9E3779B97F4A7C15ull) >> shift_);
    }

    void grow()
    {
        const std::vector<std::uint64_t> key = std::move(key_);
        const std::vector<std::uint32_t> stamp = std::move(stamp_), value = std::move(value_);
        std::size_t size = 1024, bits = 10;
        while (size < 2 * key.size())
        {
            size *= 2;
            bits++;
        }
        shift_ = 64 - static_cast<unsigned>(bits);
        key_.assign(size, 0);
        stamp_.assign(size, 0);
        value_.assign(size, 0);
        for (std::size_t j = 0; j < key.size(); j++)
            if (stamp[j] == generation_)
            {
                std::size_t i = slot(key[j]);
                while (stamp_[i] == generation_)
                    i = (i + 1) & (size - 1);
                stamp_[i] = generation_;
                key_[i] = key[j];
                value_[i] = value[j];
            }
    }

    std::vector<std::uint64_t> key_;
    std::vector<std::uint32_t> stamp_, value_;
    std::vector<EdgePair> pairs_;
    std::uint32_t generation_ = 0;
    unsigned shift_ = 54;
};

// Whether keller's point between p and q can fall within the stretch of the
// line's edge p sees, give or take `slack`, for some q on the segment whose
// Side about the line is `side`: the stationary point lies along the line
// from p at the fraction off_p / (off_p + off_q) of the way to q's foot, so it
// lies within the hull of those fractions of the feet's range -- and,
// whatever the distances, between p's foot and the feet (keller_between, the
// cheaper first test).
bool keller_between(const PathEnd::View& p, const Side& side, double slack)
{
    return std::max(p.along, side.along_max) >= p.from - slack && std::min(p.along, side.along_min) <= p.to + slack;
}

bool keller_can_reach(const PathEnd::View& p, const Side& side, double slack)
{
    if (p.off < 1e-4)
        return true;
    const double near = p.off / (p.off + side.off_min), far = p.off / (p.off + side.off_max);
    const double lo = side.along_min - p.along, hi = side.along_max - p.along;
    const double reach_min = std::min(std::min(lo * near, lo * far), std::min(hi * near, hi * far));
    const double reach_max = std::max(std::max(lo * near, lo * far), std::max(hi * near, hi * far));
    return p.along + reach_max >= p.from - slack && p.along + reach_min <= p.to + slack;
}

// The stationary point between two parallel edges' lines in closed form: the
// path unfolds about the common direction into a straight line, so its
// travel along it is shared in proportion to the distances across.  The
// second edge's parameter; false where the distances leave it undefined.
bool parallel_stationary(const Edge& g1, const Edge& g2, const Vec3& pa, const Vec3& pb, double& t2)
{
    const Edge::Projection a = g1.project(pa), c = g1.project(g2.start);
    const double zb = g1.project(pb).along;
    const double ra = a.normal.norm(), rc = c.normal.norm(), rb = g2.project(pb).normal.norm();
    const double total = ra + rc + rb;
    if (!(total > 0.0))
        return false;
    const double z2 = a.along + (zb - a.along) * (ra + rc) / total;
    t2 = (z2 - c.along) * g2.direction.dot(g1.direction);
    return std::isfinite(t2);
}

std::vector<std::size_t> joined(const std::vector<std::size_t>& a, const std::vector<std::size_t>& b)
{
    std::vector<std::size_t> out;
    out.reserve(a.size() + b.size());
    std::set_union(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));
    return out;
}

// The partners as set_union and includes need them: ascending, no repeats.
const std::vector<std::size_t>& ascending(const std::vector<std::size_t>& partners)
{
    if (std::adjacent_find(partners.begin(), partners.end(), std::greater_equal<std::size_t>()) != partners.end())
        throw std::invalid_argument("PathEnd: the partners' edges must be ascending, without repeats (edges_near)");
    return partners;
}

// The deepest images whose edges a field reads: a single diffraction's, and
// the twice-diffracted search's, which takes each end's to
// max_twice_reflections whatever max_edge_reflections is.
std::size_t edge_order(const PointPathOptions& options)
{
    const bool twice = options.diffraction && options.double_diffraction;
    return std::max({options.max_edge_reflections, twice ? options.max_twice_reflections : std::size_t {0},
            twice && options.triple_diffraction != PointPathOptions::Triple::none ? options.max_thrice_reflections
                                                                                   : std::size_t {0}});
}

// How deep an end's image tree goes: the microphone's to the specular rays'
// order and the edges'; the source's (with source_reflections) to the edges'
// order, as the microphone's side, so the ray set is the same with the ends
// swapped -- and with reliefs to the grating rays' order and one, whose
// images the diffracted searches pass over by their own orders.
std::size_t tree_order(const Buildings& buildings, bool source, const PointPathOptions& options)
{
    const std::size_t gratings = options.reliefs != nullptr ? options.grating_reflections + 1 : 0;
    if (source)
    {
        // Specular wall -> slope paths need source images even with
        // diffraction (or source-side diffraction reflections) disabled.
        const bool slopes = options.roofs && options.slope_facades && options.max_reflections > 1
                && std::any_of(buildings.footprints().begin(), buildings.footprints().end(), [](const Footprint& f) {
                    return std::any_of(f.roof.begin(), f.roof.end(), [](const RoofPlane& p) { return !p.level(); });
                });
        return std::max({options.diffraction && options.source_reflections ? edge_order(options) : 0,
                gratings, slopes ? options.max_reflections - 1 : 0});
    }
    return std::max(options.max_reflections, edge_order(options));
}

// The detour of the path from a to b by the line from p along the unit u:
// the shortest a -> q -> b, q on the line, less |b - a|; q, and its parameter
// t along the line.
double detour(const Vec3& a, const Vec3& b, const Vec3& p, const Vec3& u, Vec3& q, double& t)
{
    const Vec3 da = a - p, db = b - p;
    const double ta = da.dot(u), tb = db.dot(u);
    const double ra = (da - ta * u).norm(), rb = (db - tb * u).norm();
    t = ra + rb > 0.0 ? ta + (tb - ta) * ra / (ra + rb) : 0.5 * (ta + tb);
    q = p + t * u;
    return std::max(0.0, (a - q).norm() + (b - q).norm() - (b - a).norm());
}

constexpr std::size_t no_building = static_cast<std::size_t>(-1);

// An edge a leg passes within the fade's reach, before the ray it would
// partner is known: the edge (or its image below the floor), its building,
// its point where the leg's detour by it is least, the leg's ends (from the
// chain's start), how far along its chain the leg starts, which chain (0 the
// whole of a specular ray, 1 a diffracted ray's source side, 2 its microphone
// side), whether the edge's building stands across the leg, and the detour.
struct PendingEdge
{
    std::size_t edge = 0, building = 0;
    bool image = false;
    Vec3 point, a, b;
    double before = 0.0;
    int chain = 0;
    bool across = false;
    double detour = 0.0;
};

// The edges of the buildings near the leg from a to b (below the floor, in
// its mirror world) that a partner can diffract at: the leg's path is
// stationary within the edge, both of its ends lie in the edge's air, no
// building blocks either of the partner's two legs, and its detour by the
// edge is within the fade's reach.  Each building's are
// taken as across the leg where it blocks the leg and beside it otherwise,
// but for the buildings the leg's ends lie on (skip0, skip1), whose edges the
// leg leaves from; a leg through one of those, past a corner of a footprint
// that is not convex, is closed.  False where a building blocks the leg with
// no such edge.  Below the floor a building's image stands from -top: its
// edges are taken with their images.
bool leg_edges(const Buildings& buildings, const Vec3& a, const Vec3& b, bool below_floor, double reach,
        std::size_t skip0, std::size_t skip1, double before, int chain, std::vector<PendingEdge>& out)
{
    const double length = (b - a).norm();
    // The least detour of an edge d from the leg is at the leg's middle,
    // sqrt(length^2 + 4 d^2) - length; it reaches the fade's reach at
    // d = sqrt(reach^2 + 2 reach length) / 2, and the box is that wide.
    const double h = 0.5 * std::sqrt(reach * reach + 2.0 * reach * length) + Buildings::contact;
    thread_local std::vector<std::uint32_t> near;
    buildings.near_box(std::min(a.x, b.x) - h, std::max(a.x, b.x) + h, std::min(a.y, b.y) - h, std::max(a.y, b.y) + h,
            near);
    const auto& edges = buildings.edges();
    for (const std::uint32_t building : near)
    {
        if (building == skip0 || building == skip1)
        {
            if (buildings.blocked_by(building, a, b, below_floor))
                return false;
            continue;
        }
        const std::size_t first = out.size();
        const bool across = buildings.blocked_by(building, a, b, below_floor);
        for (std::size_t e = buildings.edge_begin(building); e < buildings.edge_begin(building + 1); e++)
        {
            const Edge& edge = edges[e];
            for (int image = 0; image < (below_floor ? 2 : 1); image++)
            {
                // The edge's image: a vertical edge from -top down to its
                // foot's image, a roof edge at -top.
                const Vec3 start = !image ? edge.start
                        : edge.vertical   ? Vec3(edge.start.x, edge.start.y, -(edge.start.z + edge.length))
                                          : Vec3(edge.start.x, edge.start.y, -edge.start.z);
                Vec3 q;
                double t;
                const double d = detour(a, b, start, edge.direction, q, t);
                if (!(d < reach) || t <= 0.0 || t >= edge.length)
                    continue;
                auto real = [&](const Vec3& v) { return image ? Vec3(v.x, v.y, -v.z) : v; };
                if (!air(edge, normal_part(edge, real(a - q)), false) || !air(edge, normal_part(edge, real(b - q)), false))
                    continue;
                // The partner's own legs, as any ray's.
                if (buildings.blocked(a, q, below_floor) || buildings.blocked(q, b, below_floor))
                    continue;
                PendingEdge p;
                p.edge = e;
                p.building = building;
                p.image = image == 1;
                p.point = q;
                p.a = a;
                p.b = b;
                p.before = before;
                p.chain = chain;
                p.across = across;
                p.detour = d;
                out.push_back(p);
            }
        }
        if (across && out.size() == first)
            return false;
    }
    return true;
}

// The partner a pending screen makes of a ray whose chains are `segments`
// long (one for a specular ray, two either side of the edge for a diffracted
// one): the screened chain is cut where the leg passes the edge, the
// spreading 1 / sqrt(s1 s2 ... (s1 + s2 + ...)) taken over the new segments
// against the old, the edge's coefficient with its angles there.
MadeEdge complete(const Buildings& buildings, const PendingEdge& p, const std::vector<double>& segments,
        const PointPathOptions& options)
{
    const Edge& edge = buildings.edges()[p.edge];
    // Directions from the edge's point to the leg's ends, in the real
    // edge's frame: an image's are mirrored back.
    auto toward = [&](const Vec3& v) {
        const Vec3 d = v - p.point;
        return p.image ? Vec3(d.x, d.y, -d.z) : d;
    };
    const Vec3 to_a = toward(p.a), to_b = toward(p.b);
    const double ra = to_a.norm(), rb = to_b.norm();
    const std::size_t chain = p.chain == 0 ? 0 : static_cast<std::size_t>(p.chain - 1);
    const double total = segments[chain];
    const double near = p.before + ra, far = std::max(total - p.before - (p.b - p.a).norm(), 0.0) + rb;
    double old_product = 1.0, old_sum = 0.0, new_product = 1.0, new_sum = 0.0;
    for (std::size_t k = 0; k < segments.size(); k++)
    {
        old_product *= segments[k];
        old_sum += segments[k];
        if (k == chain)
        {
            new_product *= near * far;
            new_sum += near + far;
        }
        else
        {
            new_product *= segments[k];
            new_sum += segments[k];
        }
    }
    MadeEdge s;
    s.sin_beta = std::max(normal_part(edge, to_b).norm() / rb, 1e-6);
    s.angles = wedge_angles(edge.n, edge.angle(normal_part(edge, to_b)), edge.angle(normal_part(edge, to_a)));
    if (options.fade_incident)
        s.angles.cot[2] = s.angles.cot[3] = 0.0;
    s.distance = near * far * s.sin_beta / (near + far);
    s.ratio = std::sqrt(old_product * old_sum / (new_product * new_sum));
    s.detour = p.detour;
    s.reach = options.fade_reach;
    return s;
}

// Pushes `ray` and the partners its pending screens make.  A building beside
// a leg leaves the ray as it was and adds a partner at each of its edges; a
// building across one stops the ray, leaving the partners.  Over several
// buildings the factors multiply: each term picks the ray or one edge per
// building.  Every building across a leg is kept, or a term would pass
// through it; beside ones fill the rest of three by least detour, and a ray
// across more than three is no ray.
void push_faded(std::vector<PointRay>& rays, PointRay ray, const std::vector<PendingEdge>& pending,
        const std::vector<double>& segments, const Buildings& buildings, const PointPathOptions& options)
{
    if (pending.empty())
    {
        rays.push_back(std::move(ray));
        return;
    }
    // The screens by building, the buildings by their least detour.
    struct Group
    {
        bool across = false;
        double least = 0.0;
        std::vector<std::size_t> members;
    };
    std::vector<Group> groups;
    for (std::size_t k = 0; k < pending.size(); k++)
    {
        if (k == 0 || pending[k].building != pending[k - 1].building || pending[k].chain != pending[k - 1].chain
                || pending[k].a.x != pending[k - 1].a.x || pending[k].a.y != pending[k - 1].a.y)
            groups.push_back({pending[k].across, pending[k].detour, {}});
        groups.back().members.push_back(k);
        groups.back().least = std::min(groups.back().least, pending[k].detour);
    }
    // Across first, then by least detour.
    std::stable_sort(groups.begin(), groups.end(),
            [](const Group& x, const Group& y) { return x.across != y.across ? x.across : x.least < y.least; });
    const std::size_t across = static_cast<std::size_t>(
            std::count_if(groups.begin(), groups.end(), [](const Group& g) { return g.across; }));
    if (across > 3)
        return;
    if (groups.size() > 3)
        groups.resize(3);
    std::vector<MadeEdge> screens(pending.size());
    for (const Group& g : groups)
        for (std::size_t k : g.members)
            screens[k] = complete(buildings, pending[k], segments, options);
    // Each group's choice: -1 the ray unscreened (beside only), else a member.
    std::vector<int> choice(groups.size());
    for (std::size_t i = 0; i < groups.size(); i++)
        choice[i] = groups[i].across ? 0 : -1;
    for (;;)
    {
        PointRay term = ray;
        for (std::size_t i = 0; i < groups.size(); i++)
            if (choice[i] >= 0)
            {
                const std::size_t k = groups[i].members[static_cast<std::size_t>(choice[i])];
                term.made_edges.push_back(screens[k]);
                term.length += pending[k].detour;
                // Its point on the edge, between the leg's ends, for drawing.
                const PendingEdge& p = pending[k];
                for (std::size_t v = 0; v + 1 < term.path.size(); v++)
                    if ((term.path[v] - p.a).norm() + (term.path[v + 1] - p.b).norm() < 1e-9
                            || (term.path[v] - p.b).norm() + (term.path[v + 1] - p.a).norm() < 1e-9)
                    {
                        term.path.insert(term.path.begin() + static_cast<std::ptrdiff_t>(v + 1), p.point);
                        break;
                    }
            }
        rays.push_back(std::move(term));
        std::size_t i = 0;
        for (; i < groups.size(); i++)
        {
            if (++choice[i] < static_cast<int>(groups[i].members.size()))
                break;
            choice[i] = groups[i].across ? 0 : -1;
        }
        if (i == groups.size())
            break;
    }
}

// Whether two sets of options trace the same image trees: every option
// PathEnd reads.
bool same_tracing(const PointPathOptions& a, const PointPathOptions& b)
{
    return a.max_reflections == b.max_reflections && a.max_edge_reflections == b.max_edge_reflections
            && a.max_twice_reflections == b.max_twice_reflections && a.diffraction == b.diffraction
            && a.triple_diffraction == b.triple_diffraction && a.max_thrice_reflections == b.max_thrice_reflections && a.double_diffraction == b.double_diffraction
            && a.source_reflections == b.source_reflections && a.radius == b.radius && a.max_images == b.max_images
            && a.prune == b.prune && (a.reliefs != nullptr) == (b.reliefs != nullptr)
            && a.grating_reflections == b.grating_reflections && a.roofs == b.roofs
            && a.slope_facades == b.slope_facades && a.twice_roofs == b.twice_roofs;
}

} // namespace

PathEnd::PathEnd(const Buildings& buildings, const Vec3& point, bool source, const PointPathOptions& options,
        const std::vector<std::size_t>& partners)
    : point_(point)
    , source_(source)
    , options_(options)
    , near_(edges_near(buildings, point, options.radius))
    , candidates_(joined(near_, ascending(partners)))
    , tree_(image_tree(buildings, point, options.radius, tree_order(buildings, source, options), options.max_images, true,
              candidates_, edge_order(options)))
{
    // A wall, a level roof or a slope the point lies on gives every ray its
    // reflection at the point; an edge would need what the model does not
    // hold.
    const Buildings::Support on = buildings.support(point);
    if (on.kind == Buildings::Support::Kind::edge)
        throw std::invalid_argument("a point lies on a building's edge, or where two of its surfaces meet");
    if (on.kind == Buildings::Support::Kind::facade)
        surface_ = on.index;
    else if (on.kind == Buildings::Support::Kind::roof || on.kind == Buildings::Support::Kind::slope)
        surface_ = buildings.roof_facade(on.index);
    if (options.roofs && options.diffraction && options.max_edge_reflections > 0
            && (!source || options.source_reflections))
    {
        std::vector<std::uint32_t> near;
        buildings.near_box(point.x - options.radius, point.x + options.radius, point.y - options.radius,
                point.y + options.radius, near);
        std::sort(near.begin(), near.end());
        for (const std::uint32_t b : near)
        {
            const Footprint& f = buildings.footprints()[b];
            for (std::size_t i = 0; i < f.roof.size(); i++)
                for (int floor = 0; floor < (options.max_edge_reflections > 1 ? 2 : 1); floor++)
                {
                    const RoofPlane& plane = f.roof[i];
                    const Vec3 from = floor ? Vec3(point.x, point.y, -point.z) : point;
                    const double d = plane.above(from);
                    if (plane.level() || !(d > Buildings::contact))
                        continue;
                    slopes_.push_back({b, i, floor, from - (2.0 * d) * plane.normal()});
                }
        }
    }
    if (options.roofs && options.diffraction && options.double_diffraction && options.twice_roofs
            && options.max_twice_reflections > 0 && (!source || options.source_reflections))
    {
        const auto& edges = buildings.edges();
        std::vector<std::uint32_t> near;
        buildings.near_box(point.x - options.radius, point.x + options.radius, point.y - options.radius,
                point.y + options.radius, near);
        std::sort(near.begin(), near.end());
        for (const std::uint32_t b : near)
        {
            const Footprint& f = buildings.footprints()[b];
            // A level roof is a flat one or a deck at the top (Buildings::on_roof).
            const bool deck = f.roof.empty()
                    || std::any_of(f.roof.begin(), f.roof.end(), [](const RoofPlane& p) { return p.level(); });
            const double* box = buildings.bounds(b);
            for (std::size_t i = 0; i <= f.roof.size(); i++)
            {
                const bool level = i == f.roof.size();
                if (level ? !deck : f.roof[i].level())
                    continue;
                // How far a point stands above the plane, along its normal.
                const auto above = [&](const Vec3& x) { return level ? x.z - f.top : f.roof[i].above(x); };
                for (int floor = 0; floor < (level || options.max_twice_reflections < 2 ? 1 : 2); floor++)
                {
                    const Vec3 from = floor ? below(point) : point;
                    const double d = above(from);
                    if (!(d > Buildings::contact))
                        continue;
                    Bounce bounce;
                    bounce.building = b;
                    bounce.plane = level ? Bounce::level : i;
                    bounce.floor = floor;
                    bounce.image = from - (2.0 * d) * (level ? Vec3(0.0, 0.0, 1.0) : f.roof[i].normal());
                    for (const std::size_t e : candidates_)
                    {
                        // The stretch of the edge above the plane, and where
                        // the lines from its ends to the image cross it.
                        const Edge& g = edges[e];
                        Vec3 lo = g.start, hi = g.start + g.length * g.direction;
                        double d_lo = above(lo), d_hi = above(hi);
                        if (!(std::max(d_lo, d_hi) > Buildings::contact))
                            continue;
                        if (d_lo < 0.0)
                        {
                            lo = lo + (d_lo / (d_lo - d_hi)) * (hi - lo);
                            d_lo = 0.0;
                        }
                        else if (d_hi < 0.0)
                        {
                            hi = hi + (d_hi / (d_hi - d_lo)) * (lo - hi);
                            d_hi = 0.0;
                        }
                        const Vec3 u = lo + (d_lo / (d_lo + d)) * (bounce.image - lo);
                        const Vec3 v = hi + (d_hi / (d_hi + d)) * (bounce.image - hi);
                        if (std::max(u.x, v.x) < box[0] || std::min(u.x, v.x) > box[1] || std::max(u.y, v.y) < box[2]
                                || std::min(u.y, v.y) > box[3])
                            continue;
                        const View view = view_from(g, bounce.image, true, 0.0, g.length);
                        if (!view.lit)
                            continue;
                        bounce.edges.push_back(e);
                        bounce.views.push_back(view);
                    }
                    if (!bounce.edges.empty())
                        bounces_.push_back(std::move(bounce));
                }
            }
        }
    }
    // Each image's places about its edges, and the stretch of each edge an
    // image to one reflection can see, once for every partner: among
    // buildings, most of what the diffracted searches would otherwise solve
    // and then refuse is a point its end cannot see.
    const bool prove = options.prune && options.diffraction;
    views_.resize(tree_.edges.size());
    for (std::size_t i = 0; i < tree_.edges.size(); i++)
    {
        const PointImage& image = tree_.images[i];
        views_[i].resize(2 * tree_.edges[i].size());
        for (std::size_t k = 0; k < tree_.edges[i].size(); k++)
        {
            const Edge& edge = buildings.edges()[tree_.edges[i][k]];
            for (int floor = 0; floor < 2; floor++)
            {
                // A source traced without its reflections is never seen off
                // the floor: nothing reads that view, so nothing proves it.
                double from = 0.0, to = edge.length;
                const bool visible = !prove || (floor == 1 && source && !options.source_reflections)
                        || visible_stretch(buildings, edge, point, image.facades, floor == 1, from, to);
                views_[i][2 * k + floor] = view_from(edge, floor ? below(image.point) : image.point, visible, from, to,
                        on_surface() && image.facades.empty());
            }
        }
    }
}

PointSourceField::PointSourceField(const Buildings& buildings, const Vec3& source, const Vec3& microphone,
        const PointPathOptions& options)
    : facade_count_(buildings.facades().size())
    , buildings_(&buildings)
    , edge_source_(options.edge_source)
    , sample_shared_faces_(options.sample_shared_faces)
    , reliefs_(options.reliefs)
    , edge_source_threads_(options.edge_source_threads)
{
    const PathEnd source_end(buildings, source, true, options, edges_near(buildings, microphone, options.radius));
    const PathEnd microphone_end(buildings, microphone, false, options, edges_near(buildings, source, options.radius));
    trace(buildings, source_end, microphone_end, options);
    for (PointRay& ray : rays_)
    {
        ray.source_surface = source_end.surface();
        ray.microphone_surface = microphone_end.surface();
    }
}

PointSourceField::PointSourceField(const Buildings& buildings, const PathEnd& source, const PathEnd& microphone,
        const PointPathOptions& options)
    : facade_count_(buildings.facades().size())
    , buildings_(&buildings)
    , edge_source_(options.edge_source)
    , sample_shared_faces_(options.sample_shared_faces)
    , reliefs_(options.reliefs)
    , edge_source_threads_(options.edge_source_threads)
{
    trace(buildings, source, microphone, options);
    for (PointRay& ray : rays_)
    {
        ray.source_surface = source.surface();
        ray.microphone_surface = microphone.surface();
    }
}

void PointSourceField::trace(const Buildings& buildings, const PathEnd& source_end, const PathEnd& microphone_end,
        const PointPathOptions& options)
{
    // The edges near this pair; each end was traced with them all, and each
    // image keeps those that cross its beam, in the order the pair's own tree
    // would list them.
    if (!source_end.source() || microphone_end.source() || !same_tracing(source_end.options(), options)
            || !same_tracing(microphone_end.options(), options))
        throw std::invalid_argument("PointSourceField: the ends must be a source and a microphone traced with the "
                                    "field's options");
    const std::vector<std::size_t> near_edges = joined(source_end.near(), microphone_end.near());
    for (const PathEnd* end : {&source_end, &microphone_end})
        if (!std::includes(end->candidates().begin(), end->candidates().end(), near_edges.begin(), near_edges.end()))
            throw std::invalid_argument("a PathEnd was traced without the edges near the end it is paired with");
    const auto& edges = buildings.edges();
    thread_local std::vector<char> is_near;
    is_near.resize(edges.size(), 0);
    for (std::size_t e : near_edges)
        is_near[e] = 1;
    auto keep = [&](const PathEnd& end, std::vector<std::vector<std::size_t>>& kept,
                        std::vector<std::vector<std::uint32_t>>& at) {
        const ImageTree& tree = end.tree();
        kept.assign(tree.edges.size(), {});
        at.assign(tree.edges.size(), {});
        for (std::size_t i = 0; i < tree.edges.size(); i++)
            for (std::size_t k = 0; k < tree.edges[i].size(); k++)
                if (is_near[tree.edges[i][k]])
                {
                    kept[i].push_back(tree.edges[i][k]);
                    at[i].push_back(static_cast<std::uint32_t>(k));
                }
    };
    std::vector<std::vector<std::size_t>> src_edges, mic_edges;
    std::vector<std::vector<std::uint32_t>> src_at, mic_at;
    keep(source_end, src_edges, src_at);
    keep(microphone_end, mic_edges, mic_at);
    // Each end's roof bounces' edges near this pair (PathEnd::bounces), by
    // their places in the bounce's list.
    auto keep_bounces = [&](const PathEnd& end, std::vector<std::vector<std::uint32_t>>& at) {
        at.assign(end.bounces().size(), {});
        for (std::size_t i = 0; i < end.bounces().size(); i++)
            for (std::size_t k = 0; k < end.bounces()[i].edges.size(); k++)
                if (is_near[end.bounces()[i].edges[k]])
                    at[i].push_back(static_cast<std::uint32_t>(k));
    };
    std::vector<std::vector<std::uint32_t>> src_bounce_at, mic_bounce_at;
    keep_bounces(source_end, src_bounce_at);
    keep_bounces(microphone_end, mic_bounce_at);
    for (std::size_t e : near_edges)
        is_near[e] = 0;
    // Specular tracing may need source images even when diffraction may not
    // use them. Keep that option's meaning for every diffracted search.
    std::vector<PointImage> direct_source;
    if (!options.source_reflections)
        direct_source.push_back(source_end.tree().images.front());
    const Ends src {options.source_reflections ? source_end.tree().images : direct_source, src_edges, src_at, source_end.views()},
            mic {microphone_end.tree().images, mic_edges, mic_at, microphone_end.views()};
    const Vec3& source = source_end.point();
    const Vec3& microphone = microphone_end.point();
    const int source_floors = options.source_reflections ? 2 : 1;
    // An end on a wall, a level roof or a slope is reached along the face of an edge
    // of that surface by a leg straight from the wedge, with no facade
    // between: end_face, that face (edge_face) or -1, its floor image's off
    // the floor (on a wall's face too).  Its leg is then on the face, every
    // other end's strictly in the air (end_air).
    const auto end_face = [](const PathEnd& end, const Edge& e, const std::vector<std::size_t>& facades, int floor) {
        return end.on_surface() && facades.empty() ? edge_face(e, floor ? below(end.point()) : end.point()) : -1;
    };
    const auto end_air = [](const Edge& e, const Vec3& v, int face) { return face >= 0 || air(e, v, false); };
    // A leg along the face an end lies on whose wedge's other leg runs along
    // that face too doubles back on it: the edge's grazing backscatter, which
    // the end's own rays over the face already carry, and no ray.
    const auto doubles_back = [](const Edge& e, int face, const Vec3& other) {
        return face >= 0 && edge_face(e, other) == face;
    };
    const bool fading = options.fade_reach > 0.0;
    // Roofs (options.roofs): the buildings whose boxes meet a box in plan,
    // ascending.
    thread_local std::vector<std::uint32_t> boxed;
    // The buildings whose roofs a folded path may bounce on: its real legs
    // each run between two consecutive members of the chain -- its two ends
    // and the facades it reflects in -- so in plan the path stays inside the
    // box of them all, wherever along a facade or an edge it meets them.  The
    // unfolded line's box can miss a roof under a reflected leg, which the
    // plan grid's coarseness then decided.
    auto buildings_along = [&](std::initializer_list<Vec3> points,
                                   const std::vector<std::size_t>& facades) -> const std::vector<std::uint32_t>& {
        double x0 = std::numeric_limits<double>::infinity(), x1 = -x0, y0 = x0, y1 = -x0;
        for (const Vec3& p : points)
        {
            x0 = std::min(x0, p.x);
            x1 = std::max(x1, p.x);
            y0 = std::min(y0, p.y);
            y1 = std::max(y1, p.y);
        }
        for (std::size_t f : facades)
        {
            const Facade& wall = buildings.facades()[f];
            x0 = std::min({x0, wall.x0, wall.x1});
            x1 = std::max({x1, wall.x0, wall.x1});
            y0 = std::min({y0, wall.y0, wall.y1});
            y1 = std::max({y1, wall.y0, wall.y1});
        }
        boxed.clear();
        buildings.near_box(x0, x1, y0, y1, boxed);
        std::sort(boxed.begin(), boxed.end());
        return boxed;
    };

    // Specular rays.
    for (const PointImage& image : mic.images)
    {
        if (image.facades.size() > options.max_reflections)
            continue;
        for (int floor = 0; floor < 2; floor++)
        {
            std::vector<Vec3> hits;
            std::vector<PendingEdge> pending;
            const std::size_t n = image.facades.size();
            double before = 0.0;
            // Leg m has m reflections before it and n - m after: its edges'
            // diffracted rays are traced with at most max_edge_reflections
            // each side.
            auto leg = [&](std::size_t m, const Vec3& p, const Vec3& q) {
                const bool traced = options.diffraction && m <= options.max_edge_reflections
                        && n - m <= options.max_edge_reflections && (m == 0 || options.source_reflections);
                const double start = before;
                before += (q - p).norm();
                if (traced && !options.fade_traced)
                    return !buildings.blocked(p, q, floor == 1);
                return leg_edges(buildings, p, q, floor == 1, options.fade_reach,
                        m == 0 ? no_building : buildings.facade_building(image.facades[n - m]),
                        m == n ? no_building : buildings.facade_building(image.facades[n - 1 - m]), start, 0,
                        pending);
            };
            if (fading ? !folded_path_open(buildings, source, microphone, image.facades, floor == 1, &hits, leg)
                       : !folded_path_open(buildings, source, microphone, image.facades, floor == 1, &hits))
                continue;
            PointRay ray;
            ray.path.push_back(source);
            ray.path.insert(ray.path.end(), hits.begin(), hits.end());
            ray.path.push_back(microphone);
            ray.length = ((floor ? below(image.point) : image.point) - source).norm();
            ray.spreading = 1.0 / ray.length;
            ray.facades = image.facades;
            if (floor)
            {
                ray.source_plane = 0.0;
                add_floor_bounce(ray, source, below(image.point), true);
            }
            if (options.reliefs != nullptr)
                meeting_toward(buildings, image.facades, (floor ? below(image.point) : image.point) - source, ray.hits);
            const std::vector<double> segments {ray.length};
            push_faded(rays_, std::move(ray), pending, segments, buildings, options);
        }
    }
    // Specular rays off a roof below both ends, in place of the floor.
    if (options.roofs)
        for (const PointImage& image : mic.images)
        {
            if (image.facades.size() + 1 > options.max_reflections)
                continue;
            for (std::uint32_t b : buildings_along({source, microphone}, image.facades))
            {
                // An end on the roof (within contact) takes its reflection at
                // the end itself (PathEnd::surface).
                const double h = buildings.footprints()[b].top;
                if (!(h < std::min(source.z, microphone.z) - Buildings::contact))
                    continue;
                std::vector<Vec3> hits;
                if (!folded_roof_path_open(buildings, source, microphone, image.facades, b, &hits))
                    continue;
                PointRay ray;
                ray.path.push_back(source);
                ray.path.insert(ray.path.end(), hits.begin(), hits.end());
                ray.path.push_back(microphone);
                ray.length = (Vec3(image.point.x, image.point.y, 2.0 * h - image.point.z) - source).norm();
                ray.spreading = 1.0 / ray.length;
                ray.facades = image.facades;
                ray.roof_facades = {buildings.roof_facade(b)};
                ray.source_plane = h;
                if (options.reliefs != nullptr)
                    meeting_toward(buildings, image.facades,
                            Vec3(image.point.x, image.point.y, 2.0 * h - image.point.z) - source, ray.hits);
                rays_.push_back(std::move(ray));
            }
        }
    // Specular rays off a pitched roof's slope: the source mirrored in the
    // slope's plane, the line from there to the microphone crossing it where
    // that plane is the roof, both ends on its air side, each with the floor
    // or without; no facade, since a slope turns a ray up or down and the
    // walls' images, vertical, then no longer hold.  An end on the slope
    // takes its reflection at the end itself (PathEnd::surface).
    if (options.roofs && options.max_reflections > 0)
    {
        thread_local std::vector<std::uint32_t> near;
        buildings.near_box(std::min(source.x, microphone.x) - options.radius,
                std::max(source.x, microphone.x) + options.radius, std::min(source.y, microphone.y) - options.radius,
                std::max(source.y, microphone.y) + options.radius, near);
        std::sort(near.begin(), near.end());
        for (const std::uint32_t b : near)
        {
            const Footprint& f = buildings.footprints()[b];
            const double* box = buildings.bounds(b);
            const auto off = [&](const Vec3& p) {
                return std::hypot(std::max({box[0] - p.x, p.x - box[1], 0.0}), std::max({box[2] - p.y, p.y - box[3], 0.0}));
            };
            if (f.roof.empty() || std::min(off(source), off(microphone)) > options.radius)
                continue;
            for (std::size_t i = 0; i < f.roof.size(); i++)
            {
                const RoofPlane& plane = f.roof[i];
                if (plane.level())
                    continue;
                const Vec3 n = plane.normal();
                for (int fa = 0; fa < 2; fa++)
                    for (int fb = 0; fb < 2; fb++)
                    {
                        const Vec3 pa = fa ? below(source) : source, pb = fb ? below(microphone) : microphone;
                        const double da = plane.above(pa), db = plane.above(pb);
                        if (!(da > Buildings::contact && db > Buildings::contact))
                            continue;
                        const Vec3 image = pa - (2.0 * da) * n;
                        const Vec3 p = image + (da / (da + db)) * (pb - image);
                        if (!buildings.on_plane(b, i, p.x, p.y))
                            continue;
                        // The floor's bounce on either side, where the
                        // unfolded leg crosses it; an end on the floor is its
                        // own bounce.
                        const Vec3 qa = fa ? pa + (-pa.z / (p.z - pa.z)) * (p - pa) : source;
                        const Vec3 qb = fb ? p + (p.z / (p.z - pb.z)) * (pb - p) : microphone;
                        if ((fa && (buildings.inside(qa.x, qa.y) < buildings.size() || buildings.blocked(source, qa)))
                                || (fb
                                        && (buildings.inside(qb.x, qb.y) < buildings.size()
                                                || buildings.blocked(qb, microphone)))
                                || buildings.blocked(qa, p) || buildings.blocked(p, qb))
                            continue;
                        PointRay ray;
                        ray.path.push_back(source);
                        if (fa && (qa - source).norm() > Buildings::contact)
                            ray.path.push_back(qa);
                        else if (fa)
                            ray.launch_mirrored = true;
                        ray.path.push_back(p);
                        if (fb && (qb - microphone).norm() > Buildings::contact)
                            ray.path.push_back(qb);
                        ray.path.push_back(microphone);
                        ray.length = (pb - image).norm();
                        ray.spreading = 1.0 / ray.length;
                        ray.roof_facades = {buildings.roof_facade(b)};
                        // The floor's spherical-wave coefficient takes the
                        // whole unfolded path as its image's: the slope is a
                        // mirror, not a new source.
                        if (fa)
                        {
                            ray.source_plane = 0.0;
                            add_floor_bounce(ray, std::abs(p.z - pa.z) / (p - pa).norm(), ray.length, false);
                        }
                        if (fb)
                            add_floor_bounce(ray, std::abs(pb.z - p.z) / (pb - p).norm(), ray.length, true);
                        rays_.push_back(std::move(ray));
                    }
            }
        }
    }

    // A slope between two chains of facades: the source's image through its
    // chain (src) and the microphone's through its (mic), each with the
    // floor or without, mirrored in the slope's plane between them -- a
    // tilted mirror does not commute with the walls' vertical ones, but with
    // the slope between the chains it need not -- each chain traced back
    // from the reflection, with the slope to max_reflections in all; the
    // chains both empty are the loop above's. Edge candidates are not a roof
    // visibility proof: a beam can reach a roof whose edges were out of range.
    if (options.roofs && options.slope_facades && options.max_reflections > 1)
    {
        const auto& source_images = source_end.tree().images;
        thread_local std::vector<std::uint32_t> near;
        buildings.near_box(std::min(source.x, microphone.x) - options.radius,
                std::max(source.x, microphone.x) + options.radius, std::min(source.y, microphone.y) - options.radius,
                std::max(source.y, microphone.y) + options.radius, near);
        std::sort(near.begin(), near.end());
        struct Candidate
        {
            std::size_t image;
            int floor;
            Vec3 point;
            double above;
        };
        std::vector<Candidate> source_side, microphone_side;
        std::vector<std::size_t> reversed_chain;
        std::vector<Vec3> from_source, to_microphone;
        for (const std::uint32_t b : near)
        {
            const Footprint& f = buildings.footprints()[b];
            if (f.roof.empty())
                continue;
            for (std::size_t i = 0; i < f.roof.size(); i++)
            {
                const RoofPlane& plane = f.roof[i];
                if (plane.level())
                    continue;
                const Vec3 n = plane.normal();
                // Each side's images on the plane's air side.
                const auto gather = [&](const std::vector<PointImage>& images, std::vector<Candidate>& out) {
                    out.clear();
                    for (std::size_t k = 0; k < images.size(); k++)
                    {
                        if (images[k].facades.size() + 1 > options.max_reflections)
                            continue;
                        for (int fl = 0; fl < 2; fl++)
                        {
                            const Vec3 pt = fl ? below(images[k].point) : images[k].point;
                            const double d = plane.above(pt);
                            if (d > Buildings::contact)
                                out.push_back({k, fl, pt, d});
                        }
                    }
                };
                gather(source_images, source_side);
                gather(mic.images, microphone_side);
                for (const Candidate& a : source_side)
                    for (const Candidate& c : microphone_side)
                    {
                        const std::vector<std::size_t>& chain_a = source_images[a.image].facades;
                        const std::vector<std::size_t>& chain_b = mic.images[c.image].facades;
                        if ((chain_a.empty() && chain_b.empty())
                                || chain_a.size() + chain_b.size() + 1 > options.max_reflections)
                            continue;
                        const Vec3 image = a.point - (2.0 * a.above) * n;
                        const Vec3 p = image + (a.above / (a.above + c.above)) * (c.point - image);
                        if (!buildings.on_plane(b, i, p.x, p.y))
                            continue;
                        // The source's side traced from the source, so that
                        // its first leg is its launch. Both chains are
                        // folded into the real half-space at the slope.
                        reversed_chain.assign(chain_a.rbegin(), chain_a.rend());
                        if (!physical_chain(buildings, source, p, reversed_chain, a.floor == 1, from_source)
                                || !physical_chain(buildings, p, microphone, chain_b, c.floor == 1, to_microphone))
                            continue;
                        PointRay ray;
                        ray.path.push_back(source);
                        ray.launch_mirrored = a.floor && source.z <= Buildings::contact;
                        ray.path.insert(ray.path.end(), from_source.begin(), from_source.end());
                        ray.path.push_back(p);
                        ray.path.insert(ray.path.end(), to_microphone.begin(), to_microphone.end());
                        ray.path.push_back(microphone);
                        ray.length = (c.point - image).norm();
                        ray.spreading = 1.0 / ray.length;
                        ray.facades = chain_a;
                        ray.facades.insert(ray.facades.end(), chain_b.begin(), chain_b.end());
                        ray.roof_facades = {buildings.roof_facade(b)};
                        if (a.floor)
                            add_floor_bounce(ray, std::abs(p.z - a.point.z) / (p - a.point).norm(), ray.length, false);
                        if (c.floor)
                            add_floor_bounce(ray, std::abs(c.point.z - p.z) / (c.point - p).norm(), ray.length, true);
                        if (options.reliefs != nullptr)
                        {
                            meeting_toward(buildings, chain_a, a.point - p, ray.hits);
                            meeting_toward(buildings, chain_b, c.point - p, ray.hits);
                        }
                        rays_.push_back(std::move(ray));
                    }
            }
        }
    }

    if (!options.diffraction)
        return;

    // Diffracted once, between an image of the source and one of the
    // microphone.
    for (std::size_t i = 0; i < src.images.size(); i++)
        for (std::size_t j = 0; j < mic.images.size(); j++)
        {
            const PointImage& a = src.images[i];
            const PointImage& b = mic.images[j];
            if (a.facades.size() > options.max_edge_reflections || b.facades.size() > options.max_edge_reflections)
                continue;
            // The edges both images keep: both lists run in near_edges' order,
            // so one pass along the two finds them.
            std::size_t next = 0;
            const auto& keep = mic.edges[j];
            for (std::size_t ki = 0; ki < src.edges[i].size(); ki++)
            {
                const std::size_t e = src.edges[i][ki];
                while (next < keep.size() && keep[next] < e)
                    next++;
                if (next == keep.size() || keep[next] != e)
                    continue;
                const Edge& edge = edges[e];
                for (int fa = 0; fa < source_floors; fa++)
                    for (int fb = 0; fb < 2; fb++)
                    {
                        // A point of the edge either end cannot see
                        // diffracts nothing between them.
                        const PathEnd::View& va = src.view(i, ki, fa);
                        const PathEnd::View& vb = mic.view(j, next, fb);
                        double t;
                        if (!keller_point(va.along, va.off, vb.along, vb.off, t))
                            continue;
                        // A Keller point past an end still has the end's wave
                        // (ended_wedge_coefficient, or the edge-source integral): the
                        // end stands in for it in every check, within the reach.
                        const bool past = (options.edge_source || options.edge_ends) && !fading
                                && past_true_end(edge, t);
                        const double tq = past ? std::clamp(t, 0.0, edge.length) : t;
                        // Fading keeps points a hair hidden, which the
                        // visibility proof would drop.
                        if (fading ? t < 0.0 || t > edge.length
                                   : !va.visible || !vb.visible || tq < std::max(va.from, vb.from)
                                        || tq > std::min(va.to, vb.to))
                            continue;
                        const Vec3 pa = fa ? below(a.point) : a.point, pb = fb ? below(b.point) : b.point;
                        const Vec3 q = edge.start + tq * edge.direction;
                        if (q.z < 0.0 || buildings.buried(q))
                            continue;
                        const Vec3 keller = edge.start + t * edge.direction;
                        if (past
                                && !end_within_reach(pa, pb, q, keller,
                                        options.edge_source ? options.edge_source_reach : options.end_reach))
                            continue;
                        const Vec3 in = pa - q, out = pb - q;
                        const int on_a = end_face(source_end, edge, a.facades, fa);
                        const int on_b = end_face(microphone_end, edge, b.facades, fb);
                        // Past an end, UTD's end wave goes wherever the end sees
                        // (onto_air); the edge-source integral keeps the wedge's.
                        const bool vertex = past && !options.edge_source;
                        if ((!vertex
                                    && (!end_air(edge, normal_part(edge, in), on_a)
                                            || !end_air(edge, normal_part(edge, out), on_b)))
                                || doubles_back(edge, on_a, pb) || doubles_back(edge, on_b, pa))
                            continue;
                        std::vector<Vec3> to_source, to_microphone;
                        std::vector<PendingEdge> pending;
                        // Leg m from the edge has m reflections between them:
                        // a second edge there is traced with
                        // max_between_reflections at most between the two
                        // and none of them the floor, and
                        // max_twice_reflections at either end.
                        double before = 0.0;
                        auto legs = [&](const std::vector<std::size_t>& facades, bool floor, int chain) {
                            const std::size_t count = facades.size();
                            before = 0.0;
                            return [&, count, floor, chain](std::size_t m, const Vec3& p, const Vec3& r) {
                                const std::size_t other = chain == 1 ? b.facades.size() : a.facades.size();
                                const bool traced = options.double_diffraction && m <= options.max_between_reflections
                                        && (m == 0 || !floor) && count - m <= options.max_twice_reflections
                                        && other <= options.max_twice_reflections;
                                const double start = before;
                                before += (r - p).norm();
                                if (traced && !options.fade_traced)
                                    return !buildings.blocked(p, r, floor);
                                return leg_edges(buildings, p, r, floor, options.fade_reach,
                                        m == 0 ? edge.building : buildings.facade_building(facades[count - m]),
                                        m == count ? no_building : buildings.facade_building(facades[count - 1 - m]),
                                        start, chain, pending);
                            };
                        };
                        if (fading)
                        {
                            if (!folded_path_open(buildings, q, source, a.facades, fa == 1, &to_source,
                                        legs(a.facades, fa == 1, 1))
                                    || !folded_path_open(buildings, q, microphone, b.facades, fb == 1, &to_microphone,
                                            legs(b.facades, fb == 1, 2)))
                                continue;
                        }
                        else if (!folded_path_open(buildings, q, source, a.facades, fa == 1, &to_source)
                                || !folded_path_open(buildings, q, microphone, b.facades, fb == 1, &to_microphone))
                            continue;
                        PointRay ray;
                        const double raw0 = end_angle(edge, in, on_a, ray.source_on_face);
                        const double raw = end_angle(edge, out, on_b, ray.microphone_on_face);
                        const double theta0 = onto_air(raw0, edge.n), theta = onto_air(raw, edge.n);
                        const double s_in = in.norm(), s_out = out.norm();
                        // The coefficient's sin(beta) and distance are the Keller point's,
                        // the same either way round (edge_ends has why).
                        const double s_in_k = (pa - keller).norm(), s_out_k = (pb - keller).norm();
                        const double sin_beta = normal_part(edge, in).norm() / s_in_k;
                        ray.length = s_in + s_out;
                        ray.spreading = 1.0 / std::sqrt(s_in * s_out * (s_in + s_out));
                        ray.facades = a.facades;
                        ray.facades.insert(ray.facades.end(), b.facades.begin(), b.facades.end());
                        if (fa)
                        {
                            ray.source_plane = 0.0;
                            add_floor_bounce(ray, pa, q, false);
                        }
                        if (fb)
                            add_floor_bounce(ray, q, pb, true);
                        ray.wedges = 1;
                        ray.wedge[0] = {wedge_angles(edge.n, theta, theta0), sin_beta,
                                s_in_k * s_out_k * sin_beta / (s_in_k + s_out_k), false};
                        if (options.edge_ends && !fading)
                        {
                            ray.wedge[0].ends = edge_ends(edge, pa, pb, t);
                            weigh_vertex(ray.wedge[0].ends, raw0, raw, edge.n, s_in, s_out);
                        }
                        ray.edge_source.edges[0] = e;
                        ray.edge_source.source = pa;
                        ray.edge_source.microphone = pb;
                        ray.edge_source.stretches[0] = {std::max(va.from, vb.from), std::min(va.to, vb.to)};
                        ray.path = {source};
                        ray.path.insert(ray.path.end(), to_source.rbegin(), to_source.rend());
                        ray.path.push_back(q);
                        ray.path.insert(ray.path.end(), to_microphone.begin(), to_microphone.end());
                        ray.path.push_back(microphone);
                        push_faded(rays_, std::move(ray), pending, {s_in, s_out}, buildings, options);
                    }
            }
        }
    // Diffracted once with a roof on one side or both, in place of the floor
    // there: a roof below that end and below the edge, between the image and
    // the edge in plan, counted as one of the side's reflections.  The roof
    // bounds which stretch of the edge its image sees, which the visibility
    // proof does not know, so these are checked whole; and they take UTD.
    if (options.roofs)
    {
        // One side's horizontal reflector: the floor (-1), a roof (its
        // building), or nothing (-2).
        constexpr long floor_side = -1, bare_side = -2;
        auto mirrored = [&](const Vec3& p, long side) {
            if (side == bare_side)
                return p;
            const double h = side == floor_side ? 0.0 : buildings.footprints()[static_cast<std::size_t>(side)].top;
            return Vec3(p.x, p.y, 2.0 * h - p.z);
        };
        std::vector<long> source_sides, microphone_sides;
        for (std::size_t i = 0; i < src.images.size(); i++)
            for (std::size_t j = 0; j < mic.images.size(); j++)
            {
                const PointImage& a = src.images[i];
                const PointImage& b = mic.images[j];
                const bool roof_a = options.source_reflections && a.facades.size() + 1 <= options.max_edge_reflections;
                const bool roof_b = b.facades.size() + 1 <= options.max_edge_reflections;
                if ((!roof_a && !roof_b) || a.facades.size() > options.max_edge_reflections
                        || b.facades.size() > options.max_edge_reflections)
                    continue;
                std::size_t next = 0;
                const auto& keep = mic.edges[j];
                for (std::size_t ki = 0; ki < src.edges[i].size(); ki++)
                {
                    const std::size_t e = src.edges[i][ki];
                    while (next < keep.size() && keep[next] < e)
                        next++;
                    if (next == keep.size() || keep[next] != e)
                        continue;
                    const Edge& edge = edges[e];
                    const Vec3 end = edge.start + edge.length * edge.direction;
                    const double edge_top = std::max(edge.start.z, end.z);
                    // Each side's reflectors at this edge: the floor and
                    // nothing as the loop above takes them, and the roofs
                    // under its real legs from the edge, by its facades, to
                    // its end.
                    auto sides = [&](const Vec3& real_end, const std::vector<std::size_t>& facades, double end_z,
                                         bool roofs_too, int floors, std::vector<long>& out) {
                        out.clear();
                        out.push_back(bare_side);
                        if (floors > 1)
                            out.push_back(floor_side);
                        if (!roofs_too)
                            return;
                        // The edge is a member of the chain: both its ends.
                        for (std::uint32_t r : buildings_along({edge.start, end, real_end}, facades))
                        {
                            // An end on the roof (within contact) takes its
                            // reflection at the end itself (PathEnd::surface).
                            const double h = buildings.footprints()[r].top;
                            if (h < end_z - Buildings::contact && h < edge_top - 1e-9)
                                out.push_back(static_cast<long>(r));
                        }
                    };
                    sides(source, a.facades, source.z, roof_a, source_floors, source_sides);
                    sides(microphone, b.facades, microphone.z, roof_b, 2, microphone_sides);
                    for (long sa : source_sides)
                        for (long sb : microphone_sides)
                        {
                            if (sa < 0 && sb < 0)
                                continue; // the loop above's
                            const Vec3 pa = mirrored(a.point, sa), pb = mirrored(b.point, sb);
                            const PathEnd::View va = sa >= 0 ? view_from(edge, pa, true, 0.0, edge.length)
                                                             : src.view(i, ki, sa == floor_side ? 1 : 0);
                            const PathEnd::View vb = sb >= 0 ? view_from(edge, pb, true, 0.0, edge.length)
                                                             : mic.view(j, next, sb == floor_side ? 1 : 0);
                            double t;
                            if (!keller_point(va.along, va.off, vb.along, vb.off, t))
                                continue;
                            // Past an end, the end stands in for the Keller
                            // point in every check (the loop above has why).
                            const bool past = options.edge_ends && past_true_end(edge, t);
                            const double tq = past ? std::clamp(t, 0.0, edge.length) : t;
                            if (!va.visible || !vb.visible || tq < std::max(va.from, vb.from)
                                    || tq > std::min(va.to, vb.to))
                                continue;
                            const Vec3 q = edge.start + tq * edge.direction, keller = edge.start + t * edge.direction;
                            if (q.z < 0.0 || buildings.buried(q)
                                    || (past && !end_within_reach(pa, pb, q, keller, options.end_reach)))
                                continue;
                            if ((sa >= 0 && !(q.z > buildings.footprints()[static_cast<std::size_t>(sa)].top + 1e-9))
                                    || (sb >= 0
                                            && !(q.z > buildings.footprints()[static_cast<std::size_t>(sb)].top + 1e-9)))
                                continue;
                            const Vec3 in = pa - q, out = pb - q;
                            const int on_a = sa < 0 ? end_face(source_end, edge, a.facades, sa == floor_side) : -1;
                            const int on_b = sb < 0 ? end_face(microphone_end, edge, b.facades, sb == floor_side) : -1;
                            // Past an end, the end wave goes wherever the end sees
                            // (onto_air).
                            if ((!past
                                        && (!end_air(edge, normal_part(edge, in), on_a)
                                                || !end_air(edge, normal_part(edge, out), on_b)))
                                    || doubles_back(edge, on_a, pb) || doubles_back(edge, on_b, pa))
                                continue;
                            std::vector<Vec3> to_source, to_microphone;
                            const bool source_open = sa >= 0
                                    ? folded_roof_path_open(buildings, q, source, a.facades,
                                              static_cast<std::size_t>(sa), &to_source)
                                    : folded_path_open(buildings, q, source, a.facades, sa == floor_side, &to_source);
                            if (!source_open)
                                continue;
                            const bool microphone_open = sb >= 0
                                    ? folded_roof_path_open(buildings, q, microphone, b.facades,
                                              static_cast<std::size_t>(sb), &to_microphone)
                                    : folded_path_open(buildings, q, microphone, b.facades, sb == floor_side,
                                              &to_microphone);
                            if (!microphone_open)
                                continue;
                            PointRay ray;
                            const double raw0 = end_angle(edge, in, on_a, ray.source_on_face);
                            const double raw = end_angle(edge, out, on_b, ray.microphone_on_face);
                            const double theta0 = onto_air(raw0, edge.n), theta = onto_air(raw, edge.n);
                            const double s_in = in.norm(), s_out = out.norm();
                            // The coefficient's sin(beta) and distance are the Keller point's,
                            // the same either way round (edge_ends has why).
                            const double s_in_k = (pa - keller).norm(), s_out_k = (pb - keller).norm();
                            const double sin_beta = normal_part(edge, in).norm() / s_in_k;
                            ray.length = s_in + s_out;
                            ray.spreading = 1.0 / std::sqrt(s_in * s_out * (s_in + s_out));
                            ray.facades = a.facades;
                            ray.facades.insert(ray.facades.end(), b.facades.begin(), b.facades.end());
                            if (sa == floor_side)
                                add_floor_bounce(ray, pa, q, false);
                            if (sb == floor_side)
                                add_floor_bounce(ray, q, pb, true);
                            for (long side : {sa, sb})
                                if (side >= 0)
                                    ray.roof_facades.push_back(buildings.roof_facade(static_cast<std::size_t>(side)));
                            if (sa != bare_side)
                                ray.source_plane
                                        = sa == floor_side ? 0.0 : buildings.footprints()[static_cast<std::size_t>(sa)].top;
                            ray.wedges = 1;
                            ray.wedge[0] = {wedge_angles(edge.n, theta, theta0), sin_beta,
                                    s_in_k * s_out_k * sin_beta / (s_in_k + s_out_k), false};
                            if (options.edge_ends)
                            {
                                ray.wedge[0].ends = edge_ends(edge, pa, pb, t);
                                weigh_vertex(ray.wedge[0].ends, raw0, raw, edge.n, s_in, s_out);
                            }
                            ray.edge_source.edges[0] = e;
                            ray.edge_source.source = pa;
                            ray.edge_source.microphone = pb;
                            ray.edge_source.stretches[0] = {std::max(va.from, vb.from), std::min(va.to, vb.to)};
                            ray.path = {source};
                            ray.path.insert(ray.path.end(), to_source.rbegin(), to_source.rend());
                            ray.path.push_back(q);
                            ray.path.insert(ray.path.end(), to_microphone.begin(), to_microphone.end());
                            ray.path.push_back(microphone);
                            rays_.push_back(std::move(ray));
                        }
                }
            }
    }

    // A roof's side of a ray, from a point q off the roof to the end: the
    // reflection where the line from q to the end's image in the roof's
    // plane -- a slope's, or a level roof's at its top (PathEnd::Bounce::
    // level) -- crosses that plane where it is the roof, then the floor before
    // the end or not; false where it misses the roof or a leg is blocked.
    // `points` gets the reflection and the floor's bounce, from q.
    const auto roof_leg = [&](std::size_t building, std::size_t plane, int floor, const Vec3& image, const Vec3& q,
                                  const Vec3& end, std::vector<Vec3>& points) {
        const Footprint& f = buildings.footprints()[building];
        const bool level = plane == PathEnd::Bounce::level;
        const double dq = level ? q.z - f.top : f.roof[plane].above(q);
        const double di = level ? f.top - image.z : -f.roof[plane].above(image);
        if (!(dq > Buildings::contact) || !(di > 0.0))
            return false;
        Vec3 p = q + (dq / (dq + di)) * (image - q);
        if (level)
            p.z = f.top;
        if (level ? !buildings.on_roof(building, p.x, p.y) : !buildings.on_plane(building, plane, p.x, p.y))
            return false;
        const Vec3 target = floor ? below(end) : end;
        const Vec3 bounce = floor ? p + (p.z / (p.z - target.z)) * (target - p) : end;
        if (buildings.blocked(q, p) || buildings.blocked(p, bounce)
                || (floor
                        && (buildings.inside(bounce.x, bounce.y) < buildings.size() || buildings.blocked(bounce, end))))
            return false;
        points = {p};
        if (floor && (bounce - end).norm() > Buildings::contact)
            points.push_back(bounce);
        return true;
    };

    // Diffracted once with a slope's reflection on one side: the partners a
    // slope's specular rays need where an edge cuts off the wave that reaches
    // the slope, or the one that leaves it.  As for those rays, no facade:
    // the end, after the floor or not, mirrored in the slope's plane, and the
    // reflection where the line from the edge's point to that image crosses
    // the plane where that plane is the roof, the edge on the plane's air
    // side and every leg open; the other side bare or off the floor.  One
    // slope a ray, counted with its floor as that side's reflections.  They
    // take UTD, as the roofs' do.
    if (options.roofs && options.max_edge_reflections > 0 && !src.images.empty() && !mic.images.empty())
    {
        using Slope = PathEnd::Slope;
        const std::vector<Slope>& source_slopes = source_end.slopes();
        const std::vector<Slope>& microphone_slopes = microphone_end.slopes();
        auto slope_leg = [&](const Slope& sl, const Vec3& q, const Vec3& end, std::vector<Vec3>& points) {
            return roof_leg(sl.building, sl.plane, sl.floor, sl.image, q, end, points);
        };
        const auto& source_keep = src.edges[0];
        const auto& microphone_keep = mic.edges[0];
        std::size_t next = 0;
        std::vector<Vec3> to_source, to_microphone;
        for (std::size_t ki = 0; ki < source_keep.size() && !(source_slopes.empty() && microphone_slopes.empty());
                ki++)
        {
            const std::size_t e = source_keep[ki];
            while (next < microphone_keep.size() && microphone_keep[next] < e)
                next++;
            if (next == microphone_keep.size() || microphone_keep[next] != e)
                continue;
            const Edge& edge = edges[e];
            // The source's side and the microphone's: a slope (its index into
            // that end's list) or the end bare (-1) or off the floor (-2).
            for (int side = 0; side < 2; side++)
            {
                const auto& slopes = side == 0 ? source_slopes : microphone_slopes;
                for (std::size_t k = 0; k < slopes.size(); k++)
                    for (int other_floor = 0; other_floor < (side == 1 ? source_floors : 2); other_floor++)
                    {
                        const Slope& sl = slopes[k];
                        const bool slope_at_source = side == 0;
                        const Vec3 other_end = slope_at_source ? microphone : source;
                        const Vec3 pa = slope_at_source ? sl.image : (other_floor ? below(source) : source);
                        const Vec3 pb = slope_at_source ? (other_floor ? below(microphone) : microphone) : sl.image;
                        const PathEnd::View va = slope_at_source ? view_from(edge, pa, true, 0.0, edge.length)
                                                                 : src.view(0, ki, other_floor);
                        const PathEnd::View vb = slope_at_source ? mic.view(0, next, other_floor)
                                                                 : view_from(edge, pb, true, 0.0, edge.length);
                        double t;
                        if (!keller_point(va.along, va.off, vb.along, vb.off, t))
                            continue;
                        // Past an end, the end stands in for the Keller point
                        // in every check (the first loop has why).
                        const bool past = options.edge_ends && past_true_end(edge, t);
                        const double tq = past ? std::clamp(t, 0.0, edge.length) : t;
                        if (!va.visible || !vb.visible || tq < std::max(va.from, vb.from)
                                || tq > std::min(va.to, vb.to))
                            continue;
                        const Vec3 q = edge.start + tq * edge.direction, keller = edge.start + t * edge.direction;
                        if (q.z < 0.0 || buildings.buried(q)
                                || (past && !end_within_reach(pa, pb, q, keller, options.end_reach)))
                            continue;
                        const Vec3 in = pa - q, out = pb - q;
                        const int on_a = slope_at_source ? -1 : end_face(source_end, edge, {}, other_floor);
                        const int on_b = slope_at_source ? end_face(microphone_end, edge, {}, other_floor) : -1;
                        // Past an end, the end wave goes wherever the end sees
                        // (onto_air).
                        if ((!past
                                    && (!end_air(edge, normal_part(edge, in), on_a)
                                            || !end_air(edge, normal_part(edge, out), on_b)))
                                || doubles_back(edge, on_a, pb) || doubles_back(edge, on_b, pa))
                            continue;
                        std::vector<Vec3>& slope_points = slope_at_source ? to_source : to_microphone;
                        std::vector<Vec3>& plain_points = slope_at_source ? to_microphone : to_source;
                        if (!slope_leg(sl, q, slope_at_source ? source : microphone, slope_points)
                                || !folded_path_open(buildings, q, other_end, {}, other_floor == 1, &plain_points))
                            continue;
                        PointRay ray;
                        const double raw0 = end_angle(edge, in, on_a, ray.source_on_face);
                        const double raw = end_angle(edge, out, on_b, ray.microphone_on_face);
                        const double theta0 = onto_air(raw0, edge.n), theta = onto_air(raw, edge.n);
                        const double s_in = in.norm(), s_out = out.norm();
                        // The coefficient's sin(beta) and distance are the Keller point's,
                        // the same either way round (edge_ends has why).
                        const double s_in_k = (pa - keller).norm(), s_out_k = (pb - keller).norm();
                        const double sin_beta = normal_part(edge, in).norm() / s_in_k;
                        ray.length = s_in + s_out;
                        ray.spreading = 1.0 / std::sqrt(s_in * s_out * (s_in + s_out));
                        // The floor's bounces, along the ray from the source:
                        // on the plain side, its leg from the edge; on the
                        // slope's, the whole unfolded leg through the slope,
                        // the slope a mirror and not a new source.
                        const auto slope_bounce = [&](bool at_microphone) {
                            const Vec3& p = slope_points.front();
                            const Vec3 target = below(at_microphone ? microphone : source);
                            add_floor_bounce(ray, std::abs(target.z - p.z) / (target - p).norm(),
                                    (sl.image - q).norm(), at_microphone);
                        };
                        if (slope_at_source)
                        {
                            if (sl.floor)
                            {
                                slope_bounce(false);
                                // A source on the floor is its own bounce,
                                // left out of the path.
                                ray.launch_mirrored = to_source.size() == 1;
                            }
                            if (other_floor)
                                add_floor_bounce(ray, q, pb, true);
                        }
                        else
                        {
                            if (other_floor)
                            {
                                ray.source_plane = 0.0;
                                add_floor_bounce(ray, pa, q, false);
                            }
                            if (sl.floor)
                                slope_bounce(true);
                        }
                        ray.roof_facades = {buildings.roof_facade(sl.building)};
                        ray.wedges = 1;
                        ray.wedge[0] = {wedge_angles(edge.n, theta, theta0), sin_beta,
                                s_in_k * s_out_k * sin_beta / (s_in_k + s_out_k), false};
                        if (options.edge_ends)
                        {
                            ray.wedge[0].ends = edge_ends(edge, pa, pb, t);
                            weigh_vertex(ray.wedge[0].ends, raw0, raw, edge.n, s_in, s_out);
                        }
                        ray.edge_source.edges[0] = e;
                        ray.edge_source.source = pa;
                        ray.edge_source.microphone = pb;
                        ray.edge_source.stretches[0] = {std::max(va.from, vb.from), std::min(va.to, vb.to)};
                        ray.path = {source};
                        ray.path.insert(ray.path.end(), to_source.rbegin(), to_source.rend());
                        ray.path.push_back(q);
                        ray.path.insert(ray.path.end(), to_microphone.begin(), to_microphone.end());
                        ray.path.push_back(microphone);
                        rays_.push_back(std::move(ray));
                    }
            }
        }
    }

    // The same with facades on either side: the slope between its side's
    // chain and the edge -- that end's image through the chain, after the
    // floor or not, mirrored in the slope's plane -- and the other side any
    // image of its end.  The slope's side to max_edge_reflections with the
    // slope and its floor, the other to max_edge_reflections; both chains
    // empty is the loop above's. The slope's box must cross the box of the
    // edge and the slope's image; the actual chains are checked at the bounce.
    if (options.roofs && options.slope_facades && options.max_edge_reflections > 0)
    {
        struct SlopeImage
        {
            std::size_t image;
            int floor;
            std::uint32_t building;
            std::size_t plane;
            Vec3 point; // the end's image through its chain, after the floor or not
            Vec3 mirrored;
        };
        const auto slope_images = [&](const Ends& end) {
            std::vector<SlopeImage> out;
            thread_local std::vector<std::uint32_t> near_b;
            buildings.near_box(std::min(source.x, microphone.x) - options.radius,
                    std::max(source.x, microphone.x) + options.radius,
                    std::min(source.y, microphone.y) - options.radius,
                    std::max(source.y, microphone.y) + options.radius, near_b);
            std::sort(near_b.begin(), near_b.end());
            for (std::size_t i = 0; i < end.images.size(); i++)
            {
                for (int fl = 0; fl < 2; fl++)
                {
                    if (end.images[i].facades.size() + 1 + fl > options.max_edge_reflections)
                        continue;
                    const Vec3 pt = fl ? below(end.images[i].point) : end.images[i].point;
                    for (const std::uint32_t b : near_b)
                    {
                        const Footprint& f = buildings.footprints()[b];
                        if (f.roof.empty())
                            continue;
                        for (std::size_t pl = 0; pl < f.roof.size(); pl++)
                        {
                            const RoofPlane& plane = f.roof[pl];
                            const double d = plane.above(pt);
                            if (plane.level() || !(d > Buildings::contact))
                                continue;
                            out.push_back({i, fl, b, pl, pt, pt - (2.0 * d) * plane.normal()});
                        }
                    }
                }
            }
            return out;
        };
        const std::vector<SlopeImage> source_slope_images = slope_images(src);
        const std::vector<SlopeImage> microphone_slope_images = slope_images(mic);
        std::vector<Vec3> to_source, to_microphone;
        std::vector<std::size_t> reversed_chain;
        for (int side = 0; side < 2; side++)
        {
            const bool slope_at_source = side == 0;
            if (slope_at_source && !options.source_reflections)
                continue;
            const std::vector<SlopeImage>& slopes = slope_at_source ? source_slope_images : microphone_slope_images;
            const Ends& slope_end = slope_at_source ? src : mic;
            const Ends& plain = slope_at_source ? mic : src;
            const PathEnd& plain_path_end = slope_at_source ? microphone_end : source_end;
            const Vec3& slope_point = slope_at_source ? source : microphone;
            const Vec3& plain_point = slope_at_source ? microphone : source;
            const int plain_floors = slope_at_source ? 2 : source_floors;
            for (const SlopeImage& sl : slopes)
            {
                const std::vector<std::size_t>& chain_a = slope_end.images[sl.image].facades;
                const RoofPlane& plane = buildings.footprints()[sl.building].roof[sl.plane];
                const double* box = buildings.bounds(sl.building);
                for (std::size_t j = 0; j < plain.images.size(); j++)
                {
                    const std::vector<std::size_t>& chain_b = plain.images[j].facades;
                    if ((chain_a.empty() && chain_b.empty()) || chain_b.size() > options.max_edge_reflections)
                        continue;
                    for (std::size_t k = 0; k < plain.edges[j].size(); k++)
                    {
                        const std::size_t e = plain.edges[j][k];
                        const Edge& edge = edges[e];
                        const Vec3 end_point = edge.start + edge.length * edge.direction;
                        // The slope's box against the edge's and its image's.
                        if (std::max({edge.start.x, end_point.x, sl.mirrored.x}) < box[0]
                                || std::min({edge.start.x, end_point.x, sl.mirrored.x}) > box[1]
                                || std::max({edge.start.y, end_point.y, sl.mirrored.y}) < box[2]
                                || std::min({edge.start.y, end_point.y, sl.mirrored.y}) > box[3]
                                || !(std::max(plane.above(edge.start), plane.above(end_point)) > Buildings::contact))
                            continue;
                        for (int fo = 0; fo < plain_floors; fo++)
                        {
                            const Vec3 plain_image = fo ? below(plain.images[j].point) : plain.images[j].point;
                            const Vec3 pa = slope_at_source ? sl.mirrored : plain_image;
                            const Vec3 pb = slope_at_source ? plain_image : sl.mirrored;
                            const PathEnd::View vs = view_from(edge, sl.mirrored, true, 0.0, edge.length);
                            const PathEnd::View& vp = plain.view(j, k, fo);
                            const PathEnd::View& va = slope_at_source ? vs : vp;
                            const PathEnd::View& vb = slope_at_source ? vp : vs;
                            double t;
                            if (!keller_point(va.along, va.off, vb.along, vb.off, t))
                                continue;
                            const bool past = options.edge_ends && past_true_end(edge, t);
                            const double tq = past ? std::clamp(t, 0.0, edge.length) : t;
                            if (!va.visible || !vb.visible || tq < std::max(va.from, vb.from)
                                    || tq > std::min(va.to, vb.to))
                                continue;
                            const Vec3 q = edge.start + tq * edge.direction, keller = edge.start + t * edge.direction;
                            if (q.z < 0.0 || buildings.buried(q)
                                    || (past && !end_within_reach(pa, pb, q, keller, options.end_reach)))
                                continue;
                            const Vec3 in = pa - q, out = pb - q;
                            const int on_plain = end_face(plain_path_end, edge, chain_b, fo);
                            const int on_a = slope_at_source ? -1 : on_plain;
                            const int on_b = slope_at_source ? on_plain : -1;
                            if ((!past
                                        && (!end_air(edge, normal_part(edge, in), on_a)
                                                || !end_air(edge, normal_part(edge, out), on_b)))
                                    || doubles_back(edge, on_a, pb) || doubles_back(edge, on_b, pa))
                                continue;
                            // The slope's side: q to the reflection, and from
                            // there by the chain to its end -- traced from the
                            // source when it is the source's, so that the first
                            // leg is the launch, as a specular ray's.
                            const double dq = plane.above(q), di = -plane.above(sl.mirrored);
                            if (!(dq > Buildings::contact) || !(di > 0.0))
                                continue;
                            const Vec3 p = q + (dq / (dq + di)) * (sl.mirrored - q);
                            if (!buildings.on_plane(sl.building, sl.plane, p.x, p.y) || buildings.blocked(q, p))
                                continue;
                            std::vector<Vec3>& slope_points = slope_at_source ? to_source : to_microphone;
                            std::vector<Vec3>& plain_points = slope_at_source ? to_microphone : to_source;
                            slope_points.clear();
                            bool open;
                            if (slope_at_source)
                            {
                                reversed_chain.assign(chain_a.rbegin(), chain_a.rend());
                                open = physical_chain(buildings, slope_point, p, reversed_chain, sl.floor == 1,
                                        slope_points);
                            }
                            else
                                open = physical_chain(buildings, p, slope_point, chain_a, sl.floor == 1, slope_points);
                            if (!open || !physical_chain(buildings, q, plain_point, chain_b, fo == 1, plain_points))
                                continue;
                            PointRay ray;
                            ray.launch_mirrored = (slope_at_source ? sl.floor : fo)
                                    && source.z <= Buildings::contact;
                            const double raw0 = end_angle(edge, in, on_a, ray.source_on_face);
                            const double raw = end_angle(edge, out, on_b, ray.microphone_on_face);
                            const double theta0 = onto_air(raw0, edge.n), theta = onto_air(raw, edge.n);
                            const double s_in = in.norm(), s_out = out.norm();
                            const double s_in_k = (pa - keller).norm(), s_out_k = (pb - keller).norm();
                            const double sin_beta = normal_part(edge, in).norm() / s_in_k;
                            ray.length = s_in + s_out;
                            ray.spreading = 1.0 / std::sqrt(s_in * s_out * (s_in + s_out));
                            // The floor's bounces, along the ray from the
                            // source: the slope side's over its whole unfolded
                            // leg, the slope a mirror; the plain side's leg from
                            // the edge.
                            const auto slope_bounce = [&](bool at_microphone) {
                                add_floor_bounce(ray, std::abs(sl.point.z - p.z) / (sl.point - p).norm(),
                                        (sl.mirrored - q).norm(), at_microphone);
                            };
                            if (slope_at_source)
                            {
                                if (sl.floor)
                                    slope_bounce(false);
                                if (fo)
                                    add_floor_bounce(ray, q, pb, true);
                            }
                            else
                            {
                                if (fo)
                                {
                                    add_floor_bounce(ray, pa, q, false);
                                }
                                if (sl.floor)
                                    slope_bounce(true);
                            }
                            ray.facades = slope_at_source ? chain_a : chain_b;
                            const std::vector<std::size_t>& later = slope_at_source ? chain_b : chain_a;
                            ray.facades.insert(ray.facades.end(), later.begin(), later.end());
                            ray.roof_facades = {buildings.roof_facade(sl.building)};
                            ray.wedges = 1;
                            ray.wedge[0] = {wedge_angles(edge.n, theta, theta0), sin_beta,
                                    s_in_k * s_out_k * sin_beta / (s_in_k + s_out_k), false};
                            if (options.edge_ends)
                            {
                                ray.wedge[0].ends = edge_ends(edge, pa, pb, t);
                                weigh_vertex(ray.wedge[0].ends, raw0, raw, edge.n, s_in, s_out);
                            }
                            ray.edge_source.edges[0] = e;
                            ray.edge_source.source = pa;
                            ray.edge_source.microphone = pb;
                            ray.edge_source.stretches[0] = {std::max(va.from, vb.from), std::min(va.to, vb.to)};
                            ray.path = {source};
                            if (slope_at_source)
                                ray.path.insert(ray.path.end(), to_source.begin(), to_source.end());
                            else
                                ray.path.insert(ray.path.end(), to_source.rbegin(), to_source.rend());
                            if (slope_at_source)
                                ray.path.push_back(p);
                            ray.path.push_back(q);
                            if (!slope_at_source)
                                ray.path.push_back(p);
                            ray.path.insert(ray.path.end(), to_microphone.begin(), to_microphone.end());
                            ray.path.push_back(microphone);
                            rays_.push_back(std::move(ray));
                        }
                    }
                }
            }
        }
    }

    if (!options.double_diffraction)
        return;

    // Diffracted twice: the source's images and the microphone's to
    // max_twice_reflections.  With t2 given, the stationary point on the first edge is
    // keller's, and with that the second's; the path is stationary where t2
    // is a fixed point of the two in turn.  Most pairs are ruled out first:
    // where the source's image does not light the first edge or the
    // microphone's does not light the second, from either side of the edge
    // whatever the point along it, or cannot see any of it, and where no
    // point of one edge can put the other's stationary point on the stretch
    // of its own segment its end sees.
    // The images to max_twice_reflections, and of each the edges it lights
    // from either floor; with pruning off, every edge.
    std::vector<std::size_t> src_first, mic_first;
    std::vector<std::vector<std::uint32_t>> src_live, mic_live;
    auto first_order = [&](const Ends& end, int floors, std::size_t depth, std::vector<std::size_t>& first,
                               std::vector<std::vector<std::uint32_t>>& live) {
        for (std::size_t i = 0; i < end.images.size(); i++)
        {
            if (end.images[i].facades.size() > depth)
                continue;
            first.push_back(i);
            live.emplace_back();
            for (std::size_t k = 0; k < end.edges[i].size(); k++)
                if (!options.prune || end.view(i, k, 0).lit || (floors == 2 && end.view(i, k, 1).lit))
                    live.back().push_back(static_cast<std::uint32_t>(k));
        }
    };
    first_order(src, source_floors, options.max_twice_reflections, src_first, src_live);
    first_order(mic, 2, options.max_twice_reflections, mic_first, mic_live);
    const double reach_slack = 1e-6;
    thread_local PairTable pairs;
    pairs.start();
    // The facades a pair's ray may reflect in between its edges
    // (max_between_reflections), each with the second edge mirrored in it and
    // the pair's geometry with that mirror.  A facade qualifies where both
    // edges stand in front of it, its heights meet the band between the
    // edges' (the path between them runs within it), and, with every end in
    // front, the plan of a segment from the first edge to the mirror can cross
    // it: the crossings fill the hull of the four from end to end.  Facades
    // farther off than the pair's own extent, at least 1 m, are not tried.
    // The floor is always a candidate, last.  Worked out once a field for each
    // pair met.
    struct Between
    {
        std::size_t facade;
        Edge mirror;
        EdgePair pair;
        bool floor = false;
    };
    std::unordered_map<std::uint64_t, std::vector<Between>> between_lists;
    std::vector<std::uint32_t> nearby;
    auto between_of = [&](std::size_t e1, std::size_t e2) -> const std::vector<Between>& {
        const std::uint64_t key = (static_cast<std::uint64_t>(e1) << 32) | static_cast<std::uint64_t>(e2);
        const auto found = between_lists.find(key);
        if (found != between_lists.end())
            return found->second;
        std::vector<Between>& list = between_lists[key];
        const Edge& g1 = edges[e1];
        const Edge& g2 = edges[e2];
        const Vec3 ends[4] = {g1.start, g1.start + g1.length * g1.direction, g2.start,
                g2.start + g2.length * g2.direction};
        double x0 = ends[0].x, x1 = x0, y0 = ends[0].y, y1 = y0, z0 = ends[0].z, z1 = z0;
        for (const Vec3& e : ends)
        {
            x0 = std::min(x0, e.x);
            x1 = std::max(x1, e.x);
            y0 = std::min(y0, e.y);
            y1 = std::max(y1, e.y);
            z0 = std::min(z0, e.z);
            z1 = std::max(z1, e.z);
        }
        const double grow = std::max({x1 - x0, y1 - y0, 1.0});
        buildings.facades_in_box(x0 - grow, x1 + grow, y0 - grow, y1 + grow, nearby);
        for (const std::uint32_t f : nearby)
        {
            const Facade& facade = buildings.facades()[f];
            if (facade.top < z0 - Buildings::contact || facade.bottom > z1 + Buildings::contact)
                continue;
            double side[4];
            bool front1 = false, front2 = false, all = true;
            for (int k = 0; k < 4; k++)
            {
                side[k] = facade_side(facade, ends[k]);
                (k < 2 ? front1 : front2) |= side[k] > Buildings::contact;
                all = all && side[k] > Buildings::contact;
            }
            if (!front1 || !front2)
                continue;
            if (all)
            {
                const double fx = facade.x1 - facade.x0, fy = facade.y1 - facade.y0, f2 = fx * fx + fy * fy;
                double lo = std::numeric_limits<double>::infinity(), hi = -lo;
                for (int k = 0; k < 2; k++)
                    for (int l = 2; l < 4; l++)
                    {
                        const Vec3 m = mirror_in(facade, ends[l]);
                        const double w = side[k] / (side[k] + side[l]);
                        const double u = ((ends[k].x + w * (m.x - ends[k].x) - facade.x0) * fx
                                                 + (ends[k].y + w * (m.y - ends[k].y) - facade.y0) * fy)
                                / f2;
                        lo = std::min(lo, u);
                        hi = std::max(hi, u);
                    }
                if (hi < -1e-9 || lo > 1.0 + 1e-9)
                    continue;
            }
            Between b {f, mirrored_edge(facade, g2), {}};
            b.pair = edge_pair(g1, b.mirror);
            list.push_back(std::move(b));
        }
        Between floor {Buildings::no_facade, floor_mirrored_edge(g2), {}, true};
        floor.pair = edge_pair(g1, floor.mirror);
        list.push_back(std::move(floor));
        return list;
    };
    std::vector<Vec3> to_source, to_microphone, across;
    // One end's side of a twice-diffracted ray as its edge sees it: an image
    // of the end through its chain, off the floor or not, or the bare end's
    // image in a level roof or a slope (PathEnd::Bounce), after the floor or
    // not for a slope.  `point` is that image, `view` its view of the edge.
    struct TwiceEnd
    {
        Vec3 point;
        const PathEnd::View* view = nullptr;
        const std::vector<std::size_t>* facades = nullptr;
        int floor = 0;
        const PathEnd::Bounce* bounce = nullptr;
    };
    const std::vector<std::size_t> bare;
    // An end's real points from its edge's point q: by its chain, or off its
    // roof.
    const auto end_leg = [&](const TwiceEnd& end, const Vec3& q, const Vec3& point, std::vector<Vec3>& points) {
        if (end.bounce == nullptr)
            return folded_path_open(buildings, q, point, *end.facades, end.floor == 1, &points);
        const PathEnd::Bounce& b = *end.bounce;
        return roof_leg(b.building, b.plane, b.floor, b.image, q, point, points);
    };
    // A roof side's floor, between the roof and the end: its grazing angle
    // from the roof's reflection, and the whole unfolded leg from the edge,
    // the roof a mirror and not a new source.  A source on the floor is its
    // own bounce, left out of the path, and launches into it.
    const auto roof_floor = [&](PointRay& ray, const TwiceEnd& end, double leg, const std::vector<Vec3>& points,
                                    bool at_microphone) {
        if (!end.bounce->floor)
            return;
        const Vec3& p = points.front();
        const Vec3 target = below(at_microphone ? microphone : source);
        add_floor_bounce(ray, std::abs(target.z - p.z) / (target - p).norm(), leg, at_microphone);
        if (!at_microphone)
            ray.launch_mirrored = points.size() == 1;
    };
    // The ray from A by the edges e1 and e2 -- the second mirrored in the
    // facade or the floor between them (via) or not -- to B, where its path
    // is stationary and open.
    const auto twice = [&](const TwiceEnd& A, const TwiceEnd& B, std::size_t e1, std::size_t e2, const Between* via,
                               const EdgePair& pair) {
        const Edge& g1 = edges[e1];
        const Edge& g2 = edges[e2];
        const Edge& h2 = via != nullptr ? via->mirror : g2;
        const PathEnd::View& from_source = *A.view;
        const PathEnd::View& from_microphone = *B.view;
        const Vec3& pa = A.point;
        const Vec3& pb_real = B.point;
        const Vec3 pb = via == nullptr ? pb_real
                : via->floor           ? below(pb_real)
                                       : mirror_in(buildings.facades()[via->facade], pb_real);
        // The fixed point of t2 -> keller(g2, q1(keller(g1, pa, q2(t2))), pb),
        // by Steffensen's method: exact in one step for parallel
        // edges, where the map is affine, and quadratic near the
        // point otherwise; a step that does not bring the
        // residual down falls back to the map itself.  `first`
        // is the first edge's parameter on the way.
        auto map = [&](double x, double& first, double& second) {
            double along, off;
            pair.on1.at(x, along, off);
            if (!keller_point(from_source.along, from_source.off, along, off, first))
                return false;
            pair.on2.at(first, along, off);
            return keller_point(along, off, from_microphone.along, from_microphone.off, second);
        };
        double t1 = 0.0, t2 = 0.5 * h2.length, first = 0.0, x1 = 0.0;
        bool converged = pair.parallel && parallel_stationary(g1, h2, pa, pb, t2);
        bool ok = converged || map(t2, first, x1);
        for (int iteration = 0; iteration < 100 && ok && !converged; iteration++)
        {
            // x1 is the map of t2.
            if (std::abs(x1 - t2) < 1e-10)
            {
                t2 = x1;
                converged = true;
                break;
            }
            double x2;
            ok = map(x1, first, x2);
            if (!ok)
                break;
            const double bend = x2 - 2.0 * x1 + t2;
            const double guess
                    = std::abs(bend) > 1e-300 ? t2 - (x1 - t2) * (x1 - t2) / bend : x2;
            double y;
            if (std::abs(bend) > 1e-300 && std::isfinite(guess) && map(guess, first, y)
                    && std::abs(y - guess) < std::abs(x2 - x1))
            {
                t2 = guess;
                x1 = y;
            }
            else
            {
                // The map of x1 is x2, already found.
                t2 = x1;
                x1 = x2;
            }
        }
        if (!ok || !converged)
            return;
        // The point's two parameters as the alternation leaves them:
        // the first edge's from t2, the second's from that.
        ok = map(t2, t1, t2);
        if (!ok || t1 < from_source.from || t1 > from_source.to || t2 < from_microphone.from
                || t2 > from_microphone.to)
            return;
        // q2 on the edge as the unfolded path meets it, and
        // where the real path does.
        const Vec3 q1 = g1.start + t1 * g1.direction, q2 = h2.start + t2 * h2.direction;
        const Vec3 q2_real = via != nullptr ? g2.start + t2 * g2.direction : q2;
        if (q1.z < 0.0 || q2_real.z < 0.0)
            return;
        if (via != nullptr && !via->floor)
        {
            const Facade& facade = buildings.facades()[via->facade];
            if (facade_side(facade, q1) <= Buildings::contact
                    || facade_side(facade, q2_real) <= Buildings::contact)
                return;
        }
        if (via != nullptr && via->floor && !(q1.z > 0.0 && q2_real.z > 0.0))
            return;
        const Vec3 in1 = pa - q1, mid = q2 - q1, out2 = pb - q2;
        if (mid.norm() < 1e-6)
            return;
        // The ray between the edges along a face both wedges share
        // (a roof between its two edges, a wall between two
        // corners): wedge 1's wave there is already the field
        // doubled by that face, and wedge 2's coefficient, its
        // incident and reflected boundaries coinciding at grazing,
        // would double it again.
        const EndFaces ends {A.bounce != nullptr ? -1 : end_face(source_end, g1, *A.facades, A.floor),
                via == nullptr && B.bounce == nullptr ? end_face(microphone_end, h2, *B.facades, B.floor) : -1};
        if (!end_air(g1, normal_part(g1, in1), ends.source)
                || !air_or_face(g1, normal_part(g1, mid))
                || !air_or_face(h2, normal_part(h2, -1.0 * mid))
                || !end_air(h2, normal_part(h2, out2), ends.microphone)
                || doubles_back(g1, ends.source, q2) || doubles_back(h2, ends.microphone, q1)
                || buildings.buried(q1) || buildings.buried(q2_real))
            return;
        if (!end_leg(A, q1, source, to_source))
            return;
        across.clear();
        if (via == nullptr ? buildings.blocked(q1, q2)
                : via->floor
                ? !folded_path_open(buildings, q1, q2_real, {}, true)
                : !folded_path_open(buildings, q1, q2_real, {via->facade}, false, &across))
            return;
        // The floor's bounce, where the unfolded leg crosses it.
        if (via != nullptr && via->floor)
            across.push_back(q1 + (q1.z / (q1.z + q2_real.z)) * (q2 - q1));
        if (!end_leg(B, q2_real, microphone, to_microphone))
            return;
        PointRay ray = twice_diffracted_ray(g1, h2, pa, pb, q1, q2, ends);
        if (options.edge_ends)
            twice_ends(g1, h2, pa, pb, ray.length, ray.wedge[0].ends, ray.wedge[1].ends);
        ray.facades = *A.facades;
        if (via != nullptr && !via->floor)
        {
            ray.facades.push_back(via->facade);
            ray.between = via->facade;
        }
        ray.between_floor = via != nullptr && via->floor;
        ray.facades.insert(ray.facades.end(), B.facades->begin(), B.facades->end());
        if (A.bounce != nullptr)
            roof_floor(ray, A, (pa - q1).norm(), to_source, false);
        else if (A.floor)
        {
            ray.source_plane = 0.0;
            add_floor_bounce(ray, pa, q1, false);
        }
        // Between the edges the second is the floor's
        // image, and so is the microphone's end.
        if (ray.between_floor)
            add_floor_bounce(ray, q1, q2, false);
        if (B.bounce != nullptr)
            roof_floor(ray, B, (pb - q2).norm(), to_microphone, true);
        else if (B.floor)
            add_floor_bounce(ray, q2, pb, true);
        for (const TwiceEnd* end : {&A, &B})
            if (end->bounce != nullptr)
                ray.roof_facades.push_back(buildings.roof_facade(end->bounce->building));
        ray.edge_source.edges[0] = e1;
        ray.edge_source.edges[1] = e2;
        ray.edge_source.source = pa;
        ray.edge_source.microphone = pb_real;
        ray.edge_source.stretches[0] = {from_source.from, from_source.to};
        ray.edge_source.stretches[1] = {from_microphone.from, from_microphone.to};
        ray.path = {source};
        ray.path.insert(ray.path.end(), to_source.rbegin(), to_source.rend());
        ray.path.push_back(q1);
        ray.path.insert(ray.path.end(), across.begin(), across.end());
        ray.path.push_back(q2_real);
        ray.path.insert(ray.path.end(), to_microphone.begin(), to_microphone.end());
        ray.path.push_back(microphone);
        if (options.shared_face && std::abs(g1.direction.dot(h2.direction)) > 1.0 - 1e-9)
        {
            // The plane normal to the edges: x from edge 1 toward
            // edge 2, z into the air off the face they share, or
            // for two that share none (a ridge over a step down to
            // a lower roof's edge), off the line between them.
            const Vec3 across_face = normal_part(g1, mid);
            const double width = across_face.norm();
            const Vec3 x = across_face / width;
            Vec3 z = g1.direction.cross(x);
            // There each wedge's angles run from the face below
            // the line; the hard wedge is symmetric, so which of
            // its faces that is does not matter.
            SharedFaceTerm term;
            bool pair = false;
            if (ray.wedge[1].grazing)
            {
                if (!in_air(g1.angle(normal_part(g1, x + 1e-3 * z)), g1.n, false))
                    z = -1.0 * z;
                term.face = {width, g1.n, h2.n};
                pair = true;
            }
            else
                pair = via == nullptr
                        && tilted_pair(buildings, g1, h2, q1, width, x, z, term.face);
            const Vec3 sa = pa - q1, sb = pb - q1;
            term.sx = sa.dot(x);
            term.sz = sa.dot(z);
            term.rx = sb.dot(x);
            term.rz = sb.dot(z);
            term.sin_beta = ray.wedge[0].sin_beta;
            if (pair && shared_face_applies(term.face, term.sx, term.sz, term.rx, term.rz))
            {
                ray.shared = true;
                ray.shared_face = term;
                ray.cascade_spreading = ray.spreading;
                ray.spreading = 1.0 / std::sqrt(ray.length * ray.wedge[0].sin_beta);
            }
        }
        rays_.push_back(std::move(ray));
    };
    // Every way between e1 and e2 for a few ends on each side: straight, and
    // by each facade or the floor between them; the ends that cannot reach
    // the pair's stretch ruled out first, for each way.
    const auto twice_pair = [&](std::size_t e1, std::size_t e2, const TwiceEnd* as, int na, const TwiceEnd* bs,
                                    int nb) {
        const EdgePair& direct = pairs.get(edges, e1, e2);
        // The ray straight between the edges, then by each facade
        // between them, the second edge mirrored in it (h2) and
        // the microphone's image with it: the path unfolds about
        // the facade into a double diffraction of its own.
        const std::vector<Between>* middle = options.max_between_reflections > 0 ? &between_of(e1, e2) : nullptr;
        const std::size_t ways = 1 + (middle != nullptr ? middle->size() : 0);
        for (std::size_t way = 0; way < ways; way++)
        {
            const Between* via = way > 0 ? &(*middle)[way - 1] : nullptr;
            const EdgePair& pair = via != nullptr ? via->pair : direct;
            // Which of the source's ends can reach the first edge's
            // segment, and the microphone's the second's: the pre-checks
            // for one end do not depend on the other's.  A mirror keeps
            // the microphone's place about its edge.
            bool reach_a[2] = {true, true}, reach_b[2] = {true, true};
            if (options.prune)
            {
                for (int k = 0; k < na; k++)
                {
                    const PathEnd::View& p = *as[k].view;
                    reach_a[k] = p.lit && keller_between(p, pair.on1, reach_slack)
                            && keller_can_reach(p, pair.on1, reach_slack);
                }
                for (int k = 0; k < nb; k++)
                {
                    const PathEnd::View& p = *bs[k].view;
                    reach_b[k] = p.lit && keller_between(p, pair.on2, reach_slack)
                            && keller_can_reach(p, pair.on2, reach_slack);
                }
            }
            for (int ka = 0; ka < na; ka++)
                for (int kb = 0; kb < nb; kb++)
                    if (reach_a[ka] && reach_b[kb])
                        twice(as[ka], bs[kb], e1, e2, via, pair);
        }
    };
    // Each image's ends at its k-th edge, off the floor or not.
    const auto image_ends = [&](const Ends& end, std::size_t i, std::size_t k, int floors, TwiceEnd* out) {
        const PointImage& image = end.images[i];
        for (int fl = 0; fl < floors; fl++)
            out[fl] = {fl ? below(image.point) : image.point, &end.view(i, k, fl), &image.facades, fl, nullptr};
    };
    TwiceEnd as[2], bs[2];
    for (std::size_t ia = 0; ia < src_first.size(); ia++)
        for (std::size_t jb = 0; jb < mic_first.size(); jb++)
        {
            const std::size_t i = src_first[ia], j = mic_first[jb];
            for (const std::uint32_t k1 : src_live[ia])
            {
                const std::size_t e1 = src.edges[i][k1];
                image_ends(src, i, k1, source_floors, as);
                for (const std::uint32_t k2 : mic_live[jb])
                {
                    const std::size_t e2 = mic.edges[j][k2];
                    if (e1 == e2)
                        continue;
                    image_ends(mic, j, k2, 2, bs);
                    twice_pair(e1, e2, as, source_floors, bs, 2);
                }
            }
        }
    // With a level roof or a slope between an end and its edge
    // (options.twice_roofs): each end's bounces (PathEnd::bounces) with the
    // other end's images, as above; not both ends off roofs.
    for (int side = 0; side < 2; side++)
    {
        const bool at_source = side == 0;
        const PathEnd& bouncing = at_source ? source_end : microphone_end;
        const auto& kept = at_source ? src_bounce_at : mic_bounce_at;
        const Ends& other = at_source ? mic : src;
        const auto& other_first = at_source ? mic_first : src_first;
        const auto& other_live = at_source ? mic_live : src_live;
        const int other_floors = at_source ? 2 : source_floors;
        TwiceEnd* roof = at_source ? as : bs;
        TwiceEnd* plain = at_source ? bs : as;
        for (std::size_t r = 0; r < bouncing.bounces().size(); r++)
        {
            const PathEnd::Bounce& bounce = bouncing.bounces()[r];
            for (const std::uint32_t kr : kept[r])
            {
                const std::size_t er = bounce.edges[kr];
                roof[0] = {bounce.image, &bounce.views[kr], &bare, 0, &bounce};
                for (std::size_t io = 0; io < other_first.size(); io++)
                {
                    const std::size_t o = other_first[io];
                    for (const std::uint32_t ko : other_live[io])
                    {
                        const std::size_t eo = other.edges[o][ko];
                        if (eo == er)
                            continue;
                        image_ends(other, o, ko, other_floors, plain);
                        if (at_source)
                            twice_pair(er, eo, roof, 1, plain, other_floors);
                        else
                            twice_pair(eo, er, plain, other_floors, roof, 1);
                    }
                }
            }
        }
    }
    // Diffracted at more than three roof edges: the taut string's corners
    // over the cross-section between each pair of the ends' points (floor
    // images or not), the path over them stationary.
    if (options.diffraction && options.double_diffraction && options.max_chain > 3)
    {
        std::vector<const Edge*> g;
        std::vector<double> t;
        std::vector<Vec3> q, to_chain_source, to_chain_microphone;
        for (int fa = 0; fa < source_floors; fa++)
            for (int fb = 0; fb < 2; fb++)
            {
                const Vec3 pa = fa ? below(source) : source, pb = fb ? below(microphone) : microphone;
                const auto chain = roof_chain(buildings, pa, pb, near_edges);
                if (chain.size() <= 3 || chain.size() > options.max_chain)
                    continue;
                g.clear();
                t.clear();
                for (const auto& [e, along] : chain)
                {
                    g.push_back(&edges[e]);
                    t.push_back(along);
                }
                if (!chain_points(g, pa, pb, t))
                    continue;
                const std::size_t n = g.size();
                q.resize(n);
                bool open = true;
                for (std::size_t i = 0; i < n && open; i++)
                {
                    q[i] = g[i]->start + t[i] * g[i]->direction;
                    open = q[i].z >= 0.0 && !buildings.buried(q[i]);
                }
                const EndFaces ends {end_face(source_end, *g.front(), {}, fa), end_face(microphone_end, *g.back(), {}, fb)};
                for (std::size_t i = 0; i < n && open; i++)
                {
                    const Edge& e = *g[i];
                    const Vec3 in = i == 0 ? pa - q[0] : q[i - 1] - q[i];
                    const Vec3 out = i + 1 == n ? pb - q[i] : q[i + 1] - q[i];
                    open = (i == 0 ? end_air(e, normal_part(e, in), ends.source) : air_or_face(e, normal_part(e, in)))
                            && (i + 1 == n ? end_air(e, normal_part(e, out), ends.microphone)
                                           : air_or_face(e, normal_part(e, out)));
                }
                if (doubles_back(*g.front(), ends.source, q[1]) || doubles_back(*g.back(), ends.microphone, q[n - 2]))
                    continue;
                for (std::size_t i = 0; i + 1 < n && open; i++)
                    open = !buildings.blocked(q[i], q[i + 1]);
                if (!open || !folded_path_open(buildings, q[0], source, {}, fa == 1, &to_chain_source)
                        || !folded_path_open(buildings, q[n - 1], microphone, {}, fb == 1, &to_chain_microphone))
                    continue;
                PointRay ray = chain_diffracted_ray(g, pa, pb, q, ends);
                if (options.edge_ends)
                    chain_ends(g, pa, pb, q, ray.length, ray);
                if (fa)
                {
                    ray.source_plane = 0.0;
                    add_floor_bounce(ray, pa, q[0], false);
                }
                if (fb)
                    add_floor_bounce(ray, q[n - 1], pb, true);
                ray.edge_source.edges[0] = chain.front().first;
                ray.edge_source.edges[1] = chain.back().first;
                ray.edge_source.source = pa;
                ray.edge_source.microphone = pb;
                ray.path = {source};
                ray.path.insert(ray.path.end(), to_chain_source.rbegin(), to_chain_source.rend());
                ray.path.insert(ray.path.end(), q.begin(), q.end());
                ray.path.insert(ray.path.end(), to_chain_microphone.begin(), to_chain_microphone.end());
                ray.path.push_back(microphone);
                rays_.push_back(std::move(ray));
            }
    }

    if (options.triple_diffraction == PointPathOptions::Triple::none)
        return;
    const bool roofs_only = options.triple_diffraction == PointPathOptions::Triple::roofs;

    // Diffracted three times: g1 lit by the source's image, g3 by the
    // microphone's, each to max_thrice_reflections, and g2 between them with
    // nothing reflecting between the edges.  With q2 given, keller's points on
    // g1 (toward the source's image) and g3 (toward the microphone's) are
    // independent of each other, and q2 is then keller's between them: the
    // path is stationary where t2 is a fixed point of that, by Steffensen's
    // method as for two edges.  The middle edges a pair (g1, g3) can take
    // depend on nothing else, and are kept per thread for the district and the
    // rule (roofs or all) from field to field: those of the buildings within the pair's own
    // extent of it (at least 1 m) that see some of each outer edge's air, and
    // whose air each outer edge sees some of -- and that meet neither at a
    // corner.  Around a corner the path's points close on the vertex, which
    // diffracts as a vertex and not as the cascade of its edges, as two edges
    // that meet have their fixed point at the corner itself.
    std::vector<std::size_t> src_third, mic_third;
    std::vector<std::vector<std::uint32_t>> src_third_live, mic_third_live;
    first_order(src, source_floors, options.max_thrice_reflections, src_third, src_third_live);
    first_order(mic, 2, options.max_thrice_reflections, mic_third, mic_third_live);
    struct MiddleCache
    {
        std::uint64_t district = 0;
        bool roofs = false;
        std::unordered_map<std::uint64_t, std::vector<std::size_t>> lists;
        std::vector<std::uint32_t> blocks;
    };
    thread_local MiddleCache cache;
    if (cache.district != buildings.id() || cache.roofs != roofs_only || cache.lists.size() > 1000000)
    {
        cache.lists.clear();
        cache.district = buildings.id();
        cache.roofs = roofs_only;
    }
    auto middle_of = [&](std::size_t e1, std::size_t e3) -> const std::vector<std::size_t>& {
        const std::uint64_t key = (static_cast<std::uint64_t>(e1) << 32) | static_cast<std::uint64_t>(e3);
        const auto found = cache.lists.find(key);
        if (found != cache.lists.end())
            return found->second;
        std::vector<std::size_t>& list = cache.lists[key];
        std::vector<std::uint32_t>& blocks = cache.blocks;
        const Edge& g1 = edges[e1];
        const Edge& g3 = edges[e3];
        const Vec3 ends[4] = {g1.start, g1.start + g1.length * g1.direction, g3.start,
                g3.start + g3.length * g3.direction};
        double x0 = ends[0].x, x1 = x0, y0 = ends[0].y, y1 = y0;
        for (const Vec3& e : ends)
        {
            x0 = std::min(x0, e.x);
            x1 = std::max(x1, e.x);
            y0 = std::min(y0, e.y);
            y1 = std::max(y1, e.y);
        }
        const double grow = std::max({x1 - x0, y1 - y0, 1.0});
        buildings.near_box(x0 - grow, x1 + grow, y0 - grow, y1 + grow, blocks);
        std::sort(blocks.begin(), blocks.end());
        for (const std::uint32_t building : blocks)
            for (std::size_t e2 = buildings.edge_begin(building); e2 < buildings.edge_begin(building + 1); e2++)
            {
                const Edge& g2 = edges[e2];
                if (e2 != e1 && e2 != e3 && !(roofs_only && g2.vertical) && !meet(g1, g2) && !meet(g2, g3)
                        && sees_air(g2, g1) && sees_air(g2, g3) && sees_air(g1, g2) && sees_air(g3, g2))
                    list.push_back(e2);
            }
        return list;
    };
    std::vector<Vec3> to_third_source, to_third_microphone;
    for (std::size_t ia = 0; ia < src_third.size(); ia++)
        for (std::size_t jb = 0; jb < mic_third.size(); jb++)
        {
            const std::size_t i = src_third[ia], j = mic_third[jb];
            const PointImage& a = src.images[i];
            const PointImage& b = mic.images[j];
            for (const std::uint32_t k1 : src_third_live[ia])
            {
                const std::size_t e1 = src.edges[i][k1];
                const Edge& g1 = edges[e1];
                if (roofs_only && g1.vertical)
                    continue;
                for (const std::uint32_t k3 : mic_third_live[jb])
                {
                    const std::size_t e3 = mic.edges[j][k3];
                    const Edge& g3 = edges[e3];
                    if (e1 == e3 || (roofs_only && g3.vertical))
                        continue;
                    for (const std::size_t e2 : middle_of(e1, e3))
                    {
                        const Edge& g2 = edges[e2];
                        // g2's points about g1's and g3's lines, and theirs
                        // about g2's (copies: the table's references last
                        // until its next get).
                        const EdgePair p12 = pairs.get(edges, e1, e2);
                        const EdgePair p32 = pairs.get(edges, e3, e2);
                        bool reach_source[2] = {true, true}, reach_microphone[2] = {true, true};
                        if (options.prune)
                        {
                            for (int fa = 0; fa < source_floors; fa++)
                            {
                                const PathEnd::View& p = src.view(i, k1, fa);
                                reach_source[fa] = p.lit && keller_between(p, p12.on1, reach_slack)
                                        && keller_can_reach(p, p12.on1, reach_slack);
                            }
                            for (int fb = 0; fb < 2; fb++)
                            {
                                const PathEnd::View& p = mic.view(j, k3, fb);
                                reach_microphone[fb] = p.lit && keller_between(p, p32.on1, reach_slack)
                                        && keller_can_reach(p, p32.on1, reach_slack);
                            }
                        }
                        for (int fa = 0; fa < source_floors; fa++)
                            for (int fb = 0; fb < 2; fb++)
                            {
                                if (!reach_source[fa] || !reach_microphone[fb])
                                    continue;
                                const PathEnd::View& from_source = src.view(i, k1, fa);
                                const PathEnd::View& from_microphone = mic.view(j, k3, fb);
                                // The map t2 -> keller(g2, q1, q3), q1 and q3 keller's
                                // toward the ends from q2(t2).
                                auto map = [&](double x, double& u1, double& u3, double& next) {
                                    double along, off, a1, o1, a3, o3;
                                    p12.on1.at(x, along, off);
                                    if (!keller_point(from_source.along, from_source.off, along, off, u1))
                                        return false;
                                    p32.on1.at(x, along, off);
                                    if (!keller_point(along, off, from_microphone.along, from_microphone.off, u3))
                                        return false;
                                    p12.on2.at(u1, a1, o1);
                                    p32.on2.at(u3, a3, o3);
                                    return keller_point(a1, o1, a3, o3, next);
                                };
                                double t1 = 0.0, t3 = 0.0, t2 = 0.5 * g2.length, x1 = 0.0;
                                bool converged = false;
                                bool ok = map(t2, t1, t3, x1);
                                for (int iteration = 0; iteration < 100 && ok && !converged; iteration++)
                                {
                                    if (std::abs(x1 - t2) < 1e-10)
                                    {
                                        t2 = x1;
                                        converged = true;
                                        break;
                                    }
                                    double x2;
                                    ok = map(x1, t1, t3, x2);
                                    if (!ok)
                                        break;
                                    const double bend = x2 - 2.0 * x1 + t2;
                                    const double guess
                                            = std::abs(bend) > 1e-300 ? t2 - (x1 - t2) * (x1 - t2) / bend : x2;
                                    double y;
                                    if (std::abs(bend) > 1e-300 && std::isfinite(guess) && map(guess, t1, t3, y)
                                            && std::abs(y - guess) < std::abs(x2 - x1))
                                    {
                                        t2 = guess;
                                        x1 = y;
                                    }
                                    else
                                    {
                                        t2 = x1;
                                        x1 = x2;
                                    }
                                }
                                if (!ok || !converged)
                                    continue;
                                double settled;
                                if (!map(t2, t1, t3, settled) || t2 < 0.0 || t2 > g2.length || t1 < from_source.from
                                        || t1 > from_source.to || t3 < from_microphone.from || t3 > from_microphone.to)
                                    continue;
                                const Vec3 pa = fa ? below(a.point) : a.point, pb = fb ? below(b.point) : b.point;
                                const Vec3 q1 = g1.start + t1 * g1.direction, q2 = g2.start + t2 * g2.direction,
                                           q3 = g3.start + t3 * g3.direction;
                                if (q1.z < 0.0 || q2.z < 0.0 || q3.z < 0.0)
                                    continue;
                                const Vec3 in1 = pa - q1, m12 = q2 - q1, m23 = q3 - q2, out3 = pb - q3;
                                if (m12.norm() < 1e-6 || m23.norm() < 1e-6)
                                    continue;
                                const EndFaces ends {end_face(source_end, g1, a.facades, fa),
                                        end_face(microphone_end, g3, b.facades, fb)};
                                if (!end_air(g1, normal_part(g1, in1), ends.source)
                                        || !air_or_face(g1, normal_part(g1, m12))
                                        || !air_or_face(g2, normal_part(g2, -1.0 * m12))
                                        || !air_or_face(g2, normal_part(g2, m23))
                                        || !air_or_face(g3, normal_part(g3, -1.0 * m23))
                                        || !end_air(g3, normal_part(g3, out3), ends.microphone)
                                        || doubles_back(g1, ends.source, q2) || doubles_back(g3, ends.microphone, q2)
                                        || buildings.buried(q1) || buildings.buried(q2) || buildings.buried(q3))
                                    continue;
                                if (!folded_path_open(buildings, q1, source, a.facades, fa == 1, &to_third_source)
                                        || buildings.blocked(q1, q2) || buildings.blocked(q2, q3)
                                        || !folded_path_open(
                                                buildings, q3, microphone, b.facades, fb == 1, &to_third_microphone))
                                    continue;
                                PointRay ray = thrice_diffracted_ray(g1, g2, g3, pa, pb, q1, q2, q3, ends);
                                if (options.edge_ends)
                                    chain_ends({&g1, &g2, &g3}, pa, pb, {q1, q2, q3}, ray.length, ray);
                                ray.facades = a.facades;
                                ray.facades.insert(ray.facades.end(), b.facades.begin(), b.facades.end());
                                if (fa)
                                {
                                    ray.source_plane = 0.0;
                                    add_floor_bounce(ray, pa, q1, false);
                                }
                                if (fb)
                                    add_floor_bounce(ray, q3, pb, true);
                                ray.edge_source.edges[0] = e1;
                                ray.edge_source.edges[1] = e3;
                                ray.edge_source.source = pa;
                                ray.edge_source.microphone = pb;
                                ray.path = {source};
                                ray.path.insert(ray.path.end(), to_third_source.rbegin(), to_third_source.rend());
                                ray.path.push_back(q1);
                                ray.path.push_back(q2);
                                ray.path.push_back(q3);
                                ray.path.insert(ray.path.end(), to_third_microphone.begin(), to_third_microphone.end());
                                ray.path.push_back(microphone);
                                rays_.push_back(std::move(ray));
                            }
                    }
                }
            }
        }
}

PointRay twice_diffracted_ray(const Edge& g1, const Edge& g2, const Vec3& pa, const Vec3& pb, const Vec3& q1,
        const Vec3& q2, EndFaces ends)
{
    const Vec3 in1 = pa - q1, mid = q2 - q1, out2 = pb - q2;
    const double s1 = in1.norm(), s2 = mid.norm(), s3 = out2.norm();
    bool leaves_on_face, arrives_on_face, source_on_face, microphone_on_face;
    const double theta0_1 = end_angle(g1, in1, ends.source, source_on_face);
    const double theta_1 = onto_face(g1.angle(normal_part(g1, mid)), g1.n, leaves_on_face);
    const double theta0_2 = onto_face(g2.angle(normal_part(g2, -1.0 * mid)), g2.n, arrives_on_face);
    const bool grazing = leaves_on_face && arrives_on_face;
    const double theta_2 = end_angle(g2, out2, ends.microphone, microphone_on_face);
    const double sb1 = normal_part(g1, in1).norm() / s1;
    const double sb2 = normal_part(g2, out2).norm() / s3;
    // Wedge 2's incident wave is curved about wedge 1 (radius s2) across that
    // edge and about the source (s1 + s2) along it.  Past wedge 2 it spreads by
    // its curvature along wedge 2's own edge (Kouyoumjian and Pathak's rho_e),
    // by Euler's formula between those two -- s1 + s2 for parallel edges, s2
    // for perpendicular ones -- as sqrt(rho_e / (s3 (rho_e + s3))).  The whole
    // spreading is then 1 / sqrt(s1 s3 ((s1 + s2)(s2 + s3) - s1 s3 c^2)), c the
    // cosine between the edges seen along the ray between them: the same from
    // either end.  The distance parameters stay the point-source ones,
    // s1 s2 / (s1 + s2) and s2 s3 / (s2 + s3), which are reciprocal too; the
    // astigmatic one that rho_e gives wedge 2 is not, and it differs only where
    // wedge 2's transition function is far from 1.
    const Vec3 ray_dir = mid / s2;
    const Vec3 along1 = g1.direction - g1.direction.dot(ray_dir) * ray_dir;
    const Vec3 along2 = g2.direction - g2.direction.dot(ray_dir) * ray_dir;
    const double n1 = along1.norm(), n2 = along2.norm();
    const double c2 = n1 > 1e-12 && n2 > 1e-12 ? std::pow(along1.dot(along2) / (n1 * n2), 2) : 1.0;
    PointRay ray;
    ray.length = s1 + s2 + s3;
    ray.spreading = c2 < 1.0 - 1e-15 ? 1.0 / std::sqrt(s1 * s3 * ((s1 + s2) * (s2 + s3) - s1 * s3 * c2))
                                     : 1.0 / std::sqrt(s1 * s2 * s3 * (s1 + s2 + s3));
    ray.wedges = 2;
    // Wedge 1's observer is wedge 2.
    ray.wedge[0] = {wedge_angles(g1.n, theta_1, theta0_1), sb1, s1 * s2 * sb1 / (s1 + s2), false};
    ray.wedge[1] = {wedge_angles(g2.n, theta_2, theta0_2), sb2, s2 * s3 * sb2 / (s2 + s3), grazing};
    ray.source_on_face = source_on_face;
    ray.microphone_on_face = microphone_on_face;
    return ray;
}

PointRay thrice_diffracted_ray(const Edge& g1, const Edge& g2, const Edge& g3, const Vec3& pa, const Vec3& pb,
        const Vec3& q1, const Vec3& q2, const Vec3& q3, EndFaces ends)
{
    return chain_diffracted_ray({&g1, &g2, &g3}, pa, pb, {q1, q2, q3}, ends);
}

PointRay chain_diffracted_ray(const std::vector<const Edge*>& g, const Vec3& pa, const Vec3& pb,
        const std::vector<Vec3>& q, EndFaces ends)
{
    const std::size_t n = g.size();
    // Leg i reaches edge i (leg n the microphone); s[i] its length.
    std::vector<Vec3> leg(n + 1);
    std::vector<double> s(n + 1);
    leg[0] = q[0] - pa;
    for (std::size_t i = 1; i < n; i++)
        leg[i] = q[i] - q[i - 1];
    leg[n] = pb - q[n - 1];
    for (std::size_t i = 0; i <= n; i++)
        s[i] = leg[i].norm();
    auto seen = [](const Edge& e, const Edge& f, const Vec3& l, double length) {
        const Vec3 d = l / length;
        const Vec3 u = e.direction - e.direction.dot(d) * d, v = f.direction - f.direction.dot(d) * d;
        const double nu = u.norm(), nv = v.norm();
        return nu > 1e-12 && nv > 1e-12 ? std::pow(u.dot(v) / (nu * nv), 2) : 1.0;
    };
    auto along_edge = [](double c2, double rho, double l) { return 1.0 / (c2 / (rho + l) + (1.0 - c2) / l); };
    PointRay ray;
    ray.length = 0.0;
    for (double l : s)
        ray.length += l;
    ray.spreading = 1.0 / s[0];
    double rho = s[0];
    bool leaves = false;
    ray.wedges = static_cast<int>(n);
    for (std::size_t i = 0; i < n; i++)
    {
        const Edge& e = *g[i];
        if (i > 0)
            rho = along_edge(seen(*g[i - 1], e, leg[i], s[i]), rho, s[i]);
        ray.spreading *= std::sqrt(rho / (s[i + 1] * (rho + s[i + 1])));
        bool arrives = false, on_face = false;
        const double theta0 = i > 0 ? onto_face(e.angle(normal_part(e, -1.0 * leg[i])), e.n, arrives)
                                    : end_angle(e, -1.0 * leg[0], ends.source, ray.source_on_face);
        const double theta = i + 1 < n ? onto_face(e.angle(normal_part(e, leg[i + 1])), e.n, on_face)
                                        : end_angle(e, leg[n], ends.microphone, ray.microphone_on_face);
        // sin(beta) by the leg on the source's side for the first edge, as
        // thrice_diffracted_ray takes it, and by the leg leaving the rest.
        const double sb = i == 0 ? normal_part(e, leg[0]).norm() / s[0] : normal_part(e, leg[i + 1]).norm() / s[i + 1];
        const WedgeTerm term {wedge_angles(e.n, theta, theta0), sb, s[i] * s[i + 1] * sb / (s[i] + s[i + 1]),
                i > 0 && leaves && arrives};
        if (i < 3)
            ray.wedge[i] = term;
        else
            ray.more_wedges.push_back(term);
        leaves = on_face;
    }
    return ray;
}

bool chain_points(const std::vector<const Edge*>& g, const Vec3& a, const Vec3& b, std::vector<double>& t,
        bool on_edges)
{
    const std::size_t n = g.size();
    if (n == 0 || t.size() != n)
        return false;
    std::vector<double> grad(n), diag(n), off(n), step(n), trial(n);
    auto point = [&](const std::vector<double>& u, std::size_t i) { return g[i]->start + u[i] * g[i]->direction; };
    auto length = [&](const std::vector<double>& u) {
        double l = (point(u, 0) - a).norm() + (b - point(u, n - 1)).norm();
        for (std::size_t i = 0; i + 1 < n; i++)
            l += (point(u, i + 1) - point(u, i)).norm();
        return l;
    };
    const double scale = std::max(1.0, (b - a).norm());
    bool settled = false;
    for (int iteration = 0; iteration < 100 && !settled; iteration++)
    {
        std::fill(grad.begin(), grad.end(), 0.0);
        std::fill(diag.begin(), diag.end(), 1e-12 / scale);
        std::fill(off.begin(), off.end(), 0.0);
        // Each leg from u (edge leg - 1's point, or a) to v (edge leg's, or b):
        // its length's gradient and Hessian in the two parameters.
        for (std::size_t l = 0; l <= n; l++)
        {
            const bool from_edge = l > 0, to_edge = l < n;
            const Vec3 u = from_edge ? point(t, l - 1) : a, v = to_edge ? point(t, l) : b;
            const Vec3 w = v - u;
            const double r = w.norm();
            if (!(r > 0.0))
                return false;
            const Vec3 h = w / r;
            const double cu = from_edge ? g[l - 1]->direction.dot(h) : 0.0;
            const double cv = to_edge ? g[l]->direction.dot(h) : 0.0;
            if (from_edge)
            {
                grad[l - 1] -= cu;
                diag[l - 1] += (1.0 - cu * cu) / r;
            }
            if (to_edge)
            {
                grad[l] += cv;
                diag[l] += (1.0 - cv * cv) / r;
            }
            if (from_edge && to_edge)
                off[l - 1] -= (g[l - 1]->direction.dot(g[l]->direction) - cu * cv) / r;
        }
        // The Newton step, H step = -grad, H tridiagonal (Thomas).
        std::vector<double> c(n), d(n);
        double pivot = diag[0];
        c[0] = n > 1 ? off[0] / pivot : 0.0;
        d[0] = -grad[0] / pivot;
        for (std::size_t i = 1; i < n; i++)
        {
            pivot = diag[i] - off[i - 1] * c[i - 1];
            if (!(pivot > 0.0))
                return false;
            c[i] = i + 1 < n ? off[i] / pivot : 0.0;
            d[i] = (-grad[i] - off[i - 1] * d[i - 1]) / pivot;
        }
        step[n - 1] = d[n - 1];
        for (std::size_t i = n - 1; i-- > 0;)
            step[i] = d[i] - c[i] * step[i + 1];
        // Damped: halve until the length does not grow.
        const double before = length(t);
        double fraction = 1.0, moved = 0.0;
        for (int halving = 0; halving < 40; halving++, fraction *= 0.5)
        {
            for (std::size_t i = 0; i < n; i++)
                trial[i] = t[i] + fraction * step[i];
            if (length(trial) <= before * (1.0 + 1e-15))
                break;
        }
        for (std::size_t i = 0; i < n; i++)
        {
            moved = std::max(moved, std::abs(trial[i] - t[i]));
            t[i] = trial[i];
        }
        settled = moved < 1e-11 * scale;
    }
    if (!settled)
        return false;
    for (std::size_t i = 0; i < n && on_edges; i++)
    {
        const double end = 1e-9 * scale;
        if (!(t[i] > end && t[i] < g[i]->length - end))
            return false;
    }
    return true;
}

Vec3 point_ray_launch(const PointRay& ray)
{
    const Vec3& source = ray.path.front();
    const Vec3& first = ray.path[1];
    Vec3 d;
    const double h = ray.source_plane;
    if (ray.launch_mirrored)
        d = Vec3(first.x - source.x, first.y - source.y, -first.z - source.z);
    else if (ray.wedges == 0)
        d = (ray.path.size() == 2 && !std::isnan(h) ? Vec3(first.x, first.y, 2.0 * h - first.z) : first) - source;
    else if (!std::isnan(h))
        d = Vec3(first.x - source.x, first.y - source.y, 2.0 * h - first.z - source.z);
    else
        d = first - source;
    const double norm = d.norm();
    return norm > 0.0 ? d / norm : d;
}

double reflected_amplitude(const PointRay& ray, const std::vector<double>& reflection, double floor_reflection)
{
    double amplitude = ray.spreading;
    for (std::size_t f : ray.facades)
        amplitude *= std::sqrt(reflection[f]);
    for (int i = 0; i < ray.floor; i++)
        amplitude *= std::sqrt(floor_reflection);
    for (std::size_t f : ray.roof_facades)
        if (f != PointRay::no_facade)
            amplitude *= std::sqrt(reflection[f]);
    return amplitude * end_surface_factor(ray.source_surface, ray.source_on_face, reflection)
            * end_surface_factor(ray.microphone_surface, ray.microphone_on_face, reflection);
}

void point_ray_arrival(const PointRay& ray, double& elevation, double& azimuth)
{
    // A leg off the floor is drawn unfolded below it, which leaves its slope
    // and its horizontal direction as they are.
    const std::size_t n = ray.path.size();
    const Vec3 d = n >= 2 ? ray.path[n - 1] - ray.path[n - 2] : Vec3(1.0, 0.0, 0.0);
    elevation = std::atan2(std::abs(d.z), std::hypot(d.x, d.y));
    azimuth = std::atan2(d.y, d.x);
}

void floor_factors(const PointRay& ray, const FloorModel& floor, const std::vector<double>& k,
        const std::complex<double>* plate, std::complex<double>* out)
{
    const std::size_t problems = k.size();
    if (!floor.impedance())
    {
        std::fill(out, out + problems, std::complex<double>(1.0));
        return;
    }
    if (floor.admittance.size() != problems)
        throw std::invalid_argument("floor_factors: the floor's admittance is not one per sub-band");
    if (!floor.surface.empty() && floor.surface.size() != problems)
        throw std::invalid_argument("floor_factors: the raised surface's reflection is not one per sub-band");
    bool at_microphone = false;
    for (int b = 0; b < ray.floor; b++)
        at_microphone = at_microphone || ray.floor_bounce[static_cast<std::size_t>(b)].microphone;
    for (std::size_t j = 0; j < problems; j++)
    {
        std::complex<double> f = plate == nullptr ? 1.0
                : !floor.surface.empty()
                ? (plate[j] + floor.surface[j] * plate[problems + j]) / (1.0 + floor.surface[j])
                : plate[at_microphone ? problems + j : j];
        for (int b = 0; b < ray.floor; b++)
        {
            const FloorBounce& bounce = ray.floor_bounce[static_cast<std::size_t>(b)];
            const std::array<double, 2> q = spherical_reflection_coefficient(bounce.sin_grazing, 0.5 * k[j] * bounce.leg,
                    floor.admittance[j], floor.boundary_loss);
            double coherent = 1.0;
            if (floor.roughness > 0.0)
            {
                const double x = k[j] * floor.roughness * bounce.sin_grazing;
                coherent = std::exp(-2.0 * x * x);
            }
            f *= coherent * std::complex<double>(q[0], q[1]);
        }
        out[j] = f;
    }
}

bool meets_relief(const PointRay& ray, const std::vector<const FacadeRelief*>* reliefs)
{
    // A ray that diffracts reflects in a relief whole: what it would scatter
    // reaches nothing in the model (the grating rays and the diffuse field's
    // lighting are specular paths), so it carries on the specular way.
    if (reliefs == nullptr || ray.wedges > 0)
        return false;
    for (std::size_t f : ray.facades)
        if (f < reliefs->size() && (*reliefs)[f] != nullptr)
            return true;
    return false;
}

void relief_factors(const PointRay& ray, const std::vector<const FacadeRelief*>& reliefs,
        const std::vector<double>& frequencies, double sound_speed, std::complex<double>* out)
{
    if (ray.hits.size() != ray.facades.size())
        throw std::invalid_argument("relief_factors: the ray was traced without its facade hits (reliefs)");
    for (std::size_t i = 0; i < ray.facades.size(); i++)
    {
        const std::size_t f = ray.facades[i];
        if (f >= reliefs.size() || reliefs[f] == nullptr)
            continue;
        const FacadeHit& hit = ray.hits[i];
        for (std::size_t j = 0; j < frequencies.size(); j++)
            out[j] *= std::sqrt(reliefs[f]->specular(hit.sin_psi, frequencies[j] * hit.cos_elevation, sound_speed));
    }
}

std::complex<double> point_ray_amplitude(const PointRay& ray, const std::vector<double>& reflection,
        double floor_reflection, double k)
{
    return point_ray_amplitude(ray, reflected_amplitude(ray, reflection, floor_reflection), k);
}

std::complex<double> point_ray_amplitude(const PointRay& ray, double reflected, double k)
{
    std::complex<double> amplitude = reflected;
    if (ray.shared)
    {
        const SharedFaceTerm& t = ray.shared_face;
        return shared_face_amplitude(
                ray, reflected, k, shared_face_field(t.face, t.sx, t.sz, t.rx, t.rz, k * t.sin_beta));
    }
    for (int w = 0; w < ray.wedges; w++)
    {
        const WedgeTerm& t = wedge_term(ray, w);
        // The wedge's coefficient in the plane normal to the edge, at the
        // wavenumber and distance projected into it: for a spherical wave the
        // distance is L / sin(beta), L = s s' sin^2(beta) / (s + s').
        amplitude *= ended_wedge_coefficient(t.angles, k * t.sin_beta, t.distance, t.ends, k) / std::sqrt(t.sin_beta)
                * (t.grazing ? 0.5 : 1.0);
    }
    for (const MadeEdge& s : ray.made_edges)
        amplitude *= made_edge_factor(s, k);
    return amplitude;
}

std::complex<double> shared_face_amplitude(const PointRay& ray, double reflected, double k, std::complex<double> field)
{
    // The 2-D field at k sin(beta), its phase along the in-plane length taken
    // out (the ray's length carries it), spread along the edges: a wave
    // exp(i k_t l_t) / sqrt(l_t) in the plane is exp(i k l) / l in space,
    // l = l_t / sin(beta), so the factor is sin(beta) / sqrt(l_t).
    const SharedFaceTerm& t = ray.shared_face;
    const double kt = k * t.sin_beta;
    // The edges' ends, as the cascade's two wedge terms take them: each one's
    // coefficient with its ends over without.
    std::complex<double> ends(1.0, 0.0);
    for (int w = 0; w < 2; w++)
    {
        const WedgeTerm& c = ray.wedge[w];
        if (!std::isfinite(c.ends.lo) && !std::isfinite(c.ends.hi))
            continue;
        const std::complex<double> bare = hard_wedge_coefficient(c.angles, k * c.sin_beta, c.distance);
        if (std::abs(bare) > 0.0)
            ends *= ended_wedge_coefficient(c.angles, k * c.sin_beta, c.distance, c.ends, k) / bare;
    }
    return std::complex<double>(reflected) * t.sin_beta * field * ends * std::polar(1.0, -kt * ray.length * t.sin_beta);
}

std::array<SubBand, 3> sub_bands(double band, double sound_speed)
{
    const double half = std::pow(2.0, 1.0 / 6.0);
    const double f1 = band / half, f2 = band * half;
    std::array<SubBand, 3> parts;
    for (int part = 0; part < 3; part++)
    {
        SubBand& s = parts[part];
        s.g1 = f1 * std::pow(f2 / f1, part / 3.0);
        s.g2 = f1 * std::pow(f2 / f1, (part + 1) / 3.0);
        s.k1 = 2.0 * pi * s.g1 / sound_speed;
        s.k2 = 2.0 * pi * s.g2 / sound_speed;
    }
    return parts;
}

double band_level(double band, const std::array<SubBand, 3>& parts, const double* energy)
{
    const double half = std::pow(2.0, 1.0 / 6.0);
    const double f1 = band / half, f2 = band * half;
    double sum = 0.0;
    for (int part = 0; part < 3; part++)
        sum += energy[part] * (parts[part].g2 - parts[part].g1);
    return 10.0 * std::log10(std::max(sum / (f2 - f1), band_level_floor));
}

double fade_reach(double k_min, double cap)
{
    // nu = sqrt(2 k d / pi) = 3.
    return std::min(4.5 * pi / k_min, cap);
}

std::complex<double> made_edge_factor(const MadeEdge& screen, double k)
{
    // Tapered off between two thirds of the reach and the reach in the
    // Fresnel parameter nu = sqrt(2 k d / pi), at most 3: a partner met beyond
    // the reach is then the same as one never made.
    const double nu = std::sqrt(2.0 * k * screen.detour / pi);
    const double end = std::min(3.0, std::sqrt(2.0 * k * screen.reach / pi)), begin = 2.0 / 3.0 * end;
    double weight = 1.0;
    if (nu >= end)
        return 0.0;
    if (nu > begin)
    {
        const double c = std::cos(0.5 * pi * (nu - begin) / (end - begin));
        weight = c * c;
    }
    return weight * screen.ratio * hard_wedge_coefficient(screen.angles, k * screen.sin_beta, screen.distance)
            / std::sqrt(screen.sin_beta);
}

namespace {

// Where the cascade is trusted: in each wedge's transition parameter
// k L a above 1 (the cascade is within a few hundredths of a dB of the
// infinite edges' integral from 0.3 up, THEORY), and every end of either
// edge and every boundary point of a partner wave more than 1000 radians of
// path excess away (an end at 30-100 radians still moves a ray by a tenth of
// a dB at the median).
constexpr double transition_lo = 0.3, transition_hi = 1.0;
constexpr double excess_lo = 100.0, excess_hi = 1000.0;

} // namespace

namespace {

// The ends' effect on a twice-diffracted ray, at each wavenumber: each edge's
// once-diffracted integral over its whole length over the infinite edge's
// UTD ray, about the two edges' stationary points q1 and q2.  The other
// edge's point is not held there: for each point z1 of the first edge the
// path on through the second is least at a point of its own, which takes the
// path length's curvature in z1 from A to A - C^2 / B (A, B and C its second
// derivatives in z1 and z2 at the stationary points; C couples the edges).
// So the first edge's integral runs from the source to a point on the leg
// toward q2 at the distance whose curvature leaves that: for parallel edges
// the second edge's distance to the receiver beyond q2, the path unfolded
// about it, and then the path through every point of the first edge is the
// one with the second's point free, exactly.  Likewise the second's, from a
// point back along the leg from q1.  Held at q2 instead, the ratio was 1-3 dB
// off the twice-diffracted integral with one end near and two parallel
// edges; this way 0.3 dB at most, with the end within a Fresnel zone.  With
// both ends near, at a corner of the edges' parameters, no product holds (to
// 11 dB off), and the caller takes the integral.  False where the stationary
// points do not settle.
bool end_ratios(const Edge& g1, const Edge& g2, const Vec3& source, const Vec3& receiver,
        const std::vector<double>& wavenumbers, std::vector<std::complex<double>>& r1,
        std::vector<std::complex<double>>& r2)
{
    const EdgeFrameEnds f2 = edge_frame_ends(g2, source, receiver);
    const double z2 = f2.z_s + (f2.z_r - f2.z_s) * f2.r_s / (f2.r_s + f2.r_r);
    const EdgeFrameEnds f1 = edge_frame_ends(g1, source, g2.start + z2 * g2.direction);
    std::vector<double> t {f1.z_s + (f1.z_r - f1.z_s) * f1.r_s / (f1.r_s + f1.r_r), z2};
    if (!chain_points({&g1, &g2}, source, receiver, t, false))
        return false;
    const Vec3 q1 = g1.start + t[0] * g1.direction, q2 = g2.start + t[1] * g2.direction;
    const double m = (q1 - source).norm(), d = (q2 - q1).norm(), l = (receiver - q2).norm();
    const Vec3 u = (q2 - q1) / d;
    const double c1m = g1.direction.dot((q1 - source) / m), c1 = g1.direction.dot(u), c2 = g2.direction.dot(u),
                 c2l = g2.direction.dot((receiver - q2) / l);
    const double s1 = 1.0 - c1 * c1, s2 = 1.0 - c2 * c2;
    const double a = (1.0 - c1m * c1m) / m + s1 / d, b = s2 / d + (1.0 - c2l * c2l) / l;
    const double c = (g1.direction.dot(g2.direction) - c1 * c2) / d;
    // The distances, d at no coupling and growing with it; a plane wave's
    // (curvature 0) is taken at 10^4 d.
    const auto distance = [&](double s, double curvature) {
        constexpr double furthest = 1e4;
        return curvature * furthest * d > s ? s / curvature : furthest * d;
    };
    const Vec3 receiver1 = q1 + distance(s1, s1 / d - c * c / b) * u;
    const Vec3 source2 = q2 - distance(s2, s2 / d - c * c / a) * u;
    r1 = edge_source_spectrum(g1.n, edge_frame_ends(g1, source, receiver1), 0.0, g1.length, wavenumbers);
    r2 = edge_source_spectrum(g2.n, edge_frame_ends(g2, source2, receiver), 0.0, g2.length, wavenumbers);
    constexpr double unbounded = 1e30;
    for (std::size_t j = 0; j < wavenumbers.size(); j++)
    {
        r1[j] /= keller_once_field(g1, source, receiver1, wavenumbers[j], -unbounded, unbounded);
        r2[j] /= keller_once_field(g2, source2, receiver, wavenumbers[j], -unbounded, unbounded);
    }
    return true;
}

} // namespace

std::vector<std::complex<double>> edge_source_factors(const Buildings& buildings, const PointRay& ray,
        const std::vector<double>& wavenumbers)
{
    std::vector<std::complex<double>> out;
    // A ray off a roof takes UTD: the roof's outline bounds the stretch of the
    // edge its image sees, which the integral over the edge does not know.
    // So does a ray that bounces on the floor between its edges, and one
    // diffracted three times; one that reflects in a facade between them takes
    // the integral unfolded about it.
    if (ray.wedges == 0 || ray.wedges > 2 || !ray.made_edges.empty() || !ray.roof_facades.empty()
            || ray.between_floor || wavenumbers.empty())
        return out;
    const EdgeSourceTerm& t = ray.edge_source;
    if (ray.wedges == 1)
    {
        const Edge g = continued(buildings.edges()[t.edges[0]]);
        const EdgeFrameEnds ends = edge_frame_ends(g, t.source, t.microphone);
        out = edge_source_spectrum(g.n, ends, 0.0, g.length, wavenumbers);
        // A ray whose Keller point lies past an end (the tracer keeps it
        // within edge_source_reach) fades out between 50 and 100 radians of
        // the path excess at the end: the tracer's cut is then no step.
        const double apex = ends.z_s + (ends.z_r - ends.z_s) * ends.r_s / (ends.r_s + ends.r_r);
        double excess = 0.0;
        if (apex < 0.0 || apex > g.length)
        {
            const double z = apex < 0.0 ? 0.0 : g.length;
            excess = std::hypot(ends.r_s, z - ends.z_s) + std::hypot(ends.r_r, z - ends.z_r)
                    - std::hypot(ends.r_s + ends.r_r, ends.z_r - ends.z_s);
        }
        for (std::size_t j = 0; j < out.size(); j++)
            out[j] *= (excess > 0.0 ? trust_taper(wavenumbers[j] * excess, 50.0, 100.0) : 1.0)
                    * std::polar(1.0 / ray.spreading, -wavenumbers[j] * ray.length);
        // A factor that is not finite falls back to the ray's UTD, which a ray
        // past an end does not have.
        for (const std::complex<double>& v : out)
            if (!std::isfinite(v.real()) || !std::isfinite(v.imag()))
                return excess > 0.0 ? std::vector<std::complex<double>>(wavenumbers.size(), 0.0)
                                    : std::vector<std::complex<double>> {};
        return out;
    }
    // A ray that reflects in a facade between its edges takes the integral
    // unfolded about it: the second edge and the microphone's image mirrored.
    const bool mirrored = ray.between != PointRay::no_facade;
    const Facade* between = mirrored ? &buildings.facades()[ray.between] : nullptr;
    const Edge g1 = continued(buildings.edges()[t.edges[0]]);
    const Edge g2 = continued(mirrored ? mirrored_edge(*between, buildings.edges()[t.edges[1]])
                                       : buildings.edges()[t.edges[1]]);
    EdgeSourceTerm t_unfolded = t;
    if (mirrored)
        t_unfolded.microphone = mirror_in(*between, t.microphone);
    const EdgeSourceTerm& u = t_unfolded;
    // The weights: the transition parameters' least, per unit wavenumber, and
    // the least path excess of each edge's ends (the other edge's point
    // free) and of a partner's boundary point.
    double transition = std::numeric_limits<double>::infinity();
    for (int w = 0; w < 2; w++)
    {
        double c2 = 1.0;
        for (double c : ray.wedge[w].angles.c)
            c2 = std::min(c2, c * c);
        transition = std::min(transition, ray.wedge[w].sin_beta * ray.wedge[w].distance * 2.0 * c2);
    }
    const auto via_first = [&](double z1) {
        const Vec3 p = g1.start + z1 * g1.direction;
        const EdgeFrameEnds e = edge_frame_ends(g2, p, u.microphone);
        return (u.source - p).norm() + std::hypot(e.r_s + e.r_r, e.z_r - e.z_s);
    };
    const auto via_second = [&](double z2) {
        const Vec3 q = g2.start + z2 * g2.direction;
        const EdgeFrameEnds e = edge_frame_ends(g1, u.source, q);
        return (u.microphone - q).norm() + std::hypot(e.r_s + e.r_r, e.z_r - e.z_s);
    };
    const double first_ends = std::max(std::min(via_first(0.0), via_first(g1.length)) - ray.length, 0.0);
    const double second_ends = std::max(std::min(via_second(0.0), via_second(g2.length)) - ray.length, 0.0);
    double partners = std::numeric_limits<double>::infinity();
    for (const PartnerWave& p : twice_diffracted_partners(g1, g2, u.source, u.microphone))
        for (double z : {p.part.lo, p.part.hi})
            if (z > 0.0 && z < (p.over_second ? g2 : g1).length)
                partners = std::min(partners, std::max((p.over_second ? via_second(z) : via_first(z)) - ray.length, 0.0));
    // Over a shared face always.  A pair that shares none has its transition
    // in the aperture (shared_face_field), which against boundary elements is
    // within 0.5 dB on EDGE's stepped terraces where the integral is 1-3 dB
    // low; it takes the integral for its ends and partners alone.
    const bool always = ray.wedge[1].grazing;
    const bool tilted = ray.shared && !always;
    // Each wavenumber on its own, so that a band's level does not depend on
    // which others the run takes.  Out of both wedges' transitions, with no
    // partner's boundary point within the excess window and either edge's
    // ends past its near side, the cascade holds but for the ends: it takes
    // each edge's ratio (end_ratios), two once-diffracted integrals in place
    // of the twice-diffracted one.  Elsewhere the integral.  Splitting the
    // wavenumbers costs nothing, the integral's cost being its geometry's.
    const std::size_t nk = wavenumbers.size();
    std::vector<double> weight(nk);
    std::vector<double> exact_k, ratio_k;
    std::vector<std::size_t> exact_j, ratio_j;
    for (std::size_t j = 0; j < nk; j++)
    {
        const double k = wavenumbers[j];
        const double in_transition = tilted ? 0.0 : trust_taper(k * transition, transition_lo, transition_hi);
        const double by_partner = trust_taper(k * partners, excess_lo, excess_hi);
        const double by_ends = trust_taper(k * std::min(first_ends, second_ends), excess_lo, excess_hi);
        weight[j] = always ? 1.0 : std::max({in_transition, by_partner, by_ends});
        if (!(weight[j] > 0.0))
            continue;
        if (!always && in_transition == 0.0 && by_partner == 0.0 && k * std::max(first_ends, second_ends) >= excess_lo)
        {
            ratio_k.push_back(k);
            ratio_j.push_back(j);
        }
        else
        {
            exact_k.push_back(k);
            exact_j.push_back(j);
        }
    }
    if (exact_k.empty() && ratio_k.empty())
        return out;
    out.resize(nk);
    for (std::size_t j = 0; j < nk; j++)
        out[j] = weight[j] < 1.0 ? point_ray_amplitude(ray, 1.0, wavenumbers[j]) : 0.0;
    std::vector<std::complex<double>> r1, r2;
    if (!ratio_k.empty() && !end_ratios(g1, g2, u.source, u.microphone, ratio_k, r1, r2))
    {
        // No stationary points to take the ratios about: the integral.
        exact_k.insert(exact_k.end(), ratio_k.begin(), ratio_k.end());
        exact_j.insert(exact_j.end(), ratio_j.begin(), ratio_j.end());
        ratio_j.clear();
    }
    // The ratios are the ends: they take the cascade without its own.
    PointRay bare;
    if (!ratio_j.empty())
    {
        bare = ray;
        bare.wedge[0].ends = bare.wedge[1].ends = EdgeEnds {};
    }
    for (std::size_t i = 0; i < ratio_j.size(); i++)
    {
        const std::size_t j = ratio_j[i];
        out[j] = weight[j] * r1[i] * r2[i] * point_ray_amplitude(bare, 1.0, wavenumbers[j]) + (1.0 - weight[j]) * out[j];
    }
    if (!exact_k.empty())
    {
        const std::vector<std::complex<double>> exact = localized_twice_diffracted_spectrum(g1, g2, u.source,
                u.microphone, ray.wedge[1].grazing ? 0.5 : 1.0, exact_k, PartnerModel::edge_source);
        for (std::size_t i = 0; i < exact_j.size(); i++)
        {
            const std::size_t j = exact_j[i];
            out[j] = weight[j] * exact[i] * std::polar(1.0 / ray.spreading, -wavenumbers[j] * ray.length)
                    + (1.0 - weight[j]) * out[j];
        }
    }
    // A factor that is not finite falls back to the cascade.
    for (const std::complex<double>& v : out)
        if (!std::isfinite(v.real()) || !std::isfinite(v.imag()))
            return {};
    return out;
}

std::vector<double> sub_band_middles(const std::vector<double>& bands, double sound_speed)
{
    std::vector<double> middles;
    for (double band : bands)
        for (const SubBand& s : sub_bands(band, sound_speed))
            middles.push_back(0.5 * (s.g1 + s.g2));
    return middles;
}

std::vector<double> sub_band_absorption(const Atmosphere& air, const std::vector<double>& bands, double sound_speed)
{
    std::vector<double> a = absorption_coefficients(air, sub_band_middles(bands, sound_speed));
    // dB of energy per metre to nepers of amplitude.
    for (double& v : a)
        v *= std::log(10.0) / 20.0;
    return a;
}

namespace {

// A shared face's transition parameter at its far edge, per unit wavenumber:
// X / k, X = k sin(beta) L (1 - cos g), g an end's angle off the face's plane,
// above or below it -- the face's reflection boundary leaves the end's edge
// as far above the face as its shadow boundary does below -- and
// L = s w / (s + w) its distance parameter there, s its distance from its
// edge and w the face's width; the smaller of the two ends'.  An end in the
// plane has X = 0 at every wavenumber, and there the cascade has the wrong
// sign.
double shared_face_transition_rate(const SharedFaceTerm& t)
{
    const double w = t.face.width;
    double rate = std::numeric_limits<double>::infinity();
    auto end = [&](double along, double off) {
        const double s = std::hypot(along, off);
        rate = std::min(rate, t.sin_beta * s * w / (s + w) * (1.0 - along / s));
    };
    if (t.sx < 0.0)
        end(-t.sx, t.sz);
    if (t.rx > w)
        end(t.rx - w, t.rz);
    return rate;
}

// Past this transition parameter the cascade is the aperture's integral to
// 0.04 dB and 0.02 rad, over the test district's shared faces and long blocks
// with grazing sources (where at 0.1 it is 16 dB low).
constexpr double cascade_stands = 30.0;
// The ratio's samples are two thirds of an octave apart, 2^(2/3): an octave
// apart they left 0.13 dB on the long blocks, two thirds 0.07.
constexpr double probe_spacing = 1.5874010519681994;

// A shared face's factor on its reflected amplitude at each wavenumber (what
// shared_face_amplitude gives at reflected = 1), sampled: the cascade of its
// two wedge terms times R, the integral's ratio to the cascade.  R is taken
// by the integral at the lowest wavenumber and then at one probe_spacing or
// more above the last, while the transition parameter is under
// cascade_stands, and at the last such; it is 1 from the first wavenumber
// past that, and in between it goes linearly in log k, in |R|'s log and R's
// unwrapped phase.  Against the integral at every wavenumber this is within
// 0.07 dB and 0.012 rad on long blocks from 20 Hz to 10 kHz (`a sampled shared
// face is its integral at every wavenumber`).  Empty where a probe's cascade
// is zero or the ratio is not finite: the caller takes the integral
// everywhere then.
} // namespace

SharedFaceSamples shared_face_samples(const SharedFaceTerm& t, const std::vector<double>& k)
{
    const std::size_t n = k.size();
    SharedFaceSamples s;
    s.free = n;
    std::vector<std::size_t> by_k(n);
    for (std::size_t j = 0; j < n; j++)
        by_k[j] = j;
    std::sort(by_k.begin(), by_k.end(), [&](std::size_t a, std::size_t b) { return k[a] < k[b]; });
    const double rate = shared_face_transition_rate(t);
    for (std::size_t r = 0; r < n; r++)
    {
        const std::size_t j = by_k[r];
        if (rate * k[j] >= cascade_stands)
        {
            s.free = j;
            if (!s.probes.empty() && s.probes.back() != by_k[r - 1])
                s.probes.push_back(by_k[r - 1]);
            break;
        }
        if (s.probes.empty() || k[j] >= probe_spacing * k[s.probes.back()])
            s.probes.push_back(j);
        if (r + 1 == n && s.probes.back() != j)
            s.probes.push_back(j);
    }
    return s;
}

namespace {

std::vector<std::complex<double>> sampled_shared_face_factors(const PointRay& ray, const std::vector<double>& k)
{
    using complex = std::complex<double>;
    const SharedFaceTerm& t = ray.shared_face;
    const std::size_t n = k.size();
    // The cascade: the ray's own wedge terms, as point_ray_amplitude takes
    // them, at the cascade's spreading.
    PointRay cascade;
    cascade.wedges = ray.wedges;
    cascade.wedge[0] = ray.wedge[0];
    cascade.wedge[1] = ray.wedge[1];
    std::vector<complex> out(n);
    const double ratio = ray.cascade_spreading / ray.spreading;
    for (std::size_t j = 0; j < n; j++)
        out[j] = point_ray_amplitude(cascade, ratio, k[j]);
    const SharedFaceSamples samples = shared_face_samples(t, k);
    const std::vector<std::size_t>& probes = samples.probes;
    const std::size_t first_free = samples.free;
    if (probes.empty())
        return out;
    std::vector<double> in_plane(probes.size());
    for (std::size_t p = 0; p < probes.size(); p++)
        in_plane[p] = k[probes[p]] * t.sin_beta;
    const std::vector<complex> field = shared_face_spectrum(t.face, t.sx, t.sz, t.rx, t.rz, in_plane);
    // The points R goes through: log k, log |R| and R's phase, unwrapped.
    std::vector<double> x, magnitude, phase;
    for (std::size_t p = 0; p < probes.size(); p++)
    {
        const std::size_t j = probes[p];
        const complex r = shared_face_amplitude(ray, 1.0, k[j], field[p]) / out[j];
        if (!(std::isfinite(r.real()) && std::isfinite(r.imag()) && std::abs(r) > 0.0))
            return {};
        double arg = std::arg(r);
        if (!phase.empty())
            arg += 2.0 * pi * std::round((phase.back() - arg) / (2.0 * pi));
        x.push_back(std::log(k[j]));
        magnitude.push_back(std::log(std::abs(r)));
        phase.push_back(arg);
    }
    if (first_free < n)
    {
        x.push_back(std::log(k[first_free]));
        magnitude.push_back(0.0);
        phase.push_back(2.0 * pi * std::round(phase.back() / (2.0 * pi)));
    }
    for (std::size_t j = 0; j < n; j++)
    {
        const double at = std::log(k[j]);
        if (at >= x.back())
        {
            if (first_free < n)
                continue;
            out[j] *= std::polar(std::exp(magnitude.back()), phase.back());
            continue;
        }
        if (at <= x.front())
        {
            out[j] *= std::polar(std::exp(magnitude.front()), phase.front());
            continue;
        }
        const std::size_t hi = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), at) - x.begin());
        const double f = (at - x[hi - 1]) / (x[hi] - x[hi - 1]);
        out[j] *= std::polar(std::exp(magnitude[hi - 1] + f * (magnitude[hi] - magnitude[hi - 1])),
                phase[hi - 1] + f * (phase[hi] - phase[hi - 1]));
    }
    return out;
}

} // namespace

void PointSourceField::sub_band_problems(const std::vector<double>& reflection, const FloorModel& floor,
        const std::vector<double>& bands, double sound_speed, const SubBandProblem& problem,
        const std::vector<double>& absorption, const double* weights, std::vector<std::size_t>* order_out) const
{
    if (!absorption.empty() && absorption.size() != 3 * bands.size())
        throw std::invalid_argument("sub_band_problems: the absorption is not three per band");
    if (floor.impedance() && floor.admittance.size() != 3 * bands.size())
        throw std::invalid_argument("sub_band_problems: the floor's admittance is not three per band");
    const double floor_reflection = floor.impedance() ? 1.0 : floor.reflection;
    // The rays in order of length, by index: band_averaged_energy wants them
    // sorted, and a ray carries its path.
    std::vector<std::size_t> order(rays_.size());
    for (std::size_t i = 0; i < order.size(); i++)
        order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return rays_[a].length < rays_[b].length; });
    if (order_out != nullptr)
        *order_out = order;
    std::vector<double> length(order.size());
    for (std::size_t i = 0; i < order.size(); i++)
        length[i] = rays_[order[i]].length;
    std::vector<std::complex<double>> amplitude(order.size());
    std::vector<double> band_reflection(facade_count_), reflected(order.size());
    // The edge-source integral's factors, every sub-band at once.
    std::vector<std::vector<std::complex<double>>> factors(order.size());
    std::vector<double> wavenumbers;
    for (double band : bands)
        for (const SubBand& part : sub_bands(band, sound_speed))
            wavenumbers.push_back(0.5 * (part.k1 + part.k2));
    if (edge_source_)
        parallel_chunks(order.size(), edge_source_threads_, [&](std::size_t begin, std::size_t end) {
            for (std::size_t i = begin; i < end; i++)
                factors[i] = edge_source_factors(*buildings_, rays_[order[i]], wavenumbers);
        });
    // A floor of finite impedance and the reliefs: each ray's factor from
    // them, every sub-band at once, for the rays they touch (a plate's
    // microphone, every ray).
    std::vector<std::vector<std::complex<double>>> floor_factor(
            floor.impedance() || reliefs_ != nullptr ? order.size() : 0);
    std::vector<std::complex<double>> plate(floor.plate ? 2 * wavenumbers.size() : 0);
    std::vector<double> frequencies(wavenumbers.size());
    for (std::size_t j = 0; j < wavenumbers.size(); j++)
        frequencies[j] = wavenumbers[j] * sound_speed / (2.0 * pi);
    for (std::size_t i = 0; i < floor_factor.size(); i++)
    {
        const PointRay& ray = rays_[order[i]];
        const bool on_floor = floor.impedance() && (ray.floor > 0 || floor.plate);
        const bool on_relief = meets_relief(ray, reliefs_);
        if (!on_floor && !on_relief)
            continue;
        floor_factor[i].assign(wavenumbers.size(), 1.0);
        if (on_floor)
        {
            if (floor.plate)
            {
                double elevation, azimuth;
                point_ray_arrival(ray, elevation, azimuth);
                floor.plate(elevation, azimuth, plate.data(), plate.data() + wavenumbers.size());
            }
            floor_factors(ray, floor, wavenumbers, floor.plate ? plate.data() : nullptr, floor_factor[i].data());
        }
        if (on_relief)
            relief_factors(ray, *reliefs_, frequencies, sound_speed, floor_factor[i].data());
    }
    // A shared face's field at every sub-band at once: sampled where that is
    // allowed and works, or its nodes' geometry shared between them
    // (shared_face_spectrum).
    std::vector<std::vector<std::complex<double>>> shared(order.size());
    std::vector<double> in_plane(wavenumbers.size());
    for (std::size_t i = 0; i < order.size(); i++)
    {
        const PointRay& ray = rays_[order[i]];
        if (!ray.shared || !factors[i].empty())
            continue;
        if (sample_shared_faces_)
        {
            factors[i] = sampled_shared_face_factors(ray, wavenumbers);
            if (!factors[i].empty())
                continue;
        }
        const SharedFaceTerm& t = ray.shared_face;
        for (std::size_t j = 0; j < wavenumbers.size(); j++)
            in_plane[j] = wavenumbers[j] * t.sin_beta;
        shared[i] = shared_face_spectrum(t.face, t.sx, t.sz, t.rx, t.rz, in_plane);
    }
    for (std::size_t band = 0; band < bands.size(); band++)
    {
        for (std::size_t f = 0; f < facade_count_; f++)
            band_reflection[f] = reflection[f * bands.size() + band];
        // What the reflections leave of each ray is the same at every
        // wavenumber of the band.
        for (std::size_t i = 0; i < order.size(); i++)
            reflected[i] = reflected_amplitude(rays_[order[i]], band_reflection, floor_reflection);
        const std::array<SubBand, 3> parts = sub_bands(bands[band], sound_speed);
        for (int part = 0; part < 3; part++)
        {
            const double kc = 0.5 * (parts[part].k1 + parts[part].k2);
            for (std::size_t i = 0; i < order.size(); i++)
            {
                const PointRay& ray = rays_[order[i]];
                if (!factors[i].empty())
                    amplitude[i] = reflected[i] * factors[i][3 * band + part];
                else if (!shared[i].empty())
                    amplitude[i] = shared_face_amplitude(ray, reflected[i], kc, shared[i][3 * band + part]);
                else
                    amplitude[i] = point_ray_amplitude(ray, reflected[i], kc);
            }
            for (std::size_t i = 0; i < floor_factor.size(); i++)
                if (!floor_factor[i].empty())
                    amplitude[i] *= floor_factor[i][3 * band + part];
            if (!absorption.empty())
                for (std::size_t i = 0; i < order.size(); i++)
                    amplitude[i] *= std::exp(-absorption[3 * band + part] * length[i]);
            if (weights != nullptr)
                for (std::size_t i = 0; i < order.size(); i++)
                    amplitude[i] *= weights[order[i] * bands.size() + band];
            problem(band, part, amplitude, length, parts[part]);
        }
    }
}

void PointSourceField::band_levels(const std::vector<double>& reflection, const FloorModel& floor,
        const std::vector<double>& bands, double sound_speed, double* out, const std::vector<double>& absorption,
        const double* weights) const
{
    std::vector<double> energy(3 * bands.size());
    sub_band_problems(
            reflection, floor, bands, sound_speed,
            [&](std::size_t band, int part, const std::vector<std::complex<double>>& amplitude,
                    const std::vector<double>& length, const SubBand& s) {
                energy[3 * band + part] = band_averaged_energy(amplitude, length, s.k1, s.k2);
            },
            absorption, weights);
    for (std::size_t band = 0; band < bands.size(); band++)
        out[band] = band_level(bands[band], sub_bands(bands[band], sound_speed), &energy[3 * band]);
}

} // namespace niceops
