#include "edge_visibility.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

#include "urban.hpp"

namespace niceops {

namespace {

// How far inside a building a passage must run to count, m: well past
// Buildings::contact, so that a segment through the shrunken solid has a
// stretch longer than contact strictly inside the real one, which is what
// blocked() refuses.
constexpr double margin = 1e-4;

// A plane n . x <= c, (nx, ny, nz, c).
using Plane = std::array<double, 4>;

struct Point2
{
    double u, w;
};

// Clips a convex polygon by a u + b w <= c (Sutherland-Hodgman).
void clip(std::vector<Point2>& polygon, double a, double b, double c, std::vector<Point2>& scratch)
{
    scratch.clear();
    const std::size_t n = polygon.size();
    for (std::size_t i = 0; i < n; i++)
    {
        const Point2& p = polygon[i];
        const Point2& q = polygon[(i + 1) % n];
        const double fp = a * p.u + b * p.w - c, fq = a * q.u + b * q.w - c;
        if (fp <= 0.0)
            scratch.push_back(p);
        if ((fp < 0.0 && fq > 0.0) || (fp > 0.0 && fq < 0.0))
        {
            const double s = fp / (fp - fq);
            scratch.push_back({p.u + s * (q.u - p.u), p.w + s * (q.w - p.w)});
        }
    }
    polygon.swap(scratch);
}

// The open interval (lo, hi) of t over which the segment from `apex` to
// q0 + t d passes through the convex solid bounded by `planes`, for t over the
// edge and a little past each end, so that an interval can reach over the
// edge's ends.  With x = apex + u (q0 - apex) + w d, u in (0, 1) along the
// segment and w = u t, the solid is a convex polygon in (u, w), and t = w / u
// ranges between its vertices' ratios.  False where it does not pass through
// (to within a sliver of area), or would only at the apex.
bool hidden_interval(const Vec3& apex, const Vec3& q0, const Vec3& d, double length, const std::vector<Plane>& planes,
        double& lo, double& hi)
{
    thread_local std::vector<Point2> polygon, scratch;
    const double past = 1e-6 * std::max(length, 1.0);
    polygon.assign({{0.0, 0.0}, {1.0, -past}, {1.0, length + past}});
    for (const Plane& p : planes)
    {
        const Vec3 n(p[0], p[1], p[2]);
        clip(polygon, n.dot(q0 - apex), n.dot(d), p[3] - n.dot(apex), scratch);
        if (polygon.size() < 3)
            return false;
    }
    double area = 0.0;
    for (std::size_t i = 0; i < polygon.size(); i++)
    {
        const Point2& p = polygon[i];
        const Point2& q = polygon[(i + 1) % polygon.size()];
        area += p.u * q.w - q.u * p.w;
    }
    if (std::abs(area) < 1e-12 * std::max(length, 1.0))
        return false;
    lo = std::numeric_limits<double>::infinity();
    hi = -lo;
    for (const Point2& p : polygon)
    {
        if (p.u < 1e-9)
            return false;
        lo = std::min(lo, p.w / p.u);
        hi = std::max(hi, p.w / p.u);
    }
    return hi > lo;
}

// The planes of a convex footprint's prism from z_low to top, each moved in
// by the margin; false for a footprint not convex or with holes.  `mirror`,
// when set, reflects the footprint in that facade's vertical plane first.  A
// pitched roof's planes stand in for the top (the solid under the least of
// them is the prism cut by each), and below the floor (z_low < 0, the
// building with its floor image) their mirrors in the floor for the bottom.
bool prism_planes(const Footprint& f, double z_low, const Facade* mirror, std::vector<Plane>& planes)
{
    if (!f.holes.empty())
        return false;
    const std::size_t n = f.x.size();
    thread_local std::vector<double> x, y;
    x.assign(f.x.begin(), f.x.end());
    y.assign(f.y.begin(), f.y.end());
    if (mirror)
    {
        const double length = std::hypot(mirror->x1 - mirror->x0, mirror->y1 - mirror->y0);
        const double nx = -(mirror->y1 - mirror->y0) / length, ny = (mirror->x1 - mirror->x0) / length;
        for (std::size_t i = 0; i < n; i++)
        {
            const double s = (x[i] - mirror->x0) * nx + (y[i] - mirror->y0) * ny;
            x[i] -= 2.0 * s * nx;
            y[i] -= 2.0 * s * ny;
        }
        // A mirror turns counterclockwise to clockwise.
        std::reverse(x.begin(), x.end());
        std::reverse(y.begin(), y.end());
    }
    for (std::size_t i = 0; i < n; i++)
    {
        const std::size_t j = (i + 1) % n, k = (i + 2) % n;
        const double cross = (x[j] - x[i]) * (y[k] - y[j]) - (y[j] - y[i]) * (x[k] - x[j]);
        if (cross < 0.0)
            return false;
    }
    for (std::size_t i = 0; i < n; i++)
    {
        const std::size_t j = (i + 1) % n;
        const double ex = x[j] - x[i], ey = y[j] - y[i], length = std::hypot(ex, ey);
        if (length == 0.0)
            continue;
        // Outward for a counterclockwise ring: the edge's right.
        const double nx = ey / length, ny = -ex / length;
        planes.push_back({nx, ny, 0.0, nx * x[i] + ny * y[i] - margin});
    }
    planes.push_back({0.0, 0.0, 1.0, f.top - margin});
    if (f.roof.empty() || z_low >= 0.0)
        planes.push_back({0.0, 0.0, -1.0, -(z_low + margin)});
    // z <= a x + b y + c, and below the floor -z <= a x + b y + c, mirrored
    // (m to m - 2 (m.n) n, the constant less 2 (f0.n)(m.n)) as the footprint.
    double nx = 0.0, ny = 0.0, offset = 0.0;
    if (mirror)
    {
        const double length = std::hypot(mirror->x1 - mirror->x0, mirror->y1 - mirror->y0);
        nx = -(mirror->y1 - mirror->y0) / length;
        ny = (mirror->x1 - mirror->x0) / length;
        offset = mirror->x0 * nx + mirror->y0 * ny;
    }
    for (const RoofPlane& r : f.roof)
        for (const double up : {1.0, -1.0})
        {
            if (up < 0.0 && z_low >= 0.0)
                continue;
            const double norm = std::sqrt(r.a * r.a + r.b * r.b + 1.0);
            const Vec3 normal = r.normal();
            double mx = normal.x, my = normal.y, c = r.c / norm;
            if (mirror)
            {
                const double k = mx * nx + my * ny;
                mx -= 2.0 * k * nx;
                my -= 2.0 * k * ny;
                c -= 2.0 * offset * k;
            }
            planes.push_back({mx, my, up / norm, c - margin});
        }
    return true;
}

// The t at which g0 + g1 t >= 0 holds, intersected into [lo, hi].
void keep_where(double g0, double g1, double& lo, double& hi)
{
    if (g1 == 0.0)
    {
        if (g0 < 0.0)
            hi = lo - 1.0;
        return;
    }
    const double root = -g0 / g1;
    if (g1 > 0.0)
        lo = std::max(lo, root);
    else
        hi = std::min(hi, root);
}

bool boxes_meet(const double* box, double x0, double x1, double y0, double y1)
{
    return box[0] < x1 && box[1] > x0 && box[2] < y1 && box[3] > y0;
}

// A facade's frame: its start, unit direction a along it, length and unit
// normal n = (-a_y, a_x) toward its front.
struct Frame
{
    double x0, y0, ax, ay, length, nx, ny;

    explicit Frame(const Facade& f)
        : x0(f.x0), y0(f.y0), length(std::hypot(f.x1 - f.x0, f.y1 - f.y0))
    {
        ax = (f.x1 - f.x0) / length;
        ay = (f.y1 - f.y0) / length;
        nx = -ay;
        ny = ax;
    }
    double side(const Vec3& p) const { return (p.x - x0) * nx + (p.y - y0) * ny; }
    Vec3 mirror(const Vec3& p) const
    {
        const double s = side(p);
        return Vec3(p.x - 2.0 * s * nx, p.y - 2.0 * s * ny, p.z);
    }
    Vec3 mirror_direction(const Vec3& d) const
    {
        const double s = d.x * nx + d.y * ny;
        return Vec3(d.x - 2.0 * s * nx, d.y - 2.0 * s * ny, d.z);
    }
    // In front of the facade's plane by the margin: n . x >= n . p0 + margin.
    Plane front() const { return {-nx, -ny, 0.0, -(nx * x0 + ny * y0) - margin}; }
};

// The edge's points t whose path can meet facade `f` at all, intersected into
// [v0, v1]: on the line from q0 + t d to the apex, the point in front of the
// facade and the apex behind it, both by more than half of contact
// (folded_path_open asks for contact), and the hit within the facade's
// length and its height, top being the gable's highest and bottom the
// lowest, and below the floor [-top, top], each with 1e-6 m of slack (it asks
// for 1e-9).  False where the apex is not behind the facade.
bool keep_window(const Facade& f, const Frame& fr, const Vec3& q0, const Vec3& d, const Vec3& apex, bool below_floor,
        double& v0, double& v1)
{
    const double s_apex = fr.side(apex);
    if (!(s_apex < -0.5 * Buildings::contact))
        return false;
    const double sq0 = fr.side(q0), nd = d.x * fr.nx + d.y * fr.ny;
    keep_where(sq0 - 0.5 * Buildings::contact, nd, v0, v1);
    // The hit, x = (-s_apex q + s_q apex) / D, D = s_q - s_apex > 0: its
    // place along the facade and its height, each a linear inequality in
    // t once multiplied through by D.
    const Vec3 n0 = -s_apex * q0 + sq0 * apex, n1 = -s_apex * d + nd * apex;
    const double d0 = sq0 - s_apex, d1 = nd;
    const double slack = 1e-6;
    const double bottom = below_floor ? -f.top : f.bottom;
    // along = ((x - p0) . a) >= -slack, <= length + slack
    const double a0 = (n0.x - fr.x0 * d0) * fr.ax + (n0.y - fr.y0 * d0) * fr.ay;
    const double a1 = (n1.x - fr.x0 * d1) * fr.ax + (n1.y - fr.y0 * d1) * fr.ay;
    keep_where(a0 + slack * d0, a1 + slack * d1, v0, v1);
    keep_where((fr.length + slack) * d0 - a0, (fr.length + slack) * d1 - a1, v0, v1);
    keep_where(n0.z - (bottom - slack) * d0, n1.z - (bottom - slack) * d1, v0, v1);
    keep_where((f.top + slack) * d0 - n0.z, (f.top + slack) * d1 - n1.z, v0, v1);
    return true;
}

// What the open intervals `hidden`, each pulled in a hair, leave of [v0, v1]:
// in order of their starts, the first point the stretch covered so far does
// not reach past, and in order of their ends, the last.  False where they
// leave nothing.
bool uncovered(std::vector<std::pair<double, double>>& hidden, double v0, double v1, double length, double& from,
        double& to)
{
    const double hair = 1e-9 * std::max(length, 1.0);
    std::sort(hidden.begin(), hidden.end());
    double first = v0;
    for (const auto& [lo, hi] : hidden)
    {
        if (lo + hair >= first)
            break;
        first = std::max(first, hi - hair);
        if (first > v1)
            return false;
    }
    std::sort(hidden.begin(), hidden.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
    double last = v1;
    for (const auto& [lo, hi] : hidden)
    {
        if (hi - hair <= last)
            break;
        last = std::min(last, lo + hair);
    }
    from = first;
    to = last;
    return true;
}

} // namespace

bool visible_stretch(const Buildings& buildings, const Edge& edge, const Vec3& end,
        const std::vector<std::size_t>& facades, bool below_floor, double& from, double& to)
{
    from = 0.0;
    to = edge.length;
    const double length = edge.length;
    const Vec3 start = below_floor ? Vec3(end.x, end.y, -end.z) : end;
    const auto& all = buildings.facades();
    const std::size_t n = facades.size();

    // The images along the chain, as folded_path_open forms them: targets[k]
    // is the end mirrored in the first k facades.
    thread_local std::vector<Frame> frames;
    thread_local std::vector<Vec3> targets;
    frames.clear();
    targets.assign(1, start);
    for (const std::size_t f : facades)
    {
        frames.emplace_back(all[f]);
        targets.push_back(frames.back().mirror(targets.back()));
    }

    // The stretch of the edge whose path can meet every facade at all.  The
    // path's leg to facade k (from the edge, or from the hit on facade k + 1)
    // lies on the line from the edge's point carried back into that leg,
    // mirrored in facades n - 1 down to k + 1, to targets[k + 1]: each hit is
    // on its facade, which its mirror fixes.
    double v0 = 0.0, v1 = length;
    Vec3 q0 = edge.start, d = edge.direction;
    for (std::size_t k = n; k-- > 0;)
    {
        if (k + 1 < n)
        {
            q0 = frames[k + 1].mirror(q0);
            d = frames[k + 1].mirror_direction(d);
        }
        if (!keep_window(all[facades[k]], frames[k], q0, d, targets[k + 1], below_floor, v0, v1))
            return false;
    }
    v0 = std::max(v0, 0.0);
    v1 = std::min(v1, length);
    if (v1 < v0)
        return false;

    thread_local std::vector<std::pair<double, double>> hidden;
    thread_local std::vector<Plane> planes;
    thread_local std::vector<std::uint32_t> near;
    hidden.clear();
    // The intervals over which the segment from `apex` to p0 + t dp passes
    // through a building near the triangle (a, b, c), each reflected in
    // `mirror` when set and cut by the half-spaces `bounds`.
    auto gather = [&](const Vec3& apex, const Vec3& p0, const Vec3& dp, const Vec3& a, const Vec3& b, const Vec3& c,
            const Facade* mirror, const Plane* bounds, std::size_t n_bounds) {
        const double x0 = std::min({a.x, b.x, c.x}) - margin, x1 = std::max({a.x, b.x, c.x}) + margin;
        const double y0 = std::min({a.y, b.y, c.y}) - margin, y1 = std::max({a.y, b.y, c.y}) + margin;
        buildings.near_box(x0, x1, y0, y1, near);
        for (const std::uint32_t k : near)
        {
            if (!boxes_meet(buildings.bounds(k), x0, x1, y0, y1))
                continue;
            const Footprint& f = buildings.footprints()[k];
            planes.clear();
            if (!prism_planes(f, below_floor ? -f.top : 0.0, mirror, planes))
                continue;
            planes.insert(planes.end(), bounds, bounds + n_bounds);
            double lo, hi;
            if (hidden_interval(apex, p0, dp, length, planes, lo, hi))
                hidden.push_back({lo, hi});
        }
    };
    const Vec3 qL = edge.start + length * edge.direction;
    if (n == 0)
        // A direct path: the buildings everywhere.
        gather(start, edge.start, edge.direction, start, edge.start, qL, nullptr, nullptr, 0);
    else if (n == 1)
    {
        // The unfolded segment from the edge's point to the end's mirror:
        // in front of the facade it passes through the buildings themselves,
        // behind it through their mirrors.  gather() looks for them in the
        // real world within the box of the triangle it is given: the
        // unfolded one in front, its mirror behind.
        const Frame& fr = frames[0];
        const Plane front = fr.front(), behind {fr.nx, fr.ny, 0.0, fr.nx * fr.x0 + fr.ny * fr.y0 - margin};
        const Facade* facade = &all[facades[0]];
        gather(targets[1], edge.start, edge.direction, targets[1], edge.start, qL, nullptr, &front, 1);
        gather(targets[1], edge.start, edge.direction, start, fr.mirror(edge.start), fr.mirror(qL), facade, &behind, 1);
    }
    else
    {
        // A chain: each leg in the real world, on the line from the edge's
        // point carried back into it to targets[k + 1], and the last from the
        // point mirrored in every facade to the end.  The leg runs between
        // the planes of the facades either side of it, in front of both, so
        // the buildings themselves, cut to the front of those, are what can
        // block it.
        Vec3 qq = edge.start, dd = edge.direction;
        for (std::size_t k = n + 1; k-- > 0;)
        {
            // Leg k - 1, between facades k - 1 and k (where they exist).
            if (k < n)
            {
                qq = frames[k].mirror(qq);
                dd = frames[k].mirror_direction(dd);
            }
            Plane bounds[2];
            std::size_t n_bounds = 0;
            if (k < n)
                bounds[n_bounds++] = frames[k].front();
            if (k >= 1)
                bounds[n_bounds++] = frames[k - 1].front();
            gather(targets[k], qq, dd, targets[k], qq, qq + length * dd, nullptr, bounds, n_bounds);
        }
    }
    return uncovered(hidden, v0, v1, length, from, to);
}

} // namespace niceops
