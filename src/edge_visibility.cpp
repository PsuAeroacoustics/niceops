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
// when set, reflects the footprint in that facade's vertical plane first.
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
    planes.push_back({0.0, 0.0, -1.0, -(z_low + margin)});
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

} // namespace

bool visible_stretch(const Buildings& buildings, const Edge& edge, const Vec3& end,
        const std::vector<std::size_t>& facades, bool below_floor, double& from, double& to)
{
    from = 0.0;
    to = edge.length;
    if (facades.size() > 1)
        return true;
    const double length = edge.length;
    const Vec3& q0 = edge.start;
    const Vec3& d = edge.direction;
    const Vec3 qL = q0 + length * d;
    const Vec3 start = below_floor ? Vec3(end.x, end.y, -end.z) : end;

    // The stretch of the edge whose path can meet the facade at all, and the
    // apex of the unfolded segments.
    double v0 = 0.0, v1 = length;
    Vec3 apex = start;
    const Facade* facade = nullptr;
    double fx0 = 0.0, fy0 = 0.0, nx = 0.0, ny = 0.0;
    if (!facades.empty())
    {
        facade = &buildings.facades()[facades[0]];
        const double fl = std::hypot(facade->x1 - facade->x0, facade->y1 - facade->y0);
        const double ax = (facade->x1 - facade->x0) / fl, ay = (facade->y1 - facade->y0) / fl;
        fx0 = facade->x0;
        fy0 = facade->y0;
        nx = -ay;
        ny = ax;
        auto side = [&](const Vec3& p) { return (p.x - fx0) * nx + (p.y - fy0) * ny; };
        const double s_end = side(start);
        apex = Vec3(start.x - 2.0 * s_end * nx, start.y - 2.0 * s_end * ny, start.z);
        const double s_apex = side(apex);
        // folded_path_open: the image behind the plane, the edge's point in
        // front, both by more than contact; taken a little wider here.
        if (!(s_apex < -0.5 * Buildings::contact))
            return false;
        const double sq0 = side(q0), nd = d.x * nx + d.y * ny;
        keep_where(sq0 - 0.5 * Buildings::contact, nd, v0, v1);
        // The hit, x = (-s_apex q + s_q apex) / D, D = s_q - s_apex > 0: its
        // place along the facade and its height, each a linear inequality in
        // t once multiplied through by D.
        const Vec3 n0 = -s_apex * q0 + sq0 * apex, n1 = -s_apex * d + nd * apex;
        const double d0 = sq0 - s_apex, d1 = nd;
        const double slack = 1e-6;
        const double bottom = below_floor ? -facade->top : facade->bottom;
        // along = ((x - p0) . a) >= -slack, <= fl + slack
        const double a0 = (n0.x - fx0 * d0) * ax + (n0.y - fy0 * d0) * ay;
        const double a1 = (n1.x - fx0 * d1) * ax + (n1.y - fy0 * d1) * ay;
        keep_where(a0 + slack * d0, a1 + slack * d1, v0, v1);
        keep_where((fl + slack) * d0 - a0, (fl + slack) * d1 - a1, v0, v1);
        keep_where(n0.z - (bottom - slack) * d0, n1.z - (bottom - slack) * d1, v0, v1);
        keep_where((facade->top + slack) * d0 - n0.z, (facade->top + slack) * d1 - n1.z, v0, v1);
    }
    v0 = std::max(v0, 0.0);
    v1 = std::min(v1, length);
    if (v1 < v0)
        return false;

    // The buildings the unfolded segments may pass through: in front of the
    // facade the buildings themselves, behind it their mirrors (for a direct
    // path, the buildings everywhere).  gather() looks for them in the real
    // world within the box of the triangle it is given: the unfolded one in
    // front, its mirror behind.
    thread_local std::vector<std::pair<double, double>> hidden;
    thread_local std::vector<Plane> planes;
    thread_local std::vector<std::uint32_t> near;
    hidden.clear();
    const double c_facade = nx * fx0 + ny * fy0;
    auto gather = [&](const Vec3& a, const Vec3& b, const Vec3& c, bool behind) {
        const double x0 = std::min({a.x, b.x, c.x}) - margin, x1 = std::max({a.x, b.x, c.x}) + margin;
        const double y0 = std::min({a.y, b.y, c.y}) - margin, y1 = std::max({a.y, b.y, c.y}) + margin;
        buildings.near_box(x0, x1, y0, y1, near);
        for (const std::uint32_t k : near)
        {
            if (!boxes_meet(buildings.bounds(k), x0, x1, y0, y1))
                continue;
            const Footprint& f = buildings.footprints()[k];
            planes.clear();
            if (!prism_planes(f, below_floor ? -f.top : 0.0, behind ? facade : nullptr, planes))
                continue;
            // In front of the facade's plane, or behind it, by the margin.
            if (facade && behind)
                planes.push_back({nx, ny, 0.0, c_facade - margin});
            else if (facade)
                planes.push_back({-nx, -ny, 0.0, -c_facade - margin});
            double lo, hi;
            if (hidden_interval(apex, q0, d, length, planes, lo, hi))
                hidden.push_back({lo, hi});
        }
    };
    gather(apex, q0, qL, false);
    if (facade)
    {
        auto mirrored = [&](const Vec3& p) {
            const double s = (p.x - fx0) * nx + (p.y - fy0) * ny;
            return Vec3(p.x - 2.0 * s * nx, p.y - 2.0 * s * ny, p.z);
        };
        gather(start, mirrored(q0), mirrored(qL), true);
    }
    // What the open intervals, each pulled in a hair, leave of [v0, v1]: in
    // order of their starts, the first point the stretch covered so far does
    // not reach past, and in order of their ends, the last.
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

} // namespace niceops
