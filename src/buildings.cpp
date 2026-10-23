#include "buildings.hpp"

#include <algorithm>
#include <atomic>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <map>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>

#include "edge_visibility.hpp"
#include "util/constants.hpp"

namespace niceops {

namespace {

double cross2(double ax, double ay, double bx, double by)
{
    return ax * by - ay * bx;
}

// Whether the closed segments p1-p2 and q1-q2 meet.
bool segments_meet(double p1x, double p1y, double p2x, double p2y, double q1x, double q1y, double q2x, double q2y)
{
    const double d1 = cross2(q2x - q1x, q2y - q1y, p1x - q1x, p1y - q1y);
    const double d2 = cross2(q2x - q1x, q2y - q1y, p2x - q1x, p2y - q1y);
    const double d3 = cross2(p2x - p1x, p2y - p1y, q1x - p1x, q1y - p1y);
    const double d4 = cross2(p2x - p1x, p2y - p1y, q2x - p1x, q2y - p1y);
    if (((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0)))
        return true;
    auto on = [](double ax, double ay, double bx, double by, double px, double py, double d) {
        return d == 0.0 && std::min(ax, bx) <= px && px <= std::max(ax, bx) && std::min(ay, by) <= py
                && py <= std::max(ay, by);
    };
    return on(q1x, q1y, q2x, q2y, p1x, p1y, d1) || on(q1x, q1y, q2x, q2y, p2x, p2y, d2)
            || on(p1x, p1y, p2x, p2y, q1x, q1y, d3) || on(p1x, p1y, p2x, p2y, q2x, q2y, d4);
}

double distance_to_side(double px, double py, double ax, double ay, double bx, double by)
{
    const double dx = bx - ax, dy = by - ay;
    const double t = std::clamp(((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy), 0.0, 1.0);
    return std::hypot(px - (ax + t * dx), py - (ay + t * dy));
}

} // namespace

// An edge's line as the ends see it.  A corner standing on the floor goes on
// below it as its mirror image: the floor's images split that one edge between
// a ray and its mirror -- a ray from the source's image to the microphone's is
// the ray from the source to the microphone by the edge's image -- and of each
// such pair the tracer keeps the one whose point is above the floor.  So the
// kept ray takes both halves, and the floor is no end; a wall that an edge
// runs square into does the same (continue_edges).
Edge continued(const Edge& g, double cap)
{
    Edge h = g;
    const double lo = std::min(g.beyond[0], cap), hi = std::min(g.beyond[1], cap);
    h.start = g.start - lo * g.direction;
    h.length = g.length + lo + hi;
    h.beyond[0] = g.beyond[0] - lo;
    h.beyond[1] = g.beyond[1] - hi;
    return h;
}

double Edge::angle(const Vec3& v) const
{
    double theta = std::atan2(v.dot(air), v.dot(face0));
    if (theta < 0.0)
        theta += 2.0 * pi;
    return theta;
}

int Edge::air_by_sign(const Vec3& v) const
{
    if (!(n > 1.00001 && n < 1.99999))
        return -1;
    const double x = v.dot(face0), y = v.dot(air);
    const double margin = 4e-7 * (std::abs(x) + std::abs(y));
    if (!(margin > 0.0))
        return -1;
    if (y > margin)
        return 1;
    if (y < -margin)
    {
        const double across = v.dot(air1);
        return across > margin ? 1 : across < -margin ? 0 : -1;
    }
    return x < -margin ? 1 : -1; // beside the half-turn, theta near pi
}

namespace {

// Every ring of a footprint, the outer one first.
template <typename F, typename Fn>
void for_each_ring(F& f, Fn fn)
{
    fn(f.x, f.y);
    for (auto& hole : f.holes)
        fn(hole.x, hole.y);
}

double signed_area(const std::vector<double>& x, const std::vector<double>& y)
{
    double area = 0.0;
    for (std::size_t i = 0, n = x.size(); i < n; i++)
        area += cross2(x[i], y[i], x[(i + 1) % n], y[(i + 1) % n]);
    return 0.5 * area;
}

// Splits each side at every vertex of another building that lies on it,
// strictly between its ends, where that building runs along it (one of the
// vertex's own sides lies on the same line): a side shared in part -- a
// T-junction, one building's wall against part of another's -- becomes
// sides shared exactly, which the constructor stands against each other.
// A building that touches a side at a vertex alone leaves it whole.  The
// footprints are valid; `bounds` and `grid` are their boxes and the grid on
// them, which a split leaves as they are.  `tolerance` (m) is how far off a
// line a point may lie and still be on it.
void split_at_junctions(std::vector<Footprint>& footprints, const std::vector<double>& bounds,
        const PlanGrid& grid, double tolerance)
{
    struct Insert
    {
        std::size_t ring, side;
        double t, x, y;
    };
    std::vector<std::vector<Insert>> inserts(footprints.size());
    // The vertices of `o` onto the sides of `b`.
    auto take = [&](std::size_t b, std::size_t o) {
        std::size_t ring = 0;
        for_each_ring(footprints[b], [&](const std::vector<double>& x, const std::vector<double>& y) {
            for (std::size_t i = 0, n = x.size(); i < n; i++)
            {
                const std::size_t j = (i + 1) % n;
                const double dx = x[j] - x[i], dy = y[j] - y[i], length = std::hypot(dx, dy);
                if (!(length > 4.0 * tolerance))
                    continue;
                auto off_line = [&](double px, double py) {
                    return std::abs(cross2(px - x[i], py - y[i], dx, dy)) / length;
                };
                for_each_ring(footprints[o], [&](const std::vector<double>& ox, const std::vector<double>& oy) {
                    for (std::size_t k = 0, m = ox.size(); k < m; k++)
                    {
                        const double t = ((ox[k] - x[i]) * dx + (oy[k] - y[i]) * dy) / (length * length);
                        if (!(t * length > 2.0 * tolerance && (1.0 - t) * length > 2.0 * tolerance)
                                || !(off_line(ox[k], oy[k]) < tolerance))
                            continue;
                        const std::size_t before = (k + m - 1) % m, after = (k + 1) % m;
                        if (off_line(ox[before], oy[before]) < tolerance || off_line(ox[after], oy[after]) < tolerance)
                            inserts[b].push_back({ring, i, t, ox[k], oy[k]});
                    }
                });
            }
            ring++;
        });
    };
    // Each pair whose boxes meet, once.
    std::vector<std::size_t> seen(footprints.size(), footprints.size());
    for (std::size_t b = 0; b < footprints.size(); b++)
    {
        const double* box = &bounds[4 * b];
        grid.box(box[0] - tolerance, box[1] + tolerance, box[2] - tolerance, box[3] + tolerance,
                [&](const std::uint32_t* first, const std::uint32_t* last) {
                    for (const std::uint32_t* o = first; o != last; o++)
                        if (*o > b && seen[*o] != b)
                        {
                            seen[*o] = b;
                            take(b, *o);
                            take(*o, b);
                        }
                });
    }
    for (std::size_t b = 0; b < footprints.size(); b++)
    {
        std::vector<Insert>& add = inserts[b];
        if (add.empty())
            continue;
        std::sort(add.begin(), add.end(), [](const Insert& u, const Insert& v) {
            return std::tie(u.ring, u.side, u.t) < std::tie(v.ring, v.side, v.t);
        });
        std::size_t ring = 0, next = 0;
        for_each_ring(footprints[b], [&](std::vector<double>& x, std::vector<double>& y) {
            std::vector<double> nx, ny;
            for (std::size_t i = 0; i < x.size(); i++)
            {
                nx.push_back(x[i]);
                ny.push_back(y[i]);
                for (; next < add.size() && add[next].ring == ring && add[next].side == i; next++)
                    // Two neighbours' vertices at one point are one vertex.
                    if (std::hypot(add[next].x - nx.back(), add[next].y - ny.back()) > 2.0 * tolerance)
                    {
                        nx.push_back(add[next].x);
                        ny.push_back(add[next].y);
                    }
            }
            x = std::move(nx);
            y = std::move(ny);
            ring++;
        });
    }
}

// Whether (px, py) is inside the ring (even-odd, boundary unspecified).
bool in_ring(const std::vector<double>& x, const std::vector<double>& y, double px, double py)
{
    bool in = false;
    for (std::size_t i = 0, n = x.size(), h = n - 1; i < n; h = i++)
        if ((y[i] > py) != (y[h] > py) && px < x[h] + (py - y[h]) * (x[i] - x[h]) / (y[i] - y[h]))
            in = !in;
    return in;
}

// Whether (px, py) is in a footprint: in its outer ring and none of its holes.
bool in_footprint(const Footprint& f, double px, double py)
{
    if (!in_ring(f.x, f.y, px, py))
        return false;
    for (const Ring& hole : f.holes)
        if (in_ring(hole.x, hole.y, px, py))
            return false;
    return true;
}

// The plane of a pitched roof that is least over (x, y).
std::size_t active_plane(const Footprint& f, double x, double y)
{
    std::size_t best = 0;
    for (std::size_t i = 1; i < f.roof.size(); i++)
        if (f.roof[i].at(x, y) < f.roof[best].at(x, y))
            best = i;
    return best;
}

// A footprint that is not convex cut into convex pieces, each with its top and
// roof, counterclockwise; empty for one that is convex, or has courtyards,
// which is left whole.  Ears are clipped (a left turn whose triangle holds no
// other vertex, its sides included) until what is left is convex, and then
// two pieces that share a diagonal are merged wherever the merged ring still
// turns left throughout.  Every piece lies within the footprint and under its
// roof, so whatever passes through a piece passes through the building; a
// clip that stalls on vertices in line keeps what it has cut, which only
// leaves part of the building out.
std::vector<Footprint> convex_pieces_of(const Footprint& f)
{
    const std::size_t n = f.x.size();
    auto turn = [&](std::size_t a, std::size_t b, std::size_t c) {
        return cross2(f.x[b] - f.x[a], f.y[b] - f.y[a], f.x[c] - f.x[b], f.y[c] - f.y[b]);
    };
    auto convex = [&](const std::vector<std::size_t>& ring) {
        const std::size_t m = ring.size();
        for (std::size_t i = 0; i < m; i++)
            if (turn(ring[i], ring[(i + 1) % m], ring[(i + 2) % m]) < 0.0)
                return false;
        return true;
    };
    std::vector<std::size_t> left(n);
    for (std::size_t i = 0; i < n; i++)
        left[i] = i;
    if (!f.holes.empty() || convex(left))
        return {};
    std::vector<std::vector<std::size_t>> pieces;
    while (left.size() > 3 && !convex(left))
    {
        const std::size_t m = left.size();
        std::size_t ear = m;
        for (std::size_t i = 0; i < m && ear == m; i++)
        {
            const std::size_t a = left[(i + m - 1) % m], b = left[i], c = left[(i + 1) % m];
            if (!(turn(a, b, c) > 0.0))
                continue;
            ear = i;
            for (const std::size_t p : left)
                if (p != a && p != b && p != c && turn(a, b, p) >= 0.0 && turn(b, c, p) >= 0.0
                        && turn(c, a, p) >= 0.0)
                {
                    ear = m;
                    break;
                }
        }
        if (ear == m)
            break;
        pieces.push_back({left[(ear + m - 1) % m], left[ear], left[(ear + 1) % m]});
        left.erase(left.begin() + static_cast<std::ptrdiff_t>(ear));
    }
    if (convex(left))
        pieces.push_back(left);
    // Merge across shared diagonals: piece p runs from i to j where piece q
    // runs from j to i.
    for (bool merged = true; merged;)
    {
        merged = false;
        for (std::size_t p = 0; p < pieces.size() && !merged; p++)
            for (std::size_t q = p + 1; q < pieces.size() && !merged; q++)
                for (std::size_t s = 0; s < pieces[p].size() && !merged; s++)
                {
                    const std::vector<std::size_t>& a = pieces[p];
                    const std::vector<std::size_t>& b = pieces[q];
                    const std::size_t i = a[s], j = a[(s + 1) % a.size()];
                    const auto at = std::find(b.begin(), b.end(), j);
                    if (at == b.end())
                        continue;
                    const std::size_t t = static_cast<std::size_t>(at - b.begin());
                    if (b[(t + 1) % b.size()] != i)
                        continue;
                    // p from j round to i, then q's vertices after i up to j.
                    std::vector<std::size_t> ring;
                    for (std::size_t k = 1; k <= a.size(); k++)
                        ring.push_back(a[(s + k) % a.size()]);
                    for (std::size_t k = 2; k < b.size(); k++)
                        ring.push_back(b[(t + k) % b.size()]);
                    if (!convex(ring))
                        continue;
                    pieces[p] = std::move(ring);
                    pieces.erase(pieces.begin() + static_cast<std::ptrdiff_t>(q));
                    merged = true;
                }
    }
    std::vector<Footprint> out;
    for (const std::vector<std::size_t>& ring : pieces)
    {
        Footprint piece;
        for (const std::size_t i : ring)
        {
            piece.x.push_back(f.x[i]);
            piece.y.push_back(f.y[i]);
        }
        piece.top = f.top;
        piece.roof = f.roof;
        out.push_back(std::move(piece));
    }
    return out;
}

} // namespace

std::vector<std::tuple<double, double, std::size_t>> roof_stretches(const Footprint& f, double ax, double ay,
        double bx, double by)
{
    std::vector<double> cuts {0.0, 1.0};
    for (std::size_t k = 0; k < f.roof.size(); k++)
        for (std::size_t l = k + 1; l < f.roof.size(); l++)
        {
            const double g0 = f.roof[k].at(ax, ay) - f.roof[l].at(ax, ay);
            const double g1 = f.roof[k].at(bx, by) - f.roof[l].at(bx, by);
            if ((g0 < 0.0) != (g1 < 0.0) && g0 != g1)
            {
                const double t = g0 / (g0 - g1);
                if (t > 1e-12 && t < 1.0 - 1e-12)
                    cuts.push_back(t);
            }
        }
    std::sort(cuts.begin(), cuts.end());
    std::vector<std::tuple<double, double, std::size_t>> out;
    for (std::size_t i = 0; i + 1 < cuts.size(); i++)
    {
        if (cuts[i + 1] - cuts[i] < 1e-12)
            continue;
        const double t = 0.5 * (cuts[i] + cuts[i + 1]);
        const std::size_t plane = active_plane(f, ax + t * (bx - ax), ay + t * (by - ay));
        if (!out.empty() && std::get<2>(out.back()) == plane)
            std::get<1>(out.back()) = cuts[i + 1];
        else
            out.emplace_back(cuts[i], cuts[i + 1], plane);
    }
    return out;
}

namespace {

// A pitched roof's ridges: where two of its planes cross and are both the
// least, over the footprint.  The roof is the least of its planes, so every
// one is convex: a ridge or a hip.
struct Ridge
{
    Vec3 a, b;
    std::size_t k = 0, l = 0;
};
std::vector<Ridge> roof_ridges(const Footprint& f)
{
    std::vector<Ridge> out;
    for (std::size_t k = 0; k < f.roof.size(); k++)
        for (std::size_t l = k + 1; l < f.roof.size(); l++)
        {
            const RoofPlane& pk = f.roof[k];
            const RoofPlane& pl = f.roof[l];
            const double nx = pk.a - pl.a, ny = pk.b - pl.b, nc = pk.c - pl.c, n2 = nx * nx + ny * ny;
            if (n2 < 1e-18)
                continue;
            // The line where they are equal: p0 + s d in plan.
            const double px = -nc * nx / n2, py = -nc * ny / n2, norm = std::sqrt(n2);
            const double dx = -ny / norm, dy = nx / norm;
            // Where both are the least: each other plane at least as high.
            double s_lo = -1e12, s_hi = 1e12;
            bool empty = false;
            for (std::size_t m = 0; m < f.roof.size() && !empty; m++)
            {
                if (m == k || m == l)
                    continue;
                // pk - pm <= 0 along the line: alpha + beta s <= 0.
                const double alpha = pk.at(px, py) - f.roof[m].at(px, py);
                const double beta = (pk.a - f.roof[m].a) * dx + (pk.b - f.roof[m].b) * dy;
                if (std::abs(beta) < 1e-15)
                    empty = alpha > 1e-12;
                else if (beta > 0.0)
                    s_hi = std::min(s_hi, -alpha / beta);
                else
                    s_lo = std::max(s_lo, -alpha / beta);
            }
            if (empty || !(s_lo < s_hi))
                continue;
            // Where the line crosses the footprint's sides; the stretches
            // between whose middles lie in the footprint.
            std::vector<double> cuts {s_lo, s_hi};
            for_each_ring(f, [&](const std::vector<double>& rx, const std::vector<double>& ry) {
                for (std::size_t i = 0, n = rx.size(), h = n - 1; i < n; h = i++)
                {
                    const double ex = rx[i] - rx[h], ey = ry[i] - ry[h];
                    const double denominator = cross2(dx, dy, ex, ey);
                    if (denominator == 0.0)
                        continue;
                    const double qx = rx[h] - px, qy = ry[h] - py;
                    const double s = cross2(qx, qy, ex, ey) / denominator;
                    const double u = cross2(qx, qy, dx, dy) / denominator;
                    if (u >= -1e-12 && u <= 1.0 + 1e-12 && s > s_lo && s < s_hi)
                        cuts.push_back(s);
                }
            });
            std::sort(cuts.begin(), cuts.end());
            // A ridge has one of its two planes the least on either side: a
            // third plane through the same line (a level one at the ridge's
            // height) makes no second ridge there.
            const auto least_is = [&](std::size_t plane, double x, double y) {
                const double z = f.roof[plane].at(x, y);
                for (const RoofPlane& other : f.roof)
                    if (other.at(x, y) < z - 1e-12)
                        return false;
                return true;
            };
            const double ex = nx / norm * Buildings::contact, ey = ny / norm * Buildings::contact;
            for (std::size_t i = 0; i + 1 < cuts.size(); i++)
            {
                const double s0 = cuts[i], s1 = cuts[i + 1], s = 0.5 * (s0 + s1);
                if (s1 - s0 < Buildings::contact || !in_footprint(f, px + s * dx, py + s * dy))
                    continue;
                // pk - pl grows along (nx, ny): l is the lower that way.
                const double mx = px + s * dx, my = py + s * dy;
                if (!least_is(l, mx + ex, my + ey) || !least_is(k, mx - ex, my - ey))
                    continue;
                Ridge r;
                r.a = Vec3(px + s0 * dx, py + s0 * dy, 0.0);
                r.b = Vec3(px + s1 * dx, py + s1 * dy, 0.0);
                r.a.z = pk.at(r.a.x, r.a.y);
                r.b.z = pk.at(r.b.x, r.b.y);
                r.k = k;
                r.l = l;
                // Joined to the one before where they meet (a side cut it).
                if (!out.empty() && out.back().k == k && out.back().l == l && (out.back().b - r.a).norm() < 1e-9)
                    out.back().b = r.b;
                else
                    out.push_back(r);
            }
        }
    return out;
}

// An edge of a pitched roof from `a` to `b`: its first face in `plane`,
// leaving the edge toward the plan direction `inward`, the air above the
// plane, and the second face leaving along `second`.
Edge roof_edge(const Vec3& a, const Vec3& b, const RoofPlane& plane, const Vec3& inward, const Vec3& second,
        std::size_t building)
{
    Edge e;
    e.start = a;
    e.length = (b - a).norm();
    e.direction = (b - a) / e.length;
    const Vec3 in(inward.x, inward.y, plane.a * inward.x + plane.b * inward.y);
    const Vec3 face0 = in - in.dot(e.direction) * e.direction;
    e.face0 = face0 / face0.norm();
    e.air = plane.normal();
    const Vec3 other = second - second.dot(e.direction) * e.direction;
    e.n = e.angle(other / other.norm()) / pi;
    e.building = building;
    return e;
}

} // namespace

// Over its outline's vertices, where its ridges end, and where two planes
// cross over a side.
double roof_top(const Footprint& f)
{
    double top = -std::numeric_limits<double>::infinity();
    for_each_ring(f, [&](const std::vector<double>& rx, const std::vector<double>& ry) {
        for (std::size_t i = 0, n = rx.size(); i < n; i++)
        {
            const std::size_t j = (i + 1) % n;
            for (const auto& [t0, t1, plane] : roof_stretches(f, rx[i], ry[i], rx[j], ry[j]))
                for (double t : {t0, t1})
                    top = std::max(top, roof_height(f, rx[i] + t * (rx[j] - rx[i]), ry[i] + t * (ry[j] - ry[i])));
        }
    });
    for (const Ridge& r : roof_ridges(f))
        top = std::max({top, r.a.z, r.b.z});
    return top;
}

struct Buildings::FrontTable
{
    FrontTable(std::size_t edges, std::size_t words)
        : formed(new std::once_flag[edges]), bits(new std::uint64_t[edges * words]), words(words)
    {
    }
    // Per edge, whether its row is formed; the rows are left unwritten
    // until they are, so that a district's untouched edges cost no memory.
    std::unique_ptr<std::once_flag[]> formed;
    std::unique_ptr<std::uint64_t[]> bits;
    std::size_t words;
};

Buildings::Buildings(std::vector<Footprint> footprints, double orientation_tolerance)
    : footprints_(std::move(footprints))
    , orientation_tolerance_(orientation_tolerance)
{
    static std::atomic<std::uint64_t> next_id {1};
    id_ = next_id++;
    ceded_.resize(footprints_.size());
    x_min_ = y_min_ = std::numeric_limits<double>::infinity();
    x_max_ = y_max_ = -std::numeric_limits<double>::infinity();
    // Each side as built, for finding the ones two buildings share: its
    // ends, its facade and roof edge, and the corner edges at its ends.
    constexpr std::size_t none = static_cast<std::size_t>(-1);
    struct Side
    {
        std::size_t building;
        double x0, y0, x1, y1;
        // Its facades and roof edges, [first, end): one each under a flat
        // roof, one a stretch of one plane under a pitched one.
        std::size_t facade, facade_end, roof, roof_end, corner0, corner1;
    };
    std::vector<Side> sides;
    for (std::size_t b = 0; b < footprints_.size(); b++)
    {
        Footprint& f = footprints_[b];
        const std::string which = "building " + std::to_string(b + 1);
        if (f.roof.empty() && (!(f.top > 0.0) || !std::isfinite(f.top)))
            throw std::invalid_argument(which + ": its roof must be above the floor");
        for (const RoofPlane& plane : f.roof)
            if (!std::isfinite(plane.a) || !std::isfinite(plane.b) || !std::isfinite(plane.c))
                throw std::invalid_argument(which + ": a roof plane is not a number");
        bool outer = true;
        auto check = [&](std::vector<double>& x, std::vector<double>& y) {
            const std::size_t n = x.size();
            if (n < 3 || y.size() != n)
                throw std::invalid_argument(which + ": a ring needs at least three vertices");
            for (std::size_t i = 0; i < n; i++)
            {
                if (!std::isfinite(x[i]) || !std::isfinite(y[i]))
                    throw std::invalid_argument(which + ": a vertex is not a number");
                if (std::hypot(x[(i + 1) % n] - x[i], y[(i + 1) % n] - y[i]) < contact)
                    throw std::invalid_argument(which + ": two consecutive vertices coincide");
            }
            const double area = signed_area(x, y);
            if (!(std::abs(area) > 0.0))
                throw std::invalid_argument(which + ": a ring has no area");
            // The solid on the left: the outer ring counterclockwise, holes
            // clockwise.
            if ((area < 0.0) == outer)
            {
                std::reverse(x.begin(), x.end());
                std::reverse(y.begin(), y.end());
            }
            outer = false;
        };
        check(f.x, f.y);
        for (Ring& hole : f.holes)
            check(hole.x, hole.y);
        // No side may meet another but its neighbours, within a ring or
        // across rings.
        std::vector<std::array<double, 4>> segments;
        std::vector<std::pair<std::size_t, std::size_t>> position; // ring, index
        std::size_t ring = 0;
        for_each_ring(f, [&](const std::vector<double>& x, const std::vector<double>& y) {
            for (std::size_t i = 0, n = x.size(); i < n; i++)
            {
                segments.push_back({x[i], y[i], x[(i + 1) % n], y[(i + 1) % n]});
                position.emplace_back(ring, i);
            }
            ring++;
        });
        std::vector<std::size_t> ring_size {f.x.size()};
        for (const Ring& hole : f.holes)
            ring_size.push_back(hole.x.size());
        for (std::size_t a = 0; a < segments.size(); a++)
            for (std::size_t c = a + 1; c < segments.size(); c++)
            {
                const auto [ra, ia] = position[a];
                const auto [rc, ic] = position[c];
                if (ra == rc && (ic == ia + 1 || (ia == 0 && ic == ring_size[ra] - 1)))
                    continue; // neighbours
                const auto& p = segments[a];
                const auto& q = segments[c];
                if (segments_meet(p[0], p[1], p[2], p[3], q[0], q[1], q[2], q[3]))
                    throw std::invalid_argument(which + (ra == rc ? ": its footprint's sides cross"
                                                                  : ": its holes and outline cross"));
            }
        for (std::size_t h = 0; h < f.holes.size(); h++)
        {
            const Ring& hole = f.holes[h];
            if (!in_ring(f.x, f.y, hole.x[0], hole.y[0]))
                throw std::invalid_argument(which + ": a hole lies outside its footprint");
            for (std::size_t k = 0; k < f.holes.size(); k++)
                if (k != h && in_ring(f.holes[k].x, f.holes[k].y, hole.x[0], hole.y[0]))
                    throw std::invalid_argument(which + ": a hole lies inside another");
        }

        // A pitched roof stands above the floor all round, and its top is its
        // highest point.
        if (!f.roof.empty())
        {
            for_each_ring(f, [&](const std::vector<double>& x, const std::vector<double>& y) {
                for (std::size_t i = 0; i < x.size(); i++)
                    if (!(roof_height(f, x[i], y[i]) > contact))
                        throw std::invalid_argument(which + ": its roof comes down to the floor");
            });
            f.top = roof_top(f);
        }
        const std::size_t n = f.x.size();
        double bx0 = f.x[0], bx1 = f.x[0], by0 = f.y[0], by1 = f.y[0];
        for (std::size_t i = 0; i < n; i++)
        {
            bx0 = std::min(bx0, f.x[i]);
            bx1 = std::max(bx1, f.x[i]);
            by0 = std::min(by0, f.y[i]);
            by1 = std::max(by1, f.y[i]);
        }
        bounds_.insert(bounds_.end(), {bx0, bx1, by0, by1});
        x_min_ = std::min(x_min_, bx0);
        x_max_ = std::max(x_max_, bx1);
        y_min_ = std::min(y_min_, by0);
        y_max_ = std::max(y_max_, by1);
        max_top_ = std::max(max_top_, f.top);
    }
    // The buildings' grid, its cells about one and a half buildings across;
    // then the junctions, which the boxes already hold.
    if (!footprints_.empty())
    {
        double building_size = 0.0;
        for (std::size_t b = 0; b < footprints_.size(); b++)
            building_size += std::max(bounds_[4 * b + 1] - bounds_[4 * b], bounds_[4 * b + 3] - bounds_[4 * b + 2]);
        building_grid_ = PlanGrid(bounds_, std::max(1.0, 1.5 * building_size / footprints_.size()));
    }
    split_at_junctions(footprints_, bounds_, building_grid_, contact);
    for (std::size_t b = 0; b < footprints_.size(); b++)
    {
        const Footprint& f = footprints_[b];
        for_each_ring(f, [&](const std::vector<double>& x, const std::vector<double>& y) {
            const std::size_t m = x.size();
            // The vertical edge at each vertex where the solid turns convex (a
            // left turn from the side before it to the side after) by more
            // than its standing within `contact` of the line between its
            // neighbours, as a vertex a junction adds does.
            std::vector<std::size_t> corner(m, none);
            for (std::size_t i = 0; i < m; i++)
            {
                const std::size_t j = (i + 1) % m, h = (i + m - 1) % m;
                const double turn = cross2(x[i] - x[h], y[i] - y[h], x[j] - x[i], y[j] - y[i]);
                const double len = std::hypot(x[j] - x[i], y[j] - y[i]);
                const double back = std::hypot(x[h] - x[i], y[h] - y[i]);
                if (turn <= contact * (len + back))
                    continue;
                const Vec3 along((x[j] - x[i]) / len, (y[j] - y[i]) / len, 0.0);
                Edge e;
                e.start = Vec3(x[i], y[i], 0.0);
                e.direction = Vec3(0.0, 0.0, 1.0);
                e.length = f.roof.empty() ? f.top : roof_height(f, x[i], y[i]);
                e.face0 = along;
                e.air = Vec3(along.y, -along.x, 0.0);
                e.n = e.angle(Vec3((x[h] - x[i]) / back, (y[h] - y[i]) / back, 0.0)) / pi;
                e.building = b;
                e.vertical = true;
                corner[i] = edges_.size();
                edges_.push_back(e);
            }
            for (std::size_t i = 0; i < m; i++)
            {
                const std::size_t j = (i + 1) % m;
                const double len = std::hypot(x[j] - x[i], y[j] - y[i]);
                const Vec3 along((x[j] - x[i]) / len, (y[j] - y[i]) / len, 0.0);
                const Vec3 outward(along.y, -along.x, 0.0);
                if (!f.roof.empty())
                {
                    // Under a pitched roof, a facade and a roof edge for each
                    // stretch of the side under one plane: an eave where the
                    // plane runs level along it, a rake where it climbs.
                    Side side {b, x[i], y[i], x[j], y[j], facades_.size(), 0, edges_.size(), 0, corner[i], corner[j]};
                    const auto stretches = roof_stretches(f, x[i], y[i], x[j], y[j]);
                    for (auto k = stretches.size(); k-- > 0;)
                    {
                        // Facades run the other way, so walk the stretches back.
                        const auto& [t0, t1, plane] = stretches[k];
                        const double xa = x[i] + t0 * (x[j] - x[i]), ya = y[i] + t0 * (y[j] - y[i]);
                        const double xb = x[i] + t1 * (x[j] - x[i]), yb = y[i] + t1 * (y[j] - y[i]);
                        Facade facade;
                        facade.x0 = xb;
                        facade.y0 = yb;
                        facade.x1 = xa;
                        facade.y1 = ya;
                        facade.bottom = 0.0;
                        facade.top_start = f.roof[plane].at(xb, yb);
                        facade.top_end = f.roof[plane].at(xa, ya);
                        facade.top = std::max(facade.top_start, facade.top_end);
                        facades_.push_back(facade);
                        facade_building_.push_back(b);
                    }
                    for (const auto& [t0, t1, plane] : stretches)
                    {
                        const double xa = x[i] + t0 * (x[j] - x[i]), ya = y[i] + t0 * (y[j] - y[i]);
                        const double xb = x[i] + t1 * (x[j] - x[i]), yb = y[i] + t1 * (y[j] - y[i]);
                        edges_.push_back(roof_edge(Vec3(xa, ya, f.roof[plane].at(xa, ya)),
                                Vec3(xb, yb, f.roof[plane].at(xb, yb)), f.roof[plane], -1.0 * outward,
                                Vec3(0.0, 0.0, -1.0), b));
                    }
                    side.facade_end = facades_.size();
                    side.roof_end = edges_.size();
                    sides.push_back(side);
                    continue;
                }
                // The side i -> j, the solid on its left; its facade runs the
                // other way, so that the outside is on the facade's left.
                Facade facade;
                facade.x0 = x[j];
                facade.y0 = y[j];
                facade.x1 = x[i];
                facade.y1 = y[i];
                facade.bottom = 0.0;
                facade.top = f.top;
                // The roof edge: the roof leaves it inward, the air is above.
                Edge roof;
                roof.start = Vec3(x[i], y[i], f.top);
                roof.direction = along;
                roof.length = len;
                roof.face0 = -1.0 * outward;
                roof.air = Vec3(0.0, 0.0, 1.0);
                roof.n = 1.5;
                roof.building = b;
                sides.push_back({b, x[i], y[i], x[j], y[j], facades_.size(), facades_.size() + 1, edges_.size(),
                        edges_.size() + 1, corner[i], corner[j]});
                facades_.push_back(facade);
                facade_building_.push_back(b);
                edges_.push_back(roof);
            }
        });
        // A pitched roof's ridges and hips, each between its two planes.
        if (!f.roof.empty())
            for (const Ridge& r : roof_ridges(f))
            {
                // Each plane leaves the ridge toward the side where it is the
                // lower: the gradient of the other less it, in plan.
                const RoofPlane& pk = f.roof[r.k];
                const RoofPlane& pl = f.roof[r.l];
                const Vec3 toward_k(pl.a - pk.a, pl.b - pk.b, 0.0), toward_l = -1.0 * toward_k;
                const Edge probe = roof_edge(r.a, r.b, pl, toward_l, Vec3(0.0, 0.0, -1.0), b);
                const Vec3 face_l = probe.face0;
                edges_.push_back(roof_edge(r.a, r.b, pk, toward_k, face_l, b));
            }
    }

    // Sides two buildings share exactly: a side's neighbour runs the other
    // way between the same two points.
    auto key = [](double ax, double ay, double bx, double by) {
        auto r = [](double v) { return static_cast<long long>(std::llround(v * 1e6)); };
        return std::array<long long, 4> {r(ax), r(ay), r(bx), r(by)};
    };
    std::map<std::array<long long, 4>, std::size_t> by_ends;
    for (std::size_t k = 0; k < sides.size(); k++)
        by_ends.emplace(key(sides[k].x0, sides[k].y0, sides[k].x1, sides[k].y1), k);
    std::vector<char> drop_facade(facades_.size(), 0), drop_edge(edges_.size(), 0);
    for (const Side& side : sides)
    {
        const auto other = by_ends.find(key(side.x1, side.y1, side.x0, side.y0));
        if (other == by_ends.end() || sides[other->second].building == side.building)
            continue;
        shared_sides_++;
        // Only what stands above the neighbour's roof is outside, judged along
        // the side: both roofs' heights over it, t from (x0, y0) to (x1, y1),
        // straight between the ends of either's stretches and where they cross.
        const Footprint& mine = footprints_[side.building];
        const Footprint& theirs = footprints_[sides[other->second].building];
        const double dx = side.x1 - side.x0, dy = side.y1 - side.y0;
        const auto height = [&](const Footprint& f, double t) {
            return roof_height(f, side.x0 + t * dx, side.y0 + t * dy);
        };
        const auto above = [&](double t) { return height(mine, t) - height(theirs, t); };
        const double level_tolerance = 1e-9;
        std::vector<double> cuts {0.0, 1.0};
        for (const Footprint* f : {&mine, &theirs})
            if (!f->roof.empty())
                for (const auto& [t0, t1, plane] : roof_stretches(*f, side.x0, side.y0, side.x1, side.y1))
                    cuts.insert(cuts.end(), {t0, t1});
        std::sort(cuts.begin(), cuts.end());
        for (std::size_t k = 0, n = cuts.size(); k + 1 < n; k++)
        {
            const double d0 = above(cuts[k]), d1 = above(cuts[k + 1]);
            if ((d0 > level_tolerance && d1 < -level_tolerance) || (d0 < -level_tolerance && d1 > level_tolerance))
                cuts.push_back(cuts[k] + (cuts[k + 1] - cuts[k]) * d0 / (d0 - d1));
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
        // The stretches where this roof stands above the neighbour's, merged.
        std::vector<std::pair<double, double>> exposed;
        bool level = true;
        for (std::size_t k = 0; k + 1 < cuts.size(); k++)
        {
            const double mid = above(0.5 * (cuts[k] + cuts[k + 1]));
            level = level && std::abs(mid) <= level_tolerance && std::abs(above(cuts[k])) <= level_tolerance;
            if (!(mid > level_tolerance))
                continue;
            if (!exposed.empty() && exposed.back().second == cuts[k])
                exposed.back().second = cuts[k + 1];
            else
                exposed.emplace_back(cuts[k], cuts[k + 1]);
        }
        // Two roofs level across the side are one roof there: the later
        // building cedes the side to the earlier, so a bounce on it is the
        // earlier's alone (on_roof), and neither keeps a wall or an edge on it.
        if (level && sides[other->second].building < side.building)
            ceded_[side.building].push_back({side.x0, side.y0, side.x1, side.y1});
        // A corner stands only above the neighbour's roof at its vertex.
        for (const auto& [c, t] : {std::pair {side.corner0, 0.0}, std::pair {side.corner1, 1.0}})
            if (c != none)
            {
                const double floor = height(theirs, t);
                if (!(above(t) > level_tolerance))
                    drop_edge[c] = 1;
                else if (edges_[c].start.z < floor)
                {
                    const double end = edges_[c].start.z + edges_[c].length;
                    edges_[c].start.z = floor;
                    edges_[c].length = end - floor;
                }
            }
        // The side's walls stand down to the neighbour's roof where it is below
        // theirs, a piece between each two cuts so that both its top and its
        // bottom run straight; its roof edges stand where its roof stands
        // above, whole across a point where the roofs only touch.  A side
        // above all along with no cut inside a wall keeps its walls and edges,
        // their bottoms raised.
        const std::size_t walls = side.facade_end - side.facade;
        const auto along = [&](double x, double y) {
            return ((x - side.x0) * dx + (y - side.y0) * dy) / (dx * dx + dy * dy);
        };
        bool whole = exposed.size() == 1 && exposed[0].first == 0.0 && exposed[0].second == 1.0;
        for (std::size_t k = side.facade; k < side.facade_end && whole; k++)
        {
            const double t0 = along(facades_[k].x1, facades_[k].y1), t1 = along(facades_[k].x0, facades_[k].y0);
            for (double c : cuts)
                whole = whole && !(c > t0 + 1e-12 && c < t1 - 1e-12);
        }
        for (std::size_t k = side.facade; k < side.facade_end; k++)
            drop_facade[k] = !whole;
        for (std::size_t k = side.roof; k < side.roof_end; k++)
            drop_edge[k] = !whole;
        for (std::size_t w = 0; w < walls; w++)
        {
            // Walls run against the side, so the side's first stretch is its
            // last wall; and each has its stretch of the side.
            const std::size_t k = side.facade + walls - 1 - w;
            const Facade original = facades_[k];
            const double t0 = along(original.x1, original.y1), t1 = along(original.x0, original.y0);
            if (!whole)
            {
                const Edge edge = edges_[side.roof + w];
                for (const auto& [u0, u1] : exposed)
                {
                    const double a = std::max(u0, t0), b = std::min(u1, t1);
                    if (!(b > a))
                        continue;
                    const Vec3 from(side.x0 + a * dx, side.y0 + a * dy, height(mine, a));
                    const Vec3 to(side.x0 + b * dx, side.y0 + b * dy, height(mine, b));
                    Edge part = edge;
                    part.start = from;
                    part.length = (to - from).norm();
                    edges_.push_back(part);
                    drop_edge.push_back(0);
                }
            }
            for (std::size_t m = 0; m + 1 < cuts.size(); m++)
            {
                const double a = std::max(cuts[m], t0), b = std::min(cuts[m + 1], t1);
                if (!(b > a) || !(above(0.5 * (a + b)) > level_tolerance))
                    continue;
                Facade piece = original;
                piece.x0 = side.x0 + b * dx;
                piece.y0 = side.y0 + b * dy;
                piece.x1 = side.x0 + a * dx;
                piece.y1 = side.y0 + a * dy;
                if (original.top_start == original.top_start)
                {
                    piece.top_start = height(mine, b);
                    piece.top_end = height(mine, a);
                    piece.top = std::max(piece.top_start, piece.top_end);
                }
                piece.bottom_start = height(theirs, b);
                piece.bottom_end = height(theirs, a);
                piece.bottom = std::min(piece.bottom_start, piece.bottom_end);
                if (piece.bottom_start == piece.bottom_end)
                    piece.bottom_start = piece.bottom_end = std::numeric_limits<double>::quiet_NaN();
                if (whole)
                {
                    facades_[k] = piece;
                    continue;
                }
                facades_.push_back(piece);
                facade_building_.push_back(side.building);
                drop_facade.push_back(0);
            }
        }
    }
    if (shared_sides_ > 0)
    {
        std::vector<Facade> kept_facades;
        std::vector<std::size_t> kept_building;
        for (std::size_t i = 0; i < facades_.size(); i++)
            if (!drop_facade[i])
            {
                kept_facades.push_back(facades_[i]);
                kept_building.push_back(facade_building_[i]);
            }
        std::vector<Edge> kept_edges;
        for (std::size_t i = 0; i < edges_.size(); i++)
            if (!drop_edge[i])
                kept_edges.push_back(edges_[i]);
        facades_ = std::move(kept_facades);
        facade_building_ = std::move(kept_building);
        edges_ = std::move(kept_edges);
    }
    for (Edge& e : edges_)
        e.air1 = std::sin(e.n * pi) * e.face0 - std::cos(e.n * pi) * e.air;
    roof_facade_.assign(footprints_.size(), no_facade);
    for (std::size_t f = facades_.size(); f-- > 0;)
        roof_facade_[facade_building_[f]] = f;

    // The grids, their cells about one and a half typical items across.
    std::vector<double> facade_boxes;
    double facade_size = 0.0;
    for (const Facade& f : facades_)
    {
        const double length = std::hypot(f.x1 - f.x0, f.y1 - f.y0);
        const double ax = (f.x1 - f.x0) / length, ay = (f.y1 - f.y0) / length;
        frames_.push_back({f.x0, f.y0, ax, ay, length, -ay, ax});
        const FacadeFrame& fi = frames_.back();
        const auto same = std::find_if(orientations_.begin(), orientations_.end(), [&](std::size_t o) {
            return std::abs(frames_[o].nx * fi.ny - frames_[o].ny * fi.nx) < std::sin(orientation_tolerance_);
        });
        if (same == orientations_.end())
        {
            orientations_.push_back(frames_.size() - 1);
            orientation_of_.push_back(frames_.size() - 1);
        }
        else
            orientation_of_.push_back(*same);
        facade_boxes.insert(facade_boxes.end(),
                {std::min(f.x0, f.x1), std::max(f.x0, f.x1), std::min(f.y0, f.y1), std::max(f.y0, f.y1)});
        facade_size += length;
    }
    if (!facades_.empty())
    {
        facade_grid_ = PlanGrid(facade_boxes, std::max(1.0, 1.5 * facade_size / facades_.size()));
        for (std::size_t i = 0; i < facades_.size(); i++)
            facade_cells_.push_back(facade_grid_.item_cells(facade_boxes[4 * i], facade_boxes[4 * i + 1],
                    facade_boxes[4 * i + 2], facade_boxes[4 * i + 3]));
    }

    // Each edge's buriers: the other buildings whose boxes reach its own.
    edge_buriers_.assign(edges_.size(), {});
    std::vector<std::uint32_t> candidates;
    for (std::size_t i = 0; i < edges_.size(); i++)
    {
        const Edge& e = edges_[i];
        const Vec3 end = e.start + e.length * e.direction;
        const double x0 = std::min(e.start.x, end.x) - contact, x1 = std::max(e.start.x, end.x) + contact;
        const double y0 = std::min(e.start.y, end.y) - contact, y1 = std::max(e.start.y, end.y) + contact;
        near_box(x0, x1, y0, y1, candidates);
        for (const std::uint32_t b : candidates)
        {
            const double* box = &bounds_[4 * b];
            if (b != e.building && box[0] <= x1 && box[1] >= x0 && box[2] <= y1 && box[3] >= y0)
                edge_buriers_[i].push_back(b);
        }
        std::sort(edge_buriers_[i].begin(), edge_buriers_[i].end());
    }
    // The edges were made building by building, and dropping the shared ones
    // kept their order.
    edge_begin_.assign(footprints_.size() + 1, edges_.size());
    for (std::size_t i = edges_.size(); i-- > 0;)
        edge_begin_[edges_[i].building] = i;
    for (std::size_t b = footprints_.size(); b-- > 0;)
        edge_begin_[b] = std::min(edge_begin_[b], edge_begin_[b + 1]);
    for (std::size_t i = 0; i + 1 < edges_.size(); i++)
        assert(edges_[i].building <= edges_[i + 1].building);
    continue_edges();
    convex_pieces_.reserve(footprints_.size());
    for (const Footprint& f : footprints_)
        convex_pieces_.push_back(convex_pieces_of(f));
    prism_begin_.reserve(footprints_.size());
    for (std::size_t b = 0; b < footprints_.size(); b++)
    {
        prism_begin_.push_back(prisms_.size() / 2);
        const std::vector<Footprint>& pieces = convex_pieces_[b];
        const Footprint* first = pieces.empty() ? &footprints_[b] : pieces.data();
        for (const Footprint* f = first; f != first + std::max<std::size_t>(pieces.size(), 1); f++)
            for (const bool below_floor : {false, true})
            {
                prisms_.emplace_back();
                if (!prism_planes(*f, below_floor ? -f->top : 0.0, nullptr, prisms_.back()))
                    prisms_.back().clear();
            }
    }
    front_ = std::make_shared<FrontTable>(edges_.size(), front_words());
}

const std::uint64_t* Buildings::front_facades(std::size_t edge) const
{
    FrontTable& table = *front_;
    std::uint64_t* row = table.bits.get() + edge * table.words;
    // call_once's completion happens before its later callers return, so
    // every thread reads the row whole.
    std::call_once(table.formed[edge], [&] {
        std::fill(row, row + table.words, 0);
        const Edge& g = edges_[edge];
        const Vec3 a = g.start, b = g.start + g.length * g.direction;
        for (std::size_t f = 0; f < frames_.size(); f++)
            if (frames_[f].side(a) > contact || frames_[f].side(b) > contact)
                row[f / 64] |= std::uint64_t(1) << (f % 64);
    });
    return row;
}

void Buildings::continue_edges()
{
    // Ends within this of each other, or of a plane, meet it (m): a
    // junction's split puts a vertex within `contact` of the line.
    const double near = 10.0 * contact;
    const auto same = [](const Vec3& u, const Vec3& v) { return (u - v).norm() < 1e-6; };
    const std::size_t none = static_cast<std::size_t>(-1);
    // The piece that goes on from each edge's far end: an edge of the same
    // wedge on the same line, starting there.
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> starts;
    const auto cell = [&](const Vec3& p, int dx, int dy, int dz) {
        const auto q = [&](double v, int d) { return static_cast<std::uint64_t>(std::llround(std::floor(v / near)) + d); };
        return (q(p.x, dx) * 73856093u) ^ (q(p.y, dy) * 19349663u) ^ (q(p.z, dz) * 83492791u);
    };
    for (std::size_t i = 0; i < edges_.size(); i++)
        starts[cell(edges_[i].start, 0, 0, 0)].push_back(i);
    std::vector<std::size_t> next(edges_.size(), none), previous(edges_.size(), none);
    for (std::size_t i = 0; i < edges_.size(); i++)
    {
        const Edge& e = edges_[i];
        const Vec3 end = e.start + e.length * e.direction;
        for (int dx = -1; dx <= 1 && next[i] == none; dx++)
            for (int dy = -1; dy <= 1 && next[i] == none; dy++)
                for (int dz = -1; dz <= 1 && next[i] == none; dz++)
                {
                    const auto found = starts.find(cell(end, dx, dy, dz));
                    if (found == starts.end())
                        continue;
                    for (std::size_t j : found->second)
                    {
                        const Edge& f = edges_[j];
                        if (j != i && previous[j] == none && (f.start - end).norm() < near
                                && same(f.direction, e.direction) && same(f.face0, e.face0) && same(f.air, e.air)
                                && std::abs(f.n - e.n) < 1e-6)
                        {
                            next[i] = j;
                            previous[j] = i;
                            break;
                        }
                    }
                }
    }
    // Whether an edge's end at p, its line running on along `out`, meets a
    // rigid plane square to it that its image is whole in: the floor, or walls
    // in one plane that cover the wedge's air about p -- a taller neighbour
    // whose wall the eave runs into, not one that stands only behind the
    // eave's own facade, as in a row whose fronts are flush, where the end is
    // a vertex.
    std::vector<std::uint32_t> walls;
    const auto mirrors = [&](const Edge& e, const Vec3& p, const Vec3& out) {
        if (std::abs(std::abs(out.z) - 1.0) < 1e-12)
            return out.z < 0.0 && std::abs(p.z) < 1e-9;
        if (std::abs(out.z) > 1e-9)
            return false;
        // Whether q lies on a wall in p's plane, square to the line.
        const auto on_wall = [&](const Vec3& q) {
            facades_in_box(q.x - near, q.x + near, q.y - near, q.y + near, walls);
            for (std::uint32_t w : walls)
            {
                const Facade& f = facades_[w];
                const FacadeFrame& a = frames_[w];
                const double along = (q.x - a.x0) * a.ax + (q.y - a.y0) * a.ay;
                if (out.x * a.nx + out.y * a.ny < -(1.0 - 1e-9) && std::abs(facade_side(f, p)) < near
                        && std::abs(facade_side(f, q)) < near && along > -near && along < a.length + near
                        && q.z > facade_bottom(f, q) - near && q.z < facade_top(f, q) + near)
                    return true;
            }
            return false;
        };
        // The wedge's air about p, face to face, a little way out.
        const double r = 0.1 * std::min(1.0, e.length);
        for (int j = 0; j <= 8; j++)
        {
            const double phi = e.n * pi * (0.02 + 0.96 * j / 8.0);
            if (!on_wall(p + r * (std::cos(phi) * e.face0 + std::sin(phi) * e.air)))
                return false;
        }
        return on_wall(p);
    };
    // Each run of pieces from its first: its length and its ends' planes.
    for (std::size_t first = 0; first < edges_.size(); first++)
    {
        if (previous[first] != none)
            continue;
        std::vector<std::size_t> run {first};
        while (next[run.back()] != none && run.size() <= edges_.size())
            run.push_back(next[run.back()]);
        const Edge& head = edges_[first];
        const Edge& tail = edges_[run.back()];
        double length = 0.0;
        for (std::size_t i : run)
            length += edges_[i].length;
        const bool lo = mirrors(head, head.start, -1.0 * head.direction);
        const bool hi = mirrors(tail, tail.start + tail.length * tail.direction, tail.direction);
        const double inf = std::numeric_limits<double>::infinity();
        double before = 0.0;
        for (std::size_t i : run)
        {
            Edge& e = edges_[i];
            const double after = length - before - e.length;
            e.beyond[0] = lo && hi ? inf : before + (lo ? length : 0.0);
            e.beyond[1] = lo && hi ? inf : after + (hi ? length : 0.0);
            before += e.length;
        }
    }
}

bool Buildings::strictly_inside(std::size_t building, double x, double y) const
{
    // The crossings first: most points tested lie outside, and only a point
    // inside needs its distance from every side.
    const Footprint& f = footprints_[building];
    bool in = false;
    for_each_ring(f, [&](const std::vector<double>& rx, const std::vector<double>& ry) {
        for (std::size_t i = 0, n = rx.size(), h = n - 1; i < n; h = i++)
            if ((ry[i] > y) != (ry[h] > y) && x < rx[h] + (y - ry[h]) * (rx[i] - rx[h]) / (ry[i] - ry[h]))
                in = !in;
    });
    if (!in)
        return false;
    bool near = false;
    for_each_ring(f, [&](const std::vector<double>& rx, const std::vector<double>& ry) {
        for (std::size_t i = 0, n = rx.size(), h = n - 1; i < n && !near; h = i++)
            near = distance_to_side(x, y, rx[h], ry[h], rx[i], ry[i]) < contact;
    });
    return !near;
}

bool Buildings::on_roof(std::size_t building, double x, double y) const
{
    // A pitched roof is a level plane to bounce on only where its least plane
    // is level: a mansard's deck, at its top.
    const Footprint& pitched = footprints_[building];
    if (!pitched.roof.empty())
    {
        const RoofPlane& p = pitched.roof[active_plane(pitched, x, y)];
        if (p.a != 0.0 || p.b != 0.0)
            return false;
    }
    return over(building, x, y);
}

long Buildings::roof_surface(std::size_t building, double x, double y) const
{
    if (!over(building, x, y))
        return roof_none;
    const Footprint& f = footprints_[building];
    if (f.roof.empty())
        return roof_level;
    const std::size_t plane = active_plane(f, x, y);
    return f.roof[plane].level() ? roof_level : static_cast<long>(plane);
}

bool Buildings::on_plane(std::size_t building, std::size_t plane, double x, double y) const
{
    const Footprint& f = footprints_[building];
    return plane < f.roof.size() && f.roof[plane].at(x, y) <= roof_height(f, x, y) + 1e-9 && over(building, x, y);
}

bool Buildings::over(std::size_t building, double x, double y) const
{
    for (const auto& s : ceded_[building])
        if (distance_to_side(x, y, s[0], s[1], s[2], s[3]) < 1e-9)
            return false;
    const Footprint& f = footprints_[building];
    bool in = false, side = false;
    for_each_ring(f, [&](const std::vector<double>& rx, const std::vector<double>& ry) {
        for (std::size_t i = 0, n = rx.size(), h = n - 1; i < n; h = i++)
        {
            if ((ry[i] > y) != (ry[h] > y) && x < rx[h] + (y - ry[h]) * (rx[i] - rx[h]) / (ry[i] - ry[h]))
                in = !in;
            side = side || distance_to_side(x, y, rx[h], ry[h], rx[i], ry[i]) < 1e-9;
        }
    });
    return in || side;
}

std::size_t Buildings::inside(double x, double y) const
{
    for (std::size_t b = 0; b < footprints_.size(); b++)
        if (x > bounds_[4 * b] && x < bounds_[4 * b + 1] && y > bounds_[4 * b + 2] && y < bounds_[4 * b + 3]
                && strictly_inside(b, x, y))
            return b;
    return footprints_.size();
}

bool Buildings::buried_in(std::size_t building, const Vec3& p) const
{
    return p.z > contact && p.z < roof_height(footprints_[building], p.x, p.y) - contact
            && strictly_inside(building, p.x, p.y);
}

bool Buildings::buried(const Vec3& p) const
{
    bool found = false;
    building_grid_.walk(p.x, p.y, 0.0, 0.0, 0.0, [&](const std::uint32_t* first, const std::uint32_t* last, double) {
        for (const std::uint32_t* i = first; i != last && !found; i++)
            found = buried_in(*i, p);
        return false;
    });
    return found;
}

bool Buildings::buried_on(std::size_t edge, const Vec3& p) const
{
    for (const std::uint32_t b : edge_buriers_[edge])
        if (buried_in(b, p))
            return true;
    return false;
}

Buildings::Support Buildings::support(const Vec3& p) const
{
    Support out;
    std::size_t surfaces = 0;
    thread_local std::vector<std::uint32_t> near;
    // Walls: within contact of the plane, along the base and up to the top.
    facades_in_box(p.x - contact, p.x + contact, p.y - contact, p.y + contact, near);
    for (const std::uint32_t f : near)
    {
        const Facade& w = facades_[f];
        const double dx = w.x1 - w.x0, dy = w.y1 - w.y0, length = std::hypot(dx, dy);
        const double along = ((p.x - w.x0) * dx + (p.y - w.y0) * dy) / length;
        if (std::abs(facade_side(w, p)) <= contact && along >= -contact && along <= length + contact
                && p.z >= facade_bottom(w, p) - contact && p.z <= facade_top(w, p) + contact)
        {
            out = {Support::Kind::facade, f};
            surfaces++;
        }
    }
    // Roofs: over the footprint (or on its outline) at the roof's height.  Two
    // level roofs across a shared side are one roof there.
    near_box(p.x - contact, p.x + contact, p.y - contact, p.y + contact, near);
    std::sort(near.begin(), near.end());
    bool roofed = false;
    for (const std::uint32_t b : near)
    {
        const Footprint& f = footprints_[b];
        if (std::abs(p.z - roof_height(f, p.x, p.y)) > contact)
            continue;
        bool on = in_footprint(f, p.x, p.y);
        for_each_ring(f, [&](const std::vector<double>& x, const std::vector<double>& y) {
            for (std::size_t i = 0, n = x.size(), h = n - 1; i < n && !on; h = i++)
                on = distance_to_side(p.x, p.y, x[h], y[h], x[i], y[i]) <= contact;
        });
        if (!on || roofed)
            continue;
        roofed = true;
        surfaces++;
        const RoofPlane* plane = f.roof.empty() ? nullptr : &f.roof[active_plane(f, p.x, p.y)];
        const bool level = plane == nullptr || (plane->a == 0.0 && plane->b == 0.0);
        if (out.kind == Support::Kind::none)
            out = {level ? Support::Kind::roof : Support::Kind::slope, b};
    }
    // An edge, or more than one surface.
    for (const std::uint32_t b : near)
        for (std::size_t e = edge_begin(b); e < edge_begin(b + 1) && surfaces < 2; e++)
        {
            const Edge& g = edges_[e];
            const double t = std::clamp((p - g.start).dot(g.direction), 0.0, g.length);
            if ((g.start + t * g.direction - p).norm() <= contact)
                surfaces = 2;
        }
    if (surfaces >= 2)
        out = {Support::Kind::edge, 0};
    return out;
}

bool Buildings::blocked_by(std::size_t building, const Vec3& a, const Vec3& b, bool below_floor) const
{
    const Footprint& f = footprints_[building];
    const double z_low = below_floor ? -f.top : 0.0;
    const double z_high = f.top;
    // Heights strictly inside, by more than round-off.
    const double z_margin = 1e-9;
    if (std::min(a.z, b.z) >= z_high - z_margin || std::max(a.z, b.z) <= z_low + z_margin)
        return false;
    const double* box = &bounds_[4 * building];
    if (std::max(a.x, b.x) <= box[0] || std::min(a.x, b.x) >= box[1] || std::max(a.y, b.y) <= box[2]
            || std::min(a.y, b.y) >= box[3])
        return false;
    const Vec3 d = b - a;
    const double length = d.norm();
    if (length < contact)
        return false;

    // Where the segment's status can change: its ends, the two heights, and
    // each side of the footprint it meets.  Between consecutive breakpoints it
    // is inside throughout or nowhere, so a midpoint decides each stretch.
    thread_local std::vector<double> t_values;
    t_values.assign({0.0, 1.0});
    auto add = [&](double t) {
        if (t > 0.0 && t < 1.0)
            t_values.push_back(t);
    };
    if (d.z != 0.0)
    {
        add((z_high - a.z) / d.z);
        add((z_low - a.z) / d.z);
    }
    // Under a pitched roof, where the segment crosses each plane (and below
    // the floor its image).
    for (const RoofPlane& plane : f.roof)
    {
        const double rate = d.z - plane.a * d.x - plane.b * d.y;
        if (rate != 0.0)
            add((plane.at(a.x, a.y) - a.z) / rate);
        if (below_floor)
        {
            const double image = d.z + plane.a * d.x + plane.b * d.y;
            if (image != 0.0)
                add((-plane.at(a.x, a.y) - a.z) / image);
        }
    }
    for_each_ring(f, [&](const std::vector<double>& rx, const std::vector<double>& ry) {
        for (std::size_t i = 0, n = rx.size(), h = n - 1; i < n; h = i++)
        {
            const double ex = rx[i] - rx[h], ey = ry[i] - ry[h];
            const double denominator = cross2(d.x, d.y, ex, ey);
            if (denominator == 0.0)
                continue;
            const double px = rx[h] - a.x, py = ry[h] - a.y;
            const double t = cross2(px, py, ex, ey) / denominator;
            const double s = cross2(px, py, d.x, d.y) / denominator;
            if (s >= -1e-12 && s <= 1.0 + 1e-12)
                add(t);
        }
    });
    std::sort(t_values.begin(), t_values.end());
    for (std::size_t k = 0; k + 1 < t_values.size(); k++)
    {
        const double t0 = t_values[k], t1 = t_values[k + 1];
        if ((t1 - t0) * length <= contact)
            continue;
        const double t = 0.5 * (t0 + t1);
        const double z = a.z + t * d.z;
        if (!(z > z_low + z_margin && z < z_high - z_margin))
            continue;
        const double x = a.x + t * d.x, y = a.y + t * d.y;
        if (!f.roof.empty())
        {
            const double h = roof_height(f, x, y);
            if (!(z > (below_floor ? -h : 0.0) + z_margin && z < h - z_margin))
                continue;
        }
        if (strictly_inside(building, x, y))
            return true;
    }
    return false;
}

namespace {

// A stamp per item, so that a query walking several cells tests each item
// once: an item is new to the query when its stamp is not the query's.
struct Stamps
{
    std::vector<std::uint32_t> stamp;
    std::uint32_t current = 0;

    void next(std::size_t items)
    {
        if (stamp.size() < items)
            stamp.resize(items, 0);
        if (++current == 0)
        {
            std::fill(stamp.begin(), stamp.end(), 0);
            current = 1;
        }
    }
    bool first(std::uint32_t item)
    {
        if (stamp[item] == current)
            return false;
        stamp[item] = current;
        return true;
    }
};

} // namespace

bool Buildings::blocked(const Vec3& a, const Vec3& b, bool below_floor) const
{
    thread_local Stamps seen;
    seen.next(footprints_.size());
    bool hit = false;
    building_grid_.walk(a.x, a.y, b.x - a.x, b.y - a.y, 1.0,
            [&](const std::uint32_t* first, const std::uint32_t* last, double) {
                for (const std::uint32_t* i = first; i != last; i++)
                    if (seen.first(*i) && blocked_by(*i, a, b, below_floor))
                    {
                        hit = true;
                        return false;
                    }
                return true;
            });
    return hit;
}

void Buildings::near_box(double x_min, double x_max, double y_min, double y_max, std::vector<std::uint32_t>& out) const
{
    thread_local Stamps seen;
    seen.next(footprints_.size());
    out.clear();
    building_grid_.box(x_min, x_max, y_min, y_max, [&](const std::uint32_t* first, const std::uint32_t* last) {
        for (const std::uint32_t* i = first; i != last; i++)
            if (seen.first(*i))
                out.push_back(*i);
    });
}

void Buildings::section(const Vec3& a, const Vec3& b, std::vector<SectionSpan>& out) const
{
    out.clear();
    const double dx = b.x - a.x, dy = b.y - a.y, span = std::hypot(dx, dy);
    if (!(span > 0.0))
        return;
    thread_local Stamps seen;
    seen.next(footprints_.size());
    thread_local std::vector<double> cuts;
    building_grid_.walk(a.x, a.y, dx, dy, 1.0, [&](const std::uint32_t* first, const std::uint32_t* last, double) {
        for (const std::uint32_t* i = first; i != last; i++)
        {
            if (!seen.first(*i))
                continue;
            const Footprint& f = footprints_[*i];
            // Where the segment crosses the rings, as fractions of it; then
            // each stretch between two crossings whose middle is inside.
            cuts.assign({0.0, 1.0});
            for_each_ring(f, [&](const std::vector<double>& rx, const std::vector<double>& ry) {
                for (std::size_t k = 0, n = rx.size(); k < n; k++)
                {
                    const std::size_t l = (k + 1) % n;
                    const double ex = rx[l] - rx[k], ey = ry[l] - ry[k], wx = rx[k] - a.x, wy = ry[k] - a.y;
                    const double det = ex * dy - dx * ey;
                    if (det == 0.0)
                        continue;
                    const double t = (ex * wy - wx * ey) / det, s = (dx * wy - dy * wx) / det;
                    if (t > 0.0 && t < 1.0 && s >= 0.0 && s <= 1.0)
                        cuts.push_back(t);
                }
            });
            std::sort(cuts.begin(), cuts.end());
            for (std::size_t k = 0; k + 1 < cuts.size(); k++)
            {
                if (!(cuts[k + 1] > cuts[k]))
                    continue;
                const double t = 0.5 * (cuts[k] + cuts[k + 1]);
                if (!in_footprint(f, a.x + t * dx, a.y + t * dy))
                    continue;
                if (!out.empty() && out.back().building == *i && out.back().d1 == cuts[k] * span)
                    out.back().d1 = cuts[k + 1] * span;
                else
                    out.push_back({cuts[k] * span, cuts[k + 1] * span, *i});
            }
        }
        return true;
    });
}

void Buildings::facades_in_box(double x_min, double x_max, double y_min, double y_max,
        std::vector<std::uint32_t>& out) const
{
    out.clear();
    PlanGrid::Cells cells;
    const FacadesInBox which = facades_in_box(x_min, x_max, y_min, y_max, cells);
    if (which == FacadesInBox::none)
        return;
    if (which == FacadesInBox::all)
    {
        out.resize(facades_.size());
        for (std::size_t i = 0; i < out.size(); i++)
            out[i] = static_cast<std::uint32_t>(i);
        return;
    }
    thread_local Stamps seen;
    seen.next(facades_.size());
    facade_grid_.box(x_min, x_max, y_min, y_max, [&](const std::uint32_t* first, const std::uint32_t* last) {
        for (const std::uint32_t* i = first; i != last; i++)
            if (seen.first(*i))
                out.push_back(*i);
    });
    std::sort(out.begin(), out.end());
}

Buildings::FacadesInBox Buildings::facades_in_box(double x_min, double x_max, double y_min, double y_max,
        PlanGrid::Cells& cells) const
{
    if (facades_.empty())
        return FacadesInBox::none;
    // Between-edge reflection searches often cover the whole district.
    // Every facade is then a candidate, already in ascending order: avoid
    // walking its grid cells, deduplicating their entries and sorting them.
    if (x_min <= x_min_ && x_max >= x_max_ && y_min <= y_min_ && y_max >= y_max_)
        return FacadesInBox::all;
    return facade_grid_.box_cells(x_min, x_max, y_min, y_max, cells) ? FacadesInBox::cells : FacadesInBox::none;
}

std::size_t Buildings::first_hit(const Vec3& pos, const Vec3& dir, bool below_floor, Vec3& hit) const
{
    std::size_t best = facades_.size();
    double best_t = std::numeric_limits<double>::infinity();
    // Above every roof nothing is in the way.
    double t_max = std::numeric_limits<double>::infinity();
    if (dir.z > 0.0)
        t_max = std::max(0.0, (max_top_ - pos.z) / dir.z) + 1.0;
    thread_local Stamps seen;
    seen.next(facades_.size());
    facade_grid_.walk(pos.x, pos.y, dir.x, dir.y, t_max,
            [&](const std::uint32_t* first, const std::uint32_t* last, double t_exit) {
                for (const std::uint32_t* item = first; item != last; item++)
                {
                    const std::uint32_t i = *item;
                    if (!seen.first(i))
                        continue;
                    const FacadeFrame& fr = frames_[i];
                    const double side = fr.side(pos);
                    const double approach = dir.x * fr.nx + dir.y * fr.ny;
                    if (side <= contact || approach >= 0.0)
                        continue;
                    const double t = side / -approach;
                    if (t >= best_t)
                        continue;
                    const Vec3 p = pos + t * dir;
                    const double along = (p.x - fr.x0) * fr.ax + (p.y - fr.y0) * fr.ay;
                    // The wall, and below the floor its image too.
                    const double top = facade_top(facades_[i], p), bottom = facade_bottom(facades_[i], p);
                    const double z = below_floor && p.z < 0.0 ? -p.z : p.z;
                    if (along >= -1e-9 && along <= fr.length + 1e-9 && z >= bottom - 1e-9 && z <= top + 1e-9)
                    {
                        best = i;
                        best_t = t;
                        hit = p;
                    }
                }
                // Every cell nearer than the best hit has been searched.
                return best_t > t_exit;
            });
    return best;
}

bool Buildings::escapes(const Vec3& a, const Vec3& direction, bool below_floor) const
{
    if (footprints_.empty())
        return true;
    const double inf = std::numeric_limits<double>::infinity();
    // Past either of these nothing can block: above every roof, or outside
    // every footprint and moving away.
    double t_rise = inf;
    if (direction.z > 0.0)
        t_rise = std::max(0.0, (max_top_ - a.z) / direction.z);
    auto exit_time = [&](double start, double speed, double low, double high) {
        if (speed > 0.0)
            return (high - start) / speed;
        if (speed < 0.0)
            return (low - start) / speed;
        return (start <= low || start >= high) ? -inf : inf;
    };
    const double t_box = std::max(0.0,
            std::min(exit_time(a.x, direction.x, x_min_ - 1.0, x_max_ + 1.0),
                    exit_time(a.y, direction.y, y_min_ - 1.0, y_max_ + 1.0)));
    const double t_far = std::min(t_rise, t_box);
    if (!std::isfinite(t_far))
        return true; // straight up or down through nothing finite: no footprint in reach
    return !blocked(a, a + (t_far + 1.0) * direction, below_floor);
}

int edge_face(const Edge& e, const Vec3& p)
{
    const Vec3 v = e.project(p).normal;
    const double c = std::cos(e.n * pi), s = std::sin(e.n * pi);
    const Vec3 other = c * e.face0 + s * e.air, other_normal = -s * e.face0 + c * e.air;
    if (std::abs(v.dot(e.air)) <= Buildings::contact && v.dot(e.face0) > 0.0)
        return 0;
    if (std::abs(v.dot(other_normal)) <= Buildings::contact && v.dot(other) > 0.0)
        return 1;
    return -1;
}

} // namespace niceops
