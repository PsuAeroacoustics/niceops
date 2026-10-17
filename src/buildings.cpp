#include "buildings.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <limits>
#include <stdexcept>
#include <string>

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

double Edge::angle(const Vec3& v) const
{
    double theta = std::atan2(v.dot(air), v.dot(face0));
    if (theta < 0.0)
        theta += 2.0 * pi;
    return theta;
}

namespace {

// Every ring of a footprint, the outer one first.
template <typename Fn>
void for_each_ring(const Footprint& f, Fn fn)
{
    fn(f.x, f.y);
    for (const Ring& hole : f.holes)
        fn(hole.x, hole.y);
}

double signed_area(const std::vector<double>& x, const std::vector<double>& y)
{
    double area = 0.0;
    for (std::size_t i = 0, n = x.size(); i < n; i++)
        area += cross2(x[i], y[i], x[(i + 1) % n], y[(i + 1) % n]);
    return 0.5 * area;
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

} // namespace

Buildings::Buildings(std::vector<Footprint> footprints, double orientation_tolerance)
    : footprints_(std::move(footprints))
    , orientation_tolerance_(orientation_tolerance)
{
    x_min_ = y_min_ = std::numeric_limits<double>::infinity();
    x_max_ = y_max_ = -std::numeric_limits<double>::infinity();
    // Each side as built, for finding the ones two buildings share: its
    // ends, its facade and roof edge, and the corner edges at its ends.
    constexpr std::size_t none = static_cast<std::size_t>(-1);
    struct Side
    {
        std::size_t building;
        double x0, y0, x1, y1;
        std::size_t facade, roof, corner0, corner1;
    };
    std::vector<Side> sides;
    for (std::size_t b = 0; b < footprints_.size(); b++)
    {
        Footprint& f = footprints_[b];
        const std::string which = "building " + std::to_string(b + 1);
        if (!(f.top > 0.0) || !std::isfinite(f.top))
            throw std::invalid_argument(which + ": its roof must be above the floor");
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

        for_each_ring(f, [&](const std::vector<double>& x, const std::vector<double>& y) {
            const std::size_t m = x.size();
            // The vertical edge at each vertex where the solid turns convex (a
            // left turn from the side before it to the side after).
            std::vector<std::size_t> corner(m, none);
            for (std::size_t i = 0; i < m; i++)
            {
                const std::size_t j = (i + 1) % m, h = (i + m - 1) % m;
                const double turn = cross2(x[i] - x[h], y[i] - y[h], x[j] - x[i], y[j] - y[i]);
                if (turn <= 0.0)
                    continue;
                const double len = std::hypot(x[j] - x[i], y[j] - y[i]);
                const Vec3 along((x[j] - x[i]) / len, (y[j] - y[i]) / len, 0.0);
                Edge e;
                e.start = Vec3(x[i], y[i], 0.0);
                e.direction = Vec3(0.0, 0.0, 1.0);
                e.length = f.top;
                e.face0 = along;
                e.air = Vec3(along.y, -along.x, 0.0);
                const double back = std::hypot(x[h] - x[i], y[h] - y[i]);
                e.n = e.angle(Vec3((x[h] - x[i]) / back, (y[h] - y[i]) / back, 0.0)) / pi;
                e.building = b;
                e.vertical = true;
                corner[i] = edges_.size();
                edges_.push_back(e);
            }
            for (std::size_t i = 0; i < m; i++)
            {
                const std::size_t j = (i + 1) % m;
                // The side i -> j, the solid on its left; its facade runs the
                // other way, so that the outside is on the facade's left.
                Facade facade;
                facade.x0 = x[j];
                facade.y0 = y[j];
                facade.x1 = x[i];
                facade.y1 = y[i];
                facade.bottom = 0.0;
                facade.top = f.top;
                const double len = std::hypot(x[j] - x[i], y[j] - y[i]);
                const Vec3 along((x[j] - x[i]) / len, (y[j] - y[i]) / len, 0.0);
                const Vec3 outward(along.y, -along.x, 0.0);
                // The roof edge: the roof leaves it inward, the air is above.
                Edge roof;
                roof.start = Vec3(x[i], y[i], f.top);
                roof.direction = along;
                roof.length = len;
                roof.face0 = -1.0 * outward;
                roof.air = Vec3(0.0, 0.0, 1.0);
                roof.n = 1.5;
                roof.building = b;
                sides.push_back({b, x[i], y[i], x[j], y[j], facades_.size(), edges_.size(), corner[i], corner[j]});
                facades_.push_back(facade);
                facade_building_.push_back(b);
                edges_.push_back(roof);
            }
        });
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
        const double top = footprints_[side.building].top;
        const double neighbour = footprints_[sides[other->second].building].top;
        if (neighbour >= top)
        {
            drop_facade[side.facade] = 1;
            drop_edge[side.roof] = 1;
            for (std::size_t c : {side.corner0, side.corner1})
                if (c != none)
                    drop_edge[c] = 1;
            continue;
        }
        // Only what stands above the neighbour's roof is outside.
        facades_[side.facade].bottom = std::max(facades_[side.facade].bottom, neighbour);
        for (std::size_t c : {side.corner0, side.corner1})
            if (c != none && edges_[c].start.z < neighbour)
            {
                edges_[c].start.z = neighbour;
                edges_[c].length = top - neighbour;
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

    // The grids, their cells about one and a half typical items across.
    std::vector<double> facade_boxes;
    double building_size = 0.0, facade_size = 0.0;
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
    for (std::size_t b = 0; b < footprints_.size(); b++)
        building_size += std::max(bounds_[4 * b + 1] - bounds_[4 * b], bounds_[4 * b + 3] - bounds_[4 * b + 2]);
    if (!footprints_.empty())
    {
        building_grid_ = PlanGrid(bounds_, std::max(1.0, 1.5 * building_size / footprints_.size()));
        facade_grid_ = PlanGrid(facade_boxes, std::max(1.0, 1.5 * facade_size / facades_.size()));
    }
}

bool Buildings::strictly_inside(std::size_t building, double x, double y) const
{
    const Footprint& f = footprints_[building];
    bool in = false, near = false;
    for_each_ring(f, [&](const std::vector<double>& rx, const std::vector<double>& ry) {
        for (std::size_t i = 0, n = rx.size(), h = n - 1; i < n; h = i++)
        {
            if (distance_to_side(x, y, rx[h], ry[h], rx[i], ry[i]) < contact)
                near = true;
            if ((ry[i] > y) != (ry[h] > y) && x < rx[h] + (y - ry[h]) * (rx[i] - rx[h]) / (ry[i] - ry[h]))
                in = !in;
        }
    });
    return in && !near;
}

std::size_t Buildings::inside(double x, double y) const
{
    for (std::size_t b = 0; b < footprints_.size(); b++)
        if (x > bounds_[4 * b] && x < bounds_[4 * b + 1] && y > bounds_[4 * b + 2] && y < bounds_[4 * b + 3]
                && strictly_inside(b, x, y))
            return b;
    return footprints_.size();
}

bool Buildings::buried(const Vec3& p) const
{
    bool found = false;
    building_grid_.walk(p.x, p.y, 0.0, 0.0, 0.0, [&](const std::uint32_t* first, const std::uint32_t* last, double) {
        for (const std::uint32_t* i = first; i != last && !found; i++)
            found = p.z > contact && p.z < footprints_[*i].top - contact && strictly_inside(*i, p.x, p.y);
        return false;
    });
    return found;
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
        if (strictly_inside(building, a.x + t * d.x, a.y + t * d.y))
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
                    const double side = (pos.x - fr.x0) * fr.nx + (pos.y - fr.y0) * fr.ny;
                    const double approach = dir.x * fr.nx + dir.y * fr.ny;
                    if (side <= contact || approach >= 0.0)
                        continue;
                    const double t = side / -approach;
                    if (t >= best_t)
                        continue;
                    const Vec3 p = pos + t * dir;
                    const double along = (p.x - fr.x0) * fr.ax + (p.y - fr.y0) * fr.ay;
                    const double bottom = below_floor ? -facades_[i].top : facades_[i].bottom;
                    if (along >= -1e-9 && along <= fr.length + 1e-9 && p.z >= bottom - 1e-9
                            && p.z <= facades_[i].top + 1e-9)
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

} // namespace niceops
