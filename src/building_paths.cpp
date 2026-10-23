#include "building_paths.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>

#include "urban.hpp"
#include "util/constants.hpp"
#include "util/parallel.hpp"
#include "wedge.hpp"

namespace niceops {

namespace {

// As the canyon: grazing incidence is degenerate for the wedges (the wave and
// its roof reflection merge), so a hair above it stands in; and a direction
// straight up keeps no round-off of a horizontal component, which would put
// the specular traces' escapes a world away.
Vec3 incident_direction(const Vec3& u)
{
    const double min_rise = std::sin(1e-3);
    const double h = std::hypot(u.x, u.y);
    if (h < 1e-9)
        return {0.0, 0.0, 1.0};
    if (u.z >= min_rise)
        return u;
    const double c = std::cos(1e-3);
    return {u.x / h * c, u.y / h * c, min_rise};
}

double horizontal_distance_to_segment(const Vec3& p, const Vec3& a, const Vec3& b)
{
    const double dx = b.x - a.x, dy = b.y - a.y;
    const double len2 = dx * dx + dy * dy;
    const double t = len2 > 0.0 ? std::clamp(((p.x - a.x) * dx + (p.y - a.y) * dy) / len2, 0.0, 1.0) : 0.0;
    return std::hypot(p.x - (a.x + t * dx), p.y - (a.y + t * dy));
}

double facade_distance(const Facade& f, const Vec3& p)
{
    return horizontal_distance_to_segment(p, Vec3(f.x0, f.y0, 0.0), Vec3(f.x1, f.y1, 0.0));
}

// A facade's base point, unit base direction and length, and the unit normal
// to its outside (the left of the base direction).
struct Frame
{
    double x0, y0, ax, ay, length, nx, ny;
};

Frame frame_of(const Facade& f)
{
    const double length = std::hypot(f.x1 - f.x0, f.y1 - f.y0);
    const double ax = (f.x1 - f.x0) / length, ay = (f.y1 - f.y0) / length;
    return {f.x0, f.y0, ax, ay, length, -ay, ax};
}

double side_of(const Frame& fr, const Vec3& p)
{
    return (p.x - fr.x0) * fr.nx + (p.y - fr.y0) * fr.ny;
}

Vec3 mirror_point(const Frame& fr, const Vec3& p)
{
    const double d = side_of(fr, p);
    return {p.x - 2.0 * d * fr.nx, p.y - 2.0 * d * fr.ny, p.z};
}

// A beam in plan view: the points from which a straight line to `apex`
// crosses the window, the segment a -> b, on the window's far side from the
// apex.  For an image of the microphone it is where a ray may come from and
// still reach the microphone by the image's reflections.  An open beam is the
// whole plane (the microphone itself).
struct Beam
{
    bool open = true;
    double px = 0.0, py = 0.0, ax = 0.0, ay = 0.0, bx = 0.0, by = 0.0;
};

double cross(double ux, double uy, double vx, double vy)
{
    return ux * vy - uy * vx;
}

// The part of the segment c -> d (plan view) inside the beam, as parameters
// [t0, t1] along it; false where it misses the beam.  A point (c = d) is
// tested as one.  Boundaries count as inside, so that facades meeting at a
// corner can pass a ray between them.
bool clip_to_beam(const Beam& beam, double cx, double cy, double dx, double dy, double& t0, double& t1)
{
    t0 = 0.0;
    t1 = 1.0;
    if (beam.open)
        return true;
    const double tolerance = 1e-9;
    // Keep where f >= 0, f linear along the segment with values f0 and f1 at
    // its ends.
    auto keep = [&](double f0, double f1) {
        if (f0 < -tolerance && f1 < -tolerance)
            return false;
        if (f0 < -tolerance)
            t0 = std::max(t0, f0 / (f0 - f1));
        else if (f1 < -tolerance)
            t1 = std::min(t1, f0 / (f0 - f1));
        return t0 <= t1;
    };
    const double wx = beam.bx - beam.ax, wy = beam.by - beam.ay;
    const double apex_side = cross(wx, wy, beam.px - beam.ax, beam.py - beam.ay) < 0.0 ? 1.0 : -1.0;
    auto beyond = [&](double x, double y) { return apex_side * cross(wx, wy, x - beam.ax, y - beam.ay); };
    const double ax = beam.ax - beam.px, ay = beam.ay - beam.py, bx = beam.bx - beam.px, by = beam.by - beam.py;
    const double turn = cross(ax, ay, bx, by) > 0.0 ? 1.0 : -1.0;
    auto after_a = [&](double x, double y) { return turn * cross(ax, ay, x - beam.px, y - beam.py); };
    auto before_b = [&](double x, double y) { return -turn * cross(bx, by, x - beam.px, y - beam.py); };
    if (!keep(beyond(cx, cy), beyond(dx, dy)) || !keep(after_a(cx, cy), after_a(dx, dy))
            || !keep(before_b(cx, cy), before_b(dx, dy)))
        return false;
    const double length = std::hypot(dx - cx, dy - cy);
    return length == 0.0 || (t1 - t0) * length > 1e-6;
}

} // namespace

std::vector<std::size_t> facing_facades(const Buildings& buildings, const Vec3& microphone, double radius)
{
    std::vector<std::size_t> out;
    const auto& facades = buildings.facades();
    for (std::size_t i = 0; i < facades.size(); i++)
        if (facade_side(facades[i], microphone) > Buildings::contact && facade_distance(facades[i], microphone) <= radius)
            out.push_back(i);
    return out;
}

std::vector<std::size_t> edges_near(const Buildings& buildings, const Vec3& point, double radius)
{
    return edges_near(buildings, std::vector<Vec3> {point}, radius);
}

std::vector<std::size_t> edges_near(const Buildings& buildings, const std::vector<Vec3>& points, double radius)
{
    std::vector<std::size_t> near;
    for (std::size_t e = 0; e < buildings.edges().size(); e++)
    {
        const Edge& edge = buildings.edges()[e];
        const Vec3 end = edge.start + edge.length * edge.direction;
        for (const Vec3& p : points)
            if (horizontal_distance_to_segment(p, edge.start, end) <= radius)
            {
                near.push_back(e);
                break;
            }
    }
    return near;
}

std::vector<std::size_t> orientations_near(const Buildings& buildings, const Vec3& point, double radius)
{
    std::vector<std::size_t> planes;
    const auto& facades = buildings.facades();
    for (std::size_t f = 0; f < facades.size(); f++)
    {
        const std::size_t o = buildings.orientation(f);
        if (facade_distance(facades[f], point) <= radius && std::find(planes.begin(), planes.end(), o) == planes.end())
            planes.push_back(o);
    }
    return planes;
}

std::vector<Vec3> specular_departures(const Buildings& buildings, const std::vector<std::size_t>& planes,
        const Vec3& u, std::size_t max_turns, std::size_t max_departures)
{
    std::vector<Vec3> departures {u};
    for (std::size_t level = 0, begin = 0; level < max_turns && begin < departures.size(); level++)
    {
        const std::size_t end = departures.size();
        for (std::size_t i = begin; i < end && departures.size() < max_departures; i++)
            for (std::size_t p : planes)
            {
                // Two orientations not quite at right angles (or at any
                // angle that divides a turn) compose into a rotation by
                // round-off, not the identity: departures within 1e-6 are
                // one.  A trace escapes "along u" within 1e-7, and a path
                // turns its departure rigidly, so no ray is found twice.
                const Vec3 w = buildings.reflect(p, departures[i]);
                const bool seen = std::any_of(departures.begin(), departures.end(), [&](const Vec3& v) {
                    return std::abs(v.x - w.x) < 1e-6 && std::abs(v.y - w.y) < 1e-6;
                });
                if (!seen)
                    departures.push_back(w);
            }
        begin = end;
    }
    return departures;
}

bool trace_specular(const Buildings& buildings, const Vec3& start, const Vec3& departure, const Vec3& u,
        bool below_floor, std::size_t max_order, SpecularTrace& trace)
{
    trace.facades.clear();
    trace.travelled = 0.0;
    Vec3 pos = start, dir = departure;
    for (std::size_t bounce = 0;; bounce++)
    {
        Vec3 hit;
        const std::size_t f = buildings.first_hit(pos, dir, below_floor, hit);
        if (f == buildings.facades().size())
        {
            if (below_floor && pos.z < 0.0 && !buildings.escapes(pos, dir, true))
                return false;
            break;
        }
        if (bounce == max_order || (below_floor && std::min(pos.z, hit.z) < 0.0 && buildings.blocked(pos, hit, true)))
            return false;
        trace.travelled += (hit - pos).norm();
        pos = hit;
        dir = buildings.reflect(f, dir);
        trace.facades.push_back(f);
    }
    trace.exit = pos;
    return std::abs(dir.x - u.x) <= 1e-7 && std::abs(dir.y - u.y) <= 1e-7;
}

bool trace_off_roof(const Buildings& buildings, const Vec3& start, const Vec3& departure, const Vec3& u,
        std::size_t max_order, SpecularTrace& trace, std::size_t& roof)
{
    trace.facades.clear();
    trace.travelled = 0.0;
    if (!(departure.z < 0.0))
        return false;
    const auto& footprints = buildings.footprints();
    thread_local std::vector<std::uint32_t> near;
    Vec3 pos = start, dir = departure;
    for (;;)
    {
        Vec3 hit;
        const std::size_t f = buildings.first_hit(pos, dir, false, hit);
        // The leg ends at the facade, or at the floor.
        const bool walled = f < buildings.facades().size();
        const Vec3 end = walled ? hit : pos + (pos.z / -dir.z) * dir;
        // The first level roof the leg comes down on, if any.
        buildings.near_box(std::min(pos.x, end.x), std::max(pos.x, end.x), std::min(pos.y, end.y),
                std::max(pos.y, end.y), near);
        double best = std::numeric_limits<double>::infinity();
        std::size_t found = buildings.size();
        for (const std::uint32_t b : near)
        {
            const double h = footprints[b].top;
            if (!(pos.z > h + Buildings::contact) || !(end.z < h))
                continue;
            const double t = (pos.z - h) / -dir.z;
            const Vec3 p = pos + t * dir;
            if (t < best && buildings.on_roof(b, p.x, p.y))
            {
                best = t;
                found = b;
            }
        }
        if (found < buildings.size())
        {
            Vec3 p = pos + best * dir;
            p.z = footprints[found].top;
            if (buildings.blocked(pos, p))
                return false;
            trace.travelled += (p - pos).norm();
            roof = found;
            SpecularTrace rest;
            if (!trace_specular(buildings, p, Vec3(dir.x, dir.y, -dir.z), u, false, max_order - trace.facades.size(),
                        rest))
                return false;
            trace.facades.insert(trace.facades.end(), rest.facades.begin(), rest.facades.end());
            trace.travelled += rest.travelled;
            trace.exit = rest.exit;
            return true;
        }
        // Falling, a leg can come down on a roof that is not level, or into a
        // building through its top.
        if (!walled || trace.facades.size() == max_order || buildings.blocked(pos, hit))
            return false;
        trace.travelled += (hit - pos).norm();
        pos = hit;
        dir = buildings.reflect(f, dir);
        trace.facades.push_back(f);
    }
}

ImageTree image_tree(const Buildings& buildings, const Vec3& point, double radius, std::size_t max_order,
        std::size_t max_images, bool beams, const std::vector<std::size_t>& candidate_edges, std::size_t edge_order)
{
    // The images a ray can reach the point through: breadth first over the
    // facades within the radius, each a mirror in a facade the current image
    // stands in front of, and only between facades that can see each other.
    // Which facades see each other is a prefilter, on a few points of each
    // face and in the floor's mirror world so that the floor images are
    // covered too; the exact test is each ray's own.
    ImageTree tree;
    tree.images.push_back({point, {}});
    tree.edges.push_back(candidate_edges);
    if (max_order == 0)
        return tree;
    const auto& facades = buildings.facades();
    std::vector<Frame> frames;
    frames.reserve(facades.size());
    for (const Facade& f : facades)
        frames.push_back(frame_of(f));
    std::vector<std::size_t> near;
    for (std::size_t i = 0; i < facades.size(); i++)
        if (facade_distance(facades[i], point) <= radius)
            near.push_back(i);
    auto side = [&](std::size_t f, const Vec3& p) { return side_of(frames[f], p); };
    std::vector<std::vector<Vec3>> samples(near.size());
    for (std::size_t a = 0; a < near.size(); a++)
    {
        const Facade& f = facades[near[a]];
        const Frame& fr = frames[near[a]];
        // Up the wall from its bottom (a neighbour's roof, where one adjoins
        // it) and down its image below the floor.
        for (double along : {0.05, 0.3, 0.5, 0.7, 0.95})
        {
            const Vec3 base(fr.x0 + along * fr.length * fr.ax + 1e-3 * fr.nx,
                    fr.y0 + along * fr.length * fr.ay + 1e-3 * fr.ny, 0.0);
            const double bottom = facade_bottom(f, base), top = facade_top(f, base);
            for (double h : {-0.5, -0.05, 0.05, 0.5})
                samples[a].emplace_back(base.x, base.y,
                        (h < 0.0 ? -1.0 : 1.0) * (bottom + std::abs(h) * (top - bottom)));
        }
    }
    const Vec3 floor_image(point.x, point.y, -point.z);
    auto seen_from = [&](const Vec3& p, std::size_t a) {
        if (side(near[a], p) <= Buildings::contact)
            return false;
        for (const auto& q : samples[a])
            if (!buildings.blocked(p, q, true))
                return true;
        return false;
    };
    // Whether two of them see each other, found when first asked: most pairs
    // never are, the beam having ruled them out (0 not yet known, 1 no, 2 yes).
    std::vector<char> known(near.size() * near.size(), 0);
    auto sees = [&](std::size_t a, std::size_t b) {
        if (a > b)
            std::swap(a, b);
        char& k = known[a * near.size() + b];
        if (k == 0)
        {
            bool found = false;
            for (std::size_t i = 0; i < samples[a].size() && !found; i++)
            {
                if (side(near[b], samples[a][i]) <= Buildings::contact)
                    continue;
                for (std::size_t j = 0; j < samples[b].size() && !found; j++)
                    found = side(near[a], samples[b][j]) > Buildings::contact
                            && !buildings.blocked(samples[a][i], samples[b][j], true);
            }
            k = found ? 2 : 1;
        }
        return k == 2;
    };
    // Each image carries its beam: a further mirror counts only in a facade,
    // or the part of one, that lies in it, and that part becomes the next
    // image's window.  Exact in plan view; heights and occlusion are left to
    // the prefilter and each ray's own test.
    struct Frontier
    {
        PointImage image;
        Beam beam;
        std::size_t last; // index into near
    };
    auto child = [&](const Frontier& parent, std::size_t b, Frontier& out) {
        const Facade& f = facades[near[b]];
        double t0 = 0.0, t1 = 1.0;
        if (side(near[b], parent.image.point) <= Buildings::contact
                || (beams && !clip_to_beam(parent.beam, f.x0, f.y0, f.x1, f.y1, t0, t1)))
            return false;
        out.image.point = mirror_point(frames[near[b]], parent.image.point);
        out.image.facades = parent.image.facades;
        out.image.facades.push_back(near[b]);
        out.beam.open = !beams;
        out.beam.px = out.image.point.x;
        out.beam.py = out.image.point.y;
        out.beam.ax = f.x0 + t0 * (f.x1 - f.x0);
        out.beam.ay = f.y0 + t0 * (f.y1 - f.y0);
        out.beam.bx = f.x0 + t1 * (f.x1 - f.x0);
        out.beam.by = f.y0 + t1 * (f.y1 - f.y0);
        out.last = b;
        return true;
    };
    // Breadth first, so that a search stopped at max_images has dropped the
    // highest orders, not whole branches.
    std::deque<Frontier> queue;
    const Frontier root {{point, {}}, Beam {}, near.size()};
    for (std::size_t a = 0; a < near.size(); a++)
    {
        Frontier next;
        if ((seen_from(point, a) || seen_from(floor_image, a)) && child(root, a, next))
            queue.push_back(std::move(next));
    }
    while (!queue.empty() && tree.images.size() < max_images)
    {
        Frontier f = std::move(queue.front());
        queue.pop_front();
        if (f.image.facades.size() < max_order)
            for (std::size_t b = 0; b < near.size(); b++)
            {
                Frontier next;
                if (b != f.last && child(f, b, next) && sees(f.last, b))
                    queue.push_back(std::move(next));
            }
        // The edges a ray could diffract from to reach this image: those
        // crossing its beam.
        static const std::vector<std::size_t> none;
        std::vector<std::size_t> edges;
        for (std::size_t e : f.image.facades.size() <= edge_order ? candidate_edges : none)
        {
            const Edge& edge = buildings.edges()[e];
            const Vec3 end = edge.start + edge.length * edge.direction;
            double t0, t1;
            if (clip_to_beam(f.beam, edge.start.x, edge.start.y, end.x, end.y, t0, t1))
                edges.push_back(e);
        }
        tree.images.push_back(std::move(f.image));
        tree.edges.push_back(std::move(edges));
    }
    tree.truncated = !queue.empty();
    return tree;
}

namespace {

// folded_path_open, each leg's openness decided by `leg_open(leg, a, b)`.
template <typename LegOpen>
bool folded_path_with(const Buildings& buildings, const Vec3& from, const Vec3& point,
        const std::vector<std::size_t>& facades, bool below_floor, std::vector<Vec3>* hits, const LegOpen& leg_open)
{
    if (hits)
        hits->clear();
    const auto& all = buildings.facades();
    // The images along the chain: targets[j] is the point mirrored in the
    // first j facades.  Per thread and reused: this runs for every candidate
    // ray of the point-source search.
    const Vec3 start = below_floor ? Vec3(point.x, point.y, -point.z) : point;
    thread_local std::vector<Vec3> targets;
    thread_local std::vector<Frame> frames;
    targets.assign(1, start);
    frames.clear();
    for (std::size_t f : facades)
    {
        frames.push_back(frame_of(all[f]));
        targets.push_back(mirror_point(frames.back(), targets.back()));
    }
    Vec3 cur = from;
    for (std::size_t j = facades.size(); j >= 1; j--)
    {
        const std::size_t f = facades[j - 1];
        const Frame& fr = frames[j - 1];
        const double s0 = side_of(fr, cur), s1 = side_of(fr, targets[j]);
        if (!(s0 > Buildings::contact && s1 < -Buildings::contact))
            return false;
        const Vec3 hit = cur + (s0 / (s0 - s1)) * (targets[j] - cur);
        const double along = (hit.x - fr.x0) * fr.ax + (hit.y - fr.y0) * fr.ay;
        // The wall, and below the floor its image too.
        const double top = facade_top(all[f], hit), bottom = facade_bottom(all[f], hit);
        const double z = below_floor && hit.z < 0.0 ? -hit.z : hit.z;
        if (along < -1e-9 || along > fr.length + 1e-9 || z < bottom - 1e-9 || z > top + 1e-9)
            return false;
        if (!leg_open(facades.size() - j, cur, hit))
            return false;
        if (hits)
            hits->push_back(hit);
        cur = hit;
    }
    return leg_open(facades.size(), cur, targets[0]);
}

} // namespace

bool folded_path_open(const Buildings& buildings, const Vec3& from, const Vec3& point,
        const std::vector<std::size_t>& facades, bool below_floor, std::vector<Vec3>* hits)
{
    return folded_path_with(buildings, from, point, facades, below_floor, hits,
            [&](std::size_t, const Vec3& a, const Vec3& b) { return !buildings.blocked(a, b, below_floor); });
}

bool folded_path_open(const Buildings& buildings, const Vec3& from, const Vec3& point,
        const std::vector<std::size_t>& facades, bool below_floor, std::vector<Vec3>* hits, const LegTest& leg_open)
{
    return folded_path_with(buildings, from, point, facades, below_floor, hits, leg_open);
}

bool folded_roof_path_open(const Buildings& buildings, const Vec3& from, const Vec3& point,
        const std::vector<std::size_t>& facades, std::size_t building, std::vector<Vec3>* hits, Vec3* bounce)
{
    if (hits)
        hits->clear();
    const double h = buildings.footprints()[building].top;
    if (!(from.z > h + 1e-9 && point.z > h + 1e-9))
        return false;
    const auto& all = buildings.facades();
    // A point of the unfolded line in the real world: past the roof, its
    // mirror in the roof's plane.
    auto real = [h](const Vec3& p) { return p.z < h ? Vec3(p.x, p.y, 2.0 * h - p.z) : p; };
    thread_local std::vector<Vec3> targets;
    thread_local std::vector<Frame> frames;
    targets.assign(1, Vec3(point.x, point.y, 2.0 * h - point.z));
    frames.clear();
    for (std::size_t f : facades)
    {
        frames.push_back(frame_of(all[f]));
        targets.push_back(mirror_point(frames.back(), targets.back()));
    }
    bool crossed = false;
    // One leg of the unfolded line: where it crosses the roof's plane, the
    // crossing is on the roof and the two real pieces are open.
    auto leg = [&](const Vec3& a, const Vec3& b) {
        if ((a.z - h) * (b.z - h) < 0.0)
        {
            Vec3 p = a + ((a.z - h) / (a.z - b.z)) * (b - a);
            p.z = h;
            if (!buildings.on_roof(building, p.x, p.y))
                return false;
            crossed = true;
            if (bounce)
                *bounce = p;
            return !buildings.blocked(real(a), p) && !buildings.blocked(p, real(b));
        }
        return !buildings.blocked(real(a), real(b));
    };
    Vec3 cur = from;
    for (std::size_t j = facades.size(); j >= 1; j--)
    {
        const std::size_t f = facades[j - 1];
        const Frame& fr = frames[j - 1];
        const double s0 = side_of(fr, cur), s1 = side_of(fr, targets[j]);
        if (!(s0 > Buildings::contact && s1 < -Buildings::contact))
            return false;
        const Vec3 hit = cur + (s0 / (s0 - s1)) * (targets[j] - cur);
        const double along = (hit.x - fr.x0) * fr.ax + (hit.y - fr.y0) * fr.ay;
        const double z = real(hit).z;
        if (along < -1e-9 || along > fr.length + 1e-9 || z < facade_bottom(all[f], hit) - 1e-9
                || z > facade_top(all[f], hit) + 1e-9)
            return false;
        if (!leg(cur, hit))
            return false;
        if (hits)
            hits->push_back(hit);
        cur = hit;
    }
    return leg(cur, targets[0]) && crossed;
}

BuildingReceiver::BuildingReceiver(const Buildings& buildings, const Vec3& microphone,
        std::vector<std::size_t> run_facades, const BuildingPathOptions& options)
    : buildings_(&buildings)
    , microphone_(microphone)
    , run_facades_(std::move(run_facades))
    , options_(options)
{
    const Buildings::Support on = buildings.support(microphone);
    if (on.kind == Buildings::Support::Kind::facade)
    {
        surface_ = on.index;
        on_surface_ = true;
        surface_facade_ = on.index;
        const Frame fr = frame_of(buildings.facades()[on.index]);
        surface_normal_ = Vec3(fr.nx, fr.ny, 0.0);
    }
    else if (options.roofs && (on.kind == Buildings::Support::Kind::roof || on.kind == Buildings::Support::Kind::slope))
    {
        on_surface_ = true;
        surface_building_ = on.index;
        surface_facade_ = buildings.roof_facade(on.index);
        surface_normal_ = Vec3(0.0, 0.0, 1.0);
        const Footprint& f = buildings.footprints()[on.index];
        if (on.kind == Buildings::Support::Kind::slope)
            for (std::size_t i = 0; i < f.roof.size(); i++)
                if (!f.roof[i].level() && std::abs(f.roof[i].above(microphone)) <= Buildings::contact
                        && buildings.on_plane(on.index, i, microphone.x, microphone.y))
                    surface_normal_ = f.roof[i].normal();
    }
    if (options.roofs)
    {
        std::vector<std::uint32_t> near;
        buildings.near_box(microphone.x - options.radius, microphone.x + options.radius,
                microphone.y - options.radius, microphone.y + options.radius, near);
        std::sort(near.begin(), near.end());
        for (const std::uint32_t b : near)
        {
            const Footprint& f = buildings.footprints()[b];
            const bool deck = f.roof.empty()
                    || std::any_of(f.roof.begin(), f.roof.end(), [](const RoofPlane& q) { return q.level(); });
            if (deck && microphone.z > f.top + Buildings::contact)
                roofs_.push_back({b, level_roof, Vec3(0.0, 0.0, 1.0), microphone.z - f.top});
            for (std::size_t i = 0; i < f.roof.size(); i++)
            {
                const double d = f.roof[i].above(microphone);
                if (!f.roof[i].level() && d > Buildings::contact)
                    roofs_.push_back({b, i, f.roof[i].normal(), d});
            }
        }
    }
    edges_ = edges_near(buildings, microphone, options.radius);
    planes_ = orientations_near(buildings, microphone, options.radius);
    ImageTree tree = image_tree(buildings, microphone, options.radius,
            options.diffraction ? options.max_diffracted_order : 0, options.max_images, options.beams, edges_);
    images_ = std::move(tree.images);
    image_edges_ = std::move(tree.edges);
    truncated_ = tree.truncated;
    image_feet_.resize(images_.size());
    for (std::size_t i = 0; i < images_.size(); i++)
    {
        const Vec3& p = images_[i].point;
        for (std::size_t e : image_edges_[i])
        {
            const Edge& edge = buildings.edges()[e];
            for (int floor = 0; floor < 2; floor++)
                image_feet_[i].push_back(foot_of(edge, floor ? Vec3(p.x, p.y, -p.z) : p));
        }
    }
}

BuildingReceiver::Foot BuildingReceiver::foot_of(const Edge& edge, const Vec3& target)
{
    const Edge::Projection place = edge.project(target);
    Foot foot;
    foot.along = place.along;
    foot.across = place.normal.norm();
    if (foot.across >= 1e-9)
    {
        // t = along - cot(beta) across within [0, length].
        const double low = (foot.along - edge.length) / foot.across, high = foot.along / foot.across;
        foot.low = low - 1e-9 * std::max(1.0, std::abs(low));
        foot.high = high + 1e-9 * std::max(1.0, std::abs(high));
        const double theta = edge.angle(place.normal);
        // Edge::angle is in [0, 2 pi), so the face-0 side's slack
        // is at the top.
        foot.air = theta < edge.n * pi + 1e-9 || theta > 2.0 * pi - 1e-9;
    }
    return foot;
}

std::vector<BuildingRay> BuildingReceiver::rays(const Vec3& u_in) const
{
    const Vec3 u = incident_direction(u_in);
    const auto& facades = buildings_->facades();
    std::vector<BuildingRay> rays;

    // Specular rays: traced back from the microphone, or from its image in the
    // floor through the floor's mirror world, along each departure up to
    // max_turns changes of orientation could turn into u.
    const auto departures = specular_departures(*buildings_, planes_, u, options_.max_turns, options_.max_departures);
    const Vec3 floor_image(microphone_.x, microphone_.y, -microphone_.z);
    // A microphone on a wall hears nothing through it: a departure into the
    // wall's own building, which the trace would start past the wall it
    // stands on, is no ray.
    const auto into_wall = [&](const Vec3& w) { return on_surface_ && w.dot(surface_normal_) < 0.0; };
    SpecularTrace trace;
    for (int floor = 0; floor < 2; floor++)
        for (const Vec3& w : departures)
        {
            if (into_wall(w))
                continue;
            if (!trace_specular(*buildings_, floor ? floor_image : microphone_, w, u, floor == 1, options_.max_order,
                        trace))
                continue;
            BuildingRay ray;
            ray.length = -u.dot(trace.exit) + trace.travelled;
            ray.facades = trace.facades;
            ray.floor = floor;
            rays.push_back(std::move(ray));
        }
    // Off a roof below the microphone: falling from it to a level roof,
    // through facades or not, and up from there; or straight down to a
    // slope, the departure mirrored in its plane, and up from there.
    const auto off_roof = [&](std::size_t building) {
        BuildingRay ray;
        ray.length = -u.dot(trace.exit) + trace.travelled;
        ray.facades = trace.facades;
        ray.roof = true;
        ray.roof_facade = buildings_->roof_facade(building);
        rays.push_back(std::move(ray));
    };
    if (std::any_of(roofs_.begin(), roofs_.end(), [](const Roof& r) { return r.plane == level_roof; }))
        for (const Vec3& w : departures)
        {
            const Vec3 down(w.x, w.y, -w.z);
            std::size_t building = 0;
            if (into_wall(down)
                    || !trace_off_roof(*buildings_, microphone_, down, u, options_.max_order, trace, building))
                continue;
            off_roof(building);
        }
    for (const Roof& roof : roofs_)
    {
        if (roof.plane == level_roof)
            continue;
        for (const Vec3& w : departures)
        {
            const double wn = w.dot(roof.normal);
            if (!(wn > 0.0))
                continue;
            const Vec3 v = w - (2.0 * wn) * roof.normal;
            const double t = roof.depth / wn;
            const Vec3 p = microphone_ + t * v;
            if (into_wall(v) || !buildings_->on_plane(roof.building, roof.plane, p.x, p.y)
                    || buildings_->blocked(microphone_, p)
                    || !trace_specular(*buildings_, p, w, u, false, options_.max_order, trace))
                continue;
            trace.travelled += t;
            off_roof(roof.building);
        }
    }
    // The run's own first-order reflections, and their floor partners where
    // the reflection from above is open too.
    for (auto& ray : rays)
        if (ray.floor == 0 && !ray.roof && ray.facades.size() == 1
                && std::find(run_facades_.begin(), run_facades_.end(), ray.facades[0]) != run_facades_.end())
        {
            ray.run_evaluates = true;
            for (auto& partner : rays)
                if (partner.floor == 1 && partner.facades == ray.facades)
                    partner.run_evaluates = true;
        }
    if (!options_.diffraction)
        return rays;

    // Diffracted rays: from each edge the wave lights, on Keller's cone, to
    // the microphone, its floor image and their images in the facades.
    struct Incidence
    {
        bool lit = false;
        double c = 0.0, sin_beta = 0.0, theta0 = 0.0, cot = 0.0;
    };
    std::vector<Incidence> incidence(buildings_->edges().size());
    for (std::size_t e : edges_)
    {
        const Edge& edge = buildings_->edges()[e];
        Incidence& in = incidence[e];
        in.c = -u.dot(edge.direction);
        in.sin_beta = std::sqrt(std::max(0.0, 1.0 - in.c * in.c));
        in.theta0 = edge.angle(u - u.dot(edge.direction) * edge.direction);
        // Not when the wave runs along the edge, or comes from behind the
        // wedge.
        in.lit = in.sin_beta >= 1e-6 && in.theta0 > 0.0 && in.theta0 < edge.n * pi;
        in.cot = in.lit ? in.c / in.sin_beta : 0.0;
    }
    // How far a point stands over a roof's plane, along its normal.
    const auto height_over = [&](const Roof& roof, const Vec3& p) {
        return roof.plane == level_roof ? p.z - buildings_->footprints()[roof.building].top
                                        : buildings_->footprints()[roof.building].roof[roof.plane].above(p);
    };
    // The microphone's image in a roof's plane: in a level roof, any of its
    // images in the facades; in a slope, the microphone's own.
    const auto mirrored = [&](const Roof& roof, const Vec3& p) {
        return p - (2.0 * height_over(roof, p)) * roof.normal;
    };
    // The real path from an edge's point q off a roof to the microphone:
    // through the image's facades off a level roof (folded_roof_path_open),
    // or straight to the slope where the line to the image `target` crosses
    // it, and up.
    const auto roof_open = [&](const Roof& roof, const Vec3& q, const Vec3& target, const Image& image) {
        if (roof.plane == level_roof)
            return folded_roof_path_open(*buildings_, q, microphone_, image.facades, roof.building);
        const double dq = height_over(roof, q), dt = -height_over(roof, target);
        const Vec3 p = q + (dq / (dq + dt)) * (target - q);
        return buildings_->on_plane(roof.building, roof.plane, p.x, p.y) && !buildings_->blocked(q, p)
                && !buildings_->blocked(p, microphone_);
    };
    // One diffracted ray from edge e on Keller's cone to `target`: the image
    // `image` off the floor or not, or mirrored in `roof` (the floor then
    // left out), standing about the edge as `foot` says.
    const auto diffract = [&](std::size_t e, const Foot& foot, const Image& image, int floor, const Vec3& target,
                                  const Roof* roof) {
        const Incidence& in = incidence[e];
        const Edge& edge = buildings_->edges()[e];
        // The bounds keep the Keller point on the edge; with its ends
        // the reach decides past them (below).
        if (foot.across < 1e-9 || (!options_.edge_ends && (!foot.air || in.cot < foot.low || in.cot > foot.high)))
            return;
        // On Keller's cone: the diffracted ray leaves at the incident
        // ray's angle to the edge.  Past an end the end's wave goes on
        // (ended_wedge_coefficient), and the end stands in for the
        // Keller point in every check, within the reach.
        const double t = foot.along - in.c * foot.across / in.sin_beta;
        // Past an end that the line goes on beyond (Edge::beyond),
        // the point is another piece's or its image's.
        const bool outside = t < 0.0 || t > edge.length;
        const bool past = outside && past_true_end(edge, t);
        if ((outside && !(past && options_.edge_ends)) || (!outside && !in.lit))
            return;
        const Vec3 keller = edge.start + t * edge.direction;
        const Vec3 q = past ? edge.start + std::clamp(t, 0.0, edge.length) * edge.direction : keller;
        // The plane wave's path through x, less the Keller path.
        const double s_keller = (target - keller).norm();
        const auto excess = [&](const Vec3& x) { return -u.dot(x - keller) + (target - x).norm() - s_keller; };
        if (past && (keller.z < 0.0 || excess(q) > options_.end_reach))
            return;
        // Off a roof, the edge's point stands over its plane.
        if (roof != nullptr && !(height_over(*roof, q) > Buildings::contact))
            return;
        const Vec3 out = target - q;
        // A microphone on one of the wedge's faces (a wall it lies
        // on) is at that face's angle exactly, not a hair to either
        // side of it.
        const int face = on_surface_ && image.facades.empty() && roof == nullptr ? edge_face(edge, target) : -1;
        double theta = face == 0 ? 0.0
                : face == 1        ? edge.n * pi
                                   : edge.angle(out - out.dot(edge.direction) * edge.direction);
        // Past an end the end wave goes wherever the end sees, at the
        // nearer face's angle beyond its plane (point_paths.cpp's
        // onto_air).
        const double raw_theta = theta;
        if (past && face < 0 && theta > edge.n * pi)
            theta = theta - edge.n * pi < 2.0 * pi - theta ? edge.n * pi : 0.0;
        if (face < 0 && !past && !(theta > 0.0 && theta < edge.n * pi))
            return;
        // A point of an edge inside a neighbour (footprints that
        // overlap, or share a side only in part) diffracts nothing.
        if (buildings_->buried(q))
            return;
        Vec3 unused;
        if (buildings_->first_hit(q, u, false, unused) != facades.size())
            return; // the wave does not reach this point of the edge
        if (roof == nullptr ? !folded_path_open(*buildings_, q, microphone_, image.facades, floor == 1)
                            : !roof_open(*roof, q, target, image))
            return;
        BuildingRay ray;
        ray.diffracted = true;
        ray.edge = e;
        ray.angles = wedge_angles(edge.n, theta,
                in.theta0 <= edge.n * pi ? in.theta0 : in.theta0 - edge.n * pi < 2.0 * pi - in.theta0 ? edge.n * pi : 0.0);
        ray.sin_beta = in.sin_beta;
        ray.s = s_keller;
        ray.length = -u.dot(q) + out.norm();
        if (options_.edge_ends)
        {
            // The edge's whole line (continued); its ends' excesses.  The
            // amplitude is the Keller ray's (s), so past an end the
            // scale is 1 and only the phase is the end's.
            const Edge g = continued(edge);
            const double tg = t + (edge.start - g.start).dot(g.direction);
            if (std::isfinite(g.beyond[0]))
                ray.ends.lo = excess(g.start);
            if (std::isfinite(g.beyond[1]))
                ray.ends.hi = excess(g.start + g.length * g.direction);
            ray.ends.past = tg < 0.0 ? -1 : tg > g.length ? 1 : 0;
            // In the wedge's solid side, the end wave weighed as the
            // point-source model weighs it (weigh_vertex), the plane
            // wave's transition distance the receiver's.
            const auto depth = [&](double a) {
                return a <= edge.n * pi ? 0.0 : std::min(a - edge.n * pi, 2.0 * pi - a);
            };
            const double deepest = std::max(depth(in.theta0), depth(raw_theta));
            if (ray.ends.past != 0 && deepest > 0.0)
            {
                const double near = std::max(ray.ends.past < 0 ? ray.ends.lo : ray.ends.hi, 1e-300);
                const double sine = std::sin(0.5 * deepest);
                ray.ends.scale /= std::sqrt(1.0 + out.norm() * 2.0 * sine * sine / near);
            }
        }
        ray.facades = image.facades;
        ray.floor = floor;
        if (roof != nullptr)
        {
            ray.roof = true;
            ray.roof_facade = buildings_->roof_facade(roof->building);
        }
        rays.push_back(std::move(ray));
    };
    for (std::size_t i = 0; i < images_.size(); i++)
    {
        const Image& image = images_[i];
        for (std::size_t j = 0; j < image_edges_[i].size(); j++)
        {
            const std::size_t e = image_edges_[i][j];
            const Incidence& in = incidence[e];
            // An edge the wave meets from inside its wedge still has its ends'
            // waves, where they see it (below).
            if (!in.lit && !(options_.edge_ends && in.sin_beta >= 1e-6))
                continue;
            for (int floor = 0; floor < 2; floor++)
                diffract(e, image_feet_[i][2 * j + floor], image, floor,
                        floor ? Vec3(image.point.x, image.point.y, -image.point.z) : image.point, nullptr);
            for (const Roof& roof : roofs_)
                if (roof.plane == level_roof || i == 0)
                {
                    const Vec3 target = mirrored(roof, image.point);
                    diffract(e, foot_of(buildings_->edges()[e], target), image, 0, target, &roof);
                }
        }
    }
    return rays;
}

double reflected_amplitude(const BuildingRay& ray, const std::vector<double>& reflection, double floor_reflection)
{
    double amplitude = 1.0;
    for (std::size_t f : ray.facades)
        amplitude *= std::sqrt(reflection[f]);
    for (int i = 0; i < ray.floor; i++)
        amplitude *= std::sqrt(floor_reflection);
    if (ray.roof && ray.roof_facade != Buildings::no_facade)
        amplitude *= std::sqrt(reflection[ray.roof_facade]);
    return amplitude;
}

std::complex<double> ray_amplitude(const BuildingRay& ray, double reflected, double k)
{
    if (!ray.diffracted)
        return reflected;
    // The wedge's coefficient in the plane normal to the edge, at the
    // wavenumber and distance projected into it: for a plane wave this is the
    // oblique-incidence form, L = s sin^2(beta).
    const double r = ray.s * ray.sin_beta;
    return reflected * ended_wedge_coefficient(ray.angles, k * ray.sin_beta, r, ray.ends, k) / std::sqrt(r);
}

std::complex<double> ray_amplitude(const BuildingRay& ray, const std::vector<double>& reflection,
        double floor_reflection, double k)
{
    return ray_amplitude(ray, reflected_amplitude(ray, reflection, floor_reflection), k);
}

double BuildingReceiver::surface_factor(const BuildingRay& ray, double r) const
{
    if (!on_surface_)
        return 1.0;
    const Vec3 end = ray.floor ? Vec3(microphone_.x, microphone_.y, -microphone_.z) : microphone_;
    const bool on_face = ray.diffracted && ray.facades.empty() && edge_face(buildings_->edges()[ray.edge], end) >= 0;
    return on_face ? 0.5 * (1.0 + r) : 1.0 + r;
}

void BuildingReceiver::ratio(const Vec3& u_in, const std::vector<double>& reflection, double floor_reflection,
        const std::vector<double>& bands, double sound_speed, double* out, const DiffuseFacades* diffuse) const
{
    const Vec3 u = incident_direction(u_in);
    // The diffuse energy is a band's mean, in the same units as a ray's
    // |amplitude|^2: it adds to the band's energy.
    std::vector<double> scattered(bands.size(), 0.0);
    if (diffuse != nullptr)
        diffuse->energy(u, scattered.data());
    auto all = rays(u);
    std::vector<BuildingRay> kept;
    for (auto& ray : all)
        if (!ray.run_evaluates)
            kept.push_back(std::move(ray));
    std::sort(kept.begin(), kept.end(), [](const BuildingRay& a, const BuildingRay& b) { return a.length < b.length; });
    const std::size_t facade_count = buildings_->facades().size();
    std::vector<double> band_reflection(facade_count);
    // On a wall, its amplitude reflection in each band; the diffuse field
    // arrives with the wall's partner too, coherently at the wall, (1 + R)^2
    // in energy.
    auto surface_r = [&](std::size_t band) {
        return !on_surface_                             ? 0.0
                : surface_facade_ == Buildings::no_facade ? 1.0
                                                          : std::sqrt(reflection[surface_facade_ * bands.size() + band]);
    };
    for (std::size_t band = 0; band < bands.size() && on_surface_; band++)
        scattered[band] *= (1.0 + surface_r(band)) * (1.0 + surface_r(band));
    if (!options_.diffraction)
    {
        for (std::size_t band = 0; band < bands.size(); band++)
        {
            for (std::size_t f = 0; f < facade_count; f++)
                band_reflection[f] = reflection[f * bands.size() + band];
            double total = scattered[band];
            for (const auto& ray : kept)
                total += std::norm(surface_factor(ray, surface_r(band))
                        * ray_amplitude(ray, band_reflection, floor_reflection, 0.0));
            out[band] = std::max(total / (1.0 + floor_reflection), 1e-30);
        }
        return;
    }
    const double half = std::pow(2.0, 1.0 / 6.0);
    // Open ground: the direct ray and the floor's.
    const std::vector<std::complex<double>> ground {1.0, std::sqrt(floor_reflection)};
    const std::vector<double> ground_length {-u.dot(microphone_), -u.dot(microphone_) + 2.0 * microphone_.z * u.z};
    std::vector<std::complex<double>> amplitude(kept.size());
    std::vector<double> length(kept.size()), reflected(kept.size());
    for (std::size_t i = 0; i < kept.size(); i++)
        length[i] = kept[i].length;
    for (std::size_t band = 0; band < bands.size(); band++)
    {
        for (std::size_t f = 0; f < facade_count; f++)
            band_reflection[f] = reflection[f * bands.size() + band];
        // What the reflections leave of each ray is the same at every
        // wavenumber of the band.
        for (std::size_t i = 0; i < kept.size(); i++)
            reflected[i] = reflected_amplitude(kept[i], band_reflection, floor_reflection)
                    * surface_factor(kept[i], surface_r(band));
        const double f1 = bands[band] / half, f2 = bands[band] * half;
        double total = scattered[band] * (f2 - f1), open = 0.0;
        for (int part = 0; part < 3; part++)
        {
            const double g1 = f1 * std::pow(f2 / f1, part / 3.0), g2 = f1 * std::pow(f2 / f1, (part + 1) / 3.0);
            const double k1 = 2.0 * pi * g1 / sound_speed, k2 = 2.0 * pi * g2 / sound_speed;
            const double kc = 0.5 * (k1 + k2);
            for (std::size_t i = 0; i < kept.size(); i++)
                amplitude[i] = ray_amplitude(kept[i], reflected[i], kc);
            total += band_averaged_energy(amplitude, length, k1, k2) * (g2 - g1);
            open += band_averaged_energy(ground, ground_length, k1, k2) * (g2 - g1);
        }
        out[band] = std::max(total / open, 1e-30);
    }
}

std::vector<double> building_gains(const BuildingReceiver& receiver, const UrbanTable::Axes& axes,
        double frame_bearing, const std::vector<double>& reflection, double floor_reflection, double sound_speed,
        unsigned threads, const DiffuseFacades* diffuse, double azimuth_offset)
{
    const double offset = azimuth_offset * 180.0 / pi;
    const std::size_t nb = axes.bands.size(), ne = axes.elevations.size();
    std::vector<double> gains(axes.azimuths.size() * ne * nb);
    parallel_chunks(axes.azimuths.size(), threads, [&](std::size_t begin, std::size_t end) {
        std::vector<double> ratio(nb);
        for (std::size_t i = begin; i < end; i++)
            for (std::size_t j = 0; j < ne; j++)
            {
                const Vec3 u = local_direction(axes.azimuths[i] + offset, axes.elevations[j], frame_bearing);
                receiver.ratio(u, reflection, floor_reflection, axes.bands, sound_speed, ratio.data(), diffuse);
                for (std::size_t band = 0; band < nb; band++)
                    gains[(i * ne + j) * nb + band] = std::max(10.0 * std::log10(ratio[band]), urban_gain_floor_db);
            }
    });
    return gains;
}

} // namespace niceops
