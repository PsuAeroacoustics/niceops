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

ImageTree image_tree(const Buildings& buildings, const Vec3& point, double radius, std::size_t max_order,
        std::size_t max_images, bool beams, const std::vector<std::size_t>& candidate_edges)
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
        for (double along : {0.05, 0.3, 0.5, 0.7, 0.95})
            for (double h : {-0.5, -0.05, 0.05, 0.5})
                samples[a].emplace_back(fr.x0 + along * fr.length * fr.ax + 1e-3 * fr.nx,
                        fr.y0 + along * fr.length * fr.ay + 1e-3 * fr.ny, h * f.top);
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
    std::vector<std::vector<char>> sees(near.size(), std::vector<char>(near.size(), 0));
    for (std::size_t a = 0; a < near.size(); a++)
        for (std::size_t b = a + 1; b < near.size(); b++)
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
            sees[a][b] = sees[b][a] = found;
        }
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
                if (b != f.last && sees[f.last][b] && child(f, b, next))
                    queue.push_back(std::move(next));
            }
        // The edges a ray could diffract from to reach this image: those
        // crossing its beam.
        std::vector<std::size_t> edges;
        for (std::size_t e : candidate_edges)
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

bool folded_path_open(const Buildings& buildings, const Vec3& from, const Vec3& point,
        const std::vector<std::size_t>& facades, bool below_floor)
{
    const auto& all = buildings.facades();
    // The images along the chain: targets[j] is the point mirrored in the
    // first j facades.
    const Vec3 start = below_floor ? Vec3(point.x, point.y, -point.z) : point;
    std::vector<Vec3> targets {start};
    std::vector<Frame> frames;
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
        const double bottom = below_floor ? -all[f].top : all[f].bottom;
        if (along < -1e-9 || along > fr.length + 1e-9 || hit.z < bottom - 1e-9 || hit.z > all[f].top + 1e-9)
            return false;
        if (buildings.blocked(cur, hit, below_floor))
            return false;
        cur = hit;
    }
    return !buildings.blocked(cur, targets[0], below_floor);
}

BuildingReceiver::BuildingReceiver(const Buildings& buildings, const Vec3& microphone,
        std::vector<std::size_t> run_facades, const BuildingPathOptions& options)
    : buildings_(&buildings)
    , microphone_(microphone)
    , run_facades_(std::move(run_facades))
    , options_(options)
{
    for (std::size_t e = 0; e < buildings.edges().size(); e++)
    {
        const Edge& edge = buildings.edges()[e];
        if (horizontal_distance_to_segment(microphone, edge.start, edge.start + edge.length * edge.direction)
                <= options.radius)
            edges_.push_back(e);
    }
    planes_ = orientations_near(buildings, microphone, options.radius);
    ImageTree tree = image_tree(buildings, microphone, options.radius,
            options.diffraction ? options.max_diffracted_order : 0, options.max_images, options.beams, edges_);
    images_ = std::move(tree.images);
    image_edges_ = std::move(tree.edges);
    truncated_ = tree.truncated;
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
    SpecularTrace trace;
    for (int floor = 0; floor < 2; floor++)
        for (const Vec3& w : departures)
        {
            if (!trace_specular(*buildings_, floor ? floor_image : microphone_, w, u, floor == 1, options_.max_order,
                        trace))
                continue;
            BuildingRay ray;
            ray.length = -u.dot(trace.exit) + trace.travelled;
            ray.facades = trace.facades;
            ray.floor = floor;
            rays.push_back(std::move(ray));
        }
    // The run's own first-order reflections, and their floor partners where
    // the reflection from above is open too.
    for (auto& ray : rays)
        if (ray.floor == 0 && ray.facades.size() == 1
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
        double c = 0.0, sin_beta = 0.0, theta0 = 0.0;
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
    }
    for (std::size_t i = 0; i < images_.size(); i++)
    {
        const Image& image = images_[i];
        for (std::size_t e : image_edges_[i])
        {
            const Incidence& in = incidence[e];
            if (!in.lit)
                continue;
            const Edge& edge = buildings_->edges()[e];
            for (int floor = 0; floor < 2; floor++)
            {
                const Vec3 target = floor ? Vec3(image.point.x, image.point.y, -image.point.z) : image.point;
                const Vec3 d = target - edge.start;
                const double along = d.dot(edge.direction);
                const double across = (d - along * edge.direction).norm();
                if (across < 1e-9)
                    continue;
                // On Keller's cone: the diffracted ray leaves at the incident
                // ray's angle to the edge.
                const double t = along - in.c * across / in.sin_beta;
                if (t < 0.0 || t > edge.length)
                    continue;
                const Vec3 q = edge.start + t * edge.direction;
                // A point of an edge inside a neighbour (footprints that
                // overlap, or share a side only in part) diffracts nothing.
                if (buildings_->buried(q))
                    continue;
                const Vec3 out = target - q;
                const double theta = edge.angle(out - out.dot(edge.direction) * edge.direction);
                if (!(theta > 0.0 && theta < edge.n * pi))
                    continue;
                Vec3 unused;
                if (buildings_->first_hit(q, u, false, unused) != facades.size())
                    continue; // the wave does not reach this point of the edge
                if (!folded_path_open(*buildings_, q, microphone_, image.facades, floor == 1))
                    continue;
                BuildingRay ray;
                ray.diffracted = true;
                ray.edge = e;
                ray.n = edge.n;
                ray.theta = theta;
                ray.theta0 = in.theta0;
                ray.sin_beta = in.sin_beta;
                ray.s = out.norm();
                ray.length = -u.dot(q) + ray.s;
                ray.facades = image.facades;
                ray.floor = floor;
                rays.push_back(std::move(ray));
            }
        }
    }
    return rays;
}

std::complex<double> ray_amplitude(const BuildingRay& ray, const std::vector<double>& reflection,
        double floor_reflection, double k)
{
    double amplitude = 1.0;
    for (std::size_t f : ray.facades)
        amplitude *= std::sqrt(reflection[f]);
    for (int i = 0; i < ray.floor; i++)
        amplitude *= std::sqrt(floor_reflection);
    if (!ray.diffracted)
        return amplitude;
    // The wedge's coefficient in the plane normal to the edge, at the
    // wavenumber and distance projected into it: for a plane wave this is the
    // oblique-incidence form, L = s sin^2(beta).
    const double r = ray.s * ray.sin_beta;
    return amplitude * hard_wedge_coefficient(ray.n, ray.theta, ray.theta0, k * ray.sin_beta, r) / std::sqrt(r);
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
    if (!options_.diffraction)
    {
        for (std::size_t band = 0; band < bands.size(); band++)
        {
            for (std::size_t f = 0; f < facade_count; f++)
                band_reflection[f] = reflection[f * bands.size() + band];
            double total = scattered[band];
            for (const auto& ray : kept)
                total += std::norm(ray_amplitude(ray, band_reflection, floor_reflection, 0.0));
            out[band] = std::max(total / (1.0 + floor_reflection), 1e-30);
        }
        return;
    }
    const double half = std::pow(2.0, 1.0 / 6.0);
    // Open ground: the direct ray and the floor's.
    const std::vector<std::complex<double>> ground {1.0, std::sqrt(floor_reflection)};
    const std::vector<double> ground_length {-u.dot(microphone_), -u.dot(microphone_) + 2.0 * microphone_.z * u.z};
    std::vector<std::complex<double>> amplitude(kept.size());
    std::vector<double> length(kept.size());
    for (std::size_t i = 0; i < kept.size(); i++)
        length[i] = kept[i].length;
    for (std::size_t band = 0; band < bands.size(); band++)
    {
        for (std::size_t f = 0; f < facade_count; f++)
            band_reflection[f] = reflection[f * bands.size() + band];
        const double f1 = bands[band] / half, f2 = bands[band] * half;
        double total = scattered[band] * (f2 - f1), open = 0.0;
        for (int part = 0; part < 3; part++)
        {
            const double g1 = f1 * std::pow(f2 / f1, part / 3.0), g2 = f1 * std::pow(f2 / f1, (part + 1) / 3.0);
            const double k1 = 2.0 * pi * g1 / sound_speed, k2 = 2.0 * pi * g2 / sound_speed;
            const double kc = 0.5 * (k1 + k2);
            for (std::size_t i = 0; i < kept.size(); i++)
                amplitude[i] = ray_amplitude(kept[i], band_reflection, floor_reflection, kc);
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
