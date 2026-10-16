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

BuildingReceiver::BuildingReceiver(const Buildings& buildings, const Vec3& microphone,
        std::vector<std::size_t> run_facades, const BuildingPathOptions& options)
    : buildings_(&buildings)
    , microphone_(microphone)
    , run_facades_(std::move(run_facades))
    , options_(options)
{
    const auto& facades = buildings.facades();
    for (std::size_t i = 0; i < facades.size(); i++)
    {
        const Facade& f = facades[i];
        const double length = std::hypot(f.x1 - f.x0, f.y1 - f.y0);
        const double ax = (f.x1 - f.x0) / length, ay = (f.y1 - f.y0) / length;
        frames_.push_back({f.x0, f.y0, ax, ay, length, -ay, ax});
        // One facade per orientation: parallel normals mirror alike.
        const Frame& fi = frames_.back();
        const bool seen = std::any_of(planes_.begin(), planes_.end(), [&](std::size_t p) {
            return std::abs(frames_[p].nx * fi.ny - frames_[p].ny * fi.nx) < 1e-12;
        });
        if (!seen)
            planes_.push_back(i);
    }
    for (std::size_t e = 0; e < buildings.edges().size(); e++)
    {
        const Edge& edge = buildings.edges()[e];
        if (horizontal_distance_to_segment(microphone, edge.start, edge.start + edge.length * edge.direction)
                <= options.radius)
            edges_.push_back(e);
    }

    // The images a diffracted ray can reach the microphone through: depth
    // first over the facades within the radius, each a mirror in a facade the
    // current image stands in front of, and only between facades that can see
    // each other.  Which facades see each other is a prefilter, on a few
    // points of each face and in the floor's mirror world so that the floor
    // images are covered too; the exact test is each ray's own.
    images_.push_back({microphone, {}});
    image_edges_.push_back(edges_);
    if (!options.diffraction || options.max_diffracted_order == 0)
        return;
    std::vector<std::size_t> near;
    for (std::size_t i = 0; i < facades.size(); i++)
        if (facade_distance(facades[i], microphone) <= options.radius)
            near.push_back(i);
    auto side = [&](std::size_t f, const Vec3& p) {
        const Frame& fr = frames_[f];
        return (p.x - fr.x0) * fr.nx + (p.y - fr.y0) * fr.ny;
    };
    std::vector<std::vector<Vec3>> samples(near.size());
    for (std::size_t a = 0; a < near.size(); a++)
    {
        const Facade& f = facades[near[a]];
        const Frame& fr = frames_[near[a]];
        for (double along : {0.05, 0.3, 0.5, 0.7, 0.95})
            for (double h : {-0.5, -0.05, 0.05, 0.5})
                samples[a].emplace_back(fr.x0 + along * fr.length * fr.ax + 1e-3 * fr.nx,
                        fr.y0 + along * fr.length * fr.ay + 1e-3 * fr.ny, h * f.top);
    }
    const Vec3 floor_image(microphone.x, microphone.y, -microphone.z);
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
        Image image;
        Beam beam;
        std::size_t last; // index into near
    };
    auto mirror = [&](std::size_t f, const Vec3& p) {
        const Frame& fr = frames_[f];
        const double d = side(f, p);
        return Vec3(p.x - 2.0 * d * fr.nx, p.y - 2.0 * d * fr.ny, p.z);
    };
    auto child = [&](const Frontier& parent, std::size_t b, Frontier& out) {
        const Facade& f = facades[near[b]];
        double t0 = 0.0, t1 = 1.0;
        if (side(near[b], parent.image.point) <= Buildings::contact
                || (options.beams && !clip_to_beam(parent.beam, f.x0, f.y0, f.x1, f.y1, t0, t1)))
            return false;
        out.image.point = mirror(near[b], parent.image.point);
        out.image.facades = parent.image.facades;
        out.image.facades.push_back(near[b]);
        out.beam.open = !options.beams;
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
    const Frontier root {{microphone, {}}, Beam {}, near.size()};
    for (std::size_t a = 0; a < near.size(); a++)
    {
        Frontier next;
        if ((seen_from(microphone, a) || seen_from(floor_image, a)) && child(root, a, next))
            queue.push_back(std::move(next));
    }
    while (!queue.empty() && images_.size() < options.max_images)
    {
        Frontier f = std::move(queue.front());
        queue.pop_front();
        if (f.image.facades.size() < options.max_diffracted_order)
            for (std::size_t b = 0; b < near.size(); b++)
            {
                Frontier next;
                if (b != f.last && sees[f.last][b] && child(f, b, next))
                    queue.push_back(std::move(next));
            }
        // The edges a ray could diffract from to reach this image: those
        // crossing its beam.
        std::vector<std::size_t> edges;
        for (std::size_t e : edges_)
        {
            const Edge& edge = buildings.edges()[e];
            const Vec3 end = edge.start + edge.length * edge.direction;
            double t0, t1;
            if (clip_to_beam(f.beam, edge.start.x, edge.start.y, end.x, end.y, t0, t1))
                edges.push_back(e);
        }
        images_.push_back(std::move(f.image));
        image_edges_.push_back(std::move(edges));
    }
    truncated_ = !queue.empty();
}

// The first facade a ray from `pos` along the unit `dir` meets from outside,
// and where; facades.size() for none.  A facade `pos` lies on (within
// contact) is behind the ray: it is where the ray reflected or diffracted.
std::size_t BuildingReceiver::first_hit(const Vec3& pos, const Vec3& dir, bool below_floor, Vec3& hit) const
{
    const auto& facades = buildings_->facades();
    std::size_t best = facades.size();
    double best_t = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < facades.size(); i++)
    {
        const Frame& fr = frames_[i];
        const double side = (pos.x - fr.x0) * fr.nx + (pos.y - fr.y0) * fr.ny;
        const double approach = dir.x * fr.nx + dir.y * fr.ny;
        if (side <= Buildings::contact || approach >= 0.0)
            continue;
        const double t = side / -approach;
        if (t >= best_t)
            continue;
        const Vec3 p = pos + t * dir;
        const double along = (p.x - fr.x0) * fr.ax + (p.y - fr.y0) * fr.ay;
        const double bottom = below_floor ? -facades[i].top : facades[i].bottom;
        if (along >= -1e-9 && along <= fr.length + 1e-9 && p.z >= bottom - 1e-9 && p.z <= facades[i].top + 1e-9)
        {
            best = i;
            best_t = t;
            hit = p;
        }
    }
    return best;
}

// Whether the diffracted ray from `q` reaches the microphone (or, below the
// floor, its floor image) as `image` says: the unfolded line crosses each of
// its facades, the last mirrored first, within the facade and from in front,
// and no leg of the folded path passes through a building.
bool BuildingReceiver::diffracted_path_open(const Vec3& q, const Image& image, bool below_floor) const
{
    const auto& facades = buildings_->facades();
    auto side = [&](std::size_t f, const Vec3& p) {
        const Frame& fr = frames_[f];
        return (p.x - fr.x0) * fr.nx + (p.y - fr.y0) * fr.ny;
    };
    const Vec3 start = below_floor ? Vec3(microphone_.x, microphone_.y, -microphone_.z) : microphone_;
    std::vector<Vec3> targets {start};
    for (std::size_t f : image.facades)
    {
        const double d = side(f, targets.back());
        targets.emplace_back(targets.back().x - 2.0 * d * frames_[f].nx, targets.back().y - 2.0 * d * frames_[f].ny,
                targets.back().z);
    }
    Vec3 cur = q;
    for (std::size_t j = image.facades.size(); j >= 1; j--)
    {
        const std::size_t f = image.facades[j - 1];
        const double s0 = side(f, cur), s1 = side(f, targets[j]);
        if (!(s0 > Buildings::contact && s1 < -Buildings::contact))
            return false;
        const Vec3 hit = cur + (s0 / (s0 - s1)) * (targets[j] - cur);
        const Frame& fr = frames_[f];
        const double along = (hit.x - fr.x0) * fr.ax + (hit.y - fr.y0) * fr.ay;
        const double bottom = below_floor ? -facades[f].top : facades[f].bottom;
        if (along < -1e-9 || along > fr.length + 1e-9 || hit.z < bottom - 1e-9 || hit.z > facades[f].top + 1e-9)
            return false;
        if (buildings_->blocked(cur, hit, below_floor))
            return false;
        cur = hit;
    }
    return !buildings_->blocked(cur, targets[0], below_floor);
}

std::vector<BuildingRay> BuildingReceiver::rays(const Vec3& u_in) const
{
    const Vec3 u = incident_direction(u_in);
    const auto& facades = buildings_->facades();
    std::vector<BuildingRay> rays;

    // Specular rays: traced back from the microphone, or from its image in the
    // floor through the floor's mirror world, along each departure up to
    // max_turns changes of orientation could turn into u.  Facades are
    // vertical, so a reflection turns only the horizontal component.
    auto mirror_direction_in = [&](std::size_t f, const Vec3& v) {
        const Frame& fr = frames_[f];
        const double d = v.x * fr.nx + v.y * fr.ny;
        return Vec3(v.x - 2.0 * d * fr.nx, v.y - 2.0 * d * fr.ny, v.z);
    };
    std::vector<Vec3> departures {u};
    for (std::size_t level = 0, begin = 0; level < options_.max_turns && begin < departures.size(); level++)
    {
        const std::size_t end = departures.size();
        for (std::size_t i = begin; i < end && departures.size() < options_.max_departures; i++)
            for (std::size_t p : planes_)
            {
                const Vec3 w = mirror_direction_in(p, departures[i]);
                const bool seen = std::any_of(departures.begin(), departures.end(), [&](const Vec3& v) {
                    return std::abs(v.x - w.x) < 1e-12 && std::abs(v.y - w.y) < 1e-12;
                });
                if (!seen)
                    departures.push_back(w);
            }
        begin = end;
    }
    const Vec3 floor_image(microphone_.x, microphone_.y, -microphone_.z);
    for (int floor = 0; floor < 2; floor++)
        for (const Vec3& w : departures)
        {
            Vec3 pos = floor ? floor_image : microphone_;
            Vec3 dir = w;
            double travelled = 0.0;
            std::vector<std::size_t> hits;
            bool escaped = false;
            for (std::size_t bounce = 0; bounce <= options_.max_order; bounce++)
            {
                Vec3 hit;
                const std::size_t f = first_hit(pos, dir, floor == 1, hit);
                if (f == facades.size())
                {
                    escaped = true;
                    break;
                }
                if (bounce == options_.max_order)
                    break;
                travelled += (hit - pos).norm();
                pos = hit;
                dir = mirror_direction_in(f, dir);
                hits.push_back(f);
            }
            if (!escaped || std::abs(dir.x - u.x) > 1e-9 || std::abs(dir.y - u.y) > 1e-9)
                continue;
            BuildingRay ray;
            ray.length = -u.dot(pos) + travelled;
            ray.facades = std::move(hits);
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
                const Vec3 out = target - q;
                const double theta = edge.angle(out - out.dot(edge.direction) * edge.direction);
                if (!(theta > 0.0 && theta < edge.n * pi))
                    continue;
                Vec3 unused;
                if (first_hit(q, u, false, unused) != facades.size())
                    continue; // the wave does not reach this point of the edge
                if (!diffracted_path_open(q, image, floor == 1))
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
        const std::vector<double>& bands, double sound_speed, double* out) const
{
    const Vec3 u = incident_direction(u_in);
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
            double total = 0.0;
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
        double total = 0.0, open = 0.0;
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
        unsigned threads)
{
    const std::size_t nb = axes.bands.size(), ne = axes.elevations.size();
    std::vector<double> gains(axes.azimuths.size() * ne * nb);
    parallel_chunks(axes.azimuths.size(), threads, [&](std::size_t begin, std::size_t end) {
        std::vector<double> ratio(nb);
        for (std::size_t i = begin; i < end; i++)
            for (std::size_t j = 0; j < ne; j++)
            {
                const Vec3 u = local_direction(axes.azimuths[i], axes.elevations[j], frame_bearing);
                receiver.ratio(u, reflection, floor_reflection, axes.bands, sound_speed, ratio.data());
                for (std::size_t band = 0; band < nb; band++)
                    gains[(i * ne + j) * nb + band] = std::max(10.0 * std::log10(ratio[band]), urban_gain_floor_db);
            }
    });
    return gains;
}

} // namespace niceops
