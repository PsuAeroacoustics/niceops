#include "diffuse_facades.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

#include "building_paths.hpp"
#include "util/constants.hpp"
#include "util/parallel.hpp"

namespace niceops {

namespace {

double horizontal_distance_to_facade(const Facade& f, const Vec3& p)
{
    const double dx = f.x1 - f.x0, dy = f.y1 - f.y0;
    const double t = std::clamp(((p.x - f.x0) * dx + (p.y - f.y0) * dy) / (dx * dx + dy * dy), 0.0, 1.0);
    return std::hypot(p.x - (f.x0 + t * dx), p.y - (f.y0 + t * dy));
}

// What a patch's unit exitance sends to the point x, as irradiance on a face
// of unit normal `facing` there (cosine-weighted), or as the sum of the
// intensities arriving from it when `facing` is zero (an omnidirectional
// microphone): its area times cos at the patch (times cos at x) over pi r^2.
// Near the patch, under three of its sizes, it is summed over a 4 by 4 grid
// of its parts, where one point would misjudge it.
double transfer(const DiffuseFacades::Patch& q, const Vec3& x, const Vec3& facing, double z_sign)
{
    // The patch's image below the floor for z_sign < 0, facing down where the
    // patch faces up.
    const Vec3 normal(q.normal.x, q.normal.y, z_sign * q.normal.z);
    auto point = [&](const Vec3& c, double area) {
        const Vec3 from(c.x, c.y, z_sign * c.z);
        const Vec3 d = x - from;
        const double r2 = d.squaredNorm();
        if (r2 == 0.0)
            return 0.0;
        const double r = std::sqrt(r2);
        const double cos_q = normal.dot(d) / r;
        if (cos_q <= 0.0)
            return 0.0;
        double cos_x = 1.0;
        if (facing.squaredNorm() > 0.0)
        {
            cos_x = -facing.dot(d) / r;
            if (cos_x <= 0.0)
                return 0.0;
        }
        return area * cos_q * cos_x / (pi * r2);
    };
    const Vec3 c(q.center.x, q.center.y, q.center.z);
    const double distance = (x - Vec3(c.x, c.y, z_sign * c.z)).norm();
    if (distance >= 3.0 * std::max(q.width, q.height))
        return point(c, q.area);
    const int m = 4;
    double sum = 0.0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
        {
            const double a = ((i + 0.5) / m - 0.5) * q.width, h = ((j + 0.5) / m - 0.5) * q.height;
            sum += point(c + a * q.along + h * q.across, q.area / (m * m));
        }
    return sum;
}

} // namespace

DiffuseFacades::DiffuseFacades(const Buildings& buildings, const Vec3& microphone,
        const std::vector<double>& reflection, const std::vector<double>& scattering, double floor_reflection,
        const std::vector<double>& bands, double sound_speed, const DiffuseOptions& options, unsigned threads)
    : buildings_(&buildings)
    , options_(options)
    , bands_(bands.size())
    , floor_reflection_(floor_reflection)
{
    const auto& facades = buildings.facades();
    const std::size_t nb = bands.size();
    specular_.resize(reflection.size());
    reflection_ = reflection;
    for (std::size_t i = 0; i < reflection.size(); i++)
        specular_[i] = reflection[i] * (1.0 - scattering[i]);
    in_reach_.assign(facades.size(), 0);
    for (std::size_t f = 0; f < facades.size(); f++)
        in_reach_[f] = horizontal_distance_to_facade(facades[f], microphone) <= options.radius;
    // With the roofs a building scatters as a whole, walls and roof, where
    // any of its walls is in reach: its roof's coefficient is one of its
    // walls', so the two must scatter together.
    if (options.roofs)
    {
        std::vector<char> building(buildings.size(), 0);
        for (std::size_t f = 0; f < facades.size(); f++)
            building[buildings.facade_building(f)] |= in_reach_[f];
        for (std::size_t f = 0; f < facades.size(); f++)
            in_reach_[f] = building[buildings.facade_building(f)];
    }
    for (std::size_t f = 0; f < facades.size(); f++)
    {
        if (!in_reach_[f])
            continue;
        const Facade& facade = facades[f];
        const double length = std::hypot(facade.x1 - facade.x0, facade.y1 - facade.y0);
        const Vec3 along((facade.x1 - facade.x0) / length, (facade.y1 - facade.y0) / length, 0.0);
        const Vec3 normal(-along.y, along.x, 0.0);
        const double tall = facade.top - facade.bottom;
        const auto na = static_cast<std::size_t>(std::max(1.0, std::ceil(length / options.patch_size)));
        const auto nh = static_cast<std::size_t>(std::max(1.0, std::ceil(tall / options.patch_size)));
        for (std::size_t i = 0; i < na; i++)
            for (std::size_t j = 0; j < nh; j++)
            {
                Patch p;
                p.width = length / na;
                p.height = tall / nh;
                p.area = p.width * p.height;
                p.along = along;
                p.normal = normal;
                p.facade = f;
                p.center = Vec3(facade.x0, facade.y0, facade.bottom) + ((i + 0.5) * p.width) * along
                        + Vec3(0.0, 0.0, (j + 0.5) * p.height) + 1e-3 * normal;
                // A patch inside a neighbour is not a wall anyone sees, and
                // one above a gable's sloping top, or below a sloping bottom
                // over a neighbour's pitched roof, is no wall.
                if (!buildings.buried(p.center) && p.center.z < facade_top(facade, p.center)
                        && p.center.z > facade_bottom(facade, p.center))
                    patches_.push_back(p);
            }
    }
    if (options.roofs)
        add_roof_patches();
    const std::size_t n = patches_.size();
    for (std::size_t p = 0; p < n; p++)
        all_.push_back(p);
    planes_ = orientations_near(buildings, microphone, options.radius);
    rho_.resize(n * nb);
    rho_s_.resize(n * nb);
    for (std::size_t p = 0; p < n; p++)
        for (std::size_t b = 0; b < nb; b++)
        {
            const std::size_t f = patches_[p].facade;
            rho_[p * nb + b] = reflection[f * nb + b];
            rho_s_[p * nb + b] = reflection[f * nb + b] * scattering[f * nb + b];
            // A relief is its walls': the roof, which takes a wall's
            // coefficients, does not scatter by it.
            if (patches_[p].roof && options.scattering_at && options.scattering_at(f, 0.0, 1.0, b) >= 0.0)
                rho_s_[p * nb + b] = 0.0;
        }

    // The form factors, through the air and off the floor (in energy), and
    // the microphone's gathering, its two paths interfering band-averaged as
    // open ground's do.
    form_factors_.resize(n);
    gather_.assign(n * nb, 0.0);
    // A microphone on a wall, a level roof or a slope hears each patch twice,
    // the second off that surface at the microphone, of one length:
    // (1 + R)^2, R the surface's specular amplitude, as the rays' partners
    // there take it.
    std::vector<double> surface(nb, 1.0);
    const Buildings::Support on = buildings.support(microphone);
    if (on.kind == Buildings::Support::Kind::facade || on.kind == Buildings::Support::Kind::roof
            || on.kind == Buildings::Support::Kind::slope)
    {
        const std::size_t f = on.kind == Buildings::Support::Kind::facade ? on.index : buildings.roof_facade(on.index);
        for (std::size_t b = 0; b < nb; b++)
        {
            const double r = f == Buildings::no_facade
                    ? 1.0
                    : std::sqrt(reflection[f * nb + b] * (1.0 - scattering[f * nb + b]));
            surface[b] = (1.0 + r) * (1.0 + r);
        }
    }
    const Vec3 none(0.0, 0.0, 0.0);
    const Vec3 image(microphone.x, microphone.y, -microphone.z);
    const double half = std::pow(2.0, 1.0 / 6.0);
    parallel_chunks(n, threads, [&](std::size_t begin, std::size_t end) {
        for (std::size_t p = begin; p < end; p++)
        {
            const Patch& pp = patches_[p];
            for (std::size_t q = 0; q < n; q++)
            {
                const Patch& qq = patches_[q];
                // The same face exchanges nothing.
                if (q == p || (qq.facade == pp.facade && qq.roof == pp.roof))
                    continue;
                // Each path's transfer first: one that faces away is zero and
                // needs no visibility test.
                double f = 0.0;
                const double direct = transfer(qq, pp.center, pp.normal, 1.0);
                if (direct > 0.0 && !buildings.blocked(pp.center, qq.center))
                    f += direct;
                if (floor_reflection > 0.0)
                {
                    const double mirrored = transfer(qq, pp.center, pp.normal, -1.0);
                    const Vec3 below(qq.center.x, qq.center.y, -qq.center.z);
                    if (mirrored > 0.0 && !buildings.blocked(pp.center, below, true))
                        f += floor_reflection * mirrored;
                }
                if (f > 0.0)
                    form_factors_[p].emplace_back(q, f);
            }
            const double t1 = transfer(pp, microphone, none, 1.0);
            const double e1 = t1 > 0.0 && !buildings.blocked(microphone, pp.center) ? t1 : 0.0;
            const double t2 = floor_reflection > 0.0 ? transfer(pp, image, none, 1.0) : 0.0;
            const double e2 = t2 > 0.0 && !buildings.blocked(pp.center, image, true) ? floor_reflection * t2 : 0.0;
            const double d1 = (microphone - pp.center).norm(), d2 = (image - pp.center).norm(), delta = d2 - d1;
            for (std::size_t b = 0; b < nb; b++)
            {
                const double k1 = 2.0 * pi * bands[b] / half / sound_speed, k2 = 2.0 * pi * bands[b] * half / sound_speed;
                const double x = 0.5 * (k2 - k1) * delta;
                const double sinc = std::abs(x) < 1e-8 ? 1.0 : std::sin(x) / x;
                const double a1 = options.air.empty() ? e1 : e1 * std::exp(-2.0 * options.air[b] * d1);
                const double a2 = options.air.empty() ? e2 : e2 * std::exp(-2.0 * options.air[b] * d2);
                gather_[p * nb + b] = surface[b]
                        * (a1 + a2 + 2.0 * std::sqrt(a1 * a2) * std::cos(0.5 * (k1 + k2) * delta) * sinc);
            }
        }
    });

    // The reverse solve, per band: y = g + F^T D y, then w = S y.
    weight_.assign(n * nb, 0.0);
    std::vector<double> y(n), next(n);
    for (std::size_t b = 0; b < nb; b++)
    {
        for (std::size_t p = 0; p < n; p++)
            y[p] = gather_[p * nb + b];
        bool done = false;
        for (std::size_t iteration = 0; iteration < options.max_iterations && !done; iteration++)
        {
            for (std::size_t p = 0; p < n; p++)
                next[p] = gather_[p * nb + b];
            for (std::size_t p = 0; p < n; p++)
            {
                const double source = rho_[p * nb + b] * y[p];
                if (source == 0.0)
                    continue;
                for (const auto& [q, f] : form_factors_[p])
                    next[q] += f * source;
            }
            double change = 0.0, largest = 0.0;
            for (std::size_t p = 0; p < n; p++)
            {
                change = std::max(change, std::abs(next[p] - y[p]));
                largest = std::max(largest, std::abs(next[p]));
            }
            y.swap(next);
            done = change <= options.tolerance * largest;
        }
        converged_ = converged_ && done;
        for (std::size_t p = 0; p < n; p++)
        {
            const bool elsewhere = options.first_scattering_elsewhere && !patches_[p].roof
                    && options.first_scattering_elsewhere(patches_[p].facade);
            weight_[p * nb + b] = rho_s_[p * nb + b] * (elsewhere ? y[p] - gather_[p * nb + b] : y[p]);
        }
    }
    double largest = 0.0;
    for (double w : weight_)
        largest = std::max(largest, w);
    for (std::size_t p = 0; p < n; p++)
    {
        double w = 0.0;
        for (std::size_t b = 0; b < nb; b++)
            w = std::max(w, weight_[p * nb + b]);
        if (w > 1e-12 * largest)
            active_.push_back(p);
    }
}

void DiffuseFacades::add_roof_patches()
{
    const Buildings& buildings = *buildings_;
    for (std::size_t b = 0; b < buildings.size(); b++)
    {
        const std::size_t f = buildings.roof_facade(b);
        if (f == Buildings::no_facade || !in_reach_[f])
            continue;
        // A flat roof, a pitched roof's level deck and each of its slopes, cut
        // alike in plan: each sample lies on one of them (roof_surface), and a
        // slope's patch is tilted with it.
        const Footprint& footprint = buildings.footprints()[b];
        const std::size_t m = footprint.x.size();
        // The patches run along the outline's longest side, so that a
        // rectangle's fit it; a cell the outline crosses keeps the part on
        // the roof, by a 4 by 4 grid of its points.
        std::size_t longest = 0;
        double best = 0.0;
        for (std::size_t i = 0; i < m; i++)
        {
            const double length = std::hypot(footprint.x[(i + 1) % m] - footprint.x[i],
                    footprint.y[(i + 1) % m] - footprint.y[i]);
            if (length > best)
            {
                best = length;
                longest = i;
            }
        }
        const double ox = footprint.x[longest], oy = footprint.y[longest];
        const Vec3 along((footprint.x[(longest + 1) % m] - ox) / best, (footprint.y[(longest + 1) % m] - oy) / best, 0.0);
        const Vec3 across(-along.y, along.x, 0.0);
        double u0 = std::numeric_limits<double>::infinity(), u1 = -u0, v0 = u0, v1 = -u0;
        for (std::size_t i = 0; i < m; i++)
        {
            const double u = (footprint.x[i] - ox) * along.x + (footprint.y[i] - oy) * along.y;
            const double v = (footprint.x[i] - ox) * across.x + (footprint.y[i] - oy) * across.y;
            u0 = std::min(u0, u);
            u1 = std::max(u1, u);
            v0 = std::min(v0, v);
            v1 = std::max(v1, v);
        }
        const auto nu = static_cast<std::size_t>(std::max(1.0, std::ceil((u1 - u0) / options_.patch_size)));
        const auto nv = static_cast<std::size_t>(std::max(1.0, std::ceil((v1 - v0) / options_.patch_size)));
        const double du = (u1 - u0) / static_cast<double>(nu), dv = (v1 - v0) / static_cast<double>(nv);
        const int k = 4;
        // Each cell's samples, each on one surface (roof_surface): the level
        // roof, or one slope -- a sample on a ridge is the first plane's, so
        // no area is counted twice.  A cell becomes a patch on each surface
        // its samples lie on, at their middle, of their share of its area.
        struct Part
        {
            int on = 0;
            double cu = 0.0, cv = 0.0;
        };
        std::map<long, Part> parts;
        for (std::size_t i = 0; i < nu; i++)
            for (std::size_t j = 0; j < nv; j++)
            {
                parts.clear();
                for (int s = 0; s < k; s++)
                    for (int t = 0; t < k; t++)
                    {
                        const double u = u0 + (static_cast<double>(i) + (s + 0.5) / k) * du;
                        const double v = v0 + (static_cast<double>(j) + (t + 0.5) / k) * dv;
                        const long surface = buildings.roof_surface(b, ox + u * along.x + v * across.x,
                                oy + u * along.y + v * across.y);
                        if (surface == Buildings::roof_none)
                            continue;
                        Part& part = parts[surface];
                        part.on++;
                        part.cu += u;
                        part.cv += v;
                    }
                for (const auto& [surface, part] : parts)
                {
                    const double cu = part.cu / part.on, cv = part.cv / part.on;
                    Patch p;
                    p.width = du;
                    p.height = dv;
                    p.area = du * dv * part.on / (k * k);
                    p.along = along;
                    p.across = across;
                    p.normal = Vec3(0.0, 0.0, 1.0);
                    p.facade = f;
                    p.roof = true;
                    const double x = ox + cu * along.x + cv * across.x, y = oy + cu * along.y + cv * across.y;
                    p.center = Vec3(x, y, footprint.top + 1e-3);
                    if (surface != Buildings::roof_level)
                    {
                        // A slope's cell is its plan cell lifted onto the
                        // plane: its sides along the lifted directions, its
                        // area the plan's over the plane's cosine.
                        const RoofPlane& r = footprint.roof[static_cast<std::size_t>(surface)];
                        const auto lift = [&](const Vec3& d) { return Vec3(d.x, d.y, r.a * d.x + r.b * d.y); };
                        const Vec3 along3 = lift(along), across3 = lift(across);
                        p.width = du * along3.norm();
                        p.height = dv * across3.norm();
                        p.area *= std::sqrt(1.0 + r.a * r.a + r.b * r.b);
                        p.along = along3 / along3.norm();
                        p.across = across3 / across3.norm();
                        p.normal = r.normal();
                        p.center = Vec3(x, y, r.at(x, y)) + 1e-3 * p.normal;
                    }
                    if (!buildings.buried(p.center))
                        patches_.push_back(p);
                }
            }
    }
}

void DiffuseFacades::irradiance(const Vec3& u, const std::vector<std::size_t>& which, std::vector<double>& e) const
{
    e.assign(patches_.size() * bands_, 0.0);
    const auto departures = specular_departures(*buildings_, planes_, u, options_.max_turns, options_.max_departures);
    SpecularTrace trace;
    for (std::size_t p : which)
    {
        const Patch& patch = patches_[p];
        for (int floor = 0; floor < 2; floor++)
        {
            // A level roof faces up: what the floor reflects rises away from
            // it, as its mirrored normal would say for every departure, so the
            // traces are saved.  A slope can face what the floor sends up.
            if (floor == 1 && (floor_reflection_ == 0.0 || patch.normal.z == 1.0))
                continue;
            const Vec3 start = floor ? Vec3(patch.center.x, patch.center.y, -patch.center.z) : patch.center;
            const Vec3 normal(patch.normal.x, patch.normal.y, floor ? -patch.normal.z : patch.normal.z);
            for (const Vec3& w : departures)
            {
                // Arriving at the face, or its image, from in front.
                const double c = w.dot(normal);
                if (c <= 0.0 || !trace_specular(*buildings_, start, w, u, floor == 1, options_.max_order, trace))
                    continue;
                for (std::size_t b = 0; b < bands_; b++)
                {
                    double carried = floor ? c * floor_reflection_ : c;
                    for (std::size_t f : trace.facades)
                        carried *= specular_[f * bands_ + b];
                    e[p * bands_ + b] += carried;
                }
            }
        }
    }
}

void DiffuseFacades::gathered(const std::vector<double>& e, double* out) const
{
    std::fill(out, out + bands_, 0.0);
    for (std::size_t p : active_)
        for (std::size_t b = 0; b < bands_; b++)
            out[b] += weight_[p * bands_ + b] * e[p * bands_ + b];
}

void DiffuseFacades::energy(const Vec3& u, double* out) const
{
    std::vector<double> e;
    irradiance(u, active_, e);
    gathered(e, out);
}

const std::vector<DiffuseFacades::PatchImage>& DiffuseFacades::patch_images(std::size_t a, std::size_t max_order,
        std::vector<PatchImage>& scratch) const
{
    if (!patch_images_.empty() && patch_image_order_ == max_order)
        return patch_images_[a];
    const ImageTree tree = image_tree(*buildings_, patches_[active_[a]].center, options_.radius, max_order,
            options_.max_images, true, {}, 0);
    scratch.clear();
    for (const PointImage& image : tree.images)
        // The patch's own face reflects what arrives, it does not light the
        // patch a second time.
        if (image.facades.size() <= max_order
                && (image.facades.empty() || patches_[active_[a]].roof
                        || image.facades[0] != patches_[active_[a]].facade))
            scratch.push_back({image.point, image.facades});
    return scratch;
}

void DiffuseFacades::prepare_point_lighting(std::size_t max_order, unsigned threads)
{
    patch_images_.clear();
    std::vector<std::vector<PatchImage>> images(active_.size());
    parallel_chunks(active_.size(), threads, [&](std::size_t begin, std::size_t end) {
        for (std::size_t a = begin; a < end; a++)
            patch_images(a, max_order, images[a]);
    });
    patch_images_ = std::move(images);
    patch_image_order_ = max_order;
}

void DiffuseFacades::point_irradiance(const Vec3& source, std::size_t max_order, const std::vector<double>& air,
        std::vector<double>& e, unsigned threads, const SourcePower* power) const
{
    e.assign(patches_.size() * bands_, 0.0);
    const auto& at = options_.scattering_at;
    const auto& facades = buildings_->facades();
    // With the roofs, the level roofs below the source and the slopes it
    // stands over, between it and the patches, each the source's mirror in
    // it: a path may reflect in one, as a ray does, counted among its
    // reflections and carrying the roof's rho (1 - s).  A level roof commutes
    // with the walls, so it may come anywhere along the patch's chain
    // (folded_roof_path_open); a slope does not, and comes first, next to
    // the source, as from a drone over the roofs onto the walls.
    constexpr std::size_t level_roof = static_cast<std::size_t>(-1);
    struct RoofMirror
    {
        std::size_t building = 0;
        std::size_t plane = level_roof; // the slope's, or a level roof at the top
        Vec3 image;
        std::size_t facade = Buildings::no_facade;
    };
    std::vector<RoofMirror> roofs;
    if (options_.roofs && options_.roof_lighting && max_order > 0 && !active_.empty())
    {
        double x0 = source.x, x1 = source.x, y0 = source.y, y1 = source.y;
        for (std::size_t p : active_)
        {
            x0 = std::min(x0, patches_[p].center.x);
            x1 = std::max(x1, patches_[p].center.x);
            y0 = std::min(y0, patches_[p].center.y);
            y1 = std::max(y1, patches_[p].center.y);
        }
        std::vector<std::uint32_t> near;
        buildings_->near_box(x0, x1, y0, y1, near);
        std::sort(near.begin(), near.end());
        for (const std::uint32_t b : near)
        {
            const Footprint& f = buildings_->footprints()[b];
            const std::size_t facade = buildings_->roof_facade(b);
            const bool deck = f.roof.empty()
                    || std::any_of(f.roof.begin(), f.roof.end(), [](const RoofPlane& q) { return q.level(); });
            if (deck && source.z > f.top + Buildings::contact)
                roofs.push_back({b, level_roof, Vec3(source.x, source.y, 2.0 * f.top - source.z), facade});
            for (std::size_t i = 0; i < f.roof.size(); i++)
            {
                const double d = f.roof[i].above(source);
                if (!f.roof[i].level() && d > Buildings::contact)
                    roofs.push_back({b, i, source - (2.0 * d) * f.roof[i].normal(), facade});
            }
        }
    }
    parallel_chunks(active_.size(), threads, [&](std::size_t begin, std::size_t end) {
        std::vector<Vec3> hits;
        std::vector<PatchImage> scratch;
        std::vector<double> sent(power != nullptr ? bands_ : 0);
        // With incidence-dependent scattering, each band's factor from the
        // walls the path reflects in and its deposit's scale at the patch.
        std::vector<double> by_incidence(at ? bands_ : 0);
        for (std::size_t a = begin; a < end; a++)
        {
            const std::size_t p = active_[a];
            const Patch& patch = patches_[p];
            // One path's deposit: `length` from the source to the patch's
            // image, arriving at the cosine `c`, launched along `launch`,
            // travelling along `d` toward the patch's image through `chain`,
            // with the floor's factor and the roof's facade (`roof`, or none).
            auto deposit = [&](double length, double c, const Vec3& launch, Vec3 d,
                                   const std::vector<std::size_t>& chain, double floor, const RoofMirror* roof) {
                if (power != nullptr)
                    (*power)(launch, length, sent.data());
                if (at)
                {
                    std::fill(by_incidence.begin(), by_incidence.end(), 1.0);
                    for (std::size_t k = chain.size(); k-- > 0;)
                    {
                        const std::size_t f = chain[k];
                        const FacadeHit hit = facade_hit(facades[f], d);
                        for (std::size_t b = 0; b < bands_; b++)
                        {
                            const double s = at(f, hit.sin_psi, hit.cos_elevation, b);
                            if (s >= 0.0 && specular_[f * bands_ + b] > 0.0)
                                by_incidence[b] *= rho_of(f, b) * (1.0 - s) / specular_[f * bands_ + b];
                        }
                        d = mirror_direction(facades[f], d);
                    }
                    if (!patch.roof)
                    {
                        const FacadeHit hit = facade_hit(facades[patch.facade], d);
                        for (std::size_t b = 0; b < bands_; b++)
                        {
                            const double s = at(patch.facade, hit.sin_psi, hit.cos_elevation, b);
                            const double nominal = rho_s_[p * bands_ + b] / rho_[p * bands_ + b];
                            if (s >= 0.0 && nominal > 0.0)
                                by_incidence[b] *= s / nominal;
                        }
                    }
                }
                for (std::size_t b = 0; b < bands_; b++)
                {
                    double carried = c / (length * length) * floor;
                    if (power != nullptr)
                        carried *= sent[b];
                    if (at)
                        carried *= by_incidence[b];
                    for (std::size_t f : chain)
                        carried *= specular_[f * bands_ + b];
                    if (roof != nullptr && roof->facade != Buildings::no_facade)
                        carried *= specular_[roof->facade * bands_ + b];
                    if (!air.empty())
                        carried *= std::exp(-2.0 * air[b] * length);
                    e[p * bands_ + b] += carried;
                }
            };
            // The cosine at the patch of the last leg, from `last` toward the
            // patch's point `end_point` (both unfolded alike), on the face
            // whose normal is `normal`.
            auto arriving = [](const Vec3& last, const Vec3& end_point, const Vec3& normal) {
                const Vec3 back = last - end_point;
                return normal.dot(back) / back.norm();
            };
            for (const PatchImage& image : patch_images(a, max_order, scratch))
            {
                for (int floor = 0; floor < 2; floor++)
                {
                    // A level roof, as above: the paths are saved.
                    if (floor == 1 && (floor_reflection_ == 0.0 || patch.normal.z == 1.0))
                        continue;
                    if (!folded_path_open(*buildings_, source, patch.center, image.facades, floor == 1, &hits))
                        continue;
                    const Vec3 target = floor ? Vec3(image.point.x, image.point.y, -image.point.z) : image.point;
                    const double length = (target - source).norm();
                    // The last leg, toward the patch: arriving from in front.
                    const Vec3 end_point = floor ? Vec3(patch.center.x, patch.center.y, -patch.center.z) : patch.center;
                    const Vec3 normal(patch.normal.x, patch.normal.y, floor ? -patch.normal.z : patch.normal.z);
                    const double c = arriving(hits.empty() ? source : hits.back(), end_point, normal);
                    if (!(c > 0.0) || !(length > 0.0))
                        continue;
                    deposit(length, c, (target - source) / length, (target - source) / length, image.facades,
                            floor ? floor_reflection_ : 1.0, nullptr);
                }
                if (image.facades.size() + 1 > max_order)
                    continue;
                for (const RoofMirror& roof : roofs)
                {
                    // Not the roof the patch lies on.
                    if (patch.roof && buildings_->facade_building(patch.facade) == roof.building)
                        continue;
                    const Footprint& f = buildings_->footprints()[roof.building];
                    if (roof.plane == level_roof)
                    {
                        const double h = f.top;
                        if (!(patch.center.z > h + Buildings::contact)
                                || !folded_roof_path_open(*buildings_, source, patch.center, image.facades,
                                        roof.building, &hits))
                            continue;
                        const Vec3 target(image.point.x, image.point.y, 2.0 * h - image.point.z);
                        const double length = (target - source).norm();
                        const double c = arriving(hits.empty() ? source : hits.back(),
                                Vec3(patch.center.x, patch.center.y, 2.0 * h - patch.center.z),
                                Vec3(patch.normal.x, patch.normal.y, -patch.normal.z));
                        if (!(c > 0.0) || !(length > 0.0))
                            continue;
                        deposit(length, c, (target - source) / length, (target - source) / length, image.facades, 1.0,
                                &roof);
                        continue;
                    }
                    const RoofPlane& plane = f.roof[roof.plane];
                    for (int floor = 0; floor < 2; floor++)
                    {
                        if (floor == 1 && (floor_reflection_ == 0.0 || patch.normal.z == 1.0))
                            continue;
                        // The slope between the source's mirror and the
                        // patch's image (off the floor or not): where the line
                        // between them crosses its plane, on the slope, the
                        // leg from the source open and the chain from there.
                        const Vec3 target = floor ? Vec3(image.point.x, image.point.y, -image.point.z) : image.point;
                        const double d_target = plane.above(target), d_image = -plane.above(roof.image);
                        if (!(d_target > Buildings::contact) || !(d_image > 0.0))
                            continue;
                        const Vec3 bounce = roof.image + (d_image / (d_image + d_target)) * (target - roof.image);
                        if (!buildings_->on_plane(roof.building, roof.plane, bounce.x, bounce.y)
                                || buildings_->blocked(source, bounce)
                                || !folded_path_open(*buildings_, bounce, patch.center, image.facades, floor == 1, &hits))
                            continue;
                        const double length = (target - roof.image).norm();
                        const Vec3 end_point = floor ? Vec3(patch.center.x, patch.center.y, -patch.center.z) : patch.center;
                        const Vec3 normal(patch.normal.x, patch.normal.y, floor ? -patch.normal.z : patch.normal.z);
                        const double c = arriving(hits.empty() ? bounce : hits.back(), end_point, normal);
                        const Vec3 launch = bounce - source;
                        if (!(c > 0.0) || !(length > 0.0) || !(launch.norm() > 0.0))
                            continue;
                        deposit(length, c, launch / launch.norm(), (target - roof.image) / length, image.facades,
                                floor ? floor_reflection_ : 1.0, &roof);
                    }
                }
            }
        }
    });
}

void DiffuseFacades::point_energy(const Vec3& source, std::size_t max_order, const std::vector<double>& air,
        double* out, unsigned threads, const SourcePower* power) const
{
    std::vector<double> e;
    point_irradiance(source, max_order, air, e, threads, power);
    gathered(e, out);
}

void DiffuseFacades::energy_forward(const Vec3& u, double* out) const
{
    std::vector<double> e;
    irradiance(u, all_, e);
    const std::size_t n = patches_.size();
    std::vector<double> exitance(n), next(n);
    for (std::size_t b = 0; b < bands_; b++)
    {
        for (std::size_t p = 0; p < n; p++)
            exitance[p] = rho_s_[p * bands_ + b] * e[p * bands_ + b];
        for (int iteration = 0; iteration < 100000; iteration++)
        {
            double change = 0.0, largest = 0.0;
            for (std::size_t p = 0; p < n; p++)
            {
                double arriving = 0.0;
                for (const auto& [q, f] : form_factors_[p])
                    arriving += f * exitance[q];
                next[p] = rho_s_[p * bands_ + b] * e[p * bands_ + b] + rho_[p * bands_ + b] * arriving;
                change = std::max(change, std::abs(next[p] - exitance[p]));
                largest = std::max(largest, next[p]);
            }
            exitance.swap(next);
            if (change <= 1e-14 * largest)
                break;
        }
        out[b] = 0.0;
        for (std::size_t p = 0; p < n; p++)
            out[b] += gather_[p * bands_ + b] * exitance[p];
    }
}

} // namespace niceops
