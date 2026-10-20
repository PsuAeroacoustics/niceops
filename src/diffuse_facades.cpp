#include "diffuse_facades.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

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
                // one above a gable's sloping top is no wall.
                if (!buildings.buried(p.center) && p.center.z < facade_top(facade, p.center))
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
        }

    // The form factors, through the air and off the floor (in energy), and
    // the microphone's gathering, its two paths interfering band-averaged as
    // open ground's do.
    form_factors_.resize(n);
    gather_.assign(n * nb, 0.0);
    // A microphone on a wall or a level roof hears each patch twice, the
    // second off that surface at the microphone, of one length: (1 + R)^2,
    // R the surface's specular amplitude, as the rays' partners there take it.
    std::vector<double> surface(nb, 1.0);
    const Buildings::Support on = buildings.support(microphone);
    if (on.kind == Buildings::Support::Kind::facade || on.kind == Buildings::Support::Kind::roof)
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
        const Footprint& footprint = buildings.footprints()[b];
        // A pitched roof does not scatter: its rays do not reflect in it.
        if (!footprint.roof.empty())
            continue;
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
        for (std::size_t i = 0; i < nu; i++)
            for (std::size_t j = 0; j < nv; j++)
            {
                int on = 0;
                double cu = 0.0, cv = 0.0;
                for (int s = 0; s < k; s++)
                    for (int t = 0; t < k; t++)
                    {
                        const double u = u0 + (static_cast<double>(i) + (s + 0.5) / k) * du;
                        const double v = v0 + (static_cast<double>(j) + (t + 0.5) / k) * dv;
                        if (buildings.on_roof(b, ox + u * along.x + v * across.x, oy + u * along.y + v * across.y))
                        {
                            on++;
                            cu += u;
                            cv += v;
                        }
                    }
                if (on == 0)
                    continue;
                cu /= on;
                cv /= on;
                Patch p;
                p.width = du;
                p.height = dv;
                p.area = du * dv * on / (k * k);
                p.along = along;
                p.across = across;
                p.normal = Vec3(0.0, 0.0, 1.0);
                p.facade = f;
                p.roof = true;
                p.center = Vec3(ox + cu * along.x + cv * across.x, oy + cu * along.y + cv * across.y, footprint.top + 1e-3);
                if (!buildings.buried(p.center))
                    patches_.push_back(p);
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
            // A roof faces up: what the floor reflects rises away from it, as
            // its mirrored normal would say for every departure, so the
            // traces are saved.
            if (floor == 1 && (floor_reflection_ == 0.0 || patch.roof))
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
    // How a wave travelling along d meets facade f: the sine of its
    // incidence and the cosine of its elevation.
    auto meet = [&](std::size_t f, const Vec3& d, double& sin_psi, double& cos_elevation) {
        const Facade& g = facades[f];
        const double ax = g.x1 - g.x0, ay = g.y1 - g.y0, length = std::hypot(ax, ay), h = std::hypot(d.x, d.y);
        sin_psi = h > 0.0 ? -(d.x * ax + d.y * ay) / (length * h) : 0.0;
        cos_elevation = h / d.norm();
    };
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
            for (const PatchImage& image : patch_images(a, max_order, scratch))
            {
                for (int floor = 0; floor < 2; floor++)
                {
                    // A roof, as above: the paths are saved.
                    if (floor == 1 && (floor_reflection_ == 0.0 || patch.roof))
                        continue;
                    if (!folded_path_open(*buildings_, source, patch.center, image.facades, floor == 1, &hits))
                        continue;
                    const Vec3 target = floor ? Vec3(image.point.x, image.point.y, -image.point.z) : image.point;
                    const double length = (target - source).norm();
                    // The last leg, toward the patch: arriving from in front.
                    const Vec3 last = hits.empty() ? source : hits.back();
                    const Vec3 end_point = floor ? Vec3(patch.center.x, patch.center.y, -patch.center.z) : patch.center;
                    const Vec3 back = last - end_point;
                    const Vec3 normal(patch.normal.x, patch.normal.y, floor ? -patch.normal.z : patch.normal.z);
                    const double c = normal.dot(back) / back.norm();
                    if (!(c > 0.0) || !(length > 0.0))
                        continue;
                    if (power != nullptr)
                        (*power)((target - source) / length, length, sent.data());
                    if (at)
                    {
                        std::fill(by_incidence.begin(), by_incidence.end(), 1.0);
                        Vec3 d = (target - source) / length;
                        double sin_psi, cos_elevation;
                        for (std::size_t k = image.facades.size(); k-- > 0;)
                        {
                            const std::size_t f = image.facades[k];
                            meet(f, d, sin_psi, cos_elevation);
                            for (std::size_t b = 0; b < bands_; b++)
                            {
                                const double s = at(f, sin_psi, cos_elevation, b);
                                if (s >= 0.0 && specular_[f * bands_ + b] > 0.0)
                                    by_incidence[b] *= rho_of(f, b) * (1.0 - s) / specular_[f * bands_ + b];
                            }
                            const Vec3 h = mirror_direction(facades[f], Vec3(d.x, d.y, 0.0));
                            d = Vec3(h.x, h.y, d.z);
                        }
                        if (!patch.roof)
                        {
                            meet(patch.facade, d, sin_psi, cos_elevation);
                            for (std::size_t b = 0; b < bands_; b++)
                            {
                                const double s = at(patch.facade, sin_psi, cos_elevation, b);
                                const double nominal = rho_s_[p * bands_ + b] / rho_[p * bands_ + b];
                                if (s >= 0.0 && nominal > 0.0)
                                    by_incidence[b] *= s / nominal;
                            }
                        }
                    }
                    for (std::size_t b = 0; b < bands_; b++)
                    {
                        double carried = c / (length * length) * (floor ? floor_reflection_ : 1.0);
                        if (power != nullptr)
                            carried *= sent[b];
                        if (at)
                            carried *= by_incidence[b];
                        for (std::size_t f : image.facades)
                            carried *= specular_[f * bands_ + b];
                        if (!air.empty())
                            carried *= std::exp(-2.0 * air[b] * length);
                        e[p * bands_ + b] += carried;
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
