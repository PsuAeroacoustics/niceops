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
    auto point = [&](const Vec3& c, double area) {
        // The patch's image below the floor for z_sign < 0.
        const Vec3 from(c.x, c.y, z_sign * c.z);
        const Vec3 d = x - from;
        const double r2 = d.squaredNorm();
        if (r2 == 0.0)
            return 0.0;
        const double r = std::sqrt(r2);
        const double cos_q = q.normal.dot(d) / r;
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
            sum += point(c + a * q.along + Vec3(0.0, 0.0, h), q.area / (m * m));
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
    for (std::size_t i = 0; i < reflection.size(); i++)
        specular_[i] = reflection[i] * (1.0 - scattering[i]);
    for (std::size_t f = 0; f < facades.size(); f++)
    {
        const Facade& facade = facades[f];
        if (horizontal_distance_to_facade(facade, microphone) > options.radius)
            continue;
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
                // A patch inside a neighbour is not a wall anyone sees.
                if (!buildings.buried(p.center))
                    patches_.push_back(p);
            }
    }
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
    const Vec3 none(0.0, 0.0, 0.0);
    const Vec3 image(microphone.x, microphone.y, -microphone.z);
    const double half = std::pow(2.0, 1.0 / 6.0);
    parallel_chunks(n, threads, [&](std::size_t begin, std::size_t end) {
        for (std::size_t p = begin; p < end; p++)
        {
            const Patch& pp = patches_[p];
            for (std::size_t q = 0; q < n; q++)
            {
                if (q == p || patches_[q].facade == pp.facade)
                    continue;
                const Patch& qq = patches_[q];
                double f = 0.0;
                if (!buildings.blocked(pp.center, qq.center))
                    f += transfer(qq, pp.center, pp.normal, 1.0);
                const Vec3 below(qq.center.x, qq.center.y, -qq.center.z);
                if (floor_reflection > 0.0 && !buildings.blocked(pp.center, below, true))
                    f += floor_reflection * transfer(qq, pp.center, pp.normal, -1.0);
                if (f > 0.0)
                    form_factors_[p].emplace_back(q, f);
            }
            const double e1 = buildings.blocked(microphone, pp.center) ? 0.0 : transfer(pp, microphone, none, 1.0);
            const double e2 = floor_reflection > 0.0 && !buildings.blocked(pp.center, image, true)
                    ? floor_reflection * transfer(pp, image, none, 1.0)
                    : 0.0;
            const double delta = (image - pp.center).norm() - (microphone - pp.center).norm();
            for (std::size_t b = 0; b < nb; b++)
            {
                const double k1 = 2.0 * pi * bands[b] / half / sound_speed, k2 = 2.0 * pi * bands[b] * half / sound_speed;
                const double x = 0.5 * (k2 - k1) * delta;
                const double sinc = std::abs(x) < 1e-8 ? 1.0 : std::sin(x) / x;
                gather_[p * nb + b] = e1 + e2 + 2.0 * std::sqrt(e1 * e2) * std::cos(0.5 * (k1 + k2) * delta) * sinc;
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
            weight_[p * nb + b] = rho_s_[p * nb + b] * y[p];
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
            if (floor == 1 && floor_reflection_ == 0.0)
                continue;
            const Vec3 start = floor ? Vec3(patch.center.x, patch.center.y, -patch.center.z) : patch.center;
            for (const Vec3& w : departures)
            {
                // Arriving at the face from in front.
                const double c = w.dot(patch.normal);
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

void DiffuseFacades::energy(const Vec3& u, double* out) const
{
    std::vector<double> e;
    irradiance(u, active_, e);
    std::fill(out, out + bands_, 0.0);
    for (std::size_t p : active_)
        for (std::size_t b = 0; b < bands_; b++)
            out[b] += weight_[p * bands_ + b] * e[p * bands_ + b];
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
