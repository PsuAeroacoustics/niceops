#include "grating_paths.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <map>

#include "building_paths.hpp"
#include "facade_relief.hpp"
#include "point_paths.hpp"
#include "urban.hpp"
#include "util/constants.hpp"

namespace niceops {

namespace {

// A facade's frame: its start, unit base direction, outward (front) normal and
// length.
struct Frame
{
    Vec3 start, along, normal;
    double length = 0.0;
};

Frame frame_of(const Facade& f)
{
    Frame out;
    const double dx = f.x1 - f.x0, dy = f.y1 - f.y0;
    out.length = std::hypot(dx, dy);
    out.start = Vec3(f.x0, f.y0, 0.0);
    out.along = Vec3(dx / out.length, dy / out.length, 0.0);
    // facade_side is positive on the left of the base direction: the front.
    out.normal = Vec3(-out.along.y, out.along.x, 0.0);
    return out;
}

Vec3 below(const Vec3& p)
{
    return Vec3(p.x, p.y, -p.z);
}

// A direction mirrored in a facade's plane, z kept.
Vec3 mirrored(const Facade& f, const Vec3& v)
{
    const Vec3 h = mirror_direction(f, Vec3(v.x, v.y, 0.0));
    return Vec3(h.x, h.y, v.z);
}

} // namespace

bool grating_point(const Facade& facade, const Vec3& a, const Vec3& b, double kick, double& x, double& z)
{
    const Frame fr = frame_of(facade);
    // The specular point: where the line from a to b's mirror crosses the plane.
    const Vec3 b_image = mirror_in(facade, b);
    const double da = facade_side(facade, a), db = facade_side(facade, b_image);
    if (!(da > 0.0 && db < 0.0))
        return false;
    const Vec3 p0 = a + (da / (da - db)) * (b_image - a);
    x = (p0 - fr.start).dot(fr.along);
    z = p0.z;
    auto point = [&](double px, double pz) { return Vec3(fr.start.x + px * fr.along.x, fr.start.y + px * fr.along.y, pz); };
    auto value = [&](double px, double pz) {
        const Vec3 p = point(px, pz);
        return (p - a).norm() + (p - b).norm() + kick * px;
    };
    if (kick == 0.0)
        return true;
    double g = value(x, z);
    for (int iteration = 0; iteration < 100; iteration++)
    {
        const Vec3 p = point(x, z);
        const Vec3 ua = p - a, ub = p - b;
        const double ra = ua.norm(), rb = ub.norm();
        const double ta = fr.along.dot(ua) / ra, tb = fr.along.dot(ub) / rb;
        const double za = ua.z / ra, zb = ub.z / rb;
        const double gx = ta + tb + kick, gz = za + zb;
        // The Hessian of |P - c| is (I - u u^T) / r, in the plane's (along, z).
        const double hxx = (1.0 - ta * ta) / ra + (1.0 - tb * tb) / rb;
        const double hxz = -ta * za / ra - tb * zb / rb;
        const double hzz = (1.0 - za * za) / ra + (1.0 - zb * zb) / rb;
        const double det = hxx * hzz - hxz * hxz;
        if (std::hypot(gx, gz) < 1e-13)
            return true;
        if (!(det > 0.0))
            return false;
        double sx = -(hzz * gx - hxz * gz) / det, sz = -(-hxz * gx + hxx * gz) / det;
        // Backtracking: g is convex, so a step that does not lower it is
        // halved.
        double step = 1.0, next = value(x + sx, z + sz);
        while (!(next <= g) && step > 1e-12)
        {
            step *= 0.5;
            next = value(x + step * sx, z + step * sz);
        }
        if (!(next <= g))
            return std::hypot(gx, gz) < 1e-9;
        x += step * sx;
        z += step * sz;
        g = next;
        // With too large a kick the least lies at infinity.
        if (!(std::abs(x) < 1e7 && std::abs(z) < 1e7))
            return false;
        if (step * std::hypot(sx, sz) < 1e-13)
            return true;
    }
    return false;
}

double grating_spreading(const Facade& facade, const Vec3& a, const Vec3& p, const Vec3& b, double kick)
{
    const Frame fr = frame_of(facade);
    const Vec3 u1 = (p - a) / (p - a).norm();
    const Vec3 u2 = (b - p) / (b - p).norm();
    // Two directions across the launch.
    const Vec3 helper = std::abs(u1.z) < 0.9 ? Vec3(0.0, 0.0, 1.0) : Vec3(1.0, 0.0, 0.0);
    Vec3 e1 = u1.cross(helper);
    e1 = e1 / e1.norm();
    const Vec3 e2 = u1.cross(e1);
    const double h = 1e-6;
    // Where the perturbed ray crosses the plane through b across the central
    // one, or false where its order does not propagate.
    auto arrival = [&](const Vec3& launch, Vec3& at) {
        const Vec3 u = launch / launch.norm();
        const double un = fr.normal.dot(u);
        if (!(un < 0.0))
            return false;
        const Vec3 hit = a + (fr.normal.dot(fr.start - a) / un) * u;
        const double t = fr.along.dot(u) + kick;
        const double n2 = 1.0 - t * t - u.z * u.z;
        if (!(n2 > 0.0))
            return false;
        const Vec3 d = t * fr.along + Vec3(0.0, 0.0, u.z) + std::sqrt(n2) * fr.normal;
        at = hit + (u2.dot(b - hit) / u2.dot(d)) * d;
        return true;
    };
    Vec3 d[2];
    const Vec3 across[2] = {e1, e2};
    for (int i = 0; i < 2; i++)
    {
        Vec3 plus, minus;
        if (!arrival(u1 + h * across[i], plus) || !arrival(u1 - h * across[i], minus))
            return 0.0;
        d[i] = (plus - minus) / (2.0 * h);
    }
    const double area = d[0].cross(d[1]).norm();
    return area > 0.0 ? 1.0 / area : 0.0;
}

void grating_band_energies(const std::vector<GratingRay>& rays, const std::vector<double>& bands, double sound_speed,
        const double* power, double* out)
{
    const double half = std::pow(2.0, 1.0 / 6.0);
    std::vector<std::complex<double>> amplitude;
    std::vector<double> length;
    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < rays.size();)
    {
        std::size_t j = i;
        while (j < rays.size() && rays[j].group == rays[i].group)
            j++;
        const std::size_t band = rays[i].band;
        // The group's rays by length, as band_averaged_energy takes them.
        order.resize(j - i);
        for (std::size_t k = 0; k < order.size(); k++)
            order[k] = i + k;
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return rays[a].length < rays[b].length; });
        amplitude.clear();
        length.clear();
        for (std::size_t k : order)
        {
            amplitude.emplace_back(std::sqrt(rays[k].energy * (power != nullptr ? power[k] : 1.0)));
            length.push_back(rays[k].length);
        }
        const double k1 = 2.0 * pi * bands[band] / half / sound_speed, k2 = 2.0 * pi * bands[band] * half / sound_speed;
        out[band] += band_averaged_energy(amplitude, length, k1, k2);
        i = j;
    }
}

void grating_rays(const Buildings& buildings, const PathEnd& source, const PathEnd& microphone,
        const GratingInput& input, std::vector<GratingRay>& out)
{
    if (input.reliefs == nullptr || input.reflection == nullptr)
        return;
    const auto& reliefs = *input.reliefs;
    const auto& facades = buildings.facades();
    const std::size_t nb = input.bands.size();
    const std::vector<PointImage>& sources = source.tree().images;
    const std::vector<PointImage>& microphones = microphone.tree().images;
    const Vec3& s = source.point();
    const Vec3& m = microphone.point();
    auto is_relief = [&](std::size_t f) { return f < reliefs.size() && reliefs[f] != nullptr; };
    // The microphone's images that reach a relief last, by that facade.
    std::map<std::size_t, std::vector<std::size_t>> reaching;
    for (std::size_t j = 0; j < microphones.size(); j++)
    {
        const auto& fs = microphones[j].facades;
        if (!fs.empty() && is_relief(fs.back()) && fs.size() - 1 <= input.max_reflections)
            reaching[fs.back()].push_back(j);
    }
    const int source_floors = input.source_reflections ? 2 : 1;
    std::size_t groups = out.empty() ? 0 : out.back().group + 1;
    std::vector<std::pair<int, double>> shares;
    std::vector<Vec3> to_source, to_microphone;
    for (const PointImage& a_image : sources)
    {
        if (a_image.facades.empty() || !is_relief(a_image.facades.back()))
            continue;
        const std::size_t f = a_image.facades.back();
        const auto found = reaching.find(f);
        if (found == reaching.end())
            continue;
        const Facade& facade = facades[f];
        const FacadeRelief& relief = *reliefs[f];
        const Frame fr = frame_of(facade);
        const std::vector<std::size_t> a_chain(a_image.facades.begin(), a_image.facades.end() - 1);
        const Vec3 a_point = mirror_in(facade, a_image.point);
        for (std::size_t j : found->second)
        {
            const PointImage& b_image = microphones[j];
            const std::vector<std::size_t> b_chain(b_image.facades.begin(), b_image.facades.end() - 1);
            if (a_chain.size() + b_chain.size() > input.max_reflections)
                continue;
            const Vec3 b_point = mirror_in(facade, b_image.point);
            for (std::size_t band = 0; band < nb; band++)
            {
                const double f_band = input.bands[band];
                const double lambda = input.sound_speed / f_band;
                const int most = static_cast<int>(std::floor(2.0 * relief.period() / lambda));
                for (int order = -most; order <= most; order++)
                {
                    if (order == 0)
                        continue;
                    const std::size_t group = groups++;
                    for (int fa = 0; fa < source_floors; fa++)
                        for (int fb = 0; fb < 2; fb++)
                        {
                            const Vec3 pa = fa ? below(a_point) : a_point, pb = fb ? below(b_point) : b_point;
                            if (!(facade_side(facade, pa) > Buildings::contact
                                        && facade_side(facade, pb) > Buildings::contact))
                                continue;
                            const double kick = order * lambda / relief.period();
                            double x, z;
                            if (!grating_point(facade, pa, pb, kick, x, z) || x < 0.0 || x > fr.length
                                    || z < std::max(facade.bottom, 0.0) || z > facade.top)
                                continue;
                            const Vec3 p(fr.start.x + x * fr.along.x, fr.start.y + x * fr.along.y, z);
                            const double r1 = (p - pa).norm(), r2 = (pb - p).norm();
                            const Vec3 u1 = (p - pa) / r1, u2 = (pb - p) / r2;
                            const double horizontal = std::hypot(u1.x, u1.y);
                            if (!(horizontal > 1e-9))
                                continue;
                            const double sin_psi = -fr.along.dot(u1) / horizontal;
                            relief.orders(sin_psi, f_band * horizontal, input.sound_speed, shares);
                            double share = 0.0;
                            for (const auto& [mm, v] : shares)
                                if (mm == order)
                                    share = v;
                            if (share < input.least_share)
                                continue;
                            if (!folded_path_open(buildings, p, s, a_chain, fa == 1, &to_source)
                                    || !folded_path_open(buildings, p, m, b_chain, fb == 1, &to_microphone))
                                continue;
                            const std::vector<double>& rho = *input.reflection;
                            double energy = share * rho[f * nb + band] * grating_spreading(facade, pa, p, pb, kick);
                            // The other reflections: the source's side met
                            // backward from P, the microphone's forward.
                            Vec3 d = u1;
                            for (std::size_t k = a_chain.size(); k-- > 0;)
                            {
                                const std::size_t g = a_chain[k];
                                d = mirrored(facades[g], d);
                                energy *= rho[g * nb + band];
                                if (is_relief(g))
                                {
                                    const Frame gf = frame_of(facades[g]);
                                    const double hz = std::hypot(d.x, d.y);
                                    energy *= reliefs[g]->specular(-gf.along.dot(d) / hz, f_band * hz, input.sound_speed);
                                }
                            }
                            d = u2;
                            for (std::size_t k = b_chain.size(); k-- > 0;)
                            {
                                const std::size_t g = b_chain[k];
                                energy *= rho[g * nb + band];
                                if (is_relief(g))
                                {
                                    const Frame gf = frame_of(facades[g]);
                                    const double hz = std::hypot(d.x, d.y);
                                    energy *= reliefs[g]->specular(-gf.along.dot(d) / hz, f_band * hz, input.sound_speed);
                                }
                                d = mirrored(facades[g], d);
                            }
                            if (input.floor)
                            {
                                if (fa)
                                    energy *= input.floor(std::abs(p.z - pa.z) / r1, r1, band);
                                if (fb)
                                    energy *= input.floor(std::abs(pb.z - p.z) / r2, r2, band);
                            }
                            const double length = r1 + r2;
                            if (!input.air.empty())
                                energy *= std::exp(-2.0 * input.air[band] * length);
                            if (!(energy > 0.0))
                                continue;
                            GratingRay ray;
                            ray.band = band;
                            ray.energy = energy;
                            ray.length = length;
                            ray.facade = f;
                            ray.order = order;
                            ray.source_floor = fa == 1;
                            ray.microphone_floor = fb == 1;
                            ray.group = group;
                            // The launch is minus the length's gradient in
                            // the source: u1 taken back through the floor and
                            // the source's facades.
                            Vec3 launch = fa ? Vec3(u1.x, u1.y, -u1.z) : u1;
                            for (std::size_t k = a_chain.size(); k-- > 0;)
                                launch = mirrored(facades[a_chain[k]], launch);
                            ray.launch = launch;
                            ray.path = {s};
                            ray.path.insert(ray.path.end(), to_source.rbegin(), to_source.rend());
                            ray.path.push_back(p);
                            ray.path.insert(ray.path.end(), to_microphone.begin(), to_microphone.end());
                            ray.path.push_back(m);
                            out.push_back(std::move(ray));
                        }
                    }
                }
        }
    }
}

} // namespace niceops
