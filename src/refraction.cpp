#include "refraction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "util/constants.hpp"

namespace niceops {

namespace {

// Gradients below this bend a ray by less than a nanoradian over 10 km: treated
// as straight, which also keeps the closed forms away from centres at infinity.
constexpr double straight_gradient = 1e-9; // 1/s

// The circles every ray lies on are centred on the line c = 0, at z = -c0/a.
// With sigma the gradient's sign and L = c0/|a|, a point at height z sits
// h(z) = L + sigma z from that line, and the sound speed there is |a| h(z).
struct Circles
{
    double sigma;
    double a;
    double L;

    explicit Circles(const LinearProfile& p)
        : sigma(p.gradient > 0.0 ? 1.0 : -1.0), a(std::fabs(p.gradient)), L(p.sound_speed / a)
    {
    }
    double h(double z) const
    {
        return L + sigma * z;
    }
    // h(z)^2 - h(0)^2, over sigma: z (2L + sigma z).  Differences of squares of
    // h are formed from this, never from the squares, which are of order L^2
    // and would cancel to nothing as the gradient shrinks.
    double kappa(double z) const
    {
        return z * (2.0 * L + sigma * z);
    }
    // Upward refraction: horizontal distance from where the grazing ray touches
    // the ground to where it passes height z.
    double grazing_half_chord(double z) const
    {
        const double clamped = std::clamp(z, 0.0, L);
        return std::sqrt(clamped * (2.0 * L - clamped));
    }
};

} // namespace

double shadow_boundary(const LinearProfile& profile, double z_s, double z_r)
{
    if (!(profile.gradient < -straight_gradient))
        return std::numeric_limits<double>::infinity();
    const Circles c(profile);
    return c.grazing_half_chord(z_s) + c.grazing_half_chord(z_r);
}

DirectEigenray direct_eigenray(const LinearProfile& profile, double z_s, double z_r, double r)
{
    const Circles c(profile);
    const double h_s = c.h(z_s);
    const double h_r = c.h(z_r);
    DirectEigenray ray;

    if (c.sigma < 0.0)
    {
        const double t_s = c.grazing_half_chord(z_s);
        const double t_r = c.grazing_half_chord(z_r);
        if (r > t_s + t_r)
        {
            // Past the limiting ray: it leaves the source at the critical angle,
            // touches the ground a distance t_s out, and would rise to z_r at
            // t_s + t_r.  What reaches further creeps along the ground.
            ray.lit = false;
            ray.shadow_depth = r - (t_s + t_r);
            ray.launch = -std::asin(t_s / c.L);
            ray.arrival = std::asin(t_r / c.L);
            ray.spreading_range = r * c.L / h_s;
            ray.path_length = c.L * (std::asin(t_s / c.L) + std::asin(t_r / c.L)) + ray.shadow_depth;
            ray.travel_time = (std::atanh(t_s / c.L) + std::atanh(t_r / c.L)) / c.a
                    + ray.shadow_depth / profile.sound_speed;
            return ray;
        }
    }

    // The circle through (0, z_s) and (r, z_r) centred on the line c = 0 at
    // horizontal position x_c.  Equal distances from the centre give
    // 2 r x_c = r^2 + h_r^2 - h_s^2.
    const double x_c = (r * r + c.sigma * (c.kappa(z_r) - c.kappa(z_s))) / (2.0 * r);
    const double rho = std::hypot(x_c, h_s);
    // Tangents, oriented from source to receiver.  On the upper half-circle
    // (downward refraction) travel is clockwise; on the lower, anticlockwise.
    ray.launch = std::atan2(c.sigma * x_c, h_s);
    ray.arrival = std::atan2(-c.sigma * (r - x_c), h_r);
    // The ray tube: r |dr/d(launch)| sin|arrival| / cos(launch) works out to
    // (r rho / h_s)^2 on any circle of the family.
    ray.spreading_range = r * rho / h_s;
    const double chord = std::hypot(r, z_s - z_r);
    ray.path_length = 2.0 * rho * std::asin(std::min(1.0, chord / (2.0 * rho)));
    // dt = ds / c = dphi / (|a| cos phi), phi the angle from the vertical through
    // the centre, so t = atanh(sin phi) / |a| between the two ends.  Written as
    // one atanh of the addition formula, whose argument simplifies to
    // r rho / (h_s^2 + x_c r) -- a difference of two atanh would cancel.
    ray.travel_time = std::atanh(r * rho / (h_s * h_s + x_c * r)) / c.a;
    return ray;
}

ReflectedEigenray reflected_eigenray(const LinearProfile& profile, double z_s, double z_r,
        double r)
{
    const Circles c(profile);
    ReflectedEigenray ray;
    if (c.sigma < 0.0 && r > c.grazing_half_chord(z_s) + c.grazing_half_chord(z_r))
        return ray;

    const double h_s = c.h(z_s);
    const double kappa_s = c.kappa(z_s);
    const double kappa_r = c.kappa(z_r);

    // Launched at elevation theta the ray's circle has radius rho = h_s/cos
    // theta, and at height z sits d(z) = sqrt(rho^2 - h(z)^2) from its centre
    // horizontally.  Down to the ground and up to the receiver on the mirrored
    // circle, the horizontal legs are differences of those d, formed without
    // cancellation as kappa / (sum of d).  One exception: a ray launched upward
    // under downward refraction passes over its top first, and the first leg
    // is then the sum d(0) + d(z_s).
    struct Legs
    {
        double to_ground, from_ground, d0, rho;
    };
    auto legs = [&](double theta) {
        const double cosine = std::cos(theta);
        const double sine2 = std::sin(theta) * std::sin(theta);
        auto d = [&](double z, double kappa_z) {
            const double v = c.sigma * (kappa_s - kappa_z) + c.h(z) * c.h(z) * sine2;
            return std::sqrt(std::max(v, 0.0)) / cosine;
        };
        Legs l;
        l.rho = h_s / cosine;
        l.d0 = d(0.0, 0.0);
        const double d_s = h_s * std::fabs(std::tan(theta));
        const double d_r = d(z_r, kappa_r);
        l.to_ground = (c.sigma > 0.0 && theta >= 0.0) ? l.d0 + d_s
                                                     : (l.d0 + d_s > 0.0 ? kappa_s / (l.d0 + d_s) : 0.0);
        l.from_ground = l.d0 + d_r > 0.0 ? kappa_r / (l.d0 + d_r) : 0.0;
        return l;
    };
    auto miss = [&](double theta) {
        const auto l = legs(theta);
        return l.to_ground + l.from_ground - r;
    };

    // The range is monotonic in the launch angle: near zero straight down, and
    // rising to the shadow boundary (upward refraction) or without limit
    // (downward).  The straight image ray's angle is a good first split:
    // downward refraction lands it short, upward long.
    constexpr double near_vertical = 1e-9;
    double lo = -0.5 * pi + near_vertical;
    double hi;
    if (c.sigma < 0.0)
        hi = -std::asin(c.grazing_half_chord(z_s) / c.L);
    else
        hi = 0.5 * pi - near_vertical;
    const double image = -std::atan2(z_s + z_r, r);
    if (image > lo && image < hi)
    {
        if (miss(image) < 0.0)
            lo = image;
        else
            hi = image;
    }
    if (c.sigma > 0.0 && lo == image)
    {
        // Walk up from the image angle rather than bisecting all the way to
        // vertical, where the range grows without limit.
        double step = 0.05;
        double probe = std::min(lo + step, 0.5 * pi - near_vertical);
        while (miss(probe) < 0.0 && probe < 0.5 * pi - near_vertical)
        {
            lo = probe;
            step *= 2.0;
            probe = std::min(probe + step, 0.5 * pi - near_vertical);
        }
        hi = probe;
    }

    // Illinois-modified regula falsi: bracketed, and superlinear on a smooth
    // monotonic function.
    double f_lo = miss(lo);
    double f_hi = miss(hi);
    if (f_lo > 0.0 || f_hi < 0.0)
        return ray; // cannot bracket: leave it to the caller's fallback
    const double tolerance = 1e-10 * std::max(r, 1.0);
    double theta = hi;
    int side = 0;
    for (int iteration = 0; iteration < 100; iteration++)
    {
        theta = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
        const double f = miss(theta);
        if (std::fabs(f) <= tolerance || hi - lo < 1e-15)
            break;
        if (f > 0.0)
        {
            hi = theta;
            f_hi = f;
            if (side == 1)
                f_lo *= 0.5;
            side = 1;
        }
        else
        {
            lo = theta;
            f_lo = f;
            if (side == -1)
                f_hi *= 0.5;
            side = -1;
        }
    }

    const auto l = legs(theta);
    ray.exists = true;
    ray.launch = theta;
    // At the ground the ray is d(0) from its centre horizontally and h(0) = L
    // vertically, so its tangent makes this angle with the ground.
    ray.grazing = std::atan2(l.d0, c.L);
    auto arc = [&](double horizontal, double vertical) {
        return 2.0 * l.rho * std::asin(std::min(1.0, std::hypot(horizontal, vertical) / (2.0 * l.rho)));
    };
    ray.path_length = arc(l.to_ground, z_s) + arc(l.from_ground, z_r);
    return ray;
}

EffectiveSoundSpeed::EffectiveSoundSpeed(const AtmosphereProfile& profile, double frame_bearing,
        double z_low, double z_high)
    : z_low_(z_low)
{
    if (!profile.has_pressure())
        throw std::logic_error("EffectiveSoundSpeed needs a profile with pressure");
    if (!(z_high >= z_low))
        throw std::logic_error("EffectiveSoundSpeed needs z_high >= z_low");

    constexpr std::size_t max_nodes = 20001;
    const double span = z_high - z_low;
    nodes_ = std::max<std::size_t>(2, static_cast<std::size_t>(std::ceil(span)) + 1);
    step_ = 1.0;
    if (nodes_ > max_nodes)
    {
        nodes_ = max_nodes;
        step_ = span / static_cast<double>(max_nodes - 1);
    }
    inv_step_ = 1.0 / step_;

    // East/north into the local axes: +x along the compass bearing, +y ninety
    // degrees to its left, as LocalFrame::track_aligned has it.
    const double bearing = frame_bearing * pi / 180.0;
    const double sin_b = std::sin(bearing);
    const double cos_b = std::cos(bearing);

    values_.resize(nodes_);
    moments_.resize(nodes_);
    for (std::size_t node = 0; node < nodes_; node++)
    {
        const double z = z_low_ + step_ * static_cast<double>(node);
        const auto enu = profile.wind_at(z);
        values_[node] = {profile.at(z).sound_speed(), enu[0] * sin_b + enu[1] * cos_b,
                -enu[0] * cos_b + enu[1] * sin_b};
    }
    // The integrals of f and z f over each cell, exact for f linear across it.
    moments_[0].fill(0.0);
    for (std::size_t node = 1; node < nodes_; node++)
    {
        const double z0 = z_low_ + step_ * static_cast<double>(node - 1);
        moments_[node] = moments_[node - 1];
        for (std::size_t q = 0; q < 3; q++)
        {
            const double f = values_[node - 1][q];
            const double df = values_[node][q] - f;
            moments_[node][2 * q] += step_ * (f + 0.5 * df);
            moments_[node][2 * q + 1] += step_ * (z0 * f + 0.5 * (z0 * df + step_ * f) + step_ * df / 3.0);
        }
    }
}

void EffectiveSoundSpeed::locate(double z, std::size_t& node, double& fraction) const
{
    const double u = (z - z_low_) * inv_step_;
    node = u <= 0.0 ? 0 : std::min(static_cast<std::size_t>(u), nodes_ - 2);
    fraction = u - static_cast<double>(node);
}

EffectiveSoundSpeed::Moments EffectiveSoundSpeed::moments(double z) const
{
    std::size_t node;
    double t;
    locate(z, node, t);
    const double z0 = z_low_ + step_ * static_cast<double>(node);
    Moments m = moments_[node];
    for (std::size_t q = 0; q < 3; q++)
    {
        const double f = values_[node][q];
        const double df = values_[node + 1][q] - f;
        m[2 * q] += step_ * t * (f + 0.5 * df * t);
        m[2 * q + 1] += step_ * t * (z0 * f + 0.5 * t * (z0 * df + step_ * f) + step_ * df * t * t / 3.0);
    }
    return m;
}

std::array<double, 2> EffectiveSoundSpeed::wind(double z) const
{
    std::size_t node;
    double t;
    locate(z, node, t);
    return {values_[node][1] + t * (values_[node + 1][1] - values_[node][1]),
            values_[node][2] + t * (values_[node + 1][2] - values_[node][2])};
}

LinearProfile EffectiveSoundSpeed::fit(double z_a, const Moments& at_a, double z_b,
        const Moments& at_b, double hx, double hy) const
{
    double z1 = z_a;
    double z2 = z_b;
    const Moments* m1 = &at_a;
    const Moments* m2 = &at_b;
    if (z2 < z1)
    {
        std::swap(z1, z2);
        std::swap(m1, m2);
    }
    Moments widened;
    if (z2 - z1 < min_fit_span)
    {
        z2 = z1 + min_fit_span;
        widened = moments(z2);
        m2 = &widened;
    }
    // Least squares over [z1, z2]: slope = (integral of z f - z_mid integral of
    // f) / (span^3 / 12), and the line passes through the mean at z_mid.
    const double span = z2 - z1;
    const double z_mid = 0.5 * (z1 + z2);
    const double weight[3] = {1.0, hx, hy};
    double slope = 0.0;
    double mean = 0.0;
    for (std::size_t q = 0; q < 3; q++)
    {
        const double s0 = (*m2)[2 * q] - (*m1)[2 * q];
        const double s1 = (*m2)[2 * q + 1] - (*m1)[2 * q + 1];
        slope += weight[q] * (s1 - z_mid * s0) * 12.0 / (span * span * span);
        mean += weight[q] * s0 / span;
    }
    return {mean - slope * z_mid, slope};
}

RefractedRay refracted_ray(const EffectiveSoundSpeed& air, const Vec3& source,
        const EffectiveSoundSpeed::Moments& at_source, const Vec3& observer,
        double receiver_height, const EffectiveSoundSpeed::Moments& at_receiver,
        GroundGeometry ground)
{
    const double z_r = ray_receiver_height(observer, receiver_height);
    const double z_s = source.z;
    const Vec3 separation = observer - source;
    const double r = separation.horizontal_norm();
    const double hx = r > 0.0 ? separation.x / r : 1.0;
    const double hy = r > 0.0 ? separation.y / r : 0.0;

    RefractedRay ray;
    ray.profile = air.fit(z_s, at_source, z_r, at_receiver, hx, hy);

    if (r < 1e-3 || std::fabs(ray.profile.gradient) < straight_gradient)
    {
        // Straight, to the same receiver the curved rays go to.
        ray.path = straight_ray(source, Vec3(observer.x, observer.y, z_r));
        ray.path.separation = separation;
        ray.travel_time = ray.path.range
                / (ray.profile.sound_speed + ray.profile.gradient * 0.5 * (z_s + z_r));
        return ray;
    }

    const auto direct = direct_eigenray(ray.profile, z_s, z_r, r);
    RayPath& path = ray.path;
    path.separation = separation;
    const double horizontal = std::cos(direct.launch);
    path.direction = Vec3(horizontal * hx, horizontal * hy, std::sin(direct.launch));
    path.range = direct.spreading_range;
    path.path_length = direct.path_length;
    ray.travel_time = direct.travel_time;
    ray.shadow_depth = direct.shadow_depth;

    // The ground models take the straight-ray geometry that has the curved
    // ray's grazing angle and length: a source height h and ground distance d
    // with (h + z_r) = L sin(grazing) and d = L cos(grazing).  Straight rays
    // give back the real h and d.
    path.source_height = source.z;
    path.ground_distance = r;
    if (ground == GroundGeometry::point && z_s > 0.0)
    {
        // Reflected off the plane z = 0 to the microphone receiver_height above
        // it, where the ground models put it.
        const auto reflected = direct.lit ? reflected_eigenray(ray.profile, z_s, receiver_height, r)
                                          : ReflectedEigenray {};
        const double grazing = reflected.exists ? reflected.grazing : 0.0;
        const double length = reflected.exists ? reflected.path_length : direct.path_length;
        path.source_height = std::max(length * std::sin(grazing) - receiver_height, 0.0);
        path.ground_distance = length * std::cos(grazing);
    }
    else if (ground == GroundGeometry::plate)
    {
        // A flush microphone sees the direct ray arrive at its own angle.
        path.source_height = direct.path_length * std::sin(std::fabs(direct.arrival));
        path.ground_distance = direct.path_length * std::cos(direct.arrival);
    }
    return ray;
}

std::vector<double> shadow_attenuation_rates(const std::vector<double>& frequencies,
        std::size_t first_band, std::size_t last_band)
{
    constexpr double tau_1 = 1.01879297; // first zero of Ai'(-x)
    const double sin_60 = std::sqrt(3.0) / 2.0;
    const double nepers_to_db = 20.0 / std::log(10.0);
    std::vector<double> rates(frequencies.size(), 0.0);
    for (std::size_t band = first_band; band < std::min(last_band, frequencies.size()); band++)
        rates[band] = nepers_to_db * tau_1 * sin_60 * std::cbrt(pi * frequencies[band]);
    return rates;
}

void shadow_loss(const RefractedRay& ray, const std::vector<double>& rates, double floor,
        std::size_t first_band, std::size_t last_band, double* loss)
{
    const double a = std::fabs(ray.profile.gradient);
    const double scale = std::cbrt(a * a) / ray.profile.sound_speed * ray.shadow_depth;
    for (std::size_t band = first_band; band < last_band; band++)
        loss[band] = std::min(rates[band] * scale, floor);
}

} // namespace niceops
