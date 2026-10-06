#include "refraction.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "exact_solutions.hpp"
#include "ground_effect.hpp"
#include "plate_receiver.hpp"
#include "util/constants.hpp"
#include "util/parallel.hpp"

namespace niceops {

namespace {

// Gradients below this bend a ray by less than a nanoradian over 10 km: treated
// as straight, which also keeps the closed forms away from centers at infinity.
constexpr double straight_gradient = 1e-9; // 1/s

// The circles every ray lies on are centered on the line c = 0, at z = -c0/a.
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
        ray.boundary = t_s + t_r;
        ray.linear_boundary = ray.boundary;
        ray.boundary_gradient = profile.gradient;
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

    // The circle through (0, z_s) and (r, z_r) centered on the line c = 0 at
    // horizontal position x_c.  Equal distances from the center give
    // 2 r x_c = r^2 + h_r^2 - h_s^2.
    const double x_c = (r * r + c.sigma * (c.kappa(z_r) - c.kappa(z_s))) / (2.0 * r);
    const double rho = std::hypot(x_c, h_s);
    // Tangents, oriented from source to receiver.  On the upper half-circle
    // (downward refraction) travel is clockwise; on the lower, counterclockwise.
    ray.launch = std::atan2(c.sigma * x_c, h_s);
    ray.arrival = std::atan2(-c.sigma * (r - x_c), h_r);
    // The ray tube: r |dr/d(launch)| sin|arrival| / cos(launch) works out to
    // (r rho / h_s)^2 on any circle of the family.
    ray.spreading_range = r * rho / h_s;
    const double chord = std::hypot(r, z_s - z_r);
    ray.path_length = 2.0 * rho * std::asin(std::min(1.0, chord / (2.0 * rho)));
    // dt = ds / c = dphi / (|a| cos phi), phi the angle from the vertical through
    // the center, so t = atanh(sin phi) / |a| between the two ends.  Written as
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
    // theta, and at height z sits d(z) = sqrt(rho^2 - h(z)^2) from its center
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
    const double tolerance = 1e-10 * std::max(r, 1.0);

    // A launch angle gives a ray only if its circle reaches the receiver's height:
    // rho = h_s / cos(theta) >= h(z_r).  Below that v is negative, d() clamps it to
    // zero, and miss() is a continuation of the range that is not monotonic, whose
    // zeros are not rays -- a receiver above the source in a gradient that bends
    // sound down (a = 0.01, z_s = 2, z_r = 3, r = 800 m) came back with a launch
    // whose traced ray ended 1.1 m short of the receiver and a grazing angle 19%
    // low.  Launches with |theta| under `window` miss it.  That binds only where
    // sound bends down and the receiver is above the source (h grows with z): where
    // it bends up, h(z_r) > h_s means the receiver is below the source, and the
    // shadow's edge on the launch angle (`hi` below) is already the stricter limit.
    const double h_r = c.h(z_r);
    const double window = c.sigma > 0.0 && h_r > h_s ? std::acos(h_s / h_r) : 0.0;
    auto reaches = [&](double theta) { return std::fabs(theta) >= window * (1.0 - 1e-9); };

    // First, from a close guess: the direct ray to the receiver's image, z_r
    // below the ground, which the reflected ray becomes as z_r shrinks.  In
    // the medium continued below ground that is the circle through the two
    // points, as in direct_eigenray.  A secant step or two from there usually
    // closes to tolerance, against a dozen for the bracketed solve below
    // started cold -- which was most of what refraction cost with a point
    // receiver.  Anything that strays out of the domain falls back to it.
    {
        const double x_c = (r * r + c.sigma * (c.kappa(-z_r) - kappa_s)) / (2.0 * r);
        double theta0 = std::atan2(c.sigma * x_c, h_s);
        const double upper = c.sigma < 0.0 ? -std::asin(c.grazing_half_chord(z_s) / c.L)
                                           : 0.5 * pi - 1e-9;
        theta0 = std::clamp(theta0, -0.5 * pi + 1e-9, upper);
        double f0 = miss(theta0);
        double theta1 = std::clamp(theta0 + (f0 > 0.0 ? -1e-6 : 1e-6), -0.5 * pi + 1e-9, upper);
        double f1 = miss(theta1);
        for (int step = 0; step < 8 && std::fabs(f1) > tolerance && f1 != f0; step++)
        {
            const double theta2 = theta1 - f1 * (theta1 - theta0) / (f1 - f0);
            if (!(theta2 > -0.5 * pi && theta2 < upper))
                break;
            theta0 = theta1;
            f0 = f1;
            theta1 = theta2;
            f1 = miss(theta1);
        }
        if (std::fabs(f1) <= tolerance && reaches(theta1))
        {
            const auto l = legs(theta1);
            ray.exists = true;
            ray.launch = theta1;
            ray.grazing = std::atan2(l.d0, c.L);
            auto arc = [&](double horizontal, double vertical) {
                return 2.0 * l.rho
                        * std::asin(std::min(1.0, std::hypot(horizontal, vertical) / (2.0 * l.rho)));
            };
            ray.path_length = arc(l.to_ground, z_s) + arc(l.from_ground, z_r);
            return ray;
        }
    }

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
    // The angles that reach the receiver are two pieces either side of the window,
    // the range monotonic along each; the root is on the lower if the range at its
    // edge is already past r, else on the upper.
    if (window > 0.0 && lo < window && hi > -window)
    {
        if (lo < -window && miss(-window) >= 0.0)
            hi = -window;
        else if (hi > window)
            lo = window;
        else
            return ray;
    }
    const double image = -std::atan2(z_s + z_r, r);
    bool walk_up = false;
    if (image > lo && image < hi)
    {
        if (miss(image) < 0.0)
        {
            lo = image;
            walk_up = true;
        }
        else
            hi = image;
    }
    else if (lo == window && window > 0.0)
        walk_up = true;
    if (c.sigma > 0.0 && walk_up)
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

    if (!reaches(theta))
        return ray;
    const auto l = legs(theta);
    ray.exists = true;
    ray.launch = theta;
    // At the ground the ray is d(0) from its center horizontally and h(0) = L
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

    const CompassToLocal to_local(frame_bearing);

    values_.resize(nodes_);
    moments_.resize(nodes_);
    for (std::size_t node = 0; node < nodes_; node++)
    {
        const double z = z_low_ + step_ * static_cast<double>(node);
        const auto wind = to_local(profile.wind_at(z));
        values_[node] = {profile.at(z).sound_speed(), wind[0], wind[1]};
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

RefractedRay direct_ray(const EffectiveSoundSpeed& air, const Vec3& source,
        const EffectiveSoundSpeed::Moments& at_source, const Vec3& receiver,
        const EffectiveSoundSpeed::Moments& at_receiver)
{
    const double z_r = receiver.z;
    const double z_s = source.z;
    const Vec3 separation = receiver - source;
    const double r = separation.horizontal_norm();
    const double hx = r > 0.0 ? separation.x / r : 1.0;
    const double hy = r > 0.0 ? separation.y / r : 0.0;

    RefractedRay ray;
    ray.profile = air.fit(z_s, at_source, z_r, at_receiver, hx, hy);
    ray.source_z = z_s;
    ray.receiver_z = z_r;
    ray.horizontal = r;

    if (r < 1e-3 || std::fabs(ray.profile.gradient) < straight_gradient)
    {
        ray.path = straight_ray(source, receiver);
        ray.travel_time = ray.path.range
                / (ray.profile.sound_speed + ray.profile.gradient * 0.5 * (z_s + z_r));
        return ray;
    }

    ray.curved = true;
    ray.direct = direct_eigenray(ray.profile, z_s, z_r, r);
    RayPath& path = ray.path;
    path.separation = separation;
    const double horizontal = std::cos(ray.direct.launch);
    path.direction = Vec3(horizontal * hx, horizontal * hy, std::sin(ray.direct.launch));
    path.range = ray.direct.spreading_range;
    path.path_length = ray.direct.path_length;
    path.source_height = source.z;
    path.ground_distance = r;
    ray.travel_time = ray.direct.travel_time;
    ray.shadow_depth = ray.direct.shadow_depth;
    return ray;
}

void finish_ground(RefractedRay& ray, GroundGeometry ground, double microphone_height)
{
    if (!ray.curved)
        return;
    // The ground models take the straight-ray geometry that has the curved
    // ray's grazing angle and length: a source height h and ground distance d
    // with (h + microphone height) = L sin(grazing) and d = L cos(grazing).
    // Straight rays give back the real h and d.
    RayPath& path = ray.path;
    const DirectEigenray& direct = ray.direct;
    if (ground == GroundGeometry::point && ray.source_z > 0.0)
    {
        // Reflected off the plane z = 0 to the microphone microphone_height
        // above it, where the ground models put it.
        const auto reflected = direct.lit
                ? reflected_eigenray(ray.profile, ray.source_z, microphone_height, ray.horizontal)
                : ReflectedEigenray {};
        const double grazing = reflected.exists ? reflected.grazing : 0.0;
        const double length = reflected.exists ? reflected.path_length : direct.path_length;
        path.source_height = std::max(length * std::sin(grazing) - microphone_height, 0.0);
        path.ground_distance = length * std::cos(grazing);
    }
    else if (ground == GroundGeometry::plate)
    {
        // A flush microphone sees the direct ray arrive at its own angle.
        path.source_height = direct.path_length * std::sin(std::fabs(direct.arrival));
        path.ground_distance = direct.path_length * std::cos(direct.arrival);
    }
}

RefractedRay refracted_ray(const EffectiveSoundSpeed& air, const Vec3& source,
        const EffectiveSoundSpeed::Moments& at_source, const Vec3& receiver,
        const EffectiveSoundSpeed::Moments& at_receiver, GroundGeometry ground,
        double microphone_height)
{
    RefractedRay ray = direct_ray(air, source, at_source, receiver, at_receiver);
    finish_ground(ray, ground, microphone_height);
    return ray;
}

namespace {

constexpr double tau_1 = 1.01879297; // first zero of Ai'(-x)
const double nepers_to_db = 20.0 / std::log(10.0);
const double shadow_log_step = std::log(10.0) / ShadowAttenuation::per_decade;

} // namespace

std::vector<double> shadow_attenuation_rates(const std::vector<double>& frequencies,
        std::size_t first_band, std::size_t last_band)
{
    const double sin_60 = std::sqrt(3.0) / 2.0;
    std::vector<double> rates(frequencies.size(), 0.0);
    for (std::size_t band = first_band; band < std::min(last_band, frequencies.size()); band++)
        rates[band] = nepers_to_db * tau_1 * sin_60 * std::cbrt(pi * frequencies[band]);
    return rates;
}

ShadowAttenuation::ShadowAttenuation(const std::vector<double>& frequencies,
        std::size_t first_band, std::size_t last_band)
    : frequencies_(frequencies), rigid_rates_(shadow_attenuation_rates(frequencies, first_band, last_band))
{
}

ShadowAttenuation::ShadowAttenuation(const std::vector<double>& frequencies,
        const std::vector<std::complex<double>>& admittance, std::size_t first_band,
        std::size_t last_band)
    : frequencies_(frequencies), rigid_rates_(shadow_attenuation_rates(frequencies, first_band, last_band)),
      log_pi_f_(frequencies.size(), 0.0), table_(frequencies.size()),
      surface_wave_(frequencies.size(), false)
{
    if (admittance.size() != frequencies.size())
        throw std::invalid_argument("ShadowAttenuation: one admittance per band");
    const int nodes = static_cast<int>(std::ceil(std::log(s_max / s_min) / shadow_log_step)) + 1;
    std::vector<double> s(nodes);
    for (int i = 0; i < nodes; i++)
        s[i] = s_min * std::exp(i * shadow_log_step);
    const std::complex<double> w = std::polar(1.0, 2.0 * pi / 3.0);
    for (std::size_t band = first_band; band < std::min(last_band, frequencies.size()); band++)
    {
        log_pi_f_[band] = std::log(pi * frequencies[band]);
        const std::complex<double> direction = std::complex<double>(0.0, -1.0) * w * w * admittance[band];
        // Re q > 0 along the whole ray q = s d, arg beta below -60 degrees: the
        // surface wave's band, which takes the rigid rate.
        surface_wave_[band] = direction.real() > 0.0;
        if (surface_wave_[band])
            continue;
        const auto zeros = first_impedance_airy_zero(direction, s);
        table_[band].resize(zeros.size());
        for (std::size_t i = 0; i < zeros.size(); i++)
            table_[band][i] = (zeros[i] / w).imag();
    }
}

double ShadowAttenuation::factor_at(std::size_t band, double log_s) const
{
    const auto& g = table_[band];
    const double x = std::clamp((log_s - std::log(s_min)) / shadow_log_step, 0.0,
            static_cast<double>(g.size() - 1));
    const std::size_t i = std::min(static_cast<std::size_t>(x), g.size() - 2);
    return g[i] + (x - static_cast<double>(i)) * (g[i + 1] - g[i]);
}

double ShadowAttenuation::decay_factor(std::size_t band, double s) const
{
    if (rigid() || surface_wave_[band])
        return tau_1 * std::sqrt(3.0) / 2.0;
    return factor_at(band, std::log(s));
}

double ShadowAttenuation::rate(std::size_t band, const LinearProfile& profile) const
{
    const double a = std::fabs(profile.gradient);
    if (rigid() || surface_wave_[band])
        return rigid_rates_[band] * std::cbrt(a * a) / profile.sound_speed;
    // alpha = Im(b / w) pi f / (c0 s^2), s^2 = (pi f / |a|)^(2/3).
    const double log_s = (log_pi_f_[band] - std::log(a)) / 3.0;
    return nepers_to_db * factor_at(band, log_s) * pi * frequencies_[band]
            / (profile.sound_speed * std::exp(2.0 * log_s));
}

void shadow_loss(const RefractedRay& ray, const ShadowAttenuation& attenuation, double floor,
        std::size_t first_band, std::size_t last_band, double* loss)
{
    if (attenuation.rigid())
    {
        // As it was before there was a table, to the bit.
        const auto& rates = attenuation.rigid_rates();
        const double a = std::fabs(ray.profile.gradient);
        const double scale = std::cbrt(a * a) / ray.profile.sound_speed * ray.shadow_depth;
        for (std::size_t band = first_band; band < last_band; band++)
            loss[band] = std::min(rates[band] * scale, floor);
        return;
    }
    for (std::size_t band = first_band; band < last_band; band++)
        loss[band] = std::min(attenuation.rate(band, ray.profile) * ray.shadow_depth, floor);
}

} // namespace niceops

namespace niceops {

ShadowTransition::ShadowTransition(const std::vector<double>& frequencies,
        const std::vector<std::complex<double>>& admittance, double receiver_height,
        double sound_speed, std::size_t first_band, std::size_t last_band,
        const PlateReceiver* plate, const GroundModel* ground, unsigned threads)
    : frequencies_(frequencies), receiver_height_(receiver_height), admittance_(admittance),
      plate_(plate), sound_speed_(sound_speed)
{
    if (plate != nullptr && ground != nullptr)
    {
        ground_ = std::make_shared<GroundModel>(*ground);
        ground_->receiver_height = receiver_height;
    }
    if (plate != nullptr && ground == nullptr)
        throw std::invalid_argument("ShadowTransition: a plate needs its ground model");
    if (admittance.size() != frequencies.size())
        throw std::invalid_argument("ShadowTransition: one admittance per band");
    xi_nodes_ = static_cast<std::size_t>(std::lround((xi_max - xi_ray) / xi_step)) + 1;
    s_nodes_ = static_cast<std::size_t>(std::ceil(std::log10(s_max / s_min) * s_per_decade)) + 1;
    table_.assign(frequencies.size(), {});
    at_boundary_.assign(frequencies.size(), {});
    max_power_ratio_.assign(frequencies.size(), 1.0);
    const std::size_t at_zero = static_cast<std::size_t>(std::lround(-xi_ray / xi_step));
    const double doubling = 20.0 * std::log10(2.0);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::size_t end = std::min(last_band, frequencies.size());
    if (first_band >= end)
        return;

    // A few seconds' worth of mode sums in all, so the bands go in parallel.
    parallel_chunks(end - first_band, threads > 0 ? threads : default_thread_count(),
            [&](std::size_t lo, std::size_t hi) {
        for (std::size_t band = first_band + lo; band < first_band + hi; band++)
        {
            const double f = frequencies[band];
            const double k0 = 2.0 * pi * f / sound_speed;
            std::vector<double> col(xi_nodes_ * s_nodes_, nan);
            std::vector<double> zero(s_nodes_, nan);
            for (std::size_t j = 0; j < s_nodes_; j++)
            {
                const double s = s_min * std::pow(10.0, static_cast<double>(j) / s_per_decade);
                const double a = pi * f / (s * s * s);
                const double l = s / k0;
                const double diffraction = 2.0 * s * s / k0;
                const double y_ref = std::min(reference_height, 0.1 * s * s);
                const double z_ref = y_ref * l;
                // k^2 linear in height has parabolic rays; the model's are circles.
                // Each is placed by its own boundary.
                const double exact_boundary
                        = diffraction * (std::sqrt(y_ref) + std::sqrt(receiver_height / l));
                const LinearProfile profile {sound_speed, -a};
                const double model_boundary = shadow_boundary(profile, z_ref, receiver_height);
                std::vector<double> ranges;
                std::vector<std::size_t> which;
                for (std::size_t i = 0; i < xi_nodes_; i++)
                {
                    const double xi = xi_ray + i * xi_step;
                    if (exact_boundary + xi * diffraction > 0.0 && model_boundary + xi * diffraction > 0.0)
                    {
                        ranges.push_back(exact_boundary + xi * diffraction);
                        which.push_back(i);
                    }
                }
                if (ranges.empty())
                    continue;
                ExactProblem problem;
                problem.frequency = f;
                problem.sound_speed = sound_speed;
                problem.admittance = admittance[band];
                problem.source_height = z_ref;
                problem.receiver_height = receiver_height;
                std::vector<std::complex<double>> field;
                try
                {
                    field = upward_pressure(problem, sound_speed / a, ranges);
                }
                catch (const std::exception&)
                {
                    continue; // no column: this band keeps the ray model at this s
                }
                double* v = col.data() + j * xi_nodes_;
                for (std::size_t k = 0; k < ranges.size(); k++)
                {
                    const std::size_t i = which[k];
                    const double r_model = model_boundary + (xi_ray + i * xi_step) * diffraction;
                    const auto ray = direct_eigenray(profile, z_ref, receiver_height, r_model);
                    const double exact_db = 20.0 * std::log10(std::abs(field[k])
                            * std::hypot(ranges[k], z_ref - receiver_height));
                    const double spreading_db = 20.0
                            * std::log10(std::hypot(r_model, z_ref - receiver_height) / ray.spreading_range);
                    const double value = exact_db - spreading_db;
                    if (std::isfinite(value))
                        v[i] = std::min(value, doubling);
                }
                // The series is trusted everywhere in the shadow; should a node
                // there still come back NaN, carry the column on at its slope.
                for (std::size_t i = at_zero; i < xi_nodes_; i++)
                    if (std::isnan(v[i]) && i >= 2 && std::isfinite(v[i - 1]) && std::isfinite(v[i - 2]))
                        v[i] = 2.0 * v[i - 1] - v[i - 2];
                zero[j] = v[at_zero];
            }
            for (const double v : col)
                if (std::isfinite(v))
                    max_power_ratio_[band] = std::max(max_power_ratio_[band], std::pow(10.0, v / 10.0));
            table_[band] = std::move(col);
            at_boundary_[band] = std::move(zero);
        }
    });

    if (plate == nullptr)
        return;
    // K: the plate over bare ground at its top, 2 degrees up, a kilometer out.
    plate_offset_.assign(frequencies.size(), std::vector<double>(plate_azimuths, 0.0));
    GroundModel bare = *ground;
    bare.receiver_height = receiver_height;
    const double elevation = 2.0 * pi / 180.0;
    const double h = 1000.0 * std::sin(elevation), d = 1000.0 * std::cos(elevation);
    std::vector<double> bare_db(frequencies.size(), 0.0), plate_db(frequencies.size(), 0.0);
    excess_ground_attenuation(bare, admittance, h, d, sound_speed, frequencies, bare_db.data(),
            first_band, end, false);
    for (int i = 0; i < plate_azimuths; i++)
    {
        plate->level(h, d, 2.0 * pi * i / plate_azimuths, plate_db.data(), first_band, end);
        for (std::size_t band = first_band; band < end; band++)
            plate_offset_[band][i] = plate_db[band] - bare_db[band];
    }
    // The bound has to cover V + K at the loudest azimuth.
    for (std::size_t band = first_band; band < end; band++)
    {
        const double k_max = *std::max_element(plate_offset_[band].begin(), plate_offset_[band].end());
        max_power_ratio_[band] = 1.0;
        for (const double v : table_[band])
            if (std::isfinite(v))
                max_power_ratio_[band] = std::max(max_power_ratio_[band], std::pow(10.0, (v + k_max) / 10.0));
    }
}

double ShadowTransition::plate_offset(std::size_t band, double azimuth) const
{
    if (plate_offset_.empty() || plate_offset_[band].empty())
        return 0.0;
    const double x = std::fmod(azimuth / (2.0 * pi) * plate_azimuths, static_cast<double>(plate_azimuths));
    const double wrapped = x < 0.0 ? x + plate_azimuths : x;
    const int i = std::min(static_cast<int>(wrapped), plate_azimuths - 1);
    const double t = wrapped - i;
    const auto& k = plate_offset_[band];
    return k[i] + t * (k[(i + 1) % plate_azimuths] - k[i]);
}

double ShadowTransition::column(std::size_t band, std::size_t j, double xi) const
{
    const double* v = table_[band].data() + j * xi_nodes_;
    const double x = (xi - xi_ray) / xi_step;
    if (x < 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    const std::size_t last = xi_nodes_ - 1;
    if (x >= static_cast<double>(last))
        return v[last] + (x - static_cast<double>(last)) * (v[last] - v[last - 1]);
    const std::size_t i = static_cast<std::size_t>(x);
    return v[i] + (x - static_cast<double>(i)) * (v[i + 1] - v[i]);
}

double ShadowTransition::level(std::size_t band, double xi, double s) const
{
    if (table_[band].empty())
        return std::numeric_limits<double>::quiet_NaN();
    const double x = std::clamp(std::log10(s / s_min) * s_per_decade, 0.0,
            static_cast<double>(s_nodes_ - 1));
    const std::size_t j = std::min(static_cast<std::size_t>(x), s_nodes_ - 2);
    const double t = x - static_cast<double>(j);
    const double v0 = column(band, j, xi);
    const double v1 = column(band, j + 1, xi);
    // Held at V(0) - floor happens in apply(); here only the table.
    return v0 + t * (v1 - v0);
}

double ShadowScattering::level(double frequency, double slant, double source_height) const
{
    const double range = std::clamp(slant, 100.0, 10000.0);
    const double height = std::clamp(source_height, 10.0, 300.0);
    return std::min(0.0, a_db + 10.0 * std::log10(gamma_t)
            + per_decade_frequency * std::log10(frequency / 1000.0)
            + per_decade_range * std::log10(range / 100.0) + per_decade_height * std::log10(height / 100.0));
}

ShadowScattering ShadowScattering::for_receiver(double gamma_t, double receiver_height, bool plate)
{
    // Fitted to the PE with gamma_t = 3e-5 and L0 = 10 m through measured
    // day-284 upwind profiles (LIDAR and balloon), 125 Hz - 2 kHz, over
    // variable porosity 200, at the points deep in each shadow where the
    // turbulent field is 10 dB clear of the still one; a third of that
    // strength moved them by 4.82 dB against single scattering's 4.77.
    ShadowScattering s;
    s.gamma_t = gamma_t;
    if (plate || receiver_height < ShadowTransition::flush_height)
    {
        // A microphone 3 mm up; sources 50 and 120 m up, 12 profiles, 9320
        // points, 1.5 dB rms about the fit.
        s.a_db = 39.18;
        s.per_decade_frequency = -1.55;
        s.per_decade_range = 2.10;
        s.per_decade_height = 5.30;
    }
    else
    {
        // 1.2 m up, sources 50 m up, 6 profiles, 2703 points where stratified
        // is 25 dB down.  The ground's interference with the scattered field
        // tilts it up with frequency.  Only one source height, so the height
        // term is the flush microphone's.
        s.a_db = 36.67;
        s.per_decade_frequency = 5.61;
        s.per_decade_range = 5.13;
        s.per_decade_height = 5.30;
    }
    return s;
}

double ShadowTransition::band_boundary(double r_s, double r_l, double z_s, double z_r,
        double frequency, double sound_speed, bool flush)
{
    if (!std::isfinite(r_l) && std::isfinite(r_s))
        return r_s * (flush ? fallback_multiple_flush : fallback_multiple_raised)
                * std::pow(frequency / 1000.0, fallback_per_decade);
    if (!(r_l > r_s))
        return r_s;
    const double root = std::sqrt(2.0 * std::max(z_s, 0.0)) + std::sqrt(2.0 * std::max(z_r, 0.0));
    const double a_linear = sound_speed * root * root / (r_l * r_l);
    // x: the band's diffraction length at the fit's curvature over the gap.
    const double s = std::cbrt(pi * frequency / a_linear);
    const double x = s * s * sound_speed / (pi * frequency) / (r_l - r_s);
    const double resolution = flush ? boundary_resolution_flush : boundary_resolution;
    return r_s + x / (x + resolution) * (r_l - r_s);
}

bool ShadowTransition::apply(const RefractedRay& ray, const ShadowAttenuation& attenuation,
        double floor, std::size_t first_band, std::size_t last_band, double* ground_weight,
        double* loss, double plate_azimuth) const
{
    if (empty() || !ray.curved || !(ray.direct.boundary_gradient < 0.0)
            || !std::isfinite(ray.direct.boundary) || first_band >= last_band)
        return false;
    const double c0 = ray.profile.sound_speed;
    const double r_s = ray.direct.boundary;
    const double r_l = ray.direct.linear_boundary;
    // Two boundaries (through layers); the layered one alone, where the fit
    // casts no shadow (a multiple of it, band_boundary); or one (a linear
    // profile).
    const bool between = std::isfinite(r_l) && r_l > r_s;
    const bool fallback = !std::isfinite(r_l);
    const bool flush = !plate_offset_.empty() || receiver_height_ < flush_height;
    const double root = std::sqrt(2.0 * std::max(ray.source_z, 0.0))
            + std::sqrt(2.0 * std::max(ray.receiver_z, 0.0));
    // The curvature that makes a boundary r: c0 root^2 / r^2.
    const double a_one = -ray.direct.boundary_gradient;
    const double r_weakest = fallback
            ? band_boundary(r_s, r_l, ray.source_z, ray.receiver_z, frequencies_[first_band], c0, flush)
            : r_l;
    const double a_linear = between || fallback ? c0 * root * root / (r_weakest * r_weakest) : a_one;
    // Every band's boundary lies in [r_s, r_l] and its diffraction length is at
    // most the lowest band's at the weakest curvature: lit by more than xi_ray
    // of that, it is lit in every band.
    const double f_lo = frequencies_[first_band];
    const double s_lo = std::cbrt(pi * f_lo / a_linear);
    if (ray.horizontal - r_s < xi_ray * s_lo * s_lo * c0 / (pi * f_lo) && !(ray.shadow_depth > 0.0))
        return false;
    const LinearProfile one_boundary {c0, -a_one};
    for (std::size_t band = first_band; band < last_band; band++)
    {
        const double f = frequencies_[band];
        double a = a_one;
        double boundary = r_s;
        double depth = ray.shadow_depth;
        if (between || fallback)
        {
            boundary = band_boundary(r_s, r_l, ray.source_z, ray.receiver_z, f, c0, flush);
            a = c0 * root * root / (boundary * boundary);
            depth = std::max(ray.horizontal - boundary, 0.0);
        }
        const LinearProfile equivalent = between || fallback ? LinearProfile {c0, -a} : one_boundary;
        const double creeping = depth > 0.0
                ? std::min(attenuation.rate(band, equivalent) * depth, floor)
                : 0.0;
        ground_weight[band] = 1.0;
        loss[band] = creeping;
        const double s = std::cbrt(pi * f / a);
        if (s < s_min)
            continue;
        const double xi = (ray.horizontal - boundary) * pi * f / (s * s * c0); // F = 2 s^2 / k0
        if (xi <= xi_ray)
            continue;
        // A boundary nearer the source than the ray model's own reach, -xi_ray
        // diffraction lengths, leaves no lit zone for the table to hand over
        // from: the source is then a few l up at most, under the 5 l the table's
        // collapse over source height was checked down to.  Such a band keeps
        // the ray model up to its boundary and blends to the table over the
        // diffraction length past it.  A source 10 m up under a weak upwind
        // gradient has it at every band to 500 Hz, and there V put the whole
        // path 2-5 dB under the parabolic equation, which matched straight
        // rays.  Keeping the ray model into the shadow as well was worse deep
        // in it (50-160 Hz 5.7 dB rms against the PE on the day-286 matrix,
        // against 5.1 this way) and through 200-630 Hz near it (4.0 against
        // 2.4); so was the ray model with only V's fall past the boundary.  A
        // 30 m source in a linear lapse comes out 0.1 dB rms worse at
        // 200-630 Hz (the old lit side was 1.3 dB low, the old shadow right).
        const bool cramped = boundary * pi * f / (s * s * c0) < -xi_ray;
        if (cramped && xi <= 0.0)
            continue;
        const double blend_from = cramped ? 0.0 : xi_ray;
        const double blend_to = cramped ? 1.0 : xi_modal;
        double v = level(band, xi, s);
        if (std::isnan(v))
            continue;
        // At worst V(0) - floor: the boundary's own level, less the floor.
        const double x = std::clamp(std::log10(s / s_min) * s_per_decade, 0.0,
                static_cast<double>(s_nodes_ - 1));
        const std::size_t j = std::min(static_cast<std::size_t>(x), s_nodes_ - 2);
        const double v_zero = at_boundary_[band][j]
                + (x - static_cast<double>(j)) * (at_boundary_[band][j + 1] - at_boundary_[band][j]);
        if (std::isfinite(v_zero) && !scattering_.on())
            v = std::max(v, v_zero - floor);
        v += plate_offset(band, plate_azimuth);
        if (scattering_.on())
        {
            // The scattered energy, from re free field at the slant range into
            // re the ray's own spreading, as V is; the plate lifts it as it
            // lifts the coherent field.  The sum is held to the table's largest,
            // which is what the peak search's bound allows for.
            const double slant = std::hypot(ray.horizontal, ray.source_z - ray.receiver_z);
            const double spreading_db = 20.0 * std::log10(slant / ray.path.range);
            const double scattered = scattering_.level(f, slant, ray.source_z)
                    + plate_offset(band, plate_azimuth) - spreading_db;
            v = std::min(10.0 * std::log10(std::pow(10.0, 0.1 * v) + std::pow(10.0, 0.1 * scattered)),
                    10.0 * std::log10(max_power_ratio_[band]));
        }
        if (xi >= blend_to)
        {
            ground_weight[band] = 0.0;
            loss[band] = -v;
        }
        else
        {
            const double t = (xi - blend_from) / (blend_to - blend_from);
            ground_weight[band] = 1.0 - t;
            // A cramped band blends from the ray model at its boundary, where
            // the creeping-wave loss, valid only deep in the shadow, has no
            // place yet.
            loss[band] = (1.0 - t) * (cramped ? 0.0 : creeping) - t * v;
        }
    }
    return true;
}

void ShadowTransition::set_cramped_ground(CrampedGround mode, const GroundModel* ground, unsigned threads)
{
    cramped_ground_ = mode;
    if (ground != nullptr && plate_ == nullptr)
    {
        ground_ = std::make_shared<GroundModel>(*ground);
        ground_->receiver_height = receiver_height_;
    }
    if (mode == CrampedGround::straight && ground_ == nullptr)
        throw std::invalid_argument("ShadowTransition: a straight cramped ground term needs the ground model");
    if (mode != CrampedGround::series || !cramped_table_.empty() || empty())
        return;
    // The residue series at the source's own height, sqrt(y_s) on its nodes,
    // over the same (xi, s) grid as V: k^2 linear in height with the
    // curvature whose model boundary is the band's.  One series per (band,
    // s, y_s), each a few milliseconds.
    cramped_table_.assign(frequencies_.size(), {});
    std::vector<std::size_t> work;
    for (std::size_t band = 0; band < frequencies_.size(); band++)
        if (!table_[band].empty())
            work.push_back(band);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::size_t per_band = s_nodes_ * cramped_root_nodes;
    cramped_xi_nodes_ = static_cast<std::size_t>(std::lround((xi_max - xi_ray) / cramped_xi_step)) + 1;
    const std::size_t xi_nodes = cramped_xi_nodes_;
    std::vector<std::vector<float>> tables(frequencies_.size());
    for (const std::size_t band : work)
        tables[band].assign(per_band * xi_nodes, static_cast<float>(nan));
    parallel_chunks(work.size() * s_nodes_, threads > 0 ? threads : default_thread_count(),
            [&](std::size_t lo, std::size_t hi) {
        for (std::size_t item = lo; item < hi; item++)
        {
            const std::size_t band = work[item / s_nodes_];
            const std::size_t j = item % s_nodes_;
            const double f = frequencies_[band];
            const double k0 = 2.0 * pi * f / sound_speed_;
            const double s = s_min * std::pow(10.0, static_cast<double>(j) / s_per_decade);
            const double a = pi * f / (s * s * s);
            const double l = s / k0;
            const double diffraction = 2.0 * s * s / k0;
            const LinearProfile profile {sound_speed_, -a};
            // Every height's ranges in one list: the modes and H0 are shared.
            std::vector<double> heights, ranges;
            std::vector<std::size_t> node_of, height_of, xi_of;
            for (int k = 0; k < cramped_root_nodes; k++)
            {
                const double root = cramped_root_max * static_cast<double>(k) / (cramped_root_nodes - 1);
                // Below the receiver too: the series is reciprocal, and the
                // nodes bracket every source above it.
                const double z_s = root * root * l;
                const double model_boundary = shadow_boundary(profile, z_s, receiver_height_);
                for (std::size_t i = 0; i < xi_nodes; i++)
                {
                    const double r = model_boundary + (xi_ray + i * cramped_xi_step) * diffraction;
                    if (r > 0.0)
                    {
                        ranges.push_back(r);
                        height_of.push_back(heights.size());
                        node_of.push_back(static_cast<std::size_t>(k));
                        xi_of.push_back(i);
                    }
                }
                heights.push_back(z_s);
            }
            if (ranges.empty())
                continue;
            ExactProblem problem;
            problem.frequency = f;
            problem.sound_speed = sound_speed_;
            problem.admittance = admittance_[band];
            problem.receiver_height = receiver_height_;
            std::vector<std::vector<std::complex<double>>> field;
            try
            {
                field = upward_pressure_heights(problem, sound_speed_ / a, heights, ranges);
            }
            catch (const std::exception&)
            {
                continue;
            }
            for (std::size_t q = 0; q < ranges.size(); q++)
            {
                const double z_s = heights[height_of[q]];
                const double db = 20.0 * std::log10(std::abs(field[height_of[q]][q])
                        * std::hypot(ranges[q], z_s - receiver_height_));
                if (std::isfinite(db))
                    tables[band][(j * cramped_root_nodes + node_of[q]) * xi_nodes + xi_of[q]] = static_cast<float>(db);
            }
        }
    });
    cramped_table_ = std::move(tables);
}

double ShadowTransition::cramped_level(std::size_t band, double xi, double s, double y_s) const
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (band >= cramped_table_.size() || cramped_table_[band].empty() || !(y_s >= 0.0))
        return nan;
    const double y = std::sqrt(y_s) / cramped_root_max * (cramped_root_nodes - 1);
    if (y > cramped_root_nodes - 1)
        return nan;
    const std::size_t k = std::min(static_cast<std::size_t>(y), static_cast<std::size_t>(cramped_root_nodes - 2));
    const double ty = y - static_cast<double>(k);
    const double x = std::clamp(std::log10(s / s_min) * s_per_decade, 0.0, static_cast<double>(s_nodes_ - 1));
    const std::size_t j = std::min(static_cast<std::size_t>(x), s_nodes_ - 2);
    const double tx = x - static_cast<double>(j);
    const double u = (xi - xi_ray) / cramped_xi_step;
    if (u < 0.0)
        return nan;
    const std::size_t last = cramped_xi_nodes_ - 1;
    auto column = [&](std::size_t jj, std::size_t kk) {
        const float* v = cramped_table_[band].data() + (jj * cramped_root_nodes + kk) * cramped_xi_nodes_;
        if (u >= static_cast<double>(last))
            return v[last] + (u - static_cast<double>(last)) * (static_cast<double>(v[last]) - v[last - 1]);
        const std::size_t i = static_cast<std::size_t>(u);
        return v[i] + (u - static_cast<double>(i)) * (static_cast<double>(v[i + 1]) - v[i]);
    };
    const double c00 = column(j, k), c01 = column(j, k + 1), c10 = column(j + 1, k), c11 = column(j + 1, k + 1);
    // A corner with no value (NaN) makes the whole lookup NaN: the caller falls back.
    return (1.0 - tx) * ((1.0 - ty) * c00 + ty * c01) + tx * ((1.0 - ty) * c10 + ty * c11);
}

bool ShadowTransition::apply_resolved(const RefractedRay& ray, const ResolvedBand* bands,
        const ShadowAttenuation& attenuation, double floor, std::size_t first_band, std::size_t last_band,
        double* ground_weight, double* loss, double plate_azimuth, const double* cramped_field) const
{
    if (empty() || !ray.curved || !(ray.direct.boundary_gradient < 0.0)
            || !std::isfinite(ray.direct.boundary) || first_band >= last_band)
        return false;
    const double c0 = ray.profile.sound_speed;
    const double r = ray.horizontal;
    const double root = std::sqrt(2.0 * std::max(ray.source_z, 0.0))
            + std::sqrt(2.0 * std::max(ray.receiver_z, 0.0));
    const double slant = std::hypot(r, ray.source_z - ray.receiver_z);
    const double spreading_db = 20.0 * std::log10(slant / ray.path.range);
    // The ray's equivalent geometry has collapsed when it is the limiting ray
    // carried into a shadow (the layered one past r_s, or the fit's past r_l).
    const bool degenerate = !ray.direct.lit
            || (ray.from_linear_fit && !(r <= ray.direct.linear_boundary));
    // The straight path's ground term, dB, computed once if any band asks.
    thread_local std::vector<double> straight_db;
    bool have_straight = false;
    auto straight = [&](std::size_t band) {
        if (!have_straight)
        {
            straight_db.assign(frequencies_.size(), 0.0);
            if (plate_ != nullptr)
                plate_->level(ray.source_z, r, plate_azimuth, straight_db.data(), first_band, last_band);
            else
                excess_ground_attenuation(*ground_, admittance_, ray.source_z, r, sound_speed_, frequencies_,
                        straight_db.data(), first_band, last_band, false);
            have_straight = true;
        }
        return std::max(straight_db[band], -200.0);
    };
    const bool use_straight = cramped_ground_ == CrampedGround::straight && ground_ != nullptr;
    for (std::size_t band = first_band; band < last_band; band++)
    {
        const double f = frequencies_[band];
        ground_weight[band] = 1.0;
        loss[band] = 0.0;
        const double boundary = bands[band].boundary;
        if (!std::isfinite(boundary) || !(boundary > 0.0))
        {
            // No shadow in this band: lit, at the ray's own level -- unless the
            // ray is a shadow ray, whose ground geometry is the grazing one.
            if (degenerate && use_straight)
            {
                ground_weight[band] = 0.0;
                loss[band] = -straight(band);
            }
            continue;
        }
        const double a = tangent_decay_ && bands[band].gradient > 0.0 ? bands[band].gradient
                                                                      : c0 * root * root / (boundary * boundary);
        const double depth = std::max(r - boundary, 0.0);
        const double creeping = depth > 0.0
                ? std::min(attenuation.rate(band, LinearProfile {c0, -a}) * depth, floor)
                : 0.0;
        loss[band] = creeping;
        const double s = std::cbrt(pi * f / a);
        const double diffraction = s * s * c0 / (pi * f); // F = 2 s^2 / k0
        const double xi = (r - boundary) / diffraction;
        const bool cramped = boundary / diffraction < -xi_ray;
        // The table's level at the boundary, for the floor.
        auto at_zero = [&]() {
            const double x = std::clamp(std::log10(s / s_min) * s_per_decade, 0.0, static_cast<double>(s_nodes_ - 1));
            const std::size_t j = std::min(static_cast<std::size_t>(x), s_nodes_ - 2);
            return at_boundary_[band][j] + (x - static_cast<double>(j)) * (at_boundary_[band][j + 1] - at_boundary_[band][j]);
        };
        const double field = !(cramped && cramped_ground_ == CrampedGround::series && s >= s_min)
                ? std::numeric_limits<double>::quiet_NaN()
                : cramped_field != nullptr ? cramped_field[band]
                                           : cramped_level(band, xi, s, ray.source_z * pi * f / (s * c0 / 2.0));
        if (std::isfinite(field))
        {
            // The field itself, re the ray's spreading, held at worst to the
            // floor below the table's boundary level and at best to the
            // table's largest, which the peak search's bound allows for.
            double v = field + spreading_db;
            const double v_zero = at_zero();
            if (std::isfinite(v_zero) && !scattering_.on())
                v = std::max(v, v_zero - floor);
            v += plate_offset(band, plate_azimuth);
            if (scattering_.on())
            {
                const double scattered = scattering_.level(f, slant, ray.source_z)
                        + plate_offset(band, plate_azimuth) - spreading_db;
                v = 10.0 * std::log10(std::pow(10.0, 0.1 * v) + std::pow(10.0, 0.1 * scattered));
            }
            ground_weight[band] = 0.0;
            loss[band] = -std::min(v, 10.0 * std::log10(max_power_ratio_[band]));
            continue;
        }
        if (s < s_min || xi <= xi_ray)
        {
            if (degenerate && use_straight && xi <= 0.0)
            {
                ground_weight[band] = 0.0;
                loss[band] = -straight(band);
            }
            continue;
        }
        if (cramped && xi <= 0.0)
        {
            if (use_straight)
            {
                ground_weight[band] = 0.0;
                loss[band] = -straight(band);
            }
            continue;
        }
        const double blend_from = cramped ? 0.0 : xi_ray;
        const double blend_to = cramped ? 1.0 : xi_modal;
        double v = level(band, xi, s);
        if (std::isnan(v))
            continue;
        const double v_zero = at_zero();
        if (std::isfinite(v_zero) && !scattering_.on())
            v = std::max(v, v_zero - floor);
        v += plate_offset(band, plate_azimuth);
        if (scattering_.on())
        {
            const double scattered = scattering_.level(f, slant, ray.source_z) + plate_offset(band, plate_azimuth)
                    - spreading_db;
            v = std::min(10.0 * std::log10(std::pow(10.0, 0.1 * v) + std::pow(10.0, 0.1 * scattered)),
                    10.0 * std::log10(max_power_ratio_[band]));
        }
        if (xi >= blend_to)
        {
            ground_weight[band] = 0.0;
            loss[band] = -v;
            continue;
        }
        const double t = (xi - blend_from) / (blend_to - blend_from);
        if (cramped && use_straight)
        {
            ground_weight[band] = 0.0;
            loss[band] = -(1.0 - t) * straight(band) - t * v;
        }
        else
        {
            ground_weight[band] = 1.0 - t;
            loss[band] = (1.0 - t) * (cramped ? 0.0 : creeping) - t * v;
        }
    }
    return true;
}

} // namespace niceops
