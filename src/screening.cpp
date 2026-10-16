#include "screening.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>

#include "util/constants.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

// Twice the signed area of (o, a, b): positive when b lies to the left of the
// directed line o -> a, which for a increasing in s is above it.
double cross(const ScreenPoint& o, const ScreenPoint& a, const ScreenPoint& b)
{
    return (a.s - o.s) * (b.z - o.z) - (a.z - o.z) * (b.s - o.s);
}

double length(const ScreenPoint& a, const ScreenPoint& b)
{
    return std::hypot(b.s - a.s, b.z - a.z);
}

// `p` mirrored in the line z = intercept + slope s.
ScreenPoint mirror(const ScreenPoint& p, const GroundLine& line)
{
    const double norm = std::sqrt(1.0 + line.slope * line.slope);
    const double above = (p.z - (line.intercept + line.slope * p.s)) / norm;
    return {p.s + 2.0 * above * line.slope / norm, p.z - 2.0 * above / norm};
}

} // namespace

std::size_t niceops::upper_hull(const ScreenPoint* points, std::size_t count, std::size_t* hull)
{
    // Andrew's monotone chain, upper half: a point that does not turn the
    // chain clockwise -- collinear included -- is dropped.
    std::size_t n = 0;
    for (std::size_t i = 0; i < count; i++)
    {
        while (n >= 2 && cross(points[hull[n - 2]], points[hull[n - 1]], points[i]) >= 0.0)
            n--;
        hull[n++] = i;
    }
    return n;
}

double niceops::path_difference(const ScreenPoint& from, const ScreenPoint* edges, std::size_t m,
        const ScreenPoint& to)
{
    assert(m >= 1);
    double path = length(from, edges[0]);
    for (std::size_t i = 0; i + 1 < m; i++)
        path += length(edges[i], edges[i + 1]);
    path += length(edges[m - 1], to);
    const double delta = path - length(from, to);
    bool above = false;
    for (std::size_t i = 0; i < m && !above; i++)
        above = cross(from, to, edges[i]) > 0.0;
    return above ? delta : -delta;
}

double niceops::cnossos_dif(double delta, double e, std::size_t edges, double wavelength)
{
    double c = 1.0;
    if (edges > 1 && e > 0.3)
    {
        const double r = 5.0 * wavelength / e;
        c = (1.0 + r * r) / (1.0 / 3.0 + r * r);
    }
    const double x = 40.0 / wavelength * c * delta;
    if (!(x >= -2.0))
        return 0.0;
    return 10.0 * std::log10(cnossos_floor + x);
}

bool niceops::rubber_band(const std::vector<double>& heights, std::size_t interior, double step, double distance,
        double microphone_height, double source_height, std::vector<ScreenPoint>& points,
        std::vector<std::size_t>& hull, std::vector<ScreenPoint>& edges, Screen& screen)
{
    assert(heights.size() >= interior + 1);
    const std::size_t count = interior + 2;
    points.resize(count);
    hull.resize(count);
    points[0] = {0.0, microphone_height};
    for (std::size_t k = 1; k <= interior; k++)
        points[k] = {static_cast<double>(k) * step, heights[k]};
    points[count - 1] = {distance, source_height};
    const std::size_t vertices = upper_hull(points.data(), count, hull.data());
    // The ends are always vertices; a blocked path has at least one between.
    if (vertices <= 2)
        return false;
    edges.clear();
    for (std::size_t i = 1; i + 1 < vertices; i++)
        edges.push_back(points[hull[i]]);
    screen = Screen {};
    screen.receiver = points[0];
    screen.source = points[count - 1];
    screen.edges = edges.size();
    screen.receiver_edge = hull[1];
    screen.source_edge = hull[vertices - 2];
    screen.e = 0.0;
    for (std::size_t i = 0; i + 1 < edges.size(); i++)
        screen.e += length(edges[i], edges[i + 1]);
    screen.delta = path_difference(screen.receiver, edges.data(), edges.size(), screen.source);
    return true;
}

void niceops::screen_sides(Screen& screen, const std::vector<ScreenPoint>& edges, const GroundLine& source_line,
        const GroundLine& receiver_line, bool flush)
{
    assert(!edges.empty() && edges.size() == screen.edges);
    // Source side: the source over the line, with the edge nearest it as the
    // receiver.  plane_geometry measures from its "microphone" at s = 0, so the
    // line is moved to start under the edge.
    const ScreenPoint& near_source = edges.back();
    const GroundLine from_edge {source_line.intercept + source_line.slope * near_source.s, source_line.slope};
    const PlaneGeometry s = plane_geometry(screen.source.s - near_source.s, screen.source.z, near_source.z, from_edge);
    screen.source_side = {std::max(0.0, s.source_height), std::max(0.0, s.receiver_height),
            std::fabs(s.ground_distance)};
    screen.delta_source_image = path_difference(screen.receiver, edges.data(), edges.size(),
            mirror(screen.source, source_line));

    // Receiver side: the edge nearest the microphone over the line, as the
    // source, and the microphone.
    const ScreenPoint& near_receiver = edges.front();
    const PlaneGeometry r = plane_geometry(near_receiver.s, near_receiver.z, screen.receiver.z, receiver_line);
    if (flush)
    {
        screen.receiver_side = {std::max(0.0, r.source_height - r.receiver_height), 0.0, std::fabs(r.ground_distance)};
        screen.delta_receiver_image = screen.delta;
    }
    else
    {
        screen.receiver_side = {std::max(0.0, r.source_height), std::max(0.0, r.receiver_height),
                std::fabs(r.ground_distance)};
        screen.delta_receiver_image = path_difference(mirror(screen.receiver, receiver_line), edges.data(),
                edges.size(), screen.source);
    }
}

void niceops::cnossos_screen_loss(const Screen& screen, const std::vector<double>& band_centers,
        double sound_speed, const double* source_gain, const double* receiver_gain, std::size_t first,
        std::size_t last, double* loss, ScreenTally* tally)
{
    assert((source_gain == nullptr) == (receiver_gain == nullptr));
    assert(source_gain == nullptr
            || (std::isfinite(screen.delta_source_image) && std::isfinite(screen.delta_receiver_image)));
    last = std::min(last, band_centers.size());
    for (std::size_t band = first; band < last; band++)
    {
        const double wavelength = sound_speed / band_centers[band];
        const double dif = cnossos_dif(screen.delta, screen.e, screen.edges, wavelength);
        double a = std::min(dif, cnossos_cap);
        if (tally != nullptr && dif > cnossos_cap)
            tally->capped = true;
        if (source_gain != nullptr)
        {
            // Each side's ground, weighted by how much deeper in the shadow
            // its image lies; the Delta_dif terms here are not capped.
            const auto side = [&](double image_delta, double gain, bool* clamped) {
                double weight = std::pow(10.0,
                        -(cnossos_dif(image_delta, screen.e, screen.edges, wavelength) - dif) / 20.0);
                if (weight > 1.0)
                {
                    weight = 1.0;
                    if (tally != nullptr)
                        *clamped = true;
                }
                return -20.0 * std::log10(1.0 + (std::sqrt(gain) - 1.0) * weight);
            };
            bool unused = false;
            a += side(screen.delta_source_image, source_gain[band],
                    tally != nullptr ? &tally->source_clamped : &unused);
            a += side(screen.delta_receiver_image, receiver_gain[band],
                    tally != nullptr ? &tally->receiver_clamped : &unused);
        }
        loss[band] = a;
    }
}

KnifeEdge niceops::knife_edge(const Screen& screen, const std::vector<ScreenPoint>& points)
{
    assert(points.size() >= 3 && screen.edges >= 1);
    const ScreenPoint& r = screen.receiver;
    const ScreenPoint& s = screen.source;
    const double distance = length(r, s);
    assert(distance > 0.0);
    // The line of sight's direction t from the microphone, and n, normal to it
    // and up (the profile's s increases toward the source).
    const double ts = (s.s - r.s) / distance;
    const double tz = (s.z - r.z) / distance;
    const double ns = -tz;
    const double nz = ts;
    const auto along = [&](const ScreenPoint& p) { return (p.s - r.s) * ts + (p.z - r.z) * tz; };
    const auto above = [&](const ScreenPoint& p) { return (p.s - r.s) * ns + (p.z - r.z) * nz; };
    // The aperture plane: through the line of sight where the profile rises
    // highest above it (the first such sample on a tie).
    std::size_t deepest = 1;
    for (std::size_t k = 2; k + 1 < points.size(); k++)
        if (above(points[k]) > above(points[deepest]))
            deepest = k;
    const double plane = along(points[deepest]);
    assert(plane > 0.0 && plane < distance);
    // Where each horizon ray crosses the plane, above the line of sight: the
    // microphone's through the vertex nearest it, the source's through the
    // vertex nearest the source.  Both vertices lie above the line, so both
    // crossings do.
    const ScreenPoint& near_receiver = points[screen.receiver_edge];
    const ScreenPoint& near_source = points[screen.source_edge];
    const double receiver_ray = above(near_receiver) * plane / along(near_receiver);
    const double source_ray = above(near_source) * (distance - plane) / (distance - along(near_source));
    KnifeEdge edge;
    edge.d1 = distance - plane;
    edge.d2 = plane;
    edge.height = std::max(receiver_ray, source_ray);
    // sqrt(d^2 + h^2) - d, written so that it keeps its digits for h << d.
    const auto extra = [&](double d) {
        const double h2 = edge.height * edge.height;
        return h2 / (std::sqrt(d * d + h2) + d);
    };
    edge.delta = extra(edge.d1) + extra(edge.d2);
    return edge;
}

void niceops::fresnel_integrals(double x, double& c, double& s)
{
    // Numerical Recipes (Press et al., 2nd ed., section 6.9): the power series
    // to |x| = 1.5, and beyond it the continued fraction for the
    // complementary error function evaluated by Lentz's method, both to
    // double precision.
    constexpr double eps = 1e-16;
    constexpr double tiny = 1e-300;
    constexpr int max_terms = 200;
    const double ax = std::fabs(x);
    if (ax < 1e-150)
    {
        c = ax;
        s = 0.0;
    }
    else if (ax <= 1.5)
    {
        // C = sum (-1)^k (pi/2)^(2k) x^(4k+1) / ((2k)! (4k+1)), S the odd terms.
        double sum = 0.0, sum_s = 0.0, sum_c = ax, sign = 1.0, term = ax;
        const double fact = 0.5 * pi * ax * ax;
        bool odd = true;
        int n = 3;
        for (int k = 1; k <= max_terms; k++)
        {
            term *= fact / k;
            sum += sign * term / n;
            const double test = std::fabs(sum) * eps;
            if (odd)
            {
                sign = -sign;
                sum_s = sum;
                sum = sum_c;
            }
            else
            {
                sum_c = sum;
                sum = sum_s;
            }
            if (term < test)
                break;
            odd = !odd;
            n += 2;
        }
        c = sum_c;
        s = sum_s;
    }
    else
    {
        const double pix2 = pi * ax * ax;
        std::complex<double> b(1.0, -pix2);
        std::complex<double> cc(1.0 / tiny, 0.0);
        std::complex<double> d = 1.0 / b;
        std::complex<double> h = d;
        int n = -1;
        for (int k = 2; k <= max_terms; k++)
        {
            n += 2;
            const double a = -static_cast<double>(n) * (n + 1);
            b += 4.0;
            d = 1.0 / (a * d + b);
            cc = b + a / cc;
            const std::complex<double> del = cc * d;
            h *= del;
            if (std::fabs(del.real() - 1.0) + std::fabs(del.imag()) < eps)
                break;
        }
        h *= std::complex<double>(ax, -ax);
        const std::complex<double> cs = std::complex<double>(0.5, 0.5)
                * (1.0 - std::complex<double>(std::cos(0.5 * pix2), std::sin(0.5 * pix2)) * h);
        c = cs.real();
        s = cs.imag();
    }
    if (x < 0.0)
    {
        c = -c;
        s = -s;
    }
}

double niceops::knife_edge_power(double v)
{
    double c = 0.0, s = 0.0;
    fresnel_integrals(v, c, s);
    return 0.5 * ((0.5 - c) * (0.5 - c) + (0.5 - s) * (0.5 - s));
}

void niceops::knife_edge_loss(const KnifeEdge& edge, const std::vector<double>& band_centers, double sound_speed,
        std::size_t first, std::size_t last, double* loss)
{
    assert(edge.delta > 0.0 && sound_speed > 0.0);
    // Eight-point Gauss-Legendre on [-1, 1]: the positive nodes and their weights.
    static constexpr double node[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267,
            0.9602898564975363};
    static constexpr double weight[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745,
            0.1012285362903763};
    const double third = std::pow(2.0, 1.0 / 6.0);
    last = std::min(last, band_centers.size());
    for (std::size_t band = first; band < last; band++)
    {
        const double lo = band_centers[band] / third;
        const double hi = band_centers[band] * third;
        const double mid = 0.5 * (hi + lo);
        const double half = 0.5 * (hi - lo);
        double mean = 0.0;
        for (int i = 0; i < 4; i++)
            for (const double side : {-1.0, 1.0})
            {
                const double f = mid + side * half * node[i];
                mean += 0.5 * weight[i] * knife_edge_power(2.0 * std::sqrt(edge.delta * f / sound_speed));
            }
        loss[band] = -10.0 * std::log10(mean);
    }
}
