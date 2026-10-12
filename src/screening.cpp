#include "screening.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

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
