#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

#include "building_paths.hpp"
#include "buildings.hpp"
#include "double_wedge.hpp"
#include "edge_source.hpp"
#include "ground_effect.hpp"
#include "point_paths.hpp"
#include "urban.hpp"
#include "util/constants.hpp"
#include "wedge.hpp"

using namespace niceops; // NOLINT

namespace {

Footprint box(double x0, double y0, double x1, double y1, double top)
{
    Footprint f;
    f.x = {x0, x1, x1, x0};
    f.y = {y0, y0, y1, y1};
    f.top = top;
    return f;
}

const std::vector<double> bands {125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0};

std::vector<double> levels(const Buildings& city, const Vec3& source, const Vec3& microphone, double rho, double floor,
        const PointPathOptions& options = {})
{
    const PointSourceField field(city, source, microphone, options);
    const std::vector<double> reflection(city.facades().size() * bands.size(), rho);
    std::vector<double> out(bands.size());
    field.band_levels(reflection, floor, bands, 343.2, out.data());
    return out;
}

// The tilted-step split is opt-in; its tests take it on.
PointPathOptions split_options()
{
    PointPathOptions options;
    options.tilted_step_split = true;
    return options;
}

// The shared face's aperture as it was first written, by the midpoint rule at
// `division` points a wavelength: the reference for its Filon rule.
std::complex<double> midpoint_aperture(const SharedFace& face, double sx, double sz, double rx, double rz, double k,
        double division)
{
    auto wrap = [](double a) { return a < 0.0 ? a + 2.0 * pi : a; };
    const double w = face.width, lambda = 2.0 * pi / k, s1 = std::hypot(sx, sz);
    const double extent = std::max(3.0 * std::sqrt(lambda * std::max(s1, w)), 30.0 * lambda);
    const double step = lambda / division;
    const double theta = wrap(std::atan2(rz, -(rx - w))), s3 = std::hypot(rx - w, rz);
    const double fx = -std::cos(face.n2 * pi), fz = std::sin(face.n2 * pi);
    const std::complex<double> launch = std::sqrt(k / (2.0 * pi)) * std::polar(1.0, -0.25 * pi);
    auto wave = [&](double x, double z) {
        const double d = std::hypot(rx - x, rz - z);
        return std::polar(1.0 / std::sqrt(d), k * d);
    };
    const bool above = rz > 0.0;
    std::complex<double> sum(0.0, 0.0), unbounded(0.0, 0.0);
    for (double t = 0.5 * step; t < extent; t += step)
    {
        const double taper = t < 0.5 * extent ? 1.0 : 0.5 * (1.0 + std::cos(pi * (t - 0.5 * extent) / (0.5 * extent)));
        const double in = std::hypot(sx, t - sz), r = std::hypot(w, t);
        const std::complex<double> amplitude
                = taper * 0.5 * (-sx / in + w / r) * launch * std::polar(1.0 / std::sqrt(in), k * in);
        const double theta0 = wrap(std::atan2(t, w));
        std::complex<double> field = std::polar(1.0 / std::sqrt(r * s3), k * (r + s3))
                * hard_wedge_coefficient(face.n2, theta, theta0, k, r * s3 / (r + s3));
        if (std::abs(theta - theta0) < pi)
            field += wave(0.0, t);
        if (theta < pi - theta0)
            field += wave(0.0, -t);
        if (theta0 > (face.n2 - 1.0) * pi && theta > (2.0 * face.n2 - 1.0) * pi - theta0)
        {
            const double along = -w * fx + t * fz;
            field += wave(w + 2.0 * along * fx + w, 2.0 * along * fz - t);
        }
        sum += amplitude * field;
        if (above)
            unbounded += amplitude * (wave(0.0, t) + wave(0.0, -t));
    }
    sum = step * (sum - unbounded);
    const double sin_g = -sz / s1, cos_g = -sx / s1;
    if (sin_g > 1e-6)
    {
        const std::complex<double> kirchhoff
                = 0.5 * (1.0 + cos_g) * std::polar(1.0, 0.25 * pi) / (std::sqrt(2.0 * pi * k) * sin_g);
        sum *= 0.5 * hard_wedge_coefficient(face.n1, 0.0, wrap(std::atan2(sz, sx)), k, 1e9) / kirchhoff;
    }
    return sum;
}

std::complex<double> midpoint_shared_face(const SharedFace& face, double sx, double sz, double rx, double rz, double k)
{
    const double forward = sz < 0.0 && sx < 0.0 ? std::atan2(-sz, -sx) : 0.0;
    const double reverse = rz < 0.0 && rx > face.width ? std::atan2(-rz, rx - face.width) : 0.0;
    std::complex<double> sum(0.0, 0.0);
    if (forward > 0.0)
        sum += forward * midpoint_aperture(face, sx, sz, rx, rz, k, 192.0);
    if (reverse > 0.0)
        sum += reverse
                * midpoint_aperture({face.width, face.n2, face.n1}, face.width - rx, rz, face.width - sx, sz, k, 192.0);
    return sum / (forward + reverse);
}

} // namespace

TEST_CASE("over open ground a point source's field is the direct ray and the floor's", "[point_paths]")
{
    const Buildings none;
    const Vec3 s(0.0, 0.0, 2.0), m(30.0, 4.0, 1.2);
    const PointSourceField field(none, s, m);
    REQUIRE(field.rays().size() == 2);
    const double r1 = (m - s).norm(), r2 = (Vec3(m.x, m.y, -m.z) - s).norm();
    const auto out = levels(none, s, m, 1.0, 0.8);
    const double half = std::pow(2.0, 1.0 / 6.0);
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        const double f1 = bands[b] / half, f2 = bands[b] * half;
        double energy = 0.0;
        for (int part = 0; part < 3; part++)
        {
            const double g1 = f1 * std::pow(f2 / f1, part / 3.0), g2 = f1 * std::pow(f2 / f1, (part + 1) / 3.0);
            energy += band_averaged_energy({1.0 / r1, std::sqrt(0.8) / r2}, {r1, r2}, 2.0 * pi * g1 / 343.2,
                              2.0 * pi * g2 / 343.2)
                    * (g2 - g1);
        }
        REQUIRE(out[b] == Catch::Approx(10.0 * std::log10(energy / (f2 - f1))).margin(1e-9));
    }
}

// Over impedance ground the floor's ray takes the spherical-wave reflection
// coefficient at its grazing angle, its whole length the image path, at each
// sub-band's middle: the two-ray sum to round-off.  With Q held at the band's
// centre that sum is the open-ground model (excess_ground_attenuation over the
// direct ray's spreading), to its band constants' six figures; Q's change
// across a band, which the model leaves out, moves it by up to 0.12 dB here,
// at a grazing source over hard ground (0.103 measured).
TEST_CASE("over impedance ground the floor's ray takes the ground's reflection coefficient", "[point_paths]")
{
    const Buildings none;
    const double c = 343.2;
    for (double flow : {200.0, 20000.0})
        for (const auto& [s, m] : {std::pair {Vec3(0.0, 0.0, 20.0), Vec3(150.0, 30.0, 1.2)},
                     std::pair {Vec3(0.0, 0.0, 3.0), Vec3(40.0, -10.0, 1.2)},
                     std::pair {Vec3(0.0, 0.0, 120.0), Vec3(60.0, 20.0, 1.2)}})
        {
            GroundModel ground;
            ground.flow_resistance = flow;
            ground.receiver_height = m.z;
            FloorModel floor;
            floor.admittance = ground_impedance(ground, sub_band_middles(bands, c), c);
            const PointSourceField field(none, s, m);
            REQUIRE(field.rays().size() == 2);
            std::vector<double> out(bands.size());
            field.band_levels({}, floor, bands, c, out.data());
            const double r1 = (m - s).norm(), r2 = (Vec3(m.x, m.y, -m.z) - s).norm();
            const double sin_grazing = (s.z + m.z) / r2;
            std::vector<double> centre(bands.size());
            const double d = std::hypot(m.x - s.x, m.y - s.y);
            excess_ground_attenuation(ground, ground_impedance(ground, bands, c), s.z, d, c, bands, centre.data());
            const std::vector<std::complex<double>> at_centre = ground_impedance(ground, bands, c);
            double worst = 0.0;
            for (std::size_t b = 0; b < bands.size(); b++)
            {
                const std::array<SubBand, 3> parts = sub_bands(bands[b], c);
                const double open = centre[b] - 20.0 * std::log10(r1);
                const auto q0 = spherical_reflection_coefficient(sin_grazing, pi * bands[b] * r2 / c, at_centre[b]);
                double energy[3], held[3];
                for (int part = 0; part < 3; part++)
                {
                    const SubBand& p = parts[part];
                    const auto q = spherical_reflection_coefficient(sin_grazing, 0.25 * (p.k1 + p.k2) * r2,
                            floor.admittance[3 * b + part]);
                    energy[part] = band_averaged_energy({1.0 / r1, std::complex<double>(q[0], q[1]) / r2}, {r1, r2},
                            p.k1, p.k2);
                    held[part] = band_averaged_energy({1.0 / r1, std::complex<double>(q0[0], q0[1]) / r2}, {r1, r2},
                            p.k1, p.k2);
                }
                REQUIRE(out[b] == Catch::Approx(band_level(bands[b], parts, energy)).margin(1e-9));
                REQUIRE(band_level(bands[b], parts, held) == Catch::Approx(open).margin(1e-5));
                worst = std::max(worst, std::abs(out[b] - open));
            }
            INFO("flow resistance " << flow << ", source " << s.z << " m up, " << d << " m off: worst "
                                    << worst << " dB from the open-ground model");
            REQUIRE(worst < 0.12);
        }
}

// Each floor bounce is its leg's: the leg between the points the ray
// diffracts at (or its ends), unfolded in the floor, its length and its
// grazing angle, ending at the microphone or not; and the ray's length is its
// legs', unfolded where they bounce.  Over the containers of EDGE's district,
// where rays bounce before, between and after their edges.
TEST_CASE("each floor bounce is its leg's, unfolded", "[point_paths]")
{
    const Buildings containers({box(4.09, 17.66, 6.53, 23.72, 7.77), box(14.81, 17.66, 17.25, 23.72, 7.77),
            box(1.6, 10.82, 4.04, 16.88, 7.77), box(1.6, 8.38, 13.72, 10.82, 5.18),
            box(20.21, 3.86, 22.65, 16.0, 7.77), box(13.5, -1.84, 19.56, 0.6, 7.77)});
    std::size_t checked = 0, kinds[4] = {0, 0, 0, 0};
    for (const Vec3& microphone : {Vec3(10.5, -1.1, 1.2), Vec3(-3.6, 6.4, 1.2), Vec3(18.0, 2.0, 1.2)})
    {
        const PointSourceField field(containers, Vec3(10.5, 20.5, 2.0), microphone);
        for (const PointRay& ray : field.rays())
        {
            if (ray.floor == 0 || !ray.facades.empty() || !ray.roof_facades.empty() || !ray.made_edges.empty())
                continue;
            // The points it turns at: its ends and edges, not the bounce
            // between two edges (on the floor).
            std::vector<Vec3> nodes;
            for (const Vec3& p : ray.path)
                if (p.z > 1e-9)
                    nodes.push_back(p);
            REQUIRE(nodes.size() == static_cast<std::size_t>(ray.wedges) + 2);
            std::vector<bool> used(nodes.size() - 1, false);
            double length = 0.0;
            for (int i = 0; i < ray.floor; i++)
            {
                const FloorBounce& b = ray.floor_bounce[static_cast<std::size_t>(i)];
                bool found = false;
                for (std::size_t k = 0; k + 1 < nodes.size() && !found; k++)
                {
                    const Vec3 &p = nodes[k], &q = nodes[k + 1];
                    const double leg = std::hypot(std::hypot(q.x - p.x, q.y - p.y), p.z + q.z);
                    if (used[k] || std::abs(leg - b.leg) > 1e-9 * leg)
                        continue;
                    REQUIRE(b.sin_grazing == Catch::Approx((p.z + q.z) / leg).epsilon(1e-12));
                    REQUIRE(b.microphone == (k + 2 == nodes.size()));
                    used[k] = found = true;
                    length += leg;
                }
                REQUIRE(found);
            }
            for (std::size_t k = 0; k + 1 < nodes.size(); k++)
                if (!used[k])
                    length += (nodes[k + 1] - nodes[k]).norm();
            REQUIRE(length == Catch::Approx(ray.length).epsilon(1e-9));
            checked++;
            kinds[std::min(ray.wedges, 3)]++;
        }
    }
    INFO(kinds[0] << " specular, " << kinds[1] << " once, " << kinds[2] << " twice, " << kinds[3] << " thrice");
    REQUIRE(kinds[0] > 0);
    REQUIRE(kinds[1] > 0);
    REQUIRE(kinds[2] > 0);
}

// A flush microphone on a plate: the direct ray takes P_d and the floor's,
// of the same length, Q P_r, so the band levels are |P_d + Q P_r|^2 over the
// spreading, from P_d and P_r at the elevation and azimuth the ray arrives
// from.
TEST_CASE("a flush microphone's plate takes P_d on the direct ray and Q P_r on the floor's", "[point_paths]")
{
    const Buildings none;
    const double c = 343.2;
    const Vec3 s(-30.0, 12.0, 25.0), m(50.0, -20.0, 0.0);
    GroundModel ground;
    FloorModel floor;
    floor.admittance = ground_impedance(ground, sub_band_middles(bands, c), c);
    const std::size_t problems = 3 * bands.size();
    std::vector<double> elevations, azimuths;
    auto direct_of = [](std::size_t j) { return std::complex<double>(0.9 + 0.01 * j, -0.05 * j); };
    auto image_of = [](std::size_t j) { return std::complex<double>(1.1 - 0.02 * j, 0.03 * j); };
    floor.plate = [&](double elevation, double azimuth, std::complex<double>* direct, std::complex<double>* image) {
        elevations.push_back(elevation);
        azimuths.push_back(azimuth);
        for (std::size_t j = 0; j < problems; j++)
        {
            direct[j] = direct_of(j);
            image[j] = image_of(j);
        }
    };
    const PointSourceField field(none, s, m);
    REQUIRE(field.rays().size() == 2);
    std::vector<double> out(bands.size());
    field.band_levels({}, floor, bands, c, out.data());
    const double r = (m - s).norm(), d = std::hypot(m.x - s.x, m.y - s.y);
    REQUIRE(elevations.size() == 2);
    for (std::size_t i = 0; i < 2; i++)
    {
        REQUIRE(elevations[i] == Catch::Approx(std::atan2(s.z, d)).margin(1e-12));
        REQUIRE(azimuths[i] == Catch::Approx(std::atan2(m.y - s.y, m.x - s.x)).margin(1e-12));
    }
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        const std::array<SubBand, 3> parts = sub_bands(bands[b], c);
        double energy[3];
        for (int part = 0; part < 3; part++)
        {
            const std::size_t j = 3 * b + static_cast<std::size_t>(part);
            const auto q = spherical_reflection_coefficient(s.z / r, 0.25 * (parts[part].k1 + parts[part].k2) * r,
                    floor.admittance[j]);
            energy[part] = std::norm(direct_of(j) + std::complex<double>(q[0], q[1]) * image_of(j)) / (r * r);
        }
        REQUIRE(out[b] == Catch::Approx(band_level(bands[b], parts, energy)).margin(1e-9));
    }
}

// The air absorbs each ray along its length, in each sub-band at that
// sub-band's middle frequency: with the floor absorbing everything, the field
// over open ground is the direct ray alone, whose band level then falls by
// exactly the sub-bands' mean of exp(-2 a L).
TEST_CASE("the air absorbs each ray along its length", "[point_paths]")
{
    const Buildings none;
    const Vec3 s(0.0, 0.0, 20.0), m(400.0, 30.0, 1.2);
    const PointSourceField field(none, s, m);
    const double r = (m - s).norm();
    Atmosphere air;
    air.relative_humidity = 40.0;
    air.temperature = 283.15;
    const std::vector<double> a = sub_band_absorption(air, bands, 343.2);
    std::vector<double> plain(bands.size()), absorbed(bands.size());
    field.band_levels({}, 0.0, bands, 343.2, plain.data());
    field.band_levels({}, 0.0, bands, 343.2, absorbed.data(), a);
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        const auto parts = sub_bands(bands[b], 343.2);
        double kept = 0.0, all = 0.0;
        for (int part = 0; part < 3; part++)
        {
            const double width = parts[part].g2 - parts[part].g1;
            kept += std::exp(-2.0 * a[3 * b + part] * r) * width;
            all += width;
            // The coefficient is ISO 9613-1's at the sub-band's middle, as
            // nepers of amplitude.
            const double middle = 0.5 * (parts[part].g1 + parts[part].g2);
            REQUIRE(a[3 * b + part] == Catch::Approx(air.attenuation_coefficient(middle) * std::log(10.0) / 20.0));
        }
        REQUIRE(absorbed[b] - plain[b] == Catch::Approx(10.0 * std::log10(kept / all)).margin(1e-9));
    }
    REQUIRE(plain.back() - absorbed.back() > 5.0);
}

TEST_CASE("the point-source field is reciprocal", "[point_paths]")
{
    // A small district: rays reflect, diffract at roof edges and corners, at
    // edges lit only by a reflection, and at two edges in turn.
    const Buildings city({box(0.0, 0.0, 6.0, 2.4, 7.8), box(8.0, 4.0, 10.4, 16.0, 7.8), box(-2.0, 8.0, 4.0, 10.4, 5.2),
            box(0.0, 14.0, 2.4, 20.0, 7.8)});
    const std::vector<std::pair<Vec3, Vec3>> pairs {{Vec3(3.0, 12.0, 2.0), Vec3(4.0, -3.0, 1.2)},
            {Vec3(6.0, 6.0, 2.0), Vec3(13.0, 10.0, 1.2)}, {Vec3(-4.0, 4.0, 1.5), Vec3(6.0, 18.0, 1.2)},
            {Vec3(-5.0, 4.0, 2.0), Vec3(6.0, 6.0, 1.2)}};
    std::size_t wedges[4] {0, 0, 0, 0}, shared = 0;
    for (const auto& [s, m] : pairs)
    {
        const auto forward = levels(city, s, m, 0.9, 1.0), back = levels(city, m, s, 0.9, 1.0);
        for (std::size_t b = 0; b < bands.size(); b++)
            REQUIRE(forward[b] == Catch::Approx(back[b]).margin(1e-9));
        const PointSourceField field(city, s, m);
        for (const auto& ray : field.rays())
        {
            wedges[std::min(ray.wedges, 3)]++;
            shared += ray.shared;
            // The path drawn is the ray: it runs from the source to the
            // microphone through its points (a bounce on the floor between two
            // edges among them), and off the floor's mirror it is as long as
            // the ray.
            REQUIRE(ray.path.size() >= 2);
            REQUIRE((ray.path.front() - s).norm() == 0.0);
            REQUIRE((ray.path.back() - m).norm() == 0.0);
            REQUIRE(ray.path.size()
                    == 2 + ray.facades.size() + static_cast<std::size_t>(ray.wedges) + (ray.between_floor ? 1 : 0));
            if (ray.floor == 0)
            {
                double length = 0.0;
                for (std::size_t v = 1; v < ray.path.size(); v++)
                    length += (ray.path[v] - ray.path[v - 1]).norm();
                REQUIRE(length == Catch::Approx(ray.length).epsilon(1e-9));
            }
        }
    }
    // Every kind was exercised, the uniform form over a shared face too.
    REQUIRE(wedges[0] > 0);
    REQUIRE(wedges[1] > 0);
    REQUIRE(wedges[2] > 0);
    REQUIRE(shared > 0);
}

// The ray set is the same with the ends swapped at every edge order: the
// twice-diffracted search takes both ends' images to one reflection even with
// none on a single diffraction's sides (--max_edge_reflections 0), where the
// source's tree once stopped at its own point and the two ways differed.
TEST_CASE("the ray set is reciprocal at every edge order", "[point_paths]")
{
    const Buildings city({box(0.0, 0.0, 6.0, 2.4, 7.8), box(8.0, 4.0, 10.4, 16.0, 7.8), box(-2.0, 8.0, 4.0, 10.4, 5.2),
            box(0.0, 14.0, 2.4, 20.0, 7.8)});
    const std::vector<std::pair<Vec3, Vec3>> pairs {{Vec3(3.0, 12.0, 2.0), Vec3(4.0, -3.0, 1.2)},
            {Vec3(6.0, 6.0, 2.0), Vec3(13.0, 10.0, 1.2)}, {Vec3(-4.0, 4.0, 1.5), Vec3(6.0, 18.0, 1.2)}};
    for (std::size_t order : {std::size_t {0}, std::size_t {1}, std::size_t {2}})
    {
        PointPathOptions options;
        options.max_edge_reflections = order;
        std::size_t twice = 0;
        for (const auto& [s, m] : pairs)
        {
            std::vector<double> lengths[2];
            for (int way = 0; way < 2; way++)
            {
                const PointSourceField field(city, way ? m : s, way ? s : m, options);
                for (const PointRay& ray : field.rays())
                    lengths[way].push_back(ray.length);
                std::sort(lengths[way].begin(), lengths[way].end());
                for (const PointRay& ray : field.rays())
                    twice += way == 0 && ray.wedges == 2;
            }
            INFO("max_edge_reflections " << order);
            REQUIRE(lengths[0].size() == lengths[1].size());
            for (std::size_t r = 0; r < lengths[0].size(); r++)
                REQUIRE(lengths[0][r] == Catch::Approx(lengths[1][r]).epsilon(1e-9));
        }
        REQUIRE(twice > 0);
    }
}

TEST_CASE("a point ray leaves the source along minus its length's gradient", "[point_paths]")
{
    // The reciprocity test's district, sources high and low, so that rays
    // leave toward facades, the floor and edges, and reflect in the floor
    // before an edge.  A ray's length is stationary in every point between
    // its ends, so moving the source by d moves it by -d . launch; rays are
    // matched across the three nudges by their facades, floors and edges.
    const Buildings city({box(0.0, 0.0, 6.0, 2.4, 7.8), box(8.0, 4.0, 10.4, 16.0, 7.8), box(-2.0, 8.0, 4.0, 10.4, 5.2),
            box(0.0, 14.0, 2.4, 20.0, 7.8)});
    const std::vector<std::pair<Vec3, Vec3>> pairs {{Vec3(3.0, 12.0, 2.0), Vec3(4.0, -3.0, 1.2)},
            {Vec3(6.0, 6.0, 2.0), Vec3(13.0, 10.0, 1.2)}, {Vec3(-4.0, 4.0, 1.5), Vec3(6.0, 18.0, 1.2)},
            {Vec3(-5.0, 4.0, 12.0), Vec3(6.0, 6.0, 1.2)}, {Vec3(5.0, 22.0, 30.0), Vec3(7.0, 6.0, 1.5)}};
    auto key = [](const PointRay& ray) {
        std::vector<std::size_t> k(ray.facades);
        k.push_back(1000000 + static_cast<std::size_t>(ray.floor));
        k.push_back(2000000 + static_cast<std::size_t>(ray.wedges));
        if (ray.wedges > 0)
        {
            k.push_back(3000000 + ray.edge_source.edges[0]);
            k.push_back(ray.edge_source.source.z < 0.0 ? 4000001 : 4000000);
        }
        if (ray.wedges > 1)
            k.push_back(5000000 + ray.edge_source.edges[1]);
        return k;
    };
    auto unique_rays = [&](const PointSourceField& field) {
        std::map<std::vector<std::size_t>, const PointRay*> found;
        std::map<std::vector<std::size_t>, int> count;
        for (const PointRay& ray : field.rays())
        {
            found[key(ray)] = &ray;
            count[key(ray)]++;
        }
        for (const auto& [k, n] : count)
            if (n > 1)
                found.erase(k);
        return found;
    };
    const double h = 1e-4;
    std::size_t checked[3] {0, 0, 0}, floor_first = 0;
    for (const auto& [s, m] : pairs)
    {
        const PointSourceField field(city, s, m);
        std::vector<PointSourceField> nudged_up, nudged_down;
        for (int axis = 0; axis < 3; axis++)
        {
            Vec3 d(axis == 0 ? h : 0.0, axis == 1 ? h : 0.0, axis == 2 ? h : 0.0);
            nudged_up.emplace_back(city, s + d, m);
            nudged_down.emplace_back(city, s - d, m);
        }
        const auto here = unique_rays(field);
        std::vector<std::map<std::vector<std::size_t>, const PointRay*>> up, down;
        for (int axis = 0; axis < 3; axis++)
        {
            up.push_back(unique_rays(nudged_up[static_cast<std::size_t>(axis)]));
            down.push_back(unique_rays(nudged_down[static_cast<std::size_t>(axis)]));
        }
        for (const auto& [k, ray] : here)
        {
            double gradient[3];
            bool matched = true;
            for (int axis = 0; axis < 3 && matched; axis++)
            {
                const auto a = up[static_cast<std::size_t>(axis)].find(k), b = down[static_cast<std::size_t>(axis)].find(k);
                matched = a != up[static_cast<std::size_t>(axis)].end() && b != down[static_cast<std::size_t>(axis)].end();
                if (matched)
                    gradient[axis] = (a->second->length - b->second->length) / (2.0 * h);
            }
            if (!matched)
                continue;
            const Vec3 launch = point_ray_launch(*ray);
            INFO("ray with " << ray->facades.size() << " reflections, floor " << ray->floor << ", " << ray->wedges
                             << " wedges: launch " << launch.x << " " << launch.y << " " << launch.z << ", -gradient "
                             << -gradient[0] << " " << -gradient[1] << " " << -gradient[2]);
            CHECK(std::abs(launch.norm() - 1.0) < 1e-12);
            CHECK(std::abs(launch.x + gradient[0]) < 1e-6);
            CHECK(std::abs(launch.y + gradient[1]) < 1e-6);
            CHECK(std::abs(launch.z + gradient[2]) < 1e-6);
            checked[std::min(ray->wedges, 2)]++;
            floor_first += ray->wedges > 0 && ray->edge_source.source.z < 0.0;
        }
    }
    INFO("checked " << checked[0] << " specular, " << checked[1] << " once and " << checked[2]
                    << " twice diffracted, " << floor_first << " off the floor before an edge");
    CHECK(checked[0] > 5);
    CHECK(checked[1] > 10);
    CHECK(checked[2] > 2);
    CHECK(floor_first > 5);
}

TEST_CASE("a ray's weight multiplies its amplitude in its band", "[point_paths]")
{
    // Unit weights are no weights, to the bit; a weight of 2 on every ray in
    // one band adds 20 lg 2 to that band alone; and a weight on one ray alone
    // is that ray's amplitude scaled, which a field with only it shows.
    const Buildings city({box(0.0, 0.0, 6.0, 2.4, 7.8), box(8.0, 4.0, 10.4, 16.0, 7.8)});
    const PointSourceField field(city, Vec3(-4.0, 1.0, 3.0), Vec3(12.0, 2.0, 1.5));
    REQUIRE(field.rays().size() > 3);
    const std::vector<double> reflection(city.facades().size() * bands.size(), 0.8);
    std::vector<double> plain(bands.size()), unit(bands.size()), doubled(bands.size());
    field.band_levels(reflection, 1.0, bands, 343.2, plain.data());
    std::vector<double> weights(field.rays().size() * bands.size(), 1.0);
    field.band_levels(reflection, 1.0, bands, 343.2, unit.data(), {}, weights.data());
    for (std::size_t r = 0; r < field.rays().size(); r++)
        weights[r * bands.size() + 2] = 2.0;
    field.band_levels(reflection, 1.0, bands, 343.2, doubled.data(), {}, weights.data());
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        CHECK(unit[b] == plain[b]);
        CHECK(doubled[b] - plain[b] == Catch::Approx(b == 2 ? 20.0 * std::log10(2.0) : 0.0).margin(1e-9));
    }
    // One ray weighted 0 is that ray left out: the energy of the rest, as a
    // field of the others by sub-band problems says.
    std::vector<double> one(field.rays().size() * bands.size(), 1.0);
    for (std::size_t b = 0; b < bands.size(); b++)
        one[0 * bands.size() + b] = 0.0;
    std::vector<double> without(bands.size());
    field.band_levels(reflection, 1.0, bands, 343.2, without.data(), {}, one.data());
    std::vector<double> energy(3 * bands.size());
    field.sub_band_problems(reflection, 1.0, bands, 343.2,
            [&](std::size_t band, int part, const std::vector<std::complex<double>>& amplitude,
                    const std::vector<double>& length, const SubBand& sub) {
                std::vector<std::complex<double>> kept;
                std::vector<double> kept_length;
                for (std::size_t i = 0; i < amplitude.size(); i++)
                    if (std::abs(amplitude[i]) > 0.0)
                    {
                        kept.push_back(amplitude[i]);
                        kept_length.push_back(length[i]);
                    }
                energy[3 * band + part] = band_averaged_energy(kept, kept_length, sub.k1, sub.k2);
            },
            {}, one.data());
    for (std::size_t b = 0; b < bands.size(); b++)
        CHECK(without[b] == Catch::Approx(band_level(bands[b], sub_bands(bands[b], 343.2), &energy[3 * b])).margin(1e-12));
}

// A roof under a reflected leg is found whatever else stands in the district:
// a low roof beside a tall wall's foot, the source's ray reflecting in the wall
// and then bouncing on the roof to a microphone above it (specular), and a
// corner's diffracted ray doing the same; with a hundred small buildings a
// kilometre off, which make the plan grid's cells small, the rays are the same.
// The unfolded line's plan box, which the search once took, misses the roof.
TEST_CASE("a roof under a reflected leg is found whatever stands far off", "[point_paths]")
{
    for (const bool diffracted : {false, true})
    {
        std::vector<Footprint> near {box(10.0, diffracted ? 5.0 : -5.0, 12.0, 35.0, 50.0), box(0.5, 25.0, 3.0, 29.0, 10.0)};
        if (diffracted)
            near.push_back(box(5.0, -10.0, 9.0, 0.0, 40.0));
        std::vector<Footprint> crowded = near;
        for (int i = 0; i < 100; i++)
            crowded.push_back(box(1000.0 + 3.0 * i, 1000.0, 1001.0 + 3.0 * i, 1001.0, 5.0));
        const Vec3 source = diffracted ? Vec3(-10.0, -20.0, 28.2) : Vec3(5.0, 0.0, 20.0);
        const Vec3 microphone(0.0, 30.0, 11.0);
        PointPathOptions options;
        options.diffraction = diffracted;
        options.double_diffraction = false;
        std::vector<double> lengths[2];
        for (int which = 0; which < 2; which++)
        {
            const Buildings city(which == 0 ? near : crowded);
            const PointSourceField field(city, source, microphone, options);
            for (const PointRay& ray : field.rays())
                if (!ray.roof_facades.empty() && ray.wedges == (diffracted ? 1 : 0))
                    lengths[which].push_back(ray.length);
            std::sort(lengths[which].begin(), lengths[which].end());
        }
        INFO((diffracted ? "diffracted" : "specular"));
        REQUIRE(!lengths[0].empty());
        REQUIRE(lengths[1] == lengths[0]);
    }
}

// Two roofs level across a side they share are one roof there: a bounce
// exactly on the join is one ray, not one per building, and the field 1 mm off
// the join is the same to round-off (two rays there put it 1.4 dB high).
TEST_CASE("a bounce on two level roofs' join is one ray", "[point_paths]")
{
    const Buildings city({box(0.0, 0.0, 10.0, 10.0, 6.0), box(10.0, 0.0, 20.0, 10.0, 6.0)});
    const std::vector<double> bands {500.0, 1000.0};
    const std::vector<double> reflection(city.facades().size() * bands.size(), 1.0);
    PointPathOptions options;
    options.diffraction = false;
    double level[2][2];
    for (int which = 0; which < 2; which++)
    {
        const Vec3 microphone(15.0 + which * 1e-3, 5.0, 20.0);
        const PointSourceField field(city, Vec3(5.0, 5.0, 20.0), microphone, options);
        std::size_t roofs = 0;
        for (const PointRay& ray : field.rays())
            roofs += ray.roof_facades.size();
        REQUIRE(roofs == 1);
        field.band_levels(reflection, 1.0, bands, 343.2, level[which]);
    }
    REQUIRE(level[0][0] == Catch::Approx(level[1][0]).margin(0.01));
    REQUIRE(level[0][1] == Catch::Approx(level[1][1]).margin(0.01));
}

// A courtyard closed by four buildings meeting at T-junctions, each side
// wall's ends against the long walls' middles: from a low source outside,
// every ray reaches the microphone over a roof.  Two buildings' walls along
// one line were two separate solids, and rays diffracted at the junctions'
// corners and ran along the seam between them into the courtyard.
TEST_CASE("a courtyard closed at T-junctions is reached over its roofs alone", "[point_paths]")
{
    const double top = 10.0;
    const Buildings well({box(-12.0, -12.0, 12.0, -9.0, top), box(-12.0, 9.0, 12.0, 12.0, top),
            box(-12.0, -9.0, -9.0, 9.0, top), box(9.0, -9.0, 12.0, 9.0, top)});
    REQUIRE(well.shared_sides() == 8);
    for (const Vec3& source : {Vec3(-200.0, 3.0, 2.0), Vec3(150.0, -120.0, 5.0), Vec3(-12.0, -40.0, 1.0)})
    {
        INFO("source " << source.x << ", " << source.y << ", " << source.z);
        const PointSourceField field(well, source, Vec3(0.0, 0.0, 1.2));
        REQUIRE(!field.rays().empty());
        for (const PointRay& ray : field.rays())
        {
            double highest = 0.0;
            for (const Vec3& p : ray.path)
                highest = std::max(highest, p.z);
            REQUIRE(highest >= top - 1e-9);
        }
    }
}

// --fade adds partners and nothing else: the rays it leaves without a made
// edge are the rays without it.  Two ways that failed: a reflected leg back
// through the other wing of the L-shaped building it leaves (its own
// building's blocking went untested), and a ray across one building with three
// beside the leg nearer than it, where keeping the three nearest dropped the
// one across and let the bare ray through it.
TEST_CASE("the fade adds partners and nothing else", "[point_paths]")
{
    Footprint ell;
    ell.x = {150.0, 190.0, 190.0, 170.0, 170.0, 150.0};
    ell.y = {0.0, 0.0, 20.0, 20.0, 50.0, 50.0};
    ell.top = 25.0;
    std::vector<Footprint> screened {box(100.0, -50.0, 102.0, 50.0, 50.0), box(50.0, 20.0, 52.0, 29.7, 10.0)};
    for (double x : {20.0, 40.0, 60.0})
        screened.push_back(box(x - 2.0, 0.2 * x + 0.001, x, 0.2 * x + 2.001, 10.0));
    struct Scene
    {
        std::vector<Footprint> footprints;
        Vec3 source, microphone;
    };
    const std::vector<Scene> scenes {{{ell}, Vec3(180.0, 26.0, 5.0), Vec3(200.0, 10.0, 5.0)},
            {screened, Vec3(0.0, 0.0, 2.0), Vec3(0.0, 40.0, 2.0)}};
    std::size_t made = 0;
    for (std::size_t s = 0; s < scenes.size(); s++)
    {
        const Buildings city(scenes[s].footprints);
        std::vector<double> lengths[2];
        for (int fade = 0; fade < 2; fade++)
        {
            PointPathOptions options;
            options.diffraction = false;
            options.fade_reach = fade ? 0.5 : 0.0;
            const PointSourceField field(city, scenes[s].source, scenes[s].microphone, options);
            for (const PointRay& ray : field.rays())
            {
                if (ray.made_edges.empty())
                    lengths[fade].push_back(ray.length);
                made += !ray.made_edges.empty();
            }
            std::sort(lengths[fade].begin(), lengths[fade].end());
        }
        INFO("scene " << s);
        REQUIRE(lengths[1].size() == lengths[0].size());
        for (std::size_t r = 0; r < lengths[0].size(); r++)
            REQUIRE(lengths[1][r] == Catch::Approx(lengths[0][r]).epsilon(1e-12));
    }
    // The second scene's partners exist: the fade did something.
    REQUIRE(made > 0);
}

// A ray's amplitudes across the wavenumbers at once are point_ray_amplitude's
// bits at each: a made edge's partner, an end on either side or both, past
// an end, a grazing wedge, and a reflected amplitude other than 1.
TEST_CASE("a ray's amplitudes across the wavenumbers are its amplitude at each", "[point_paths]")
{
    std::vector<Footprint> screened {box(100.0, -50.0, 102.0, 50.0, 50.0), box(50.0, 20.0, 52.0, 29.7, 10.0)};
    for (double x : {20.0, 40.0, 60.0})
        screened.push_back(box(x - 2.0, 0.2 * x + 0.001, x, 0.2 * x + 2.001, 10.0));
    const Buildings city(screened);
    PointPathOptions options;
    options.fade_reach = 0.5;
    const PointSourceField field(city, Vec3(0.0, 0.0, 2.0), Vec3(0.0, 40.0, 2.0), options);
    std::vector<double> k;
    for (int j = 0; j < 51; j++)
        k.push_back(2.0 * pi * 100.0 * std::pow(2.0, j / 9.0) / 343.2);
    std::vector<PointRay> rays;
    for (const PointRay& ray : field.rays())
    {
        if (ray.shared || ray.profile)
            continue;
        rays.push_back(ray);
        if (ray.wedges == 0)
            continue;
        // Its first wedge's ends made each other way.
        const EdgeEnds as_traced = ray.wedge[0].ends;
        for (int variant = 0; variant < 5; variant++)
        {
            PointRay v = ray;
            EdgeEnds& e = v.wedge[0].ends;
            e = as_traced;
            e.lo = variant == 1 ? std::numeric_limits<double>::infinity() : 0.3 + 0.1 * variant;
            e.hi = variant == 2 ? std::numeric_limits<double>::infinity() : 0.8;
            e.past = variant == 3 ? -1 : variant == 4 ? 1 : 0;
            e.scale = variant >= 3 ? 0.7 : 1.0;
            v.wedge[0].grazing = variant % 2 == 0;
            rays.push_back(v);
        }
    }
    std::size_t made = 0, compared = 0;
    std::vector<std::complex<double>> out(k.size()), scratch;
    for (const PointRay& ray : rays)
    {
        made += !ray.made_edges.empty();
        std::fill(out.begin(), out.end(), std::complex<double>(0.7));
        point_ray_amplitudes(ray, k.data(), k.size(), out.data(), scratch);
        for (std::size_t j = 0; j < k.size(); j++)
        {
            const std::complex<double> one = point_ray_amplitude(ray, 0.7, k[j]);
            REQUIRE(std::memcmp(&one, &out[j], sizeof one) == 0);
            compared++;
        }
    }
    REQUIRE(made > 0);
    REQUIRE(compared > 51 * 20);
}

TEST_CASE("a roof reflects as the floor does", "[point_paths]")
{
    // A 6 m block, the source and the microphone above it: the roof's ray is
    // the microphone's image in the roof's plane, and it is the only new one;
    // with roofs off it is gone.
    const Buildings city({box(0.0, 0.0, 10.0, 10.0, 6.0)});
    const Vec3 s(2.0, 5.0, 20.0), m(8.0, 5.0, 10.0);
    PointPathOptions without;
    without.roofs = false;
    const PointSourceField on(city, s, m), off(city, s, m, without);
    std::size_t roof_rays = 0;
    for (const PointRay& ray : on.rays())
        if (!ray.roof_facades.empty())
        {
            roof_rays++;
            CHECK(ray.wedges == 0);
            CHECK(ray.length == Catch::Approx((Vec3(m.x, m.y, 12.0 - m.z) - s).norm()).epsilon(1e-12));
            CHECK(ray.source_plane == 6.0);
        }
    CHECK(roof_rays == 1);
    CHECK(on.rays().size() == off.rays().size() + 1);
    for (const PointRay& ray : off.rays())
        CHECK(ray.roof_facades.empty());
    // Below the roof's height either end leaves it out.
    for (const Vec3& low : {Vec3(8.0, 5.0, 5.0), Vec3(14.0, 5.0, 3.0)})
    {
        const PointSourceField field(city, s, low);
        for (const PointRay& ray : field.rays())
            CHECK((ray.roof_facades.empty() || ray.wedges > 0));
    }
}

TEST_CASE("the field is continuous where a roof's edge hides its reflection", "[point_paths]")
{
    // A microphone rising past the roof's reflection boundary at its edge:
    // the edge's diffracted ray changes sign there by the roof's ray, so the
    // two together are continuous.  Source (3, 20, 15) above a 6 m block whose
    // east edge is at x = 10; the boundary at x = 14 is the line from the
    // source's image (3, 20, -3) through (10, 20, 6).  (Without the roof's ray
    // the band-averaged energy hardly steps there -- the sign flip is
    // incoherent with rays several wavelengths shorter -- but misses the
    // reflection's energy throughout the zone it lights.)
    const Buildings city({box(0.0, 0.0, 10.0, 40.0, 6.0)});
    const Vec3 s(3.0, 20.0, 15.0);
    const double boundary = -3.0 + 9.0 * 11.0 / 7.0;
    for (const double eps : {1e-6, 1e-4})
    {
        const auto below = levels(city, s, Vec3(14.0, 20.0, boundary - eps), 1.0, 1.0),
                   above = levels(city, s, Vec3(14.0, 20.0, boundary + eps), 1.0, 1.0);
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            INFO(bands[b] << " Hz, " << eps << " m either side: " << above[b] - below[b] << " dB");
            CHECK(std::abs(above[b] - below[b]) < 0.01);
        }
    }
}

TEST_CASE("on a wide roof the field is the direct ray and the roof's", "[point_paths]")
{
    // A microphone on a 200 m roof, its edges 90 m off: the direct ray and
    // its reflection in the roof, as two rays over a rigid floor; without
    // roofs, the reflection's 3 dB or so is missing.
    const Buildings city({box(0.0, 0.0, 200.0, 200.0, 10.0)});
    const Vec3 s(100.0, 100.0, 40.0), m(110.0, 100.0, 11.5);
    const PointSourceField field(city, s, m);
    PointPathOptions without;
    without.roofs = false;
    const auto on = levels(city, s, m, 1.0, 1.0), off = levels(city, s, m, 1.0, 1.0, without);
    const Vec3 image(m.x, m.y, 20.0 - m.z);
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        const auto parts = sub_bands(bands[b], 343.2);
        double energy[3];
        for (int part = 0; part < 3; part++)
        {
            const std::vector<std::complex<double>> amplitude {1.0 / (m - s).norm(), 1.0 / (image - s).norm()};
            const std::vector<double> length {(m - s).norm(), (image - s).norm()};
            energy[part] = band_averaged_energy(amplitude, length, parts[static_cast<std::size_t>(part)].k1,
                    parts[static_cast<std::size_t>(part)].k2);
        }
        const double exact = band_level(bands[b], parts, energy);
        INFO(bands[b] << " Hz: " << on[b] - exact << " dB off the two rays, " << on[b] - off[b] << " over no roofs");
        CHECK(std::abs(on[b] - exact) < 0.05);
        CHECK(on[b] - off[b] > 1.5);
    }
}

TEST_CASE("rays off a roof are reciprocal and leave along minus their length's gradient", "[point_paths]")
{
    // A low block, a street and a tall one: from above the low roof the
    // sound reflects up to the tall roof's near edge and diffracts down into
    // the street, with and without the floor and a facade.
    const Buildings city({box(0.0, 0.0, 10.0, 40.0, 6.0), box(20.0, 0.0, 30.0, 40.0, 20.0)});
    const Vec3 s(0.0, 20.0, 15.0), m(15.0, 20.0, 1.5);
    const auto forward = levels(city, s, m, 0.8, 0.9), back = levels(city, m, s, 0.8, 0.9);
    for (std::size_t b = 0; b < bands.size(); b++)
        CHECK(forward[b] == Catch::Approx(back[b]).margin(1e-9));
    const PointSourceField field(city, s, m);
    const double h = 1e-4;
    std::size_t checked = 0;
    for (const PointRay& ray : field.rays())
    {
        if (ray.roof_facades.empty())
            continue;
        // A twice-diffracted ray off a roof draws its bounce and has no
        // source_plane (NaN).
        auto same = [&](const PointRay& other) {
            return other.facades == ray.facades && other.floor == ray.floor && other.wedges == ray.wedges
                    && other.roof_facades == ray.roof_facades && other.between_floor == ray.between_floor
                    && (other.source_plane == ray.source_plane
                            || (std::isnan(other.source_plane) && std::isnan(ray.source_plane)))
                    && (ray.wedges == 0 || other.edge_source.edges[0] == ray.edge_source.edges[0])
                    && (ray.wedges < 2 || other.edge_source.edges[1] == ray.edge_source.edges[1])
                    && other.path.size() == ray.path.size();
        };
        double gradient[3];
        bool matched = true;
        for (int axis = 0; axis < 3 && matched; axis++)
        {
            const Vec3 d(axis == 0 ? h : 0.0, axis == 1 ? h : 0.0, axis == 2 ? h : 0.0);
            const PointSourceField up(city, s + d, m), down(city, s - d, m);
            const PointRay *u = nullptr, *w = nullptr;
            for (const PointRay& r : up.rays())
                if (same(r))
                    u = &r;
            for (const PointRay& r : down.rays())
                if (same(r))
                    w = &r;
            matched = u != nullptr && w != nullptr;
            if (matched)
                gradient[axis] = (u->length - w->length) / (2.0 * h);
        }
        REQUIRE(matched);
        const Vec3 launch = point_ray_launch(ray);
        CHECK(std::abs(launch.x + gradient[0]) < 1e-6);
        CHECK(std::abs(launch.y + gradient[1]) < 1e-6);
        CHECK(std::abs(launch.z + gradient[2]) < 1e-6);
        checked++;
    }
    CHECK(checked >= 3);
}

TEST_CASE("a far point source is the plane wave", "[point_paths]")
{
    // One building, the microphone near it; the source a thousand kilometres
    // off in the direction u.  Scaled by the range, the point source's field
    // is the plane-wave builder's, every ray of it.
    const Buildings city({box(0.0, 0.0, 30.0, 12.0, 20.0)});
    const Vec3 m(-6.0, 5.0, 1.5);
    const double range = 1e6;
    BuildingPathOptions plane;
    plane.max_diffracted_order = 2;
    const BuildingReceiver receiver(city, m, {}, plane);
    PointPathOptions point;
    point.double_diffraction = false;
    point.source_reflections = false;
    point.max_edge_reflections = 2;
    point.radius = 2e6;
    const std::vector<double> reflection(city.facades().size() * bands.size(), 0.85);
    for (const Vec3 u : {local_direction(250.0, 10.0, 90.0), local_direction(80.0, 30.0, 90.0),
                 local_direction(120.0, 5.0, 90.0)})
    {
        std::vector<double> ratio(bands.size());
        receiver.ratio(u, reflection, 0.9, bands, 343.2, ratio.data());
        const auto far = levels(city, m + range * u, m, 0.85, 0.9, point);
        // The plane wave's total is its ratio times open ground's band energy.
        const auto open = levels(Buildings {}, m + range * u, m, 0.85, 0.9, point);
        for (std::size_t b = 0; b < bands.size(); b++)
            REQUIRE(far[b] - open[b] == Catch::Approx(10.0 * std::log10(ratio[b])).margin(2e-3));
    }
}

// The same over the roofs (BuildingPathOptions::roofs): a microphone 3 m over
// a flat roof and 2 m over a slope, and one lying on each.  The plane wave's
// rays off the roof -- specular, and diffracted to the microphone's mirror in
// it -- are the far point source's, and a microphone on the roof takes the
// roof's reflection as one on a wall does.  The source is 1e8 m off: the
// floor's reflection reaches a microphone 23 m up from 130 m away, where
// 1e6 m leaves a millimetre of curvature in its path.
TEST_CASE("a far point source is the plane wave over a roof", "[point_paths]")
{
    Footprint shed = box(0.0, 0.0, 30.0, 12.0, 40.0);
    shed.roof = {{0.0, 0.4, 14.0}};
    const RoofPlane& slope = shed.roof.front();
    const Vec3 on_slope(15.0, 6.0, slope.at(15.0, 6.0));
    struct Case
    {
        Buildings city;
        Vec3 microphone;
        bool over; // over the roof, not on it
    };
    const std::vector<Case> cases {{Buildings({box(0.0, 0.0, 30.0, 12.0, 20.0)}), Vec3(15.0, 6.0, 23.0), true},
            {Buildings({box(0.0, 0.0, 30.0, 12.0, 20.0)}), Vec3(15.0, 6.0, 20.0), false},
            {Buildings({shed}), on_slope + 2.0 * slope.normal(), true}, {Buildings({shed}), on_slope, false}};
    const double range = 1e8;
    BuildingPathOptions plane;
    plane.max_diffracted_order = 2;
    PointPathOptions point;
    point.double_diffraction = false;
    point.source_reflections = false;
    point.max_edge_reflections = 2;
    point.radius = 2e6;
    for (const Case& c : cases)
    {
        const BuildingReceiver receiver(c.city, c.microphone, {}, plane);
        const std::vector<double> reflection(c.city.facades().size() * bands.size(), 0.85);
        std::size_t off_roof = 0;
        for (const Vec3 u : {local_direction(250.0, 10.0, 90.0), local_direction(80.0, 30.0, 90.0),
                     local_direction(120.0, 5.0, 90.0), local_direction(10.0, 60.0, 90.0)})
        {
            for (const BuildingRay& ray : receiver.rays(u))
                off_roof += ray.roof;
            std::vector<double> ratio(bands.size());
            receiver.ratio(u, reflection, 0.9, bands, 343.2, ratio.data());
            const auto far = levels(c.city, c.microphone + range * u, c.microphone, 0.85, 0.9, point);
            const auto open = levels(Buildings {}, c.microphone + range * u, c.microphone, 0.85, 0.9, point);
            for (std::size_t b = 0; b < bands.size(); b++)
            {
                INFO("microphone " << c.microphone.x << ", " << c.microphone.y << ", " << c.microphone.z << "; u "
                                   << u.x << ", " << u.y << ", " << u.z << "; band " << bands[b]);
                CHECK(far[b] - open[b] == Catch::Approx(10.0 * std::log10(ratio[b])).margin(2e-3));
            }
        }
        // Over a roof some rays are off it; on one, none.
        CHECK((off_roof > 0) == c.over);
    }
}

// A microphone lying on a roof, flat or a slope, and one 0.1 mm over it: the
// roof's reflection taken at the microphone (1 + R) is the reflected ray the
// other hears from 0.1 mm below, so over rigid surfaces their tables are one.
// (An absorbing roof's own edges diffract as hard wedges, whose face doubles
// the wave a hair above it where the roof gives 1 + R on it: 0.1 dB at
// R^2 = 0.85.  A wall corner's end wave from its top, on the roof's plane, is
// doubled on the roof and not over it, the point-source model's choice too:
// 0.09 dB at 125 Hz here, so the edges' ends are left out.  Within contact,
// 1e-6 m, a point is on the roof.)
TEST_CASE("a table on a roof is the table a hair above it", "[point_paths]")
{
    Footprint shed = box(0.0, 0.0, 30.0, 12.0, 40.0);
    shed.roof = {{0.0, 0.4, 14.0}};
    const RoofPlane& slope = shed.roof.front();
    const Vec3 on_slope(15.0, 6.0, slope.at(15.0, 6.0));
    const Buildings flat({box(0.0, 0.0, 30.0, 12.0, 20.0)}), pitched({shed});
    for (const auto& [city, on, normal] : {std::tuple {&flat, Vec3(15.0, 6.0, 20.0), Vec3(0.0, 0.0, 1.0)},
                 std::tuple {&pitched, on_slope, slope.normal()}})
    {
        BuildingPathOptions options;
        options.max_diffracted_order = 2;
        options.edge_ends = false;
        const BuildingReceiver lying(*city, on, {}, options), hair(*city, on + 1e-4 * normal, {}, options);
        const std::vector<double> reflection(city->facades().size() * bands.size(), 1.0);
        for (const Vec3 u : {local_direction(250.0, 10.0, 90.0), local_direction(80.0, 30.0, 90.0),
                     local_direction(10.0, 60.0, 90.0)})
        {
            std::vector<double> a(bands.size()), b(bands.size());
            lying.ratio(u, reflection, 0.9, bands, 343.2, a.data());
            hair.ratio(u, reflection, 0.9, bands, 343.2, b.data());
            for (std::size_t k = 0; k < bands.size(); k++)
            {
                INFO("u " << u.x << ", " << u.y << ", " << u.z << "; band " << bands[k]);
                CHECK(10.0 * std::log10(a[k] / b[k]) == Catch::Approx(0.0).margin(0.01));
            }
        }
    }
}

// Where every boundary's partner is traced the field is continuous: a ray
// that appears or vanishes there hands its share to an edge's diffracted
// ray.  West of a tall box, a source; north of it, lines of microphones 2.5 cm
// apart crossing the reflection boundary of the box's west wall and the
// shadow boundary of its north-west corner, directly and off the floor, near
// the corner and further off.  The box is tall enough that its roof's edges,
// and the corner's end, are far from every path, and the near line stops
// short of where the north-east corner's twice-diffracted ray reaches the
// floor: carrying it on would need the floor between the two edges, which the
// model leaves out (a 1.3 dB step at 125 Hz, 22.5 m along).  A receiver
// exactly on a boundary, as the lines' round numbers put some, takes the side
// its ray is kept on; it took the other, 1 + 1/2 for 1/2, 9 dB high.  A step
// between neighbours above a fraction of a dB in any band is a boundary
// handled wrong; the largest the field takes on its own is under 0.35 dB.
// (urban/continuity.py in the validation harness sweeps whole cities the same
// way, and finds where the model's ray set is not complete.)
TEST_CASE("the field is continuous across a corner's shadow and reflection boundaries", "[point_paths]")
{
    const Buildings tower({box(0.0, 0.0, 10.0, 10.0, 200.0)});
    const Vec3 source(-10.0, 5.0, 2.0);
    const std::vector<double> continuity_bands {125.0, 500.0, 2000.0};
    const std::vector<double> reflection(tower.facades().size() * continuity_bands.size(), 0.9);
    double worst = 0.0;
    std::size_t lit = 0, reflected = 0, shadowed = 0;
    for (const auto& [y, last] : {std::pair {12.0, 15.0}, std::pair {20.0, 40.0}})
    {
        std::vector<double> previous;
        for (double x = -30.0; x <= last + 1e-9; x += 0.025)
        {
            const PointSourceField field(tower, source, Vec3(x, y, 1.2));
            bool direct = false, specular = false;
            for (const PointRay& ray : field.rays())
            {
                direct |= ray.wedges == 0 && ray.facades.empty();
                specular |= ray.wedges == 0 && !ray.facades.empty();
            }
            lit += direct;
            reflected += specular;
            shadowed += !direct;
            std::vector<double> out(continuity_bands.size());
            field.band_levels(reflection, 0.9, continuity_bands, 343.2, out.data());
            if (!previous.empty())
                for (std::size_t b = 0; b < out.size(); b++)
                    worst = std::max(worst, std::abs(out[b] - previous[b]));
            previous = out;
        }
    }
    INFO("largest step between neighbours " << worst << " dB; " << lit << " lit, " << reflected << " reflected, "
                                            << shadowed << " in the corner's shadow");
    // Both boundaries are crossed on both lines.
    REQUIRE(reflected > 0);
    REQUIRE(shadowed > 0);
    REQUIRE(worst < 0.5);
}

TEST_CASE("the field is continuous where a thick building's roof plane hides its far edge", "[point_paths]")
{
    // A long building 10 m thick and 8 m high, the source low before it, the
    // microphone behind it crossing the roof plane: above, the near edge's ray
    // grazes the roof and reaches it; below, only the ray diffracted at both
    // edges does.  The uniform form over the roof (double_wedge.hpp) carries
    // the one into the other.
    const Buildings city({box(-200.0, 0.0, 200.0, 10.0, 8.0)});
    const Vec3 s(0.0, -15.0, 2.0);
    PointPathOptions options;
    options.max_reflections = 0;
    options.max_edge_reflections = 0;
    for (const double y : {14.0, 25.0})
    {
        const auto above = levels(city, s, Vec3(0.0, y, 8.0 + 1e-3), 1.0, 0.0, options);
        const auto under = levels(city, s, Vec3(0.0, y, 8.0 - 1e-3), 1.0, 0.0, options);
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            INFO("microphone " << y << " m behind, band " << bands[b]);
            REQUIRE(std::abs(above[b] - under[b]) < 0.5);
        }
    }
}

namespace {

// The long block's boundary-element solutions (tools/thick_barrier_bem.py):
// frequency, the source's y and height, the microphone's distance behind the
// block and height, and the level, dB re the free field.
struct BlockCase
{
    double f, source_y, source_z, behind, z, exact;
};

const std::vector<BlockCase>& long_block_cases()
{
    static const std::vector<BlockCase> cases {
        {125, -25, 4.5, 5, 1.2, -17.01},
        {125, -25, 4.5, 5, 4.8, -7.37},
        {125, -25, 4.5, 5, 5.26, -7.20},
        {125, -25, 4.5, 5, 6.5, -0.48},
        {125, -25, 4.5, 25, 1.2, -6.88},
        {125, -25, 4.5, 25, 4.8, -5.72},
        {125, -25, 4.5, 25, 5.26, -4.29},
        {125, -25, 4.5, 25, 6.5, -2.38},
        {125, -15, 6.5, 5, 1.2, -19.05},
        {125, -15, 6.5, 5, 4.8, -9.04},
        {125, -15, 6.5, 25, 1.2, -6.60},
        {125, -15, 6.5, 25, 4.8, -4.40},
        {125, -10, 2, 5, 1.2, -19.70},
        {125, -10, 2, 5, 4.8, -11.46},
        {125, -10, 2, 5, 5.26, -12.13},
        {125, -10, 2, 5, 6.5, -5.35},
        {125, -10, 2, 25, 1.2, -10.42},
        {125, -10, 2, 25, 4.8, -10.90},
        {125, -10, 2, 25, 5.26, -9.35},
        {125, -10, 2, 25, 6.5, -7.42},
        {250, -25, 4.5, 5, 1.2, -18.15},
        {250, -25, 4.5, 5, 4.8, -11.96},
        {250, -25, 4.5, 5, 5.26, -7.76},
        {250, -25, 4.5, 5, 6.5, -5.11},
        {250, -25, 4.5, 25, 1.2, -17.58},
        {250, -25, 4.5, 25, 4.8, -12.95},
        {250, -25, 4.5, 25, 5.26, -10.00},
        {250, -25, 4.5, 25, 6.5, -5.81},
        {250, -15, 6.5, 5, 1.2, -10.90},
        {250, -15, 6.5, 5, 4.8, -4.96},
        {250, -15, 6.5, 25, 1.2, -9.85},
        {250, -15, 6.5, 25, 4.8, -5.93},
        {250, -10, 2, 5, 1.2, -20.69},
        {250, -10, 2, 5, 4.8, -16.45},
        {250, -10, 2, 5, 5.26, -11.72},
        {250, -10, 2, 5, 6.5, -9.54},
        {250, -10, 2, 25, 1.2, -20.80},
        {250, -10, 2, 25, 4.8, -20.37},
        {250, -10, 2, 25, 5.26, -17.72},
        {250, -10, 2, 25, 6.5, -11.11},
        {500, -25, 4.5, 5, 1.2, -32.26},
        {500, -25, 4.5, 5, 4.8, -8.09},
        {500, -25, 4.5, 5, 5.26, -6.71},
        {500, -25, 4.5, 5, 6.5, 1.20},
        {500, -25, 4.5, 25, 1.2, -14.40},
        {500, -25, 4.5, 25, 4.8, -7.17},
        {500, -25, 4.5, 25, 5.26, -6.10},
        {500, -25, 4.5, 25, 6.5, -3.73},
        {500, -15, 6.5, 5, 1.2, -26.98},
        {500, -15, 6.5, 5, 4.8, -2.92},
        {500, -15, 6.5, 25, 1.2, -6.93},
        {500, -15, 6.5, 25, 4.8, -0.34},
        {500, -10, 2, 5, 1.2, -37.68},
        {500, -10, 2, 5, 4.8, -15.04},
        {500, -10, 2, 5, 5.26, -14.99},
        {500, -10, 2, 5, 6.5, -8.12},
        {500, -10, 2, 25, 1.2, -19.90},
        {500, -10, 2, 25, 4.8, -15.47},
        {500, -10, 2, 25, 5.26, -14.75},
        {500, -10, 2, 25, 6.5, -14.35},
        {1000, -25, 4.5, 5, 1.2, -19.87},
        {1000, -25, 4.5, 5, 4.8, -10.00},
        {1000, -25, 4.5, 5, 5.26, -7.61},
        {1000, -25, 4.5, 5, 6.5, -1.61},
        {1000, -25, 4.5, 25, 1.2, -24.09},
        {1000, -25, 4.5, 25, 4.8, -10.02},
        {1000, -25, 4.5, 25, 5.26, -8.09},
        {1000, -25, 4.5, 25, 6.5, -5.26},
        {1000, -15, 6.5, 5, 1.2, -15.40},
        {1000, -15, 6.5, 5, 4.8, -4.28},
        {1000, -15, 6.5, 25, 1.2, -16.63},
        {1000, -15, 6.5, 25, 4.8, -1.34},
        {1000, -10, 2, 5, 1.2, -25.89},
        {1000, -10, 2, 5, 4.8, -17.77},
        {1000, -10, 2, 5, 5.26, -16.37},
        {1000, -10, 2, 5, 6.5, -9.36},
        {1000, -10, 2, 25, 1.2, -31.45},
        {1000, -10, 2, 25, 4.8, -20.50},
        {1000, -10, 2, 25, 5.26, -17.71},
        {1000, -10, 2, 25, 6.5, -18.45},
    };
    return cases;
}

} // namespace

TEST_CASE("over a long block the point-source field is the boundary-element solution's", "[point_paths]")
{
    // A block 2.44 m thick and 5.18 m high -- two stacked shipping containers,
    // DigiCity's (FIDELITY) -- on rigid ground, long enough to be two-
    // dimensional, the source and the microphone in one plane across it.
    // There every ray is a 2-D ray times 1/sqrt(its length), and the field in
    // that plane is the line source's: the exact values are a boundary-element
    // solution for one (tools/thick_barrier_bem.py, 60 elements a wavelength,
    // converged to 0.05 dB here), in dB re the free field.  Sources low, near
    // grazing and above the roof; microphones deep in the shadow, at the roof
    // plane and above it.  The uniform form is within 1 dB everywhere; the
    // cascade misses by up to 10.
    const std::vector<BlockCase>& cases = long_block_cases();
    const double width = 2.44, height = 5.18;
    Footprint block;
    block.x = {-2000.0, 2000.0, 2000.0, -2000.0};
    block.y = {0.0, 0.0, width, width};
    block.top = height;
    const Buildings city({block});
    const std::vector<double> reflection(city.facades().size(), 1.0);
    double worst[2] {0.0, 0.0}, square = 0.0;
    for (const BlockCase& c : cases)
    {
        const double k = 2.0 * pi * c.f / 343.2;
        const Vec3 s(0.0, c.source_y, c.source_z), m(0.0, width + c.behind, c.z);
        const double free = 1.0 / std::sqrt((m - s).norm());
        for (int uniform = 0; uniform < 2; uniform++)
        {
            PointPathOptions options;
            options.shared_face = uniform == 1;
            const PointSourceField field(city, s, m, options);
            std::complex<double> sum(0.0, 0.0);
            for (const PointRay& ray : field.rays())
                sum += point_ray_amplitude(ray, reflection, 1.0, k) * std::sqrt(ray.length)
                        * std::polar(1.0, k * ray.length);
            const double error = 20.0 * std::log10(std::abs(sum) / free) - c.exact;
            worst[uniform] = std::max(worst[uniform], std::abs(error));
            if (uniform == 1)
            {
                square += error * error;
                INFO(c.f << " Hz, source " << c.source_y << " m, " << c.source_z << " m up, microphone " << c.behind
                          << " m behind, " << c.z << " m up");
                REQUIRE(std::abs(error) < 1.0);
            }
        }
    }
    REQUIRE(std::sqrt(square / cases.size()) < 0.3);
    REQUIRE(worst[0] > 5.0);
}

// The same block with the edge-source integral (edge_source.hpp): the field
// that crosses the roof is its second order over the shared face, which holds
// the transition the uniform form approximates and the cascade misses.  Every
// other case.  Within 0.28 dB, rms 0.10, where the uniform form is within
// 0.94 dB, rms 0.21; the BEM's own convergence is 0.05 dB and the point
// source's 2-D equivalent (each ray times the root of its length) is good to
// O(1/kL), most of what is left at 125 Hz.
TEST_CASE("over a long block the edge-source field is the boundary-element solution's", "[point_paths][edge_source]")
{
    const double width = 2.44, height = 5.18;
    Footprint block;
    block.x = {-2000.0, 2000.0, 2000.0, -2000.0};
    block.y = {0.0, 0.0, width, width};
    block.top = height;
    const Buildings city({block});
    const std::vector<double> reflection(city.facades().size(), 1.0);
    double worst = 0.0, square = 0.0;
    int count = 0;
    const std::vector<BlockCase>& cases = long_block_cases();
    for (std::size_t i = 0; i < cases.size(); i += 2)
    {
        const BlockCase& c = cases[i];
        const double k = 2.0 * pi * c.f / 343.2;
        const Vec3 s(0.0, c.source_y, c.source_z), m(0.0, width + c.behind, c.z);
        const double free = 1.0 / std::sqrt((m - s).norm());
        PointPathOptions options;
        options.edge_source = true;
        const PointSourceField field(city, s, m, options);
        std::complex<double> sum(0.0, 0.0);
        for (const PointRay& ray : field.rays())
        {
            const std::vector<std::complex<double>> f = edge_source_factors(city, ray, {k});
            const std::complex<double> a = f.empty() ? point_ray_amplitude(ray, reflection, 1.0, k)
                                                     : reflected_amplitude(ray, reflection, 1.0) * f[0];
            sum += a * std::sqrt(ray.length) * std::polar(1.0, k * ray.length);
        }
        const double error = 20.0 * std::log10(std::abs(sum) / free) - c.exact;
        INFO(c.f << " Hz, source " << c.source_y << " m, " << c.source_z << " m up, microphone " << c.behind
                  << " m behind, " << c.z << " m up");
        REQUIRE(std::abs(error) < 0.35);
        worst = std::max(worst, std::abs(error));
        square += error * error;
        count++;
    }
    INFO("worst " << worst << " dB");
    REQUIRE(std::sqrt(square / count) < 0.15);
}

namespace {

// A step: a long block 10 m across and 10 m high, then one 10 m across and
// 7 m high against it.  The boundary-element solutions (tools/thick_barrier_bem.py
// --width 20 --height 7 --roof '0,10;10,10;10,7', 45 elements a wavelength,
// converged to 0.15 dB against 30), BlockCase's columns, the microphone's
// distance behind the lower block; none on the line through the two blocks'
// far edges, where the uniform form's level stepped before the jump split
// (FIDELITY, "Across a step").
const std::vector<BlockCase>& step_cases()
{
    static const std::vector<BlockCase> cases {
        {125, -100, 10.5, 2, 1.2, -15.69},
        {250, -100, 10.5, 2, 1.2, -31.37},
        {125, -100, 10.5, 2, 4, -15.20},
        {250, -100, 10.5, 2, 4, -33.29},
        {500, -100, 10.5, 2, 4, -24.19},
        {125, -100, 10.5, 5, 1.2, -17.64},
        {250, -100, 10.5, 5, 1.2, -36.33},
        {500, -100, 10.5, 5, 1.2, -25.87},
        {125, -100, 10.5, 5, 2, -14.84},
        {250, -100, 10.5, 5, 2, -31.09},
        {500, -100, 10.5, 5, 2, -28.46},
        {125, -100, 10.5, 5, 3, -15.08},
        {250, -100, 10.5, 5, 3, -33.67},
        {500, -100, 10.5, 5, 3, -23.36},
        {125, -100, 10.5, 5, 4, -15.39},
        {250, -100, 10.5, 5, 4, -28.41},
        {500, -100, 10.5, 5, 4, -22.33},
        {125, -100, 10.5, 5, 5, -11.03},
        {250, -100, 10.5, 5, 5, -24.25},
        {500, -100, 10.5, 5, 5, -19.56},
        {125, -100, 10.5, 5, 6, -10.68},
        {250, -100, 10.5, 5, 6, -20.51},
        {500, -100, 10.5, 5, 6, -15.14},
        {125, -100, 10.5, 5, 7, -7.24},
        {250, -100, 10.5, 5, 7, -16.53},
        {500, -100, 10.5, 5, 7, -12.58},
        {125, -100, 10.5, 10, 1.2, -28.02},
        {250, -100, 10.5, 10, 1.2, -22.62},
        {500, -100, 10.5, 10, 1.2, -21.94},
        {125, -100, 10.5, 10, 2, -12.14},
        {250, -100, 10.5, 10, 2, -26.72},
        {500, -100, 10.5, 10, 2, -26.24},
        {125, -100, 10.5, 10, 3, -13.65},
        {250, -100, 10.5, 10, 3, -28.73},
        {500, -100, 10.5, 10, 3, -19.90},
        {125, -100, 10.5, 10, 5, -8.63},
        {250, -100, 10.5, 10, 5, -17.94},
        {500, -100, 10.5, 10, 5, -15.30},
        {125, -100, 10.5, 10, 6, -10.81},
        {250, -100, 10.5, 10, 6, -16.46},
        {500, -100, 10.5, 10, 6, -14.55},
        {125, -100, 10.5, 10, 7, -6.35},
        {250, -100, 10.5, 10, 7, -15.76},
        {500, -100, 10.5, 10, 7, -10.88},
        {125, -60, 12, 5, 1.2, -23.26},
        {250, -60, 12, 5, 1.2, -29.28},
        {500, -60, 12, 5, 1.2, -26.15},
        {125, -60, 12, 5, 4, -20.89},
        {250, -60, 12, 5, 4, -21.35},
        {500, -60, 12, 5, 4, -22.60},
        {125, -60, 12, 10, 1.2, -33.14},
        {250, -60, 12, 10, 1.2, -15.50},
        {500, -60, 12, 10, 1.2, -22.13},
        {125, -30, 15, 2, 1.2, -14.40},
        {250, -30, 15, 2, 1.2, -22.77},
        {500, -30, 15, 2, 1.2, -27.48},
        {125, -30, 15, 2, 4, -13.97},
        {250, -30, 15, 2, 4, -24.82},
        {500, -30, 15, 2, 4, -18.74},
        {125, -30, 15, 5, 1.2, -16.19},
        {250, -30, 15, 5, 1.2, -27.71},
        {500, -30, 15, 5, 1.2, -20.24},
        {125, -30, 15, 5, 2, -13.40},
        {250, -30, 15, 5, 2, -22.41},
        {500, -30, 15, 5, 2, -22.82},
        {125, -30, 15, 5, 3, -13.59},
        {250, -30, 15, 5, 3, -24.89},
        {500, -30, 15, 5, 3, -17.70},
        {125, -30, 15, 5, 4, -13.92},
        {250, -30, 15, 5, 4, -19.53},
        {500, -30, 15, 5, 4, -16.63},
        {125, -30, 15, 5, 5, -9.48},
        {250, -30, 15, 5, 5, -15.34},
        {500, -30, 15, 5, 5, -13.78},
        {125, -30, 15, 5, 6, -9.20},
        {250, -30, 15, 5, 6, -11.57},
        {500, -30, 15, 5, 6, -9.20},
        {125, -30, 15, 5, 7, -5.69},
        {250, -30, 15, 5, 7, -7.82},
        {500, -30, 15, 5, 7, -6.33},
        {125, -30, 15, 10, 1.2, -25.72},
        {250, -30, 15, 10, 1.2, -13.60},
        {500, -30, 15, 10, 1.2, -16.02},
        {125, -30, 15, 10, 2, -10.26},
        {250, -30, 15, 10, 2, -17.54},
        {500, -30, 15, 10, 2, -20.19},
        {125, -30, 15, 10, 3, -11.88},
        {250, -30, 15, 10, 3, -19.31},
        {500, -30, 15, 10, 3, -13.91},
        {125, -30, 15, 10, 5, -6.78},
        {250, -30, 15, 10, 5, -8.77},
        {500, -30, 15, 10, 5, -8.95},
        {125, -30, 15, 10, 6, -8.87},
        {250, -30, 15, 10, 6, -7.11},
        {500, -30, 15, 10, 6, -7.76},
        {125, -30, 15, 10, 7, -4.41},
        {250, -30, 15, 10, 7, -6.64},
        {500, -30, 15, 10, 7, -4.01},
    };
    return cases;
}

// The line source's level, dB re the free field, over long buildings: every
// ray times the root of its length.
double line_source_level(const Buildings& city, const Vec3& s, const Vec3& m, double f, const PointPathOptions& options)
{
    const double k = 2.0 * pi * f / 343.2;
    const std::vector<double> reflection(city.facades().size(), 1.0);
    const PointSourceField field(city, s, m, options);
    std::complex<double> sum(0.0, 0.0);
    for (const PointRay& ray : field.rays())
        sum += point_ray_amplitude(ray, reflection, 1.0, k) * std::sqrt(ray.length) * std::polar(1.0, k * ray.length);
    return 20.0 * std::log10(std::abs(sum) * std::sqrt((m - s).norm()));
}

Footprint long_block(double y0, double y1, double top)
{
    Footprint b;
    b.x = {-4000.0, 4000.0, 4000.0, -4000.0};
    b.y = {y0, y0, y1, y1};
    b.top = top;
    return b;
}

} // namespace

// Over the step the sound from above crosses the taller block's far edge and
// then the lower one's, which share no face: the line between them touches
// both wedges from above, the lower roof tilted 16.7 degrees below it, and the
// aperture over the pair (double_wedge.hpp) takes the transition the cascade
// loses.  98 cases, sources 10.5-15 m up and 30-100 m off, microphones 2-10 m
// behind and 1.2-7 m up, 125-500 Hz: within 1.8 dB, rms 0.48, with the jump
// split at the line through the far edges (0.53 and 2.2 without), where the
// cascade is rms 0.80 and its mean 0.37 dB low.
TEST_CASE("over a step down to a lower roof the pair of edges is uniform", "[point_paths]")
{
    const Buildings city({long_block(0.0, 10.0, 10.0), long_block(10.0, 20.0, 7.0)});
    double square[2] {0.0, 0.0}, mean[2] {0.0, 0.0}, worst = 0.0;
    const std::vector<BlockCase>& cases = step_cases();
    for (const BlockCase& c : cases)
    {
        const Vec3 s(0.0, c.source_y, c.source_z), m(0.0, 20.0 + c.behind, c.z);
        for (int uniform = 0; uniform < 2; uniform++)
        {
            PointPathOptions options;
            options.shared_face = uniform == 1;
            const double error = line_source_level(city, s, m, c.f, options) - c.exact;
            square[uniform] += error * error;
            mean[uniform] += error;
            if (uniform == 1)
            {
                INFO(c.f << " Hz, source " << c.source_y << " m, " << c.source_z << " m up, microphone " << c.behind
                          << " m behind, " << c.z << " m up");
                REQUIRE(std::abs(error) < 2.5);
                worst = std::max(worst, std::abs(error));
            }
        }
    }
    const double n = static_cast<double>(cases.size());
    INFO("uniform rms " << std::sqrt(square[1] / n) << " mean " << mean[1] / n << " worst " << worst << "; cascade rms "
                        << std::sqrt(square[0] / n) << " mean " << mean[0] / n);
    REQUIRE(std::sqrt(square[1] / n) < 0.65);
    REQUIRE(std::abs(mean[1] / n) < 0.25);
    REQUIRE(std::sqrt(square[0] / n) > std::sqrt(square[1] / n) + 0.2);
}

// The same step, microphones above the lower roof (8-11 m up, 5 and 10 m
// behind; 45 elements a wavelength, converged to 0.01 dB against 30), among
// them two boundaries: the lower roof's reflection of the first edge's wave
// (8.5 m up at 5 m behind, 10 m at 10) and the taller roof's plane (10 m up).
// A microphone exactly on either took the waves of both sides at once: in
// the taller roof's plane the pair took the cascade, as above it, while the
// ray over the near edge alone stayed blocked, as below, 3-9.6 dB high; on
// the reflection boundary the roof's reflected ray and the pair's aperture
// both, 13 dB at a point.  Each boundary now goes to the side the ray set
// gives it.  Within 1 dB, rms 0.26.
TEST_CASE("above a step's lower roof the field holds to its boundaries", "[point_paths]")
{
    const Buildings city({long_block(0.0, 10.0, 10.0), long_block(10.0, 20.0, 7.0)});
    static const std::vector<BlockCase> cases {
        {125, -100, 10.5, 5, 8.00, -4.87}, {125, -100, 10.5, 5, 8.40, -5.24}, {125, -100, 10.5, 5, 8.60, -5.41}, {125, -100, 10.5, 5, 9.00, -5.29},
        {125, -100, 10.5, 5, 10.00, -3.89}, {125, -100, 10.5, 5, 11.00, -3.58}, {125, -100, 10.5, 10, 8.00, -4.97}, {125, -100, 10.5, 10, 8.40, -5.37},
        {125, -100, 10.5, 10, 8.60, -5.52}, {125, -100, 10.5, 10, 9.00, -5.30}, {125, -100, 10.5, 10, 10.00, -2.88}, {125, -100, 10.5, 10, 11.00, -2.38},
        {125, -30, 15, 5, 8.00, -3.75}, {125, -30, 15, 5, 8.40, -4.33}, {125, -30, 15, 5, 8.60, -4.58}, {125, -30, 15, 5, 9.00, -4.51},
        {125, -30, 15, 5, 10.00, -3.50}, {125, -30, 15, 5, 11.00, -3.01}, {125, -30, 15, 10, 8.00, -3.63}, {125, -30, 15, 10, 8.40, -4.22},
        {125, -30, 15, 10, 8.60, -4.41}, {125, -30, 15, 10, 9.00, -4.14}, {125, -30, 15, 10, 10.00, -2.01}, {125, -30, 15, 10, 11.00, -2.50},
        {250, -100, 10.5, 5, 8.00, -13.71}, {250, -100, 10.5, 5, 8.40, -14.43}, {250, -100, 10.5, 5, 8.60, -13.76}, {250, -100, 10.5, 5, 9.00, -12.15},
        {250, -100, 10.5, 5, 10.00, -6.75}, {250, -100, 10.5, 5, 11.00, -4.89}, {250, -100, 10.5, 10, 8.00, -12.84}, {250, -100, 10.5, 10, 8.40, -12.96},
        {250, -100, 10.5, 10, 8.60, -12.27}, {250, -100, 10.5, 10, 9.00, -10.63}, {250, -100, 10.5, 10, 10.00, -8.44}, {250, -100, 10.5, 10, 11.00, -4.68},
        {250, -30, 15, 5, 8.00, -5.51}, {250, -30, 15, 5, 8.40, -6.18}, {250, -30, 15, 5, 8.60, -5.36}, {250, -30, 15, 5, 9.00, -3.84},
        {250, -30, 15, 5, 10.00, 1.11}, {250, -30, 15, 5, 11.00, 0.14}, {250, -30, 15, 10, 8.00, -4.43}, {250, -30, 15, 10, 8.40, -4.48},
        {250, -30, 15, 10, 8.60, -3.70}, {250, -30, 15, 10, 9.00, -2.35}, {250, -30, 15, 10, 10.00, -0.53}, {250, -30, 15, 10, 11.00, 1.58},
        {500, -100, 10.5, 5, 8.00, -11.14}, {500, -100, 10.5, 5, 8.40, -7.97}, {500, -100, 10.5, 5, 8.60, -7.25}, {500, -100, 10.5, 5, 9.00, -5.80},
        {500, -100, 10.5, 5, 10.00, -3.96}, {500, -100, 10.5, 5, 11.00, -2.24}, {500, -100, 10.5, 10, 8.00, -9.85}, {500, -100, 10.5, 10, 8.40, -7.90},
        {500, -100, 10.5, 10, 8.60, -7.38}, {500, -100, 10.5, 10, 9.00, -6.27}, {500, -100, 10.5, 10, 10.00, -4.38}, {500, -100, 10.5, 10, 11.00, -1.22},
        {500, -30, 15, 5, 8.00, -4.46}, {500, -30, 15, 5, 8.40, -1.39}, {500, -30, 15, 5, 8.60, -0.85}, {500, -30, 15, 5, 9.00, 0.13},
        {500, -30, 15, 5, 10.00, 2.43}, {500, -30, 15, 5, 11.00, 0.58}, {500, -30, 15, 10, 8.00, -2.71}, {500, -30, 15, 10, 8.40, -0.79},
        {500, -30, 15, 10, 8.60, -0.50}, {500, -30, 15, 10, 9.00, 0.52}, {500, -30, 15, 10, 10.00, 0.98}, {500, -30, 15, 10, 11.00, 2.79},
    };
    double square = 0.0;
    for (const BlockCase& c : cases)
    {
        const Vec3 s(0.0, c.source_y, c.source_z), m(0.0, 20.0 + c.behind, c.z);
        const double error = line_source_level(city, s, m, c.f, PointPathOptions()) - c.exact;
        INFO(c.f << " Hz, source " << c.source_y << " m, " << c.source_z << " m up, microphone " << c.behind
                  << " m behind, " << c.z << " m up");
        REQUIRE(std::abs(error) < 1.1);
        square += error * error;
    }
    REQUIRE(std::sqrt(square / static_cast<double>(cases.size())) < 0.35);
    // On each boundary the level is its side's: the plane's below it, the
    // reflection boundary's on the reflected side.
    for (const double f : {125.0, 250.0, 500.0})
    {
        const Vec3 s(0.0, -100.0, 10.5);
        auto at = [&](double behind, double z) {
            return line_source_level(city, s, Vec3(0.0, 20.0 + behind, z), f, PointPathOptions());
        };
        INFO(f << " Hz");
        REQUIRE(at(5.0, 10.0) == Catch::Approx(at(5.0, 10.0 - 1e-6)).margin(0.01));
        REQUIRE(at(5.0, 8.5) == Catch::Approx(at(5.0, 8.5 + 1e-6)).margin(0.01));
    }
}

// The same step, the source and the microphone both above the taller roof's
// plane: the source 0.3 and 1.9 degrees above it, the microphone 0.1-2 m up
// from it, 5 and 15 m behind (boundary elements at 45 elements a wavelength,
// within 0.09 dB of 30).  Neither end is below the line between the roof's
// edges, and the pair takes each end's mirror below it, negated
// (double_wedge.cpp, mirrored_weights): within 0.3 dB, rms 0.10, where the
// cascade is 1.6 dB rms and up to 4.3 dB low.  Through the plane the level
// holds to 0.11 dB, where the cascade stepped up to 8 dB.
TEST_CASE("over a roof with both ends above its plane the field holds through it", "[point_paths]")
{
    const Buildings city({long_block(0.0, 10.0, 10.0), long_block(10.0, 20.0, 7.0)});
    static const std::vector<BlockCase> cases {
            {125, -100, 10.5, 5, 10.1, -3.84}, {125, -100, 10.5, 5, 10.25, -3.86}, {125, -100, 10.5, 5, 10.5, -3.94}, {125, -100, 10.5, 5, 11, -3.52},
            {125, -100, 10.5, 5, 12, 1.10}, {125, -100, 10.5, 15, 10.1, -2.57}, {125, -100, 10.5, 15, 10.25, -2.63}, {125, -100, 10.5, 15, 10.5, -2.77},
            {125, -100, 10.5, 15, 11, -2.94}, {125, -100, 10.5, 15, 12, -1.36}, {125, -60, 12, 5, 10.1, -8.69}, {125, -60, 12, 5, 10.25, -8.65},
            {125, -60, 12, 5, 10.5, -8.60}, {125, -60, 12, 5, 11, -7.76}, {125, -60, 12, 5, 12, -3.30}, {125, -60, 12, 15, 10.1, -7.06},
            {125, -60, 12, 15, 10.25, -7.11}, {125, -60, 12, 15, 10.5, -7.20}, {125, -60, 12, 15, 11, -7.11}, {125, -60, 12, 15, 12, -5.17},
            {250, -100, 10.5, 5, 10.1, -6.18}, {250, -100, 10.5, 5, 10.25, -5.52}, {250, -100, 10.5, 5, 10.5, -4.92}, {250, -100, 10.5, 5, 11, -4.89},
            {250, -100, 10.5, 5, 12, -3.23}, {250, -100, 10.5, 15, 10.1, -7.65}, {250, -100, 10.5, 15, 10.25, -7.47}, {250, -100, 10.5, 15, 10.5, -7.06},
            {250, -100, 10.5, 15, 11, -5.43}, {250, -100, 10.5, 15, 12, -2.80}, {250, -60, 12, 5, 10.1, -0.50}, {250, -60, 12, 5, 10.25, 0.08},
            {250, -60, 12, 5, 10.5, 0.47}, {250, -60, 12, 5, 11, -0.19}, {250, -60, 12, 5, 12, 0.09}, {250, -60, 12, 15, 10.1, -2.76},
            {250, -60, 12, 15, 10.25, -2.71}, {250, -60, 12, 15, 10.5, -2.48}, {250, -60, 12, 15, 11, -0.94}, {250, -60, 12, 15, 12, 1.06},
            {500, -100, 10.5, 5, 10.1, -2.90}, {500, -100, 10.5, 5, 10.25, -1.56}, {500, -100, 10.5, 5, 10.5, -0.77}, {500, -100, 10.5, 5, 11, -2.23},
            {500, -100, 10.5, 5, 12, 1.71}, {500, -100, 10.5, 15, 10.1, -3.75}, {500, -100, 10.5, 15, 10.25, -3.81}, {500, -100, 10.5, 15, 10.5, -3.50},
            {500, -100, 10.5, 15, 11, -2.27}, {500, -100, 10.5, 15, 12, 0.29}, {500, -60, 12, 5, 10.1, -2.36}, {500, -60, 12, 5, 10.25, -1.04},
            {500, -60, 12, 5, 10.5, -0.19}, {500, -60, 12, 5, 11, -1.15}, {500, -60, 12, 5, 12, 2.85}, {500, -60, 12, 15, 10.1, -2.76},
            {500, -60, 12, 15, 10.25, -2.74}, {500, -60, 12, 15, 10.5, -2.34}, {500, -60, 12, 15, 11, -1.03}, {500, -60, 12, 15, 12, 1.64},
            {1000, -100, 10.5, 5, 10.1, -5.21}, {1000, -100, 10.5, 5, 10.25, -2.97}, {1000, -100, 10.5, 5, 10.5, -1.54}, {1000, -100, 10.5, 5, 11, -0.75},
            {1000, -100, 10.5, 5, 12, 1.11}, {1000, -100, 10.5, 15, 10.1, -4.49}, {1000, -100, 10.5, 15, 10.25, -4.30}, {1000, -100, 10.5, 15, 10.5, -3.14},
            {1000, -100, 10.5, 15, 11, -1.18}, {1000, -100, 10.5, 15, 12, 0.25}, {1000, -60, 12, 5, 10.1, -4.21}, {1000, -60, 12, 5, 10.25, -2.20},
            {1000, -60, 12, 5, 10.5, -1.07}, {1000, -60, 12, 5, 11, 0.02}, {1000, -60, 12, 5, 12, 1.45}, {1000, -60, 12, 15, 10.1, -3.30},
            {1000, -60, 12, 15, 10.25, -3.02}, {1000, -60, 12, 15, 10.5, -1.89}, {1000, -60, 12, 15, 11, -0.15}, {1000, -60, 12, 15, 12, 1.18},
    };
    double square = 0.0, cascade_square = 0.0;
    for (const BlockCase& c : cases)
    {
        const Vec3 s(0.0, c.source_y, c.source_z), m(0.0, 20.0 + c.behind, c.z);
        const double error = line_source_level(city, s, m, c.f, PointPathOptions()) - c.exact;
        PointPathOptions cascade;
        cascade.shared_face = false;
        const double cascade_error = line_source_level(city, s, m, c.f, cascade) - c.exact;
        INFO(c.f << " Hz, source " << c.source_y << " m, " << c.source_z << " m up, microphone " << c.behind
                  << " m behind, " << c.z << " m up");
        REQUIRE(std::abs(error) < 0.5);
        square += error * error;
        cascade_square += cascade_error * cascade_error;
    }
    const double n = static_cast<double>(cases.size());
    INFO("rms " << std::sqrt(square / n) << ", cascade " << std::sqrt(cascade_square / n));
    REQUIRE(std::sqrt(square / n) < 0.15);
    REQUIRE(std::sqrt(cascade_square / n) > 1.0);
    for (const double f : {125.0, 250.0, 500.0})
        for (const double source_z : {10.5, 12.0, 15.0})
            for (const double behind : {5.0, 15.0})
            {
                const Vec3 s(0.0, source_z == 10.5 ? -100.0 : -60.0, source_z);
                auto at = [&](double z) {
                    return line_source_level(city, s, Vec3(0.0, 20.0 + behind, z), f, PointPathOptions());
                };
                INFO(f << " Hz, source " << source_z << " m up, " << behind << " m behind");
                REQUIRE(at(10.0 + 1e-6) == Catch::Approx(at(10.0)).margin(0.2));
            }
}

// EDGE's pent-roofed terrace (UE02) at model scale: a ridge 0.732 m up over a
// slope from the eaves, a step down to a deck 0.612 m up, the deck's edge,
// and a 45-degree slope to the far eave.  The drone's sound (1.5 m up, 2 m
// before the terrace) and the floor source's (2.9 cm up, 1.5 m before)
// cross the ridge and then the deck's edge, 1.2 degrees inside the ridge's
// shadow boundary from the drone; the microphone on the floor 1.25 m behind.
// Boundary elements (45 elements a wavelength, converged to 0.1 dB): the
// uniform form is within 0.6 dB, where the cascade is 4-5 dB low from the
// drone -- the bias EDGE's measurements show there (FIDELITY, "Pitched
// roofs").
TEST_CASE("EDGE's pent-roofed terrace takes its ridge and deck edge uniformly", "[point_paths]")
{
    Footprint ridge = long_block(0.0, 0.48, 0.732);
    ridge.roof = {{0.0, 0.5, 0.492}};
    Footprint deck = long_block(0.48, 0.96, 0.612);
    deck.roof = {{0.0, -1.0, 1.452}, {0.0, 0.0, 0.612}};
    const Buildings city({ridge, deck});
    struct Case
    {
        double f, sy, sz, exact;
    };
    const Case cases[] {{1600, -2.0, 1.5, -2.66}, {3150, -2.0, 1.5, -3.54}, {6300, -2.0, 1.5, -5.35},
            {1600, -1.5, 0.02925, -7.78}, {3150, -1.5, 0.02925, -11.81}, {6300, -1.5, 0.02925, -20.95}};
    double cascade_worst = 0.0;
    for (const Case& c : cases)
    {
        const Vec3 s(0.0, c.sy, c.sz), m(0.0, 0.96 + 1.25, 0.00318);
        INFO(c.f << " Hz, source " << c.sz << " m up");
        REQUIRE(std::abs(line_source_level(city, s, m, c.f, PointPathOptions()) - c.exact) < 0.6);
        PointPathOptions cascade;
        cascade.shared_face = false;
        cascade_worst = std::max(cascade_worst, std::abs(line_source_level(city, s, m, c.f, cascade) - c.exact));
    }
    REQUIRE(cascade_worst > 4.0);
}

// A microphone on the line through the step's far edges is seen from the
// first past the second -- the path search traces the first's ray grazing
// past -- so the aperture takes that ray's field back out there as it does
// just above.  Taking the microphone as below the line there counted the
// field twice, 9-13 dB high.
TEST_CASE("on the line through a stepped pair the field is the field just above it", "[point_paths]")
{
    const Buildings city({long_block(0.0, 10.0, 10.0), long_block(10.0, 20.0, 7.0)});
    const Vec3 s(0.0, -100.0, 10.5);
    for (double f : {125.0, 250.0, 500.0})
    {
        INFO(f << " Hz");
        const double on = line_source_level(city, s, Vec3(0.0, 30.0, 4.0), f, PointPathOptions());
        const double above = line_source_level(city, s, Vec3(0.0, 30.0, 4.001), f, PointPathOptions());
        REQUIRE(std::abs(on - above) < 0.05);
    }
}

namespace {

// A source `above` degrees above the line through the step's far edges,
// `range` m from the taller block's, and a point `off` m off the line along
// its normal, at `behind` m past the lower block.
Vec3 above_step_line(double range, double above)
{
    const double a = std::atan(0.3) + above * pi / 180.0;
    return Vec3(0.0, 10.0 - range * std::cos(a), 10.0 + range * std::sin(a));
}
Vec3 off_step_line(double behind, double off)
{
    const double t = std::atan(0.3), y = 20.0 + behind;
    return Vec3(0.0, y + off * std::sin(t), 10.0 - 0.3 * (y - 10.0) + off * std::cos(t));
}

} // namespace

// Where the microphone crosses the line through the step's far edges the
// rays change: below it the pair is the receiver's aperture and the first
// edge's ray is blocked, above it the first edge's ray appears beside the
// cascade, and the level stepped by up to 6 dB with the source near the
// line.  The jump split (TiltedStepSplit) shares the jump between the sides:
// sources 0.5, 3.3 and 10 degrees above the line and 16 below, 20 and 60 m
// from the taller block's edge, microphones 2, 10 and 20 m past the lower
// block, 0.1 mm either side, 125 Hz to 1 kHz: within 0.001 dB, and the band
// levels as close.  At 0.1 mm the microphone is on the foot fields' own
// points, where the carry and the fade are 1 and the sides agree for any
// beta: that pins the sides' shares and signs.  At 1 mm, where the carry and
// the fade are the microphone's own, the field changes across the line by no
// more than it does over 2 mm on either side alone.  Neither says what beta
// or the fade should be: that is the boundary-element tables' to say.
TEST_CASE("across the line through a stepped pair the field does not step", "[point_paths]")
{
    const Buildings city({long_block(0.0, 10.0, 10.0), long_block(10.0, 20.0, 7.0)});
    PointPathOptions unsplit;
    unsplit.tilted_step_split = false;
    double worst = 0.0, unsplit_worst = 0.0;
    for (double range : {60.0, 20.0})
        for (double above : {0.5, 3.3, 10.0, -16.0})
            for (double behind : {2.0, 10.0, 20.0})
                for (double f : {125.0, 250.0, 500.0, 1000.0})
                {
                    const Vec3 s = above_step_line(range, above);
                    const Vec3 lo = off_step_line(behind, -1e-4), hi = off_step_line(behind, 1e-4);
                    const double step = line_source_level(city, s, hi, f, split_options())
                            - line_source_level(city, s, lo, f, split_options());
                    INFO(f << " Hz, source " << above << " degrees above the line, " << range << " m, microphone "
                           << behind << " m behind");
                    REQUIRE(std::abs(step) < 0.02);
                    worst = std::max(worst, std::abs(step));
                    // At 1 mm either side, no more than the field changes
                    // over 2 mm on either side alone (0.9 of it, up to 0.16
                    // dB at 1 kHz).
                    double at[4];
                    for (int i = 0; i < 4; i++)
                        at[i] = line_source_level(city, s, off_step_line(behind, (2 * i - 3) * 1e-3), f, split_options());
                    const double alone = 0.5 * (std::abs(at[1] - at[0]) + std::abs(at[3] - at[2]));
                    REQUIRE(std::abs(at[2] - at[1]) < alone + 0.002);
                    unsplit_worst = std::max(unsplit_worst, std::abs(line_source_level(city, s, hi, f, unsplit)
                                                                    - line_source_level(city, s, lo, f, unsplit)));
                }
    INFO("worst " << worst << " dB, unsplit " << unsplit_worst);
    REQUIRE(unsplit_worst > 2.0);
    // The band levels too: a point source, a rigid floor and walls absorbing.
    const Vec3 s = above_step_line(60.0, 3.3);
    for (double behind : {2.0, 10.0})
    {
        const std::vector<double> below = levels(city, s, off_step_line(behind, -1e-4), 0.9, 1.0, split_options());
        const std::vector<double> over = levels(city, s, off_step_line(behind, 1e-4), 0.9, 1.0, split_options());
        for (std::size_t b = 0; b < bands.size(); b++)
        {
            INFO(bands[b] << " Hz, " << behind << " m behind");
            REQUIRE(std::abs(over[b] - below[b]) < 0.02);
        }
    }
}

// Blocks 40 m long and the source 15 m along them, so the pair's rays reach
// microphones past the blocks' ends: the split goes on as far as they do,
// not only where the microphone's foot is within the edges.  It stopped at
// the ends, and the field stepped 1.0-2.1 dB across the line 1-2 m past them.
TEST_CASE("past a stepped pair's ends the split goes on with its rays", "[point_paths]")
{
    const Buildings city({box(-20.0, 0.0, 20.0, 10.0, 10.0), box(-20.0, 10.0, 20.0, 20.0, 7.0)});
    Vec3 s = above_step_line(60.0, 3.3);
    s.x = -15.0;
    for (double x : {19.0, 21.0, 22.0})
    {
        Vec3 lo = off_step_line(5.0, -1e-3), hi = off_step_line(5.0, 1e-3);
        lo.x = hi.x = x;
        const std::vector<double> below = levels(city, s, lo, 0.9, 1.0, split_options()), over = levels(city, s, hi, 0.9, 1.0, split_options());
        for (std::size_t b = 0; b < 3; b++)
        {
            INFO(bands[b] << " Hz, microphone at x = " << x);
            REQUIRE(std::abs(over[b] - below[b]) < 0.05);
        }
    }
}

// Where the line through the step's far edges goes into the floor, 23.3 m
// past the lower block, a microphone 1.2 m up still stands by it, and its
// foot slides back along the line to where the line leaves the air.  The
// split was dropped there, and the level stepped by up to 1.8 dB onto the
// rays alone.  Boundary elements, 45 elements a wavelength: the split within
// 0.8 dB on both sides of that point; the rays alone 1.1-2.4 dB out.
TEST_CASE("where a stepped pair's line meets the floor the split goes on", "[point_paths]")
{
    const Buildings city({long_block(0.0, 10.0, 10.0), long_block(10.0, 20.0, 7.0)});
    static const std::vector<BlockCase> cases {
            {125, -46.3818, 30.5205, 22.8, 1.2, -3.48}, {125, -46.3818, 30.5205, 23.6, 1.2, -2.99},
            {125, -46.3818, 30.5205, 23.7, 1.2, -2.93}, {125, -46.3818, 30.5205, 24.4, 1.2, -2.53},
            {250, -46.3818, 30.5205, 22.8, 1.2, -6.25}, {250, -46.3818, 30.5205, 23.6, 1.2, -6.59},
            {250, -46.3818, 30.5205, 23.7, 1.2, -6.63}, {250, -46.3818, 30.5205, 24.4, 1.2, -6.96},
    };
    PointPathOptions unsplit;
    unsplit.tilted_step_split = false;
    for (const BlockCase& c : cases)
    {
        const Vec3 s(0.0, c.source_y, c.source_z), m(0.0, 20.0 + c.behind, c.z);
        INFO(c.f << " Hz, microphone " << c.behind << " m behind");
        REQUIRE(std::abs(line_source_level(city, s, m, c.f, split_options()) - c.exact) < 0.8);
        REQUIRE(std::abs(line_source_level(city, s, m, c.f, unsplit) - c.exact) > 1.0);
    }
}

// The same step near the line through its far edges (boundary elements, 45
// elements a wavelength, within 0.15 dB of 30): sources 1 and 3.3 degrees
// above it, 60 and 20 m from the taller block's edge, microphones 2, 10 and
// 20 m past the lower block, 0.02 to 1 m either side, 125 and 250 Hz.  Split,
// within 1.2 dB, rms 0.47; the rays alone stepped there and were rms 1.21,
// up to 4.1 dB low just above the line.
TEST_CASE("by the line through a stepped pair the split field is the boundary-element solution's", "[point_paths]")
{
    const Buildings city({long_block(0.0, 10.0, 10.0), long_block(10.0, 20.0, 7.0)});
    static const std::vector<BlockCase> cases {
            {125, -47.1599, 28.2412, 2, 5.4, -7.22}, {125, -47.1599, 28.2412, 2, 6.1, -2.61}, {125, -47.1599, 28.2412, 2, 6.3, -2.85},
            {125, -47.1599, 28.2412, 2, 6.38, -3.03}, {125, -47.1599, 28.2412, 2, 6.42, -3.13}, {125, -47.1599, 28.2412, 2, 6.5, -3.29},
            {125, -47.1599, 28.2412, 2, 6.7, -3.35}, {125, -47.1599, 28.2412, 2, 7.4, 0.78}, {125, -47.1599, 28.2412, 10, 3, -5.84},
            {125, -47.1599, 28.2412, 10, 3.7, -7.43}, {125, -47.1599, 28.2412, 10, 3.9, -5.25}, {125, -47.1599, 28.2412, 10, 3.98, -4.46},
            {125, -47.1599, 28.2412, 10, 4.02, -4.09}, {125, -47.1599, 28.2412, 10, 4.1, -3.40}, {125, -47.1599, 28.2412, 10, 4.3, -2.06},
            {125, -47.1599, 28.2412, 10, 5, -0.77}, {125, -47.1599, 28.2412, 20, 0.7, -1.11}, {125, -47.1599, 28.2412, 20, 0.9, -2.36},
            {125, -47.1599, 28.2412, 20, 0.98, -2.99}, {125, -47.1599, 28.2412, 20, 1.02, -3.34}, {125, -47.1599, 28.2412, 20, 1.1, -4.10},
            {125, -47.1599, 28.2412, 20, 1.3, -6.46}, {125, -47.1599, 28.2412, 20, 2, -8.96}, {125, -46.3818, 30.5205, 2, 5.4, -7.71},
            {125, -46.3818, 30.5205, 2, 6.1, -3.18}, {125, -46.3818, 30.5205, 2, 6.3, -3.45}, {125, -46.3818, 30.5205, 2, 6.38, -3.64},
            {125, -46.3818, 30.5205, 2, 6.42, -3.74}, {125, -46.3818, 30.5205, 2, 6.5, -3.91}, {125, -46.3818, 30.5205, 2, 6.7, -3.98},
            {125, -46.3818, 30.5205, 2, 7.4, 0.10}, {125, -46.3818, 30.5205, 10, 3, -6.26}, {125, -46.3818, 30.5205, 10, 3.7, -7.34},
            {125, -46.3818, 30.5205, 10, 3.9, -5.20}, {125, -46.3818, 30.5205, 10, 3.98, -4.44}, {125, -46.3818, 30.5205, 10, 4.02, -4.09},
            {125, -46.3818, 30.5205, 10, 4.1, -3.44}, {125, -46.3818, 30.5205, 10, 4.3, -2.19}, {125, -46.3818, 30.5205, 10, 5, -1.21},
            {125, -46.3818, 30.5205, 20, 0.7, -1.14}, {125, -46.3818, 30.5205, 20, 0.9, -2.47}, {125, -46.3818, 30.5205, 20, 0.98, -3.15},
            {125, -46.3818, 30.5205, 20, 1.02, -3.52}, {125, -46.3818, 30.5205, 20, 1.1, -4.34}, {125, -46.3818, 30.5205, 20, 1.3, -6.90},
            {125, -46.3818, 30.5205, 20, 2, -8.20}, {125, -8.7939, 16.8402, 2, 5.4, -6.89}, {125, -8.7939, 16.8402, 2, 6.1, -2.34},
            {125, -8.7939, 16.8402, 2, 6.3, -2.62}, {125, -8.7939, 16.8402, 2, 6.38, -2.82}, {125, -8.7939, 16.8402, 2, 6.42, -2.92},
            {125, -8.7939, 16.8402, 2, 6.5, -3.10}, {125, -8.7939, 16.8402, 2, 6.7, -3.21}, {125, -8.7939, 16.8402, 2, 7.4, 0.83},
            {125, -8.7939, 16.8402, 10, 3, -5.34}, {125, -8.7939, 16.8402, 10, 3.7, -7.18}, {125, -8.7939, 16.8402, 10, 3.9, -4.91},
            {125, -8.7939, 16.8402, 10, 3.98, -4.09}, {125, -8.7939, 16.8402, 10, 4.02, -3.72}, {125, -8.7939, 16.8402, 10, 4.1, -3.02},
            {125, -8.7939, 16.8402, 10, 4.3, -1.67}, {125, -8.7939, 16.8402, 10, 5, -0.49}, {125, -8.7939, 16.8402, 20, 0.7, -0.51},
            {125, -8.7939, 16.8402, 20, 0.9, -1.76}, {125, -8.7939, 16.8402, 20, 0.98, -2.39}, {125, -8.7939, 16.8402, 20, 1.02, -2.73},
            {125, -8.7939, 16.8402, 20, 1.1, -3.50}, {125, -8.7939, 16.8402, 20, 1.3, -5.90}, {125, -8.7939, 16.8402, 20, 2, -9.43},
            {250, -47.1599, 28.2412, 2, 5.4, -10.03}, {250, -47.1599, 28.2412, 2, 6.1, -5.02}, {250, -47.1599, 28.2412, 2, 6.3, -5.92},
            {250, -47.1599, 28.2412, 2, 6.38, -6.09}, {250, -47.1599, 28.2412, 2, 6.42, -6.00}, {250, -47.1599, 28.2412, 2, 6.5, -5.43},
            {250, -47.1599, 28.2412, 2, 6.7, -2.96}, {250, -47.1599, 28.2412, 2, 7.4, -1.35}, {250, -47.1599, 28.2412, 10, 3.7, -4.17},
            {250, -47.1599, 28.2412, 10, 4.3, -6.37}, {250, -47.1599, 28.2412, 10, 5, -1.92}, {250, -47.1599, 28.2412, 20, 0.7, -9.01},
            {250, -47.1599, 28.2412, 20, 0.98, -14.52}, {250, -47.1599, 28.2412, 20, 1.02, -12.83}, {250, -47.1599, 28.2412, 20, 1.1, -9.71},
            {250, -47.1599, 28.2412, 20, 1.3, -4.65}, {250, -47.1599, 28.2412, 20, 2, -0.51}, {250, -46.3818, 30.5205, 2, 5.4, -10.40},
            {250, -46.3818, 30.5205, 2, 6.1, -5.37}, {250, -46.3818, 30.5205, 2, 6.3, -6.28}, {250, -46.3818, 30.5205, 2, 6.38, -6.42},
            {250, -46.3818, 30.5205, 2, 6.42, -6.31}, {250, -46.3818, 30.5205, 2, 6.5, -5.71}, {250, -46.3818, 30.5205, 2, 6.7, -3.25},
            {250, -46.3818, 30.5205, 2, 7.4, -2.13}, {250, -46.3818, 30.5205, 10, 3, -10.71}, {250, -46.3818, 30.5205, 10, 3.7, -4.10},
            {250, -46.3818, 30.5205, 10, 3.9, -5.38}, {250, -46.3818, 30.5205, 10, 3.98, -6.06}, {250, -46.3818, 30.5205, 10, 4.02, -6.37},
            {250, -46.3818, 30.5205, 10, 4.1, -6.75}, {250, -46.3818, 30.5205, 10, 4.3, -5.69}, {250, -46.3818, 30.5205, 10, 5, -1.92},
            {250, -46.3818, 30.5205, 20, 0.7, -8.83}, {250, -46.3818, 30.5205, 20, 0.9, -15.29}, {250, -46.3818, 30.5205, 20, 0.98, -12.75},
            {250, -46.3818, 30.5205, 20, 1.02, -11.05}, {250, -46.3818, 30.5205, 20, 1.1, -8.10}, {250, -46.3818, 30.5205, 20, 1.3, -3.42},
            {250, -46.3818, 30.5205, 20, 2, -0.04}, {250, -8.7939, 16.8402, 2, 5.4, -11.15}, {250, -8.7939, 16.8402, 2, 6.1, -6.15},
            {250, -8.7939, 16.8402, 2, 6.3, -7.11}, {250, -8.7939, 16.8402, 2, 6.38, -7.28}, {250, -8.7939, 16.8402, 2, 6.42, -7.19},
            {250, -8.7939, 16.8402, 2, 6.5, -6.61}, {250, -8.7939, 16.8402, 2, 6.7, -4.13}, {250, -8.7939, 16.8402, 2, 7.4, -2.84},
            {250, -8.7939, 16.8402, 10, 3, -11.85}, {250, -8.7939, 16.8402, 10, 3.7, -4.84}, {250, -8.7939, 16.8402, 10, 4.1, -7.86},
            {250, -8.7939, 16.8402, 10, 4.3, -7.01}, {250, -8.7939, 16.8402, 10, 5, -2.78}, {250, -8.7939, 16.8402, 20, 0.7, -9.38},
            {250, -8.7939, 16.8402, 20, 0.9, -17.08}, {250, -8.7939, 16.8402, 20, 0.98, -15.29}, {250, -8.7939, 16.8402, 20, 1.02, -13.45},
            {250, -8.7939, 16.8402, 20, 1.1, -10.16}, {250, -8.7939, 16.8402, 20, 1.3, -5.00}, {250, -8.7939, 16.8402, 20, 2, -0.92},
    };
    PointPathOptions unsplit;
    unsplit.tilted_step_split = false;
    double square = 0.0, unsplit_square = 0.0;
    for (const BlockCase& c : cases)
    {
        const Vec3 s(0.0, c.source_y, c.source_z), m(0.0, 20.0 + c.behind, c.z);
        const double error = line_source_level(city, s, m, c.f, split_options()) - c.exact;
        INFO(c.f << " Hz, source " << c.source_y << " m, " << c.source_z << " m up, microphone " << c.behind
                  << " m behind, " << c.z << " m up");
        REQUIRE(std::abs(error) < 1.5);
        square += error * error;
        const double unsplit_error = line_source_level(city, s, m, c.f, unsplit) - c.exact;
        unsplit_square += unsplit_error * unsplit_error;
    }
    const double n = static_cast<double>(cases.size());
    INFO("rms " << std::sqrt(square / n) << ", unsplit " << std::sqrt(unsplit_square / n));
    REQUIRE(std::sqrt(square / n) < 0.6);
    REQUIRE(std::sqrt(unsplit_square / n) > 1.0);
}

// Further from the line (boundary elements, 45 elements a wavelength,
// within 0.15 dB of 30): sources 0.5, 2, 5 and 10 degrees above it, 60 m
// from the taller block's edge, microphones 5 and 10 m past the lower block,
// 2.4 m below the line to 3.4 m above, 125 Hz to 1 kHz; 68 of 352 cases.
// Split, rms 0.27 (0.26 over all 352, within 0.9 dB); the rays alone 0.57
// (0.67, up to 4.2 dB low).
TEST_CASE("past a stepped pair the split field is the boundary-element solution's", "[point_paths]")
{
    const Buildings city({long_block(0.0, 10.0, 10.0), long_block(10.0, 20.0, 7.0)});
    static const std::vector<BlockCase> cases {
            {125, -47.3169, 27.7417, 5, 3, -7.65}, {125, -47.3169, 27.7417, 5, 5, -3.50}, {125, -47.3169, 27.7417, 5, 7.5, 1.73},
            {125, -47.3169, 27.7417, 10, 4.5, -1.28}, {125, -47.3169, 27.7417, 10, 5, -0.81}, {125, -47.3169, 27.7417, 10, 7.5, 1.98},
            {125, -46.8329, 29.2360, 5, 4.5, -8.81}, {125, -46.8329, 29.2360, 5, 6.5, -2.91}, {125, -46.8329, 29.2360, 5, 7.5, 1.33},
            {125, -46.8329, 29.2360, 10, 4.5, -1.27}, {125, -46.8329, 29.2360, 10, 6.5, -0.56}, {125, -45.7482, 32.1841, 5, 4.1, -8.80},
            {125, -45.7482, 32.1841, 10, 6.5, 0.25}, {125, -43.6026, 36.9584, 5, 4.1, -8.33}, {125, -43.6026, 36.9584, 5, 6, -3.36},
            {125, -43.6026, 36.9584, 10, 3.9, -2.83}, {125, -43.6026, 36.9584, 10, 4.1, -1.31}, {125, -43.6026, 36.9584, 10, 6, -2.17},
            {250, -47.3169, 27.7417, 5, 3.9, -14.39}, {250, -47.3169, 27.7417, 5, 5.6, -4.52}, {250, -47.3169, 27.7417, 5, 6, -4.85},
            {250, -47.3169, 27.7417, 10, 3.9, -5.65}, {250, -47.3169, 27.7417, 10, 5.6, -2.99}, {250, -46.8329, 29.2360, 5, 3.5, -8.46},
            {250, -46.8329, 29.2360, 10, 5.6, -3.06}, {250, -45.7482, 32.1841, 5, 3.5, -7.58}, {250, -45.7482, 32.1841, 5, 5.4, -3.26},
            {250, -45.7482, 32.1841, 10, 3, -9.60}, {250, -45.7482, 32.1841, 10, 3.5, -3.18}, {250, -45.7482, 32.1841, 10, 5.4, -2.91},
            {250, -43.6026, 36.9584, 5, 3, -14.41}, {250, -43.6026, 36.9584, 5, 5, -4.53}, {250, -43.6026, 36.9584, 5, 5.4, -1.15},
            {250, -43.6026, 36.9584, 10, 3, -6.06}, {250, -43.6026, 36.9584, 10, 5, 0.35}, {250, -43.6026, 36.9584, 10, 7.5, -0.11},
            {500, -47.3169, 27.7417, 10, 5, -3.15}, {500, -47.3169, 27.7417, 10, 7.5, 1.45}, {500, -46.8329, 29.2360, 5, 4.5, -10.11},
            {500, -46.8329, 29.2360, 5, 6.5, -1.02}, {500, -46.8329, 29.2360, 5, 7.5, 0.14}, {500, -46.8329, 29.2360, 10, 4.5, -4.66},
            {500, -46.8329, 29.2360, 10, 6.5, 0.10}, {500, -45.7482, 32.1841, 5, 4.1, -11.06}, {500, -45.7482, 32.1841, 5, 4.5, -7.18},
            {500, -45.7482, 32.1841, 5, 6.5, 0.96}, {500, -45.7482, 32.1841, 10, 4.1, -2.04}, {500, -45.7482, 32.1841, 10, 6, 1.17},
            {500, -43.6026, 36.9584, 10, 4.1, 1.38}, {500, -43.6026, 36.9584, 10, 6, 0.35}, {1000, -47.3169, 27.7417, 5, 3.9, -15.89},
            {1000, -47.3169, 27.7417, 5, 5.6, -7.29}, {1000, -47.3169, 27.7417, 5, 6, -3.10}, {1000, -47.3169, 27.7417, 10, 3.9, -7.77},
            {1000, -47.3169, 27.7417, 10, 5.6, 0.45}, {1000, -46.8329, 29.2360, 5, 3.5, -15.09}, {1000, -46.8329, 29.2360, 5, 3.9, -14.61},
            {1000, -46.8329, 29.2360, 5, 5.6, -5.68}, {1000, -46.8329, 29.2360, 10, 3.5, -9.35}, {1000, -46.8329, 29.2360, 10, 5.4, 0.55},
            {1000, -45.7482, 32.1841, 10, 3.5, -5.11}, {1000, -45.7482, 32.1841, 10, 5.4, 0.41}, {1000, -43.6026, 36.9584, 5, 3, -9.91},
            {1000, -43.6026, 36.9584, 5, 5, -2.96}, {1000, -43.6026, 36.9584, 5, 5.4, 0.36}, {1000, -43.6026, 36.9584, 10, 3, -0.61},
            {1000, -43.6026, 36.9584, 10, 5, -1.12}, {1000, -43.6026, 36.9584, 10, 7.5, 0.31},
    };
    PointPathOptions unsplit;
    unsplit.tilted_step_split = false;
    double square = 0.0, unsplit_square = 0.0;
    for (const BlockCase& c : cases)
    {
        const Vec3 s(0.0, c.source_y, c.source_z), m(0.0, 20.0 + c.behind, c.z);
        const double error = line_source_level(city, s, m, c.f, split_options()) - c.exact;
        INFO(c.f << " Hz, source " << c.source_y << " m, " << c.source_z << " m up, microphone " << c.behind
                  << " m behind, " << c.z << " m up");
        REQUIRE(std::abs(error) < 1.0);
        square += error * error;
        const double unsplit_error = line_source_level(city, s, m, c.f, unsplit) - c.exact;
        unsplit_square += unsplit_error * unsplit_error;
    }
    const double n = static_cast<double>(cases.size());
    INFO("rms " << std::sqrt(square / n) << ", unsplit " << std::sqrt(unsplit_square / n));
    REQUIRE(std::sqrt(square / n) < 0.35);
    REQUIRE(std::sqrt(unsplit_square / n) > std::sqrt(square / n) + 0.2);
}

// The aperture's Filon rule against the midpoint rule at a 192nd of a
// wavelength, on shared faces from DigiCity and a district: both ends below
// the face's plane, the receiver above it, a source in the plane (the reverse
// aperture's receiver at theta = pi, both boundaries at the edge), a source
// above it.  Over 4000 such faces it is within 0.016 dB and 1.5e-3 rad; a node
// on a cut, where the wedge's coefficient jumps, was 6-28 dB off.
TEST_CASE("the aperture's Filon rule is the fine midpoint rule's", "[point_paths]")
{
    struct Face
    {
        double width, sx, sz, rx, rz, n;
    };
    // The last two are obtuse corners (n = 1.25), where the wave reflected
    // off the second wedge's other face reaches the receiver: right-angled
    // ones never take that path.
    const std::vector<Face> faces {
        {2.438, -3.9627, -3.2527, 14.3812, -29.5687, 1.5},
        {2.438, -4.3169, -21.6163, 6.1392, 6.4313, 1.5},
        {12.116, -4.9019, -20.0857, 19.3141, 0.9958, 1.5},
        {60.0, -90.0, 0.0, 90.0, -28.5, 1.5},
        {60.0, -90.0, 10.0, 70.0, -90.0, 1.5},
        {6.0, -10.0, -5.0, 16.0, -3.64, 1.25},
        {4.0, -12.0, -4.0, 9.0, -2.5, 1.25},
    };
    double worst = 0.0, worst_phase = 0.0;
    for (const Face& f : faces)
        for (double frequency : {125.0, 1000.0, 5000.0})
        {
            const SharedFace face {f.width, f.n, f.n};
            const double k = 2.0 * pi * frequency / 343.2;
            const std::complex<double> filon = shared_face_field(face, f.sx, f.sz, f.rx, f.rz, k);
            const std::complex<double> midpoint = midpoint_shared_face(face, f.sx, f.sz, f.rx, f.rz, k);
            worst = std::max(worst, std::abs(20.0 * std::log10(std::abs(filon) / std::abs(midpoint))));
            worst_phase = std::max(worst_phase, std::abs(std::arg(filon / midpoint)));
        }
    INFO("worst " << worst << " dB, " << worst_phase << " rad");
    REQUIRE(worst < 0.03);
    REQUIRE(worst_phase < 0.005);
}

// The aperture worked out once for every wavenumber from the lowest band's
// (shared_face_geometry), its nodes no longer ending at each k's taper and
// extent, still gives each k's field as the fine midpoint rule does, on the
// Filon test's faces from 50 Hz to 10 kHz.
TEST_CASE("the aperture worked out once is every wavenumber's", "[point_paths]")
{
    struct Face
    {
        double width, sx, sz, rx, rz, n;
    };
    const std::vector<Face> faces {
        {2.438, -3.9627, -3.2527, 14.3812, -29.5687, 1.5},
        {2.438, -4.3169, -21.6163, 6.1392, 6.4313, 1.5},
        {12.116, -4.9019, -20.0857, 19.3141, 0.9958, 1.5},
        {60.0, -90.0, 0.0, 90.0, -28.5, 1.5},
        {60.0, -90.0, 10.0, 70.0, -90.0, 1.5},
        {6.0, -10.0, -5.0, 16.0, -3.64, 1.25},
        {4.0, -12.0, -4.0, 9.0, -2.5, 1.25},
    };
    const double c = 343.2;
    const std::vector<double> frequencies {50.0, 125.0, 400.0, 1000.0, 2500.0, 5000.0, 10000.0};
    std::vector<double> wavenumbers;
    for (double frequency : frequencies)
        wavenumbers.push_back(2.0 * pi * frequency / c);
    double worst = 0.0, worst_phase = 0.0, worst_cpu = 0.0;
    for (const Face& f : faces)
    {
        const SharedFace face {f.width, f.n, f.n};
        const SharedFaceGeometry geometry = shared_face_geometry(face, f.sx, f.sz, f.rx, f.rz, wavenumbers);
        const double reference = std::hypot(f.sx, f.sz) + f.width + std::hypot(f.rx - f.width, f.rz);
        for (double frequency : frequencies)
        {
            const double k = 2.0 * pi * frequency / c;
            const std::complex<double> once = shared_face_field_at(geometry, k) * std::polar(1.0, k * reference);
            const std::complex<double> midpoint = midpoint_shared_face(face, f.sx, f.sz, f.rx, f.rz, k);
            const std::complex<double> cpu = shared_face_field(face, f.sx, f.sz, f.rx, f.rz, k);
            worst = std::max(worst, std::abs(20.0 * std::log10(std::abs(once) / std::abs(midpoint))));
            worst_phase = std::max(worst_phase, std::abs(std::arg(once / midpoint)));
            worst_cpu = std::max(worst_cpu, std::abs(20.0 * std::log10(std::abs(once) / std::abs(cpu))));
        }
    }
    INFO("worst " << worst << " dB, " << worst_phase << " rad against the midpoint rule; " << worst_cpu
                  << " dB against shared_face_field");
    REQUIRE(worst < 0.03);
    REQUIRE(worst_phase < 0.005);
}

// Sampled (PointPathOptions::sample_shared_faces), a shared face's amplitude is
// the cascade's times the integral's ratio to it, taken two thirds of an octave apart while
// the first wedge's transition parameter is under 30: against the integral at
// every one of 84 sub-bands from 20 Hz to 10 kHz, over long blocks 1 to 10 m
// thick with sources low, grazing, in the roof's plane and just above it and
// microphones behind at heights through its plane.  There the cascade alone
// is up to 16 dB off, and of the wrong sign with an end in the plane.
TEST_CASE("a sampled shared face is its integral at every wavenumber", "[point_paths]")
{
    std::vector<double> bands;
    for (int b = 13; b <= 40; b++)
        bands.push_back(std::pow(10.0, b / 10.0));
    const double c = 343.2;
    double worst = 0.0, worst_phase = 0.0;
    std::size_t shared = 0;
    for (double width : {1.0, 2.44, 5.0, 10.0})
    {
        Footprint block;
        block.x = {-2000.0, 2000.0, 2000.0, -2000.0};
        block.y = {0.0, 0.0, width, width};
        block.top = 5.18;
        const Buildings city({block});
        const std::vector<double> reflection(city.facades().size() * bands.size(), 1.0);
        for (const Vec3& s : {Vec3(0.0, -25.0, 4.5), Vec3(0.0, -10.0, 2.0), Vec3(0.0, -200.0, 4.0),
                     Vec3(0.0, -15.0, 6.5), Vec3(30.0, -60.0, 5.1), Vec3(0.0, -140.0, 5.18), Vec3(10.0, -80.0, 5.3)})
            for (double behind : {2.0, 25.0})
                for (double z : {1.2, 5.26, 6.5})
                {
                    const Vec3 m(0.0, width + behind, z);
                    PointPathOptions exact;
                    exact.sample_shared_faces = false;
                    const PointSourceField sampled_field(city, s, m), exact_field(city, s, m, exact);
                    std::vector<std::vector<std::complex<double>>> amplitude[2];
                    std::vector<std::size_t> order;
                    for (int which = 0; which < 2; which++)
                        (which ? exact_field : sampled_field)
                                .sub_band_problems(
                                        reflection, 1.0, bands, c,
                                        [&](std::size_t, int, const std::vector<std::complex<double>>& a,
                                                const std::vector<double>&, const SubBand&) { amplitude[which].push_back(a); },
                                        {}, nullptr, &order);
                    REQUIRE(amplitude[0].size() == amplitude[1].size());
                    for (std::size_t r = 0; r < order.size(); r++)
                    {
                        if (!sampled_field.rays()[order[r]].shared)
                            continue;
                        shared++;
                        for (std::size_t q = 0; q < amplitude[0].size(); q++)
                        {
                            const std::complex<double> ratio = amplitude[0][q][r] / amplitude[1][q][r];
                            worst = std::max(worst, std::abs(20.0 * std::log10(std::abs(ratio))));
                            worst_phase = std::max(worst_phase, std::abs(std::arg(ratio)));
                        }
                    }
                }
    }
    INFO(shared << " rays over shared faces; worst " << worst << " dB, " << worst_phase << " rad");
    REQUIRE(shared > 50);
    REQUIRE(worst < 0.1);
    REQUIRE(worst_phase < 0.03);
}

// shared_face_spectrum shares its nodes' geometry between wavenumbers but not
// its nodes: at every one of a run's 72 sub-bands (24 bands, 50 Hz to 10 kHz)
// it is shared_face_field to the bit, on the Filon test's faces -- which take
// the forward aperture alone, the reverse alone, and both.
TEST_CASE("the shared face's spectrum is its field at each wavenumber", "[point_paths]")
{
    struct Face
    {
        double width, sx, sz, rx, rz, n;
    };
    const std::vector<Face> faces {
        {2.438, -3.9627, -3.2527, 14.3812, -29.5687, 1.5},
        {2.438, -4.3169, -21.6163, 6.1392, 6.4313, 1.5},
        {12.116, -4.9019, -20.0857, 19.3141, 0.9958, 1.5},
        {60.0, -90.0, 0.0, 90.0, -28.5, 1.5},
        {60.0, -90.0, 10.0, 70.0, -90.0, 1.5},
        {6.0, -10.0, -5.0, 16.0, -3.64, 1.25},
        {4.0, -12.0, -4.0, 9.0, -2.5, 1.25},
    };
    std::vector<double> wavenumbers;
    for (int b = 17; b <= 40; b++)
        for (const SubBand& part : sub_bands(std::pow(10.0, b / 10.0), 343.2))
            wavenumbers.push_back(0.5 * (part.k1 + part.k2));
    std::size_t compared = 0;
    for (const Face& f : faces)
    {
        const SharedFace face {f.width, f.n, f.n};
        const std::vector<std::complex<double>> spectrum
                = shared_face_spectrum(face, f.sx, f.sz, f.rx, f.rz, wavenumbers);
        REQUIRE(spectrum.size() == wavenumbers.size());
        for (std::size_t j = 0; j < wavenumbers.size(); j++)
        {
            const std::complex<double> field = shared_face_field(face, f.sx, f.sz, f.rx, f.rz, wavenumbers[j]);
            INFO("width " << f.width << ", k " << wavenumbers[j]);
            REQUIRE(spectrum[j].real() == field.real());
            REQUIRE(spectrum[j].imag() == field.imag());
            compared++;
        }
    }
    REQUIRE(compared == 7 * 72);
}

// An end within a hair of the line between the edges, the other above it:
// the aperture's end-point terms, each as 1/g, once failed to cancel there,
// so the field strayed from 3e-5 rad and grew without bound below 1e-6 (to
// 13 dB over the direct wave on a suburban street, where a microphone lay in a
// slope's plane).  It now holds the limit it reaches by 1e-3 rad, through the
// field, the spectrum and the sampled geometry alike.
TEST_CASE("an end at the line between the edges keeps the shared face's limit", "[point_paths]")
{
    struct Face
    {
        double width, d, r, a, n;
    };
    const std::vector<Face> faces {{15.7, 31.4, 100.0, 0.5, 1.5}, {10.0, 20.0, 30.0, 0.05, 1.5},
            {5.0, 50.0, 10.0, 0.3, 1.5}, {15.7, 31.4, 100.0, 0.5, 1.9}};
    const std::vector<double> wavenumbers {2.3, 9.2, 36.6};
    for (const Face& f : faces)
    {
        const SharedFace face {f.width, f.n, f.n};
        const double rx = f.width + f.r * std::cos(f.a), rz = f.r * std::sin(f.a);
        auto level = [&](double g, double k) {
            return 20.0 * std::log10(std::abs(shared_face_field(face, -f.d * std::cos(g), -f.d * std::sin(g), rx, rz, k)));
        };
        for (const double g : {1e-4, 3e-5, 1e-6, 1e-8, 1e-12})
        {
            const double sx = -f.d * std::cos(g), sz = -f.d * std::sin(g);
            const auto spectrum = shared_face_spectrum(face, sx, sz, rx, rz, wavenumbers);
            const SharedFaceGeometry geometry = shared_face_geometry(face, sx, sz, rx, rz, wavenumbers);
            for (std::size_t j = 0; j < wavenumbers.size(); j++)
            {
                const double k = wavenumbers[j];
                INFO("width " << f.width << ", g " << g << ", k " << k);
                CHECK(std::abs(level(g, k) - level(3e-3, k)) < 0.3);
                CHECK(std::abs(20.0 * std::log10(std::abs(spectrum[j])) - level(g, k)) < 1e-9);
                CHECK(std::abs(20.0 * std::log10(std::abs(shared_face_field_at(geometry, k))) - level(g, k)) < 0.05);
            }
        }
    }
}

// The diffracted paths ruled out before they are solved or checked are only
// those that could not have been rays: among DigiCity's containers, a block
// of a district, and a grid of blocks of mixed heights over several of
// Buildings' grid cells, with the source low and high, the rays are the same
// with the pruning off, to the bit.
TEST_CASE("pruning the diffracted paths loses no ray", "[point_paths]")
{
    const Buildings containers({box(4.09, 17.66, 6.53, 23.72, 7.77), box(14.81, 17.66, 17.25, 23.72, 7.77),
            box(1.6, 10.82, 4.04, 16.88, 7.77), box(1.6, 8.38, 13.72, 10.82, 5.18),
            box(20.21, 3.86, 22.65, 16.0, 7.77), box(13.5, -1.84, 19.56, 0.6, 7.77)});
    const Buildings blocks({box(0.0, 0.0, 60.0, 60.0, 20.0), box(80.0, 0.0, 140.0, 60.0, 35.0),
            box(0.0, 80.0, 60.0, 140.0, 10.0), box(80.0, 80.0, 140.0, 140.0, 40.0)});
    std::vector<Footprint> grid_blocks;
    const double heights[] = {6.0, 15.0, 9.0, 24.0, 12.0, 30.0, 9.0};
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            grid_blocks.push_back(box(36.0 * i, 36.0 * j, 36.0 * i + 22.0 + 2.0 * (j % 3), 36.0 * j + 24.0 - 3.0 * (i % 2),
                    heights[(3 * i + j) % 7]));
    const Buildings grid(grid_blocks);
    struct Case
    {
        const Buildings* city;
        Vec3 source, microphone;
    };
    const std::vector<Case> cases {
        {&containers, {10.5, 20.5, 2.0}, {10.5, -1.1, 1.2}},
        {&containers, {10.5, 20.5, 2.0}, {-3.6, 6.4, 1.2}},
        {&containers, {8.0, 12.0, 15.0}, {18.0, 2.0, 1.2}},
        {&blocks, {70.0, 70.0, 30.0}, {30.0, 70.0, 1.5}},
        {&blocks, {70.0, 30.0, 5.0}, {110.0, 150.0, 1.5}},
        {&grid, {101.0, 31.0, 40.0}, {29.0, 140.0, 1.5}},
        {&grid, {66.0, 101.0, 2.0}, {140.0, 29.0, 1.2}},
    };
    for (const Case& c : cases)
    {
        PointPathOptions all;
        all.prune = false;
        const PointSourceField pruned(*c.city, c.source, c.microphone);
        const PointSourceField every(*c.city, c.source, c.microphone, all);
        REQUIRE(pruned.rays().size() == every.rays().size());
        std::size_t twice = 0;
        for (std::size_t i = 0; i < every.rays().size(); i++)
        {
            REQUIRE(pruned.rays()[i].length == every.rays()[i].length);
            REQUIRE(pruned.rays()[i].wedges == every.rays()[i].wedges);
            twice += every.rays()[i].wedges == 2;
        }
        REQUIRE(twice > 0);
    }
}

// A run over many sources traces each end's image tree once, its images
// keeping the edges near any end it will meet, and each pair keeps those near
// it (PathEnd): the field must be the one traced for the pair alone, ray for
// ray.  Leaving the partners' edges out of an end's candidates gives
// different rays, and is refused.
TEST_CASE("a field from ends traced once is the pair's own", "[point_paths]")
{
    const Buildings blocks({box(0.0, 0.0, 60.0, 60.0, 20.0), box(80.0, 0.0, 140.0, 60.0, 35.0),
            box(0.0, 80.0, 60.0, 140.0, 10.0), box(80.0, 80.0, 140.0, 140.0, 40.0),
            box(400.0, 400.0, 460.0, 460.0, 25.0)});
    PointPathOptions options;
    options.radius = 120.0;
    // The last source and microphone are far apart: each pair keeps a
    // different set of edges.
    const std::vector<Vec3> sources {{70.0, 70.0, 30.0}, {70.0, 30.0, 5.0}, {300.0, 300.0, 40.0}};
    const std::vector<Vec3> microphones {{30.0, 70.0, 1.5}, {110.0, 150.0, 1.5}, {350.0, 380.0, 1.2}};
    const std::vector<std::size_t> near_sources = edges_near(blocks, sources, options.radius),
                                   near_microphones = edges_near(blocks, microphones, options.radius);
    std::size_t compared = 0;
    for (const Vec3& s : sources)
    {
        const PathEnd source(blocks, s, true, options, near_microphones);
        for (const Vec3& m : microphones)
        {
            const PathEnd microphone(blocks, m, false, options, near_sources);
            const PointSourceField alone(blocks, s, m, options);
            const PointSourceField shared(blocks, source, microphone, options);
            REQUIRE(alone.rays().size() == shared.rays().size());
            for (std::size_t i = 0; i < alone.rays().size(); i++)
            {
                REQUIRE(alone.rays()[i].length == shared.rays()[i].length);
                REQUIRE(alone.rays()[i].spreading == shared.rays()[i].spreading);
                REQUIRE(alone.rays()[i].facades == shared.rays()[i].facades);
            }
            compared += alone.rays().size();
        }
    }
    REQUIRE(compared > 1000);
    // An end traced without the edges near its partner is refused.
    const PathEnd lonely(blocks, sources[2], true, options, {});
    const PathEnd far_microphone(blocks, microphones[0], false, options, near_sources);
    REQUIRE_THROWS_AS(PointSourceField(blocks, lonely, far_microphone, options), std::invalid_argument);
    // So are ends in the wrong roles, ends traced with other options, and
    // partners out of order.
    const PathEnd source(blocks, sources[0], true, options, near_microphones);
    const PathEnd microphone(blocks, microphones[0], false, options, near_sources);
    REQUIRE_THROWS_AS(PointSourceField(blocks, microphone, source, options), std::invalid_argument);
    PointPathOptions other = options;
    other.max_edge_reflections = 1;
    REQUIRE_THROWS_AS(PointSourceField(blocks, source, microphone, other), std::invalid_argument);
    // Without roofs an end finds no slopes (PathEnd::slopes).
    PointPathOptions flat = options;
    flat.roofs = false;
    REQUIRE_THROWS_AS(PointSourceField(blocks, source, microphone, flat), std::invalid_argument);
    std::vector<std::size_t> reversed(near_sources.rbegin(), near_sources.rend());
    REQUIRE_THROWS_AS(PathEnd(blocks, microphones[0], false, options, reversed), std::invalid_argument);
}

// The twice-diffracted search takes the microphone's images to one
// reflection whatever max_edge_reflections is: at 0 the rays diffracted
// twice are the default's, with the source's side the same in both.
TEST_CASE("twice-diffracted rays reflect at the microphone whatever the edge order", "[point_paths]")
{
    const Buildings containers({box(4.09, 17.66, 6.53, 23.72, 7.77), box(14.81, 17.66, 17.25, 23.72, 7.77),
            box(1.6, 10.82, 4.04, 16.88, 7.77), box(1.6, 8.38, 13.72, 10.82, 5.18),
            box(20.21, 3.86, 22.65, 16.0, 7.77), box(13.5, -1.84, 19.56, 0.6, 7.77)});
    PointPathOptions deep;
    deep.source_reflections = false;
    PointPathOptions none = deep;
    none.max_edge_reflections = 0;
    std::size_t reflected = 0;
    for (const Vec3& microphone : {Vec3(10.5, -1.1, 1.2), Vec3(-3.6, 6.4, 1.2), Vec3(18.0, 2.0, 1.2)})
    {
        const Vec3 source(10.5, 20.5, 2.0);
        const PointSourceField a(containers, source, microphone, deep), b(containers, source, microphone, none);
        std::vector<const PointRay*> twice_a, twice_b;
        for (const PointRay& ray : a.rays())
            if (ray.wedges == 2)
                twice_a.push_back(&ray);
        for (const PointRay& ray : b.rays())
            if (ray.wedges == 2)
                twice_b.push_back(&ray);
        REQUIRE(twice_a.size() == twice_b.size());
        for (std::size_t i = 0; i < twice_a.size(); i++)
        {
            REQUIRE(twice_a[i]->length == twice_b[i]->length);
            REQUIRE(twice_a[i]->facades == twice_b[i]->facades);
            reflected += !twice_b[i]->facades.empty();
        }
    }
    REQUIRE(reflected > 0);
}

// A microphone in a street behind a row, which the source reaches over the
// row's roof: what reaches it reflects between the street's walls after the
// second edge.  Each end's images go to max_twice_reflections, so the rays at
// one are the rays at two that reflect once at most, to the bit, and two adds
// rays that reflect twice.  Swapping the ends gives the same rays.
TEST_CASE("twice-diffracted rays reflect to max_twice_reflections at either end", "[point_paths]")
{
    const Buildings rows({box(0.0, -20.0, 10.0, 20.0, 8.0), box(16.0, -20.0, 26.0, 20.0, 6.0)});
    const Vec3 source(-15.0, 0.0, 1.5), microphone(13.0, 0.0, 1.2);
    PointPathOptions once;
    once.max_twice_reflections = 1;
    const PointPathOptions twice;
    REQUIRE(twice.max_twice_reflections == 2);
    auto lengths = [](const PointSourceField& field, std::size_t& deep) {
        std::vector<double> out;
        deep = 0;
        for (const PointRay& ray : field.rays())
            if (ray.wedges == 2)
            {
                out.push_back(ray.length);
                deep += ray.facades.size() >= 2;
            }
        std::sort(out.begin(), out.end());
        return out;
    };
    std::size_t deep_once = 0, deep_twice = 0, deep_swapped = 0;
    const std::vector<double> a = lengths(PointSourceField(rows, source, microphone, once), deep_once);
    const std::vector<double> b = lengths(PointSourceField(rows, source, microphone, twice), deep_twice);
    REQUIRE(!a.empty());
    REQUIRE(b.size() > a.size());
    REQUIRE(std::includes(b.begin(), b.end(), a.begin(), a.end()));
    REQUIRE(deep_twice > 0);
    const std::vector<double> c = lengths(PointSourceField(rows, microphone, source, twice), deep_swapped);
    REQUIRE(c.size() == b.size());
    for (std::size_t i = 0; i < b.size(); i++)
        REQUIRE(c[i] == Catch::Approx(b[i]).epsilon(1e-12));
}


// The sorted lengths of a field's rays with `wedges` edges, and the count of
// those that reflect between two edges.
std::vector<double> lengths_of(const PointSourceField& field, int wedges, std::size_t* between = nullptr)
{
    std::vector<double> out;
    std::size_t count = 0;
    for (const PointRay& ray : field.rays())
        if (ray.wedges == wedges)
        {
            out.push_back(ray.length);
            count += ray.between != PointRay::no_facade;
        }
    std::sort(out.begin(), out.end());
    if (between != nullptr)
        *between = count;
    return out;
}

// A notch: a low house between two tall ones, the source in front of it and
// the microphone behind, on its middle line.  Over the low roof the ray may
// reflect in either tall house's wall between the roof's two edges: unfolded
// about the wall it is the double diffraction to the microphone's mirror, the
// two edges parallel, so its length is the transverse path over the roof with
// the mirror's offset along the edges added in quadrature.  Swapping the ends
// gives the same rays, and so does tracing them with pruning off.
TEST_CASE("a ray reflects between two edges as the double diffraction to the mirror", "[point_paths]")
{
    const Buildings notch({box(0.0, 10.0, 10.0, 20.0, 12.0), box(0.0, 0.0, 10.0, 10.0, 8.0),
            box(0.0, -10.0, 10.0, 0.0, 12.0)});
    const Vec3 source(-15.0, 5.0, 1.5), microphone(25.0, 5.0, 1.2);
    PointPathOptions straight;
    straight.max_between_reflections = 0;
    const PointPathOptions between;
    REQUIRE(between.max_between_reflections == 1);
    std::size_t none = 0, some = 0, swapped = 0, unpruned = 0;
    const PointSourceField field(notch, source, microphone, between);
    const std::vector<double> plain = lengths_of(PointSourceField(notch, source, microphone, straight), 2, &none);
    const std::vector<double> with = lengths_of(field, 2, &some);
    REQUIRE(none == 0);
    REQUIRE(some > 0);
    REQUIRE(std::includes(with.begin(), with.end(), plain.begin(), plain.end()));
    const double transverse = std::hypot(15.0, 6.5) + 10.0 + std::hypot(15.0, 6.8);
    const double expected = std::hypot(transverse, 10.0);
    bool found = false;
    for (const PointRay& ray : field.rays())
        if (ray.between != PointRay::no_facade && ray.floor == 0 && ray.facades.size() == 1)
        {
            found = found || std::abs(ray.length - expected) < 1e-9;
            // It reflects in a tall house's wall, at the roof's height.
            const std::size_t house = notch.facade_building(ray.between);
            REQUIRE((house == 0 || house == 2));
            REQUIRE(ray.path.size() == 5);
            REQUIRE(ray.path[2].z == Catch::Approx(8.0).margin(1e-9));
        }
    REQUIRE(found);
    const std::vector<double> back = lengths_of(PointSourceField(notch, microphone, source, between), 2, &swapped);
    REQUIRE(swapped == some);
    REQUIRE(back.size() == with.size());
    for (std::size_t i = 0; i < with.size(); i++)
        REQUIRE(back[i] == Catch::Approx(with[i]).epsilon(1e-12));
    PointPathOptions every = between;
    every.prune = false;
    REQUIRE(lengths_of(PointSourceField(notch, source, microphone, every), 2, &unpruned) == with);
    REQUIRE(unpruned == some);
}

// Among DigiCity-like containers some twice-diffracted rays reflect in the
// floor between their edges: each bounces on the floor in the air, its
// length is its real legs' sum, and the rays are the same with the ends
// swapped and with pruning off.
TEST_CASE("a ray reflects in the floor between two edges", "[point_paths]")
{
    const Buildings containers({box(4.09, 17.66, 6.53, 23.72, 7.77), box(14.81, 17.66, 17.25, 23.72, 7.77),
            box(1.6, 10.82, 4.04, 16.88, 7.77), box(1.6, 8.38, 13.72, 10.82, 5.18),
            box(20.21, 3.86, 22.65, 16.0, 7.77), box(13.5, -1.84, 19.56, 0.6, 7.77)});
    std::size_t bounced = 0, checked = 0;
    for (const Vec3& microphone : {Vec3(10.5, -1.1, 1.2), Vec3(-3.6, 6.4, 1.2), Vec3(18.0, 2.0, 1.2)})
    {
        const Vec3 source(10.5, 20.5, 2.0);
        const PointSourceField field(containers, source, microphone);
        for (const PointRay& ray : field.rays())
        {
            if (!ray.between_floor)
                continue;
            bounced++;
            REQUIRE(ray.wedges == 2);
            REQUIRE(ray.floor >= 1);
            // A ray that also bounces on the floor at an end draws that leg
            // unfolded below it; with the one bounce between its edges, the
            // path is the real one.
            if (ray.floor == 1 && ray.facades.empty())
            {
                double length = 0.0, lowest = 1e300;
                for (std::size_t k = 1; k < ray.path.size(); k++)
                    lowest = std::min(lowest, ray.path[k].z);
                REQUIRE(lowest == Catch::Approx(0.0).margin(1e-9));
                checked++;
                for (std::size_t k = 1; k < ray.path.size(); k++)
                    length += (ray.path[k] - ray.path[k - 1]).norm();
                REQUIRE(length == Catch::Approx(ray.length).epsilon(1e-9));
            }
        }
        const std::vector<double> there = lengths_of(field, 2);
        const std::vector<double> back = lengths_of(PointSourceField(containers, microphone, source), 2);
        REQUIRE(back.size() == there.size());
        for (std::size_t i = 0; i < there.size(); i++)
            REQUIRE(back[i] == Catch::Approx(there[i]).epsilon(1e-12));
        PointPathOptions every;
        every.prune = false;
        REQUIRE(lengths_of(PointSourceField(containers, source, microphone, every), 2) == there);
    }
    REQUIRE(bounced > 0);
    REQUIRE(checked > 0);
}

// Two rows, the first the taller, the microphone behind the second: the path
// crosses the first row's two roof edges and the second's far one, the second's
// near edge lit below it.  The three edges are parallel, so the ray's length is
// the transverse path and its spreading 1 / sqrt(s1 s2 s3 s4 (s1 + s2 + s3 + s4)).
// None without three edges; the same rays with the ends swapped and with
// pruning off.
TEST_CASE("a ray crosses three roof edges as the unfolded path does", "[point_paths]")
{
    const Buildings rows({box(0.0, -20.0, 10.0, 20.0, 9.0), box(16.0, -20.0, 26.0, 20.0, 6.0)});
    const Vec3 source(-15.0, 0.0, 1.5), microphone(31.0, 0.0, 1.2);
    const double s1 = std::hypot(15.0, 7.5), s2 = 10.0, s3 = std::hypot(16.0, 3.0), s4 = std::hypot(5.0, 4.8);
    const double total = s1 + s2 + s3 + s4;
    const PointPathOptions options;
    REQUIRE(options.triple_diffraction == PointPathOptions::Triple::roofs);
    const PointSourceField field(rows, source, microphone, options);
    bool found = false;
    for (const PointRay& ray : field.rays())
        if (ray.wedges == 3 && ray.floor == 0 && ray.facades.empty() && std::abs(ray.length - total) < 1e-9)
        {
            found = true;
            REQUIRE(ray.spreading == Catch::Approx(1.0 / std::sqrt(s1 * s2 * s3 * s4 * total)).epsilon(1e-12));
            REQUIRE(ray.path.size() == 5);
            REQUIRE(ray.path[3].x == Catch::Approx(26.0).margin(1e-9));
        }
    REQUIRE(found);
    PointPathOptions without = options;
    without.triple_diffraction = PointPathOptions::Triple::none;
    REQUIRE(lengths_of(PointSourceField(rows, source, microphone, without), 3).empty());
    const std::vector<double> there = lengths_of(field, 3);
    const std::vector<double> back = lengths_of(PointSourceField(rows, microphone, source, options), 3);
    REQUIRE(back.size() == there.size());
    for (std::size_t i = 0; i < there.size(); i++)
        REQUIRE(back[i] == Catch::Approx(there[i]).epsilon(1e-12));
    PointPathOptions every = options;
    every.prune = false;
    REQUIRE(lengths_of(PointSourceField(rows, source, microphone, every), 3) == there);
}

// Three tall wedge-shaped blocks, apexes down, between a source and a
// microphone low on either side: the sound goes round the apexes, one
// vertical edge each.  As the middle apex moves across the line between the
// outer two, the ray round those is cut off and the ray round all three takes
// over, its middle coefficient's transition making up the step.  Every edge
// (Triple::all), the corners being vertical; without three edges the field
// steps there.  The outer blocks run twice as deep as the middle one, so that
// no ray round their far corners meets the middle one's far side as it moves;
// and the sweep steps past the apex exactly on the line, where three rays
// graze at once.  The edges' tops are ends, the pair's and the three's alike
// (twice_ends, chain_ends): 0.063 dB at worst (125 Hz, at the line) against
// 0.046 with no ends, and 0.15 with ends on the pair alone.
TEST_CASE("the field is continuous where a third edge cuts a twice-diffracted ray", "[point_paths]")
{
    PointPathOptions options;
    options.triple_diffraction = PointPathOptions::Triple::all;
    options.max_reflections = 0;
    options.max_edge_reflections = 0;
    options.max_twice_reflections = 0;
    options.max_between_reflections = 0;
    PointPathOptions without = options;
    without.triple_diffraction = PointPathOptions::Triple::none;
    const Vec3 source(-25.0, 40.0, 1.5), microphone(25.0, 40.0, 1.5);
    auto wedge = [](double x, double y, double depth) {
        Footprint f;
        f.x = {x, x + 0.08 * depth, x - 0.08 * depth};
        f.y = {y, y + depth, y + depth};
        f.top = 60.0;
        return f;
    };
    // The outer apexes at (-12, 10) and (12, 10), the middle one at (0, c).
    auto field_at = [&](double c, const PointPathOptions& o, std::size_t* thrice = nullptr) {
        const Buildings city({wedge(-12.0, 10.0, 200.0), wedge(0.0, c, 100.0), wedge(12.0, 10.0, 200.0)});
        if (thrice != nullptr)
            *thrice = lengths_of(PointSourceField(city, source, microphone, o), 3).size();
        return levels(city, source, microphone, 1.0, 0.0, o);
    };
    std::size_t below = 0, above = 0;
    field_at(9.6, options, &below);
    field_at(10.4, options, &above);
    REQUIRE(below > 0);
    REQUIRE(above > 0);
    double worst = 0.0, worst_without = 0.0;
    std::vector<double> previous = field_at(9.801, options), previous_without = field_at(9.801, without);
    for (int step = 1; step <= 50; step++)
    {
        const double c = 9.801 + 0.008 * step;
        const std::vector<double> now = field_at(c, options), now_without = field_at(c, without);
        for (std::size_t b = 0; b < 3; b++)
        {
            worst = std::max(worst, std::abs(now[b] - previous[b]));
            worst_without = std::max(worst_without, std::abs(now_without[b] - previous_without[b]));
        }
        previous = now;
        previous_without = now_without;
    }
    REQUIRE(worst_without > 1.0);
    REQUIRE(worst < 0.08);
}


// An edge's end makes no step (THEORY 6.5).  As a microphone rises past where
// the corner's Keller point leaves its top, the ray inside the edge hands over
// to the one past its end; and past the end, a microphone crossing the plane
// of the corner's face into the wedge's solid side, where the end still sees
// it, takes the face's coefficient weighed by its depth (weigh_vertex).  A
// micrometre either side, each within 1e-3 dB; with the edges unending the
// first is a step of decibels.
TEST_CASE("an edge's end is no step in the field", "[point_paths]")
{
    const Buildings city({box(0.0, 0.0, 20.0, 20.0, 6.0)});
    PointPathOptions options, unending;
    unending.edge_ends = false;
    const double hair = 1e-6;
    // The corner at the origin: its Keller point at z = (3 + z_m) / 2, the
    // source and microphone equally far from it, crosses its top at z_m = 9.
    const Vec3 source(-10.0, 15.0, 3.0);
    const auto step = [&](const Vec3& s, const Vec3& below, const Vec3& above, const PointPathOptions& o) {
        const std::vector<double> a = levels(city, s, below, 1.0, 0.0, o), b = levels(city, s, above, 1.0, 0.0, o);
        double worst = 0.0;
        for (std::size_t i = 0; i < bands.size(); i++)
        {
            REQUIRE(std::isfinite(a[i]));
            worst = std::max(worst, std::abs(a[i] - b[i]));
        }
        return worst;
    };
    const double past_top = step(source, Vec3(15.0, -10.0, 9.0 - hair), Vec3(15.0, -10.0, 9.0 + hair), options);
    const double without = step(source, Vec3(15.0, -10.0, 9.0 - hair), Vec3(15.0, -10.0, 9.0 + hair), unending);
    // Past the top, across the plane x = 0 of the corner's face.
    const double across_face
            = step(Vec3(-10.0, -8.0, 3.0), Vec3(-hair, 10.0, 14.0), Vec3(hair, 10.0, 14.0), options);
    INFO("past the top " << past_top << " dB (unending " << without << "), across the face's plane " << across_face);
    REQUIRE(past_top < 1e-3);
    REQUIRE(across_face < 1e-3);
    REQUIRE(without > 1.0);
}

// UTD with its edges' ends against the edge-source integral, which takes each
// edge's ends exactly: a building's corner, the microphone swept up past its
// top and across past its end.  Away from the planes of its faces through the
// vertex, where two ends' transitions meet and the per-term weights are not
// the vertex's own coefficient (FIDELITY), within 0.3 dB rms and 1 dB at
// worst (0.18 and 0.86 measured); with the edges unending, 1.4 dB rms and
// 7 dB at worst.
TEST_CASE("UTD with its edges' ends is the edge-source integral past a corner", "[point_paths][edge_source]")
{
    const Buildings city({box(0.0, 0.0, 10.0, 20.0, 6.0)});
    const Vec3 source(-15.0, 5.0, 3.0);
    const std::vector<double> swept {125.0, 250.0, 500.0, 1000.0, 2000.0};
    const double k = 2.0 * pi * swept.front() / 343.2;
    const auto field = [&](const Vec3& m, int model) {
        PointPathOptions o;
        o.edge_ends = model != 1;
        o.edge_source = model == 2;
        o.edge_source_reach = 100.0 / k;
        o.end_reach = end_fade_hi / k;
        const PointSourceField f(city, source, m, o);
        std::vector<double> out(swept.size());
        f.band_levels(std::vector<double>(city.facades().size() * swept.size(), 1.0), FloorModel(1.0), swept, 343.2,
                out.data());
        return out;
    };
    std::vector<Vec3> points;
    for (double z = 1.0; z <= 14.0 + 1e-9; z += 0.25)
        if (std::abs(z - 6.0) >= 0.5)
            points.push_back(Vec3(15.0, -5.0, z));
    for (double y = 8.0; y >= -8.0 - 1e-9; y -= 0.25)
        if (std::abs(y) >= 0.5)
            points.push_back(Vec3(15.0, y, 9.0));
    double sum = 0.0, sum_unending = 0.0, worst = 0.0, worst_unending = 0.0;
    std::size_t count = 0;
    for (const Vec3& m : points)
    {
        const std::vector<double> ends = field(m, 0), bare = field(m, 1), exact = field(m, 2);
        for (std::size_t b = 0; b < swept.size(); b++)
        {
            const double d = ends[b] - exact[b], u = bare[b] - exact[b];
            sum += d * d;
            sum_unending += u * u;
            worst = std::max(worst, std::abs(d));
            worst_unending = std::max(worst_unending, std::abs(u));
            count++;
        }
    }
    const double rms = std::sqrt(sum / count), rms_unending = std::sqrt(sum_unending / count);
    INFO("with ends " << rms << " dB rms, " << worst << " at worst; unending " << rms_unending << ", "
                      << worst_unending);
    REQUIRE(rms < 0.3);
    REQUIRE(worst < 1.0);
    REQUIRE(rms_unending > 1.0);
    REQUIRE(worst_unending > 5.0);
}

// A vertical edge standing on the floor goes on below it as its image, and the
// floor's images split that doubled edge between rays and their mirrors: of
// the four once-diffracted rays around a corner (the source or its image, the
// microphone or its), the tracer keeps the two whose points stand above the
// floor, and with the edge-source integral each takes the doubled edge.  Their
// sum is the doubled edge's field from the source and from its image.
TEST_CASE("with the edge source a corner on the floor is one edge with its image", "[point_paths][edge_source]")
{
    Footprint box;
    box.x = {0.0, 10.0, 10.0, 0.0};
    box.y = {0.0, 0.0, 10.0, 10.0};
    box.top = 8.0;
    const Buildings city({box});
    const Vec3 s(-6.0, 8.0, 1.5), m(8.0, -6.0, 1.0);
    PointPathOptions options;
    options.edge_source = true;
    options.double_diffraction = false;
    const PointSourceField field(city, s, m, options);
    const std::vector<double> reflection(city.facades().size(), 1.0);
    // The corner at the origin.
    std::size_t corner = city.edges().size();
    for (std::size_t e = 0; e < city.edges().size(); e++)
        if (city.edges()[e].vertical && std::hypot(city.edges()[e].start.x, city.edges()[e].start.y) < 1e-9)
            corner = e;
    REQUIRE(corner < city.edges().size());
    const Edge& g = city.edges()[corner];
    Edge doubled = g;
    doubled.start = g.start - g.length * g.direction;
    doubled.length = 2.0 * g.length;
    const Vec3 image(s.x, s.y, -s.z);
    for (double k : {3.0, 30.0})
    {
        std::complex<double> sum = 0.0;
        int count = 0;
        for (const PointRay& ray : field.rays())
            if (ray.wedges == 1 && ray.edge_source.edges[0] == corner && ray.facades.empty())
            {
                sum += reflected_amplitude(ray, reflection, 1.0) * edge_source_factors(city, ray, {k})[0]
                        * std::polar(1.0, k * ray.length);
                count++;
            }
        REQUIRE(count == 2);
        const std::complex<double> exact = edge_source_field(g.n, edge_frame_ends(doubled, s, m), k, 0.0, doubled.length)
                + edge_source_field(g.n, edge_frame_ends(doubled, image, m), k, 0.0, doubled.length);
        INFO("k " << k << ": " << std::abs(sum) << " against " << std::abs(exact));
        REQUIRE(std::abs(sum - exact) < 1e-4 * std::abs(exact));
    }
}

// Where a roof edge meets the vertical corner below it, the two edges share a
// vertex, and a twice-diffracted ray down the wall between them has a point of
// one edge on the other's line, which has no angle about it.  The partner
// search once took that point for a wave's boundary there, a stretch of no
// length whose once-diffracted field was not a number, and the whole band
// with it.  (RWTH EDGE's detached house, UE14, at model scale.)
TEST_CASE("with the edge source a corner's adjacent edges give a finite field", "[point_paths][edge_source]")
{
    Footprint house;
    house.x = {0.0, 0.96, 0.96, 0.0};
    house.y = {0.0, 0.0, 1.0, 1.0};
    house.top = 0.492;
    const Buildings city({house});
    PointPathOptions options;
    options.edge_source = true;
    const PointSourceField field(city, Vec3(-2.0, 0.5, 1.5), Vec3(1.46, 0.5, 0.0032), options);
    const std::vector<double> ks {2.0 * pi * 2000.0 / 343.2, 2.0 * pi * 10000.0 / 343.2};
    // The field's own factors fall back to the cascade where not finite, so
    // the integral is held directly, its vertical edge through the floor as
    // the factors take it (continued).
    int adjacent = 0;
    for (const PointRay& ray : field.rays())
    {
        if (ray.wedges != 2)
            continue;
        const Edge g1 = continued(city.edges()[ray.edge_source.edges[0]]);
        const Edge g2 = continued(city.edges()[ray.edge_source.edges[1]]);
        const Vec3 a1 = g1.start, b1 = g1.start + g1.length * g1.direction;
        const Vec3 a2 = g2.start, b2 = g2.start + g2.length * g2.direction;
        if (std::min({(a1 - a2).norm(), (a1 - b2).norm(), (b1 - a2).norm(), (b1 - b2).norm()}) > 1e-9)
            continue;
        adjacent++;
        for (const std::complex<double>& v : localized_twice_diffracted_spectrum(g1, g2, ray.edge_source.source,
                     ray.edge_source.microphone, ray.wedge[1].grazing ? 0.5 : 1.0, ks, PartnerModel::edge_source))
            REQUIRE((std::isfinite(v.real()) && std::isfinite(v.imag())));
    }
    REQUIRE(adjacent > 0);
}

// A path that bends round one block's corner and then the next block's
// opposite one (an S-bend) is a twice-diffracted ray of two parallel edges,
// and they are coupled: as the first edge's point moves, the path through the
// second is least at a point of the second's own.  With the low block's
// corner ending a metre or two above the ray and the tall one's ends far off,
// out of both transitions, the ray takes the cascade times each edge's end
// ratio, and these hold the integral to 0.3 dB (0.27 measured, at 4 kHz with
// the end 1.5 m off and within a Fresnel zone; 0.04 with it 2.5 m off);
// taken about the other edge's fixed stationary point they were 0.94 dB off.
// And each wavenumber's factor is its own, whatever others it is asked with:
// the cascade's place was once decided by the least, so a band moved with
// the run's band limits.
TEST_CASE("the edge source's end ratios hold for coupled edges, and each band is its own",
        "[point_paths][edge_source]")
{
    const Buildings city({box(0.0, 0.0, 10.0, 20.0, 40.0), box(12.0, -20.0, 22.0, 3.46, 200.0)});
    std::vector<double> ks;
    for (double f : {500.0, 1000.0, 2000.0, 4000.0})
        ks.push_back(2.0 * pi * f / 343.2);
    PointPathOptions options;
    options.edge_source = true;
    std::size_t checked = 0;
    for (double z : {37.5, 38.5})
    {
        const PointSourceField field(city, Vec3(-20.0, -1.0, z), Vec3(40.0, 10.0, z), options);
        for (const PointRay& ray : field.rays())
        {
            if (ray.wedges != 2 || ray.shared || ray.floor > 0 || ray.between != PointRay::no_facade
                    || !ray.facades.empty())
                continue;
            const Edge g1 = continued(city.edges()[ray.edge_source.edges[0]]);
            const Edge g2 = continued(city.edges()[ray.edge_source.edges[1]]);
            const std::vector<std::complex<double>> factors = edge_source_factors(city, ray, ks);
            const std::vector<std::complex<double>> exact = localized_twice_diffracted_spectrum(g1, g2,
                    ray.edge_source.source, ray.edge_source.microphone, 1.0, ks, PartnerModel::edge_source);
            REQUIRE(factors.size() == ks.size());
            for (std::size_t j = 0; j < ks.size(); j++)
            {
                const std::complex<double> reference = exact[j] * std::polar(1.0 / ray.spreading, -ks[j] * ray.length);
                INFO(ks[j] * 343.2 / (2.0 * pi) << " Hz, the ends " << z << " m up");
                REQUIRE(std::abs(20.0 * std::log10(std::abs(factors[j] / reference))) < 0.3);
                REQUIRE(std::abs(std::arg(factors[j] / reference)) < 0.03);
                // Asked alone, and with the lowest band of a run beside it.
                for (const std::vector<double>& with : {std::vector<double> {ks[j]},
                             std::vector<double> {2.0 * pi * 125.0 / 343.2, ks[j]}})
                {
                    const std::complex<double> v = edge_source_factors(city, ray, with).back();
                    REQUIRE(std::abs(20.0 * std::log10(std::abs(v / factors[j]))) < 0.01);
                }
            }
            checked++;
        }
    }
    REQUIRE(checked >= 2);
}
