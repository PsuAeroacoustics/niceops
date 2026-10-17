#include "point_paths.hpp"

#include <algorithm>
#include <cmath>

#include "urban.hpp"
#include "util/constants.hpp"
#include "wedge.hpp"

namespace niceops {

namespace {

Vec3 below(const Vec3& p)
{
    return {p.x, p.y, -p.z};
}

double horizontal_distance(const Vec3& p, const Edge& e)
{
    const Vec3 end = e.start + e.length * e.direction;
    const double dx = end.x - e.start.x, dy = end.y - e.start.y;
    const double len2 = dx * dx + dy * dy;
    const double t = len2 > 0.0 ? std::clamp(((p.x - e.start.x) * dx + (p.y - e.start.y) * dy) / len2, 0.0, 1.0) : 0.0;
    return std::hypot(p.x - (e.start.x + t * dx), p.y - (e.start.y + t * dy));
}

// Where on an edge's line the path from a to b is stationary: the point that
// leaves the edge toward b at the angle a arrives at.  Parameter along the
// edge from its start; false where either point lies on the line.
bool keller(const Edge& e, const Vec3& a, const Vec3& b, double& t)
{
    const Vec3 da = a - e.start, db = b - e.start;
    const double ta = da.dot(e.direction), tb = db.dot(e.direction);
    const double ra = (da - ta * e.direction).norm(), rb = (db - tb * e.direction).norm();
    if (ra < 1e-9 || rb < 1e-9)
        return false;
    t = ta + (tb - ta) * ra / (ra + rb);
    return true;
}

// The component of v normal to the edge.
Vec3 normal_part(const Edge& e, const Vec3& v)
{
    return v - v.dot(e.direction) * e.direction;
}

// Whether an angle lies in the wedge's air, strictly or with its faces.
bool in_air(double theta, double n, bool faces)
{
    const double slack = 1e-9;
    return faces ? (theta > -slack && theta < n * pi + slack) || theta > 2.0 * pi - slack
                 : theta > slack && theta < n * pi - slack;
}

// An angle on a face (within round-off) taken exactly on it.
double onto_face(double theta, double n, bool& on_face)
{
    const double slack = 1e-7;
    on_face = true;
    if (theta < slack || theta > 2.0 * pi - slack)
        return 0.0;
    if (std::abs(theta - n * pi) < slack)
        return n * pi;
    on_face = false;
    return theta;
}

} // namespace

PointSourceField::PointSourceField(const Buildings& buildings, const Vec3& source, const Vec3& microphone,
        const PointPathOptions& options)
    : facade_count_(buildings.facades().size())
{
    const auto& edges = buildings.edges();
    std::vector<std::size_t> near_edges;
    for (std::size_t e = 0; e < edges.size(); e++)
        if (horizontal_distance(source, edges[e]) <= options.radius
                || horizontal_distance(microphone, edges[e]) <= options.radius)
            near_edges.push_back(e);
    const std::size_t receiver_order = std::max(options.max_reflections, options.max_edge_reflections);
    const ImageTree mic = image_tree(buildings, microphone, options.radius, receiver_order, options.max_images, true,
            near_edges);
    const ImageTree src = image_tree(buildings, source, options.radius,
            options.diffraction && options.source_reflections ? options.max_edge_reflections : 0, options.max_images,
            true, near_edges);
    const int source_floors = options.source_reflections ? 2 : 1;

    // Specular rays.
    for (const PointImage& image : mic.images)
    {
        if (image.facades.size() > options.max_reflections)
            continue;
        for (int floor = 0; floor < 2; floor++)
        {
            std::vector<Vec3> hits;
            if (!folded_path_open(buildings, source, microphone, image.facades, floor == 1, &hits))
                continue;
            PointRay ray;
            ray.path.push_back(source);
            ray.path.insert(ray.path.end(), hits.begin(), hits.end());
            ray.path.push_back(microphone);
            ray.length = ((floor ? below(image.point) : image.point) - source).norm();
            ray.spreading = 1.0 / ray.length;
            ray.facades = image.facades;
            ray.floor = floor;
            rays_.push_back(std::move(ray));
        }
    }
    if (!options.diffraction)
        return;

    // Diffracted once, between an image of the source and one of the
    // microphone.
    for (std::size_t i = 0; i < src.images.size(); i++)
        for (std::size_t j = 0; j < mic.images.size(); j++)
        {
            const PointImage& a = src.images[i];
            const PointImage& b = mic.images[j];
            if (b.facades.size() > options.max_edge_reflections)
                continue;
            for (std::size_t e : src.edges[i])
            {
                if (std::find(mic.edges[j].begin(), mic.edges[j].end(), e) == mic.edges[j].end())
                    continue;
                const Edge& edge = edges[e];
                for (int fa = 0; fa < source_floors; fa++)
                    for (int fb = 0; fb < 2; fb++)
                    {
                        const Vec3 pa = fa ? below(a.point) : a.point, pb = fb ? below(b.point) : b.point;
                        double t;
                        if (!keller(edge, pa, pb, t) || t < 0.0 || t > edge.length)
                            continue;
                        const Vec3 q = edge.start + t * edge.direction;
                        if (q.z < 0.0 || buildings.buried(q))
                            continue;
                        const Vec3 in = pa - q, out = pb - q;
                        const double theta0 = edge.angle(normal_part(edge, in));
                        const double theta = edge.angle(normal_part(edge, out));
                        if (!in_air(theta0, edge.n, false) || !in_air(theta, edge.n, false))
                            continue;
                        std::vector<Vec3> to_source, to_microphone;
                        if (!folded_path_open(buildings, q, source, a.facades, fa == 1, &to_source)
                                || !folded_path_open(buildings, q, microphone, b.facades, fb == 1, &to_microphone))
                            continue;
                        const double s_in = in.norm(), s_out = out.norm();
                        const double sin_beta = normal_part(edge, in).norm() / s_in;
                        PointRay ray;
                        ray.length = s_in + s_out;
                        ray.spreading = 1.0 / std::sqrt(s_in * s_out * (s_in + s_out));
                        ray.facades = a.facades;
                        ray.facades.insert(ray.facades.end(), b.facades.begin(), b.facades.end());
                        ray.floor = fa + fb;
                        ray.wedges = 1;
                        ray.wedge[0] = {edge.n, theta, theta0, sin_beta, s_in * s_out * sin_beta / (s_in + s_out), false};
                        ray.path = {source};
                        ray.path.insert(ray.path.end(), to_source.rbegin(), to_source.rend());
                        ray.path.push_back(q);
                        ray.path.insert(ray.path.end(), to_microphone.begin(), to_microphone.end());
                        ray.path.push_back(microphone);
                        rays_.push_back(std::move(ray));
                    }
            }
        }
    if (!options.double_diffraction)
        return;

    // Diffracted twice: the source's images and the microphone's to one
    // reflection, the two points found in turn until the path is stationary.
    for (std::size_t i = 0; i < src.images.size(); i++)
        for (std::size_t j = 0; j < mic.images.size(); j++)
        {
            const PointImage& a = src.images[i];
            const PointImage& b = mic.images[j];
            if (a.facades.size() > 1 || b.facades.size() > 1)
                continue;
            for (std::size_t e1 : src.edges[i])
                for (std::size_t e2 : mic.edges[j])
                {
                    if (e1 == e2)
                        continue;
                    const Edge& g1 = edges[e1];
                    const Edge& g2 = edges[e2];
                    for (int fa = 0; fa < source_floors; fa++)
                        for (int fb = 0; fb < 2; fb++)
                        {
                            const Vec3 pa = fa ? below(a.point) : a.point, pb = fb ? below(b.point) : b.point;
                            double t1 = 0.5 * g1.length, t2 = 0.5 * g2.length;
                            bool ok = true;
                            for (int iteration = 0; iteration < 200 && ok; iteration++)
                            {
                                const Vec3 q2 = g2.start + t2 * g2.direction;
                                double n1, n2;
                                ok = keller(g1, pa, q2, n1);
                                if (!ok)
                                    break;
                                const Vec3 q1 = g1.start + n1 * g1.direction;
                                ok = keller(g2, q1, pb, n2);
                                if (!ok)
                                    break;
                                const double change = std::max(std::abs(n1 - t1), std::abs(n2 - t2));
                                t1 = n1;
                                t2 = n2;
                                if (change < 1e-10)
                                    break;
                            }
                            if (!ok || t1 < 0.0 || t1 > g1.length || t2 < 0.0 || t2 > g2.length)
                                continue;
                            const Vec3 q1 = g1.start + t1 * g1.direction, q2 = g2.start + t2 * g2.direction;
                            if (q1.z < 0.0 || q2.z < 0.0 || buildings.buried(q1) || buildings.buried(q2))
                                continue;
                            const Vec3 in1 = pa - q1, mid = q2 - q1, out2 = pb - q2;
                            const double s1 = in1.norm(), s2 = mid.norm(), s3 = out2.norm();
                            if (s2 < 1e-6)
                                continue;
                            // The ray between the edges along a face both wedges share
                            // (a roof between its two edges, a wall between two
                            // corners): wedge 1's wave there is already the field
                            // doubled by that face, and wedge 2's coefficient, its
                            // incident and reflected boundaries coinciding at grazing,
                            // would double it again.
                            bool leaves_on_face, arrives_on_face;
                            const double theta0_1 = g1.angle(normal_part(g1, in1));
                            const double theta_1 = onto_face(g1.angle(normal_part(g1, mid)), g1.n, leaves_on_face);
                            const double theta0_2 = onto_face(g2.angle(normal_part(g2, -1.0 * mid)), g2.n, arrives_on_face);
                            const bool grazing = leaves_on_face && arrives_on_face;
                            const double theta_2 = g2.angle(normal_part(g2, out2));
                            if (!in_air(theta0_1, g1.n, false) || !in_air(theta_1, g1.n, true)
                                    || !in_air(theta0_2, g2.n, true) || !in_air(theta_2, g2.n, false))
                                continue;
                            std::vector<Vec3> to_source, to_microphone;
                            if (!folded_path_open(buildings, q1, source, a.facades, fa == 1, &to_source)
                                    || buildings.blocked(q1, q2)
                                    || !folded_path_open(buildings, q2, microphone, b.facades, fb == 1, &to_microphone))
                                continue;
                            const double sb1 = normal_part(g1, in1).norm() / s1;
                            const double sb2 = normal_part(g2, out2).norm() / s3;
                            PointRay ray;
                            ray.length = s1 + s2 + s3;
                            ray.spreading = 1.0 / std::sqrt(s1 * s2 * s3 * (s1 + s2 + s3));
                            ray.facades = a.facades;
                            ray.facades.insert(ray.facades.end(), b.facades.begin(), b.facades.end());
                            ray.floor = fa + fb;
                            ray.wedges = 2;
                            // Wedge 2's incident wave is curved about wedge 1 (radius s2)
                            // across the edge and about the source (s1 + s2) along it:
                            // its distance parameter takes the first, its spreading the
                            // second.  Wedge 1's observer is wedge 2.
                            ray.wedge[0] = {g1.n, theta_1, theta0_1, sb1, s1 * s2 * sb1 / (s1 + s2), false};
                            ray.wedge[1] = {g2.n, theta_2, theta0_2, sb2, s2 * s3 * sb2 / (s2 + s3), grazing};
                            ray.path = {source};
                            ray.path.insert(ray.path.end(), to_source.rbegin(), to_source.rend());
                            ray.path.push_back(q1);
                            ray.path.push_back(q2);
                            ray.path.insert(ray.path.end(), to_microphone.begin(), to_microphone.end());
                            ray.path.push_back(microphone);
                            if (grazing && options.shared_face && std::abs(g1.direction.dot(g2.direction)) > 1.0 - 1e-9)
                            {
                                // The plane normal to the edges: x from edge 1 toward
                                // edge 2, z into the air off the face.
                                const Vec3 across = normal_part(g1, mid);
                                const double width = across.norm();
                                const Vec3 x = across / width;
                                Vec3 z = g1.direction.cross(x);
                                if (!in_air(g1.angle(normal_part(g1, x + 1e-3 * z)), g1.n, false))
                                    z = -1.0 * z;
                                const Vec3 a = pa - q1, b = pb - q1;
                                // There each wedge's angles run from the shared face;
                                // the hard wedge is symmetric, so which of its faces
                                // that is does not matter.
                                SharedFaceTerm term;
                                term.face = {width, g1.n, g2.n};
                                term.sx = a.dot(x);
                                term.sz = a.dot(z);
                                term.rx = b.dot(x);
                                term.rz = b.dot(z);
                                term.sin_beta = sb1;
                                if (shared_face_applies(term.face, term.sx, term.sz, term.rx, term.rz))
                                {
                                    ray.shared = true;
                                    ray.shared_face = term;
                                    ray.spreading = 1.0 / std::sqrt(ray.length * sb1);
                                }
                            }
                            rays_.push_back(std::move(ray));
                        }
                }
        }
}

std::complex<double> point_ray_amplitude(const PointRay& ray, const std::vector<double>& reflection,
        double floor_reflection, double k)
{
    std::complex<double> amplitude = ray.spreading;
    for (std::size_t f : ray.facades)
        amplitude *= std::sqrt(reflection[f]);
    for (int i = 0; i < ray.floor; i++)
        amplitude *= std::sqrt(floor_reflection);
    if (ray.shared)
    {
        // The 2-D field at k sin(beta), its phase along the in-plane length
        // taken out (the ray's length carries it), spread along the edges:
        // a wave exp(i k_t l_t) / sqrt(l_t) in the plane is exp(i k l) / l in
        // space, l = l_t / sin(beta), so the factor is sin(beta) / sqrt(l_t).
        const SharedFaceTerm& t = ray.shared_face;
        const double kt = k * t.sin_beta;
        return amplitude * t.sin_beta
                * shared_face_field(t.face, t.sx, t.sz, t.rx, t.rz, kt)
                * std::polar(1.0, -kt * ray.length * t.sin_beta);
    }
    for (int w = 0; w < ray.wedges; w++)
    {
        const WedgeTerm& t = ray.wedge[w];
        // The wedge's coefficient in the plane normal to the edge, at the
        // wavenumber and distance projected into it: for a spherical wave the
        // distance is L / sin(beta), L = s s' sin^2(beta) / (s + s').
        amplitude *= hard_wedge_coefficient(t.n, t.theta, t.theta0, k * t.sin_beta, t.distance) / std::sqrt(t.sin_beta)
                * (t.grazing ? 0.5 : 1.0);
    }
    return amplitude;
}

void PointSourceField::band_levels(const std::vector<double>& reflection, double floor_reflection,
        const std::vector<double>& bands, double sound_speed, double* out) const
{
    std::vector<PointRay> sorted = rays_;
    std::sort(sorted.begin(), sorted.end(), [](const PointRay& a, const PointRay& b) { return a.length < b.length; });
    std::vector<double> length(sorted.size());
    for (std::size_t i = 0; i < sorted.size(); i++)
        length[i] = sorted[i].length;
    std::vector<std::complex<double>> amplitude(sorted.size());
    std::vector<double> band_reflection(facade_count_);
    const double half = std::pow(2.0, 1.0 / 6.0);
    for (std::size_t band = 0; band < bands.size(); band++)
    {
        for (std::size_t f = 0; f < facade_count_; f++)
            band_reflection[f] = reflection[f * bands.size() + band];
        const double f1 = bands[band] / half, f2 = bands[band] * half;
        double energy = 0.0;
        for (int part = 0; part < 3; part++)
        {
            const double g1 = f1 * std::pow(f2 / f1, part / 3.0), g2 = f1 * std::pow(f2 / f1, (part + 1) / 3.0);
            const double k1 = 2.0 * pi * g1 / sound_speed, k2 = 2.0 * pi * g2 / sound_speed;
            const double kc = 0.5 * (k1 + k2);
            for (std::size_t i = 0; i < sorted.size(); i++)
                amplitude[i] = point_ray_amplitude(sorted[i], band_reflection, floor_reflection, kc);
            energy += band_averaged_energy(amplitude, length, k1, k2) * (g2 - g1);
        }
        out[band] = 10.0 * std::log10(std::max(energy / (f2 - f1), 1e-30));
    }
}

} // namespace niceops
