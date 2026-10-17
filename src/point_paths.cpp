#include "point_paths.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iterator>
#include <stdexcept>

#include "urban.hpp"
#include "util/constants.hpp"
#include "wedge.hpp"

namespace niceops {

namespace {

Vec3 below(const Vec3& p)
{
    return {p.x, p.y, -p.z};
}

// Where the path through a point on a line is stationary, the point leaving
// toward b at the angle a arrives at: from a's and b's places about the line
// (along it, and their distances off it), the parameter along it.  False
// where either lies on the line.
bool keller_point(double ta, double ra, double tb, double rb, double& t)
{
    if (ra < 1e-9 || rb < 1e-9)
        return false;
    t = ta + (tb - ta) * ra / (ra + rb);
    return true;
}

// The same for the points themselves: where on an edge's line the path from
// a to b is stationary.
bool keller(const Edge& e, const Vec3& a, const Vec3& b, double& t)
{
    const Edge::Projection pa = e.project(a), pb = e.project(b);
    return keller_point(pa.along, pa.normal.norm(), pb.along, pb.normal.norm(), t);
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

// in_air(e.angle(v), e.n, faces), by the signs where they decide it.
bool air(const Edge& e, const Vec3& v, bool faces)
{
    const int sign = e.air_by_sign(v);
    return sign >= 0 ? sign == 1 : in_air(e.angle(v), e.n, faces);
}

// in_air(onto_face(e.angle(v), e.n, on_face), e.n, true): a direction on a
// face, within onto_face's 1e-7, counts as in the air.
bool air_or_face(const Edge& e, const Vec3& v)
{
    const int sign = e.air_by_sign(v);
    if (sign >= 0)
        return sign == 1;
    bool on_face;
    return in_air(onto_face(e.angle(v), e.n, on_face), e.n, true);
}

// Where a point stands about an edge, which does not depend on where along
// the edge it is seen from: its parameter along the edge's line, its distance
// from that line, and whether it lies in the wedge's air, strictly.  The air
// test is a hair wider than in_air's, for the angle taken from the edge's
// start rather than from a point along it.
struct Seen
{
    double along = 0.0;
    double off = 0.0;
    bool lit = false;
};

Seen seen_from(const Edge& e, const Vec3& p)
{
    const Edge::Projection place = e.project(p);
    Seen s;
    s.along = place.along;
    s.off = place.normal.norm();
    const int sign = e.air_by_sign(place.normal);
    if (s.off < 1e-4 || sign == 1)
        s.lit = true;
    else if (sign == -1)
    {
        const double theta = e.angle(place.normal);
        const double slack = 1e-9 - 1e-10;
        s.lit = theta > slack && theta < e.n * pi - slack;
    }
    return s;
}

// How one edge's points stand about another edge's line.  The point at
// parameter t of the segment lies along the line at c0 + c1 t and off it by
// |n0 + t nd|; over the segment those range over [along_min, along_max] and
// [off_min, off_max].
struct Side
{
    double c0 = 0.0, c1 = 0.0;
    Vec3 n0, nd;
    double along_min = 0.0, along_max = 0.0;
    double off_min = 0.0, off_max = 0.0;

    // keller_point's places for the point at parameter x.
    void at(double x, double& along, double& off) const
    {
        along = c0 + c1 * x;
        off = (n0 + x * nd).norm();
    }
};

Side side_of(const Edge& line, const Edge& segment)
{
    const Edge::Projection start = line.project(segment.start);
    Side s;
    s.c0 = start.along;
    s.c1 = segment.direction.dot(line.direction);
    s.n0 = start.normal;
    s.nd = normal_part(line, segment.direction);
    const double end = s.c0 + segment.length * s.c1;
    s.along_min = std::min(s.c0, end);
    s.along_max = std::max(s.c0, end);
    const double nn = s.nd.dot(s.nd);
    const double closest = nn > 0.0 ? std::clamp(-s.n0.dot(s.nd) / nn, 0.0, segment.length) : 0.0;
    s.off_min = (s.n0 + closest * s.nd).norm();
    s.off_max = std::max(s.n0.norm(), (s.n0 + segment.length * s.nd).norm());
    return s;
}

// Two edges about each other, for the twice-diffracted search: the second's
// points about the first's line and the first's about the second's.
struct EdgePair
{
    Side on1, on2;
    bool parallel = false;
};

EdgePair edge_pair(const Edge& g1, const Edge& g2)
{
    EdgePair p;
    p.on1 = side_of(g1, g2);
    p.on2 = side_of(g2, g1);
    p.parallel = std::abs(p.on2.c1) > 1.0 - 1e-12;
    return p;
}

// The edge pairs a field has met, by their two edges: worked out once per
// field however many images meet them.  Open addressing, and a generation per
// slot so that the next field starts empty without clearing; its memory
// follows the pairs met, not the square of the edges near.  A reference it
// returns lasts until the next get().
class PairTable
{
public:
    void start()
    {
        pairs_.clear();
        if (++generation_ == 0)
        {
            std::fill(stamp_.begin(), stamp_.end(), 0u);
            generation_ = 1;
        }
    }

    const EdgePair& get(const std::vector<Edge>& edges, std::size_t e1, std::size_t e2)
    {
        if (2 * (pairs_.size() + 1) > key_.size())
            grow();
        const std::uint64_t key = (static_cast<std::uint64_t>(e1) << 32) | static_cast<std::uint64_t>(e2);
        std::size_t i = slot(key);
        while (stamp_[i] == generation_)
        {
            if (key_[i] == key)
                return pairs_[value_[i]];
            i = (i + 1) & (key_.size() - 1);
        }
        stamp_[i] = generation_;
        key_[i] = key;
        value_[i] = static_cast<std::uint32_t>(pairs_.size());
        pairs_.push_back(edge_pair(edges[e1], edges[e2]));
        return pairs_.back();
    }

private:
    std::size_t slot(std::uint64_t key) const
    {
        return static_cast<std::size_t>((key * 0x9E3779B97F4A7C15ull) >> shift_);
    }

    void grow()
    {
        const std::vector<std::uint64_t> key = std::move(key_);
        const std::vector<std::uint32_t> stamp = std::move(stamp_), value = std::move(value_);
        std::size_t size = 1024, bits = 10;
        while (size < 2 * key.size())
        {
            size *= 2;
            bits++;
        }
        shift_ = 64 - static_cast<unsigned>(bits);
        key_.assign(size, 0);
        stamp_.assign(size, 0);
        value_.assign(size, 0);
        for (std::size_t j = 0; j < key.size(); j++)
            if (stamp[j] == generation_)
            {
                std::size_t i = slot(key[j]);
                while (stamp_[i] == generation_)
                    i = (i + 1) & (size - 1);
                stamp_[i] = generation_;
                key_[i] = key[j];
                value_[i] = value[j];
            }
    }

    std::vector<std::uint64_t> key_;
    std::vector<std::uint32_t> stamp_, value_;
    std::vector<EdgePair> pairs_;
    std::uint32_t generation_ = 0;
    unsigned shift_ = 54;
};

// Whether keller(line, p, q) can fall within the line's edge, give or take
// `slack`, for some q on the segment whose Side about the line is `side`: the
// stationary point lies along the line from p at the fraction
// off_p / (off_p + off_q) of the way to q's foot, so it lies within the hull of
// those fractions of the feet's range -- and, whatever the distances, between
// p's foot and the feet (keller_between, the cheaper first test).
bool keller_between(const Edge& line, const Seen& p, const Side& side, double slack)
{
    return std::max(p.along, side.along_max) >= -slack && std::min(p.along, side.along_min) <= line.length + slack;
}

bool keller_can_reach(const Edge& line, const Seen& p, const Side& side, double slack)
{
    if (p.off < 1e-4)
        return true;
    const double near = p.off / (p.off + side.off_min), far = p.off / (p.off + side.off_max);
    const double lo = side.along_min - p.along, hi = side.along_max - p.along;
    const double reach_min = std::min(std::min(lo * near, lo * far), std::min(hi * near, hi * far));
    const double reach_max = std::max(std::max(lo * near, lo * far), std::max(hi * near, hi * far));
    return p.along + reach_max >= -slack && p.along + reach_min <= line.length + slack;
}

// The stationary point between two parallel edges' lines in closed form: the
// path unfolds about the common direction into a straight line, so its
// travel along it is shared in proportion to the distances across.  The
// second edge's parameter; false where the distances leave it undefined.
bool parallel_stationary(const Edge& g1, const Edge& g2, const Vec3& pa, const Vec3& pb, double& t2)
{
    const Edge::Projection a = g1.project(pa), c = g1.project(g2.start);
    const double zb = g1.project(pb).along;
    const double ra = a.normal.norm(), rc = c.normal.norm(), rb = g2.project(pb).normal.norm();
    const double total = ra + rc + rb;
    if (!(total > 0.0))
        return false;
    const double z2 = a.along + (zb - a.along) * (ra + rc) / total;
    t2 = (z2 - c.along) * g2.direction.dot(g1.direction);
    return std::isfinite(t2);
}

std::size_t tree_order(bool source, const PointPathOptions& options)
{
    if (source)
        return options.diffraction && options.source_reflections ? options.max_edge_reflections : 0;
    return std::max(options.max_reflections, options.max_edge_reflections);
}

std::vector<std::size_t> joined(const std::vector<std::size_t>& a, const std::vector<std::size_t>& b)
{
    std::vector<std::size_t> out;
    out.reserve(a.size() + b.size());
    std::set_union(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));
    return out;
}

// The partners as set_union and includes need them: ascending, no repeats.
const std::vector<std::size_t>& ascending(const std::vector<std::size_t>& partners)
{
    if (std::adjacent_find(partners.begin(), partners.end(), std::greater_equal<std::size_t>()) != partners.end())
        throw std::invalid_argument("PathEnd: the partners' edges must be ascending, without repeats (edges_near)");
    return partners;
}

// The deepest images whose edges a field reads: a single diffraction's, and
// the twice-diffracted search's, which takes the microphone's to one
// reflection whatever max_edge_reflections is.
std::size_t edge_order(const PointPathOptions& options)
{
    return std::max(options.max_edge_reflections,
            options.diffraction && options.double_diffraction ? std::size_t {1} : std::size_t {0});
}

// Whether two sets of options trace the same image trees: every option
// PathEnd reads.
bool same_tracing(const PointPathOptions& a, const PointPathOptions& b)
{
    return a.max_reflections == b.max_reflections && a.max_edge_reflections == b.max_edge_reflections
            && a.diffraction == b.diffraction && a.double_diffraction == b.double_diffraction
            && a.source_reflections == b.source_reflections && a.radius == b.radius && a.max_images == b.max_images;
}

} // namespace

PathEnd::PathEnd(const Buildings& buildings, const Vec3& point, bool source, const PointPathOptions& options,
        const std::vector<std::size_t>& partners)
    : point_(point)
    , source_(source)
    , options_(options)
    , near_(edges_near(buildings, point, options.radius))
    , candidates_(joined(near_, ascending(partners)))
    , tree_(image_tree(buildings, point, options.radius, tree_order(source, options), options.max_images, true,
              candidates_, edge_order(options)))
{
}

PointSourceField::PointSourceField(const Buildings& buildings, const Vec3& source, const Vec3& microphone,
        const PointPathOptions& options)
    : facade_count_(buildings.facades().size())
{
    trace(buildings, PathEnd(buildings, source, true, options, edges_near(buildings, microphone, options.radius)),
            PathEnd(buildings, microphone, false, options, edges_near(buildings, source, options.radius)), options);
}

PointSourceField::PointSourceField(const Buildings& buildings, const PathEnd& source, const PathEnd& microphone,
        const PointPathOptions& options)
    : facade_count_(buildings.facades().size())
{
    trace(buildings, source, microphone, options);
}

void PointSourceField::trace(const Buildings& buildings, const PathEnd& source_end, const PathEnd& microphone_end,
        const PointPathOptions& options)
{
    // The edges near this pair; each end was traced with them all, and each
    // image keeps those that cross its beam, in the order the pair's own tree
    // would list them.
    if (!source_end.source() || microphone_end.source() || !same_tracing(source_end.options(), options)
            || !same_tracing(microphone_end.options(), options))
        throw std::invalid_argument("PointSourceField: the ends must be a source and a microphone traced with the "
                                    "field's options");
    const std::vector<std::size_t> near_edges = joined(source_end.near(), microphone_end.near());
    for (const PathEnd* end : {&source_end, &microphone_end})
        if (!std::includes(end->candidates().begin(), end->candidates().end(), near_edges.begin(), near_edges.end()))
            throw std::invalid_argument("a PathEnd was traced without the edges near the end it is paired with");
    const auto& edges = buildings.edges();
    thread_local std::vector<char> is_near;
    is_near.resize(edges.size(), 0);
    for (std::size_t e : near_edges)
        is_near[e] = 1;
    auto keep = [&](const ImageTree& tree) {
        std::vector<std::vector<std::size_t>> kept(tree.edges.size());
        for (std::size_t i = 0; i < tree.edges.size(); i++)
            for (std::size_t e : tree.edges[i])
                if (is_near[e])
                    kept[i].push_back(e);
        return kept;
    };
    const std::vector<std::vector<std::size_t>> src_edges = keep(source_end.tree()),
                                                mic_edges = keep(microphone_end.tree());
    for (std::size_t e : near_edges)
        is_near[e] = 0;
    const Ends src {source_end.tree().images, src_edges}, mic {microphone_end.tree().images, mic_edges};
    const Vec3& source = source_end.point();
    const Vec3& microphone = microphone_end.point();
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
            // The edges both images keep: both lists run in near_edges' order,
            // so one pass along the two finds them.
            std::size_t next = 0;
            const auto& keep = mic.edges[j];
            for (std::size_t e : src.edges[i])
            {
                while (next < keep.size() && keep[next] < e)
                    next++;
                if (next == keep.size() || keep[next] != e)
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
                        if (!air(edge, normal_part(edge, in), false) || !air(edge, normal_part(edge, out), false))
                            continue;
                        std::vector<Vec3> to_source, to_microphone;
                        if (!folded_path_open(buildings, q, source, a.facades, fa == 1, &to_source)
                                || !folded_path_open(buildings, q, microphone, b.facades, fb == 1, &to_microphone))
                            continue;
                        const double theta0 = edge.angle(normal_part(edge, in));
                        const double theta = edge.angle(normal_part(edge, out));
                        const double s_in = in.norm(), s_out = out.norm();
                        const double sin_beta = normal_part(edge, in).norm() / s_in;
                        PointRay ray;
                        ray.length = s_in + s_out;
                        ray.spreading = 1.0 / std::sqrt(s_in * s_out * (s_in + s_out));
                        ray.facades = a.facades;
                        ray.facades.insert(ray.facades.end(), b.facades.begin(), b.facades.end());
                        ray.floor = fa + fb;
                        ray.wedges = 1;
                        ray.wedge[0] = {wedge_angles(edge.n, theta, theta0), sin_beta,
                                s_in * s_out * sin_beta / (s_in + s_out), false};
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
    // reflection.  With t2 given, the stationary point on the first edge is
    // keller's, and with that the second's; the path is stationary where t2
    // is a fixed point of the two in turn.  Most pairs are ruled out first:
    // where the source's image does not light the first edge or the
    // microphone's does not see the second, from either side of the edge
    // whatever the point along it, and where no point of one edge can put the
    // other's stationary point on its own segment.
    // Per image and floor (2 image + floor), each of its edges' Seen.
    auto seen_all = [&](const Ends& tree, std::size_t image, const Vec3& p) {
        std::vector<Seen> seen;
        for (std::size_t e : tree.edges[image])
            seen.push_back(seen_from(edges[e], p));
        return seen;
    };
    std::vector<std::size_t> src_first, mic_first;
    for (std::size_t i = 0; i < src.images.size(); i++)
        if (src.images[i].facades.size() <= 1)
            src_first.push_back(i);
    for (std::size_t j = 0; j < mic.images.size(); j++)
        if (mic.images[j].facades.size() <= 1)
            mic_first.push_back(j);
    std::vector<std::vector<Seen>> src_lit(src_first.size() * 2), mic_lit(mic_first.size() * 2);
    for (std::size_t a = 0; a < src_first.size(); a++)
        for (int fa = 0; fa < source_floors; fa++)
        {
            const Vec3& p = src.images[src_first[a]].point;
            src_lit[a * 2 + fa] = seen_all(src, src_first[a], fa ? below(p) : p);
        }
    for (std::size_t b = 0; b < mic_first.size(); b++)
        for (int fb = 0; fb < 2; fb++)
        {
            const Vec3& p = mic.images[mic_first[b]].point;
            mic_lit[b * 2 + fb] = seen_all(mic, mic_first[b], fb ? below(p) : p);
        }
    const double reach_slack = 1e-6;
    thread_local PairTable pairs;
    pairs.start();
    std::vector<Vec3> to_source, to_microphone;
    for (std::size_t ia = 0; ia < src_first.size(); ia++)
        for (std::size_t jb = 0; jb < mic_first.size(); jb++)
        {
            const std::size_t i = src_first[ia], j = mic_first[jb];
            const PointImage& a = src.images[i];
            const PointImage& b = mic.images[j];
            for (std::size_t k1 = 0; k1 < src.edges[i].size(); k1++)
            {
                const std::size_t e1 = src.edges[i][k1];
                if (options.prune && !src_lit[ia * 2][k1].lit && (source_floors == 1 || !src_lit[ia * 2 + 1][k1].lit))
                    continue;
                for (std::size_t k2 = 0; k2 < mic.edges[j].size(); k2++)
                {
                    const std::size_t e2 = mic.edges[j][k2];
                    if (e1 == e2
                            || (options.prune && !mic_lit[jb * 2][k2].lit && !mic_lit[jb * 2 + 1][k2].lit))
                        continue;
                    const Edge& g1 = edges[e1];
                    const Edge& g2 = edges[e2];
                    const EdgePair& pair = pairs.get(edges, e1, e2);
                    // Which of the source's images can reach the first edge's
                    // segment, and the microphone's the second's, floor by
                    // floor: the pre-checks for one end do not depend on the
                    // other's.
                    bool reach_source[2] = {true, true}, reach_microphone[2] = {true, true};
                    if (options.prune)
                    {
                        for (int fa = 0; fa < source_floors; fa++)
                        {
                            const Seen& p = src_lit[ia * 2 + fa][k1];
                            reach_source[fa] = p.lit && keller_between(g1, p, pair.on1, reach_slack)
                                    && keller_can_reach(g1, p, pair.on1, reach_slack);
                        }
                        for (int fb = 0; fb < 2; fb++)
                        {
                            const Seen& p = mic_lit[jb * 2 + fb][k2];
                            reach_microphone[fb] = p.lit && keller_between(g2, p, pair.on2, reach_slack)
                                    && keller_can_reach(g2, p, pair.on2, reach_slack);
                        }
                    }
                    for (int fa = 0; fa < source_floors; fa++)
                        for (int fb = 0; fb < 2; fb++)
                        {
                            if (!reach_source[fa] || !reach_microphone[fb])
                                continue;
                            const Seen& from_source = src_lit[ia * 2 + fa][k1];
                            const Seen& from_microphone = mic_lit[jb * 2 + fb][k2];
                            const Vec3 pa = fa ? below(a.point) : a.point, pb = fb ? below(b.point) : b.point;
                            // The fixed point of t2 -> keller(g2, q1(keller(g1, pa, q2(t2))), pb),
                            // by Steffensen's method: exact in one step for parallel
                            // edges, where the map is affine, and quadratic near the
                            // point otherwise; a step that does not bring the
                            // residual down falls back to the map itself.  `first`
                            // is the first edge's parameter on the way.
                            auto map = [&](double x, double& first, double& second) {
                                double along, off;
                                pair.on1.at(x, along, off);
                                if (!keller_point(from_source.along, from_source.off, along, off, first))
                                    return false;
                                pair.on2.at(first, along, off);
                                return keller_point(along, off, from_microphone.along, from_microphone.off, second);
                            };
                            double t1 = 0.0, t2 = 0.5 * g2.length, first = 0.0, x1 = 0.0;
                            bool converged = pair.parallel && parallel_stationary(g1, g2, pa, pb, t2);
                            bool ok = converged || map(t2, first, x1);
                            for (int iteration = 0; iteration < 100 && ok && !converged; iteration++)
                            {
                                // x1 is the map of t2.
                                if (std::abs(x1 - t2) < 1e-10)
                                {
                                    t2 = x1;
                                    converged = true;
                                    break;
                                }
                                double x2;
                                ok = map(x1, first, x2);
                                if (!ok)
                                    break;
                                const double bend = x2 - 2.0 * x1 + t2;
                                const double guess = std::abs(bend) > 1e-300 ? t2 - (x1 - t2) * (x1 - t2) / bend : x2;
                                double y;
                                if (std::abs(bend) > 1e-300 && std::isfinite(guess) && map(guess, first, y)
                                        && std::abs(y - guess) < std::abs(x2 - x1))
                                {
                                    t2 = guess;
                                    x1 = y;
                                }
                                else
                                {
                                    // The map of x1 is x2, already found.
                                    t2 = x1;
                                    x1 = x2;
                                }
                            }
                            if (!ok || !converged)
                                continue;
                            // The point's two parameters as the alternation leaves them:
                            // the first edge's from t2, the second's from that.
                            ok = map(t2, t1, t2);
                            if (!ok || t1 < 0.0 || t1 > g1.length || t2 < 0.0 || t2 > g2.length)
                                continue;
                            const Vec3 q1 = g1.start + t1 * g1.direction, q2 = g2.start + t2 * g2.direction;
                            if (q1.z < 0.0 || q2.z < 0.0)
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
                            if (!air(g1, normal_part(g1, in1), false) || !air_or_face(g1, normal_part(g1, mid))
                                    || !air_or_face(g2, normal_part(g2, -1.0 * mid)) || !air(g2, normal_part(g2, out2), false)
                                    || buildings.buried(q1) || buildings.buried(q2))
                                continue;
                            if (!folded_path_open(buildings, q1, source, a.facades, fa == 1, &to_source)
                                    || buildings.blocked(q1, q2)
                                    || !folded_path_open(buildings, q2, microphone, b.facades, fb == 1, &to_microphone))
                                continue;
                            bool leaves_on_face, arrives_on_face;
                            const double theta0_1 = g1.angle(normal_part(g1, in1));
                            const double theta_1 = onto_face(g1.angle(normal_part(g1, mid)), g1.n, leaves_on_face);
                            const double theta0_2 = onto_face(g2.angle(normal_part(g2, -1.0 * mid)), g2.n, arrives_on_face);
                            const bool grazing = leaves_on_face && arrives_on_face;
                            const double theta_2 = g2.angle(normal_part(g2, out2));
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
                            ray.wedge[0] = {wedge_angles(g1.n, theta_1, theta0_1), sb1, s1 * s2 * sb1 / (s1 + s2), false};
                            ray.wedge[1] = {wedge_angles(g2.n, theta_2, theta0_2), sb2, s2 * s3 * sb2 / (s2 + s3), grazing};
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
}

double reflected_amplitude(const PointRay& ray, const std::vector<double>& reflection, double floor_reflection)
{
    double amplitude = ray.spreading;
    for (std::size_t f : ray.facades)
        amplitude *= std::sqrt(reflection[f]);
    for (int i = 0; i < ray.floor; i++)
        amplitude *= std::sqrt(floor_reflection);
    return amplitude;
}

std::complex<double> point_ray_amplitude(const PointRay& ray, const std::vector<double>& reflection,
        double floor_reflection, double k)
{
    return point_ray_amplitude(ray, reflected_amplitude(ray, reflection, floor_reflection), k);
}

std::complex<double> point_ray_amplitude(const PointRay& ray, double reflected, double k)
{
    std::complex<double> amplitude = reflected;
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
        amplitude *= hard_wedge_coefficient(t.angles, k * t.sin_beta, t.distance) / std::sqrt(t.sin_beta)
                * (t.grazing ? 0.5 : 1.0);
    }
    return amplitude;
}

void PointSourceField::band_levels(const std::vector<double>& reflection, double floor_reflection,
        const std::vector<double>& bands, double sound_speed, double* out) const
{
    // The rays in order of length, by index: band_averaged_energy wants them
    // sorted, and a ray carries its path.
    std::vector<std::size_t> order(rays_.size());
    for (std::size_t i = 0; i < order.size(); i++)
        order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return rays_[a].length < rays_[b].length; });
    std::vector<double> length(order.size());
    for (std::size_t i = 0; i < order.size(); i++)
        length[i] = rays_[order[i]].length;
    std::vector<std::complex<double>> amplitude(order.size());
    std::vector<double> band_reflection(facade_count_), reflected(order.size());
    const double half = std::pow(2.0, 1.0 / 6.0);
    for (std::size_t band = 0; band < bands.size(); band++)
    {
        for (std::size_t f = 0; f < facade_count_; f++)
            band_reflection[f] = reflection[f * bands.size() + band];
        // What the reflections leave of each ray is the same at every
        // wavenumber of the band.
        for (std::size_t i = 0; i < order.size(); i++)
            reflected[i] = reflected_amplitude(rays_[order[i]], band_reflection, floor_reflection);
        const double f1 = bands[band] / half, f2 = bands[band] * half;
        double energy = 0.0;
        for (int part = 0; part < 3; part++)
        {
            const double g1 = f1 * std::pow(f2 / f1, part / 3.0), g2 = f1 * std::pow(f2 / f1, (part + 1) / 3.0);
            const double k1 = 2.0 * pi * g1 / sound_speed, k2 = 2.0 * pi * g2 / sound_speed;
            const double kc = 0.5 * (k1 + k2);
            for (std::size_t i = 0; i < order.size(); i++)
                amplitude[i] = point_ray_amplitude(rays_[order[i]], reflected[i], kc);
            energy += band_averaged_energy(amplitude, length, k1, k2) * (g2 - g1);
        }
        out[band] = 10.0 * std::log10(std::max(energy / (f2 - f1), 1e-30));
    }
}

} // namespace niceops
