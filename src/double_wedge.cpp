#include "double_wedge.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

#include "util/constants.hpp"
#include "wedge.hpp"

namespace niceops {

namespace {

using complex = std::complex<double>;

double wrap(double angle)
{
    return angle < 0.0 ? angle + 2.0 * pi : angle;
}

// From the first edge, the line to the second along +x; from the second,
// along -x.  Each wedge's own angles add its face's `drop` (from_face).
double angle1(double x, double z)
{
    return wrap(std::atan2(z, x));
}

double angle2(double x, double z)
{
    return wrap(std::atan2(z, -x));
}

// Buildings::contact: a path within this of an edge grazes past it.
constexpr double contact = 1e-6;

// An angle from the line between the edges turned to one from a wedge's face
// `drop` below it.
double from_face(double angle, double drop)
{
    angle += drop;
    return angle >= 2.0 * pi ? angle - 2.0 * pi : angle;
}

// The end correction's share of the first wedge's Keller coefficient toward
// the second edge: half on a face the wedges share, where the second wedge's
// grazing coefficient doubles the wave that the first's, on its own face,
// has doubled already; all of it where the second edge stands off the first
// wedge's faces.
double keller_scale(const SharedFace& face)
{
    return face.drop1 == 0.0 ? 0.5 : 1.0;
}

// The pair seen from the second edge, for the aperture there.
SharedFace reversed(const SharedFace& face)
{
    return {face.width, face.n2, face.n1, face.drop2, face.drop1};
}

// Whether the aperture at the first edge applies, and its weight: the source
// below the line between the edges and before the first, weighted by its angle
// below the line, the receiver past the aperture.  Where the second wedge's
// face toward the first is tilted
// below the line, the source can see that face's plane; with the receiver in
// the reflection of the first edge's wave off it, the unbounded field the
// aperture takes back out would hold the source's own reflection off the
// face's plane past the second edge, which no ray carries, and the aperture
// at that end does not apply.
double forward_weight(const SharedFace& face, double sx, double sz, double rx, double rz)
{
    if (!(sz < 0.0 && sx < 0.0))
        return 0.0;
    if (face.drop2 > 0.0)
    {
        // The aperture parts the source from the receiver.  On a shared face
        // the receiver is always past it; tilted, a reflection on the
        // receiver's side in the first wedge's face toward the second (a step's
        // wall) puts its image behind the first edge.
        if (!(rx > 0.0))
            return 0.0;
        const double drop = face.drop2;
        const double theta = from_face(angle2(rx - face.width, rz), drop);
        const bool reflected = theta + drop < pi;
        const bool sees_face = sz * std::cos(drop) - (sx - face.width) * std::sin(drop) > 0.0;
        if (reflected && sees_face)
            return 0.0;
    }
    return std::atan2(-sz, -sx);
}

// The same from the second edge's end.
double reverse_weight(const SharedFace& face, double sx, double sz, double rx, double rz)
{
    return forward_weight(reversed(face), face.width - rx, rz, face.width - sx, sz);
}

// The moments' power series below a radian of phase: m_j, the integral of
// v^j exp(i k v) over [0, s], is s^(j+1) times the sum over n of
// (i phase)^n / (n! (n + j + 1)), phase = k s.  Its even terms are the real
// part and its odd ones the imaginary, each a polynomial in phase^2 whose
// coefficients, (-1)^(n/2) / (n! (n + j + 1)), are fixed; to n = 21 they leave
// under 1e-19 of the sum.
struct FilonSeries
{
    double c[3][22];
};

constexpr FilonSeries filon_series()
{
    FilonSeries series {};
    double factorial = 1.0;
    for (int n = 0; n < 22; n++)
    {
        if (n > 0)
            factorial *= n;
        for (int j = 0; j < 3; j++)
            series.c[j][n] = ((n / 2) % 2 == 0 ? 1.0 : -1.0) / (factorial * (n + j + 1));
    }
    return series;
}

constexpr FilonSeries filon_coefficients = filon_series();

// m_j / s^(j+1) for phase under a radian, by Horner in phase^2.
complex filon_moment_series(int j, double phase)
{
    const double u = phase * phase;
    const double* c = filon_coefficients.c[j];
    double even = c[20], odd = c[21];
    for (int n = 18; n >= 0; n -= 2)
    {
        even = even * u + c[n];
        odd = odd * u + c[n + 1];
    }
    return {even, phase * odd};
}

// The integral over [u0, u2] of g(u) exp(i k u), g the quadratic through
// (u0, g0), (u1, g1) and (u2, g2): Filon's rule, exact in the oscillation
// however many wavelengths the panel spans.  Its moments are taken from their
// recurrence where the panel spans a radian or more of phase, and from their
// power series below, where the recurrence cancels.
complex filon_panel(double k, double u0, double u1, double u2, complex g0, complex g1, complex g2)
{
    const double span = u2 - u0, d1 = u1 - u0;
    const complex c1 = (g1 - g0) / d1;
    const complex c2 = ((g2 - g0) / span - c1) / (u2 - u1);
    // m0, m1, m2: the integrals of v^m exp(i k v) over [0, span].
    const double phase = k * span;
    complex m0, m1, m2;
    if (std::abs(phase) < 1.0)
    {
        m0 = span * filon_moment_series(0, phase);
        m1 = span * span * filon_moment_series(1, phase);
        m2 = span * span * span * filon_moment_series(2, phase);
    }
    else
    {
        // Dividing by i k is multiplying by -i / k.
        const complex e = std::polar(1.0, phase), by_ik(0.0, -1.0 / k);
        m0 = (e - 1.0) * by_ik;
        m1 = (span * e - m0) * by_ik;
        m2 = (span * span * e - 2.0 * m1) * by_ik;
    }
    return std::polar(1.0, k * u0) * (g0 * m0 + c1 * m1 + c2 * (m2 - d1 * m1));
}

// One of the aperture integrand's waves at a point of the aperture: its path
// length (the phase per unit k), that length's rate along the aperture, and
// its amplitude.
struct Wave
{
    double length = 0.0;
    double rate = 0.0;
    complex amplitude;
};

// The taper of the aperture's far half at t, for an aperture of `extent`.
double taper(double t, double extent)
{
    if (t < 0.5 * extent)
        return 1.0;
    if (t >= extent)
        return 0.0;
    return 0.5 * (1.0 + std::cos(pi * (t - 0.5 * extent) / (0.5 * extent)));
}

// What the aperture's integrand has at its point t that no wavenumber
// changes: the incident leg `in`, the leg r to the second edge and its
// obliquity factor, the angle theta0 at the second wedge, and for each wave
// its path length, that length's rate along the aperture and the square root
// its amplitude is divided by (`root`).  Only the waves `count` holds are
// filled.  aperture() and aperture_geometry() both take their nodes from
// here, so the CPU's and the GPU's integrands cannot part.
struct NodeGeometry
{
    double in = 0.0, r = 0.0, obliquity = 0.0, theta0 = 0.0;
    double length[4] = {}, rate[4] = {}, root[4] = {};
};

// The aperture's fixed geometry: the second wedge's receiver angle theta and
// distance s3, the direction of its other face, and its face toward the first
// edge's tilt below the line between them.  `above`: the receiver sees the
// first edge past the second; `above_face`: the reflection of the first edge's
// wave off that face reaches it.
struct ApertureFrame
{
    double w = 0.0, sx = 0.0, sz = 0.0, rx = 0.0, rz = 0.0;
    double s1 = 0.0, theta = 0.0, s3 = 0.0, fx = 0.0, fz = 0.0;
    double drop = 0.0, cos_drop = 1.0, sin_drop = 0.0;
    bool above = false, above_face = false;
};

ApertureFrame aperture_frame(const SharedFace& face, double sx, double sz, double rx, double rz)
{
    ApertureFrame f;
    f.w = face.width;
    f.sx = sx;
    f.sz = sz;
    f.rx = rx;
    f.rz = rz;
    f.s1 = std::hypot(sx, sz);
    f.drop = face.drop2;
    f.cos_drop = std::cos(face.drop2);
    f.sin_drop = std::sin(face.drop2);
    f.theta = from_face(angle2(rx - f.w, rz), f.drop);
    f.s3 = std::hypot(rx - f.w, rz);
    // The second wedge's other face, for the reflection off it.
    f.fx = -std::cos(face.n2 * pi - f.drop);
    f.fz = std::sin(face.n2 * pi - f.drop);
    // On a shared face the receiver sees the first edge, and the face's
    // reflection of its wave reaches it, above the face's plane.  Tilted, the
    // line between the edges touches both wedges from above, so the receiver
    // sees the first edge on or above that line -- on it within the path
    // search's contact tolerance, which traces the first edge's ray grazing
    // past the second there -- and the face's reflection reaches it where
    // hard_wedge_field's test says, the first edge at drop from the face.
    f.above = face.drop2 == 0.0 ? rz > 0.0 : rz > -contact;
    f.above_face = face.drop2 == 0.0 ? rz > 0.0 : f.theta + f.drop < pi;
    return f;
}

NodeGeometry node_geometry(const ApertureFrame& f, double t, const int* count)
{
    NodeGeometry n;
    n.in = std::hypot(f.sx, t - f.sz);
    n.r = std::hypot(f.w, t);
    n.obliquity = 0.5 * (-f.sx / n.in + f.w / n.r);
    const double rate = (t - f.sz) / n.in;
    n.theta0 = std::atan2(t, f.w) + f.drop;
    n.length[0] = n.in + n.r + f.s3;
    n.rate[0] = rate + t / n.r;
    n.root[0] = std::sqrt(n.r * f.s3);
    if (count[1])
    {
        const double d = std::hypot(f.rx, f.rz - t);
        n.length[1] = n.in + d;
        n.rate[1] = rate + (t - f.rz) / d;
        n.root[1] = std::sqrt(d);
    }
    if (count[2])
    {
        // The aperture point mirrored in the plane of the second wedge's face
        // toward the first edge: through (w, 0), tilted `drop` below the line.
        // Along the aperture the image moves by (2 sin cos, 2 sin^2 - 1).
        const double along = f.w * f.cos_drop - t * f.sin_drop;
        const double ix = 2.0 * f.w - 2.0 * along * f.cos_drop, iz = -2.0 * along * f.sin_drop - t;
        const double d = std::hypot(f.rx - ix, f.rz - iz);
        n.length[2] = n.in + d;
        n.rate[2] = rate
                - ((f.rx - ix) * (2.0 * f.sin_drop * f.cos_drop) + (f.rz - iz) * (2.0 * f.sin_drop * f.sin_drop - 1.0))
                        / d;
        n.root[2] = std::sqrt(d);
    }
    if (count[3])
    {
        // The aperture point mirrored in the other face's plane.
        const double along = -f.w * f.fx + t * f.fz;
        const double ix = f.w + 2.0 * along * f.fx + f.w, iz = 2.0 * along * f.fz - t;
        const double d = std::hypot(f.rx - ix, f.rz - iz);
        const double rate_d = -((f.rx - ix) * 2.0 * f.fz * f.fx + (f.rz - iz) * (2.0 * f.fz * f.fz - 1.0)) / d;
        n.length[3] = n.in + d;
        n.rate[3] = rate + rate_d;
        n.root[3] = std::sqrt(d);
    }
    return n;
}

// The cuts, where theta0 crosses one of hard_wedge_field's boundaries, out to
// `extent`, added to `cuts` and `features`.
void boundary_cuts(const SharedFace& face, const ApertureFrame& f, double extent, std::vector<double>& cuts,
        std::vector<double>& features)
{
    const double theta = f.theta;
    for (double boundary : {theta - pi, pi - theta, (face.n2 - 1.0) * pi, (2.0 * face.n2 - 1.0) * pi - theta})
    {
        // theta0 = atan(t / w) + drop.
        const double off_line = boundary - f.drop;
        if (off_line > 0.0 && off_line < 0.5 * pi && f.w * std::tan(off_line) < extent)
        {
            cuts.push_back(f.w * std::tan(off_line));
            features.push_back(cuts.back());
        }
    }
}

// Which waves reach the receiver between the cuts a and b (hard_wedge_field's
// tests), each +1, -1 or 0: where the receiver sees the first edge the
// incident wave counts one less, and where the second wedge's face toward it
// reflects the first edge's wave to the receiver the reflection off that face
// does (aperture() has why).
void wave_counts(const SharedFace& face, const ApertureFrame& f, double a, double b, int* count)
{
    const double middle = std::atan2(0.5 * (a + b), f.w) + f.drop;
    const double theta = f.theta;
    count[0] = 1;
    count[1] = (std::abs(theta - middle) < pi ? 1 : 0) - (f.above ? 1 : 0);
    count[2] = (theta < pi - middle ? 1 : 0) - (f.above_face ? 1 : 0);
    count[3] = middle > (face.n2 - 1.0) * pi && theta > (2.0 * face.n2 - 1.0) * pi - middle ? 1 : 0;
}

// A panel's length from t: a sixth to a seventh of the distance from the
// nearest feature, and at most the distance from the second edge.
double panel_step(const std::vector<double>& features, double t, double smallest, double w)
{
    double near = std::abs(t - features[0]);
    for (double f : features)
        near = std::min(near, std::abs(t - f));
    return std::min(0.15 * (near + smallest), t + w);
}

// What an aperture's node at t has that no wavenumber changes: node_geometry's,
// and the second wedge's angles and distance parameter there.
struct NodeAt
{
    NodeGeometry geometry;
    WedgeAngles angles;
    double distance = 0.0;
};

NodeAt node_at(const ApertureFrame& f, const WedgeSweep& sweep, double t, const int* count)
{
    NodeAt n;
    n.geometry = node_geometry(f, t, count);
    n.angles = sweep.at(n.geometry.theta0);
    n.distance = n.geometry.r * f.s3 / (n.geometry.r + f.s3);
    return n;
}

// node_at's results by t, for one aperture taken at many wavenumbers.  Each
// wavenumber's nodes are graded from the same edge and cuts, so most of them
// land, to the bit, on points another has taken already; only those near its
// taper's start and its extent are its own.  A node's wave counts follow from
// the cuts that bound it, so they agree wherever its t does; they are compared
// anyway, a mismatch taking the node afresh.  Open addressing on t's bits,
// emptied by a generation count rather than cleared.
class NodeTable
{
public:
    void clear()
    {
        if (++generation_ == 0)
        {
            for (Slot& s : slots_)
                s.generation = 0;
            generation_ = 1;
        }
        used_ = 0;
        if (slots_.empty())
            slots_.resize(1024);
    }

    // The node at t with these wave counts: the table's, or make()'s, kept.
    template <typename Make>
    const NodeAt& find(double t, const int* count, Make make)
    {
        if (2 * (used_ + 1) > slots_.size())
            grow();
        Slot& s = probe(t);
        if (s.generation == generation_ && std::equal(count, count + 4, s.count))
            return s.node;
        if (s.generation != generation_)
            used_++;
        s.generation = generation_;
        s.t = t;
        std::copy(count, count + 4, s.count);
        s.node = make();
        return s.node;
    }

private:
    struct Slot
    {
        unsigned generation = 0;
        int count[4] = {0, 0, 0, 0};
        double t = 0.0;
        NodeAt node;
    };

    Slot& probe(double t)
    {
        std::uint64_t bits;
        std::memcpy(&bits, &t, sizeof bits);
        const std::size_t mask = slots_.size() - 1;
        for (std::size_t i = (bits * 0x9E3779B97F4A7C15ull) >> 40 & mask;; i = (i + 1) & mask)
            if (slots_[i].generation != generation_ || slots_[i].t == t)
                return slots_[i];
    }

    void grow()
    {
        std::vector<Slot> old(2 * slots_.size());
        old.swap(slots_);
        for (Slot& s : old)
            if (s.generation == generation_)
                probe(s.t) = s;
    }

    std::vector<Slot> slots_;
    std::size_t used_ = 0;
    unsigned generation_ = 0;
};

// The source below the line between the edges and before the first.  The aperture
// is the half-line x = 0, z > 0, out to 30 wavelengths or three Fresnel zones,
// its far half tapered.  Each of its points sends a line source's wave to the
// second wedge, whose field -- the geometrical-optics waves past it and its
// diffraction -- the receiver sums.
//
// The integrand is four waves, each a smooth amplitude times exp(i k l), l a
// path length that grows along the aperture: the second wedge's diffracted
// wave, the incident wave, and its reflections off the wedge's face toward
// the first edge and off its other face.  Each geometrical-optics wave reaches the receiver
// or not, and the coefficient jumps, only where the second wedge's shadow and
// reflection boundaries cross the aperture, so it is cut there and between the
// cuts each wave is integrated by Filon's rule in its own path length, the
// amplitude per unit length quadratic over each panel.  The panels grow
// geometrically away from the edge and the cuts, where the amplitude changes
// fastest.  Against the midpoint rule at a 192nd of a wavelength this is
// within 0.016 dB (0.0015 dB rms) over the shared faces of DigiCity and of a
// district; the midpoint rule at a 24th of a wavelength, which it replaced,
// was within 0.027 dB (0.014-0.018 dB rms), for 6 to 24 times the evaluations.
//
// What does not depend on the wavenumber is worked out on construction; at()
// takes its nodes' geometry from `nodes(t, count)`, which either works it out
// (aperture()) or finds it in a NodeTable (aperture_spectrum()), the same
// arithmetic either way.
class Aperture
{
public:
    Aperture(const SharedFace& face, double sx, double sz, double rx, double rz)
        : face_(face), f_(aperture_frame(face, sx, sz, rx, rz)), sweep_(face.n2, f_.theta)
    {
        sin_g_ = -sz / f_.s1;
        cos_g_ = -sx / f_.s1;
        if (sin_g_ > 1e-6)
            keller_ = wedge_angles(face.n1, face.drop1, from_face(angle1(sx, sz), face.drop1));
        keller_scale_ = keller_scale(face);
    }

    NodeAt node(double t, const int* count) const
    {
        return node_at(f_, sweep_, t, count);
    }

    template <typename Nodes>
    complex at(double k, Nodes&& nodes) const
    {
        const ApertureFrame& f = f_;
        const double w = f.w, s1 = f.s1;
        const double extent = aperture_extent(std::max(s1, w), k);
        const complex launch = std::sqrt(k / (2.0 * pi)) * std::polar(1.0, -0.25 * pi);
        // Above the line between the edges the receiver sees the first edge,
        // and the aperture's waves past the second edge hold that wedge's whole
        // field -- the source's own wave where the wedge does not shadow it, and
        // the edge's diffracted wave with the face's reflection of it -- which
        // the rays of the first edge carry already.  The same aperture with the
        // face going on past the second edge is that field alone, in the
        // aperture's own approximation, and is taken back out: what remains is
        // what the second edge does to it.  So there the incident wave is
        // counted less one, the unbounded field's, and where the face's
        // reflection of the first edge's wave reaches the receiver, that
        // reflection is too (wave_counts).

        // The cuts: where theta0, the aperture point's angle at the second wedge,
        // crosses one of hard_wedge_field's boundaries.  The taper's start is a
        // panel end too, though the amplitude changes smoothly there.
        std::vector<double> cuts {0.0, 0.5 * extent, extent};
        std::vector<double> features {0.0};
        boundary_cuts(face_, f, extent, cuts, features);
        std::sort(cuts.begin(), cuts.end());

        // The waves at t, those of `count` (+1, -1 or 0) only.  The second
        // wedge's coefficient is wanted at one theta and many theta0.
        auto waves = [&](double t, const int* count, Wave* out) {
            const NodeAt& at = nodes(t, count);
            const NodeGeometry& n = at.geometry;
            const complex amplitude = taper(t, extent) * n.obliquity * launch / std::sqrt(n.in);
            out[0].length = n.length[0];
            out[0].rate = n.rate[0];
            out[0].amplitude = amplitude * hard_wedge_coefficient(at.angles, k, at.distance) / n.root[0];
            for (int j = 1; j < 4; j++)
                if (count[j])
                    out[j] = {n.length[j], n.rate[j], double(count[j]) * amplitude / n.root[j]};
        };

        const double smallest = std::max(0.01 * std::min(w, s1), 1e-6 * extent);
        complex sum(0.0, 0.0);
        std::vector<double> ends;
        for (std::size_t c = 0; c + 1 < cuts.size(); c++)
        {
            const double a = cuts[c], b = cuts[c + 1];
            if (b - a < 1e-9 * extent)
                continue;
            int count[4];
            wave_counts(face_, f, a, b, count);
            ends.assign(1, a);
            for (double t = a; t < b;)
            {
                const double h = panel_step(features, t, smallest, w);
                double next = t + h;
                if (next > b || b - next < 0.3 * h)
                    next = b;
                ends.push_back(next);
                t = next;
            }
            // Each segment's ends are taken a hair inside it, on the side of a cut
            // the segment's waves belong to.
            const double inside = 1e-7 * (b - a);
            Wave first[4], middle_wave[4], last[4];
            waves(a + inside, count, first);
            for (std::size_t p = 0; p + 1 < ends.size(); p++)
            {
                waves(0.5 * (ends[p] + ends[p + 1]), count, middle_wave);
                waves(p + 2 == ends.size() ? b - inside : ends[p + 1], count, last);
                for (int j = 0; j < 4; j++)
                    if (count[j])
                    {
                        // Per unit path length the amplitude is A / (dl/dt).  The
                        // lengths grow along the aperture wherever a wave counts:
                        // its stationary point, the straight path through the
                        // aperture, lies where the second wedge shadows it, or
                        // below the line between the edges.
                        sum += filon_panel(k, first[j].length, middle_wave[j].length, last[j].length,
                                first[j].amplitude / first[j].rate, middle_wave[j].amplitude / middle_wave[j].rate,
                                last[j].amplitude / last[j].rate);
                    }
                std::copy(last, last + 4, first);
            }
        }
        // The end point at the edge gives the aperture's own coefficient,
        // (1 + cos g) / 2 exp(i pi/4) / (sqrt(2 pi k) sin g), g the source's angle
        // below the line between the edges; the wedge's is its Keller
        // coefficient toward the second edge, half of it on a shared face, the
        // face doubling it (keller_scale).  Both grow as 1/g, and their ratio
        // tends to 1.
        if (sin_g_ > 1e-6)
        {
            const complex kirchhoff
                    = 0.5 * (1.0 + cos_g_) * std::polar(1.0, 0.25 * pi) / (std::sqrt(2.0 * pi * k) * sin_g_);
            const complex keller = hard_wedge_coefficient(keller_, k, 1e9);
            sum *= keller_scale_ * keller / kirchhoff;
        }
        return sum;
    }

private:
    SharedFace face_;
    ApertureFrame f_;
    WedgeSweep sweep_;
    double sin_g_ = 0.0, cos_g_ = 0.0;
    WedgeAngles keller_;
    double keller_scale_ = 0.5;
};

// Where the receiver is above the line between the edges and the
// aperture's source within about 1e-5 rad below it, the aperture's end-point
// terms, each growing as 1/g, no longer cancel in double precision: the field
// strays by decibels and then grows without bound, to 6-13 dB over the direct
// wave on a suburban street.  Above about 1e-4 rad it has reached its limit
// as g goes to zero (within 0.2 dB from 1e-3 rad), so a source nearer the line
// than nearest_angle is taken at that angle, at its own distance.  Its weight
// stays its own angle.
constexpr double nearest_angle = 1e-3;

void off_the_line(double& sx, double& sz)
{
    const double s1 = std::hypot(sx, sz);
    if (-sz < s1 * std::sin(nearest_angle))
    {
        sx = -s1 * std::cos(nearest_angle);
        sz = -s1 * std::sin(nearest_angle);
    }
}

complex aperture(const SharedFace& face, double sx, double sz, double rx, double rz, double k)
{
    off_the_line(sx, sz);
    const Aperture a(face, sx, sz, rx, rz);
    return a.at(k, [&](double t, const int* count) { return a.node(t, count); });
}

// aperture() at each of the wavenumbers, into out, its nodes' geometry shared
// through `table`.
void aperture_spectrum(const SharedFace& face, double sx, double sz, double rx, double rz,
        const std::vector<double>& wavenumbers, NodeTable& table, complex* out)
{
    off_the_line(sx, sz);
    const Aperture a(face, sx, sz, rx, rz);
    table.clear();
    auto nodes = [&](double t, const int* count) -> const NodeAt& {
        return table.find(t, count, [&] { return a.node(t, count); });
    };
    for (std::size_t j = 0; j < wavenumbers.size(); j++)
        out[j] = a.at(wavenumbers[j], nodes);
}

// aperture()'s integrand and grading, worked out once for the wavenumbers
// given (none of them zero).
ApertureGeometry aperture_geometry(const SharedFace& face, double sx, double sz, double rx, double rz,
        const std::vector<double>& wavenumbers)
{
    off_the_line(sx, sz);
    ApertureGeometry g;
    const ApertureFrame f = aperture_frame(face, sx, sz, rx, rz);
    const double w = f.w, s1 = f.s1;
    g.extent_scale = std::max(s1, w);
    const double extent = aperture_extent(g.extent_scale, *std::min_element(wavenumbers.begin(), wavenumbers.end()));
    const double shortest = aperture_extent(g.extent_scale, *std::max_element(wavenumbers.begin(), wavenumbers.end()));
    g.reference = s1 + w + f.s3;
    std::vector<double> cuts {0.0, extent};
    std::vector<double> features {0.0};
    boundary_cuts(face, f, extent, cuts, features);
    // Each wavenumber's taper start and extent end a panel; they are not
    // features, which would grade the panels around them.
    std::vector<double> ends_fixed;
    for (double k : wavenumbers)
        for (double e : {0.5 * aperture_extent(g.extent_scale, k), aperture_extent(g.extent_scale, k)})
            if (e > 0.0 && e < extent)
                ends_fixed.push_back(e);
    std::sort(cuts.begin(), cuts.end());
    std::sort(ends_fixed.begin(), ends_fixed.end());
    const WedgeSweep sweep(face.n2, f.theta);
    auto node = [&](double t, const int* count) {
        const NodeGeometry geometry = node_geometry(f, t, count);
        ApertureNode n;
        n.t = t;
        const double base = geometry.obliquity / std::sqrt(geometry.in);
        for (int j = 0; j < 4; j++)
            if (count[j])
            {
                n.length[j] = geometry.length[j] - g.reference;
                n.rate[j] = geometry.rate[j];
                n.weight[j] = count[j] * base / geometry.root[j];
            }
        n.angles = sweep.at(geometry.theta0);
        n.distance = geometry.r * f.s3 / (geometry.r + f.s3);
        g.nodes.push_back(n);
    };
    const double smallest = std::max(0.01 * std::min(w, s1), 1e-6 * extent);
    for (std::size_t c = 0; c + 1 < cuts.size(); c++)
    {
        const double a = cuts[c], b = cuts[c + 1];
        if (b - a < 1e-9 * extent)
            continue;
        int count[4];
        wave_counts(face, f, a, b, count);
        std::vector<double> ends {a};
        for (double t = a; t < b;)
        {
            const double h = panel_step(features, t, smallest, w);
            double next = t + h;
            if (next > b || b - next < 0.3 * h)
                next = b;
            // A taper's start or an extent before the step's end ends the
            // panel there instead.
            const auto fixed = std::upper_bound(ends_fixed.begin(), ends_fixed.end(), t + 1e-9 * (b - a));
            if (fixed != ends_fixed.end() && *fixed < next - 1e-9 * (b - a))
                next = *fixed;
            ends.push_back(next);
            t = next;
        }
        // A hair inside, as aperture() takes it, but never more than a hair of
        // the shortest aperture: a geometrical-optics wave's stationary point
        // is on its cut, where its amplitude per unit length goes as
        // 1 / sqrt(u - u_c), and how much of that the first panel misses
        // depends on how far in the node sits.  1e-7 of a 200 m segment moved
        // a 10 kHz field by 0.06 dB.
        const double inside = 1e-7 * std::min(b - a, shortest);
        g.segments.push_back(g.nodes.size());
        node(a + inside, count);
        for (std::size_t p = 0; p + 1 < ends.size(); p++)
        {
            node(0.5 * (ends[p] + ends[p + 1]), count);
            node(p + 2 == ends.size() ? b - inside : ends[p + 1], count);
        }
    }
    g.segments.push_back(g.nodes.size());
    g.sin_g = -sz / s1;
    g.cos_g = -sx / s1;
    g.keller = wedge_angles(face.n1, face.drop1, from_face(angle1(sx, sz), face.drop1));
    g.keller_scale = keller_scale(face);
    return g;
}

complex aperture_at(const ApertureGeometry& g, double k)
{
    const double extent = aperture_extent(g.extent_scale, k);
    complex sum(0.0, 0.0);
    for (std::size_t s = 0; s + 1 < g.segments.size(); s++)
        for (std::size_t i = g.segments[s]; i + 2 < g.segments[s + 1]; i += 2)
        {
            const ApertureNode* n = &g.nodes[i];
            if (n[0].t >= extent)
                break;
            complex d2[3];
            for (int m = 0; m < 3; m++)
                d2[m] = hard_wedge_coefficient(n[m].angles, k, n[m].distance);
            for (int j = 0; j < 4; j++)
            {
                if (n[0].weight[j] == 0.0 && n[2].weight[j] == 0.0)
                    continue;
                complex v[3];
                for (int m = 0; m < 3; m++)
                    v[m] = taper(n[m].t, extent) * n[m].weight[j] / n[m].rate[j] * (j == 0 ? d2[m] : complex(1.0, 0.0));
                sum += filon_panel(k, n[0].length[j], n[1].length[j], n[2].length[j], v[0], v[1], v[2]);
            }
        }
    sum *= std::sqrt(k / (2.0 * pi)) * std::polar(1.0, -0.25 * pi);
    if (g.sin_g > 1e-6)
    {
        const complex kirchhoff = 0.5 * (1.0 + g.cos_g) * std::polar(1.0, 0.25 * pi) / (std::sqrt(2.0 * pi * k) * g.sin_g);
        sum *= g.keller_scale * hard_wedge_coefficient(g.keller, k, 1e9) / kirchhoff;
    }
    return sum;
}

} // namespace

bool shared_face_applies(const SharedFace& face, double sx, double sz, double rx, double rz)
{
    return forward_weight(face, sx, sz, rx, rz) > 0.0 || reverse_weight(face, sx, sz, rx, rz) > 0.0;
}

std::complex<double> shared_face_field(const SharedFace& face, double sx, double sz, double rx, double rz, double k)
{
    // Each end below the line between the edges can be the aperture's source;
    // each evaluation is weighted by its source's angle below the line, so that
    // the one whose source nears it gives way to the other smoothly, and the
    // weights swap with the ends.
    const double forward = forward_weight(face, sx, sz, rx, rz);
    const double reverse = reverse_weight(face, sx, sz, rx, rz);
    if (forward + reverse <= 0.0)
        return {std::numeric_limits<double>::quiet_NaN(), 0.0};
    complex sum(0.0, 0.0);
    if (forward > 0.0)
        sum += forward * aperture(face, sx, sz, rx, rz, k);
    if (reverse > 0.0)
        sum += reverse * aperture(reversed(face), face.width - rx, rz, face.width - sx, sz, k);
    return sum / (forward + reverse);
}

std::vector<std::complex<double>> shared_face_spectrum(const SharedFace& face, double sx, double sz, double rx,
        double rz, const std::vector<double>& wavenumbers)
{
    // As shared_face_field, wavenumber by wavenumber.
    const double forward = forward_weight(face, sx, sz, rx, rz);
    const double reverse = reverse_weight(face, sx, sz, rx, rz);
    std::vector<complex> out(wavenumbers.size(), {std::numeric_limits<double>::quiet_NaN(), 0.0});
    if (forward + reverse <= 0.0)
        return out;
    thread_local NodeTable table_scratch;
    thread_local std::vector<complex> forward_scratch, reverse_scratch;
    NodeTable& table = table_scratch;
    std::vector<complex>& forward_field = forward_scratch;
    std::vector<complex>& reverse_field = reverse_scratch;
    forward_field.resize(wavenumbers.size());
    reverse_field.resize(wavenumbers.size());
    if (forward > 0.0)
        aperture_spectrum(face, sx, sz, rx, rz, wavenumbers, table, forward_field.data());
    if (reverse > 0.0)
        aperture_spectrum(reversed(face), face.width - rx, rz, face.width - sx, sz, wavenumbers, table,
                reverse_field.data());
    for (std::size_t j = 0; j < wavenumbers.size(); j++)
    {
        complex sum(0.0, 0.0);
        if (forward > 0.0)
            sum += forward * forward_field[j];
        if (reverse > 0.0)
            sum += reverse * reverse_field[j];
        out[j] = sum / (forward + reverse);
    }
    return out;
}

double aperture_extent(double extent_scale, double k)
{
    const double lambda = 2.0 * pi / k;
    return std::max(3.0 * std::sqrt(lambda * extent_scale), 30.0 * lambda);
}

SharedFaceGeometry shared_face_geometry(const SharedFace& face, double sx, double sz, double rx, double rz,
        const std::vector<double>& wavenumbers)
{
    if (wavenumbers.empty() || *std::min_element(wavenumbers.begin(), wavenumbers.end()) <= 0.0)
        throw std::invalid_argument("shared_face_geometry: the wavenumbers must be positive, and there must be some");
    SharedFaceGeometry g;
    g.forward_weight = forward_weight(face, sx, sz, rx, rz);
    g.reverse_weight = reverse_weight(face, sx, sz, rx, rz);
    if (g.forward_weight > 0.0)
        g.forward = aperture_geometry(face, sx, sz, rx, rz, wavenumbers);
    if (g.reverse_weight > 0.0)
        g.reverse = aperture_geometry(reversed(face), face.width - rx, rz, face.width - sx, sz, wavenumbers);
    return g;
}

std::complex<double> shared_face_field_at(const SharedFaceGeometry& g, double k)
{
    if (g.forward_weight + g.reverse_weight <= 0.0)
        return {std::numeric_limits<double>::quiet_NaN(), 0.0};
    complex sum(0.0, 0.0);
    if (g.forward_weight > 0.0)
        sum += g.forward_weight * aperture_at(g.forward, k);
    if (g.reverse_weight > 0.0)
        sum += g.reverse_weight * aperture_at(g.reverse, k);
    return sum / (g.forward_weight + g.reverse_weight);
}

} // namespace niceops
