#include "stratified.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>

#include "airy.hpp"
#include "util/constants.hpp"
#include "util/parallel.hpp"

namespace niceops {

namespace {

constexpr double nan = std::numeric_limits<double>::quiet_NaN();

// sin of the ray's angle to the horizontal where the sound speed is c, for
// Snell invariant p: 0 at a turning point.
inline double sine(double p, double c)
{
    const double v = 1.0 - p * p * c * c;
    return v > 0.0 ? std::sqrt(v) : 0.0;
}

struct Crossing
{
    double dx = 0.0; // horizontal distance
    double dt = 0.0; // travel time
    double ds = 0.0; // arc length
};

// What crossing a layer from height z0 (speed c0) to z1 (speed c1), c linear
// between, costs a ray with invariant p that does not turn inside it.  Written
// so that neither a vanishing gradient nor a vertical ray divides by zero: with
// q the sines at the two ends,
//   dx = p (c0 + c1) dz / (q0 + q1),
//   dt = log1p(g dz K / ((1 + q1) c0)) / g,
//   ds = atan2(p g dz A, q0 q1 + p^2 c0 c1) / (g p),
// where K and A are the factors that remain once the difference of the two
// ends' q is written as p^2 (c1^2 - c0^2) / (q0 + q1).  Each tends to the
// straight-ray value as g -> 0 and to the vertical one as p -> 0.
//
// At a turning point the sine is zero by definition; evaluated from p and c it
// is instead the square root of a rounding error, ~1e-8, which costs the
// layer's crossing that much accuracy, so the caller passes q = 0 there.  The
// sines q0 and q1 come in formed (Path::between forms each level's once).
template <bool full>
inline Crossing cross_sines(double p, double z0, double c0, double q0, double z1, double c1, double q1)
{
    Crossing out;
    const double dz = z1 - z0;
    if (!(dz > 0.0))
        return out;
    const double sum = q0 + q1;
    out.dx = sum > 0.0 ? p * (c0 + c1) * dz / sum : 0.0;
    if constexpr (full)
    {
        const double g = (c1 - c0) / dz;
        const double pp = p * p * (c0 + c1) / (sum > 0.0 ? sum : 1.0);
        // Time.
        const double K = 1.0 + q0 + c0 * pp;
        const double base = dz * K / ((1.0 + q1) * c0);
        const double x = g * base;
        out.dt = std::fabs(x) < 1e-8 ? base * (1.0 - 0.5 * x) : base * std::log1p(x) / x;
        // Arc.
        const double A = q0 + c0 * pp;
        const double D = q0 * q1 + p * p * c0 * c1;
        const double y = p * g * dz * A;
        out.ds = std::fabs(y) < 1e-8 * std::fabs(D) ? dz * A / D : std::atan2(y, D) / (g * p);
    }
    return out;
}


// One path's layered speed, with what a ray needs precomputed.
struct Path
{
    const Layers& L;

    // A pair's solve walks between the same few heights -- the source, the
    // receiver, the ground -- for every p it tries: each height's first level
    // above it and its speed, looked up once (`where`) and kept.
    struct Height
    {
        double z = std::numeric_limits<double>::quiet_NaN();
        const double* upper = nullptr;
        double c = 0.0;
    };
    mutable Height remembered[4];
    mutable unsigned next = 0;

    const Height& where(double z) const
    {
        for (const Height& h : remembered)
            if (h.z == z)
                return h;
        Height& h = remembered[next++ % 4];
        h.z = z;
        h.upper = std::upper_bound(L.z, L.z + L.levels, z);
        h.c = at_below(z, h.upper);
        return h;
    }

    // Layers::at(height) given `upper`, the first level above `height`
    // (std::upper_bound's answer), so that a walk that has found it does not
    // search again: the same branches and the same arithmetic.
    double at_below(double height, const double* upper) const
    {
        if (L.levels == 1 || height <= L.z[0])
            return L.c[0];
        if (height >= L.z[L.levels - 1])
            return L.c[L.levels - 1];
        const std::size_t i = static_cast<std::size_t>(upper - L.z) - 1;
        return L.c[i] + (height - L.z[i]) / (L.z[i + 1] - L.z[i]) * (L.c[i + 1] - L.c[i]);
    }

    // Sum the crossings of every layer between heights za < zb; `turn_a` /
    // `turn_b` say the ray turns at that end.
    template <bool full>
    Crossing between(double p, double za, double zb, bool turn_a = false, bool turn_b = false) const
    {
        // Each level's sine is formed once, as the top of one layer and the
        // bottom of the next.  A layer of zero thickness contributes nothing,
        // whatever its sines.
        Crossing total;
        if (!(zb > za))
            return total;
        const double* begin = L.z;
        const double* end = L.z + L.levels;
        const Height& a = where(za);
        const double* level = a.upper;
        double z0 = za;
        double c0 = a.c;
        double q0 = turn_a ? 0.0 : sine(p, c0);
        for (; level != end && *level < zb; ++level)
        {
            const double c1 = L.c[level - begin];
            const double q1 = sine(p, c1);
            const auto piece = cross_sines<full>(p, z0, c0, q0, *level, c1, q1);
            total.dx += piece.dx;
            total.dt += piece.dt;
            total.ds += piece.ds;
            z0 = *level;
            c0 = c1;
            q0 = q1;
        }
        const double cb = where(zb).c;
        const auto piece = cross_sines<full>(p, z0, c0, q0, zb, cb, turn_b ? 0.0 : sine(p, cb));
        total.dx += piece.dx;
        total.dt += piece.dt;
        total.ds += piece.ds;
        return total;
    }

    // The largest sound speed on [za, zb]: at the ends or at a level between.
    double fastest(double za, double zb) const
    {
        double best = std::max(L.at(za), L.at(zb));
        for (std::size_t k = 0; k < L.levels; k++)
        {
            if (L.z[k] > za && L.z[k] < zb)
                best = std::max(best, L.c[k]);
        }
        return best;
    }

    // The largest sound speed anywhere at or above z (the profile is constant
    // above its top level).
    double fastest_above(double z) const
    {
        double best = L.at(z);
        for (std::size_t k = 0; k < L.levels; k++)
        {
            if (L.z[k] > z)
                best = std::max(best, L.c[k]);
        }
        return best;
    }

    // Where a ray with invariant p, going down from z, first turns: the
    // height at which c reaches 1/p.  NaN if it gets to the ground first.
    double turning_below(double p, double z) const
    {
        const double target = 1.0 / p;
        double z1 = z;
        double c1 = L.at(z);
        if (c1 >= target)
            return z;
        const double* begin = L.z;
        const double* end = L.z + L.levels;
        const double* level = std::lower_bound(begin, end, z);
        while (level != begin)
        {
            --level;
            if (*level < 0.0)
                break;
            const double c0 = L.c[level - begin];
            if (c0 >= target)
                return z1 - (target - c1) / (c0 - c1) * (z1 - *level);
            z1 = *level;
            c1 = c0;
        }
        const double c0 = L.at(0.0);
        if (c0 >= target && z1 > 0.0)
            return z1 - (target - c1) / (c0 - c1) * z1;
        return nan;
    }

    // Where a ray with invariant p, going up from z, first turns: the height
    // at which c reaches 1/p.  NaN if it never does.
    double turning_above(double p, double z) const
    {
        const double target = 1.0 / p;
        const double* begin = L.z;
        const double* end = L.z + L.levels;
        const Height& h = where(z);
        const double* level = h.upper;
        double z0 = z;
        double c0 = h.c;
        if (c0 >= target)
            return z;
        for (; level != end; ++level)
        {
            const double c1 = L.c[level - begin];
            if (c1 >= target)
                return z0 + (target - c0) / (c1 - c0) * (*level - z0);
            z0 = *level;
            c0 = c1;
        }
        return nan;
    }
};

// The invariants at which a ray's range jumps.  A ray going up from z_s turns
// at the first height where c reaches 1/p, so as 1/p passes the speed at a
// local maximum of c that is the fastest air yet met above z_s, the turning
// point leaps from just under that maximum to wherever the air is next as
// fast -- and the range with it.  The same holds going down from z_r for the
// rays that dip.  A family is continuous between these, not across them, and
// a bracket straddling one can change sign with no root inside.
//
// Listed ascending.  Only a performance matter: an extra entry costs solve()
// two evaluations, a missing one (past the capacity) leaves solve() to find
// out the slow way.
struct Jumps
{
    static constexpr int capacity = 32;
    double p[capacity];
    int count = 0;

    Jumps(const Layers& L, double z_s, double z_r)
    {
        const double* begin = L.z;
        const double* end = L.z + L.levels;
        double fastest = L.at(z_s);
        for (const double* level = std::upper_bound(begin, end, z_s); level + 1 < end; ++level)
        {
            const std::size_t k = level - begin;
            if (L.c[k] > fastest)
            {
                fastest = L.c[k];
                if (L.c[k + 1] <= L.c[k])
                    add(1.0 / L.c[k]);
            }
        }
        fastest = L.at(z_r);
        for (const double* level = std::lower_bound(begin, end, z_r); level != begin;)
        {
            --level;
            if (*level <= 0.0)
                break;
            const std::size_t k = level - begin;
            if (L.c[k] > fastest)
            {
                fastest = L.c[k];
                if ((level == begin ? L.at(0.0) : L.c[k - 1]) <= L.c[k])
                    add(1.0 / L.c[k]);
            }
        }
        std::sort(p, p + count);
    }

    void add(double value)
    {
        if (count < capacity)
            p[count++] = value;
    }
};

// A ray with Snell invariant p from z_s: launched up (turning aloft first) or
// down, and arriving at z_r either on its first crossing, from above, or on
// its second -- after dipping below z_r and turning back up without touching
// the ground.  With `reflected` it bounces off z = 0 instead and rises to z_r.
// NaN range when the ray never gets there.
struct Arm
{
    bool up = false;
    bool second = false;
};

template <bool full>
Crossing trace_p(const Path& path, double z_s, double z_r, double p, Arm arm, bool reflected)
{
    Crossing out;
    auto add = [&out](const Crossing& c, double times) {
        out.dx += times * c.dx;
        out.dt += times * c.dt;
        out.ds += times * c.ds;
    };
    if (arm.up)
    {
        const double top = path.turning_above(p, z_s);
        if (std::isnan(top))
        {
            out.dx = nan;
            return out;
        }
        add(path.between<full>(p, z_s, top, false, true), 2.0);
    }
    if (reflected)
    {
        add(path.between<full>(p, 0.0, z_s), 1.0);
        add(path.between<full>(p, 0.0, z_r), 1.0);
        return out;
    }
    add(path.between<full>(p, z_r, z_s), 1.0);
    if (arm.second)
    {
        const double bottom = path.turning_below(p, z_r);
        if (std::isnan(bottom))
        {
            out.dx = nan;
            return out;
        }
        add(path.between<full>(p, bottom, z_r, true, false), 2.0);
    }
    return out;
}

// One family of rays: an arm, the interval of p it spans, and which way the
// range runs along it -- up with p when launched down to a first crossing,
// down with p otherwise (a steeper launch dips deeper, or turns higher).
struct Family
{
    Arm arm;
    double p_lo = 0.0;
    double p_hi = 0.0;
    bool rising = true;
};

// The families that reach z_r, in the order they are tried.  With the fastest
// air between receiver and source c_rs, below the receiver c_0r, and above the
// source c_up: a ray passes z_r only if 1/p > c_rs, dips back up through it
// only if 1/p < c_0r, and turns aloft only if 1/p <= c_up.  `reflected` asks
// for the ground-bouncing families instead, which need 1/p above the fastest
// air all the way down.
struct Families
{
    Family list[4];
    int count = 0;
    // The family whose end is the shadow boundary, and that end's p.
    int limiting = -1;
    double limiting_p = 0.0;
};

// int_{z0}^{z1} sin(angle to the horizontal) dz for a ray of invariant p that
// does not turn on [z0, z1], c linear between c0 and c1: the vertical phase per
// unit wavenumber.  With u = p c the integrand is sqrt(1 - u^2) du / (p g), whose
// antiderivative is (u q + asin u) / 2; the difference of the two asin is taken
// as one atan2, and a layer too flat for the quotient is the trapezoid.
double vertical_phase(double p, double z0, double c0, double z1, double c1)
{
    const double dz = z1 - z0;
    if (!(dz > 0.0))
        return 0.0;
    const double u0 = p * c0, u1 = p * c1;
    const double q0 = sine(p, c0), q1 = sine(p, c1);
    const double du = u1 - u0;
    if (std::fabs(du) <= 1e-9 * std::max(std::fabs(u0), std::fabs(u1)))
        return 0.5 * (q0 + q1) * dz;
    return dz * (u1 * q1 - u0 * q0 + std::atan2(u1 * q0 - u0 * q1, u0 * u1 + q0 * q1)) / (2.0 * du);
}

// The largest invariant, at most `pass`, whose ray the wave can follow past the
// fastest air of [z_r, z_s] (c_rs), for families().
//
// A ray of invariant 1/c_rs is horizontal where the air is that fast.  If the
// air around that height has any gradient, the ray turns there over a vertical
// extent the wave resolves at high enough frequency, and the range of the rays
// approaching it converges.  If it is flat -- a calm, isothermal floor, a
// profile written to a few decimals -- the ray runs along it without limit, and
// the downward family's last rays (pass is 1/c_rs less 1e-12) reached 20 km
// across a 3 cm floor, a figure set by that inset alone.  The shadow boundary
// went with them, and with it the hand-over and ShadowTransition; the pairs
// short of it took near-grazing rays whose tubes were tens of times their
// range, 22-28 dB low, where the same profile without the floor was not.
//
// No wave can be steered along a stretch thinner than its vertical wavelength
// at that angle: one confined to a height h diffracts out of it within about
// k h^2.  So around each height where the air is this fast, the turning region
// is taken to run to where the limiting ray has gathered one radian of
// vertical phase, k int sin dz = 1, at the highest frequency the run evaluates
// (Layers::resolution_frequency), and the ray's horizontal travel across it,
// over its height h, is capped at 1.5 k h^2.  In a linear gradient ray theory
// itself travels (4/3) k h^2 there (sqrt(2 R h) against an Airy height h with
// k^2 h^3 = 9 R / 8), 8/9 of the cap, and less on exact arcs; so the cap
// leaves every resolved turning alone and binds only on flat or nearly flat
// air, which gathers no phase.  The margin is what a floor may still add: the
// cap's slack over (4/3) k h^2.  At 2 k h^2 that was 0.26 dB on the demo
// grid's plates under a 0.1 ft calm floor, and 0.67 dB with the run's top band
// at 200 Hz, for a floor a fiftieth of its wavelength; at 1.5, 0.12 and 0.19.  Where it
// binds, the family stops at the invariant whose travel across the region is
// the cap, whose ray still crosses the region at a small but finite angle; the
// boundary is where that ray lands, and past it is shadow -- what the profile
// without the floor gives, to within the cap's slack over what ray theory
// travels in a resolved turning.
//
// The region is a property of the air, not of the path: where the path ends
// before the radian is gathered, the walk goes on past it, up through the
// profile or down to the ground, for as far as the limiting ray could go (not
// into air faster than 1/pass, where it would have turned).  Only the travel
// inside the path, which is all the ray has, is held to the cap.  Cut off at
// the path's ends, a source within a radian of the fastest height -- 1.5 m
// over a 1.2 m microphone at 10 kHz under a plain lapse, 5 m at 200 Hz in a
// light wind -- had a region as short as the path, and the cap bound on
// smooth gradients: one such boundary came in from 789 m to 34 m.
//
// `pass` is returned unchanged, bit for bit, wherever the cap does not bind.
double resolved_pass(const Path& path, double z_s, double z_r, double c_rs, double pass)
{
    const Layers& L = path.L;
    if (!(L.resolution_frequency > 0.0) || !std::isfinite(L.resolution_frequency) || !(z_s > z_r))
        return pass;
    const double k = 2.0 * pi * L.resolution_frequency / c_rs;
    // The heights the path crosses, indexed 0 to n - 1: z_r, the levels
    // strictly between, z_s.  Segment i runs from height i to i + 1.
    const std::size_t first = static_cast<std::size_t>(std::upper_bound(L.z, L.z + L.levels, z_r) - L.z);
    const std::size_t past = static_cast<std::size_t>(std::lower_bound(L.z, L.z + L.levels, z_s) - L.z);
    const int n = 2 + static_cast<int>(past > first ? past - first : 0);
    const double c_r = L.at(z_r), c_s = L.at(z_s);
    auto zb = [&](int i) { return i == 0 ? z_r : i == n - 1 ? z_s : L.z[first + i - 1]; };
    auto cb = [&](int i) { return i == 0 ? c_r : i == n - 1 ? c_s : L.c[first + i - 1]; };
    // Segment i's phase, and a floor under it for one square root: on a layer
    // q^2 is a concave quadratic in z, so q is concave and its trapezoid never
    // exceeds its integral.
    auto phase_of = [&](int i) { return k * vertical_phase(pass, zb(i), cb(i), zb(i + 1), cb(i + 1)); };
    auto phase_floor = [&](int i) {
        return 0.5 * k * (sine(pass, cb(i)) + sine(pass, cb(i + 1))) * (zb(i + 1) - zb(i));
    };
    // Where in segment i the phase gathered from its end nearer the fastest air
    // (the top, if `from_top`) reaches `target`, which is less than the
    // segment's `total`.  Safeguarded Newton on the distance from that end; the
    // first guess takes q to grow as the square root of it, as from a turning
    // height, so the phase as its 3/2 power.
    auto inside_of = [&](double z0, double c0, double z1, double c1, bool from_top, double target,
            double total) {
        const double length = z1 - z0;
        auto speed = [&](double d) {
            return from_top ? c1 + d / length * (c0 - c1) : c0 + d / length * (c1 - c0);
        };
        auto gathered = [&](double d) {
            return from_top ? k * vertical_phase(pass, z1 - d, speed(d), z1, c1)
                            : k * vertical_phase(pass, z0, c0, z0 + d, speed(d));
        };
        double lo = 0.0, hi = length;
        double d = length * std::pow(std::clamp(target / total, 0.0, 1.0), 2.0 / 3.0);
        for (int iteration = 0; iteration < 60 && hi - lo > 1e-12 * length; iteration++)
        {
            const double miss = gathered(d) - target;
            if (std::fabs(miss) <= 1e-10)
                break;
            (miss < 0.0 ? lo : hi) = d;
            const double slope = k * sine(pass, speed(d));
            const double next = slope > 0.0 ? d - miss / slope : 0.5 * (lo + hi);
            d = next > lo && next < hi ? next : 0.5 * (lo + hi);
        }
        return from_top ? z1 - d : z0 + d;
    };
    auto inside = [&](int i, bool from_top, double target, double total) {
        return inside_of(zb(i), cb(i), zb(i + 1), cb(i + 1), from_top, target, total);
    };
    // Past the path's end at z_end (speed c_end), up through the profile or
    // down to the ground: the height at which `remaining` more phase is
    // gathered, or the first height the limiting ray cannot pass -- air faster
    // than 1/pass, where it turns, or the ground -- if that comes first.  Above
    // the top level the air is constant; a ray horizontal there (the profile's
    // fastest air extending to its top) never gathers it, and the region is
    // unbounded.
    const double turning_speed = 1.0 / pass;
    auto beyond = [&](double z_end, double c_end, bool up, double remaining) {
        double z0 = z_end, c0 = c_end;
        // One piece from (z0, c0) toward (z1, c1), cut short where the ray
        // would turn; true when the walk ends in it.
        double edge = z_end;
        auto piece = [&](double z1, double c1) {
            bool turns = false;
            if (c1 > turning_speed)
            {
                z1 = z0 + (z1 - z0) * (turning_speed - c0) / (c1 - c0);
                c1 = turning_speed;
                turns = true;
            }
            const bool upward = z1 > z0;
            const double lo = upward ? z0 : z1, hi = upward ? z1 : z0;
            const double c_lo = upward ? c0 : c1, c_hi = upward ? c1 : c0;
            const double gathered = k * vertical_phase(pass, lo, c_lo, hi, c_hi);
            if (gathered >= remaining)
            {
                edge = inside_of(lo, c_lo, hi, c_hi, !upward, remaining, gathered);
                return true;
            }
            remaining -= gathered;
            z0 = z1;
            c0 = c1;
            edge = z1;
            return turns;
        };
        if (up)
        {
            for (const double* level = std::upper_bound(L.z, L.z + L.levels, z_end); level != L.z + L.levels;
                    ++level)
            {
                if (piece(*level, L.c[level - L.z]))
                    return edge;
            }
            const double q = sine(pass, c0);
            return q > 0.0 ? z0 + remaining / (k * q) : std::numeric_limits<double>::infinity();
        }
        for (const double* level = std::lower_bound(L.z, L.z + L.levels, z_end); level != L.z;)
        {
            --level;
            if (!(*level > 0.0))
                break;
            if (piece(*level, L.c[level - L.z]))
                return edge;
        }
        if (z0 > 0.0)
            piece(0.0, L.at(0.0));
        return edge;
    };
    // The walk from a run of heights at the fastest air, one way: the segment
    // in which one radian is reached and what remains of it there, or none
    // (-1) if the path ends first, with what remains of it then (`beyond`
    // carries on from there).  `next_to_run`: that segment adjoins the run.
    // When the trapezoid alone shows it (total 0), nothing was integrated.
    struct Reach
    {
        int segment = -1;
        double target = 0.0, total = 0.0;
        bool next_to_run = false;
    };
    auto walk = [&](int from, int step, bool single) {
        Reach out;
        double phase = 0.0;
        for (int i = from; step < 0 ? i > 0 : i + 1 < n; i += step)
        {
            const int segment = step < 0 ? i - 1 : i;
            if (single && i == from && phase_floor(segment) >= 1.0)
            {
                out.segment = segment;
                out.next_to_run = true;
                return out;
            }
            const double piece = phase_of(segment);
            if (phase + piece >= 1.0)
            {
                out = {segment, 1.0 - phase, piece, i == from};
                return out;
            }
            phase += piece;
        }
        out.target = 1.0 - phase;
        return out;
    };
    const double fastest = c_rs * (1.0 - 1e-12);
    double result = pass;
    for (int a = 0; a < n; a++)
    {
        if (!(cb(a) >= fastest))
            continue;
        int b = a;
        while (b + 1 < n && cb(b + 1) >= fastest)
            b++;
        const Reach below = walk(a, -1, a == b);
        const Reach above = walk(b, +1, a == b);
        // A single height at the fastest air whose radian is gathered within
        // the segment on each side of it (or that ends the path there) is a
        // resolved turning in a linear gradient: the cap does not bind.
        if (a == b && (a == 0 || below.next_to_run) && (b == n - 1 || above.next_to_run))
        {
            a = b;
            continue;
        }
        auto edge = [&](const Reach& reach, bool from_top, double end, double c_end) {
            if (reach.segment < 0)
                return beyond(end, c_end, !from_top, reach.target);
            if (reach.total == 0.0) // the floor alone placed it: integrate now
                return inside(reach.segment, from_top, 1.0, phase_of(reach.segment));
            return inside(reach.segment, from_top, reach.target, reach.total);
        };
        const double lo = edge(below, true, z_r, c_r);
        const double hi = edge(above, false, z_s, c_s);
        const double h = hi - lo;
        const double cap = 1.5 * k * h * h;
        // The ray's own travel: inside the path.
        const double from = std::max(lo, z_r), to = std::min(hi, z_s);
        if (!(cap < std::numeric_limits<double>::infinity()) || path.between<false>(pass, from, to).dx <= cap)
        {
            a = b;
            continue;
        }
        // The travel across the region rises with p: bisect on s = 1 - p c_rs,
        // in log s, between pass and a vertical ray, keeping the side within
        // the cap.
        double s_lo = std::log(std::max(1.0 - pass * c_rs, 1e-300)), s_hi = 0.0;
        for (int iteration = 0; iteration < 80; iteration++)
        {
            const double mid = 0.5 * (s_lo + s_hi);
            const double p = (1.0 - std::exp(mid)) / c_rs;
            (path.between<false>(p, from, to).dx > cap ? s_lo : s_hi) = mid;
        }
        result = std::min(result, (1.0 - std::exp(s_hi)) / c_rs);
        a = b;
    }
    return result;
}

Families families(const Path& path, double z_s, double z_r, bool reflected)
{
    constexpr double inset = 1e-12;
    const double c_s = path.L.at(z_s);
    const double c_up = path.fastest_above(z_s);
    Families f;
    auto push = [&f](Arm arm, double lo, double hi, bool rising) {
        if (hi > lo)
            f.list[f.count++] = Family {arm, lo, hi, rising};
    };
    if (reflected)
    {
        const double c_0s = path.fastest(0.0, z_s);
        const double hi = std::min(1.0 / c_s, 1.0 / c_0s) * (1.0 - inset);
        push({false, false}, 0.0, hi, true);
        if (c_up > c_s)
            push({true, false}, 1.0 / c_up * (1.0 + inset), hi, false);
        return f;
    }
    const double c_rs = path.fastest(z_r, z_s);
    const double c_0r = path.fastest(0.0, z_r);
    // Not past what the wave resolves at the fastest air (resolved_pass) --
    // unless that air is the source's and there is faster air above it.  Then
    // pass is the horizontal launch, where the downward family meets the
    // families launched up to turn aloft, and is no end of the rays: capping it
    // opened a gap between them that no ray reached, a false shadow inside the
    // lit zone of a downward-refracting profile.
    const double unresolved = std::min(1.0 / c_s, 1.0 / c_rs) * (1.0 - inset);
    const bool launch_junction = c_up > c_s && !(c_rs > c_s);
    const double pass = launch_junction ? unresolved : resolved_pass(path, z_s, z_r, c_rs, unresolved);
    const double dip = 1.0 / c_0r * (1.0 + inset);
    const double turn = 1.0 / c_up * (1.0 + inset);
    push({false, false}, 0.0, pass, true);
    const int down_dip = f.count;
    push({false, true}, dip, pass, false);
    const bool has_down_dip = f.count > down_dip;
    if (c_up > c_s)
    {
        push({true, false}, turn, pass, false);
        push({true, true}, std::max(turn, dip), pass, false);
    }
    // The shadow, when there is one, starts where the downward rays run out:
    // at the end of the dipping family if it exists, else of the first.
    f.limiting = has_down_dip ? down_dip : 0;
    f.limiting_p = has_down_dip ? dip : pass;
    return f;
}

// Solve range(p) = r along a family, starting from p_guess when it lies in
// it.  NaN if r is outside the family's reach.  `jumps`, when given, are where
// the range may be discontinuous (Jumps).
template <typename Range>
double solve(Range&& range, double r, const Family& family, double p_guess,
        const Jumps* jumps = nullptr)
{
    const double sign = family.rising ? 1.0 : -1.0;
    auto f = [&](double p) {
        const double x = range(p);
        return std::isnan(x) ? sign * std::numeric_limits<double>::infinity() : sign * (x - r);
    };
    // f rises with p on the family, whichever way the range runs.
    double lo = family.p_lo, hi = family.p_hi;
    double f_lo = f(lo);
    double f_hi = f(hi);
    if (!(f_lo <= 0.0) || !(f_hi >= 0.0))
        return nan;
    if (p_guess > lo && p_guess < hi)
    {
        const double f_guess = f(p_guess);
        const double width = hi - lo;
        double step = 1e-4 * width;
        if (f_guess < 0.0)
        {
            lo = p_guess;
            f_lo = f_guess;
            double probe = std::min(p_guess + step, hi);
            double f_probe = f(probe);
            while (f_probe < 0.0 && probe < hi)
            {
                lo = probe;
                f_lo = f_probe;
                step *= 4.0;
                probe = std::min(probe + step, hi);
                f_probe = f(probe);
            }
            hi = probe;
            f_hi = f_probe;
        }
        else
        {
            hi = p_guess;
            f_hi = f_guess;
            double probe = std::max(p_guess - step, lo);
            double f_probe = f(probe);
            while (f_probe > 0.0 && probe > lo)
            {
                hi = probe;
                f_hi = f_probe;
                step *= 4.0;
                probe = std::max(probe - step, lo);
                f_probe = f(probe);
            }
            lo = probe;
            f_lo = f_probe;
        }
    }
    // Narrow the bracket to the first continuous piece of it across which f
    // changes sign.  Without this a sign change at a jump alone took the full
    // 200 iterations to be refused: a third of the solves in a surface duct,
    // and nine-tenths of their cost.
    if (jumps != nullptr)
    {
        for (const double* j = std::upper_bound(jumps->p, jumps->p + jumps->count, lo);
                j != jumps->p + jumps->count && *j < hi; ++j)
        {
            const double below = *j * (1.0 - 1e-12);
            const double above = *j * (1.0 + 1e-12);
            if (!(below > lo && above < hi))
                continue;
            const double f_below = f(below);
            if (f_lo <= 0.0 && f_below >= 0.0)
            {
                hi = below;
                f_hi = f_below;
                break;
            }
            lo = above;
            f_lo = f(above);
        }
        if (!(f_lo <= 0.0) || !(f_hi >= 0.0))
            return nan;
    }
    const double tolerance = 1e-10 * std::max(r, 1.0);
    double p = 0.5 * (lo + hi);
    int side = 0;
    double value = nan;
    for (int iteration = 0; iteration < 200; iteration++)
    {
        p = std::isfinite(f_lo) && std::isfinite(f_hi) ? (lo * f_hi - hi * f_lo) / (f_hi - f_lo)
                                                       : 0.5 * (lo + hi);
        value = f(p);
        // Stop, too, once the bracket has closed to adjacent doubles and the
        // next estimate cannot fall strictly inside it.
        if (std::fabs(value) <= tolerance || !(p > lo && p < hi))
            break;
        if (value > 0.0)
        {
            hi = p;
            f_hi = value;
            if (side == 1)
                f_lo *= 0.5;
            side = 1;
        }
        else
        {
            lo = p;
            f_lo = value;
            if (side == -1)
                f_hi *= 0.5;
            side = -1;
        }
    }
    // A bracket can be false.  A ray that never gets there stands in as an
    // infinity of the right sign, so a family whose end is such a ray passes
    // the test above even when its reach stops short of r, and the iteration
    // then closes on that end instead of a root.  In a surface duct that
    // handed back the ray grazing the duct's top as a lit arrival 575 times
    // its range away, 55 dB low.  Only a ray that lands is a root.
    if (!(std::fabs(value) <= 1e-6 * std::max(r, 1.0)))
        return nan;
    return p;
}

// Launch elevation of a ray of invariant p on an arm, and back.
inline double launch_of(double p, double c_s, bool up)
{
    const double a = std::acos(std::min(1.0, p * c_s));
    return up ? a : -a;
}

// r |dr/dtheta| sin|arrival| / cos(launch), by central differences in launch
// angle, kept inside the family.
template <typename Range>
double tube_range(Range&& range, const Family& family, double p, double c_s, double r,
        double arrival)
{
    const double theta = launch_of(p, c_s, family.arm.up);
    auto at = [&](double t) { return range(std::cos(t) / c_s); };
    // The family's ends in launch angle.
    double t_a = launch_of(family.p_lo, c_s, family.arm.up);
    double t_b = launch_of(family.p_hi, c_s, family.arm.up);
    if (t_a > t_b)
        std::swap(t_a, t_b);
    const double delta = std::min(1e-6, 0.25 * std::min(theta - t_a, t_b - theta));
    const double a = std::max(theta - delta, t_a);
    const double b = std::min(theta + delta, t_b);
    // A root that has closed on the end of its family leaves no window: delta is
    // 0, a == b, and the slope below is 0/0.  The tube is then the straight one,
    // r / cos(launch), as it is where the dipping family is a sliver.
    if (!(b > a))
        return r / std::cos(theta);
    const double slope = (at(b) - at(a)) / (b - a);
    return std::sqrt(r * std::fabs(slope) * std::sin(std::fabs(arrival)) / std::cos(theta));
}

} // namespace

double Layers::at(double height) const
{
    if (levels == 1 || height <= z[0])
        return c[0];
    if (height >= z[levels - 1])
        return c[levels - 1];
    const double* upper = std::upper_bound(z, z + levels, height);
    const std::size_t i = static_cast<std::size_t>(upper - z) - 1;
    return c[i] + (height - z[i]) / (z[i + 1] - z[i]) * (c[i + 1] - c[i]);
}

bool Layers::uniform() const
{
    for (std::size_t k = 1; k < levels; k++)
    {
        if (std::fabs(c[k] - c[0]) > 1e-12 * c[0])
            return false;
    }
    return true;
}

LayeredSoundSpeed::LayeredSoundSpeed(const AtmosphereProfile& profile, double frame_bearing,
        double resolution_frequency)
    : heights_(profile.height), resolution_frequency_(resolution_frequency)
{
    const CompassToLocal to_local(frame_bearing);
    for (std::size_t k = 0; k < heights_.size(); k++)
    {
        sound_speed_.push_back(profile.at(heights_[k]).sound_speed());
        const auto wind = to_local(profile.wind_at(heights_[k]));
        wind_x_.push_back(wind[0]);
        wind_y_.push_back(wind[1]);
    }
}

void LayeredSoundSpeed::along(double hx, double hy, double* out) const
{
    for (std::size_t k = 0; k < heights_.size(); k++)
        out[k] = sound_speed_[k] + hx * wind_x_[k] + hy * wind_y_[k];
}

DirectEigenray layered_direct_eigenray(const Layers& layers, double z_s, double z_r, double r,
        double guess)
{
    const Path path {layers};
    const double c_s = layers.at(z_s);
    const double c_r = layers.at(z_r);
    const auto fams = families(path, z_s, z_r, false);
    const Jumps jumps(layers, z_s, z_r);
    const double p_guess = std::cos(guess) / c_s;
    const double chord = std::hypot(r, z_s - z_r);
    const double focus_floor = chord * std::pow(10.0, -max_focusing_gain_db / 20.0);

    DirectEigenray out;
    // A ground-grazing shadow: the fastest air below the source is at the
    // ground, and no ray from aloft reaches the receiver.  Its boundary is where
    // the limiting ray reaches the receiver's height.
    //
    // Faster air above the source is not enough to rule that out: a ray turned
    // back down by it turns up again wherever the air is faster still, so it
    // reaches the receiver only if the air above the source outruns everything
    // between the two -- exactly when families() gives the upward arms.  A
    // morning inversion over an upwind shear is the common case: 330.5 m/s aloft
    // against 339 at the ground, where the shear, not the inversion, casts the
    // shadow.  Testing against the source's own speed instead left those paths
    // with no boundary, so no transition, and the shadow came in as a cliff.
    //
    // ShadowTransition does not take it at face value: a log law's steep
    // lowest meters set it, and the wave resolves them only at high frequency
    // (see there, and layered_direct_ray's linear_boundary).
    //
    // And when the air aloft does outrun it, the rays it turns back come down
    // only past their skip distance, kilometers out for an inversion a few
    // hundred feet up.  Short of that the receiver is in the same ground-grazing
    // shadow, so the boundary is set there too (boundary_under_aloft), and
    // taken back by layered_direct_ray where a ray from aloft or along a duct
    // does arrive.  Without it a source 10 m up, upwind under a morning
    // inversion, had no boundary 200 m out: no transition and no hand-over, but
    // the layered shadow ray's tube, focused by the log law's lowest meters,
    // 6.5 dB over its range, plus the whole ground term -- 12 dB over free
    // field at 63 Hz where the parabolic equation gives 4.
    const bool aloft = path.fastest_above(z_s) > path.fastest(z_r, z_s);
    if (fams.count > 0 && !(path.fastest(0.0, z_s) > layers.at(0.0)))
    {
        const Family& limit = fams.list[fams.limiting];
        out.boundary = trace_p<false>(path, z_s, z_r, fams.limiting_p, limit.arm, false).dx;
        const double root = std::sqrt(2.0 * std::max(z_s, 0.0)) + std::sqrt(2.0 * std::max(z_r, 0.0));
        if (out.boundary > 0.0)
            out.boundary_gradient = -layers.at(0.0) * root * root / (out.boundary * out.boundary);
        out.boundary_under_aloft = aloft;
    }
    for (int k = 0; k < fams.count; k++)
    {
        const Family& family = fams.list[k];
        auto range = [&](double p) { return trace_p<false>(path, z_s, z_r, p, family.arm, false).dx; };
        const double p = solve(range, r, family, p_guess, &jumps);
        if (std::isnan(p))
            continue;
        if (family.arm.up)
        {
            // Turned back from aloft: no ground-grazing shadow about it.
            out.boundary = std::numeric_limits<double>::infinity();
            out.boundary_gradient = 0.0;
            out.boundary_under_aloft = false;
        }
        const auto full = trace_p<true>(path, z_s, z_r, p, family.arm, false);
        const double rising = std::asin(sine(p, c_r));
        out.launch = launch_of(p, c_s, family.arm.up);
        out.arrival = family.arm.second ? rising : -rising;
        out.path_length = full.ds;
        out.travel_time = full.dt;
        out.spreading_range
                = std::max(focus_floor, tube_range(range, family, p, c_s, r, out.arrival));
        return out;
    }

    // Beyond the downward rays' reach, and no ray from aloft arrives: the
    // shadow.  Its limiting ray is the last downward one -- the one that
    // turns at the fastest air below the receiver, or, with no such dip, the
    // one that turns at the fastest air between receiver and source -- and
    // what arrives further creeps along the ground from where it gets to.
    out.lit = false;
    const Family& family = fams.list[fams.limiting];
    auto range = [&](double p) { return trace_p<false>(path, z_s, z_r, p, family.arm, false).dx; };
    const double p_limit = fams.limiting_p;
    const auto full = trace_p<true>(path, z_s, z_r, p_limit, family.arm, false);
    out.launch = launch_of(p_limit, c_s, false);
    out.arrival = 0.0;
    // The tube a little inside the limit, carried out in proportion to range.
    const double inside_theta = out.launch + (family.rising ? -1e-4 : 1e-4);
    const double p_inside = std::clamp(std::cos(inside_theta) / c_s, family.p_lo, family.p_hi);
    const double reach = range(p_inside);
    const double rising = std::asin(sine(p_inside, c_r));
    const double tube = tube_range(range, family, p_inside, c_s, reach,
            family.arm.second ? rising : -rising);
    out.shadow_depth = std::max(r - full.dx, 0.0);
    // A receiver a few millimeters up leaves the dipping family a sliver of
    // heights that deep, whose tube comes out 0 / 0; the limiting ray is then
    // spread as a linear profile's, r / cos(launch).
    const double carried = std::isfinite(tube) && tube > 0.0 && reach > 0.0
            ? tube * r / reach
            : r / std::cos(out.launch);
    out.spreading_range = std::max(focus_floor, carried);
    out.path_length = full.ds + out.shadow_depth;
    out.travel_time = full.dt + out.shadow_depth / layers.at(0.0);
    return out;
}

ReflectedEigenray layered_reflected_eigenray(const Layers& layers, double z_s, double z_r,
        double r, double guess)
{
    const Path path {layers};
    const double c_s = layers.at(z_s);
    const auto fams = families(path, z_s, z_r, true);
    const Jumps jumps(layers, z_s, z_r);
    const double p_guess = std::cos(guess) / c_s;
    ReflectedEigenray out;
    for (int k = 0; k < fams.count; k++)
    {
        const Family& family = fams.list[k];
        auto range = [&](double p) { return trace_p<false>(path, z_s, z_r, p, family.arm, true).dx; };
        const double p = solve(range, r, family, p_guess, &jumps);
        if (std::isnan(p))
            continue;
        const auto full = trace_p<true>(path, z_s, z_r, p, family.arm, true);
        out.exists = true;
        out.launch = launch_of(p, c_s, family.arm.up);
        out.grazing = std::asin(sine(p, layers.at(0.0)));
        out.path_length = full.ds;
        return out;
    }
    return out;
}

namespace {

// One hop of a ray with invariant p that reaches the ground from z_s: up from
// the ground to its turning height and back down.  NaN range if it never turns.
template <bool full>
Crossing hop(const Path& path, double z_s, double p)
{
    Crossing out;
    // Nothing between the ground and z_s turns it, so it first turns above z_s.
    const double top = path.turning_above(p, z_s);
    if (std::isnan(top))
    {
        out.dx = nan;
        return out;
    }
    const auto up = path.between<full>(p, 0.0, top, false, true);
    out.dx = 2.0 * up.dx;
    out.dt = 2.0 * up.dt;
    out.ds = 2.0 * up.ds;
    return out;
}


} // namespace

std::vector<Arrival> layered_arrivals(const Layers& layers, double z_s, double z_r, double r,
        const DirectEigenray& primary, bool bounced_only, int sample_count)
{
    std::vector<Arrival> out;
    const Path path {layers};
    const Jumps jumps(layers, z_s, z_r);
    const double c_s = layers.at(z_s);
    const double c_r = layers.at(z_r);
    const double c_0 = layers.at(0.0);
    const double chord = std::hypot(r, z_s - z_r);
    const double focus_floor = chord * std::pow(10.0, -max_focusing_gain_db / 20.0);
    // A ray going down past the receiver meets the ground only if nothing below
    // the receiver turns it first.
    const double c_0r = path.fastest(0.0, z_r);

    // The duct: rays that reach the ground from the source and turn aloft.
    const double c_0s = path.fastest(0.0, z_s);
    const double c_up = path.fastest_above(z_s);
    constexpr double inset = 1e-12;
    const double duct_lo = 1.0 / c_up * (1.0 + inset);
    const double duct_hi = 1.0 / c_0s * (1.0 - inset);
    const bool duct = c_up > c_0s && duct_hi > duct_lo;
    // dx/dz = p c / q grows with p, so the duct's steepest ray crosses any
    // height interval in the least distance: every hop is at least twice its
    // crossing of [0, z_s], every first descent at least its crossing of
    // [z_r, z_s].  A path shorter than that has no bounced arrival.
    int m_max = 0;
    if (duct)
    {
        const double leg = path.between<false>(duct_lo, 0.0, z_s).dx;
        const double approach = path.between<false>(duct_lo, z_r, z_s).dx;
        if (leg > 0.0)
            m_max = static_cast<int>(std::min((r - approach) / (2.0 * leg), 64.0));
    }

    auto add = [&](const Arm& arm, int m, double p, const Family& whole, auto&& range) {
        const auto first = trace_p<true>(path, z_s, z_r, p, arm, false);
        Arrival a;
        a.bounces = m;
        a.launch = launch_of(p, c_s, arm.up);
        const double rising = std::asin(sine(p, c_r));
        a.arrival = arm.second ? rising : -rising;
        a.grazing = !arm.second && 1.0 / p > c_0r ? std::asin(sine(p, c_0)) : 0.0;
        a.path_length = first.ds;
        a.travel_time = first.dt;
        if (m > 0)
        {
            const auto hops = hop<true>(path, z_s, p);
            a.path_length += m * hops.ds;
            a.travel_time += m * hops.dt;
        }
        a.spreading_range = std::max(focus_floor, tube_range(range, whole, p, c_s, r, a.arrival));
        out.push_back(a);
    };

    // Every family that reaches the receiver directly, and each first-crossing
    // arm's bounced rays.  The downward first-crossing family alone is
    // monotonic (its range is one crossing of [z_r, z_s]); the rest are
    // sampled in launch angle, every crossing of r solved for.  A
    // family turning aloft is not monotonic when the source sits at a local
    // maximum of the sound speed with faster air above: its shallowest rays
    // skim the near-isovelocity layer and travel far, its steepest turn high
    // and travel far, and the arrivals lie between -- a bracket on the ends
    // misses them all.
    //
    // Every family's range is assembled from the same few pieces of a ray of
    // invariant p -- down from the source to the receiver, up from the source
    // to its turning height, the dip below the receiver, and the ground to the
    // receiver for a hop -- so the samples are taken once, on one grid in
    // launch angle over all the families, and shared.  Sampled family by
    // family, a duct's walks through the layers were repeated 3-5 times over.
    //
    // How many samples is a matter of how fast the range oscillates in p: a ray
    // that bounces m times travels X_0(p) + m L(p), and the m L(p) term turns
    // the same sampling interval through m times as many wiggles, so a
    // multi-hop root is missed or, worse, bracketed twice over and lost.  The
    // count grows with the most hops the path can hold; a path with no duct
    // keeps 32.  Against a 2048-sample reference (12 bearings, sources 5-300 m,
    // receivers 0.3 and 1.2 m, 50 m to 10 km, example_inversion.csv), the pairs
    // more than 1 dB out fell from 11.3% to 2.9%, and example_windy.csv's from
    // 1.05% to 0.01%.  The grid is also bunched toward the shallow end, where
    // the rays travel far and the range turns fastest in launch angle: free,
    // and 0.05% instead of 1.05% on the windy profile with 32 samples.
    constexpr int base_samples = 32;
    constexpr int samples_per_hop = 8;
    constexpr int max_samples = 512;
    const int samples = sample_count > 0 ? std::clamp(sample_count, 2, max_samples)
                                         : std::clamp(base_samples + samples_per_hop * m_max, base_samples, max_samples);
    const auto fams = families(path, z_s, z_r, false);
    struct Sampled
    {
        Family family; // over the p it is sampled on
        bool hops;
    };
    Sampled sampled[4];
    int n_sampled = 0;
    bool need_up = false, need_dip = false, need_ground = false;
    double p_min = std::numeric_limits<double>::infinity(), p_max = 0.0;
    for (int k = 0; k < fams.count; k++)
    {
        const Family& family = fams.list[k];
        const Arm arm = family.arm;
        const bool monotonic = !arm.up && !arm.second;
        const bool hops = duct && !arm.second && m_max >= 1;
        // The downward family's own direct ray is never an extra:
        // layered_direct_eigenray tries that family first, so its root, if it
        // has one, is the primary.  `bounced_only` drops the dipping family's
        // direct rays too, which the linear fit stands for, but not those
        // turned back from aloft, which it does not.
        if ((monotonic || (bounced_only && !arm.up)) && !hops)
            continue;
        // The downward family is sampled only for its hops, so over the duct;
        // the others over all of themselves.
        const double p_a = monotonic ? std::max(family.p_lo, duct_lo) : family.p_lo;
        const double p_b = monotonic ? std::min(family.p_hi, duct_hi) : family.p_hi;
        if (!(p_b > p_a))
            continue;
        sampled[n_sampled++] = {Family {arm, p_a, p_b, family.rising}, hops};
        need_up = need_up || arm.up || hops;
        need_dip = need_dip || arm.second;
        need_ground = need_ground || hops;
        p_min = std::min(p_min, p_a);
        p_max = std::max(p_max, p_b);
    }
    if (n_sampled == 0)
        return out;

    // The pieces at one p; NaN where the ray does not turn there.
    struct Row
    {
        double p, down, up, dip, ground;
    };
    auto row = [&](double p) {
        Row w {p, path.between<false>(p, z_r, z_s).dx, nan, nan, nan};
        if (need_up)
        {
            const double top = path.turning_above(p, z_s);
            if (!std::isnan(top))
                w.up = path.between<false>(p, z_s, top, false, true).dx;
        }
        if (need_dip)
        {
            const double bottom = path.turning_below(p, z_r);
            if (!std::isnan(bottom))
                w.dip = path.between<false>(p, bottom, z_r, true, false).dx;
        }
        if (need_ground && p > duct_lo && p < duct_hi)
            w.ground = path.between<false>(p, 0.0, z_r).dx;
        return w;
    };
    // Shallowest first, as each family's samples run.
    Row grid[max_samples];
    {
        const double t_a = std::acos(std::min(1.0, p_max * c_s));
        const double t_b = std::acos(std::min(1.0, p_min * c_s));
        for (int j = 0; j < samples; j++)
        {
            const double u = static_cast<double>(j) / (samples - 1);
            const double t = t_a + (t_b - t_a) * (u * u);
            grid[j] = row(std::clamp(std::cos(t) / c_s, p_min, p_max));
        }
    }

    std::vector<Row> s;
    for (int k = 0; k < n_sampled; k++)
    {
        const Family& whole = sampled[k].family;
        const Arm arm = whole.arm;
        const bool monotonic = !arm.up && !arm.second;
        const bool hops = sampled[k].hops;
        // Whether this arm's own direct rays are listed.
        const bool direct = !monotonic && (!bounced_only || arm.up);
        // Its own ends, and the shared samples strictly between.
        s.clear();
        s.push_back(whole.p_hi == p_max ? grid[0] : row(whole.p_hi));
        for (int j = 0; j < samples; j++)
        {
            const Row& w = grid[j];
            if (w.p < whole.p_hi && w.p > whole.p_lo)
                s.push_back(w);
        }
        s.push_back(whole.p_lo == p_min ? grid[samples - 1] : row(whole.p_lo));
        auto direct_of = [arm](const Row& w) {
            return w.down + (arm.up ? 2.0 * w.up : 0.0) + (arm.second ? 2.0 * w.dip : 0.0);
        };
        auto hop_of = [](const Row& w) { return 2.0 * (w.ground + w.down + w.up); };
        auto direct_range = [&](double p) { return trace_p<false>(path, z_s, z_r, p, arm, false).dx; };
        // The primary's own bracket is split around it (below): this is the
        // primary's arm if it was launched and arrives the same way.
        const bool primary_arm = primary.lit && (primary.launch > 0.0) == arm.up
                && (primary.arrival > 0.0) == arm.second;
        const double p_primary = std::cos(primary.launch) / c_s;
        for (int m = direct ? 0 : 1; m <= (hops ? m_max : 0); m++)
        {
            auto range = [&](double p) {
                return m == 0 ? direct_range(p) : direct_range(p) + m * hop<false>(path, z_s, p).dx;
            };
            auto offset = [&](const Row& w) { return direct_of(w) + (m == 0 ? 0.0 : m * hop_of(w)) - r; };
            // Solve between lo < hi, where range - r is f_lo and f_hi.
            auto solve_between = [&](double lo, double f_lo, double hi, double f_hi) {
                if (!(std::isfinite(f_lo) && std::isfinite(f_hi) && (f_lo < 0.0) != (f_hi < 0.0)))
                    return;
                const double p = solve(range, r, Family {arm, lo, hi, f_hi > f_lo}, nan, &jumps);
                if (!std::isnan(p)
                        && !(m == 0 && primary.lit && std::fabs(launch_of(p, c_s, arm.up) - primary.launch) < 1e-7))
                    add(arm, m, p, whole, range);
            };
            double before = offset(s[0]);
            for (std::size_t j = 1; j < s.size(); j++)
            {
                // p falls along the samples; each bracket runs up in p.
                const double here = offset(s[j]);
                const double lo = s[j].p, hi = s[j - 1].p;
                if (m == 0 && primary_arm && p_primary > lo && p_primary < hi)
                {
                    // The primary is a root here already.  Either side of it,
                    // f changes sign again only if there are others: near a
                    // caustic, two more.  Solving the whole bracket returned
                    // the primary, or one of them, at random.
                    const double gap = 1e-6 * (hi - lo);
                    const double left = std::max(p_primary - gap, lo), right = std::min(p_primary + gap, hi);
                    if (left > lo)
                        solve_between(lo, here, left, range(left) - r);
                    if (right < hi)
                        solve_between(right, range(right) - r, hi, before);
                }
                else
                {
                    solve_between(lo, here, hi, before);
                }
                before = here;
            }
        }
    }
    return out;
}

namespace {

// The path's layers along its bearing, in `scratch`.
Layers path_layers(const LayeredSoundSpeed& air, double hx, double hy, std::vector<double>& scratch)
{
    scratch.resize(air.levels());
    air.along(hx, hy, scratch.data());
    return Layers {air.heights().data(), scratch.data(), air.levels(), air.resolution_frequency()};
}

} // namespace

namespace {

// int_{za}^{zb} sqrt(max(0, p^2 - 1/c(z)^2)) dz through the layers, in seconds:
// the decay exponent, per unit angular frequency, of the field where a ray of
// invariant p cannot go.  Within a layer c is linear, and with u = p c the
// integrand is sqrt(1 - 1/u^2) du / g, whose antiderivative is
// sqrt(u^2 - 1) - acos(1/u) for u >= 1.
double evanescent_integral(const Layers& L, double p, double za, double zb)
{
    auto antiderivative = [](double u) {
        const double e = u - 1.0;
        if (e < 1e-6) // both terms tend to sqrt(2 e); their difference to (2/3) sqrt(2) e^(3/2)
            return e > 0.0 ? (2.0 / 3.0) * std::sqrt(2.0) * e * std::sqrt(e) : 0.0;
        return std::sqrt(u * u - 1.0) - std::acos(1.0 / u);
    };
    auto piece = [&](double z0, double c0, double z1, double c1) {
        const double dz = z1 - z0;
        if (!(dz > 0.0))
            return 0.0;
        const double u0 = p * c0, u1 = p * c1;
        if (u0 <= 1.0 && u1 <= 1.0)
            return 0.0;
        if (std::fabs(c1 - c0) < 1e-9 * c0)
        {
            const double u = 0.5 * (u0 + u1);
            return u > 1.0 ? std::sqrt(p * p - 1.0 / (u / p * u / p)) * dz : 0.0;
        }
        const double g = (c1 - c0) / dz;
        // Only where u > 1: the evanescent end of the layer, if it is not all.
        return (antiderivative(std::max(u1, 1.0)) - antiderivative(std::max(u0, 1.0))) / g;
    };
    double total = 0.0;
    if (!(zb > za))
        return total;
    double z0 = za;
    double c0 = L.at(za);
    const double* begin = L.z;
    const double* end = L.z + L.levels;
    for (const double* level = std::upper_bound(begin, end, za); level != end && *level < zb; ++level)
    {
        const double c1 = L.c[level - begin];
        total += piece(z0, c0, *level, c1);
        z0 = *level;
        c0 = c1;
    }
    return total + piece(z0, c0, zb, L.at(zb));
}

// Ai(x)^2 for real x >= 0.
double airy_squared(double x)
{
    std::complex<double> ai, ai_prime;
    airy(std::complex<double>(x, 0.0), ai, ai_prime);
    return ai.real() * ai.real();
}

constexpr double airy_zero_squared = 0.355028053887817239 * 0.355028053887817239; // Ai(0)^2
// Ai's largest value, at its first maximum x = -1.0188, squared: what a caustic's
// tube is normalized by, lit side included.
constexpr double airy_peak_squared = 0.535656656015700 * 0.535656656015700;

// Ai^2 averaged over a third-octave band, seven frequencies f 2^((k - 3)/18)
// about its center: on a caustic's lit side Ai oscillates, and a band averages
// its fringes.  Both of add_wave_loss's arguments go as omega^(2/3), so the
// mean depends on the band only through x at its center: the k-th frequency's
// is x rho_k, rho_k = 2^((k - 3)/27).  With it, the mean's slope in x.
struct BandAiry
{
    double mean = 0.0;
    double slope = 0.0;
};

BandAiry band_airy_direct(double x)
{
    constexpr int n = 7;
    BandAiry b;
    for (int k = 0; k < n; k++)
    {
        const double rho = std::pow(2.0, 2.0 * (k - 0.5 * (n - 1)) / (9.0 * (n - 1)));
        std::complex<double> ai, ai_prime;
        airy(std::complex<double>(x * rho, 0.0), ai, ai_prime);
        b.mean += ai.real() * ai.real();
        b.slope += 2.0 * ai.real() * ai_prime.real() * rho;
    }
    b.mean /= n;
    b.slope /= n;
    return b;
}

// The mean's logarithm and its slope every 1/128 in x from -60 to 64, for
// cubic Hermite interpolation; built on first use (a function-local static,
// so thread-safe).  Beyond 72 every sub-band's Ai^2 is under e^-726, below
// add_wave_loss's floor of 1e-300 in every band, so the mean is taken as
// zero; elsewhere outside the table it is evaluated directly.  The arrival
// loop asks for it per band per arrival, seven complex Airy evaluations each
// directly -- half a run with --tunneling and --skip_fill in a duct.
struct BandAiryTable
{
    static constexpr double lowest = -60.0;
    static constexpr double highest = 64.0;
    static constexpr double step = 1.0 / 128.0;
    static constexpr double vanishing = 72.0;
    std::vector<double> log_mean, slope;

    BandAiryTable()
    {
        const auto count = static_cast<std::size_t>(std::lround((highest - lowest) / step)) + 1;
        log_mean.resize(count);
        slope.resize(count);
        for (std::size_t k = 0; k < count; k++)
        {
            const BandAiry b = band_airy_direct(lowest + step * static_cast<double>(k));
            log_mean[k] = std::log(b.mean);
            slope[k] = b.slope / b.mean;
        }
    }

    double operator()(double x) const
    {
        if (x >= vanishing)
            return 0.0;
        const double u = (x - lowest) / step;
        if (!(u >= 0.0 && u < static_cast<double>(log_mean.size() - 1)))
            return band_airy_direct(x).mean;
        const auto k = static_cast<std::size_t>(u);
        const double t = u - static_cast<double>(k);
        const double r = 1.0 - t;
        return std::exp((1.0 + 2.0 * t) * r * r * log_mean[k] + t * r * r * step * slope[k]
                + t * t * (3.0 - 2.0 * t) * log_mean[k + 1] - t * t * r * step * slope[k + 1]);
    }
};

double band_airy_squared(double x)
{
    static const BandAiryTable table;
    return table(x);
}

// The Airy length at a turning point where the sound speed is c and falls at
// |g| per meter going up, at wave_reference_frequency: (c^3 / (2 |g| omega^2))^(1/3).
double reference_airy_length(double c, double g)
{
    const double omega = 2.0 * pi * wave_reference_frequency;
    return std::cbrt(c * c * c / (2.0 * std::max(std::fabs(g), 1e-9) * omega * omega));
}

// The sound speed's gradient at z, from the layer it lies in (above z when it
// is a level).
double gradient_at(const Layers& L, double z)
{
    if (L.levels < 2)
        return 0.0;
    const double* upper = std::upper_bound(L.z, L.z + L.levels, z);
    if (upper == L.z || upper == L.z + L.levels)
        return 0.0;
    const std::size_t i = static_cast<std::size_t>(upper - L.z) - 1;
    return (L.c[i + 1] - L.c[i]) / (L.z[i + 1] - L.z[i]);
}

} // namespace

void add_wave_loss(const WaveCorrection& wave, double chord, double spreading_range,
        const std::vector<double>& frequencies, std::size_t first_band, std::size_t last_band,
        double* loss)
{
    if (!wave.active())
        return;
    const double tube_gain_db = 20.0 * std::log10(chord / spreading_range);
    for (std::size_t band = first_band; band < last_band; band++)
    {
        const double f = frequencies[band];
        const double omega = 2.0 * pi * f;
        // The tube is the reference frequency's: the Airy length, and with it
        // the energy near a turning point or caustic, scales as omega^(1/3).
        const double scale = std::cbrt(f / wave_reference_frequency);
        double db = 0.0;
        if (wave.evanescent > 0.0)
        {
            const double a = band_airy_squared(std::pow(1.5 * omega * wave.evanescent, 2.0 / 3.0));
            db -= 10.0 * std::log10(std::max(scale * a / airy_zero_squared, 1e-300));
        }
        if (wave.caustic)
        {
            const double a = band_airy_squared(std::cbrt(omega * omega) * wave.caustic_depth);
            db -= 10.0 * std::log10(std::max(scale * a / airy_peak_squared, 1e-300));
        }
        loss[band] += std::max(db, tube_gain_db - max_focusing_gain_db);
    }
}

std::vector<Arrival> layered_wave_arrivals(const Layers& layers, double z_s, double z_r, double r,
        const WaveOptions& options)
{
    std::vector<Arrival> out;
    if (!options.any() || !(z_s > z_r))
        return out;
    const Path path {layers};
    const double c_s = layers.at(z_s);
    const double c_r = layers.at(z_r);
    const double c_0 = layers.at(0.0);
    const double c_up = path.fastest_above(z_s);
    const double c_rs = path.fastest(z_r, z_s);
    const double c_0r = path.fastest(0.0, z_r);
    constexpr double inset = 1e-12;
    // Nothing turns back from aloft: no family here comes from there.
    if (!(c_up > c_s))
        return out;
    // A tail this deep is under 60 dB at 20 Hz, and smaller at every band.
    constexpr double lowest_frequency = 20.0;
    const double omega_low = 2.0 * pi * lowest_frequency;
    const double max_evanescent = 60.0 / (10.0 / std::log(10.0) * omega_low);
    const double max_caustic_depth = 10.0 / std::cbrt(omega_low * omega_low);
    constexpr int samples = 48;
    constexpr int max_cycles = 16;

    // A family: its rays' range to the receiver (or, tunneled, to their last
    // turning point above it) as a function of p, over [p_lo, p_hi].
    enum class Kind
    {
        direct,
        tunneled
    };
    struct Spec
    {
        Kind kind;
        Arm arm;
        int m; // hops (direct) or cycles (tunneled)
        double p_lo, p_hi;
    };
    // The ray's lowest turning point above the receiver, or NaN.
    auto tunnel_height = [&](double p) {
        const double z_t = path.turning_below(p, z_s);
        return std::isfinite(z_t) && z_t > z_r ? z_t : nan;
    };
    // Range (dx), time and arc of a family's ray of invariant p.
    auto trace = [&](const Spec& s, double p, auto full_tag) {
        constexpr bool full = decltype(full_tag)::value;
        Crossing total;
        auto add = [&total](const Crossing& c, double times) {
            total.dx += times * c.dx;
            total.dt += times * c.dt;
            total.ds += times * c.ds;
        };
        if (s.kind == Kind::direct)
        {
            total = trace_p<full>(path, z_s, z_r, p, s.arm, false);
            if (s.m > 0 && std::isfinite(total.dx))
                add(hop<full>(path, z_s, p), s.m);
            return total;
        }
        const double z_t = tunnel_height(p);
        const double top = path.turning_above(p, z_s);
        if (std::isnan(z_t) || std::isnan(top))
        {
            total.dx = nan;
            return total;
        }
        if (s.arm.up)
            add(path.between<full>(p, z_s, top, false, true), 2.0);
        add(path.between<full>(p, z_t, z_s, true, false), 1.0);
        if (s.m > 0)
            add(path.between<full>(p, z_t, top, true, true), 2.0 * s.m);
        return total;
    };
    auto range_of = [&](const Spec& s) {
        return [&, s](double p) { return trace(s, p, std::false_type {}).dx; };
    };

    // Each family once, its rays' range written as base(p) + m cycle(p): for
    // the direct families the hops of a duct, for the tunneled ones the cycles
    // between their two turning heights.  Sampled once, shared over m.
    struct Group
    {
        Kind kind;
        Arm arm;
        int m_lo, m_hi;
        double p_lo, p_hi;
    };
    std::vector<Group> groups;
    const double turn = 1.0 / c_up * (1.0 + inset);
    const double pass = std::min(1.0 / c_s, 1.0 / c_rs) * (1.0 - inset);
    if (options.skip_fill)
    {
        // The families from aloft that reach the receiver, as families() has
        // them, and in a duct each first-crossing family's bounced rays.
        const auto fams = families(path, z_s, z_r, false);
        for (int k = 0; k < fams.count; k++)
            if (fams.list[k].arm.up)
                groups.push_back({Kind::direct, fams.list[k].arm, 0, 0, fams.list[k].p_lo, fams.list[k].p_hi});
        const double c_0s = path.fastest(0.0, z_s);
        const double duct_hi = std::min(1.0 / c_0s * (1.0 - inset), pass);
        if (c_up > c_0s && duct_hi > turn)
        {
            const double leg = path.between<false>(turn, 0.0, z_s).dx;
            const int m_max = leg > 0.0 ? static_cast<int>(std::min(r / (2.0 * leg) + 1.0, 64.0)) : 0;
            if (m_max >= 1)
            {
                groups.push_back({Kind::direct, {false, false}, 1, m_max, turn, duct_hi});
                groups.push_back({Kind::direct, {true, false}, 1, m_max, turn, duct_hi});
            }
        }
    }
    if (options.tunneling)
    {
        // Rays that pass the source, turn aloft, and turn again between the
        // receiver and the source: 1/p over c_s, at most c_up, under c_rs.
        const double p_lo = std::max(turn, 1.0 / c_rs * (1.0 + inset));
        const double p_hi = (1.0 / c_s) * (1.0 - inset);
        if (p_hi > p_lo)
        {
            groups.push_back({Kind::tunneled, {true, false}, 0, max_cycles, p_lo, p_hi});
            groups.push_back({Kind::tunneled, {false, false}, 1, max_cycles, p_lo, p_hi});
        }
    }

    // The arrival of a family's ray of invariant p.
    auto make = [&](const Spec& s, double p) {
        Arrival a;
        const auto full = trace(s, p, std::true_type {});
        a.bounces = s.kind == Kind::direct ? s.m : 0;
        a.launch = launch_of(p, c_s, s.arm.up);
        a.path_length = full.ds;
        a.travel_time = full.dt;
        if (s.kind == Kind::tunneled)
        {
            // Its tube a reference Airy length above the turning point, where
            // the ray is still a ray.
            const double z_t = tunnel_height(p);
            const double z_a = z_t + reference_airy_length(layers.at(z_t), gradient_at(layers, z_t));
            a.arrival = -std::asin(sine(p, layers.at(z_a)));
            a.grazing = 0.0;
            a.wave.evanescent = std::max(evanescent_integral(layers, p, z_r, z_t), 1e-300);
        }
        else
        {
            const double rising = std::asin(sine(p, c_r));
            a.arrival = s.arm.second ? rising : -rising;
            a.grazing = !s.arm.second && 1.0 / p > c_0r ? std::asin(sine(p, c_0)) : 0.0;
        }
        return a;
    };
    const double omega_ref = 2.0 * pi * wave_reference_frequency;

    for (const Group& group : groups)
    {
        if (!(group.p_hi > group.p_lo))
            continue;
        // Samples uniform in launch angle, shallowest (largest p) first: the
        // m = 0 range and one cycle.
        double p[samples], base[samples], cycle[samples];
        const double t_a = std::acos(std::min(1.0, group.p_hi * c_s));
        const double t_b = std::acos(std::min(1.0, group.p_lo * c_s));
        for (int j = 0; j < samples; j++)
        {
            p[j] = std::clamp(std::cos(t_a + (t_b - t_a) * j / (samples - 1)) / c_s, group.p_lo, group.p_hi);
            const Spec zero {group.kind, group.arm, 0, group.p_lo, group.p_hi};
            base[j] = trace(zero, p[j], std::false_type {}).dx;
            cycle[j] = nan;
            if (group.m_hi > 0)
            {
                if (group.kind == Kind::direct)
                {
                    cycle[j] = hop<false>(path, z_s, p[j]).dx;
                }
                else
                {
                    const double z_t = tunnel_height(p[j]);
                    const double top = path.turning_above(p[j], z_s);
                    if (std::isfinite(z_t) && std::isfinite(top))
                        cycle[j] = 2.0 * path.between<false>(p[j], z_t, top, true, true).dx;
                }
            }
        }
        // A cycle never shorter than the shortest sampled: past r / that, no m reaches.
        double shortest = std::numeric_limits<double>::infinity();
        for (int j = 0; j < samples; j++)
            if (cycle[j] > 0.0)
                shortest = std::min(shortest, cycle[j]);
        const int m_top = group.m_hi == 0 ? 0
                : std::isfinite(shortest) ? std::min(group.m_hi, static_cast<int>(r / shortest) + 1)
                                          : 0;
        for (int m = group.m_lo; m <= m_top; m++)
        {
            const Spec s {group.kind, group.arm, m, group.p_lo, group.p_hi};
            const auto range = range_of(s);
            const Family whole {s.arm, s.p_lo, s.p_hi, false};
            const bool tunneled = s.kind == Kind::tunneled;
            double x[samples];
            for (int j = 0; j < samples; j++)
                x[j] = m == 0 ? base[j] : base[j] + m * cycle[j];

            // The caustics: interior extrema of the range, refined by
            // golden-section search over the samples either side.
            struct Extremum
            {
                int j;
                double p, x, second;
            };
            std::vector<Extremum> extrema;
            for (int j = 1; j + 1 < samples; j++)
            {
                if (!(std::isfinite(x[j - 1]) && std::isfinite(x[j]) && std::isfinite(x[j + 1])))
                    continue;
                const bool minimum = x[j] < x[j - 1] && x[j] <= x[j + 1];
                const bool maximum = x[j] > x[j - 1] && x[j] >= x[j + 1];
                if (!minimum && !maximum)
                    continue;
                // Near enough to matter: the receiver within a twentieth of
                // the range of it on the dark side, or beside it on the lit.
                const double sign = minimum ? 1.0 : -1.0;
                if (!((x[j] - r) * sign > -std::max(std::fabs(x[j - 1] - x[j]), std::fabs(x[j + 1] - x[j]))))
                    continue;
                if (!((x[j] - r) * sign < 0.05 * r) && !options.skip_fill)
                    continue;
                double lo = p[j + 1], hi = p[j - 1];
                constexpr double phi = 0.6180339887498949;
                double c = hi - phi * (hi - lo), d = lo + phi * (hi - lo);
                double fc = sign * range(c), fd = sign * range(d);
                for (int it = 0; it < 60 && hi - lo > 1e-13 * hi; it++)
                {
                    if (fc < fd)
                    {
                        hi = d;
                        d = c;
                        fd = fc;
                        c = hi - phi * (hi - lo);
                        fc = sign * range(c);
                    }
                    else
                    {
                        lo = c;
                        c = d;
                        fc = fd;
                        d = lo + phi * (hi - lo);
                        fd = sign * range(d);
                    }
                }
                const double p_c = 0.5 * (lo + hi);
                const double x_c = range(p_c);
                const double h = 1e-3 * (p[j - 1] - p[j + 1]);
                const double second = (range(p_c + h) - 2.0 * x_c + range(p_c - h)) / (h * h);
                if (std::isfinite(x_c) && std::isfinite(second) && second != 0.0 && (second > 0.0) == minimum)
                    extrema.push_back({j, p_c, x_c, second});
            }

            // A caustic's Airy arrival, `depth` = (X_c - r)(2 / X'')^(1/3):
            // positive on the dark side, negative on the lit.
            auto caustic_arrival = [&](const Extremum& e, double depth) {
                Arrival arrival = make(s, e.p);
                arrival.spreading_range = 0.0;
                if (tunneled && !(arrival.wave.evanescent <= max_evanescent))
                    return arrival;
                const double theta = arrival.launch;
                const double g = e.x * std::fabs(std::sin(theta)) / c_s * std::sin(std::fabs(arrival.arrival))
                        / std::cos(theta);
                if (!(g > 0.0))
                    return arrival;
                // C Ai_peak^2 at the reference frequency, so that the band's
                // Airy factor is at most 1; a tunneled caustic also takes the
                // turning point's 4 pi Ai(0)^2.
                double intensity = 4.0 * pi * std::cbrt(omega_ref)
                        / (std::cbrt(2.0) * g * std::cbrt(e.second * e.second)) * airy_peak_squared;
                if (tunneled)
                    intensity *= 4.0 * pi * airy_zero_squared;
                arrival.spreading_range = 1.0 / std::sqrt(intensity);
                arrival.wave.caustic = true;
                arrival.wave.caustic_depth = depth;
                return arrival;
            };

            // Root intervals [j - 1, j] taken by a caustic's lit-side arrival.
            bool consumed[samples] = {};
            for (const Extremum& e : extrema)
            {
                const double depth = (e.x - r) * std::cbrt(2.0 / e.second);
                if (depth > 0.0)
                {
                    if (!options.skip_fill || !(depth < max_caustic_depth))
                        continue;
                }
                else
                {
                    // Lit side: only the tunneled pair beside its caustic,
                    // whose tubes coalesce there and would take the cap.
                    if (!tunneled || !options.tunneling)
                        continue;
                    const double reach = std::max(std::fabs(x[e.j - 1] - e.x), std::fabs(x[e.j + 1] - e.x));
                    if (!(std::fabs(r - e.x) <= reach))
                        continue;
                    consumed[e.j] = consumed[e.j + 1] = true;
                }
                Arrival arrival = caustic_arrival(e, depth);
                if (arrival.spreading_range > 0.0)
                    out.push_back(arrival);
            }

            if (tunneled && options.tunneling)
            {
                // Every other crossing of r.
                for (int j = 1; j < samples; j++)
                {
                    if (consumed[j])
                        continue;
                    const double f_lo = x[j] - r, f_hi = x[j - 1] - r;
                    if (!(std::isfinite(f_lo) && std::isfinite(f_hi) && (f_lo < 0.0) != (f_hi < 0.0)))
                        continue;
                    const double root = solve(range, r, Family {s.arm, p[j], p[j - 1], f_hi > f_lo}, nan);
                    if (std::isnan(root))
                        continue;
                    Arrival a = make(s, root);
                    if (!(a.wave.evanescent <= max_evanescent))
                        continue;
                    const double tube = tube_range(range, whole, root, c_s, r, a.arrival);
                    if (!(std::isfinite(tube) && tube > 0.0))
                        continue;
                    // Normalized so that the Airy factor's largest value is 1
                    // (add_wave_loss): at the turning point the field is
                    // 4 pi Ai(0)^2 times one ray's.
                    a.spreading_range = tube / std::sqrt(4.0 * pi * airy_zero_squared);
                    out.push_back(a);
                }
            }
        }
    }
    return out;
}

RefractedRay layered_direct_ray(const LayeredSoundSpeed& air, const EffectiveSoundSpeed& linear,
        const Vec3& source, const EffectiveSoundSpeed::Moments& at_source, const Vec3& receiver,
        const EffectiveSoundSpeed::Moments& at_receiver, std::vector<double>& scratch,
        bool linear_past_shadow, const WaveOptions& waves, const LayeredSolution* given,
        LayeredSolution* solved)
{
    // The linear model's ray: the warm start, and the answer wherever the
    // layered solve does not apply -- straight paths, and a receiver at or
    // above the source, which only the linear model's closed form covers.
    RefractedRay ray = direct_ray(linear, source, at_source, receiver, at_receiver);
    const double r = ray.horizontal;
    const double z_s = source.z;
    // A microphone a few centimeters below the plane -- a ground board on
    // ground just under the survey datum -- is solved for on the plane.  The
    // layers hold their lowest speed below z = 0, and the direct-ray solve
    // knows no ground, so it would otherwise send a ray grazing along under
    // the plane to it: on a log-law wind 2.7 cm down gave a "lit" arrival of
    // -0.003 degrees and a ray tube 46 times its range, 33 dB low, where on
    // the plane the pair is in shadow.  The ray still ends at the microphone.
    const double z_r = std::max(receiver.z, 0.0);
    if (r < 1e-3 || !(z_s > z_r))
        return ray;
    const double hx = (receiver.x - source.x) / r;
    const double hy = (receiver.y - source.y) / r;
    const auto layers = path_layers(air, hx, hy, scratch);
    if (layers.uniform())
    {
        // Every ray straight: exactly what straight_ray() gives.
        RefractedRay straight;
        straight.path = straight_ray(source, receiver);
        straight.profile = {layers.c[0], 0.0};
        straight.travel_time = straight.path.range / layers.c[0];
        straight.source_z = z_s;
        straight.receiver_z = z_r;
        straight.horizontal = r;
        return straight;
    }
    const double linear_boundary = ray.direct.boundary;
    DirectEigenray direct;
    bool handover = false;
    std::vector<Arrival> arrivals;
    std::vector<Arrival> arrivals_wave;
    // Given the solves' results (interpolated between two solved states),
    // they stand if this pair hands over as they did.
    bool taken = false;
    if (given != nullptr && given->solved)
    {
        direct = given->direct;
        handover = linear_past_shadow && r > handover_fraction * direct.boundary && ray.curved
                && linear_boundary > direct.boundary;
        if (handover == given->handover)
        {
            arrivals = given->arrivals;
            arrivals_wave = given->waves;
            taken = true;
        }
    }
    if (!taken)
    {
        const double guess = ray.curved ? ray.direct.launch : std::atan2(z_r - z_s, r);
        direct = layered_direct_eigenray(layers, z_s, z_r, r, guess);
        // See below.  The fit's ray then stands for every downward direct ray, so
        // only bounced rays and those turned back from aloft are looked for; the
        // latter exist only under a boundary_under_aloft.
        handover = linear_past_shadow && r > handover_fraction * direct.boundary && ray.curved
                && linear_boundary > direct.boundary;
        arrivals = layered_arrivals(layers, z_s, z_r, r, direct, handover);
        // A boundary set under faster air aloft stands only where nothing from
        // aloft or along a duct arrives (layered_direct_eigenray).  Where something
        // does, past the direct rays' reach, this is the duct's gap and is handled
        // as one, below, exactly as before such boundaries were set.
        if (direct.boundary_under_aloft && !direct.lit && !arrivals.empty())
        {
            direct.boundary = std::numeric_limits<double>::infinity();
            direct.boundary_gradient = 0.0;
            direct.boundary_under_aloft = false;
            if (handover)
            {
                handover = false;
                arrivals = layered_arrivals(layers, z_s, z_r, r, direct, false);
            }
        }
        // --tunneling, --skip_fill: arrivals ray theory does not have.  Kept out
        // of `arrivals` until now, so they take no part in the boundary's removal
        // or the promotion to primary below: each carries a per-band loss that only
        // an extra arrival is given.
        if (waves.any())
        {
            auto added = layered_wave_arrivals(layers, z_s, z_r, r, waves);
            for (Arrival& a : added)
                arrivals_wave.push_back(a);
        }
    }
    if (solved != nullptr)
    {
        solved->solved = true;
        solved->chord = std::hypot(r, z_s - z_r);
        solved->direct = direct;
        solved->handover = handover;
        solved->arrivals = arrivals;
        solved->waves = arrivals_wave;
    }
    // The pair's path, for the arrivals to start from: its separation is theirs.
    auto form_arrivals = [&](RefractedRay& into) {
        arrivals.insert(arrivals.end(), arrivals_wave.begin(), arrivals_wave.end());
        for (Arrival& a : arrivals)
        {
            a.path = into.path;
            const double horizontal = std::cos(a.launch);
            a.path.direction = Vec3(horizontal * hx, horizontal * hy, std::sin(a.launch));
            a.path.range = a.spreading_range;
            a.path.path_length = a.path_length;
            a.path.source_height = z_s;
            a.path.ground_distance = r;
        }
        into.extra = std::move(arrivals);
    };
    // From handover_fraction of the layered boundary on, not only where the
    // layered solve fails.  On a log law the dipping family's range is not
    // monotonic, and reaches some tens of meters past its end, where its
    // near-grazing rays have no reflected partner and the ground term would
    // take the grazing geometry.  Short of the boundary the layered ray is
    // worse still: the steep lowest layer that sets the boundary also focuses
    // the rays that just clear it, and the parabolic equation, which does not
    // resolve that layer, shows neither (see handover_fraction).
    //
    // And the fit's ray stays on past the fit's own boundary, in its shadow.
    // Its limiting ray carries its lit spreading on continuously, where going
    // back to the layered shadow ray there stepped by its different spreading,
    // ~2 dB on a log law.
    if (handover)
    {
        ray.direct.boundary = direct.boundary;
        ray.direct.linear_boundary = linear_boundary;
        ray.direct.boundary_gradient = direct.boundary_gradient;
        ray.from_linear_fit = true;
        ray.interpolated = taken;
        form_arrivals(ray);
        return ray;
    }
    // In a duct, past the reach of every direct ray, the "shadow" is only the
    // gap before the bounced rays: it has no boundary (none of the ground-grazing
    // kind), and what arrives is those rays.  The strongest becomes the primary.
    std::size_t promoted = arrivals.size();
    if (!direct.lit && !std::isfinite(direct.boundary) && !arrivals.empty())
    {
        promoted = 0;
        for (std::size_t k = 1; k < arrivals.size(); k++)
        {
            if (arrivals[k].spreading_range < arrivals[promoted].spreading_range)
                promoted = k;
        }
        const Arrival& a = arrivals[promoted];
        direct.lit = true;
        direct.shadow_depth = 0.0;
        direct.launch = a.launch;
        direct.arrival = a.arrival;
        direct.spreading_range = a.spreading_range;
        direct.path_length = a.path_length;
        direct.travel_time = a.travel_time;
        ray.bounces = a.bounces;
        ray.ground_grazing = a.grazing;
        arrivals.erase(arrivals.begin() + static_cast<std::ptrdiff_t>(promoted));
    }
    ray.curved = true;
    ray.direct = direct;
    ray.direct.linear_boundary = linear_boundary;
    RayPath& path = ray.path;
    const double horizontal = std::cos(direct.launch);
    path.direction = Vec3(horizontal * hx, horizontal * hy, std::sin(direct.launch));
    path.range = direct.spreading_range;
    path.path_length = direct.path_length;
    path.source_height = source.z;
    path.ground_distance = r;
    ray.travel_time = direct.travel_time;
    ray.shadow_depth = direct.shadow_depth;
    // ray.profile stays the path's linear fit: the shadow's attenuation rate
    // reads its gradient as the path's overall refraction.
    ray.interpolated = taken;
    form_arrivals(ray);
    return ray;
}

void layered_finish_ground(const LayeredSoundSpeed& air, RefractedRay& ray, GroundGeometry ground,
        double microphone_height, std::vector<double>& scratch, const LayeredSolution* given,
        LayeredSolution* solved)
{
    if (!ray.curved)
        return;
    // The arrivals' ground geometry is the primary's construction, taken from
    // the ray itself: its partner rising from the next bounce meets the ground
    // at the same invariant, so at the same grazing angle, and its length
    // differs by the receiver's few meters.
    for (Arrival& a : ray.extra)
    {
        if (ground == GroundGeometry::point && a.grazing > 0.0)
        {
            a.path.source_height = std::max(a.path_length * std::sin(a.grazing) - microphone_height, 0.0);
            a.path.ground_distance = a.path_length * std::cos(a.grazing);
        }
        else if (ground == GroundGeometry::plate)
        {
            a.path.source_height = a.path_length * std::sin(std::fabs(a.arrival));
            a.path.ground_distance = a.path_length * std::cos(a.arrival);
        }
    }
    if (ray.from_linear_fit)
    {
        finish_ground(ray, ground, microphone_height);
        return;
    }
    RayPath& path = ray.path;
    const DirectEigenray& direct = ray.direct;
    if (ground == GroundGeometry::point && ray.ground_grazing > 0.0)
    {
        // A primary promoted from the arrivals: their construction.
        path.source_height = std::max(direct.path_length * std::sin(ray.ground_grazing) - microphone_height, 0.0);
        path.ground_distance = direct.path_length * std::cos(ray.ground_grazing);
    }
    else if (ground == GroundGeometry::point && ray.source_z > microphone_height)
    {
        double grazing = 0.0, length = 0.0;
        if (given != nullptr && given->ground_set)
        {
            grazing = given->grazing;
            length = given->length;
        }
        else
        {
            const double r = ray.horizontal;
            const double hx = path.separation.x / r;
            const double hy = path.separation.y / r;
            const auto layers = path_layers(air, hx, hy, scratch);
            const auto reflected = direct.lit
                    ? layered_reflected_eigenray(layers, ray.source_z, microphone_height, r,
                              -std::atan2(ray.source_z + microphone_height, r))
                    : ReflectedEigenray {};
            grazing = reflected.exists ? reflected.grazing : 0.0;
            length = reflected.exists ? reflected.path_length : direct.path_length;
            if (!reflected.exists && direct.lit && direct.arrival < 0.0)
            {
                // The reflected solve takes its families to be monotonic, so it
                // misses its ray wherever the direct one is not (a source at a
                // local maximum of the sound speed), and grazing 0 put the pair in
                // the ground's null.  The partner of a ray arriving from above
                // meets the ground at its own invariant's angle, as the arrivals'
                // do, if nothing below the receiver turns it first.
                const Path through {layers};
                const double p = std::cos(direct.launch) / layers.at(ray.source_z);
                if (1.0 / p > through.fastest(0.0, microphone_height))
                    grazing = std::asin(sine(p, layers.at(0.0)));
            }
        }
        if (solved != nullptr)
        {
            solved->ground_set = true;
            solved->grazing = grazing;
            solved->length = length;
        }
        path.source_height = std::max(length * std::sin(grazing) - microphone_height, 0.0);
        path.ground_distance = length * std::cos(grazing);
    }
    else
    {
        finish_ground(ray, ground, microphone_height);
    }
}

namespace {

double lerp(double a, double b, double t)
{
    return a + t * (b - a);
}

bool close_angles(double a, double b)
{
    return std::fabs(a - b) <= ray_interpolation_angle;
}

bool close_tubes(double a, double b)
{
    return a > 0.0 && b > 0.0 && std::fabs(b / a - 1.0) <= ray_interpolation_spreading;
}

// A length at a and b, each over its own chord, at t, over this chord.
struct Chords
{
    double a, b, t, here;

    bool close(double x, double y) const
    {
        return close_tubes(x / a, y / b);
    }
    double lerp_scaled(double x, double y) const
    {
        return here * lerp(x / a, y / b, t);
    }
};

// Two arrival lists alike, and `out` the one between.
bool interpolate_arrivals(const std::vector<Arrival>& a, const std::vector<Arrival>& b, const Chords& c,
        std::vector<Arrival>& out)
{
    if (a.size() != b.size())
        return false;
    out = a;
    for (std::size_t k = 0; k < a.size(); k++)
    {
        const Arrival& x = a[k];
        const Arrival& y = b[k];
        if (x.bounces != y.bounces || (x.grazing > 0.0) != (y.grazing > 0.0)
                || (x.wave.evanescent > 0.0) != (y.wave.evanescent > 0.0) || x.wave.caustic != y.wave.caustic
                || !close_angles(x.launch, y.launch) || !close_angles(x.arrival, y.arrival)
                || !close_angles(x.grazing, y.grazing) || !c.close(x.spreading_range, y.spreading_range)
                || !c.close(x.path_length, y.path_length))
            return false;
        Arrival& z = out[k];
        z.launch = lerp(x.launch, y.launch, c.t);
        z.arrival = lerp(x.arrival, y.arrival, c.t);
        z.grazing = lerp(x.grazing, y.grazing, c.t);
        z.spreading_range = c.lerp_scaled(x.spreading_range, y.spreading_range);
        z.path_length = c.lerp_scaled(x.path_length, y.path_length);
        z.travel_time = c.lerp_scaled(x.travel_time, y.travel_time);
        z.wave.evanescent = lerp(x.wave.evanescent, y.wave.evanescent, c.t);
        z.wave.caustic_depth = lerp(x.wave.caustic_depth, y.wave.caustic_depth, c.t);
    }
    return true;
}

} // namespace

bool interpolate_layered(const LayeredSolution& a, const LayeredSolution& b, double t, double chord,
        LayeredSolution& out)
{
    if (!a.solved || !b.solved || a.handover != b.handover || !(a.chord > 0.0) || !(b.chord > 0.0)
            || !(chord > 0.0))
        return false;
    const Chords c {a.chord, b.chord, t, chord};
    const DirectEigenray& x = a.direct;
    const DirectEigenray& y = b.direct;
    if (x.lit != y.lit || x.boundary_under_aloft != y.boundary_under_aloft
            || std::isfinite(x.boundary) != std::isfinite(y.boundary)
            || std::isfinite(x.linear_boundary) != std::isfinite(y.linear_boundary)
            || (x.shadow_depth > 0.0) != (y.shadow_depth > 0.0) || !close_angles(x.launch, y.launch)
            || !close_angles(x.arrival, y.arrival) || !c.close(x.spreading_range, y.spreading_range)
            || !c.close(x.path_length, y.path_length))
        return false;
    if (std::isfinite(x.boundary) && !close_tubes(x.boundary, y.boundary))
        return false;
    LayeredSolution between;
    if (!interpolate_arrivals(a.arrivals, b.arrivals, c, between.arrivals)
            || !interpolate_arrivals(a.waves, b.waves, c, between.waves))
        return false;
    out = std::move(between);
    out.solved = true;
    out.chord = chord;
    out.handover = a.handover;
    out.direct = x;
    DirectEigenray& z = out.direct;
    z.launch = lerp(x.launch, y.launch, t);
    z.arrival = lerp(x.arrival, y.arrival, t);
    z.spreading_range = c.lerp_scaled(x.spreading_range, y.spreading_range);
    z.path_length = c.lerp_scaled(x.path_length, y.path_length);
    z.travel_time = c.lerp_scaled(x.travel_time, y.travel_time);
    z.shadow_depth = lerp(x.shadow_depth, y.shadow_depth, t);
    if (std::isfinite(x.boundary))
    {
        z.boundary = lerp(x.boundary, y.boundary, t);
        z.boundary_gradient = lerp(x.boundary_gradient, y.boundary_gradient, t);
    }
    if (std::isfinite(x.linear_boundary))
        z.linear_boundary = lerp(x.linear_boundary, y.linear_boundary, t);
    out.ground_set = false;
    return true;
}

bool interpolate_ground(const LayeredSolution& a, const LayeredSolution& b, double t, double chord,
        LayeredSolution& out)
{
    if (!a.ground_set || !b.ground_set || (a.grazing > 0.0) != (b.grazing > 0.0) || !(a.chord > 0.0)
            || !(b.chord > 0.0) || !(chord > 0.0))
        return false;
    const Chords c {a.chord, b.chord, t, chord};
    if (!close_angles(a.grazing, b.grazing) || !c.close(a.length, b.length))
        return false;
    out.ground_set = true;
    out.grazing = lerp(a.grazing, b.grazing, t);
    out.length = c.lerp_scaled(a.length, b.length);
    return true;
}

} // namespace niceops

namespace niceops {

double resolved_height(const Layers& layers, double frequency, double phase)
{
    const double c_ground = layers.at(0.0);
    const double k = 2.0 * pi * frequency / c_ground;
    // Per unit wavenumber, the phase still to gather.
    double remaining = phase / k;
    // int sqrt(2 (c(0) - c(z)) / c(0)) dz with the shortfall d(z) = c(0) - c(z)
    // linear across each piece: sqrt(2 / c0) int sqrt(d) dz, and over a piece
    // from d0 to d1 (clamped at 0) that is (2/3) (d1^1.5 - d0^1.5) dz / (d1 - d0).
    const double scale = std::sqrt(2.0 / c_ground);
    auto piece = [&](double z0, double c0, double z1, double c1, double& out) {
        // Where along [z0, z1] the remaining phase is gathered; false if not.
        const double dz = z1 - z0;
        if (!(dz > 0.0))
            return false;
        const double d0 = c_ground - c0, d1 = c_ground - c1;
        // The part of the piece where d > 0.
        double a = z0, b = z1, da = d0, db = d1;
        if (d0 <= 0.0 && d1 <= 0.0)
            return false;
        if (d0 < 0.0)
        {
            a = z0 + dz * (-d0) / (d1 - d0);
            da = 0.0;
        }
        if (d1 < 0.0)
        {
            b = z0 + dz * d0 / (d0 - d1);
            db = 0.0;
        }
        const double len = b - a;
        const double slope = (db - da) / len;
        // Integral from a to a + x of scale sqrt(da + slope u) du.
        auto gathered = [&](double x) {
            if (std::fabs(slope) * len <= 1e-12 * std::max(da, db))
                return scale * std::sqrt(da) * x;
            const double e = std::max(da + slope * x, 0.0);
            return scale * 2.0 / 3.0 * (e * std::sqrt(e) - da * std::sqrt(da)) / slope;
        };
        const double total = gathered(len);
        if (total < remaining)
        {
            remaining -= total;
            return false;
        }
        if (std::fabs(slope) * len <= 1e-12 * std::max(da, db))
        {
            out = a + remaining / (scale * std::sqrt(da));
            return true;
        }
        // (da + slope x)^1.5 = da^1.5 + 1.5 slope remaining / scale.
        const double target = da * std::sqrt(da) + 1.5 * slope * remaining / scale;
        const double e = std::cbrt(std::max(target, 0.0) * std::max(target, 0.0));
        out = std::clamp(a + (e - da) / slope, a, b);
        return true;
    };
    double z0 = 0.0, c0 = c_ground, out = 0.0;
    for (std::size_t k_level = 0; k_level < layers.levels; k_level++)
    {
        const double z1 = layers.z[k_level];
        if (!(z1 > z0))
            continue;
        const double c1 = layers.c[k_level];
        if (piece(z0, c0, z1, c1, out))
            return out;
        z0 = z1;
        c0 = c1;
    }
    // Above the top the air is constant.
    const double d = c_ground - c0;
    if (d > 0.0)
        return z0 + remaining / (scale * std::sqrt(d));
    return std::numeric_limits<double>::infinity();
}

Layers tangent_below(const Layers& layers, double height, std::vector<double>& z, std::vector<double>& c,
        double& gradient)
{
    z.clear();
    c.clear();
    gradient = 0.0;
    if (layers.levels > 1 && height >= layers.z[0] && height < layers.z[layers.levels - 1])
    {
        const double* upper = std::upper_bound(layers.z, layers.z + layers.levels, height);
        const std::size_t i = static_cast<std::size_t>(upper - layers.z) - 1;
        gradient = (layers.c[i + 1] - layers.c[i]) / (layers.z[i + 1] - layers.z[i]);
    }
    const double c_h = layers.at(height);
    if (height > 0.0)
    {
        z.push_back(0.0);
        c.push_back(c_h - gradient * height);
    }
    z.push_back(height);
    c.push_back(c_h);
    for (std::size_t k = 0; k < layers.levels; k++)
    {
        if (layers.z[k] > height)
        {
            z.push_back(layers.z[k]);
            c.push_back(layers.c[k]);
        }
    }
    return Layers {z.data(), c.data(), z.size(), layers.resolution_frequency};
}

double layered_shadow_boundary(const Layers& layers, double z_s, double z_r)
{
    const double inf = std::numeric_limits<double>::infinity();
    if (!(z_s > z_r) || layers.uniform())
        return inf;
    const Path path {layers};
    const auto fams = families(path, z_s, z_r, false);
    // As layered_direct_eigenray: the fastest air below the source is at the ground.
    if (!(fams.count > 0) || path.fastest(0.0, z_s) > layers.at(0.0))
        return inf;
    const Family& limit = fams.list[fams.limiting];
    const double boundary = trace_p<false>(path, z_s, z_r, fams.limiting_p, limit.arm, false).dx;
    return boundary > 0.0 && std::isfinite(boundary) ? boundary : inf;
}

void resolved_shadow_bands(const Layers& layers, const LinearProfile& fit, double z_s, double z_r,
        const std::vector<double>& frequencies, std::size_t first_band, std::size_t last_band,
        double phase, bool flush, ResolvedBand* out, bool resolve_linear, bool boundary_law)
{
    const double inf = std::numeric_limits<double>::infinity();
    // The fit's span, as EffectiveSoundSpeed::fit takes it.
    const double z1 = std::min(z_s, z_r);
    const double z2 = std::max(std::max(z_s, z_r), z1 + EffectiveSoundSpeed::min_fit_span);
    const double span = z2 - z1;
    const double z_mid = 0.5 * (z1 + z2);
    std::vector<double> z, c;
    // Every h_f within one layer gives the same resolved profile -- the
    // layer's own line, carried down to the ground -- so its two boundaries
    // are worked out once per layer and shared by the bands in it.
    long cached_layer = -2;
    double cached_layered = inf, cached_linear = inf, cached_gradient = 0.0;
    for (std::size_t band = first_band; band < std::min(last_band, frequencies.size()); band++)
    {
        ResolvedBand& b = out[band];
        b = ResolvedBand {};
        const double f = frequencies[band];
        b.height = resolved_height(layers, f, phase);
        if (!std::isfinite(b.height))
            continue;
        const long layer = layers.levels > 1 && b.height >= layers.z[0] && b.height < layers.z[layers.levels - 1]
                ? static_cast<long>(std::upper_bound(layers.z, layers.z + layers.levels, b.height) - layers.z) - 1
                : b.height < layers.z[0] ? -1 : static_cast<long>(layers.levels);
        if (layer == cached_layer && layer >= 0 && layer < static_cast<long>(layers.levels))
        {
            b.layered_boundary = cached_layered;
            b.linear_boundary = cached_linear;
            b.boundary = std::isfinite(b.layered_boundary) && !boundary_law ? b.layered_boundary
                    : std::isfinite(b.layered_boundary)
                    ? ShadowTransition::band_boundary(b.layered_boundary, b.linear_boundary, z_s, z_r, f,
                              fit.sound_speed, flush)
                    : inf;
            b.gradient = cached_gradient;
            continue;
        }
        double gradient = 0.0;
        const Layers tangent = tangent_below(layers, b.height, z, c, gradient);
        b.layered_boundary = layered_shadow_boundary(tangent, z_s, z_r);
        // The fit's change: the least-squares line through the tangent's change
        // to the profile, d(z) = tangent - c(z) below h_f, over the fit's span.
        double integral = 0.0, moment = 0.0;
        const double top = std::min(z2, b.height);
        if (top > z1)
        {
            std::vector<double> cuts {z1};
            for (std::size_t k = 0; k < layers.levels; k++)
                if (layers.z[k] > z1 && layers.z[k] < top)
                    cuts.push_back(layers.z[k]);
            cuts.push_back(top);
            const double c_h = layers.at(b.height);
            auto change = [&](double height) { return c_h + gradient * (height - b.height) - layers.at(height); };
            for (std::size_t i = 0; i + 1 < cuts.size(); i++)
            {
                const double a0 = cuts[i], dz = cuts[i + 1] - cuts[i];
                const double d0 = change(a0), d1 = change(cuts[i + 1]);
                integral += 0.5 * dz * (d0 + d1);
                moment += dz * (0.5 * a0 * (d0 + d1) + dz * (d0 / 6.0 + d1 / 3.0));
            }
        }
        const double slope = 12.0 * (moment - z_mid * integral) / (span * span * span);
        const double mean = integral / span;
        const LinearProfile resolved_fit {fit.sound_speed + mean - slope * z_mid, fit.gradient + slope};
        b.linear_boundary = shadow_boundary(resolve_linear ? resolved_fit : fit, z_s, z_r);
        b.boundary = std::isfinite(b.layered_boundary) && !boundary_law ? b.layered_boundary
                : std::isfinite(b.layered_boundary)
                ? ShadowTransition::band_boundary(b.layered_boundary, b.linear_boundary, z_s, z_r, f,
                          fit.sound_speed, flush)
                : inf;
        b.gradient = gradient < 0.0 ? -gradient : 0.0;
        cached_layer = layer;
        cached_layered = b.layered_boundary;
        cached_linear = b.linear_boundary;
        cached_gradient = b.gradient;
    }
}

ResolvedShadowTable::ResolvedShadowTable(const LayeredSoundSpeed& air, const EffectiveSoundSpeed& linear,
        const std::vector<double>& frequencies, std::size_t first_band, std::size_t last_band, double phase,
        double z_low, double z_high, const std::vector<double>& receiver_heights, const std::vector<bool>& flush,
        unsigned threads)
{
    last_band = std::min(last_band, frequencies.size());
    if (first_band >= last_band || receiver_heights.empty())
        return;
    first_band_ = first_band;
    bands_ = last_band - first_band;
    bearings_ = static_cast<std::size_t>(std::lround(360.0 / bearing_step));
    root_low_ = std::floor(std::sqrt(std::max(z_low, 0.0)) * height_nodes_per_meter_root)
            / height_nodes_per_meter_root;
    root_step_ = 1.0 / height_nodes_per_meter_root;
    const double root_high = std::sqrt(std::max(z_high, z_low));
    heights_ = static_cast<std::size_t>(std::ceil((root_high - root_low_) / root_step_)) + 2;
    inverse_boundary_.assign(receiver_heights.size(), std::vector<float>(bearings_ * heights_ * bands_, 0.0f));
    gradient_.assign(receiver_heights.size(), std::vector<float>(bearings_ * heights_ * bands_, 0.0f));
    parallel_chunks(bearings_, threads > 0 ? threads : default_thread_count(), [&](std::size_t lo, std::size_t hi) {
        std::vector<double> scratch;
        std::vector<ResolvedBand> bands(frequencies.size());
        for (std::size_t i = lo; i < hi; i++)
        {
            const double bearing = static_cast<double>(i) * bearing_step * pi / 180.0;
            // Compass bearing in the local frame of `air` and `linear` is the
            // caller's: (hx, hy) here are the local-frame direction's x and y.
            const double hx = std::cos(bearing), hy = std::sin(bearing);
            scratch.resize(air.levels());
            air.along(hx, hy, scratch.data());
            const Layers layers {air.heights().data(), scratch.data(), air.levels(), air.resolution_frequency()};
            for (std::size_t j = 0; j < heights_; j++)
            {
                const double root = root_low_ + static_cast<double>(j) * root_step_;
                const double z_s = root * root;
                for (std::size_t q = 0; q < receiver_heights.size(); q++)
                {
                    const double z_r = std::max(receiver_heights[q], 0.0);
                    const std::size_t base = (i * heights_ + j) * bands_;
                    if (!(z_s > z_r))
                        continue;
                    const auto fit = linear.fit(z_s, linear.moments(z_s), z_r, linear.moments(z_r), hx, hy);
                    resolved_shadow_bands(layers, fit, z_s, z_r, frequencies, first_band, last_band, phase,
                            flush[q], bands.data());
                    for (std::size_t b = 0; b < bands_; b++)
                    {
                        const ResolvedBand& rb = bands[first_band + b];
                        inverse_boundary_[q][base + b] = std::isfinite(rb.boundary)
                                ? static_cast<float>(1.0 / rb.boundary)
                                : 0.0f;
                        gradient_[q][base + b] = static_cast<float>(rb.gradient);
                    }
                }
            }
        }
    });
}

void ResolvedShadowTable::lookup(double hx, double hy, double z_s, std::size_t receiver, ResolvedBand* out) const
{
    double angle = std::atan2(hy, hx) * 180.0 / pi;
    if (angle < 0.0)
        angle += 360.0;
    const double x = angle / bearing_step;
    std::size_t i0 = static_cast<std::size_t>(x);
    const double tx = x - static_cast<double>(i0);
    i0 %= bearings_;
    const std::size_t i1 = (i0 + 1) % bearings_;
    const double y = std::clamp((std::sqrt(std::max(z_s, 0.0)) - root_low_) / root_step_, 0.0,
            static_cast<double>(heights_ - 1));
    const std::size_t j0 = std::min(static_cast<std::size_t>(y), heights_ - 2);
    const double ty = y - static_cast<double>(j0);
    const double w00 = (1.0 - tx) * (1.0 - ty), w01 = (1.0 - tx) * ty, w10 = tx * (1.0 - ty), w11 = tx * ty;
    const float* ib = inverse_boundary_[receiver].data();
    const float* g = gradient_[receiver].data();
    const std::size_t c00 = (i0 * heights_ + j0) * bands_, c01 = c00 + bands_;
    const std::size_t c10 = (i1 * heights_ + j0) * bands_, c11 = c10 + bands_;
    for (std::size_t b = 0; b < bands_; b++)
    {
        ResolvedBand& r = out[first_band_ + b];
        const double inverse = w00 * ib[c00 + b] + w01 * ib[c01 + b] + w10 * ib[c10 + b] + w11 * ib[c11 + b];
        r.boundary = inverse > 0.0 ? 1.0 / inverse : std::numeric_limits<double>::infinity();
        r.gradient = w00 * g[c00 + b] + w01 * g[c01 + b] + w10 * g[c10 + b] + w11 * g[c11 + b];
    }
}

} // namespace niceops
