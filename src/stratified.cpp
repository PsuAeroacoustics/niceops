#include "stratified.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "util/constants.hpp"

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
// layer's crossing that much accuracy.  `turn0` / `turn1` pin it to zero.
template <bool full>
inline Crossing cross(double p, double z0, double c0, double z1, double c1, bool turn0 = false,
        bool turn1 = false)
{
    Crossing out;
    const double dz = z1 - z0;
    if (!(dz > 0.0))
        return out;
    const double q0 = turn0 ? 0.0 : sine(p, c0);
    const double q1 = turn1 ? 0.0 : sine(p, c1);
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

    // Sum the crossings of every layer between heights za < zb; `turn_a` /
    // `turn_b` say the ray turns at that end.
    template <bool full>
    Crossing between(double p, double za, double zb, bool turn_a = false, bool turn_b = false) const
    {
        Crossing total;
        if (!(zb > za))
            return total;
        double z0 = za;
        double c0 = L.at(za);
        bool first = true;
        const double* begin = L.z;
        const double* end = L.z + L.levels;
        for (const double* level = std::upper_bound(begin, end, za); level != end && *level < zb;
                ++level)
        {
            const double c1 = L.c[level - begin];
            const auto piece = cross<full>(p, z0, c0, *level, c1, first && turn_a, false);
            total.dx += piece.dx;
            total.dt += piece.dt;
            total.ds += piece.ds;
            z0 = *level;
            c0 = c1;
            first = false;
        }
        const auto piece = cross<full>(p, z0, c0, zb, L.at(zb), first && turn_a, turn_b);
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
        double z0 = z;
        double c0 = L.at(z);
        if (c0 >= target)
            return z;
        const double* begin = L.z;
        const double* end = L.z + L.levels;
        for (const double* level = std::upper_bound(begin, end, z); level != end; ++level)
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
    const double pass = std::min(1.0 / c_s, 1.0 / c_rs) * (1.0 - inset);
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

LayeredSoundSpeed::LayeredSoundSpeed(const AtmosphereProfile& profile, double frame_bearing)
    : heights_(profile.height)
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
    // lowest metres set it, and the wave resolves them only at high frequency
    // (see there, and layered_direct_ray's linear_boundary).
    if (fams.count > 0 && !(path.fastest_above(z_s) > path.fastest(z_r, z_s))
            && !(path.fastest(0.0, z_s) > layers.at(0.0)))
    {
        const Family& limit = fams.list[fams.limiting];
        out.boundary = trace_p<false>(path, z_s, z_r, fams.limiting_p, limit.arm, false).dx;
        const double root = std::sqrt(2.0 * std::max(z_s, 0.0)) + std::sqrt(2.0 * std::max(z_r, 0.0));
        if (out.boundary > 0.0)
            out.boundary_gradient = -layers.at(0.0) * root * root / (out.boundary * out.boundary);
    }
    for (int k = 0; k < fams.count; k++)
    {
        const Family& family = fams.list[k];
        auto range = [&](double p) { return trace_p<false>(path, z_s, z_r, p, family.arm, false).dx; };
        const double p = solve(range, r, family, p_guess, &jumps);
        if (std::isnan(p))
            continue;
        const auto full = trace_p<true>(path, z_s, z_r, p, family.arm, false);
        const double rising = std::asin(sine(p, c_r));
        out.launch = launch_of(p, c_s, family.arm.up);
        out.arrival = family.arm.second ? rising : -rising;
        out.path_length = full.ds;
        out.travel_time = full.dt;
        out.spreading_range
                = std::max(tube_range(range, family, p, c_s, r, out.arrival), focus_floor);
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
    // A receiver a few millimetres up leaves the dipping family a sliver of
    // heights that deep, whose tube comes out 0 / 0; the limiting ray is then
    // spread as a linear profile's, r / cos(launch).
    const double carried = std::isfinite(tube) && tube > 0.0 && reach > 0.0
            ? tube * r / reach
            : r / std::cos(out.launch);
    out.spreading_range = std::max(carried, focus_floor);
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

// The most arrivals kept per pair, nearest-spread first.  A source a few metres
// up in a deep duct has hops of a hundred metres and could otherwise list
// dozens of near-grazing rays at long range, each a sphere lookup.
constexpr std::size_t max_arrivals = 24;

} // namespace

std::vector<Arrival> layered_arrivals(const Layers& layers, double z_s, double z_r, double r,
        const DirectEigenray& primary, bool bounced_only)
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
        a.spreading_range = std::max(tube_range(range, whole, p, c_s, r, a.arrival), focus_floor);
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
    constexpr int samples = 32;
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
        // has one, is the primary.
        if ((monotonic || bounced_only) && !hops)
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
    Row grid[samples];
    {
        const double t_a = std::acos(std::min(1.0, p_max * c_s));
        const double t_b = std::acos(std::min(1.0, p_min * c_s));
        for (int j = 0; j < samples; j++)
        {
            const double t = t_a + (t_b - t_a) * j / (samples - 1);
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
        const bool direct = !monotonic && !bounced_only;
        // Its own ends, and the shared samples strictly between.
        s.clear();
        s.push_back(whole.p_hi == p_max ? grid[0] : row(whole.p_hi));
        for (const Row& w : grid)
        {
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
    if (out.size() > max_arrivals)
    {
        std::partial_sort(out.begin(), out.begin() + max_arrivals, out.end(),
                [](const Arrival& a, const Arrival& b) { return a.spreading_range < b.spreading_range; });
        out.resize(max_arrivals);
    }
    return out;
}

namespace {

// The path's layers along its bearing, in `scratch`.
Layers path_layers(const LayeredSoundSpeed& air, double hx, double hy, std::vector<double>& scratch)
{
    scratch.resize(air.levels());
    air.along(hx, hy, scratch.data());
    return Layers {air.heights().data(), scratch.data(), air.levels()};
}

} // namespace

RefractedRay layered_direct_ray(const LayeredSoundSpeed& air, const EffectiveSoundSpeed& linear,
        const Vec3& source, const EffectiveSoundSpeed::Moments& at_source, const Vec3& receiver,
        const EffectiveSoundSpeed::Moments& at_receiver, std::vector<double>& scratch,
        bool linear_past_shadow)
{
    // The linear model's ray: the warm start, and the answer wherever the
    // layered solve does not apply -- straight paths, and a receiver at or
    // above the source, which only the linear model's closed form covers.
    RefractedRay ray = direct_ray(linear, source, at_source, receiver, at_receiver);
    const double r = ray.horizontal;
    const double z_s = source.z;
    // A microphone a few centimetres below the plane -- a ground board on
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
    const double guess = ray.curved ? ray.direct.launch : std::atan2(z_r - z_s, r);
    auto direct = layered_direct_eigenray(layers, z_s, z_r, r, guess);
    const double linear_boundary = ray.direct.boundary;
    // See below.  The fit's ray then stands for every direct ray: a finite
    // ground-grazing boundary means no ray turned back from aloft reaches the
    // receiver (layered_direct_eigenray), so the direct rays left are the
    // downward ones it replaces, and only bounced rays are looked for.
    const bool handover = linear_past_shadow && r > handover_fraction * direct.boundary && ray.curved
            && linear_boundary > direct.boundary;
    auto arrivals = layered_arrivals(layers, z_s, z_r, r, direct, handover);
    // The pair's path, for the arrivals to start from: its separation is theirs.
    auto form_arrivals = [&](RefractedRay& into) {
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
    // monotonic, and reaches some tens of metres past its end, where its
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
    form_arrivals(ray);
    return ray;
}

void layered_finish_ground(const LayeredSoundSpeed& air, RefractedRay& ray, GroundGeometry ground,
        double microphone_height, std::vector<double>& scratch)
{
    if (!ray.curved)
        return;
    // The arrivals' ground geometry is the primary's construction, taken from
    // the ray itself: its partner rising from the next bounce meets the ground
    // at the same invariant, so at the same grazing angle, and its length
    // differs by the receiver's few metres.
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
        const double r = ray.horizontal;
        const double hx = path.separation.x / r;
        const double hy = path.separation.y / r;
        const auto layers = path_layers(air, hx, hy, scratch);
        const auto reflected = direct.lit
                ? layered_reflected_eigenray(layers, ray.source_z, microphone_height, r,
                          -std::atan2(ray.source_z + microphone_height, r))
                : ReflectedEigenray {};
        double grazing = reflected.exists ? reflected.grazing : 0.0;
        const double length = reflected.exists ? reflected.path_length : direct.path_length;
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
        path.source_height = std::max(length * std::sin(grazing) - microphone_height, 0.0);
        path.ground_distance = length * std::cos(grazing);
    }
    else
    {
        finish_ground(ray, ground, microphone_height);
    }
}

} // namespace niceops
