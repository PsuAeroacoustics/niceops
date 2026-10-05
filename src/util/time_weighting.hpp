#pragma once

// Time-weighted maximum levels and time above a threshold (--lamax_weighting,
// --lamax_file, --nax_file, --ta_file).
//
// A sound level meter does not report the instantaneous mean square: it reports
// it through a time weighting.  IEC 61672-1 defines two, F (fast) and S (slow),
// as first-order exponential averages of the squared, frequency-weighted
// pressure with time constants of 0.125 s and 1 s, so LAFmax is the largest
// value of
//
//     y(t) = (1/tau) integral_{-inf}^{t} p^2(s) exp(-(t - s)/tau) ds.
//
// ECAC Doc 32 (and NORAH2 after it) report a "LASmax" that is really the largest
// 0.5 s equivalent level, which Doc 32 itself calls a surrogate for S; that is
// leq_half_second here, the largest running 0.5 s average of the mean square.
// `instantaneous` is the largest sample of the history, which is what
// --lamax_file wrote before this option existed and is not a meter quantity: on
// the 2017 calm-day approaches it sat +0.95 dB (median; p95 2.29) above LAFmax.
//
// The history is a level per *emission* time, but a meter hears it in
// *reception* time, where each sample lasts its own share (reception_factor in
// sel.hpp, the same weights the SEL integral uses).  Every function here takes
// the samples as a piecewise-constant mean square over those reception
// durations, back to back, and evaluates the weighting on that signal exactly:
// no resampling, no blocking, so the answer does not depend on a phase or a grid.
//
// Before the first sample the signal is taken to have held the first sample's
// level (the detector starts settled on it), and nothing is counted after the
// last sample ends.  An event that starts quiet, as every footprint history
// does, does not notice either convention.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "sel.hpp"

namespace niceops {

enum class LevelWeighting
{
    instantaneous, // the largest sample: no time weighting
    fast,          // IEC 61672-1 F, tau = 0.125 s
    slow,          // IEC 61672-1 S, tau = 1 s
    leq_half_second, // largest running 0.5 s Leq (ECAC Doc 32's LASmax surrogate)
};

// The weighting's time constant, s: the exponential's tau, or the averaging
// window for leq_half_second.  Zero for instantaneous.
inline double weighting_time_constant(LevelWeighting weighting)
{
    switch (weighting)
    {
    case LevelWeighting::fast:
        return 0.125;
    case LevelWeighting::slow:
        return 1.0;
    case LevelWeighting::leq_half_second:
        return 0.5;
    case LevelWeighting::instantaneous:
        break;
    }
    return 0.0;
}

// The coarsest history step a weighting is meaningful on: half its time
// constant.  Coarser, and the exponential sees fewer than two samples per time
// constant -- its output is then a property of the sampling, not of the signal
// -- so a run asking for it is refused rather than given a number.
inline double coarsest_weighting_step(LevelWeighting weighting)
{
    if (weighting == LevelWeighting::instantaneous)
        return std::numeric_limits<double>::infinity();
    return 0.5 * weighting_time_constant(weighting);
}

// The command-line spelling.
inline const char* to_string(LevelWeighting weighting)
{
    switch (weighting)
    {
    case LevelWeighting::fast:
        return "fast";
    case LevelWeighting::slow:
        return "slow";
    case LevelWeighting::leq_half_second:
        return "leq0.5";
    case LevelWeighting::instantaneous:
        break;
    }
    return "instantaneous";
}

// What the maximum is called in a log line: LAFmax, LASmax, ...
inline const char* maximum_name(LevelWeighting weighting)
{
    switch (weighting)
    {
    case LevelWeighting::fast:
        return "LAFmax";
    case LevelWeighting::slow:
        return "LASmax";
    case LevelWeighting::leq_half_second:
        return "LAeq,0.5s max";
    case LevelWeighting::instantaneous:
        break;
    }
    return "instantaneous LAmax";
}

// Throws nothing: the command line restricts the spelling (CLI::IsMember), and
// anything else falls through to fast, the default.
inline LevelWeighting parse_level_weighting(const std::string& name)
{
    if (name == "instantaneous")
        return LevelWeighting::instantaneous;
    if (name == "slow")
        return LevelWeighting::slow;
    if (name == "leq0.5")
        return LevelWeighting::leq_half_second;
    return LevelWeighting::fast;
}

namespace detail {

// Mean square re (20 uPa)^2 of a level; a culled (-infinity) sample is silence.
inline double mean_square(double level_db)
{
    return std::pow(10.0, 0.1 * level_db);
}

inline double level_of(double mean_square)
{
    return mean_square > 0.0 ? 10.0 * std::log10(mean_square)
                             : -std::numeric_limits<double>::infinity();
}

// The running integral of a piecewise-constant mean square, for the running
// Leq: sample k holds over [start[k], start[k + 1]).  Before 0 the first
// sample's mean square is taken to have held, so E(t) = p0 t for t < 0.
// `at` walks a hint forward, so a monotonic sequence of queries is O(n) overall.
struct RunningEnergy
{
    std::vector<double> start; // n + 1 boundaries, start[0] = 0
    std::vector<double> energy; // E at each boundary
    std::vector<double> power;

    RunningEnergy(const std::vector<double>& levels, const std::vector<double>& durations)
    {
        const std::size_t n = levels.size();
        start.resize(n + 1);
        energy.resize(n + 1);
        power.resize(n);
        start[0] = 0.0;
        energy[0] = 0.0;
        for (std::size_t k = 0; k < n; k++)
        {
            power[k] = mean_square(levels[k]);
            start[k + 1] = start[k] + durations[k];
            energy[k + 1] = energy[k] + power[k] * durations[k];
        }
    }

    double end() const { return start.back(); }

    double at(double t, std::size_t& hint) const
    {
        if (t <= 0.0)
            return power.empty() ? 0.0 : power[0] * t;
        const std::size_t n = power.size();
        if (t >= start[n])
            return energy[n];
        while (hint + 1 < n && start[hint + 1] <= t)
            hint++;
        return energy[hint] + power[hint] * (t - start[hint]);
    }
};

// Calls visit(t_a, W_a, t_b, W_b) for each piece of the running window energy
// W(t) = E(t) - E(t - window), t over [0, end], between consecutive
// breakpoints.  W is linear between them: its breakpoints are where t or
// t - window crosses a sample boundary.
template <typename Visit>
inline void for_each_window_piece(const RunningEnergy& e, double window, Visit&& visit)
{
    const std::size_t n = e.power.size();
    std::size_t hint_now = 0;
    std::size_t hint_back = 0;
    auto window_energy = [&](double t) { return e.at(t, hint_now) - e.at(t - window, hint_back); };
    // Merge the boundaries {start[k]} and {start[k] + window} up to the end.
    std::size_t i = 0;
    std::size_t j = 0;
    double t_prev = 0.0;
    double w_prev = window_energy(0.0);
    while (true)
    {
        const double a = i <= n ? e.start[i] : std::numeric_limits<double>::infinity();
        const double b = j <= n ? e.start[j] + window : std::numeric_limits<double>::infinity();
        double t = std::min(a, b);
        if (t > e.end())
            t = e.end();
        if (t > t_prev)
        {
            const double w = window_energy(t);
            visit(t_prev, w_prev, t, w);
            t_prev = t;
            w_prev = w;
        }
        if (t >= e.end())
            break;
        if (a <= b)
            i++;
        else
            j++;
    }
}

} // namespace detail

// Each sample's reception-time duration, s: the emission step T times its
// reception_factor, exactly the weights the SEL integral applies, so a maximum,
// a time above and an exposure taken from one history all hear it on the same
// clock.  `ranges` and `c` as for sound_exposure_level; empty ranges (or no
// sound speed) give every sample T.
inline void reception_durations(std::size_t count, double T, const std::vector<double>& ranges,
        double c, std::vector<double>& durations)
{
    durations.resize(count);
    const bool reception = ranges.size() == count && c > 0.0;
    for (std::size_t i = 0; i < count; i++)
    {
        durations[i] = T
                * (reception ? reception_factor(i, count, T, c,
                                       [&ranges](std::size_t k) { return ranges[k]; })
                             : 1.0);
    }
}

// The largest time-weighted level of a history, dB.  `durations[i]` is how
// long sample i lasts in reception time (the same length as `levels`).
// Returns -infinity for an empty or silent history.
inline double weighted_maximum(const std::vector<double>& levels,
        const std::vector<double>& durations, LevelWeighting weighting)
{
    if (levels.empty())
        return -std::numeric_limits<double>::infinity();
    switch (weighting)
    {
    case LevelWeighting::instantaneous:
        // The level itself rather than a round trip through the mean square,
        // so this is the number --lamax_file wrote before, to the bit.
        return *std::max_element(levels.begin(), levels.end());
    case LevelWeighting::fast:
    case LevelWeighting::slow:
    {
        // Exact for a piecewise-constant input: across sample i the detector
        // relaxes from y toward p_i as exp(-s / tau), monotonically, so its
        // largest value is at a sample's end (or the settled start).
        const double tau = weighting_time_constant(weighting);
        double y = detail::mean_square(levels.front());
        double best = y;
        for (std::size_t i = 0; i < levels.size(); i++)
        {
            const double p = detail::mean_square(levels[i]);
            const double decay = std::exp(-durations[i] / tau);
            y = p + (y - p) * decay;
            best = std::max(best, y);
        }
        return detail::level_of(best);
    }
    case LevelWeighting::leq_half_second:
    {
        const double window = weighting_time_constant(weighting);
        const detail::RunningEnergy energy(levels, durations);
        double best = 0.0;
        // A linear piece peaks at one of its ends.
        detail::for_each_window_piece(energy, window, [&](double, double w_a, double, double w_b) {
            best = std::max(best, std::max(w_a, w_b));
        });
        if (energy.end() <= 0.0)
            best = energy.power.empty() ? 0.0 : energy.power[0] * window;
        return detail::level_of(best / window);
    }
    }
    return -std::numeric_limits<double>::infinity();
}

// Seconds of reception time during which the time-weighted level exceeds
// `threshold_db` (strictly), over the history's own span: the TA metric.  Each
// weighting is followed exactly between samples, so a level that crosses the
// threshold inside a sample is credited to the instant it crosses.
inline double time_above(const std::vector<double>& levels, const std::vector<double>& durations,
        LevelWeighting weighting, double threshold_db)
{
    if (levels.empty())
        return 0.0;
    const double threshold = detail::mean_square(threshold_db);
    double total = 0.0;
    switch (weighting)
    {
    case LevelWeighting::instantaneous:
        for (std::size_t i = 0; i < levels.size(); i++)
        {
            if (levels[i] > threshold_db)
                total += durations[i];
        }
        return total;
    case LevelWeighting::fast:
    case LevelWeighting::slow:
    {
        const double tau = weighting_time_constant(weighting);
        double y0 = detail::mean_square(levels.front());
        for (std::size_t i = 0; i < levels.size(); i++)
        {
            const double p = detail::mean_square(levels[i]);
            const double d = durations[i];
            const double y1 = p + (y0 - p) * std::exp(-d / tau);
            const bool above0 = y0 > threshold;
            const bool above1 = y1 > threshold;
            if (above0 && above1)
                total += d;
            else if (above0)
            {
                // Falling toward p < threshold: y(s) = p + (y0 - p) e^{-s/tau}
                // reaches the threshold at s = tau ln((y0 - p) / (threshold - p)).
                total += std::clamp(tau * std::log((y0 - p) / (threshold - p)), 0.0, d);
            }
            else if (above1)
            {
                // Rising toward p > threshold: above from that crossing on.
                total += d - std::clamp(tau * std::log((p - y0) / (p - threshold)), 0.0, d);
            }
            y0 = y1;
        }
        return total;
    }
    case LevelWeighting::leq_half_second:
    {
        const double window = weighting_time_constant(weighting);
        const detail::RunningEnergy energy(levels, durations);
        const double limit = threshold * window;
        detail::for_each_window_piece(energy, window, [&](double t_a, double w_a, double t_b, double w_b) {
            const bool above_a = w_a > limit;
            const bool above_b = w_b > limit;
            if (above_a && above_b)
                total += t_b - t_a;
            else if (above_a != above_b)
            {
                // Linear in between: the crossing's share of the piece.
                const double f = (limit - w_a) / (w_b - w_a);
                total += above_a ? f * (t_b - t_a) : (1.0 - f) * (t_b - t_a);
            }
        });
        return total;
    }
    }
    return total;
}

} // namespace niceops
