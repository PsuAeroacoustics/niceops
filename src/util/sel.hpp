#pragma once

// Sound exposure level integration.  ten_db_down_interval() finds the event's
// 10 dB-down duration; sound_exposure_level() integrates it into a SEL:
// 10*log10(sum 10^(0.1*L) dtau), dtau being each sample's share of *reception*
// time.  Kept in a header so the unit tests can reach them.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "logger.hpp"

namespace niceops {

// The sound exposure level integrates the event's 10 dB-down duration, which is
// also what --db_cutoff's margin has to clear.
inline constexpr double sel_window_db = 10.0;

// The 10 dB-down duration: from the first to the last sample at or above
// `threshold` (the peak less sel_window_db).  Samples that dip below the
// threshold between them are inside the interval.  This is the duration of the
// aircraft noise metrics -- the EPNL duration of 14 CFR 36 A36.4.5, and the SEL
// convention that follows it -- and it is what panam's sound_exposure_level
// uses, so a prediction and a measurement are integrated alike.
//
// It replaces the D original's walk out from the peak, which stopped at the
// first sample below the threshold on either side.  That cut an event off at
// its first dip: an approach heard from upstream (the pass, then the hover at
// the approach point) lost its second rise, and on the 2017 approach data the
// walk and this rule differed by up to 0.7 dB at a single microphone.
//
// `level_at(i)` reads levels through an accessor so the peak-expanding path
// can evaluate lazily; `certainly_below(i)`, when it returns true, promises
// level_at(i) < threshold without the caller evaluating it.
template <typename LevelAt, typename CertainlyBelow>
inline std::array<std::size_t, 2> ten_db_down_interval(std::size_t count, double threshold,
        LevelAt&& level_at, CertainlyBelow&& certainly_below)
{
    std::array<std::size_t, 2> index = {0, count - 1};
    for (std::size_t i = 0; i < count; i++)
    {
        if (!certainly_below(i) && level_at(i) >= threshold)
        {
            index[0] = i;
            break;
        }
    }
    for (std::size_t i = count; i-- > index[0];)
    {
        if (!certainly_below(i) && level_at(i) >= threshold)
        {
            index[1] = i;
            break;
        }
    }
    return index;
}

inline std::array<std::size_t, 2> peak_range(const std::vector<double>& levels,
        double down = sel_window_db)
{
    // An empty time history has no peak.  Guarding here rather than letting
    // levels.size() - 1 wrap and max_element dereference end(): that crashed
    // with a segfault well downstream of the real cause, an unparseable
    // trajectory that yielded no states.
    if (levels.empty())
        fatal("peak_range: empty level history; the trajectory produced no states.");
    const double peak = *std::max_element(levels.begin(), levels.end());
    return ten_db_down_interval(
            levels.size(), peak - down, [&levels](std::size_t i) { return levels[i]; },
            [](std::size_t) { return false; });
}

// Each sample's share of reception time, as a multiple of the emission step.
// Sound emitted at t reaches the observer at t + r(t)/c, so an emission step T
// arrives spread over T (1 + (dr/dt)/c): compressed while the source closes,
// stretched while it recedes.  The level history is sampled in emission time;
// the exposure a microphone records accumulates in reception time, so without
// this an approaching source -- all an approach-point observer ever sees -- is
// over-counted by 1 / (1 - M_r).  `range_at(i)` is the source-observer range at
// sample i and `c` the sound speed, in consistent units.
template <typename RangeAt>
inline double reception_factor(std::size_t i, std::size_t count, double T, double c,
        RangeAt&& range_at)
{
    if (count < 2 || !(c > 0.0) || !(T > 0.0))
        return 1.0;
    const std::size_t lo = i > 0 ? i - 1 : i;
    const std::size_t hi = i + 1 < count ? i + 1 : i;
    const double drdt = (range_at(hi) - range_at(lo)) / (static_cast<double>(hi - lo) * T);
    // A helicopter is subsonic, so this stays positive; the floor only keeps a
    // corrupt track from producing a negative exposure.
    return std::max(1.0 + drdt / c, 1e-3);
}

// Mean spacing of a monotonically increasing time series, and how far the worst
// interval strays from it (relative).  The SEL integral weights every sample by
// this, so it must come from the data rather than from an assumed step.
inline double uniform_sample_interval(const std::vector<double>& times, double* relative_spread = nullptr)
{
    if (times.size() < 2)
    {
        if (relative_spread != nullptr)
            *relative_spread = 0.0;
        return 0.0;
    }
    const double mean = (times.back() - times.front()) / static_cast<double>(times.size() - 1);
    if (relative_spread != nullptr)
    {
        double worst = 0.0;
        if (mean > 0.0)
        {
            for (std::size_t i = 1; i < times.size(); i++)
                worst = std::max(worst, std::fabs((times[i] - times[i - 1]) - mean) / mean);
        }
        *relative_spread = worst;
    }
    return mean;
}

// A level with an ambient level mixed in, both in dB, as the sum of their
// energies: 10 log10(10^(L/10) + 10^(A/10)).  Aircraft and ambient noise are
// uncorrelated, so it is their mean squares that add, and two equal levels give
// 3.01 dB more, not 6.02.  The D original summed the pressure amplitudes
// (20 log10(10^(L/20) + 10^(A/20))), which adds the two coherently and read up to
// 3 dB high -- 6.02 dB over two equal levels, 2.4 dB over one 10 dB below the
// ambient where the energy sum gives 0.4.  `ambient_energy` is 10^(A/10), formed
// once per run.  A -infinity level (a culled state) gives the ambient alone.
inline double level_with_ambient(double level_db, double ambient_energy)
{
    return 10.0 * std::log10(std::pow(10.0, 0.1 * level_db) + ambient_energy);
}

// Integrate an already-located window into a sound exposure level, each sample
// weighted by its share of reception time (`weight_at`, see reception_factor;
// pass one that returns 1 for plain emission-time integration).  Ascending
// order, because the accumulation order is part of the answer at 1e-9.
template <typename LevelAt, typename WeightAt>
inline double sel_over_window(std::size_t first, std::size_t last, double T, LevelAt&& level_at,
        WeightAt&& weight_at)
{
    double sum = 0.0;
    for (std::size_t i = first; i <= last; i++)
    {
        sum += std::pow(10.0, 0.1 * level_at(i)) * weight_at(i);
    }
    return 10.0 * std::log10(T * sum);
}

// SEL of a level history sampled every T in emission time.  `ranges`, when
// given (same length, same units as c), integrates in reception time; empty,
// every sample counts T.
inline double sound_exposure_level(const std::vector<double>& levels, double T,
        const std::vector<double>& ranges = {}, double c = 0.0)
{
    auto index = peak_range(levels);
    const bool reception = ranges.size() == levels.size() && c > 0.0;
    return sel_over_window(
            index[0], index[1], T, [&levels](std::size_t i) { return levels[i]; },
            [&](std::size_t i) {
                return reception ? reception_factor(i, ranges.size(), T, c,
                                           [&ranges](std::size_t k) { return ranges[k]; })
                                 : 1.0;
            });
}

// The same integral over the whole history, rather than the event's 10 dB-down
// window (--full_history).  Every sample counts, dips included.
inline double full_sound_exposure_level(const std::vector<double>& levels, double T,
        const std::vector<double>& ranges = {}, double c = 0.0)
{
    if (levels.empty())
        return -std::numeric_limits<double>::infinity();
    const bool reception = ranges.size() == levels.size() && c > 0.0;
    return sel_over_window(
            0, levels.size() - 1, T, [&levels](std::size_t i) { return levels[i]; },
            [&](std::size_t i) {
                return reception ? reception_factor(i, ranges.size(), T, c,
                                           [&ranges](std::size_t k) { return ranges[k]; })
                                 : 1.0;
            });
}

// The day-evening-night metrics (-D, --lden_file, --cnel_file).  Each is the
// energy average over a 24-hour day of every operation's exposure, with a
// penalty added to the level of every sample heard in a penalized period:
//
//   DNL   night 22:00-07:00 +10 dB (14 CFR 150; the D original's metric).
//   Lden  evening 19:00-23:00 +5 dB, night 23:00-07:00 +10 dB, day 07:00-19:00
//         (EU Directive 2002/49/EC Annex I, its default periods).
//   CNEL  evening 19:00-22:00 +10 lg 3 = 4.77 dB (each evening operation weighted
//         three times), night 22:00-07:00 +10 dB (California Code of
//         Regulations Title 21, section 5001).
//
// Written out, Lden = 10 lg[(12 10^(Ld/10) + 4 10^((Le+5)/10) + 8 10^((Ln+10)/10))
// / 24] with Ld, Le, Ln the equivalent levels over each period's own length,
// which is the same number as the penalized energies summed and divided by
// 86 400 s; DNL and CNEL likewise.  The period boundaries belong to the period
// they open: 07:00 is day, 19:00 evening, 22:00 (23:00 for Lden) night.
enum class DayMetric
{
    dnl,
    lden,
    cnel,
};

inline constexpr std::size_t day_metric_count = 3;

// A penalized stretch of the local day, hours, [start, end); start > end wraps
// through midnight.
struct PenaltyPeriod
{
    double start_hour;
    double end_hour;
    double penalty_db;
};

inline std::vector<PenaltyPeriod> penalty_periods(DayMetric metric)
{
    // 10 lg 3, spelled so the constant is visibly the weighting by three.
    const double cnel_evening = 10.0 * std::log10(3.0);
    switch (metric)
    {
    case DayMetric::lden:
        return {{19.0, 23.0, 5.0}, {23.0, 7.0, 10.0}};
    case DayMetric::cnel:
        return {{19.0, 22.0, cnel_evening}, {22.0, 7.0, 10.0}};
    case DayMetric::dnl:
        break;
    }
    return {{22.0, 7.0, 10.0}};
}

// The option that asks for a metric, for messages.
inline const char* day_metric_option(DayMetric metric)
{
    switch (metric)
    {
    case DayMetric::lden:
        return "--lden_file";
    case DayMetric::cnel:
        return "--cnel_file";
    case DayMetric::dnl:
        break;
    }
    return "-D";
}

inline const char* day_metric_name(DayMetric metric)
{
    switch (metric)
    {
    case DayMetric::lden:
        return "day-evening-night level (Lden)";
    case DayMetric::cnel:
        return "community noise equivalent level (CNEL)";
    case DayMetric::dnl:
        break;
    }
    return "day-night level (DNL)";
}

// Seconds since local midnight, in [0, 86400).
inline double local_seconds_of_day(double utcsec, double utc_offset_hours)
{
    double local = std::fmod(utcsec + utc_offset_hours * 3600.0, 86400.0);
    if (local < 0.0)
        local += 86400.0;
    return local;
}

// Whether a sample falls in DNL's night.  It runs 22:00 to 07:00 local, so the
// UTC seconds from the trajectory need the local offset applied before the
// window is tested.
inline bool is_night(double utcsec, double utc_offset_hours)
{
    const double local = local_seconds_of_day(utcsec, utc_offset_hours);
    return local < 7.0 * 3600.0 || local >= 22.0 * 3600.0;
}

// The penalty, dB, a sample heard at `utcsec` takes under `metric`.
inline double schedule_penalty(double utcsec, double utc_offset_hours, DayMetric metric)
{
    const double local = local_seconds_of_day(utcsec, utc_offset_hours);
    for (const auto& period : penalty_periods(metric))
    {
        const double start = period.start_hour * 3600.0;
        const double end = period.end_hour * 3600.0;
        const bool inside = start < end ? (local >= start && local < end)
                                        : (local >= start || local < end);
        if (inside)
            return period.penalty_db;
    }
    return 0.0;
}

// One operation's penalized exposure: the sound exposure level of the history
// with `penalty_db[i]` added to sample i.  The penalty goes on before the
// window is found, so a night operation's window is located on the penalized
// history exactly as SEL's is on the raw one; `full` integrates the whole
// history instead (--full_history), as AAM and NORAH2 do.  An unpenalized
// sample adds exactly 0.0, so a daytime DNL is SEL's own number less the day.
//
// `scratch` holds the penalized copy; it is a parameter so a caller looping over
// observers can reuse one buffer instead of allocating per observer.
inline double penalized_exposure_level(const std::vector<double>& levels, double T,
        const std::vector<double>& penalty_db, std::vector<double>& scratch,
        const std::vector<double>& ranges = {}, double c = 0.0, bool full = false)
{
    scratch.resize(levels.size());
    for (std::size_t i = 0; i < levels.size(); i++)
        scratch[i] = levels[i] + penalty_db[i];
    return full ? full_sound_exposure_level(scratch, T, ranges, c)
                : sound_exposure_level(scratch, T, ranges, c);
}

// Day-night level of one operation: the sound exposure level with a 10 dB
// penalty applied to every sample taken at night, as the D original computes it
// (penalized_exposure_level with DNL's schedule, from night flags).
inline double day_night_level(const std::vector<double>& levels, double T,
        const std::vector<char>& night, std::vector<double>& scratch,
        const std::vector<double>& ranges = {}, double c = 0.0, bool full = false)
{
    std::vector<double> penalty(levels.size());
    for (std::size_t i = 0; i < levels.size(); i++)
        penalty[i] = night[i] ? 10.0 : 0.0;
    return penalized_exposure_level(levels, T, penalty, scratch, ranges, c, full);
}

} // namespace niceops
