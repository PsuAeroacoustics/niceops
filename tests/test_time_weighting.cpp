// Time-weighted maxima, time above, and the day-evening-night schedules
// (util/time_weighting.hpp, util/sel.hpp).
//
// The weightings are checked against their closed forms: a first-order
// exponential's step and tone-burst responses (IEC 61672-1 specifies the tone
// burst response as exactly these), and a running average's ramp.  Every
// function evaluates the weighting exactly on a piecewise-constant input, so
// the closed forms hold to rounding, at any sampling, and not only in a limit.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <random>
#include <vector>

#include "util/sel.hpp"
#include "util/time_weighting.hpp"

using namespace niceops;

namespace {

constexpr double silence = -std::numeric_limits<double>::infinity();

// A history of `count` samples each `dt` long: silence, then `level` from
// sample `on` to sample `off` (exclusive), then silence.
struct History
{
    std::vector<double> levels;
    std::vector<double> durations;
};

History burst(double dt, std::size_t count, std::size_t on, std::size_t off, double level)
{
    History h;
    h.levels.assign(count, silence);
    h.durations.assign(count, dt);
    for (std::size_t i = on; i < off; i++)
        h.levels[i] = level;
    return h;
}

bool close(double a, double b, double tolerance)
{
    return std::fabs(a - b) <= tolerance;
}

} // namespace

TEST_CASE("the exponential weightings' step response is the closed form", "[time_weighting]")
{
    // Silence for 1 s, 80 dB for 2 s, silence for 1 s.  y(t) = P (1 - e^{-t/tau})
    // after the onset, so the maximum is at the step's end, and the level is
    // more than 3 dB down for tau ln 2 after the onset and tau ln(2 y_end / P)
    // after the end.
    for (const double dt : {0.001, 0.02, 0.05})
    {
        const auto n = static_cast<std::size_t>(std::lround(4.0 / dt));
        const auto h = burst(dt, n, n / 4, 3 * n / 4, 80.0);
        for (const auto weighting : {LevelWeighting::fast, LevelWeighting::slow})
        {
            const double tau = weighting_time_constant(weighting);
            const double settled = 1.0 - std::exp(-2.0 / tau);
            REQUIRE(close(weighted_maximum(h.levels, h.durations, weighting),
                    80.0 + 10.0 * std::log10(settled), 1e-9));
            const double expected_above = 2.0 - tau * std::log(2.0) + tau * std::log(2.0 * settled);
            REQUIRE(close(time_above(h.levels, h.durations, weighting, 80.0 - 10.0 * std::log10(2.0)),
                    expected_above, 1e-9));
        }
    }
}

TEST_CASE("the tone-burst response is IEC 61672-1's", "[time_weighting]")
{
    // A burst of length T_b from silence reads 10 lg(1 - e^{-T_b/tau}) under the
    // steady level: IEC 61672-1:2013 Table 4 lists -1.0 dB (F, 200 ms),
    // -18.0 (F, 2 ms), -27.0 (F, 0.25 ms) and -7.4 (S, 200 ms), and the meter's
    // tolerance is on those values.
    struct Case
    {
        LevelWeighting weighting;
        double burst;
        double table;
    };
    for (const Case c : {Case {LevelWeighting::fast, 0.2, -1.0}, Case {LevelWeighting::fast, 0.002, -18.0},
                 Case {LevelWeighting::fast, 0.00025, -27.0}, Case {LevelWeighting::slow, 0.2, -7.4}})
    {
        // One sample of the burst's own length between two of silence: exact
        // whatever the sampling, since the input is piecewise constant.
        History h;
        h.levels = {silence, 94.0, silence};
        h.durations = {0.5, c.burst, 0.5};
        const double tau = weighting_time_constant(c.weighting);
        const double response = weighted_maximum(h.levels, h.durations, c.weighting) - 94.0;
        REQUIRE(close(response, 10.0 * std::log10(1.0 - std::exp(-c.burst / tau)), 1e-9));
        REQUIRE(close(response, c.table, 0.05));
    }
}

TEST_CASE("the weighting follows each sample's own duration", "[time_weighting]")
{
    // The same burst cut into unequal pieces is the same signal.  An
    // implementation that took one step for every sample (the emission step)
    // instead of each sample's reception duration fails this.
    History even;
    even.levels = {silence, 90.0, 90.0, 90.0, 90.0, silence};
    even.durations = {0.3, 0.05, 0.05, 0.05, 0.05, 0.3};
    History uneven;
    uneven.levels = {silence, 90.0, 90.0, 90.0, silence};
    uneven.durations = {0.3, 0.01, 0.13, 0.06, 0.3};
    for (const auto weighting : {LevelWeighting::fast, LevelWeighting::slow, LevelWeighting::leq_half_second})
    {
        REQUIRE(close(weighted_maximum(even.levels, even.durations, weighting),
                weighted_maximum(uneven.levels, uneven.durations, weighting), 1e-9));
        REQUIRE(close(time_above(even.levels, even.durations, weighting, 85.0),
                time_above(uneven.levels, uneven.durations, weighting, 85.0), 1e-9));
    }
}

TEST_CASE("reception time compresses an approaching burst", "[time_weighting]")
{
    // A source closing at half the speed of sound: each 0.02 s of emission is
    // heard in 0.01 s (reception_factor 1 - 1/2).  A 0.2 s burst in emission is
    // a 0.1 s burst at the microphone, and the fast detector reads it so.
    const double c = 340.0;
    const double T = 0.02;
    const std::size_t n = 200;
    std::vector<double> ranges(n);
    for (std::size_t i = 0; i < n; i++)
        ranges[i] = 2000.0 - 0.5 * c * T * static_cast<double>(i);
    std::vector<double> levels(n, silence);
    for (std::size_t i = 100; i < 110; i++)
        levels[i] = 85.0;
    std::vector<double> durations;
    reception_durations(n, T, ranges, c, durations);
    REQUIRE(close(durations[100], 0.01, 1e-12));
    REQUIRE(close(weighted_maximum(levels, durations, LevelWeighting::fast),
            85.0 + 10.0 * std::log10(1.0 - std::exp(-0.1 / 0.125)), 1e-9));
    // Without ranges every sample lasts T.
    reception_durations(n, T, {}, c, durations);
    REQUIRE(close(weighted_maximum(levels, durations, LevelWeighting::fast),
            85.0 + 10.0 * std::log10(1.0 - std::exp(-0.2 / 0.125)), 1e-9));
}

TEST_CASE("the running 0.5 s Leq is the closed form", "[time_weighting]")
{
    // A 0.2 s burst: the largest 0.5 s window holds all of it.
    History h;
    h.levels = {silence, 80.0, silence};
    h.durations = {1.0, 0.2, 1.0};
    REQUIRE(close(weighted_maximum(h.levels, h.durations, LevelWeighting::leq_half_second),
            80.0 + 10.0 * std::log10(0.2 / 0.5), 1e-9));
    // A burst longer than the window reads its own level, whatever the phase of
    // the samples (no blocks).
    const auto long_burst = burst(0.013, 300, 77, 231, 80.0);
    REQUIRE(close(weighted_maximum(long_burst.levels, long_burst.durations, LevelWeighting::leq_half_second),
            80.0, 1e-9));
    // The window's energy ramps up over 0.5 s after the onset and down over
    // 0.5 s after the end, crossing half power half a window in on each side,
    // so the time above half power is the burst's own length.
    const auto two_seconds = burst(0.01, 400, 100, 300, 80.0);
    REQUIRE(close(time_above(two_seconds.levels, two_seconds.durations, LevelWeighting::leq_half_second,
                          80.0 - 10.0 * std::log10(2.0)),
            2.0, 1e-9));
    // A quarter of the power is crossed a quarter-window in on the way up and
    // three quarters in on the way down: 2 + 0.5/2 s.
    REQUIRE(close(time_above(two_seconds.levels, two_seconds.durations, LevelWeighting::leq_half_second,
                          80.0 - 10.0 * std::log10(4.0)),
            2.0 + 2.0 * 0.125, 1e-9));
}

TEST_CASE("the running Leq turns where the window's trailing edge crosses a sample", "[time_weighting]")
{
    // Samples of 0.4, 0.3, 0.4 and 0.4 s, the 0.3 s one loud.  The 0.5 s
    // window's energy rises over [0.4, 0.7] to 0.3 P, holds until its trailing
    // edge reaches the burst at t = 0.9 -- not a sample boundary -- and falls to
    // zero at 1.2.  Above 0.3 P (W > 0.15 P) from 0.55 to 1.05.  Taking W as
    // linear between sample ends alone would put the fall's crossing at 1.0.
    History h;
    h.levels = {silence, 80.0, silence, silence};
    h.durations = {0.4, 0.3, 0.4, 0.4};
    const auto w = LevelWeighting::leq_half_second;
    REQUIRE(close(weighted_maximum(h.levels, h.durations, w), 80.0 + 10.0 * std::log10(0.6), 1e-9));
    REQUIRE(close(time_above(h.levels, h.durations, w, 80.0 + 10.0 * std::log10(0.3)), 0.5, 1e-9));
}

TEST_CASE("the instantaneous maximum is the largest sample, to the bit", "[time_weighting]")
{
    const std::vector<double> levels = {61.3, 74.123456789, 70.0, silence, 74.1234567};
    const std::vector<double> durations = {0.02, 0.02, 0.02, 0.02, 0.02};
    REQUIRE(weighted_maximum(levels, durations, LevelWeighting::instantaneous) == 74.123456789);
    // Time above counts whole samples strictly over the threshold.
    REQUIRE(close(time_above(levels, durations, LevelWeighting::instantaneous, 70.0), 0.04, 1e-15));
    REQUIRE(close(time_above(levels, durations, LevelWeighting::instantaneous, 61.0), 0.08, 1e-15));
}

TEST_CASE("slower weightings never read higher", "[time_weighting]")
{
    // Each weighting is an average of the one before it in this list, so on
    // any history instantaneous >= fast >= slow, and instantaneous >= leq0.5.
    std::mt19937 random(20261005);
    std::uniform_real_distribution<double> level(40.0, 95.0);
    std::uniform_real_distribution<double> step(0.005, 0.06);
    for (int trial = 0; trial < 50; trial++)
    {
        std::vector<double> levels(300);
        std::vector<double> durations(300);
        for (std::size_t i = 0; i < levels.size(); i++)
        {
            levels[i] = level(random);
            durations[i] = step(random);
        }
        const double inst = weighted_maximum(levels, durations, LevelWeighting::instantaneous);
        const double fast = weighted_maximum(levels, durations, LevelWeighting::fast);
        const double slow = weighted_maximum(levels, durations, LevelWeighting::slow);
        const double leq = weighted_maximum(levels, durations, LevelWeighting::leq_half_second);
        REQUIRE(inst >= fast);
        REQUIRE(fast >= slow);
        REQUIRE(inst >= leq);
        // Time above at the maximum itself is zero; just under it, positive.
        for (const auto w : {LevelWeighting::fast, LevelWeighting::slow, LevelWeighting::leq_half_second})
        {
            const double peak = weighted_maximum(levels, durations, w);
            REQUIRE(time_above(levels, durations, w, peak + 1e-9) == 0.0);
            REQUIRE(time_above(levels, durations, w, peak - 0.5) > 0.0);
        }
    }
}

TEST_CASE("the one-pass maximum does not depend on the time-above threshold", "[time_weighting]")
{
    // One pass computes both metrics; the maximum must be the same, to the
    // bit, whatever threshold rides along (+infinity is the LAmax-only call),
    // and +infinity is never exceeded.  Histories with culled stretches, single
    // samples and an empty one, for every weighting.
    std::mt19937 random(20261010);
    std::uniform_real_distribution<double> level(40.0, 95.0);
    std::uniform_real_distribution<double> step(0.005, 0.06);
    std::uniform_int_distribution<int> length(0, 400);
    std::bernoulli_distribution culled(0.15);
    const double never = std::numeric_limits<double>::infinity();
    for (int trial = 0; trial < 60; trial++)
    {
        const std::size_t n = trial < 2 ? static_cast<std::size_t>(trial) : static_cast<std::size_t>(length(random));
        std::vector<double> levels(n);
        std::vector<double> durations(n);
        for (std::size_t i = 0; i < n; i++)
        {
            levels[i] = trial % 7 == 3 || culled(random) ? silence : level(random);
            durations[i] = step(random);
        }
        for (const auto w : {LevelWeighting::instantaneous, LevelWeighting::fast, LevelWeighting::slow,
                     LevelWeighting::leq_half_second})
        {
            const auto alone = weighted_maximum_and_time_above(levels, durations, w, never);
            REQUIRE(alone.time_above == 0.0);
            for (const double threshold : {30.0, 60.0, 70.0, 80.0, 90.0, 120.0})
            {
                INFO("trial " << trial << ", n " << n << ", " << to_string(w) << ", threshold " << threshold);
                REQUIRE(weighted_maximum_and_time_above(levels, durations, w, threshold).maximum == alone.maximum);
            }
        }
    }
}

TEST_CASE("a weighting is refused on a history coarser than half its time constant", "[time_weighting]")
{
    REQUIRE(coarsest_weighting_step(LevelWeighting::fast) == 0.0625);
    REQUIRE(coarsest_weighting_step(LevelWeighting::slow) == 0.5);
    REQUIRE(coarsest_weighting_step(LevelWeighting::leq_half_second) == 0.25);
    REQUIRE(std::isinf(coarsest_weighting_step(LevelWeighting::instantaneous)));
    for (const auto w : {LevelWeighting::fast, LevelWeighting::slow, LevelWeighting::leq_half_second,
                 LevelWeighting::instantaneous})
        REQUIRE(parse_level_weighting(to_string(w)) == w);
}

TEST_CASE("the day-evening-night schedules penalize their own periods", "[sel]")
{
    const double hour = 3600.0;
    const double cnel_evening = 10.0 * std::log10(3.0);
    struct Probe
    {
        double local_hour;
        double dnl;
        double lden;
        double cnel;
    };
    for (const Probe p : {Probe {0.0, 10, 10, 10}, Probe {6.999, 10, 10, 10}, Probe {7.0, 0, 0, 0},
                 Probe {12.0, 0, 0, 0}, Probe {18.999, 0, 0, 0}, Probe {19.0, 0, 5, cnel_evening},
                 Probe {21.999, 0, 5, cnel_evening}, Probe {22.0, 10, 5, 10}, Probe {22.999, 10, 5, 10},
                 Probe {23.0, 10, 10, 10}, Probe {23.999, 10, 10, 10}})
    {
        const double utc = p.local_hour * hour;
        REQUIRE(schedule_penalty(utc, 0.0, DayMetric::dnl) == p.dnl);
        REQUIRE(schedule_penalty(utc, 0.0, DayMetric::lden) == p.lden);
        REQUIRE(schedule_penalty(utc, 0.0, DayMetric::cnel) == p.cnel);
        // The DNL schedule and is_night agree.
        REQUIRE((schedule_penalty(utc, 0.0, DayMetric::dnl) > 0.0) == is_night(utc, 0.0));
    }
    // The offset shifts UTC into local time: 20:00 UTC is 15:00 at UTC-5.
    REQUIRE(schedule_penalty(20.0 * hour, -5.0, DayMetric::lden) == 0.0);
    REQUIRE(schedule_penalty(20.0 * hour, 0.0, DayMetric::lden) == 5.0);
    // 02:00 UTC at UTC-5 is 21:00 the day before: evening under both.
    REQUIRE(schedule_penalty(2.0 * hour, -5.0, DayMetric::lden) == 5.0);
    REQUIRE(schedule_penalty(2.0 * hour, -5.0, DayMetric::cnel) == cnel_evening);
    REQUIRE(schedule_penalty(2.0 * hour, -5.0, DayMetric::dnl) == 0.0);
    // 04:00 UTC is 23:00: night under all three.
    REQUIRE(schedule_penalty(4.0 * hour, -5.0, DayMetric::lden) == 10.0);
}

TEST_CASE("a day of steady noise gives each metric's period weighting", "[sel]")
{
    // A steady L for 24 hours, integrated whole: DNL = L + 10 lg((15 + 9*10)/24),
    // Lden = L + 10 lg((12 + 4*10^0.5 + 8*10)/24), CNEL = L + 10 lg((12 + 3*3 +
    // 9*10)/24) -- the directive's and Title 21's own formulas.
    const double L = 55.0;
    const double T = 10.0;
    const std::size_t n = 8640;
    std::vector<double> levels(n, L);
    std::vector<double> scratch;
    struct Expected
    {
        DayMetric metric;
        double offset;
    };
    for (const Expected e : {Expected {DayMetric::dnl, 10.0 * std::log10((15.0 + 90.0) / 24.0)},
                 Expected {DayMetric::lden, 10.0 * std::log10((12.0 + 4.0 * std::sqrt(10.0) + 80.0) / 24.0)},
                 Expected {DayMetric::cnel, 10.0 * std::log10((12.0 + 9.0 + 90.0) / 24.0)}})
    {
        std::vector<double> penalty(n);
        for (std::size_t i = 0; i < n; i++)
            penalty[i] = schedule_penalty((static_cast<double>(i) + 0.5) * T, 0.0, e.metric);
        const double exposure = penalized_exposure_level(levels, T, penalty, scratch, {}, 0.0, true);
        REQUIRE(close(exposure - 10.0 * std::log10(86400.0), L + e.offset, 1e-9));
    }
}

TEST_CASE("--full_history governs the penalized exposure", "[sel]")
{
    // An event and a second one 15 dB down: the 10 dB-down window leaves the
    // second out, the whole history keeps it.
    const double T = 0.5;
    std::vector<double> levels = {40, 60, 80, 60, 40, 40, 65, 40};
    const std::vector<double> zero(levels.size(), 0.0);
    std::vector<double> scratch;
    REQUIRE(penalized_exposure_level(levels, T, zero, scratch) == sound_exposure_level(levels, T));
    REQUIRE(penalized_exposure_level(levels, T, zero, scratch, {}, 0.0, true)
            == full_sound_exposure_level(levels, T));
    REQUIRE(full_sound_exposure_level(levels, T) > sound_exposure_level(levels, T) + 0.1);
    // day_night_level, the D-era spelling, is the DNL schedule's case of it.
    const std::vector<char> night(levels.size(), 1);
    const std::vector<double> ten(levels.size(), 10.0);
    REQUIRE(day_night_level(levels, T, night, scratch, {}, 0.0, true)
            == penalized_exposure_level(levels, T, ten, scratch, {}, 0.0, true));
}

TEST_CASE("the period levels count their own period's samples, over its own length", "[sel]")
{
    const double hour = 3600.0;
    const double outside = -std::numeric_limits<double>::infinity();
    struct Probe
    {
        double local_hour;
        double day, evening, night;
    };
    // Each boundary belongs to the period it opens.
    for (const Probe p : {Probe {0.0, outside, outside, 0}, Probe {6.999, outside, outside, 0},
                 Probe {7.0, 0, outside, outside}, Probe {18.999, 0, outside, outside},
                 Probe {19.0, outside, 0, outside}, Probe {22.999, outside, 0, outside},
                 Probe {23.0, outside, outside, 0}})
    {
        const double utc = p.local_hour * hour;
        REQUIRE(schedule_penalty(utc, 0.0, DayMetric::lday) == p.day);
        REQUIRE(schedule_penalty(utc, 0.0, DayMetric::levening) == p.evening);
        REQUIRE(schedule_penalty(utc, 0.0, DayMetric::lnight) == p.night);
    }
    const ClockWindow sixteen {7.0, 23.0};
    REQUIRE(schedule_penalty(6.5 * hour, 0.0, DayMetric::laeq, sixteen) == outside);
    REQUIRE(schedule_penalty(7.0 * hour, 0.0, DayMetric::laeq, sixteen) == 0.0);
    REQUIRE(schedule_penalty(23.0 * hour, 0.0, DayMetric::laeq, sixteen) == outside);
    // A window through midnight.
    const ClockWindow late {22.0, 6.0};
    REQUIRE(schedule_penalty(23.5 * hour, 0.0, DayMetric::laeq, late) == 0.0);
    REQUIRE(schedule_penalty(5.5 * hour, 0.0, DayMetric::laeq, late) == 0.0);
    REQUIRE(schedule_penalty(12.0 * hour, 0.0, DayMetric::laeq, late) == outside);

    REQUIRE(averaging_seconds(DayMetric::lday) == 12.0 * hour);
    REQUIRE(averaging_seconds(DayMetric::levening) == 4.0 * hour);
    REQUIRE(averaging_seconds(DayMetric::lnight) == 8.0 * hour);
    REQUIRE(averaging_seconds(DayMetric::laeq, sixteen) == 16.0 * hour);
    REQUIRE(averaging_seconds(DayMetric::laeq, late) == 8.0 * hour);
    REQUIRE(averaging_seconds(DayMetric::laeq, ClockWindow {0.0, 24.0}) == 24.0 * hour);
    REQUIRE(averaging_seconds(DayMetric::lden) == 24.0 * hour);
}

TEST_CASE("Lden is its period levels recombined", "[sel]")
{
    // EU Directive 2002/49/EC Annex I: Lden = 10 lg[(12 10^(Ld/10) + 4 10^((Le+5)/10)
    // + 8 10^((Ln+10)/10)) / 24].  One level history per period, a sample each
    // hour of the day with a different level, so every period has energy.
    std::vector<double> levels, lden_penalty, day, evening, night;
    for (int h = 0; h < 24; h++)
    {
        levels.push_back(60.0 + h);
        const double utc = (h + 0.5) * 3600.0;
        lden_penalty.push_back(schedule_penalty(utc, 0.0, DayMetric::lden));
        day.push_back(schedule_penalty(utc, 0.0, DayMetric::lday));
        evening.push_back(schedule_penalty(utc, 0.0, DayMetric::levening));
        night.push_back(schedule_penalty(utc, 0.0, DayMetric::lnight));
    }
    std::vector<double> scratch;
    const double T = 3600.0;
    const auto level = [&](const std::vector<double>& penalty, DayMetric metric) {
        return 10.0 * std::log10(std::pow(10.0, 0.1 * penalized_exposure_level(levels, T, penalty, scratch, {}, 0.0, true))
                / averaging_seconds(metric));
    };
    const double lden = level(lden_penalty, DayMetric::lden);
    const double ld = level(day, DayMetric::lday);
    const double le = level(evening, DayMetric::levening);
    const double ln = level(night, DayMetric::lnight);
    const double recombined = 10.0 * std::log10((12.0 * std::pow(10.0, 0.1 * ld) + 4.0 * std::pow(10.0, 0.1 * (le + 5.0))
                                                        + 8.0 * std::pow(10.0, 0.1 * (ln + 10.0)))
            / 24.0);
    REQUIRE(std::fabs(lden - recombined) < 1e-9);
}
