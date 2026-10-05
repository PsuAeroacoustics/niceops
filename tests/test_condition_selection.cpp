// --blend's sphere selection on small synthetic condition sets
// (src/condition_selection.hpp): AAM's bracketing search and its gates, ties
// and fallbacks, and the nearest-condition search.  The Database-level cases --
// the policies wired into blend_at(), the two-sided load-factor law -- are in
// test_sphere_database.cpp, next to the fixture that writes databases.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

#include "condition_selection.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace)

namespace {

// select()'s output as {id: weight}, so a case can say which spheres it
// expects without depending on the slot order.
std::map<std::size_t, double> picked(const AamSelector& selector, double speed, double fpa,
        double speed_gate = 20.0, const AamLimits& limits = {})
{
    std::array<double, 4> weights {};
    std::array<std::size_t, 4> ids {};
    const std::size_t filled = selector.select(speed, fpa, speed_gate, limits, weights, ids);
    REQUIRE(filled >= 1);
    REQUIRE(filled <= 4);
    std::map<std::size_t, double> out;
    double total = 0.0;
    for (std::size_t k = 0; k < filled; k++)
    {
        REQUIRE(out.count(ids[k]) == 0);
        out[ids[k]] = weights[k];
        total += weights[k];
    }
    for (std::size_t k = filled; k < 4; k++)
        REQUIRE(weights[k] == 0.0);
    REQUIRE(total == Catch::Approx(1.0).epsilon(1e-15));
    return out;
}

AamSelector selector_of(const std::vector<std::array<double, 2>>& speed_fpa)
{
    std::vector<double> speeds, fpas;
    std::vector<std::size_t> ids;
    for (std::size_t i = 0; i < speed_fpa.size(); i++)
    {
        speeds.push_back(speed_fpa[i][0]);
        fpas.push_back(speed_fpa[i][1]);
        ids.push_back(i);
    }
    AamSelector selector;
    selector.build(speeds, fpas, ids);
    return selector;
}

} // namespace

TEST_CASE("--blend and --load_factor_correction names round-trip", "[condition_selection]")
{
    for (const auto policy : {BlendPolicy::barycentric, BlendPolicy::nearest, BlendPolicy::aam})
        REQUIRE(parse_blend_policy(to_string(policy)) == policy);
    REQUIRE_FALSE(parse_blend_policy("AAM").has_value());
    for (const auto correction : {LoadFactorCorrection::one_sided, LoadFactorCorrection::two_sided})
        REQUIRE(parse_load_factor_correction(to_string(correction)) == correction);
    REQUIRE_FALSE(parse_load_factor_correction("both").has_value());
}

TEST_CASE("AAM rounds sphere path angles to the half degree, halves away from zero",
        "[condition_selection]")
{
    REQUIRE(aam_rounded_flight_path_angle(-6.2) == -6.0);
    REQUIRE(aam_rounded_flight_path_angle(-6.3) == -6.5);
    REQUIRE(aam_rounded_flight_path_angle(-6.25) == -6.5);
    REQUIRE(aam_rounded_flight_path_angle(6.25) == 6.5);
    REQUIRE(aam_rounded_flight_path_angle(0.24) == 0.0);
    REQUIRE(aam_rounded_flight_path_angle(3.0) == 3.0);
}

TEST_CASE("one AAM axis: exact, gated bracket, too wide, past the ends, ties", "[condition_selection]")
{
    const std::vector<double> knots = {0.0, 10.0, 40.0};

    // Exactly on a knot: that knot alone.
    auto p = aam_axis_pick(knots, 10.0, 20.0, true);
    REQUIRE((p.lower == 1 && p.upper == 1));

    // Within the gate: a bracket and the fraction across it.
    p = aam_axis_pick(knots, 2.5, 20.0, true);
    REQUIRE((p.lower == 0 && p.upper == 1));
    REQUIRE(p.fraction == 0.25);

    // A pair exactly as wide as the gate still interpolates: the limit is a
    // maximum (Table 2-6, "dFPA deg max", "dSPD kts max").
    p = aam_axis_pick(knots, 5.0, 10.0, false);
    REQUIRE((p.lower == 0 && p.upper == 1));
    REQUIRE(p.fraction == 0.5);
    // ... and a hair narrower gate does not.
    p = aam_axis_pick(knots, 4.0, 9.999, false);
    REQUIRE((p.lower == 0 && p.upper == 0));

    // Too wide: the nearer knot, whichever side.
    p = aam_axis_pick(knots, 20.0, 20.0, true);
    REQUIRE((p.lower == 1 && p.upper == 1));
    p = aam_axis_pick(knots, 30.0, 20.0, false);
    REQUIRE((p.lower == 2 && p.upper == 2));

    // Exactly halfway across a pair too wide to interpolate: the tie rule.
    p = aam_axis_pick(knots, 25.0, 20.0, true);
    REQUIRE((p.lower == 2 && p.upper == 2));
    p = aam_axis_pick(knots, 25.0, 20.0, false);
    REQUIRE((p.lower == 1 && p.upper == 1));

    // Past either end: the end knot, never an extrapolation.
    p = aam_axis_pick(knots, -5.0, 1e9, true);
    REQUIRE((p.lower == 0 && p.upper == 0));
    p = aam_axis_pick(knots, 45.0, 1e9, false);
    REQUIRE((p.lower == 2 && p.upper == 2));

    // One knot answers everything.
    p = aam_axis_pick({7.0}, -100.0, 0.0, false);
    REQUIRE((p.lower == 0 && p.upper == 0));

    REQUIRE_THROWS(aam_axis_pick({}, 1.0, 1.0, true));
    REQUIRE_THROWS(aam_axis_pick(knots, std::numeric_limits<double>::quiet_NaN(), 1.0, true));
}

TEST_CASE("AAM brackets in speed within each path-angle row, then in path angle",
        "[condition_selection]")
{
    // Two rows with different speeds, as a scattered database has them: each
    // row brackets the state's speed on its own (Figure 2-20, I1/I2 at FPAng1
    // and I3/I4 at FPAng2), and the rows are then blended in path angle.
    const auto selector = selector_of({
            {50.0, -6.0}, // 0
            {70.0, -6.0}, // 1
            {60.0, -3.0}, // 2
            {80.0, -3.0}, // 3
            {65.0, 0.0},  // 4
    });
    REQUIRE(selector.row_count() == 3);
    const auto w = picked(selector, 65.0, -4.5);
    REQUIRE(w.size() == 4);
    REQUIRE(w.at(0) == Catch::Approx(0.5 * 0.25));
    REQUIRE(w.at(1) == Catch::Approx(0.5 * 0.75));
    REQUIRE(w.at(2) == Catch::Approx(0.5 * 0.75));
    REQUIRE(w.at(3) == Catch::Approx(0.5 * 0.25));

    // On a row's angle: no path-angle interpolation, two spheres.
    const auto on_row = picked(selector, 55.0, -6.0);
    REQUIRE(on_row.size() == 2);
    REQUIRE(on_row.at(0) == Catch::Approx(0.75));
    REQUIRE(on_row.at(1) == Catch::Approx(0.25));

    // On a sphere: that sphere whole.
    const auto exact = picked(selector, 80.0, -3.0);
    REQUIRE(exact.size() == 1);
    REQUIRE(exact.at(3) == 1.0);
}

TEST_CASE("AAM's path-angle gate is 3 degrees descending and 5 otherwise", "[condition_selection]")
{
    // Rows 4 degrees apart on each side of level: inside the departure gate,
    // outside the approach one.
    const auto selector = selector_of({
            {60.0, -5.0}, // 0
            {60.0, -1.0}, // 1
            {60.0, 1.0},  // 2
            {60.0, 5.0},  // 3
    });
    // Climbing at 3 degrees: interpolated between +1 and +5.
    const auto climb = picked(selector, 60.0, 3.0);
    REQUIRE(climb.size() == 2);
    REQUIRE(climb.at(2) == Catch::Approx(0.5));
    REQUIRE(climb.at(3) == Catch::Approx(0.5));

    // Descending at 2 degrees: the -5/-1 pair is wider than 3, so the nearer
    // row alone.
    const auto descent = picked(selector, 60.0, -2.0);
    REQUIRE(descent.size() == 1);
    REQUIRE(descent.count(1) == 1);
    const auto steeper = picked(selector, 60.0, -4.0);
    REQUIRE(steeper.size() == 1);
    REQUIRE(steeper.count(0) == 1);

    // Level flight is a departure: between rows 2 degrees either side of it
    // (4 apart, past the approach gate) a level state interpolates, and a
    // state a hair into the descent takes the nearer row alone.
    const auto straddle = selector_of({{60.0, -2.0}, {60.0, 2.0}});
    const auto level = picked(straddle, 60.0, 0.0);
    REQUIRE(level.size() == 2);
    REQUIRE(level.at(0) == Catch::Approx(0.5));
    REQUIRE(level.at(1) == Catch::Approx(0.5));
    const auto just_down = picked(straddle, 60.0, -1e-6);
    REQUIRE(just_down.size() == 1);
    REQUIRE(just_down.count(0) == 1);

    // The limits are parameters; widened, the descent interpolates too.
    AamLimits wide;
    wide.approach_fpa_deg = 4.0;
    const auto widened = picked(selector, 60.0, -2.0, 20.0, wide);
    REQUIRE(widened.size() == 2);
    REQUIRE(widened.at(0) == Catch::Approx(0.25));
    REQUIRE(widened.at(1) == Catch::Approx(0.75));
}

TEST_CASE("AAM path-angle ties go toward level flight, by Table 2-8's labels", "[condition_selection]")
{
    // Exactly halfway between two rows too far apart to interpolate.  Table
    // 2-8: approaches take the higher angle, departures the lower.
    const auto descending = selector_of({{60.0, -8.0}, {60.0, -2.0}});
    const auto approach = picked(descending, 60.0, -5.0);
    REQUIRE(approach.size() == 1);
    REQUIRE(approach.count(1) == 1);

    const auto climbing = selector_of({{60.0, 2.0}, {60.0, 10.0}});
    const auto departure = picked(climbing, 60.0, 6.0);
    REQUIRE(departure.size() == 1);
    REQUIRE(departure.count(0) == 1);
}

TEST_CASE("AAM's speed never disqualifies: the nearest speed, ties to the faster", "[condition_selection]")
{
    const auto selector = selector_of({{40.0, 0.0}, {100.0, 0.0}});
    // 60 kt: the 40/100 pair is 60 kt wide, past the 20 kt gate, so the nearer
    // speed alone -- however far it is.
    auto w = picked(selector, 60.0, 0.0);
    REQUIRE(w.size() == 1);
    REQUIRE(w.count(0) == 1);
    // Exactly halfway: the faster sphere, "generally a conservative result".
    w = picked(selector, 70.0, 0.0);
    REQUIRE(w.size() == 1);
    REQUIRE(w.count(1) == 1);
    // Past the fastest: the fastest.
    w = picked(selector, 150.0, 0.0);
    REQUIRE(w.size() == 1);
    REQUIRE(w.count(1) == 1);
    // The gate is in the caller's unit: as advance ratio for a 200 m/s tip
    // speed, 20 kt is 0.0514 -- a pair 0.06 apart is past it.
    const auto mu = selector_of({{0.20, 0.0}, {0.26, 0.0}});
    const double gate = 20.0 * 1852.0 / 3600.0 / 200.0;
    REQUIRE(picked(mu, 0.22, 0.0, gate).size() == 1);
    REQUIRE(picked(mu, 0.22, 0.0, 0.07).size() == 2);
}

TEST_CASE("AAM falls back to the nearest row, then the nearest speed in it", "[condition_selection]")
{
    // A database too sparse to bracket anything within the gates: a state
    // takes the row nearest in path angle (path angle selection comes first
    // and may disqualify), then the speed nearest in that row -- not the
    // sphere nearest in some combined distance.
    const auto selector = selector_of({
            {60.0, -12.0}, // 0
            {120.0, -12.0}, // 1
            {90.0, 0.0},   // 2
    });
    // -9 degrees is nearer -12 than 0; in that row 80 kt is nearer 60 than
    // 120, though sphere 2 at 90 kt is nearer in speed.
    const auto w = picked(selector, 80.0, -9.0);
    REQUIRE(w.size() == 1);
    REQUIRE(w.count(0) == 1);
    // Below the lowest row and past every speed: still one sphere.
    const auto far = picked(selector, 500.0, -40.0);
    REQUIRE(far.size() == 1);
    REQUIRE(far.count(1) == 1);
}

TEST_CASE("AAM rows are rounded angles; duplicates in a row take the nearer stored angle",
        "[condition_selection]")
{
    // -6.2 and -5.8 both round to -6: one row, and at one speed two spheres
    // AAM would stop on.  The one whose own angle is nearer the state wins.
    const auto selector = selector_of({
            {60.0, -6.2}, // 0
            {60.0, -5.8}, // 1
            {80.0, -6.1}, // 2
    });
    REQUIRE(selector.row_count() == 1);
    auto w = picked(selector, 60.0, -5.9);
    REQUIRE(w.size() == 1);
    REQUIRE(w.count(1) == 1);
    w = picked(selector, 60.0, -6.3);
    REQUIRE(w.count(0) == 1);
    // The row's angle is the rounded one: a state at -6 is on the row, a
    // state at -6.4 below every row, and both take it unblended in angle.
    w = picked(selector, 70.0, -6.4);
    REQUIRE(w.size() == 2);
    REQUIRE(w.at(0) == Catch::Approx(0.5));
    REQUIRE(w.at(2) == Catch::Approx(0.5));
}

TEST_CASE("AamSelector refuses mismatched or non-finite input", "[condition_selection]")
{
    AamSelector selector;
    REQUIRE_THROWS(selector.build({1.0, 2.0}, {0.0}, {0, 1}));
    REQUIRE_THROWS(selector.build({1.0}, {std::numeric_limits<double>::quiet_NaN()}, {0}));
    std::array<double, 4> weights {};
    std::array<std::size_t, 4> ids {};
    REQUIRE_THROWS(selector.select(1.0, 0.0, 20.0, {}, weights, ids));
}

TEST_CASE("nearest_point takes the nearest, the lowest index on a tie", "[condition_selection]")
{
    const std::vector<double> points = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0};
    const double near_second[2] = {0.9, 0.2};
    REQUIRE(nearest_point(2, 4, points.data(), near_second) == 1);
    const double halfway[2] = {0.5, 0.0}; // equidistant from 0, 1 and 3
    REQUIRE(nearest_point(2, 4, points.data(), halfway) == 0);
    const double near_third[2] = {0.1, 0.8};
    REQUIRE(nearest_point(2, 4, points.data(), near_third) == 2);
}
