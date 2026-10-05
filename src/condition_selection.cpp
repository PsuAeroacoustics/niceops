#include "condition_selection.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace niceops {

std::optional<BlendPolicy> parse_blend_policy(const std::string& name)
{
    if (name == "barycentric")
        return BlendPolicy::barycentric;
    if (name == "nearest")
        return BlendPolicy::nearest;
    if (name == "aam")
        return BlendPolicy::aam;
    return std::nullopt;
}

const char* to_string(BlendPolicy policy)
{
    switch (policy)
    {
    case BlendPolicy::nearest:
        return "nearest";
    case BlendPolicy::aam:
        return "aam";
    case BlendPolicy::barycentric:
        break;
    }
    return "barycentric";
}

std::optional<LoadFactorCorrection> parse_load_factor_correction(const std::string& name)
{
    if (name == "one_sided")
        return LoadFactorCorrection::one_sided;
    if (name == "two_sided")
        return LoadFactorCorrection::two_sided;
    return std::nullopt;
}

const char* to_string(LoadFactorCorrection correction)
{
    return correction == LoadFactorCorrection::two_sided ? "two_sided" : "one_sided";
}

double aam_rounded_flight_path_angle(double degrees)
{
    // std::round takes halves away from zero, which is NINT's rule.
    return std::round(2.0 * degrees) / 2.0;
}

AxisPick aam_axis_pick(const std::vector<double>& knots, double x, double gate, bool tie_to_upper)
{
    if (knots.empty())
        input_error("AAM sphere search: an axis with no knots.");
    if (!std::isfinite(x))
        input_error("AAM sphere search: the state's condition is ", x, ", not a number it can bracket.");

    // First knot at or above x.  At or below the first knot, and past the last,
    // there is nothing to bracket with and the end knot is the nearest.
    const auto above = std::lower_bound(knots.begin(), knots.end(), x);
    if (above == knots.begin())
        return {0, 0, 0.0};
    if (above == knots.end())
        return {knots.size() - 1, knots.size() - 1, 0.0};
    const std::size_t upper = static_cast<std::size_t>(above - knots.begin());
    if (knots[upper] == x)
        return {upper, upper, 0.0};
    const std::size_t lower = upper - 1;
    const double width = knots[upper] - knots[lower];
    if (width <= gate)
        return {lower, upper, (x - knots[lower]) / width};

    // Too far apart to interpolate: the nearer one alone.
    const double to_lower = x - knots[lower];
    const double to_upper = knots[upper] - x;
    if (to_upper < to_lower || (to_upper == to_lower && tie_to_upper))
        return {upper, upper, 0.0};
    return {lower, lower, 0.0};
}

void AamSelector::build(const std::vector<double>& speeds, const std::vector<double>& flight_path_angles,
        const std::vector<std::size_t>& ids)
{
    if (speeds.size() != flight_path_angles.size() || speeds.size() != ids.size())
        input_error("AAM sphere search: ", speeds.size(), " speeds, ", flight_path_angles.size(),
                " flight path angles and ", ids.size(), " sphere ids.");
    for (std::size_t i = 0; i < speeds.size(); i++)
    {
        if (!std::isfinite(speeds[i]) || !std::isfinite(flight_path_angles[i]))
            input_error("AAM sphere search: sphere ", ids[i], " is at speed ", speeds[i],
                    " and flight path angle ", flight_path_angles[i], ".");
    }

    std::vector<double> rounded(speeds.size());
    for (std::size_t i = 0; i < speeds.size(); i++)
        rounded[i] = aam_rounded_flight_path_angle(flight_path_angles[i]);

    // Table 2-5's sort, by (rounded) path angle and then speed; the id last
    // only makes the order total.
    std::vector<std::size_t> order(speeds.size());
    std::iota(order.begin(), order.end(), std::size_t {0});
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        if (rounded[a] != rounded[b])
            return rounded[a] < rounded[b];
        if (speeds[a] != speeds[b])
            return speeds[a] < speeds[b];
        return ids[a] < ids[b];
    });

    rows_.clear();
    row_angles_.clear();
    for (const std::size_t i : order)
    {
        if (rows_.empty() || rows_.back().flight_path_angle != rounded[i])
        {
            rows_.push_back({rounded[i], {}, {}});
            row_angles_.push_back(rounded[i]);
        }
        Row& row = rows_.back();
        if (row.speeds.empty() || row.speeds.back() != speeds[i])
        {
            row.speeds.push_back(speeds[i]);
            row.entries.push_back({speeds[i], {}});
        }
        row.entries.back().spheres.push_back({flight_path_angles[i], ids[i]});
    }
}

std::size_t AamSelector::resolve(const Speed& speed, double flight_path_angle) const
{
    // AAM itself stops on two spheres at one profile ("Same profile data found
    // in multiple spheres").  Two whose angles only round alike are distinct
    // conditions, so the nearer one's is the answer; ties keep the sort order.
    std::size_t best = speed.spheres.front().id;
    double best_distance = std::fabs(speed.spheres.front().flight_path_angle - flight_path_angle);
    for (std::size_t k = 1; k < speed.spheres.size(); k++)
    {
        const double distance = std::fabs(speed.spheres[k].flight_path_angle - flight_path_angle);
        if (distance < best_distance)
        {
            best_distance = distance;
            best = speed.spheres[k].id;
        }
    }
    return best;
}

std::size_t AamSelector::select(double speed, double flight_path_angle, double speed_gate,
        const AamLimits& limits, std::array<double, 4>& weights, std::array<std::size_t, 4>& ids) const
{
    weights.fill(0.0);
    ids.fill(0);
    if (rows_.empty())
        input_error("AAM sphere search over no spheres.");

    // Path angle first: it may disqualify spheres.  A descending state is an
    // approach (3 degrees), anything else a departure (5).  The tie rule of
    // Table 2-8 is taken by its labels -- approaches the higher angle,
    // departures the lower, i.e. toward level flight in both -- because its
    // signs contradict the document's own convention (positive = departure).
    const bool approach = flight_path_angle < 0.0;
    const AxisPick rows = aam_axis_pick(row_angles_, flight_path_angle,
            approach ? limits.approach_fpa_deg : limits.departure_fpa_deg, approach);

    std::size_t filled = 0;
    const auto add_row = [&](std::size_t r, double row_weight) {
        // Speed second: it guides the choice but never disqualifies, so the
        // nearest speed in the row stands in when nothing brackets the state
        // within the gate.  A tie goes to the faster sphere, "generally a
        // conservative result".
        const Row& row = rows_[r];
        const AxisPick pick = aam_axis_pick(row.speeds, speed, speed_gate, true);
        weights[filled] = row_weight * (1.0 - pick.fraction);
        ids[filled++] = resolve(row.entries[pick.lower], flight_path_angle);
        if (pick.upper != pick.lower)
        {
            weights[filled] = row_weight * pick.fraction;
            ids[filled++] = resolve(row.entries[pick.upper], flight_path_angle);
        }
    };
    add_row(rows.lower, 1.0 - rows.fraction);
    if (rows.upper != rows.lower)
        add_row(rows.upper, rows.fraction);
    return filled;
}

std::size_t nearest_point(std::size_t dim, std::size_t count, const double* points, const double* query)
{
    std::size_t best = 0;
    double best_distance = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < count; i++)
    {
        double distance = 0.0;
        for (std::size_t j = 0; j < dim; j++)
        {
            const double delta = points[i * dim + j] - query[j];
            distance += delta * delta;
        }
        if (distance < best_distance)
        {
            best_distance = distance;
            best = i;
        }
    }
    return best;
}

} // namespace niceops
