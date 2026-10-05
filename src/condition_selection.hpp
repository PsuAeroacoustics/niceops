#pragma once

// How a trajectory state's operating condition picks the database spheres it
// is evaluated from (--blend), and how a fixed_load_factor database's
// load-factor law is applied (--load_factor_correction).
//
// The default, barycentric, is the Delaunay blend in interpolator.hpp and is
// not here.  The two alternatives exist for comparison with other codes and
// with regulatory practice, not because they predict better:
//
//   nearest  the single stored condition nearest the state, in the same
//            normalized condition space the triangulation is built in.  ECAC
//            Doc 32 assigns one hemisphere to each constant segment of a
//            procedure; a trajectory tool has no segments, and taking the
//            nearest sphere at every state is the closest it can come.
//
//   aam      the rotorcraft sphere search of AAM v3.1 (Technical Reference and
//            User's Guide, May 2021, section 2.5.1.1, pp. 31-34): rows of
//            spheres at one flight path angle (each rounded to 0.5 degrees),
//            the two rows bracketing the state, within each the two speeds
//            bracketing it, up to four spheres blended linearly in dB, speed
//            first and then path angle.  A row pair wider than 3 degrees (the
//            state descending) or 5 (climbing or level) is not interpolated,
//            and the nearer row is taken alone; a speed pair wider than 20 kt
//            likewise gives the nearer speed.  Path angle may disqualify a
//            sphere, speed may not: the nearest row always supplies one.
//
// Everything here is a pure function of the conditions it is given, so the
// unit tests drive it on small synthetic sets without writing a database.

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace niceops {

enum class BlendPolicy
{
    barycentric,
    nearest,
    aam,
};

std::optional<BlendPolicy> parse_blend_policy(const std::string& name);
const char* to_string(BlendPolicy policy);

// Whether a fixed_load_factor database's law, 20 lg(C_T / C_T,ref), is applied
// only above the reference (one_sided, panam's convention and the default) or
// on both sides of it (two_sided, floored at a minimum load factor).
enum class LoadFactorCorrection
{
    one_sided,
    two_sided,
};

std::optional<LoadFactorCorrection> parse_load_factor_correction(const std::string& name);
const char* to_string(LoadFactorCorrection correction);

// The interpolation limits of AAM's Table 2-7, rotary wing.  Approach is a
// descending state (flight path angle below zero); level and climbing states
// take the departure limit.  Not user options: they are AAM's, and the policy
// exists to reproduce AAM.
struct AamLimits
{
    double approach_fpa_deg = 3.0;
    double departure_fpa_deg = 5.0;
    double speed_kt = 20.0;
};

// How Database::blend_at() picks and weights its spheres, and how it applies a
// fixed_load_factor database's law.  The defaults are the behavior every
// caller had before these existed.
struct BlendOptions
{
    BlendPolicy policy = BlendPolicy::barycentric;
    LoadFactorCorrection load_factor_correction = LoadFactorCorrection::one_sided;
    // The floor of the two-sided law, in g: below it the correction is held.
    // Ignored one-sided.
    double min_load_factor = 0.5;
    // Rotor tip speed, m/s, which turns AAM's 20 kt speed gate into advance
    // ratio when the database records none of its own (its own is what its
    // spheres' mu were formed with, so it is preferred).  BlendPolicy::aam
    // refuses to run with neither.
    double tip_speed = std::numeric_limits<double>::quiet_NaN();
    AamLimits aam_limits;
};

// AAM sorts and brackets on flight path angles rounded to the nearest
// 0.5 degree (Technical Reference p. 31, since RNM 7).  Halves round away from
// zero, as Fortran's NINT does.
double aam_rounded_flight_path_angle(double degrees);

// Which of an ascending list of distinct knots a query at `x` takes, by AAM's
// rule for one axis.  An exact match takes that knot.  A query between two
// knots no more than `gate` apart interpolates between them, `fraction` of the
// way from `lower` to `upper`.  Otherwise -- the pair is wider than the gate,
// or the query is past an end -- the nearest knot is taken alone; a query
// exactly halfway takes `upper` when `tie_to_upper`, else `lower`.  A single
// knot is lower == upper.
struct AxisPick
{
    std::size_t lower = 0;
    std::size_t upper = 0;
    double fraction = 0.0;
};
AxisPick aam_axis_pick(const std::vector<double>& knots, double x, double gate, bool tie_to_upper);

// AAM's search over one set of spheres, each with a speed (in any unit
// proportional to airspeed -- the database's advance ratio serves, and the
// gate is converted to it) and a flight path angle in degrees.
class AamSelector
{
public:
    // `ids` names each sphere in the caller's numbering; it is what select()
    // hands back.
    void build(const std::vector<double>& speeds, const std::vector<double>& flight_path_angles,
            const std::vector<std::size_t>& ids);

    bool empty() const
    {
        return rows_.empty();
    }
    std::size_t row_count() const
    {
        return rows_.size();
    }

    // Up to four spheres and their weights (summing to one) for a state at
    // `speed` (same unit as build()'s) and `flight_path_angle` (degrees,
    // positive climbing), with `speed_gate` the 20 kt limit in that unit.
    // Returns how many slots it filled; the rest of `weights` is zero.
    std::size_t select(double speed, double flight_path_angle, double speed_gate,
            const AamLimits& limits, std::array<double, 4>& weights,
            std::array<std::size_t, 4>& ids) const;

private:
    struct Sphere
    {
        double flight_path_angle; // as stored, not rounded
        std::size_t id;
    };
    // One distinct speed within a row.  Nearly always one sphere; two spheres
    // whose path angles round to the same row at the same speed are both
    // kept, and the one whose own path angle is nearer the state is taken.
    struct Speed
    {
        double value;
        std::vector<Sphere> spheres;
    };
    struct Row
    {
        double flight_path_angle; // rounded to 0.5 degrees
        std::vector<double> speeds;
        std::vector<Speed> entries;
    };
    std::size_t resolve(const Speed& speed, double flight_path_angle) const;

    std::vector<double> row_angles_;
    std::vector<Row> rows_;
};

// The index of the point nearest `query` among `count` points of `dim`
// coordinates, row-major, by Euclidean distance; the lowest index wins a tie.
std::size_t nearest_point(std::size_t dim, std::size_t count, const double* points,
        const double* query);

} // namespace niceops
