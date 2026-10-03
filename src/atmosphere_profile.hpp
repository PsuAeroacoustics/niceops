#pragma once

// A horizontally stratified atmosphere -- temperature, humidity and pressure
// against height -- and the layered absorption it implies.
//
// ISO 9613-1's coefficient depends on all three, and between a helicopter at a
// few hundred metres and the ground they are not one value.  Integrating the
// coefficient along the path rather than taking it at the ground moved LAE by
// -0.6 to +1.3 dB(A) at 500 m in sonAIR's comparison (Zellmann & Wunderli,
// Inter-Noise 2014).  It is the first step of the propagation roadmap because
// it needs no ray model: along a straight path, what a band loses is the
// coefficient's mean over the heights the path spans, times the path length.
//
// Heights are above the observer plane -- the same z the trajectory uses -- so
// the ground is z = 0 whatever --ground_elevation says.
//
// Wind, when the file gives it, only matters to refraction (--refraction):
// absorption does not depend on it.

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "atmosphere.hpp"
#include "util/constants.hpp"

namespace niceops {

struct AtmosphereProfile
{
    // The file it came from, for messages.
    std::string source;
    // Metres above the observer plane, strictly ascending.
    std::vector<double> height;
    std::vector<double> temperature; // K
    std::vector<double> relative_humidity; // percent
    // kPa.  Empty when the file gave no pressure; with_hydrostatic_pressure()
    // fills it.
    std::vector<double> pressure;
    // m/s toward the east and toward the north (the meteorological u and v).
    // Empty when the file gives no wind: still air.
    std::vector<double> wind_east;
    std::vector<double> wind_north;

    bool has_pressure() const
    {
        return !pressure.empty();
    }
    bool has_wind() const
    {
        return !wind_east.empty();
    }
    double bottom() const
    {
        return height.front();
    }
    double top() const
    {
        return height.back();
    }

    // The atmosphere at height z, in metres.  Temperature and humidity are
    // linear between levels.  Pressure is linear in its logarithm, since it
    // falls off exponentially.  Every quantity is held at its end value outside
    // the levels.  Needs has_pressure().
    Atmosphere at(double z) const;

    // The wind at height z, m/s east and north, interpolated like temperature.
    // Zero without wind columns.
    std::array<double, 2> wind_at(double z) const;

    // A copy with the pressure filled in: the hydrostatic equation integrated
    // through the temperature profile from `surface_pressure` kPa at z = 0,
    // treating the air as dry.  Humidity would raise the virtual temperature,
    // which moves the pressure at 1 km by well under 0.1%, and absorption is
    // the least sensitive of the three inputs to pressure.
    AtmosphereProfile with_hydrostatic_pressure(double surface_pressure) const;
};

// Turns a wind's east and north components into the local frame, whose +x
// axis points along compass bearing `frame_bearing` degrees and +y ninety
// degrees to its left, as LocalFrame::track_aligned has it.  The ray models and
// the air-referenced flight conditions all take the wind through this one
// rotation.
struct CompassToLocal
{
    explicit CompassToLocal(double frame_bearing)
        : sin_b(std::sin(frame_bearing * pi / 180.0))
        , cos_b(std::cos(frame_bearing * pi / 180.0))
    {
    }

    std::array<double, 2> operator()(const std::array<double, 2>& east_north) const
    {
        return {east_north[0] * sin_b + east_north[1] * cos_b,
                -east_north[0] * cos_b + east_north[1] * sin_b};
    }

    // The inverse: a vector in the local frame as (east, north).
    std::array<double, 2> east_north(const std::array<double, 2>& local) const
    {
        return {local[0] * sin_b - local[1] * cos_b, local[0] * cos_b + local[1] * sin_b};
    }

    double sin_b;
    double cos_b;
};

// Read a profile.  Columns are found by name (read_labelled_csv), in any order,
// one row per level:
//   height above the observers   z, agl                   in `unit`; z_ft / z_m
//                                                          are checked against it
//   or altitude above sea level  alt, altitude, height    in `unit`, less
//                                                          `ground_elevation` (m)
//   temperature                  T_K, T_C                 or T / temperature,
//                                                          unit told from values
//   relative humidity, percent   RH, humidity
//   pressure (optional)          p_kPa, p_hPa, p_Pa       or p / pressure,
//                                                          unit told from values
//   wind (optional), either      wind_speed_mps / _kt     and wind_direction, the
//                                                          compass bearing it
//                                                          blows FROM, degrees
//                                or u_mps / u_kt, v_...   toward east and north
// Wind columns must name their unit: nothing about the values tells m/s from
// knots.  Other columns are ignored and named in the log.  A single row is a uniform
// atmosphere.  Rows may come in either height order.
AtmosphereProfile read_atmosphere_profile(const std::string& path, LengthUnit unit,
        double ground_elevation);

// Per-band absorption integrated through a profile, tabulated in height so the
// observer loop never evaluates ISO 9613-1.
//
// The table holds alpha_b(z) at nodes a metre apart and its running integral
// F_b(z).  Between nodes alpha is linear and F its exact integral, a
// quadratic.  That is why (F_b(z_a) - F_b(z_b)) / (z_a - z_b) tends smoothly to
// alpha_b as the two heights meet, rather than jumping at a node.
//
// The observer loop needs F at a state's height and at the observer's.  Both
// are fixed while the other varies, so they are looked up once per state and
// once per observer.  Per pair, a band then costs a subtraction and a multiply.
// On a grid, where every observer shares a height, the caller forms each
// state's mean once and pays nothing per pair.
//
// Straight paths only.  For a curved ray the same mean is a first-order
// approximation to the integral along the arc.  The mean is taken over the
// whole height span, including the stretch inside the sphere radius that the
// path length excludes.  For the shipped databases that stretch is a few
// metres against hundreds.
class LayeredAbsorption
{
public:
    LayeredAbsorption() = default;
    // Tabulate over heights [z_low, z_high] m, for bands [first_band, last_band)
    // of `frequencies`.  The profile must carry pressure.
    LayeredAbsorption(const AtmosphereProfile& profile, const std::vector<double>& frequencies,
            std::size_t first_band, std::size_t last_band, double z_low, double z_high);

    // F_b(z), dB: the absorption integral from z_low up to z.  Written to
    // out[band] for each band in range; entries outside it are left alone.
    void integral(double z, double* out) const;

    // alpha_b(z), dB/m, written like integral().
    void coefficients(double z, double* out) const;

    // The mean coefficient over the heights between z_a and z_b, given
    // integral() at each: (F_a - F_b) / (z_a - z_b).  Within a millimetre of
    // each other it takes the coefficient at their midpoint instead, where the
    // difference quotient would lose its digits.  Written like integral().
    void mean_coefficients(double z_a, const double* integral_a, double z_b,
            const double* integral_b, double* out) const;

    bool empty() const
    {
        return alpha_.empty();
    }
    std::size_t first_band() const
    {
        return first_band_;
    }
    std::size_t last_band() const
    {
        return last_band_;
    }
    double z_low() const
    {
        return z_low_;
    }
    double z_high() const
    {
        return z_low_ + step_ * static_cast<double>(nodes_ - 1);
    }

private:
    // The node below z, and how far toward the next one z lies (0 to 1 inside
    // the table; outside it the end cell extrapolates, which covers round-off
    // at the edges).
    void locate(double z, std::size_t& node, double& fraction) const;

    double z_low_ = 0.0;
    double step_ = 1.0;
    double inv_step_ = 1.0;
    std::size_t nodes_ = 0;
    std::size_t first_band_ = 0;
    std::size_t last_band_ = 0;
    std::size_t bands_ = 0;
    // [node][band - first_band_]
    std::vector<double> alpha_;
    std::vector<double> integral_;
};

} // namespace niceops
