#include "geodetic.hpp"

#include <cmath>

#include "util/constants.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

constexpr double HALF_PI = pi / 2.0;
// Farther from the frame's origin than this is not a flight; see to_local.
constexpr double max_distance_from_origin = 1.0e6; // m
constexpr double EPSILON10 = 1.0e-10;
constexpr double EPSILON12 = 1.0e-12;

// WGS84.
constexpr double EQUATOR_RADIUS = 6378137.0;
constexpr double FLATTENING = 1.0 / 298.257223563;
constexpr double SQUARED_ECCENTRICITY = 2.0 * FLATTENING - FLATTENING * FLATTENING;

// Series coefficients for the meridional arc, as in the original.
constexpr double C00 = 1.0;
constexpr double C02 = 0.25;
constexpr double C04 = 0.046875;
constexpr double C06 = 0.01953125;
constexpr double C08 = 0.01068115234375;
constexpr double C22 = 0.75;
constexpr double C44 = 0.46875;
constexpr double C46 = 0.01302083333333333333;
constexpr double C48 = 0.00712076822916666666;
constexpr double C66 = 0.36458333333333333333;
constexpr double C68 = 0.00569661458333333333;
constexpr double C88 = 0.3076171875;

constexpr double E = SQUARED_ECCENTRICITY;
constexpr double R0 = C00 - E * (C02 + E * (C04 + E * (C06 + E * C08)));
constexpr double R1 = E * (C22 - E * (C04 + E * (C06 + E * C08)));
constexpr double R2T = E * E;
constexpr double R2 = R2T * (C44 - E * (C46 + E * C48));
constexpr double R3T = R2T * E;
constexpr double R3 = R3T * (C66 - E * C68);
constexpr double R4 = R3T * E * C88;

// Meridional arc length from the equator to `lat` (radians), in units of the
// equatorial radius.
double meridian_arc(double lat)
{
    const double sin_phi = std::sin(lat);
    const double cos_phi = std::cos(lat) * sin_phi;
    const double squared_sin_phi = sin_phi * sin_phi;
    return R0 * lat
            - cos_phi * (R1 + squared_sin_phi * (R2 + squared_sin_phi * (R3 + squared_sin_phi * R4)));
}

double parallel_radius(double sin_phi, double cos_phi)
{
    return cos_phi / std::sqrt(1.0 - SQUARED_ECCENTRICITY * sin_phi * sin_phi);
}

// American polyconic forward, in units of the equatorial radius, with `lon`
// already relative to the reference meridian.  Returns (easting, northing).
struct Projected
{
    double east;
    double north;
};

Projected project(double lat, double lon, double meridian_at_reference)
{
    if (std::fabs(lat) < EPSILON10)
    {
        // On the equator the projection degenerates to the plate carree.
        return {lon, -meridian_at_reference};
    }
    const double sin_lat = std::sin(lat);
    const double cos_lat = std::cos(lat);
    const double ms
            = std::fabs(sin_lat) > EPSILON10 ? parallel_radius(sin_lat, cos_lat) / sin_lat : 0.0;
    // NOTE: the scaled longitude is what both components use -- this is the
    // original's `sin(lon *= sin(lat))`, whose side effect carries into the
    // line below it.
    const double scaled_lon = lon * sin_lat;
    return {ms * std::sin(scaled_lon),
            (meridian_arc(lat) - meridian_at_reference) + ms * (1.0 - std::cos(scaled_lon))};
}

} // namespace

namespace niceops {

LocalFrame::LocalFrame(const Geodetic& reference)
    : reference_(reference)
    , meridian_at_reference_(meridian_arc(reference.latitude * deg_to_rad))
{
}

LocalFrame LocalFrame::track_aligned(const Geodetic& reference, double bearing)
{
    LocalFrame frame(reference);
    // A compass bearing runs clockwise from north; the rotation that carries
    // east/north onto it is its complement.
    const double rotation = (90.0 - bearing) * deg_to_rad;
    frame.cos_rotation_ = std::cos(rotation);
    frame.sin_rotation_ = std::sin(rotation);
    return frame;
}

Vec3 LocalFrame::to_local(const Geodetic& position) const
{
    double lat = position.latitude * deg_to_rad;
    double lon = position.longitude * deg_to_rad;

    const double from_pole = std::fabs(lat) - HALF_PI;
    if (from_pole > EPSILON12 || std::fabs(lon) > 10.0)
    {
        input_error("Position is outside the WGS84 projection's domain: latitude ",
                position.latitude, ", longitude ", position.longitude, ".");
    }
    if (std::fabs(from_pole) < EPSILON12)
    {
        lat = lat < 0.0 ? -HALF_PI : HALF_PI;
    }
    lon -= reference_.longitude * deg_to_rad;
    // The difference of two longitudes is an angle, not a number to project: a
    // track across the antimeridian (reference 179.99 E, point 179.995 W) came out
    // 32,000 km away.  remainder() leaves a difference already within half a turn
    // exactly as it was, so nothing that worked before moves.
    lon = std::remainder(lon, 2.0 * pi);

    const auto projected = project(lat, lon, meridian_at_reference_);
    const double east = EQUATOR_RADIUS * projected.east;
    const double north = EQUATOR_RADIUS * projected.north;
    // This is a local frame.  A position a thousand kilometers from its origin is
    // a wrong coordinate or a sign-flipped --origin, not a flight, and projecting
    // it anyway (it used to) gave a plausible-looking footprint for nowhere.
    if (std::hypot(east, north) > max_distance_from_origin)
    {
        input_error("Position (latitude ", position.latitude, ", longitude ", position.longitude,
                ") is ", std::hypot(east, north) / 1000.0, " km from the frame's origin (",
                reference_.latitude, ", ", reference_.longitude, "), too far for a local frame. "
                "Check the coordinate, and the sign of --origin.");
    }
    return {east * cos_rotation_ + north * sin_rotation_,
            -east * sin_rotation_ + north * cos_rotation_, position.altitude};
}

Vec3 LocalFrame::from_east_north_up(const Vec3& enu) const
{
    // The same rotation to_local applies to a projected position, without the
    // projection: over a footprint's few kilometers the convergence of the
    // meridians is far below the accuracy of any velocity a trajectory file
    // carries, so one rotation serves the whole frame.
    return {enu(0) * cos_rotation_ + enu(1) * sin_rotation_,
            -enu(0) * sin_rotation_ + enu(1) * cos_rotation_, enu(2)};
}

Geodetic LocalFrame::to_geodetic(const Vec3& local) const
{
    // Newton on the 2x2 forward map, from the reference point.  The projection
    // is smooth and very nearly an identity near the reference, so this lands
    // on the answer in three or four steps.
    const double target_east = local(0);
    const double target_north = local(1);
    Geodetic guess = reference_;
    guess.altitude = local(2);

    // About a millimeter at the surface, big enough to keep the difference
    // quotient well clear of round-off and small enough to be accurate.
    const double step = 1.0e-8;
    for (int iteration = 0; iteration < 50; iteration++)
    {
        const Vec3 at = to_local(guess);
        const double residual_east = at(0) - target_east;
        const double residual_north = at(1) - target_north;
        if (std::hypot(residual_east, residual_north) < 1.0e-9)
            break;

        Geodetic d_lat = guess;
        d_lat.latitude += step;
        Geodetic d_lon = guess;
        d_lon.longitude += step;
        const Vec3 at_lat = to_local(d_lat);
        const Vec3 at_lon = to_local(d_lon);

        // The 2x2 Jacobian of (easting, northing) against (latitude, longitude),
        // by finite differences, solved by hand -- Cramer's rule on two rows
        // needs no linear algebra library.
        const double d_east_d_lat = (at_lat(0) - at(0)) / step;
        const double d_north_d_lat = (at_lat(1) - at(1)) / step;
        const double d_east_d_lon = (at_lon(0) - at(0)) / step;
        const double d_north_d_lon = (at_lon(1) - at(1)) / step;

        const double determinant = d_east_d_lat * d_north_d_lon - d_east_d_lon * d_north_d_lat;
        if (std::fabs(determinant) < 1.0e-20)
        {
            // Degenerate at the poles, where every longitude maps to one point.
            break;
        }
        const double correction_lat
                = (d_north_d_lon * residual_east - d_east_d_lon * residual_north) / determinant;
        const double correction_lon
                = (d_east_d_lat * residual_north - d_north_d_lat * residual_east) / determinant;
        guess.latitude -= correction_lat;
        guess.longitude -= correction_lon;
    }
    return guess;
}

} // namespace niceops
