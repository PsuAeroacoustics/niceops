#include "geodetic.hpp"

#include <cmath>

#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double DEG_TO_RAD = PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / PI;
constexpr double HALF_PI = PI / 2.0;
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

// Transverse Mercator forward, in units of the equatorial radius, with `lon`
// already relative to the reference meridian.
Eigen::Vector2d project(double lat, double lon, double meridian_at_reference)
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
    , meridian_at_reference_(meridian_arc(reference.latitude * DEG_TO_RAD))
{
}

LocalFrame LocalFrame::track_aligned(const Geodetic& reference, double bearing)
{
    LocalFrame frame(reference);
    // A compass bearing runs clockwise from north; the rotation that carries
    // east/north onto it is its complement.
    const double rotation = (90.0 - bearing) * DEG_TO_RAD;
    frame.cos_rotation_ = std::cos(rotation);
    frame.sin_rotation_ = std::sin(rotation);
    return frame;
}

Eigen::Vector3d LocalFrame::to_local(const Geodetic& position) const
{
    double lat = position.latitude * DEG_TO_RAD;
    double lon = position.longitude * DEG_TO_RAD;

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
    lon -= reference_.longitude * DEG_TO_RAD;

    const auto projected = project(lat, lon, meridian_at_reference_);
    const double east = EQUATOR_RADIUS * projected(0);
    const double north = EQUATOR_RADIUS * projected(1);
    return {east * cos_rotation_ + north * sin_rotation_,
            -east * sin_rotation_ + north * cos_rotation_, position.altitude};
}

Geodetic LocalFrame::to_geodetic(const Eigen::Vector3d& local) const
{
    // Newton on the 2x2 forward map, from the reference point.  The projection
    // is smooth and very nearly an identity near the reference, so this lands
    // on the answer in three or four steps.
    const Eigen::Vector2d target(local(0), local(1));
    Geodetic guess = reference_;
    guess.altitude = local(2);

    // About a millimetre at the surface, big enough to keep the difference
    // quotient well clear of round-off and small enough to be accurate.
    const double step = 1.0e-8;
    for (int iteration = 0; iteration < 50; iteration++)
    {
        const Eigen::Vector3d at = to_local(guess);
        const Eigen::Vector2d residual = Eigen::Vector2d(at(0), at(1)) - target;
        if (residual.norm() < 1.0e-9)
            break;

        Geodetic d_lat = guess;
        d_lat.latitude += step;
        Geodetic d_lon = guess;
        d_lon.longitude += step;
        const Eigen::Vector3d at_lat = to_local(d_lat);
        const Eigen::Vector3d at_lon = to_local(d_lon);

        Eigen::Matrix2d jacobian;
        jacobian(0, 0) = (at_lat(0) - at(0)) / step;
        jacobian(1, 0) = (at_lat(1) - at(1)) / step;
        jacobian(0, 1) = (at_lon(0) - at(0)) / step;
        jacobian(1, 1) = (at_lon(1) - at(1)) / step;

        const double determinant = jacobian.determinant();
        if (std::fabs(determinant) < 1.0e-20)
        {
            // Degenerate at the poles, where every longitude maps to one point.
            break;
        }
        const Eigen::Vector2d correction = jacobian.inverse() * residual;
        guess.latitude -= correction(0);
        guess.longitude -= correction(1);
    }
    return guess;
}

} // namespace niceops
