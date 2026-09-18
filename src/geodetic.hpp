#pragma once

// WGS84 geodetic coordinates and their projection onto the local Cartesian
// frame the rest of NICE-OPS works in.  Ported from the D original's
// worldtransforms.d and wgs84tocartesian.d, the latter itself a port of
// https://github.com/chrberger/WGS84toCartesian.
//
// The projection is a transverse Mercator about a reference point, with
// altitude carried straight through as z.  That makes it a map projection with
// a height attached rather than a true ENU frame: it does not rotate the
// horizontal axes with height, and it does not account for the Earth curving
// away beneath a distant observer.  Over the few kilometres a noise footprint
// spans, both effects are far below the modelling error, and it is the
// projection the shipped trajectory files were built with -- their x/y columns
// agree with it to a few centimetres.

#include <string>

#include <Eigen/Dense>

namespace niceops {

// A WGS84 position.  Angles in degrees, altitude in metres.
struct Geodetic
{
    double latitude = 0.0;
    double longitude = 0.0;
    double altitude = 0.0;
};

// Projects WGS84 positions onto a local Cartesian frame centred on a reference
// point: x east, y north, z up, all in metres.
class LocalFrame
{
public:
    LocalFrame() = default;

    // x east, y north.
    explicit LocalFrame(const Geodetic& reference);

    // A frame whose +x axis points along `bearing` (degrees clockwise from
    // north), +y ninety degrees to its left.  The shipped AC trajectories use
    // one of these rather than east/north: their x/y run along the flight
    // track, rotated by about the heading in the file's own heading column.
    // Geodetic input has to be projected into the same frame as everything
    // else in a run, so the rotation is part of the transform, not an
    // afterthought -- the D original's WorldTransform had no notion of it and
    // produced east/north only.
    static LocalFrame track_aligned(const Geodetic& reference, double bearing);

    const Geodetic& reference() const
    {
        return reference_;
    }

    // Geodetic -> local Cartesian metres.  Throws for a position outside the
    // projection's domain; the D original quietly returned the origin instead,
    // which silently teleports a bad coordinate to the reference point.
    Eigen::Vector3d to_local(const Geodetic& position) const;

    // ...and back.  Newton iteration on the forward projection, which converges
    // to the last bit in a handful of steps.  The D original walked latitude
    // and then longitude in fixed 1e-5 degree increments until the error
    // stopped shrinking, stopping a centimetre out after thousands of forward
    // evaluations.
    Geodetic to_geodetic(const Eigen::Vector3d& local) const;

private:
    Geodetic reference_;
    // Rotation from east/north into the local axes: cos and sin of the angle
    // the frame is turned through.  Identity for a plain east/north frame.
    double cos_rotation_ = 1.0;
    double sin_rotation_ = 0.0;
    // Meridional arc length at the reference latitude, the projection's origin
    // offset.  Depends only on the reference, so it is computed once.
    double meridian_at_reference_ = 0.0;
};

} // namespace niceops
