#pragma once

#include <string>
#include <utility>
#include <vector>

#include "geodetic.hpp"
#include "vec3.hpp"

namespace niceops {

// A bare-earth elevation grid on geographic axes, such as a USGS 3DEP
// 1/3 arc-second terrain model converted to netCDF.
//
// The file holds 1-D `latitude` and `longitude` coordinate variables (degrees,
// strictly increasing) and a 2-D `elevation(latitude, longitude)` in meters.
// An optional text attribute `vertical_datum` records the datum the heights are
// on.  NICE-OPS does not convert between datums: the trajectory's altitude and
// the terrain must already share one, and a mismatch is the caller's error.
class Terrain
{
public:
    explicit Terrain(const std::string& path);

    // Height (m) at a geographic position, bilinear in the four surrounding
    // nodes.  NaN outside the grid, and NaN where a node is missing.
    double elevation(double latitude, double longitude) const;

    bool covers(double latitude, double longitude) const;

    const std::string& vertical_datum() const
    {
        return vertical_datum_;
    }

private:
    std::vector<double> latitudes_;
    std::vector<double> longitudes_;
    // Row-major: elevations_[i * longitudes_.size() + j] is (latitude i, longitude j).
    std::vector<double> elevations_;
    std::string vertical_datum_;
};

// A line fitted by least squares to terrain heights along a path:
// z = intercept + slope * s, with s the distance along the path from its start.
// `rms` is the root-mean-square residual, which is how far the real ground
// departs from its mean line.
struct GroundPlane
{
    double intercept = 0.0;
    double slope = 0.0;
    double rms = 0.0;
};

GroundPlane fit_ground_plane(const std::vector<double>& distances, const std::vector<double>& heights);

// One terrain height along a path: its distance from the start, and the height.
struct ProfileSample
{
    double distance;
    double height;
};

// Terrain heights along the straight path from `from` to `to` (local meters in
// `frame`), at most `spacing` apart, both ends included.  Positions are
// projected back to geographic coordinates to read the grid.
std::vector<ProfileSample> terrain_profile(const Terrain& terrain, const LocalFrame& frame,
        const Vec3& from, const Vec3& to, double spacing);

} // namespace niceops
