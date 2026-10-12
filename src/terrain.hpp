#pragma once

#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "geodetic.hpp"
#include "vec3.hpp"

namespace niceops {

// A bare-earth elevation grid on geographic axes, such as a USGS 3DEP
// 1/3 arc-second terrain model converted to netCDF (the validation harness's
// terrain/dem_to_netcdf.py writes it).
//
// The file holds 1-D `latitude` and `longitude` coordinate variables (degrees,
// strictly increasing) and a 2-D `elevation(latitude, longitude)` in meters.
// An optional text attribute `vertical_datum` records the datum the heights are
// on.  NICE-OPS does not convert between datums: --terrain_offset is what puts
// the heights on the trajectory's datum (at Amedee, 3DEP's NAVD88 against the
// 2017 tracks' ellipsoidal heights, -23.6 m), and a mismatch is the caller's
// error -- which run() and process_trajectory() check for as far as the inputs
// allow.
class Terrain
{
public:
    explicit Terrain(const std::string& path);

    // The same grid from memory, for tests: `elevations` row-major,
    // [i * longitudes.size() + j] at (latitude i, longitude j).  Checked as the
    // file is.
    Terrain(std::vector<double> latitudes, std::vector<double> longitudes,
            std::vector<double> elevations, std::string vertical_datum = "unstated");

    // Height (m) at a geographic position, bilinear in the four surrounding
    // nodes.  NaN outside the grid, and NaN where a node is missing.
    double elevation(double latitude, double longitude) const;

    bool covers(double latitude, double longitude) const;

    const std::string& vertical_datum() const
    {
        return vertical_datum_;
    }

    // The grid's extent, degrees.
    double south() const
    {
        return latitudes_.front();
    }
    double north() const
    {
        return latitudes_.back();
    }
    double west() const
    {
        return longitudes_.front();
    }
    double east() const
    {
        return longitudes_.back();
    }

    // The grid's node spacing in meters at its center: the smaller of the
    // mean north-south and the mean east-west spacing.  What --terrain_step
    // defaults to (floored at min_terrain_step).
    double spacing_meters() const;

private:
    void validate();

    std::vector<double> latitudes_;
    std::vector<double> longitudes_;
    // Row-major: elevations_[i * longitudes_.size() + j] is (latitude i, longitude j).
    std::vector<double> elevations_;
    std::string vertical_datum_;
};

// --terrain_step's floor when it defaults to the terrain's own spacing, m: a
// 1 m lidar grid would otherwise put three thousand samples on every 3 km ray
// for no change in a mean plane fitted over the whole of it.
constexpr double min_terrain_step = 5.0;

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
// projected back to geographic coordinates to read the grid.  The reference the
// run's faster path (TerrainRaster, GroundPlaneCache) is tested against.
std::vector<ProfileSample> terrain_profile(const Terrain& terrain, const LocalFrame& frame,
        const Vec3& from, const Vec3& to, double spacing);

// The terrain resampled onto the run's local frame: z (m, the frame's own z,
// whose zero is --ground_elevation on the trajectory's datum) at the nodes of a
// square grid, bilinear between them.  The observer loop reads terrain millions
// of times, and each geographic read costs an inverse projection (a Newton
// iteration); this costs one per node, once per trajectory.
class TerrainRaster
{
public:
    TerrainRaster() = default;

    // Samples `terrain` every `step` meters over [x_min, x_max] x [y_min, y_max]
    // of `frame`, adding `height_offset` (m) to every height: --terrain_offset
    // less --ground_elevation.  Refuses (input_error) a node the terrain does not
    // cover or has no height for, naming the extent it would need.
    TerrainRaster(const Terrain& terrain, const LocalFrame& frame, double x_min, double x_max,
            double y_min, double y_max, double step, double height_offset);

    // A raster from its heights directly, for tests: node (i, j) at
    // (x0 + i step, y0 + j step), heights[j * nx + i].
    TerrainRaster(double x0, double y0, double step, std::size_t nx, std::size_t ny,
            std::vector<double> heights);

    // Bilinear in the four surrounding nodes; a point past the edge takes the
    // edge's height (the raster is built with a margin for that, see
    // GroundPlaneCache).
    double height(double x, double y) const;

    double step() const
    {
        return step_;
    }
    bool empty() const
    {
        return heights_.empty();
    }
    std::size_t nodes() const
    {
        return heights_.size();
    }

private:
    double x0_ = 0.0;
    double y0_ = 0.0;
    double step_ = 1.0;
    double inverse_step_ = 1.0;
    std::size_t nx_ = 0;
    std::size_t ny_ = 0;
    // Row-major in y: heights_[j * nx_ + i] is (x0 + i step, y0 + j step).
    std::vector<double> heights_;
};

// The mean ground under a path, in the vertical plane through it:
// z = intercept + slope * s, s the horizontal distance from the observer toward
// the source's footprint.  `intercept` is the line's height at the observer.
struct GroundLine
{
    double intercept = 0.0;
    double slope = 0.0;
};

// The least-squares line through the raster's heights at n + 1 equally spaced
// points from (x0, y0) to (x1, y1), both ends included: the direct, per-path
// fit GroundPlaneCache approximates.  n >= 1.
GroundLine fit_ground_line(const TerrainRaster& raster, double x0, double y0, double x1, double y1,
        std::size_t n);

// The mean ground line under every path from one observer, for the observer
// loop (--terrain).  A direct fit (fit_ground_line) samples the whole profile
// for every pair -- 300 samples for a 3 km path, a microsecond, 94 million
// times on the demo grid.  Here the profiles are taken along `bearings` fixed
// rays from the observer, each sampled every raster step and kept as running
// sums of z and k z, from which the least-squares line over any whole number of
// steps is a few flops.  A ray is sampled only when, and only as far as, a pair
// first needs it, so an observer pays for the bearings and ranges its track
// actually covers.  A path between two rays, or ending between two samples,
// takes the line interpolated linearly in bearing and in length between the
// four fits around it -- continuous everywhere, and exact on a uniformly
// sloping plane up to the cosine's curvature between rays (a slope of 0.05
// errs by 1e-7, 0.4 mm at 3 km).  A path shorter than one step takes the fit
// over the first step along its bearing.
//
// One per thread: set_observer() rewinds it, and it holds no locks.
//
// With `screening` (--terrain_screening) each ray also keeps its raw heights
// and its horizon from the microphone, the running maximum of
// (z_k - h_mic) / (k step): about 16 bytes more per sample, and nothing at all
// without it.  blocked() asks whether any sample between a pair's ends rises
// above its line of sight, on the profile blended between the two rays in
// bearing as the line is; fit_range() is the least-squares line over part of
// that profile, from the same running sums.
class GroundPlaneCache
{
public:
    // 0.25 degrees: rays 13 m apart at 3 km, about a 3DEP 1/3 arc-second cell.
    static constexpr int bearings = 1440;

    explicit GroundPlaneCache(const TerrainRaster& raster, bool screening = false);

    // Start on a new observer, whose ground point is (x, y).
    void set_observer(double x, double y);
    // The same, with the microphone `microphone_height` above the raster's
    // ground at (x, y), which the horizon is taken from (screening only).
    void set_observer(double x, double y, double microphone_height);
    // The same, with the microphone at `microphone_z` in the raster's frame
    // (an observer file's own z), held no lower than the ground at its foot.
    void set_observer_at(double x, double y, double microphone_z);

    // Where a path to the footprint (dx, dy) away horizontally falls among the
    // rays: the ray at or counterclockwise below its bearing, the next one,
    // the fraction of the way between them, and the whole steps and fraction
    // of a step along them.  `distance` = hypot(dx, dy).
    struct RayPair
    {
        int ray = 0;
        int next = 0;
        double turn = 0.0;
        std::size_t n = 1;
        double fraction = 0.0;
        double distance = 0.0;
    };
    RayPair locate(double dx, double dy, double distance) const;

    // The line under the path from the observer to the footprint (dx, dy) away
    // horizontally, `distance` = hypot(dx, dy).  line_from(locate(...)) is the
    // same line to the bit.
    GroundLine line_to(double dx, double dy, double distance);
    GroundLine line_from(const RayPair& pair);

    // A pair's terrain profile, blended between its two rays: heights relative
    // to the observer's ground at samples k = 0 .. interior + 1, k steps out,
    // where 1 .. interior are the samples strictly between the microphone and
    // the source's footprint, and interior + 1 is the first at or past it.
    // `source_height` is the source's height on the same base, raised onto
    // the profile at the footprint if it is below it.
    struct Profile
    {
        std::vector<double> heights;
        std::size_t interior = 0;
        double step = 0.0;
        double distance = 0.0;
        double microphone_height = 0.0;
        double source_height = 0.0;
    };

    // Screening only.  The quick test: false when neither ray's horizon over
    // the interior samples rises above the pair's line of sight to a source at
    // local z `source_z`, which proves the blended profile does not either.
    bool may_block(const RayPair& pair, double source_z);
    // Whether some interior sample of the blended profile lies strictly above
    // the line of sight; `profile` is filled when it does.  A sample exactly on
    // it does not block.
    bool blocked(const RayPair& pair, double source_z, Profile& profile);

    // Screening only: the least-squares line through the blended profile's
    // samples k0 .. k1 (k0 < k1, both already sampled), heights relative to
    // the observer's ground, s from the observer.
    GroundLine fit_range(const RayPair& pair, std::size_t k0, std::size_t k1) const;

    // The raster's height at the observer's ground point, which the profile's
    // heights are relative to.
    double ground_z() const
    {
        return z_;
    }

private:
    // The fit over the first n steps of ray `ray` (n + 1 samples), relative to
    // the observer's height.
    GroundLine fit(int ray, std::size_t n);
    // The fit over samples k0 .. k1 of ray `ray`, relative to the observer's
    // height, s from the observer.
    GroundLine fit_ray_range(int ray, std::size_t k0, std::size_t k1) const;
    // The ray's running sums through sample n, sampling it further if needed.
    void ensure(int ray, std::size_t n);
    GroundLine along(int ray, std::size_t n, double fraction);
    // Screening: the samples strictly inside a pair's distance, and the
    // source's height raised onto the blended profile at the footprint.
    std::size_t interior_samples(const RayPair& pair) const;
    double source_height(const RayPair& pair, double source_z) const;
    // The slope of a sample's height `z` (relative to the observer's ground)
    // as seen from the microphone, `distance` away: the one spelling the
    // horizon, the line of sight and the scan all use.
    double sight_slope(double z, double distance) const
    {
        return (z - microphone_) / distance;
    }

    const TerrainRaster* raster_;
    double step_;
    bool screening_ = false;
    double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;
    double microphone_ = 0.0;
    unsigned generation_ = 0;
    std::vector<double> cos_;
    std::vector<double> sin_;
    // Per ray: the generation its sums belong to, and the sums themselves --
    // sum_z[k] = z_0 + ... + z_k, sum_kz[k] = 0 z_0 + ... + k z_k, heights
    // relative to the observer's.
    std::vector<unsigned> stamp_;
    std::vector<std::vector<double>> sum_z_;
    std::vector<std::vector<double>> sum_kz_;
    // Screening only: each sample's height relative to the observer's, and the
    // horizon, the largest sight_slope of samples 1 .. k (-infinity at k = 0).
    std::vector<std::vector<double>> z_raw_;
    std::vector<std::vector<double>> horizon_;
};

// What the ground models take, over a sloping mean ground line, for a source
// `distance` meters from the observer horizontally at height `source_z` and a
// microphone over the observer at `microphone_z` (both the local frame's z):
// the source's and the microphone's heights above the line, measured
// perpendicular to it, and the distance between their feet along it.  A
// rotation of the vertical plane, so the source-microphone distance is
// unchanged; on a level line at height zero it is exactly (source_z,
// microphone_z, distance), bit for bit.  Heights may come out negative (below
// the line); the caller decides what that means.
struct PlaneGeometry
{
    double source_height = 0.0;
    double receiver_height = 0.0;
    double ground_distance = 0.0;
};

inline PlaneGeometry plane_geometry(double distance, double source_z, double microphone_z,
        const GroundLine& line)
{
    // Division by the norm, not multiplication by its reciprocal: a level line
    // divides by exactly 1, so the flat case reproduces the plain geometry.
    const double norm = std::sqrt(1.0 + line.slope * line.slope);
    PlaneGeometry g;
    g.source_height = (source_z - (line.intercept + line.slope * distance)) / norm;
    g.receiver_height = (microphone_z - line.intercept) / norm;
    g.ground_distance = (distance + line.slope * (source_z - microphone_z)) / norm;
    return g;
}

} // namespace niceops
