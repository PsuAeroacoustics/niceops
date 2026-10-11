#pragma once

// Contour polygons of a gridded level field (--contour_file): the region at or
// above a level, as polygons with holes, by marching squares.  Plotting draws
// contours too, but only as figures, only in a build with plotting, and not as
// anything a GIS can read or measure.

#include <array>
#include <string>
#include <vector>

namespace niceops {

class LocalFrame;

// One closed ring, its first point not repeated at the end.
using Ring = std::vector<std::array<double, 2>>;

struct ContourPolygon
{
    Ring outer;              // counterclockwise
    std::vector<Ring> holes; // clockwise
};

// The region where values >= level, on the grid x (nx, ascending) by y (ny,
// ascending), values row-major with x fastest (values[iy * nx + ix]), as -o
// writes it.  A non-finite value is below every level.  The grid's edge closes
// the region: a polygon reaching it runs along it.  A cell whose four corners
// straddle the level two and two (a saddle) is decided by its center, the mean
// of the corners.
std::vector<ContourPolygon> contour_polygons(const std::vector<double>& x, const std::vector<double>& y,
        const std::vector<double>& values, double level);

// Signed area of a ring (positive counterclockwise), and a polygon's area, its
// outer ring less its holes.
double ring_area(const Ring& ring);
double polygon_area(const ContourPolygon& polygon);

// Writes the contours of `values` at each of `levels` as a GeoJSON
// FeatureCollection (RFC 7946: WGS84 longitude and latitude, exterior rings
// counterclockwise, holes clockwise), one MultiPolygon feature per level with
// its level, the metric's name and its enclosed area.  Returns the areas, m^2,
// in the order of `levels`.  Throws (input_error) if the file cannot be written.
std::vector<double> write_contour_geojson(const std::string& path, const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& values, const std::vector<double>& levels,
        const LocalFrame& frame, const std::string& metric);

} // namespace niceops
