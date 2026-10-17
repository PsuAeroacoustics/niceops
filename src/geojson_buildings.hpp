#pragma once

// LoD1 buildings from GeoJSON (RFC 7946), for niceops_urban_table --buildings:
// each Polygon or MultiPolygon feature a footprint in WGS84 longitude and
// latitude, projected into the run's local frame by the same transform the run
// gives a geodetic track (geodetic.hpp), its holes kept as courtyards.  The
// height comes from the feature's properties, as OpenStreetMap tags them:
//
//   height, building:height   metres, a number or text ("12", "12.5 m",
//                             "40 ft"; ft, feet and ' mean feet)
//   building:levels           storeys, times `level_height`
//
// in that order of preference, or `default_height` where none gives one.  A
// feature raised off the ground (min_height or building:min_level above
// zero: a bridge, an overhang, a building part on another) is no extruded
// footprint, and is left out and counted.  Optional numeric properties
// `absorption` and `scattering` set the walls' coefficients.

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "buildings.hpp"
#include "geodetic.hpp"

namespace niceops {

struct GeojsonOptions
{
    double level_height = 3.0; // m per storey
    // For a feature that gives no height; NaN leaves it out, and counts it.
    double default_height = std::numeric_limits<double>::quiet_NaN();
};

struct GeojsonBuilding
{
    Footprint footprint;
    std::string id;      // the feature's id, or its position in the file
    double absorption = std::numeric_limits<double>::quiet_NaN();
    double scattering = std::numeric_limits<double>::quiet_NaN();
};

struct GeojsonReport
{
    std::size_t features = 0;
    std::size_t footprints = 0;
    std::size_t defaulted = 0;      // took default_height
    std::size_t no_height = 0;      // left out: no height and no default
    std::size_t raised = 0;         // left out: off the ground
    std::size_t not_polygons = 0;   // left out: points, lines, no geometry
};

// Throws std::runtime_error for a file that is not GeoJSON, a coordinate
// reference system other than WGS84 longitude and latitude, or a malformed
// polygon.
std::vector<GeojsonBuilding> read_geojson_buildings(const std::string& path, const LocalFrame& frame,
        const GeojsonOptions& options, GeojsonReport& report);

// A height property's value in metres: a number, or text with an optional
// unit; NaN for anything else.
double geojson_height(const std::string& text);

} // namespace niceops
