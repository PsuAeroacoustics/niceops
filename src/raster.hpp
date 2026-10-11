#pragma once

// The gridded levels of a run as one georeferenced raster (--raster_file): a
// CF-1.8 netCDF that a GIS (GDAL, QGIS) or xarray opens as it is, where the
// text grids -o writes need parsing and carry no position on the Earth.

#include <string>
#include <vector>

namespace niceops {

class LocalFrame;

struct RasterField
{
    std::string name;      // the variable's name
    std::string long_name; // what it is, for the long_name attribute
    std::string units;     // "dB" or "s", ...
    const std::vector<double>* values = nullptr; // row-major, x fastest, as -o writes it
};

// Writes the fields on the grid x (nx) by y (ny), local meters, each a
// variable (y, x) with the local axes as coordinates.  With a frame, the
// grid's WGS84 latitude and longitude go beside them as two-dimensional
// auxiliary coordinates, and the frame's reference and bearing as global
// attributes.  A non-finite value (an observer that heard nothing) is written
// as the fill value, NaN.  Throws (input_error) if the file cannot be written.
void write_raster(const std::string& path, const std::vector<double>& x, const std::vector<double>& y,
        const std::vector<RasterField>& fields, const LocalFrame* frame, double frame_bearing);

} // namespace niceops
