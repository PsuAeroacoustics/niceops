#include "raster.hpp"

#include <netcdf.h>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "geodetic.hpp"
#include "util/logger.hpp"
#include "util/netcdf_write.hpp"
#include "util/paths.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

constexpr const char* option = "--raster_file";

void nc_check(int status, const std::string& what, const std::string& path)
{
    netcdf::check_write(status, option, what, path);
}

void put_text(int ncid, int var, const char* name, const std::string& value, const std::string& path)
{
    netcdf::put_text(ncid, var, name, value, option, path);
}

void put_double(int ncid, int var, const char* name, double value, const std::string& path)
{
    nc_check(nc_put_att_double(ncid, var, name, NC_DOUBLE, 1, &value),
            std::string("writing attribute ") + name, path);
}

} // namespace

namespace niceops {

void write_raster(const std::string& path, const std::vector<double>& x, const std::vector<double>& y,
        const std::vector<RasterField>& fields, const LocalFrame* frame, double frame_bearing)
{
    const std::size_t nx = x.size();
    const std::size_t ny = y.size();
    ensure_parent_directory(path);
    netcdf::File file;
    nc_check(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file.id), "creating it", path);
    const int ncid = file.id;

    int x_dim = -1, y_dim = -1;
    nc_check(nc_def_dim(ncid, "x", nx, &x_dim), "defining x", path);
    nc_check(nc_def_dim(ncid, "y", ny, &y_dim), "defining y", path);
    const int grid[2] = {y_dim, x_dim};

    int x_var = -1, y_var = -1;
    nc_check(nc_def_var(ncid, "x", NC_DOUBLE, 1, &x_dim, &x_var), "defining x", path);
    nc_check(nc_def_var(ncid, "y", NC_DOUBLE, 1, &y_dim, &y_var), "defining y", path);
    put_text(ncid, x_var, "units", "m", path);
    put_text(ncid, y_var, "units", "m", path);
    put_text(ncid, x_var, "long_name", "x in the run's local frame", path);
    put_text(ncid, y_var, "long_name", "y in the run's local frame", path);
    put_text(ncid, x_var, "axis", "X", path);
    put_text(ncid, y_var, "axis", "Y", path);

    int lat_var = -1, lon_var = -1;
    if (frame)
    {
        nc_check(nc_def_var(ncid, "latitude", NC_DOUBLE, 2, grid, &lat_var), "defining latitude", path);
        nc_check(nc_def_var(ncid, "longitude", NC_DOUBLE, 2, grid, &lon_var), "defining longitude", path);
        put_text(ncid, lat_var, "standard_name", "latitude", path);
        put_text(ncid, lat_var, "units", "degrees_north", path);
        put_text(ncid, lon_var, "standard_name", "longitude", path);
        put_text(ncid, lon_var, "units", "degrees_east", path);
    }

    const double fill = std::numeric_limits<double>::quiet_NaN();
    std::vector<int> field_vars;
    for (const auto& field : fields)
    {
        assert(field.values && field.values->size() == nx * ny);
        int var = -1;
        nc_check(nc_def_var(ncid, field.name.c_str(), NC_DOUBLE, 2, grid, &var), "defining " + field.name, path);
        nc_check(nc_def_var_fill(ncid, var, 0, &fill), "setting the fill value of " + field.name, path);
        put_text(ncid, var, "long_name", field.long_name, path);
        put_text(ncid, var, "units", field.units, path);
        if (frame)
            put_text(ncid, var, "coordinates", "latitude longitude", path);
        field_vars.push_back(var);
    }

    put_text(ncid, NC_GLOBAL, "Conventions", "CF-1.8", path);
    put_text(ncid, NC_GLOBAL, "source", "NICE-OPS", path);
    if (frame)
    {
        put_double(ncid, NC_GLOBAL, "frame_reference_latitude_deg", frame->reference().latitude, path);
        put_double(ncid, NC_GLOBAL, "frame_reference_longitude_deg", frame->reference().longitude, path);
        put_double(ncid, NC_GLOBAL, "frame_reference_altitude_m", frame->reference().altitude, path);
        // NaN for a plain east/north frame: x east and y north.
        put_double(ncid, NC_GLOBAL, "frame_bearing_deg", frame_bearing, path);
    }
    nc_check(nc_enddef(ncid), "ending its definition", path);

    nc_check(nc_put_var_double(ncid, x_var, x.data()), "writing x", path);
    nc_check(nc_put_var_double(ncid, y_var, y.data()), "writing y", path);
    if (frame)
    {
        // Observers sit on the ground plane, as write_geodetic_levels places them.
        std::vector<double> latitude(nx * ny), longitude(nx * ny);
        for (std::size_t iy = 0; iy < ny; iy++)
            for (std::size_t ix = 0; ix < nx; ix++)
            {
                const auto position = frame->to_geodetic(Vec3(x[ix], y[iy], 0.0));
                latitude[iy * nx + ix] = position.latitude;
                longitude[iy * nx + ix] = position.longitude;
            }
        nc_check(nc_put_var_double(ncid, lat_var, latitude.data()), "writing latitude", path);
        nc_check(nc_put_var_double(ncid, lon_var, longitude.data()), "writing longitude", path);
    }
    std::vector<double> values(nx * ny);
    for (std::size_t k = 0; k < fields.size(); k++)
    {
        const auto& source = *fields[k].values;
        for (std::size_t i = 0; i < values.size(); i++)
            values[i] = std::isfinite(source[i]) ? source[i] : fill;
        nc_check(nc_put_var_double(ncid, field_vars[k], values.data()), "writing " + fields[k].name, path);
    }
    nc_check(file.close(), "closing it", path);
}

} // namespace niceops
