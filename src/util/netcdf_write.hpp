#pragma once

// The netCDF writers' shared pieces (--band_frames, --raster_file,
// --provenance_spheres): a failed call refuses with the option and file named,
// and a file handle that closes on scope exit.

#include <netcdf.h>

#include <string>

#include "util/logger.hpp"

namespace niceops::netcdf {

inline void check_write(int status, const char* option, const std::string& what, const std::string& path)
{
    if (status != NC_NOERR)
        input_error(option, " ", path, ": netCDF error while ", what, ": ", nc_strerror(status), ".");
}

// Closes on scope exit unless closed explicitly, which is the only close whose
// status is looked at.
struct File
{
    int id = -1;
    File() = default;
    File(const File&) = delete;
    File& operator=(const File&) = delete;
    ~File()
    {
        if (id >= 0)
            nc_close(id);
    }
    int close()
    {
        const int status = id >= 0 ? nc_close(id) : NC_NOERR;
        id = -1;
        return status;
    }
};

inline void put_text(int ncid, int var, const char* name, const std::string& value, const char* option,
        const std::string& path)
{
    check_write(nc_put_att_text(ncid, var, name, value.size(), value.c_str()),
            option, std::string("writing attribute ") + name, path);
}

} // namespace niceops::netcdf
