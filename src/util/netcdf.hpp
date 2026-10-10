#pragma once

// Shape checks for netCDF reads.  nc_get_var_* and nc_get_att_* write every
// value the file holds, whatever the buffer they are handed, so a reader that
// sizes its buffer from a dimension (or from one value) must first check that
// the variable or attribute is that shape.  A malformed file would otherwise
// overrun the buffer before any later check could refuse it.

#include <netcdf.h>

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace niceops::netcdf {

// The names of a variable's dimensions, outermost first; nullopt if they
// cannot be read.
inline std::optional<std::vector<std::string>> dimension_names(int ncid, int varid)
{
    int ndims = 0;
    if (nc_inq_varndims(ncid, varid, &ndims) != NC_NOERR || ndims < 0 || ndims > NC_MAX_VAR_DIMS)
        return std::nullopt;
    std::vector<int> dims(static_cast<std::size_t>(ndims));
    if (ndims > 0 && nc_inq_vardimid(ncid, varid, dims.data()) != NC_NOERR)
        return std::nullopt;
    std::vector<std::string> names;
    for (const int dim : dims)
    {
        char name[NC_MAX_NAME + 1] = {};
        if (nc_inq_dimname(ncid, dim, name) != NC_NOERR)
            return std::nullopt;
        names.emplace_back(name);
    }
    return names;
}

// a * b, or nullopt where it would not fit in a size_t.
inline std::optional<std::size_t> checked_product(std::size_t a, std::size_t b)
{
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a)
        return std::nullopt;
    return a * b;
}

// How many values a variable holds; nullopt if its shape cannot be read or the
// count does not fit in a size_t.
inline std::optional<std::size_t> element_count(int ncid, int varid)
{
    int ndims = 0;
    if (nc_inq_varndims(ncid, varid, &ndims) != NC_NOERR || ndims < 0 || ndims > NC_MAX_VAR_DIMS)
        return std::nullopt;
    std::vector<int> dims(static_cast<std::size_t>(ndims));
    if (ndims > 0 && nc_inq_vardimid(ncid, varid, dims.data()) != NC_NOERR)
        return std::nullopt;
    std::optional<std::size_t> count = 1;
    for (const int dim : dims)
    {
        std::size_t length = 0;
        if (nc_inq_dimlen(ncid, dim, &length) != NC_NOERR)
            return std::nullopt;
        count = checked_product(*count, length);
        if (!count)
            return std::nullopt;
    }
    return count;
}

// How many values an attribute holds; nullopt if there is no such attribute.
inline std::optional<std::size_t> attribute_length(int ncid, int varid, const char* name)
{
    std::size_t length = 0;
    if (nc_inq_attlen(ncid, varid, name, &length) != NC_NOERR)
        return std::nullopt;
    return length;
}

} // namespace niceops::netcdf
