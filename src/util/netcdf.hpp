#pragma once

// Shape-checked netCDF reads.  nc_get_var_* and nc_get_att_* write every value
// the file holds, whatever the buffer they are handed, so a buffer sized from
// anything but the variable's (or attribute's) own shape can be overrun by a
// malformed file before any later check could refuse it.  These size it from
// the shape; the caller then checks the shape is what it expects.

#include <netcdf.h>

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace niceops::netcdf {

// A variable's dimension ids, outermost first; nullopt if they cannot be read.
inline std::optional<std::vector<int>> dimension_ids(int ncid, int varid)
{
    int ndims = 0;
    if (nc_inq_varndims(ncid, varid, &ndims) != NC_NOERR || ndims < 0 || ndims > NC_MAX_VAR_DIMS)
        return std::nullopt;
    std::vector<int> dims(static_cast<std::size_t>(ndims));
    if (ndims > 0 && nc_inq_vardimid(ncid, varid, dims.data()) != NC_NOERR)
        return std::nullopt;
    return dims;
}

// The names of a variable's dimensions, outermost first; nullopt if they
// cannot be read.
inline std::optional<std::vector<std::string>> dimension_names(int ncid, int varid)
{
    const auto dims = dimension_ids(ncid, varid);
    if (!dims)
        return std::nullopt;
    std::vector<std::string> names;
    for (const int dim : *dims)
    {
        char name[NC_MAX_NAME + 1] = {};
        if (nc_inq_dimname(ncid, dim, name) != NC_NOERR)
            return std::nullopt;
        names.emplace_back(name);
    }
    return names;
}

// How many values a variable holds; nullopt if its shape cannot be read or the
// count does not fit in a size_t.
inline std::optional<std::size_t> element_count(int ncid, int varid)
{
    const auto dims = dimension_ids(ncid, varid);
    if (!dims)
        return std::nullopt;
    std::size_t count = 1;
    for (const int dim : *dims)
    {
        std::size_t length = 0;
        if (nc_inq_dimlen(ncid, dim, &length) != NC_NOERR
                || (length != 0 && count > std::numeric_limits<std::size_t>::max() / length))
            return std::nullopt;
        count *= length;
    }
    return count;
}

// Every value of a variable, as double, in its own (row-major) order; nullopt
// if it cannot be read.
inline std::optional<std::vector<double>> read_doubles(int ncid, int varid)
{
    const auto count = element_count(ncid, varid);
    if (!count)
        return std::nullopt;
    std::vector<double> values(*count);
    if (*count > 0 && nc_get_var_double(ncid, varid, values.data()) != NC_NOERR)
        return std::nullopt;
    return values;
}

enum class Single
{
    ok,
    missing,
    not_one_number, // more or fewer than one value, or not numeric
};

// An attribute that holds one number.
inline Single single_value(int ncid, int varid, const char* name, double& value)
{
    std::size_t length = 0;
    if (nc_inq_attlen(ncid, varid, name, &length) != NC_NOERR)
        return Single::missing;
    if (length != 1 || nc_get_att_double(ncid, varid, name, &value) != NC_NOERR)
        return Single::not_one_number;
    return Single::ok;
}

inline Single single_value(int ncid, int varid, const char* name, int& value)
{
    std::size_t length = 0;
    if (nc_inq_attlen(ncid, varid, name, &length) != NC_NOERR)
        return Single::missing;
    if (length != 1 || nc_get_att_int(ncid, varid, name, &value) != NC_NOERR)
        return Single::not_one_number;
    return Single::ok;
}

} // namespace niceops::netcdf
