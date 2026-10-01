#include "propagation_table.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <netcdf.h>

#include "util/logger.hpp"
#include "util/paths.hpp"

namespace niceops {

namespace {

constexpr const char* axis_names[5] = {"bearing", "source_height", "receiver_height", "range", "band"};

void check_axes(const PropagationTable::Axes& axes, std::size_t values)
{
    auto ascending = [](const std::vector<double>& v) {
        for (std::size_t k = 1; k < v.size(); k++)
            if (!(v[k] > v[k - 1]))
                return false;
        return !v.empty() && std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x); });
    };
    if (!ascending(axes.source_heights) || !ascending(axes.receiver_heights) || !ascending(axes.ranges)
            || !ascending(axes.bands))
        throw std::invalid_argument("propagation table: the height, range and band axes must be finite "
                                    "and ascend");
    const std::size_t n = axes.bearings.size();
    if (n == 0)
        throw std::invalid_argument("propagation table: no bearings");
    for (std::size_t k = 0; k < n; k++)
        if (std::abs(axes.bearings[k] - 360.0 * static_cast<double>(k) / static_cast<double>(n)) > 1e-9)
            throw std::invalid_argument("propagation table: bearings must be uniform from 0 over 360 degrees");
    const std::size_t expected = n * axes.source_heights.size() * axes.receiver_heights.size()
            * axes.ranges.size() * axes.bands.size();
    if (values != expected)
        throw std::invalid_argument("propagation table: the values do not fill the axes");
}

// The cell of an ascending axis that holds x, and the weight of its upper
// end; false outside the axis.  A one-point axis holds only its point.
bool locate(const std::vector<double>& axis, double x, std::size_t& lower, double& fraction)
{
    if (axis.size() == 1)
    {
        lower = 0;
        fraction = 0.0;
        return std::abs(x - axis[0]) <= 1e-6 * std::max(1.0, std::abs(axis[0]));
    }
    if (!(x >= axis.front()) || !(x <= axis.back()))
        return false;
    const auto upper = std::upper_bound(axis.begin(), axis.end(), x);
    lower = std::min(static_cast<std::size_t>(upper - axis.begin()), axis.size() - 1) - 1;
    fraction = (x - axis[lower]) / (axis[lower + 1] - axis[lower]);
    return true;
}

void nc_check(int status, const std::string& what, const std::string& path)
{
    if (status != NC_NOERR)
        input_error("Propagation table ", path, ": ", what, " (", nc_strerror(status), ").");
}

} // namespace

PropagationTable::PropagationTable(Axes axes, std::vector<double> values, std::string method,
        bool includes_absorption, std::map<std::string, std::string> provenance)
    : axes_(std::move(axes))
    , values_(std::move(values))
    , method_(std::move(method))
    , includes_absorption_(includes_absorption)
    , provenance_(std::move(provenance))
{
    check_axes(axes_, values_.size());
}

void PropagationTable::save(const std::string& path) const
{
    ensure_parent_directory(path);
    int ncid = 0;
    nc_check(nc_create(path.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid), "cannot create it", path);
    const std::vector<double>* axes[5] = {&axes_.bearings, &axes_.source_heights, &axes_.receiver_heights,
            &axes_.ranges, &axes_.bands};
    const char* units[5] = {"degrees", "m", "m", "m", "Hz"};
    int dims[5] = {};
    int coordinates[5] = {};
    for (int k = 0; k < 5; k++)
    {
        nc_check(nc_def_dim(ncid, axis_names[k], axes[k]->size(), &dims[k]), "defining a dimension", path);
        nc_check(nc_def_var(ncid, axis_names[k], NC_DOUBLE, 1, &dims[k], &coordinates[k]), "defining an axis",
                path);
        nc_put_att_text(ncid, coordinates[k], "units", std::string(units[k]).size(), units[k]);
    }
    int level = 0;
    nc_check(nc_def_var(ncid, "relative_level", NC_DOUBLE, 5, dims, &level), "defining the values", path);
    const std::string level_units = "dB re free field at the straight-line distance; NaN where there is no value";
    nc_put_att_text(ncid, level, "units", level_units.size(), level_units.c_str());
    const int version = format_version;
    nc_put_att_int(ncid, NC_GLOBAL, "propagation_table_version", NC_INT, 1, &version);
    const int absorption = includes_absorption_ ? 1 : 0;
    nc_put_att_int(ncid, NC_GLOBAL, "includes_absorption", NC_INT, 1, &absorption);
    nc_put_att_text(ncid, NC_GLOBAL, "method", method_.size(), method_.c_str());
    const std::string bearing_note = "compass bearing the sound travels toward, source to observer";
    nc_put_att_text(ncid, NC_GLOBAL, "bearing_convention", bearing_note.size(), bearing_note.c_str());
    for (const auto& [key, value] : provenance_)
        nc_put_att_text(ncid, NC_GLOBAL, ("provenance_" + key).c_str(), value.size(), value.c_str());
    nc_check(nc_enddef(ncid), "ending its definition", path);
    for (int k = 0; k < 5; k++)
        nc_check(nc_put_var_double(ncid, coordinates[k], axes[k]->data()), "writing an axis", path);
    nc_check(nc_put_var_double(ncid, level, values_.data()), "writing the values", path);
    nc_check(nc_close(ncid), "closing it", path);
}

PropagationTable PropagationTable::load(const std::string& path)
{
    int ncid = 0;
    if (nc_open(path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR)
        input_error("Cannot open propagation table ", path, ".");
    struct Closer
    {
        int id;
        ~Closer()
        {
            nc_close(id);
        }
    } closer {ncid};

    int version = 0;
    if (nc_get_att_int(ncid, NC_GLOBAL, "propagation_table_version", &version) != NC_NOERR)
        input_error("Propagation table ", path, " has no propagation_table_version attribute: it is not "
                                                "one, or it predates the format.");
    if (version != format_version)
        input_error("Propagation table ", path, " is format version ", version, "; this build reads ",
                format_version, ".");

    Axes axes;
    std::vector<double>* targets[5] = {&axes.bearings, &axes.source_heights, &axes.receiver_heights,
            &axes.ranges, &axes.bands};
    std::size_t count = 1;
    for (int k = 0; k < 5; k++)
    {
        int dim = 0, var = 0;
        std::size_t length = 0;
        if (nc_inq_dimid(ncid, axis_names[k], &dim) != NC_NOERR || nc_inq_dimlen(ncid, dim, &length) != NC_NOERR
                || nc_inq_varid(ncid, axis_names[k], &var) != NC_NOERR)
            input_error("Propagation table ", path, " has no '", axis_names[k], "' axis.");
        targets[k]->resize(length);
        nc_check(nc_get_var_double(ncid, var, targets[k]->data()), std::string("reading ") + axis_names[k], path);
        count *= length;
    }
    int level = 0;
    if (nc_inq_varid(ncid, "relative_level", &level) != NC_NOERR)
        input_error("Propagation table ", path, " has no 'relative_level' variable.");
    int ndims = 0;
    int dims[NC_MAX_VAR_DIMS] = {};
    nc_inq_varndims(ncid, level, &ndims);
    nc_inq_vardimid(ncid, level, dims);
    for (int k = 0; k < 5 && ndims == 5; k++)
    {
        char name[NC_MAX_NAME + 1] = {};
        nc_inq_dimname(ncid, dims[k], name);
        if (std::string(name) != axis_names[k])
            ndims = -1;
    }
    if (ndims != 5)
        input_error("Propagation table ", path, ": relative_level must be dimensioned (bearing, source_height, "
                                                "receiver_height, range, band).");
    std::vector<double> values(count);
    nc_check(nc_get_var_double(ncid, level, values.data()), "reading relative_level", path);

    auto text = [&](const std::string& name) {
        std::size_t length = 0;
        if (nc_inq_attlen(ncid, NC_GLOBAL, name.c_str(), &length) != NC_NOERR)
            return std::string();
        std::string value(length, '\0');
        nc_get_att_text(ncid, NC_GLOBAL, name.c_str(), value.data());
        return value;
    };
    int absorption = 0;
    nc_get_att_int(ncid, NC_GLOBAL, "includes_absorption", &absorption);
    std::map<std::string, std::string> provenance;
    int attributes = 0;
    nc_inq_natts(ncid, &attributes);
    for (int k = 0; k < attributes; k++)
    {
        char name[NC_MAX_NAME + 1] = {};
        nc_inq_attname(ncid, NC_GLOBAL, k, name);
        const std::string key(name);
        if (key.rfind("provenance_", 0) == 0)
            provenance[key.substr(11)] = text(key);
    }
    try
    {
        return PropagationTable(std::move(axes), std::move(values), text("method"), absorption != 0,
                std::move(provenance));
    }
    catch (const std::invalid_argument& e)
    {
        input_error("Propagation table ", path, ": ", e.what(), ".");
    }
}

PropagationTable::BandMap PropagationTable::match_bands(const std::vector<double>& run_bands,
        std::size_t first_band, std::size_t last_band) const
{
    last_band = std::min(last_band, run_bands.size());
    BandMap map;
    map.first = first_band;
    map.last = last_band;
    map.column.assign(run_bands.size(), 0);
    std::string missing;
    for (std::size_t band = first_band; band < last_band; band++)
    {
        const auto& table = axes_.bands;
        const auto nearest = std::min_element(table.begin(), table.end(), [&](double a, double b) {
            return std::abs(std::log(a / run_bands[band])) < std::abs(std::log(b / run_bands[band]));
        });
        if (std::abs(*nearest / run_bands[band] - 1.0) > 0.005)
            missing += (missing.empty() ? "" : ", ") + std::to_string(run_bands[band]);
        map.column[band] = static_cast<std::size_t>(nearest - table.begin());
    }
    if (!missing.empty())
        throw std::runtime_error("the propagation table has no band at " + missing
                + " Hz, which the run's band limits include");
    return map;
}

bool PropagationTable::covers_receiver(double height) const
{
    std::size_t lower = 0;
    double fraction = 0.0;
    return locate(axes_.receiver_heights, height, lower, fraction);
}

bool PropagationTable::lookup(double bearing_degrees, double source_height, double receiver_height, double range,
        const BandMap& bands, double* out) const
{
    std::size_t s0 = 0, h0 = 0, r0 = 0;
    double fs = 0.0, fh = 0.0, fr = 0.0;
    if (!locate(axes_.source_heights, source_height, s0, fs) || !locate(axes_.receiver_heights, receiver_height, h0, fh)
            || !locate(axes_.ranges, range, r0, fr))
        return false;
    // Bearing: periodic, uniform.
    const std::size_t nb = axes_.bearings.size();
    std::size_t b0 = 0, b1 = 0;
    double fb = 0.0;
    if (nb > 1)
    {
        double u = std::fmod(bearing_degrees, 360.0);
        if (u < 0.0)
            u += 360.0;
        u *= static_cast<double>(nb) / 360.0;
        b0 = std::min(static_cast<std::size_t>(u), nb - 1);
        fb = u - static_cast<double>(b0);
        b1 = (b0 + 1) % nb;
    }
    // The sixteen corners and their weights; a corner of zero weight is not
    // read, so a NaN just past a point that lies on the grid does not spoil it.
    struct Corner
    {
        std::size_t offset;
        double weight;
    };
    Corner corners[16];
    int used = 0;
    const std::size_t bearing_index[2] = {b0, b1};
    const double bearing_weight[2] = {1.0 - fb, fb};
    const double source_weight[2] = {1.0 - fs, fs};
    const double receiver_weight[2] = {1.0 - fh, fh};
    const double range_weight[2] = {1.0 - fr, fr};
    const std::size_t nbands = axes_.bands.size();
    for (int b = 0; b < 2; b++)
        for (int s = 0; s < 2; s++)
            for (int h = 0; h < 2; h++)
                for (int r = 0; r < 2; r++)
                {
                    const double weight = bearing_weight[b] * source_weight[s] * receiver_weight[h] * range_weight[r];
                    if (weight == 0.0)
                        continue;
                    corners[used++] = {index(bearing_index[b], s0 + s, h0 + h, r0 + r) * nbands, weight};
                }
    for (std::size_t band = bands.first; band < bands.last; band++)
    {
        const std::size_t column = bands.column[band];
        double sum = 0.0;
        for (int c = 0; c < used; c++)
            sum += corners[c].weight * values_[corners[c].offset + column];
        if (!std::isfinite(sum))
            return false;
        out[band] = sum;
    }
    return true;
}

std::vector<double> PropagationTable::max_power_ratio(const BandMap& bands) const
{
    std::vector<double> ratio(bands.column.size(), 1.0);
    const std::size_t nbands = axes_.bands.size();
    const std::size_t cells = values_.size() / nbands;
    for (std::size_t band = bands.first; band < bands.last; band++)
    {
        double largest = -std::numeric_limits<double>::infinity();
        for (std::size_t cell = 0; cell < cells; cell++)
        {
            const double v = values_[cell * nbands + bands.column[band]];
            if (std::isfinite(v))
                largest = std::max(largest, v);
        }
        ratio[band] = std::max(1.0, std::pow(10.0, 0.1 * largest));
    }
    return ratio;
}

} // namespace niceops
