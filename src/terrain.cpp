#include "terrain.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <netcdf.h>

#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

// Closes the file however the reader leaves.  input_error() aborts, so the
// guard matters on the error paths as well as the normal one.
struct CloseOnExit
{
    int id;
    ~CloseOnExit()
    {
        nc_close(id);
    }
};

void check(int status, const std::string& what)
{
    if (status != NC_NOERR)
        input_error("Terrain file: ", what, ": ", nc_strerror(status));
}

// Length of the dimension a coordinate variable runs along.
std::size_t axis_length(int ncid, const char* name, std::vector<double>& values)
{
    int varid = 0;
    check(nc_inq_varid(ncid, name, &varid), std::string("no coordinate variable ") + name);
    int dimid = 0;
    std::size_t length = 0;
    check(nc_inq_vardimid(ncid, varid, &dimid), std::string("reading the dimension of ") + name);
    check(nc_inq_dimlen(ncid, dimid, &length), std::string("reading the length of ") + name);
    if (length < 2)
        input_error("Terrain file: ", name, " needs at least two nodes.");
    values.resize(length);
    check(nc_get_var_double(ncid, varid, values.data()), std::string("reading ") + name);
    for (std::size_t i = 1; i < length; i++)
    {
        if (!(values[i] > values[i - 1]))
            input_error("Terrain file: ", name, " must increase strictly.");
    }
    return length;
}

// The cell along a sorted axis that holds v, and how far across it v lies.
// Positions outside the axis are the caller's to exclude before this is used.
std::pair<std::size_t, double> cell_on(const std::vector<double>& axis, double v)
{
    std::size_t hi = static_cast<std::size_t>(std::upper_bound(axis.begin(), axis.end(), v) - axis.begin());
    hi = std::clamp<std::size_t>(hi, 1, axis.size() - 1);
    const std::size_t lo = hi - 1;
    return {lo, (v - axis[lo]) / (axis[hi] - axis[lo])};
}

} // namespace

Terrain::Terrain(const std::string& path)
{
    int ncid = 0;
    check(nc_open(path.c_str(), NC_NOWRITE, &ncid), "opening " + path);
    const CloseOnExit closer {ncid};

    const std::size_t n_latitude = axis_length(ncid, "latitude", latitudes_);
    const std::size_t n_longitude = axis_length(ncid, "longitude", longitudes_);

    int varid = 0;
    check(nc_inq_varid(ncid, "elevation", &varid), "no elevation variable in " + path);
    int dimids[2] = {0, 0};
    check(nc_inq_vardimid(ncid, varid, dimids), "reading the dimensions of elevation");
    std::size_t lengths[2] = {0, 0};
    check(nc_inq_dimlen(ncid, dimids[0], &lengths[0]), "reading the length of elevation");
    check(nc_inq_dimlen(ncid, dimids[1], &lengths[1]), "reading the length of elevation");
    if (lengths[0] != n_latitude || lengths[1] != n_longitude)
        input_error("Terrain file: elevation must run along latitude then longitude, and match them.");
    elevations_.resize(n_latitude * n_longitude);
    check(nc_get_var_double(ncid, varid, elevations_.data()), "reading elevation");

    // Optional: say which vertical datum the heights are on.
    std::size_t datum_length = 0;
    if (NC_NOERR == nc_inq_attlen(ncid, NC_GLOBAL, "vertical_datum", &datum_length) && datum_length > 0)
    {
        std::vector<char> text(datum_length + 1, '\0');
        check(nc_get_att_text(ncid, NC_GLOBAL, "vertical_datum", text.data()), "reading vertical_datum");
        vertical_datum_.assign(text.data(), datum_length);
    }
    else
    {
        vertical_datum_ = "unstated";
    }
}

bool Terrain::covers(double latitude, double longitude) const
{
    return latitude >= latitudes_.front() && latitude <= latitudes_.back()
            && longitude >= longitudes_.front() && longitude <= longitudes_.back();
}

double Terrain::elevation(double latitude, double longitude) const
{
    if (!covers(latitude, longitude))
        return std::numeric_limits<double>::quiet_NaN();
    const auto [i, fi] = cell_on(latitudes_, latitude);
    const auto [j, fj] = cell_on(longitudes_, longitude);
    const std::size_t width = longitudes_.size();
    const double z00 = elevations_[i * width + j];
    const double z01 = elevations_[i * width + j + 1];
    const double z10 = elevations_[(i + 1) * width + j];
    const double z11 = elevations_[(i + 1) * width + j + 1];
    // A missing node makes the whole cell missing, not just the part it weighs.
    return (1.0 - fi) * (1.0 - fj) * z00 + (1.0 - fi) * fj * z01 + fi * (1.0 - fj) * z10
            + fi * fj * z11;
}

GroundPlane niceops::fit_ground_plane(const std::vector<double>& distances, const std::vector<double>& heights)
{
    // Least squares on centered sums, which stays exact for a path many
    // kilometers long where the raw sums of squares would cancel.
    GroundPlane plane;
    const std::size_t n = std::min(distances.size(), heights.size());
    if (n == 0)
        return plane;
    double mean_s = 0.0;
    double mean_z = 0.0;
    for (std::size_t k = 0; k < n; k++)
    {
        mean_s += distances[k];
        mean_z += heights[k];
    }
    mean_s /= static_cast<double>(n);
    mean_z /= static_cast<double>(n);
    double sxx = 0.0;
    double sxz = 0.0;
    for (std::size_t k = 0; k < n; k++)
    {
        const double ds = distances[k] - mean_s;
        sxx += ds * ds;
        sxz += ds * (heights[k] - mean_z);
    }
    plane.slope = sxx > 0.0 ? sxz / sxx : 0.0;
    plane.intercept = mean_z - plane.slope * mean_s;
    double squares = 0.0;
    for (std::size_t k = 0; k < n; k++)
    {
        const double residual = heights[k] - (plane.intercept + plane.slope * distances[k]);
        squares += residual * residual;
    }
    plane.rms = std::sqrt(squares / static_cast<double>(n));
    return plane;
}

std::vector<ProfileSample> niceops::terrain_profile(const Terrain& terrain, const LocalFrame& frame,
        const Vec3& from, const Vec3& to, double spacing)
{
    if (!(spacing > 0.0) || !std::isfinite(spacing))
        input_error("A terrain profile needs a positive spacing, not ", spacing, " m.");
    const Vec3 path = to - from;
    const double length = path.norm();
    const std::size_t steps = std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(length / spacing)));
    std::vector<ProfileSample> profile;
    profile.reserve(steps + 1);
    for (std::size_t k = 0; k <= steps; k++)
    {
        const double t = static_cast<double>(k) / static_cast<double>(steps);
        const Geodetic position = frame.to_geodetic(from + t * path);
        profile.push_back({t * length, terrain.elevation(position.latitude, position.longitude)});
    }
    return profile;
}
