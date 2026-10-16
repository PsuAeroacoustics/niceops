#include "terrain.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <netcdf.h>

#include "util/logger.hpp"
#include "util/netcdf.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- matches the rest of src/

namespace {

constexpr double pi = 3.14159265358979323846;

// Closes the file however the reader leaves.  input_error() throws, so the
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

// A coordinate variable's values, along whatever dimension it runs.
void read_axis(int ncid, const char* name, std::vector<double>& values)
{
    int varid = 0;
    check(nc_inq_varid(ncid, name, &varid), std::string("no coordinate variable ") + name);
    int dimension_count = 0;
    check(nc_inq_varndims(ncid, varid, &dimension_count), std::string("reading the rank of ") + name);
    if (dimension_count != 1)
        input_error("Terrain file: ", name, " must be 1-D, not ", dimension_count, "-D.");
    auto read = netcdf::read_doubles(ncid, varid);
    if (!read)
        input_error("Terrain file: cannot read ", name, ".");
    values = std::move(*read);
}

// A numeric attribute of a variable, if it has one.
bool double_attribute(int ncid, int varid, const char* name, double& value)
{
    std::size_t length = 0;
    if (nc_inq_attlen(ncid, varid, name, &length) != NC_NOERR || length == 0)
        return false;
    std::vector<double> values(length);
    check(nc_get_att_double(ncid, varid, name, values.data()), std::string("reading elevation:") + name);
    value = values.front();
    return true;
}

std::string text_attribute(int ncid, int varid, const char* name)
{
    std::size_t length = 0;
    if (nc_inq_attlen(ncid, varid, name, &length) != NC_NOERR || length == 0)
        return {};
    nc_type type = NC_NAT;
    check(nc_inq_atttype(ncid, varid, name, &type), std::string("reading the type of ") + name);
    if (type != NC_CHAR)
        return {};
    std::vector<char> text(length + 1, '\0');
    check(nc_get_att_text(ncid, varid, name, text.data()), std::string("reading ") + name);
    return std::string(text.data(), length);
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

// Meters per degree of latitude and of longitude at `latitude`, on the WGS84
// ellipsoid (the usual series), for stating extents and spacings in the log.
double meters_per_degree_latitude(double latitude)
{
    const double phi = latitude * pi / 180.0;
    return 111132.954 - 559.822 * std::cos(2.0 * phi) + 1.175 * std::cos(4.0 * phi);
}

double meters_per_degree_longitude(double latitude)
{
    const double phi = latitude * pi / 180.0;
    return 111412.84 * std::cos(phi) - 93.5 * std::cos(3.0 * phi);
}

} // namespace

Terrain::Terrain(const std::string& path)
{
    int ncid = 0;
    check(nc_open(path.c_str(), NC_NOWRITE, &ncid), "opening " + path);
    const CloseOnExit closer {ncid};

    read_axis(ncid, "latitude", latitudes_);
    read_axis(ncid, "longitude", longitudes_);

    int varid = 0;
    check(nc_inq_varid(ncid, "elevation", &varid), "no elevation variable in " + path);
    int dimension_count = 0;
    check(nc_inq_varndims(ncid, varid, &dimension_count), "reading the rank of elevation");
    if (dimension_count != 2)
        input_error("Terrain file ", path, ": elevation must be 2-D (latitude, longitude), not ",
                dimension_count, "-D.");
    int dimids[2] = {0, 0};
    check(nc_inq_vardimid(ncid, varid, dimids), "reading the dimensions of elevation");
    std::size_t lengths[2] = {0, 0};
    check(nc_inq_dimlen(ncid, dimids[0], &lengths[0]), "reading the length of elevation");
    check(nc_inq_dimlen(ncid, dimids[1], &lengths[1]), "reading the length of elevation");
    if (lengths[0] != latitudes_.size() || lengths[1] != longitudes_.size())
        input_error("Terrain file ", path, ": elevation must run along latitude then longitude, "
                "and match them.");
    // The heights must be meters: a grid in feet read as meters puts the ground
    // three times too high or low with nothing downstream able to tell.
    const std::string units = text_attribute(ncid, varid, "units");
    if (!units.empty() && units != "m" && units != "meter" && units != "meters" && units != "metre"
            && units != "metres")
    {
        input_error("Terrain file ", path, ": elevation is in '", units, "'; NICEOPS reads "
                "terrain heights in meters. Convert the grid to meters.");
    }
    elevations_.resize(lengths[0] * lengths[1]);
    if (!elevations_.empty())
        check(nc_get_var_double(ncid, varid, elevations_.data()), "reading elevation");
    // Missing nodes become NaN, whatever the file marked them with, so that
    // elevation() reports them rather than reading -9999 m or 1e20 m as ground.
    // A packed grid (CF scale_factor/add_offset) is unpacked after the test,
    // which is on the stored value.
    double fill = 0.0;
    const bool has_fill = double_attribute(ncid, varid, "_FillValue", fill);
    double missing = 0.0;
    const bool has_missing = double_attribute(ncid, varid, "missing_value", missing);
    double scale = 1.0;
    double offset = 0.0;
    // `|`, not `||`: both attributes have to be read.
    const bool packed = double_attribute(ncid, varid, "scale_factor", scale)
            | double_attribute(ncid, varid, "add_offset", offset);
    for (double& z : elevations_)
    {
        if ((has_fill && z == fill) || (has_missing && z == missing))
            z = std::numeric_limits<double>::quiet_NaN();
        else if (packed)
            z = z * scale + offset;
    }

    vertical_datum_ = text_attribute(ncid, NC_GLOBAL, "vertical_datum");
    if (vertical_datum_.empty())
        vertical_datum_ = "unstated";
    try
    {
        validate();
    }
    catch (const std::runtime_error& problem)
    {
        input_error("Terrain file ", path, ": ", problem.what());
    }
}

Terrain::Terrain(std::vector<double> latitudes, std::vector<double> longitudes,
        std::vector<double> elevations, std::string vertical_datum)
    : latitudes_(std::move(latitudes))
    , longitudes_(std::move(longitudes))
    , elevations_(std::move(elevations))
    , vertical_datum_(std::move(vertical_datum))
{
    validate();
}

void Terrain::validate()
{
    for (const auto* axis : {&latitudes_, &longitudes_})
    {
        const char* name = axis == &latitudes_ ? "latitude" : "longitude";
        const double limit = axis == &latitudes_ ? 90.0 : 180.0;
        if (axis->size() < 2)
            throw std::runtime_error(std::string(name) + " needs at least two nodes.");
        for (std::size_t i = 0; i < axis->size(); i++)
        {
            if (!(std::fabs((*axis)[i]) <= limit))
                throw std::runtime_error(std::string(name) + " holds " + std::to_string((*axis)[i])
                        + ", which is not a WGS84 " + name + " in degrees.");
            if (i > 0 && !((*axis)[i] > (*axis)[i - 1]))
                throw std::runtime_error(std::string(name) + " must increase strictly.");
        }
    }
    if (elevations_.size() != latitudes_.size() * longitudes_.size())
        throw std::runtime_error("elevation must hold one height per (latitude, longitude) node.");
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
    // Written as nested linear interpolations from a corner, which returns a
    // constant grid's height exactly: the four weights of the product form do
    // not sum to exactly 1, and level ground has to come out level to the bit.
    const double south = z00 + fj * (z01 - z00);
    const double north = z10 + fj * (z11 - z10);
    return south + fi * (north - south);
}

double Terrain::spacing_meters() const
{
    const double center = 0.5 * (latitudes_.front() + latitudes_.back());
    const double north_south = (latitudes_.back() - latitudes_.front())
            / static_cast<double>(latitudes_.size() - 1) * meters_per_degree_latitude(center);
    const double east_west = (longitudes_.back() - longitudes_.front())
            / static_cast<double>(longitudes_.size() - 1) * meters_per_degree_longitude(center);
    return std::min(north_south, east_west);
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

TerrainRaster::TerrainRaster(const Terrain& terrain, const LocalFrame& frame, double x_min,
        double x_max, double y_min, double y_max, double step, double height_offset)
{
    if (!(step > 0.0) || !std::isfinite(step))
        input_error("--terrain_step must be a positive length, not ", step, " m.");
    if (!(x_max >= x_min) || !(y_max >= y_min) || !std::isfinite(x_min + x_max + y_min + y_max))
        input_error("Terrain raster: no extent to cover (", x_min, " to ", x_max, " by ", y_min,
                " to ", y_max, " m).");
    step_ = step;
    inverse_step_ = 1.0 / step;
    x0_ = x_min;
    y0_ = y_min;
    nx_ = static_cast<std::size_t>(std::ceil((x_max - x_min) * inverse_step_)) + 1;
    ny_ = static_cast<std::size_t>(std::ceil((y_max - y_min) * inverse_step_)) + 1;
    nx_ = std::max<std::size_t>(nx_, 2);
    ny_ = std::max<std::size_t>(ny_, 2);
    // 1e8 nodes is 800 MB; a run that big wants a coarser step, not an
    // out-of-memory kill.
    constexpr double most_nodes = 1e8;
    if (static_cast<double>(nx_) * static_cast<double>(ny_) > most_nodes)
    {
        input_error("--terrain: sampling ", x_max - x_min, " by ", y_max - y_min, " m every ", step,
                " m would take ", static_cast<double>(nx_) * static_cast<double>(ny_),
                " nodes. Give a coarser --terrain_step.");
    }
    heights_.resize(nx_ * ny_);
    std::size_t uncovered = 0;
    std::size_t missing = 0;
    double lat_lo = 90.0, lat_hi = -90.0, lon_lo = 180.0, lon_hi = -180.0;
    double void_lat_lo = 90.0, void_lat_hi = -90.0, void_lon_lo = 180.0, void_lon_hi = -180.0;
    for (std::size_t j = 0; j < ny_; j++)
    {
        for (std::size_t i = 0; i < nx_; i++)
        {
            const Vec3 local(x0_ + static_cast<double>(i) * step_, y0_ + static_cast<double>(j) * step_, 0.0);
            const Geodetic where = frame.to_geodetic(local);
            lat_lo = std::min(lat_lo, where.latitude);
            lat_hi = std::max(lat_hi, where.latitude);
            lon_lo = std::min(lon_lo, where.longitude);
            lon_hi = std::max(lon_hi, where.longitude);
            const double z = terrain.elevation(where.latitude, where.longitude);
            if (!terrain.covers(where.latitude, where.longitude))
                uncovered++;
            else if (!std::isfinite(z))
            {
                missing++;
                void_lat_lo = std::min(void_lat_lo, where.latitude);
                void_lat_hi = std::max(void_lat_hi, where.latitude);
                void_lon_lo = std::min(void_lon_lo, where.longitude);
                void_lon_hi = std::max(void_lon_hi, where.longitude);
            }
            heights_[j * nx_ + i] = z + height_offset;
        }
    }
    if (uncovered > 0)
    {
        input_error("--terrain covers latitudes ", terrain.south(), " to ", terrain.north(),
                " and longitudes ", terrain.west(), " to ", terrain.east(), ", but this run needs ",
                lat_lo, " to ", lat_hi, " and ", lon_lo, " to ", lon_hi, " (every observer and every "
                "state's footprint, with a margin); ", uncovered, " of ", heights_.size(),
                " sample points fall outside it. Give a terrain grid that covers the whole scene "
                "(check --origin and --frame_bearing if it should).");
    }
    if (missing > 0)
    {
        input_error("--terrain has no height at ", missing, " of the ", heights_.size(),
                " points this run samples (missing nodes between latitudes ", void_lat_lo, " and ",
                void_lat_hi, ", longitudes ", void_lon_lo, " and ", void_lon_hi, "). Fill the grid's "
                "voids first.");
    }
}

TerrainRaster::TerrainRaster(double x0, double y0, double step, std::size_t nx, std::size_t ny,
        std::vector<double> heights)
    : x0_(x0)
    , y0_(y0)
    , step_(step)
    , inverse_step_(1.0 / step)
    , nx_(nx)
    , ny_(ny)
    , heights_(std::move(heights))
{
    if (!(step > 0.0) || nx < 2 || ny < 2 || heights_.size() != nx * ny)
        input_error("Terrain raster: needs a positive step, at least 2 x 2 nodes, and a height per node.");
}

double TerrainRaster::height(double x, double y) const
{
    // Clamped into the grid: a point past the edge reads the edge.
    const double u = std::clamp((x - x0_) * inverse_step_, 0.0, static_cast<double>(nx_ - 1));
    const double v = std::clamp((y - y0_) * inverse_step_, 0.0, static_cast<double>(ny_ - 1));
    const std::size_t i = std::min(static_cast<std::size_t>(u), nx_ - 2);
    const std::size_t j = std::min(static_cast<std::size_t>(v), ny_ - 2);
    const double fu = u - static_cast<double>(i);
    const double fv = v - static_cast<double>(j);
    const double* row = heights_.data() + j * nx_ + i;
    const double* next = row + nx_;
    // From a corner, like Terrain::elevation, so level ground stays level.
    const double low = row[0] + fu * (row[1] - row[0]);
    const double high = next[0] + fu * (next[1] - next[0]);
    return low + fv * (high - low);
}

namespace {

// The least-squares line through n + 1 heights equally spaced `step` apart,
// given as their sum and their index-weighted sum (heights relative to
// whatever `base` they were taken from).  Closed form: the indices' mean is
// n/2 and their centered sum of squares n(n+1)(n+2)/12.
GroundLine line_from_sums(double sum_z, double sum_kz, std::size_t n, double step, double base)
{
    const double count = static_cast<double>(n);
    const double sxx = count * (count + 1.0) * (count + 2.0) / 12.0;
    const double per_sample = (sum_kz - 0.5 * count * sum_z) / sxx;
    GroundLine line;
    line.slope = per_sample / step;
    line.intercept = base + sum_z / (count + 1.0) - per_sample * 0.5 * count;
    return line;
}

} // namespace

GroundLine niceops::fit_ground_line(const TerrainRaster& raster, double x0, double y0, double x1,
        double y1, std::size_t n)
{
    n = std::max<std::size_t>(n, 1);
    const double base = raster.height(x0, y0);
    double sum_z = 0.0;
    double sum_kz = 0.0;
    for (std::size_t k = 0; k <= n; k++)
    {
        const double t = static_cast<double>(k) / static_cast<double>(n);
        const double z = raster.height(x0 + t * (x1 - x0), y0 + t * (y1 - y0)) - base;
        sum_z += z;
        sum_kz += static_cast<double>(k) * z;
    }
    const double length = std::hypot(x1 - x0, y1 - y0);
    return line_from_sums(sum_z, sum_kz, n, length / static_cast<double>(n), base);
}

GroundPlaneCache::GroundPlaneCache(const TerrainRaster& raster, bool screening, double min_protrusion)
    : raster_(&raster)
    , step_(raster.step())
    , screening_(screening)
    , min_protrusion_(min_protrusion)
    , cos_(bearings)
    , sin_(bearings)
    , stamp_(bearings, 0u)
    , sum_z_(bearings)
    , sum_kz_(bearings)
    , z_raw_(screening ? bearings : 0)
    , horizon_(screening ? bearings : 0)
{
    for (int ray = 0; ray < bearings; ray++)
    {
        const double angle = 2.0 * pi * ray / bearings;
        cos_[ray] = std::cos(angle);
        sin_[ray] = std::sin(angle);
    }
}

void GroundPlaneCache::set_observer(double x, double y)
{
    x_ = x;
    y_ = y;
    z_ = raster_->height(x, y);
    if (++generation_ == 0)
    {
        // Wrapped: an old stamp could be taken for this observer's.
        std::fill(stamp_.begin(), stamp_.end(), 0u);
        generation_ = 1;
    }
}

void GroundPlaneCache::set_observer(double x, double y, double microphone_height)
{
    // The horizon depends on the microphone's height, so it belongs to the
    // observer, and the generation stamp that rewinds the sums rewinds it.
    microphone_ = microphone_height;
    set_observer(x, y);
}

void GroundPlaneCache::set_observer_at(double x, double y, double microphone_z)
{
    set_observer(x, y);
    // A microphone the file puts below the raster's ground would sit under its
    // own profile and see every path blocked; it is taken at the ground.
    microphone_ = std::max(0.0, microphone_z - z_);
}

void GroundPlaneCache::ensure(int ray, std::size_t n)
{
    auto& sum_z = sum_z_[ray];
    auto& sum_kz = sum_kz_[ray];
    if (stamp_[ray] != generation_)
    {
        stamp_[ray] = generation_;
        sum_z.clear();
        sum_kz.clear();
        if (screening_)
        {
            z_raw_[ray].clear();
            horizon_[ray].clear();
        }
    }
    if (sum_z.size() > n)
        return;
    // Sample ahead a little, so a track walking outward along a ray does not
    // come back for one sample at a time.
    const std::size_t want = n + 1 + std::min<std::size_t>(n / 4, 64);
    const double dx = step_ * cos_[ray];
    const double dy = step_ * sin_[ray];
    for (std::size_t k = sum_z.size(); k < want; k++)
    {
        const double s = static_cast<double>(k);
        const double z = raster_->height(x_ + s * dx, y_ + s * dy) - z_;
        sum_z.push_back((k > 0 ? sum_z[k - 1] : 0.0) + z);
        sum_kz.push_back((k > 0 ? sum_kz[k - 1] : 0.0) + s * z);
        if (screening_)
        {
            auto& horizon = horizon_[ray];
            z_raw_[ray].push_back(z);
            horizon.push_back(k > 0 ? std::max(horizon[k - 1], sight_slope(z, s * step_))
                                    : -std::numeric_limits<double>::infinity());
        }
    }
}

GroundLine GroundPlaneCache::fit(int ray, std::size_t n)
{
    return line_from_sums(sum_z_[ray][n], sum_kz_[ray][n], n, step_, z_);
}

GroundLine GroundPlaneCache::along(int ray, std::size_t n, double fraction)
{
    ensure(ray, fraction > 0.0 ? n + 1 : n);
    const GroundLine a = fit(ray, n);
    if (!(fraction > 0.0))
        return a;
    const GroundLine b = fit(ray, n + 1);
    return {(1.0 - fraction) * a.intercept + fraction * b.intercept,
            (1.0 - fraction) * a.slope + fraction * b.slope};
}

GroundPlaneCache::RayPair GroundPlaneCache::locate(double dx, double dy, double distance) const
{
    RayPair pair;
    double angle = std::atan2(dy, dx);
    if (angle < 0.0)
        angle += 2.0 * pi;
    const double u = angle * (bearings / (2.0 * pi));
    pair.ray = std::min(static_cast<int>(u), bearings - 1);
    pair.turn = std::clamp(u - pair.ray, 0.0, 1.0);
    pair.next = pair.ray + 1 == bearings ? 0 : pair.ray + 1;

    const double v = distance / step_;
    pair.n = static_cast<std::size_t>(v);
    pair.fraction = v - static_cast<double>(pair.n);
    if (pair.n < 1)
    {
        // Shorter than one step: the line over the first step along its bearing.
        pair.n = 1;
        pair.fraction = 0.0;
    }
    pair.distance = distance;
    return pair;
}

GroundLine GroundPlaneCache::line_from(const RayPair& pair)
{
    const GroundLine a = along(pair.ray, pair.n, pair.fraction);
    if (!(pair.turn > 0.0))
        return a;
    const GroundLine b = along(pair.next, pair.n, pair.fraction);
    const double turn = pair.turn;
    return {(1.0 - turn) * a.intercept + turn * b.intercept, (1.0 - turn) * a.slope + turn * b.slope};
}

GroundLine GroundPlaneCache::line_to(double dx, double dy, double distance)
{
    return line_from(locate(dx, dy, distance));
}

std::size_t GroundPlaneCache::interior_samples(const RayPair& pair) const
{
    // ceil(D / step) - 1: the samples strictly between the microphone and the
    // footprint.  line_from() has already, or will have, sampled each ray to
    // floor(D / step), plus one past a fraction, which covers these and the
    // one after.
    const double v = pair.distance / step_;
    if (!(v > 1.0))
        return 0;
    return static_cast<std::size_t>(std::ceil(v)) - 1;
}

double GroundPlaneCache::source_height(const RayPair& pair, double source_z) const
{
    // The blended profile at the footprint, between the samples around it.  A
    // state a few meters under the terrain is let through upstream
    // (process_trajectory), and is taken onto it, as terrain_ray() takes a
    // height below the mean line onto the line.
    const auto at = [&](std::size_t k) {
        const double a = z_raw_[pair.ray][k];
        if (!(pair.turn > 0.0))
            return a;
        return (1.0 - pair.turn) * a + pair.turn * z_raw_[pair.next][k];
    };
    double foot = at(pair.n);
    if (pair.fraction > 0.0)
        foot = (1.0 - pair.fraction) * foot + pair.fraction * at(pair.n + 1);
    return std::max(source_z - z_, foot);
}

bool GroundPlaneCache::may_block(const RayPair& pair, double source_z)
{
    assert(screening_);
    const std::size_t interior = interior_samples(pair);
    if (interior == 0)
        return false;
    ensure(pair.ray, interior + 1);
    if (pair.turn > 0.0)
        ensure(pair.next, interior + 1);
    const double sight = sight_slope(source_height(pair, source_z), pair.distance);
    // A blended sample's slope is never above the larger of its two rays'.
    if (horizon_[pair.ray][interior] > sight)
        return true;
    return pair.turn > 0.0 && horizon_[pair.next][interior] > sight;
}

bool GroundPlaneCache::blocked(const RayPair& pair, double source_z, Profile& profile)
{
    if (!may_block(pair, source_z))
        return false;
    const std::size_t interior = interior_samples(pair);
    const double source = source_height(pair, source_z);
    const double sight = sight_slope(source, pair.distance);
    const double* a = z_raw_[pair.ray].data();
    const double* b = z_raw_[pair.next].data();
    const bool turned = pair.turn > 0.0;
    const double turn = pair.turn;
    bool above = false;
    for (std::size_t k = 1; k <= interior && !above; k++)
    {
        const double z = turned ? (1.0 - turn) * a[k] + turn * b[k] : a[k];
        const double s = static_cast<double>(k) * step_;
        above = sight_slope(z, s) > sight && (min_protrusion_ == 0.0 || z - (microphone_ + sight * s) > min_protrusion_);
    }
    if (!above)
        return false;
    profile.heights.resize(interior + 2);
    for (std::size_t k = 0; k <= interior + 1; k++)
        profile.heights[k] = turned ? (1.0 - turn) * a[k] + turn * b[k] : a[k];
    profile.interior = interior;
    profile.step = step_;
    profile.distance = pair.distance;
    profile.microphone_height = microphone_;
    profile.source_height = source;
    return true;
}

GroundLine GroundPlaneCache::fit_ray_range(int ray, std::size_t k0, std::size_t k1) const
{
    // Re-based to k0, so the index-weighted sum is over (k - k0) and the
    // closed form of line_from_sums applies; then moved back to s from the
    // observer.
    const auto& sum_z = sum_z_[ray];
    const auto& sum_kz = sum_kz_[ray];
    const double z = sum_z[k1] - (k0 > 0 ? sum_z[k0 - 1] : 0.0);
    const double kz = sum_kz[k1] - (k0 > 0 ? sum_kz[k0 - 1] : 0.0) - static_cast<double>(k0) * z;
    GroundLine line = line_from_sums(z, kz, k1 - k0, step_, 0.0);
    line.intercept -= line.slope * (static_cast<double>(k0) * step_);
    return line;
}

GroundLine GroundPlaneCache::fit_range(const RayPair& pair, std::size_t k0, std::size_t k1) const
{
    assert(screening_ && k0 < k1 && k1 < sum_z_[pair.ray].size());
    const GroundLine a = fit_ray_range(pair.ray, k0, k1);
    if (!(pair.turn > 0.0))
        return a;
    assert(k1 < sum_z_[pair.next].size());
    const GroundLine b = fit_ray_range(pair.next, k0, k1);
    const double turn = pair.turn;
    return {(1.0 - turn) * a.intercept + turn * b.intercept, (1.0 - turn) * a.slope + turn * b.slope};
}
