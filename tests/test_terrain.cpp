#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <random>
#include <string>
#include <vector>

#include <netcdf.h>

#include "propagation.hpp"
#include "terrain.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace)

namespace {

// A grid whose heights are linear in latitude and longitude, so bilinear
// interpolation has to reproduce them exactly: 1000 m per degree of latitude,
// 500 m per degree of longitude, 300 m at the south-west corner.
std::string write_linear_dem(const std::string& name)
{
    const std::string path = (std::filesystem::temp_directory_path() / name).string();
    const std::size_t n_lat = 6;
    const std::size_t n_lon = 6;
    std::vector<double> latitude(n_lat), longitude(n_lon), elevation(n_lat * n_lon);
    for (std::size_t i = 0; i < n_lat; i++)
        latitude[i] = 40.0 + 0.1 * static_cast<double>(i);
    for (std::size_t j = 0; j < n_lon; j++)
        longitude[j] = -77.5 + 0.1 * static_cast<double>(j);
    for (std::size_t i = 0; i < n_lat; i++)
        for (std::size_t j = 0; j < n_lon; j++)
            elevation[i * n_lon + j] = 300.0 + 1000.0 * (latitude[i] - 40.0) + 500.0 * (longitude[j] + 77.5);

    int ncid = 0;
    REQUIRE(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid) == NC_NOERR);
    int lat_dim = 0;
    int lon_dim = 0;
    REQUIRE(nc_def_dim(ncid, "latitude", n_lat, &lat_dim) == NC_NOERR);
    REQUIRE(nc_def_dim(ncid, "longitude", n_lon, &lon_dim) == NC_NOERR);
    int lat_var = 0;
    int lon_var = 0;
    int elev_var = 0;
    REQUIRE(nc_def_var(ncid, "latitude", NC_DOUBLE, 1, &lat_dim, &lat_var) == NC_NOERR);
    REQUIRE(nc_def_var(ncid, "longitude", NC_DOUBLE, 1, &lon_dim, &lon_var) == NC_NOERR);
    const int elev_dims[2] = {lat_dim, lon_dim};
    REQUIRE(nc_def_var(ncid, "elevation", NC_DOUBLE, 2, elev_dims, &elev_var) == NC_NOERR);
    const std::string datum = "NAVD88";
    REQUIRE(nc_put_att_text(ncid, NC_GLOBAL, "vertical_datum", datum.size(), datum.c_str()) == NC_NOERR);
    REQUIRE(nc_enddef(ncid) == NC_NOERR);
    REQUIRE(nc_put_var_double(ncid, lat_var, latitude.data()) == NC_NOERR);
    REQUIRE(nc_put_var_double(ncid, lon_var, longitude.data()) == NC_NOERR);
    REQUIRE(nc_put_var_double(ncid, elev_var, elevation.data()) == NC_NOERR);
    REQUIRE(nc_close(ncid) == NC_NOERR);
    return path;
}

// Write a terrain file: `height(latitude, longitude)` at every node of the
// given axes, with optional attributes on elevation.  Fill values are written
// as `fill` where `height` returns NaN.
struct TerrainFileOptions
{
    std::string units = "m";
    bool with_fill = false;
    double fill = -9999.0;
    double scale = 1.0;
    double offset = 0.0;
    bool packed = false;
    std::string datum = "NAVD88";
};

void write_terrain_file(const std::string& path, const std::vector<double>& latitude,
        const std::vector<double>& longitude, const std::function<double(double, double)>& height,
        const TerrainFileOptions& options = {})
{
    std::vector<double> elevation(latitude.size() * longitude.size());
    for (std::size_t i = 0; i < latitude.size(); i++)
        for (std::size_t j = 0; j < longitude.size(); j++)
        {
            const double z = height(latitude[i], longitude[j]);
            elevation[i * longitude.size() + j] = std::isnan(z) && options.with_fill
                    ? options.fill
                    : (options.packed ? (z - options.offset) / options.scale : z);
        }
    int ncid = 0;
    REQUIRE(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid) == NC_NOERR);
    int lat_dim = 0, lon_dim = 0, lat_var = 0, lon_var = 0, elev_var = 0;
    REQUIRE(nc_def_dim(ncid, "latitude", latitude.size(), &lat_dim) == NC_NOERR);
    REQUIRE(nc_def_dim(ncid, "longitude", longitude.size(), &lon_dim) == NC_NOERR);
    REQUIRE(nc_def_var(ncid, "latitude", NC_DOUBLE, 1, &lat_dim, &lat_var) == NC_NOERR);
    REQUIRE(nc_def_var(ncid, "longitude", NC_DOUBLE, 1, &lon_dim, &lon_var) == NC_NOERR);
    const int dims[2] = {lat_dim, lon_dim};
    REQUIRE(nc_def_var(ncid, "elevation", NC_DOUBLE, 2, dims, &elev_var) == NC_NOERR);
    if (!options.units.empty())
        REQUIRE(nc_put_att_text(ncid, elev_var, "units", options.units.size(), options.units.c_str()) == NC_NOERR);
    if (options.with_fill)
        REQUIRE(nc_put_att_double(ncid, elev_var, "_FillValue", NC_DOUBLE, 1, &options.fill) == NC_NOERR);
    if (options.packed)
    {
        REQUIRE(nc_put_att_double(ncid, elev_var, "scale_factor", NC_DOUBLE, 1, &options.scale) == NC_NOERR);
        REQUIRE(nc_put_att_double(ncid, elev_var, "add_offset", NC_DOUBLE, 1, &options.offset) == NC_NOERR);
    }
    REQUIRE(nc_put_att_text(ncid, NC_GLOBAL, "vertical_datum", options.datum.size(), options.datum.c_str())
            == NC_NOERR);
    REQUIRE(nc_enddef(ncid) == NC_NOERR);
    REQUIRE(nc_put_var_double(ncid, lat_var, latitude.data()) == NC_NOERR);
    REQUIRE(nc_put_var_double(ncid, lon_var, longitude.data()) == NC_NOERR);
    REQUIRE(nc_put_var_double(ncid, elev_var, elevation.data()) == NC_NOERR);
    REQUIRE(nc_close(ncid) == NC_NOERR);
}

std::vector<double> axis(double from, double to, std::size_t n)
{
    std::vector<double> values(n);
    for (std::size_t k = 0; k < n; k++)
        values[k] = from + (to - from) * static_cast<double>(k) / static_cast<double>(n - 1);
    return values;
}

// A raster from a height function of local x, y, every `step` meters over
// [-half, half] in both.
TerrainRaster raster_of(const std::function<double(double, double)>& height, double half, double step)
{
    const std::size_t n = static_cast<std::size_t>(std::llround(2.0 * half / step)) + 1;
    std::vector<double> heights(n * n);
    for (std::size_t j = 0; j < n; j++)
        for (std::size_t i = 0; i < n; i++)
            heights[j * n + i] = height(-half + static_cast<double>(i) * step, -half + static_cast<double>(j) * step);
    return TerrainRaster(-half, -half, step, n, n, std::move(heights));
}

// Rolling ground: a few hills and a tilt, smooth on the scale of a step.
double rolling(double x, double y)
{
    return 12.0 + 0.02 * x - 0.035 * y + 15.0 * std::exp(-((x - 400.0) * (x - 400.0) + (y + 250.0) * (y + 250.0)) / 1.2e5)
            - 9.0 * std::exp(-((x + 700.0) * (x + 700.0) + (y - 600.0) * (y - 600.0)) / 2.0e5)
            + 2.0 * std::sin(x / 130.0) * std::cos(y / 170.0);
}

} // namespace

TEST_CASE("Terrain reads a grid and interpolates it bilinearly", "[terrain]")
{
    const auto path = write_linear_dem("niceops_linear_dem.nc");
    const Terrain terrain(path);
    CHECK(terrain.vertical_datum() == "NAVD88");

    // Between nodes a linear grid is reproduced exactly.
    CHECK(terrain.elevation(40.25, -77.25) == Catch::Approx(300.0 + 250.0 + 125.0).margin(1e-9));
    CHECK(terrain.elevation(40.0, -77.5) == Catch::Approx(300.0).margin(1e-9));
    CHECK(terrain.elevation(40.5, -77.0) == Catch::Approx(300.0 + 500.0 + 250.0).margin(1e-9));

    // Outside the grid there is no height, and the grid says so.
    CHECK(!terrain.covers(39.9, -77.25));
    CHECK(std::isnan(terrain.elevation(39.9, -77.25)));
    CHECK(std::isnan(terrain.elevation(40.25, -76.9)));

    std::filesystem::remove(path);
}

TEST_CASE("A ground plane is the least-squares line through the heights", "[terrain]")
{
    std::vector<double> s;
    std::vector<double> exact;
    std::vector<double> bumpy;
    for (int k = 0; k <= 100; k++)
    {
        const double distance = 25.0 * static_cast<double>(k);
        s.push_back(distance);
        exact.push_back(120.0 + 0.04 * distance);
        bumpy.push_back(120.0 + 0.04 * distance + 3.0 * std::sin(distance / 150.0));
    }

    const auto flat_line = fit_ground_plane(s, exact);
    CHECK(flat_line.slope == Catch::Approx(0.04).margin(1e-12));
    CHECK(flat_line.intercept == Catch::Approx(120.0).margin(1e-9));
    CHECK(flat_line.rms == Catch::Approx(0.0).margin(1e-9));

    // Bumps leave a residual: the fit is the mean line, and the rms says how
    // far the ground departs from it.
    const auto rough = fit_ground_plane(s, bumpy);
    CHECK(rough.slope == Catch::Approx(0.04).margin(1e-3));
    CHECK(rough.rms > 1.0);
    CHECK(rough.rms < 3.0);
}

TEST_CASE("A terrain profile follows the grid along a path", "[terrain]")
{
    const auto path = write_linear_dem("niceops_profile_dem.nc");
    const Terrain terrain(path);
    const LocalFrame frame(Geodetic {40.25, -77.25, 0.0});
    const Vec3 from = frame.to_local(Geodetic {40.1, -77.4, 0.0});
    const Vec3 to = frame.to_local(Geodetic {40.4, -77.1, 0.0});
    const auto profile = terrain_profile(terrain, frame, from, to, 1000.0);

    REQUIRE(profile.size() >= 2);
    CHECK(profile.front().distance == Catch::Approx(0.0));
    CHECK(profile.back().distance == Catch::Approx((to - from).norm()));
    CHECK(profile.front().height == Catch::Approx(terrain.elevation(40.1, -77.4)).margin(1e-6));
    CHECK(profile.back().height == Catch::Approx(terrain.elevation(40.4, -77.1)).margin(1e-6));

    // A linear grid over a path of about 40 km is close to a straight line in
    // local meters, so the profile fits one to well under a meter.
    std::vector<double> distances;
    std::vector<double> heights;
    for (const auto& sample : profile)
    {
        distances.push_back(sample.distance);
        heights.push_back(sample.height);
    }
    CHECK(fit_ground_plane(distances, heights).rms < 1.0);

    std::filesystem::remove(path);
}

TEST_CASE("The terrain reader honors fill values, units and packing", "[terrain]")
{
    const auto dir = std::filesystem::temp_directory_path();
    const auto lat = axis(40.0, 40.1, 5);
    const auto lon = axis(-77.1, -77.0, 5);
    const auto plane = [](double la, double lo) { return 300.0 + 1000.0 * (la - 40.0) + 500.0 * (lo + 77.1); };

    SECTION("a fill value is a missing node, not a height")
    {
        const auto path = (dir / "niceops_terrain_fill.nc").string();
        TerrainFileOptions options;
        options.with_fill = true;
        write_terrain_file(path, lat, lon,
                [&](double la, double lo) { return la > 40.07 && lo > -77.03 ? std::nan("") : plane(la, lo); },
                options);
        const Terrain terrain(path);
        CHECK(terrain.elevation(40.01, -77.09) == Catch::Approx(plane(40.01, -77.09)).margin(1e-9));
        CHECK(std::isnan(terrain.elevation(40.09, -77.01)));
        std::filesystem::remove(path);
    }

    SECTION("heights in feet are refused, not read as meters")
    {
        const auto path = (dir / "niceops_terrain_feet.nc").string();
        TerrainFileOptions options;
        options.units = "ft";
        write_terrain_file(path, lat, lon, plane, options);
        CHECK_THROWS_WITH(Terrain(path), Catch::Matchers::ContainsSubstring("meters"));
        std::filesystem::remove(path);
    }

    SECTION("a packed grid is unpacked")
    {
        const auto path = (dir / "niceops_terrain_packed.nc").string();
        TerrainFileOptions options;
        options.packed = true;
        options.scale = 0.5;
        options.offset = 200.0;
        write_terrain_file(path, lat, lon, plane, options);
        const Terrain terrain(path);
        CHECK(terrain.elevation(40.05, -77.05) == Catch::Approx(plane(40.05, -77.05)).margin(1e-9));
        std::filesystem::remove(path);
    }

    SECTION("a decreasing axis and a one-node axis are refused")
    {
        CHECK_THROWS_WITH(Terrain({40.1, 40.0}, {-77.1, -77.0}, {1.0, 2.0, 3.0, 4.0}),
                Catch::Matchers::ContainsSubstring("increase"));
        CHECK_THROWS_WITH(Terrain({40.0}, {-77.1, -77.0}, {1.0, 2.0}),
                Catch::Matchers::ContainsSubstring("two nodes"));
        CHECK_THROWS(Terrain({40.0, 40.1}, {-77.1, -77.0}, {1.0, 2.0, 3.0}));
    }

    SECTION("the spacing is stated in meters")
    {
        // 0.025 degrees: 2.78 km north-south, 2.13 km east-west at 40 N.
        const Terrain terrain(lat, lon, std::vector<double>(25, 0.0));
        CHECK(terrain.spacing_meters() == Catch::Approx(0.025 * 111412.84 * std::cos(40.05 * 3.14159265358979 / 180.0)).epsilon(2e-3));
    }
}

TEST_CASE("The raster resamples the terrain onto the local frame", "[terrain]")
{
    // Heights defined as a plane in the local frame and written at the grid's
    // nodes: bilinear between nodes 0.001 degrees apart follows the plane to
    // well under a millimeter, so the raster must reproduce it.
    const LocalFrame frame = LocalFrame::track_aligned(Geodetic {40.25, -77.25, 0.0}, 210.0);
    const auto in_frame = [&](double la, double lo) {
        const Vec3 p = frame.to_local(Geodetic {la, lo, 0.0});
        return 150.0 + 0.04 * p.x - 0.03 * p.y;
    };
    const Terrain terrain(axis(40.22, 40.28, 61), axis(-77.29, -77.21, 81), [&] {
        std::vector<double> z;
        for (const double la : axis(40.22, 40.28, 61))
            for (const double lo : axis(-77.29, -77.21, 81))
                z.push_back(in_frame(la, lo));
        return z;
    }());
    const TerrainRaster raster(terrain, frame, -1500.0, 1500.0, -1000.0, 1000.0, 10.0, -100.0);
    for (const auto& [x, y] : std::vector<std::pair<double, double>> {{0.0, 0.0}, {-1234.5, 777.7}, {1499.0, -999.0}, {33.3, 3.3}})
        CHECK(raster.height(x, y) == Catch::Approx(50.0 + 0.04 * x - 0.03 * y).margin(1e-3));

    // Asking for more than the terrain covers is refused with the extent named.
    CHECK_THROWS_WITH(TerrainRaster(terrain, frame, -9000.0, 1500.0, -1000.0, 1000.0, 10.0, 0.0),
            Catch::Matchers::ContainsSubstring("covers latitudes"));
    // So is a void inside it.
    std::vector<double> holed;
    for (const double la : axis(40.22, 40.28, 61))
        for (const double lo : axis(-77.29, -77.21, 81))
            holed.push_back(std::fabs(la - 40.25) < 0.002 && std::fabs(lo + 77.25) < 0.002 ? std::nan("") : 0.0);
    const Terrain void_terrain(axis(40.22, 40.28, 61), axis(-77.29, -77.21, 81), holed);
    CHECK_THROWS_WITH(TerrainRaster(void_terrain, frame, -500.0, 500.0, -500.0, 500.0, 10.0, 0.0),
            Catch::Matchers::ContainsSubstring("no height"));
}

TEST_CASE("The cached mean ground line is the direct least-squares fit", "[terrain]")
{
    const TerrainRaster raster = raster_of(rolling, 2500.0, 10.0);
    GroundPlaneCache cache(raster);

    // Along a ray's own bearing, a whole number of steps out, the cache's line
    // is the direct fit over the same samples.
    const double ox = -137.0, oy = 211.0;
    cache.set_observer(ox, oy);
    for (const int ray : {0, 1, 359, 720, 1000, 1439})
    {
        const double angle = 2.0 * 3.14159265358979323846 * ray / GroundPlaneCache::bearings;
        for (const std::size_t n : {1, 2, 7, 150})
        {
            const double distance = 10.0 * static_cast<double>(n);
            const double dx = distance * std::cos(angle), dy = distance * std::sin(angle);
            const GroundLine cached = cache.line_to(dx, dy, distance);
            const GroundLine direct = fit_ground_line(raster, ox, oy, ox + dx, oy + dy, n);
            CHECK(cached.intercept == Catch::Approx(direct.intercept).margin(1e-7));
            CHECK(cached.slope == Catch::Approx(direct.slope).margin(1e-9));
        }
    }

    // Anywhere else it is the interpolation of fits around the path, which on
    // rolling ground stays within centimeters of the direct fit at both ends.
    std::mt19937 random(7);
    std::uniform_real_distribution<double> where(-1500.0, 1500.0);
    double worst = 0.0;
    for (int observer = 0; observer < 20; observer++)
    {
        const double x0 = where(random), y0 = where(random);
        cache.set_observer(x0, y0);
        for (int pair = 0; pair < 50; pair++)
        {
            const double x1 = where(random), y1 = where(random);
            const double distance = std::hypot(x1 - x0, y1 - y0);
            const GroundLine cached = cache.line_to(x1 - x0, y1 - y0, distance);
            const auto n = std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(distance / 10.0)));
            const GroundLine direct = fit_ground_line(raster, x0, y0, x1, y1, n);
            worst = std::max(worst, std::fabs(cached.intercept - direct.intercept));
            worst = std::max(worst, std::fabs(cached.intercept + cached.slope * distance
                                            - (direct.intercept + direct.slope * distance)));
        }
    }
    CHECK(worst < 0.15);
}

TEST_CASE("Over a uniformly tilted plane the ground models get the plane's own geometry", "[terrain]")
{
    // z = 40 + 0.06 x - 0.025 y, a 3.7 degree slope.  Along any path the mean
    // ground line is the plane's own section, so the heights above it and the
    // distance along it are known in closed form.
    const double gx = 0.06, gy = -0.025;
    const auto plane = [&](double x, double y) { return 40.0 + gx * x + gy * y; };
    const TerrainRaster raster = raster_of(plane, 3000.0, 10.0);
    GroundPlaneCache cache(raster);
    struct Pair
    {
        Vec3 observer, source;
        double microphone;
    };
    const std::vector<Pair> pairs = {
            {{0.0, 0.0, 0.0}, {800.0, 300.0, 400.0}, 1.2},
            {{-500.0, 250.0, 0.0}, {-500.0, -1500.0, 150.0}, 1.2},
            {{1000.0, -900.0, 0.0}, {-1200.0, 1300.0, 600.0}, 4.0},
            {{300.0, 300.0, 0.0}, {303.0, 304.0, 120.0}, 1.2}, // nearly overhead
    };
    for (const auto& p : pairs)
    {
        const Vec3 ground(p.observer.x, p.observer.y, plane(p.observer.x, p.observer.y));
        const Vec3 microphone = ground + Vec3(0.0, 0.0, p.microphone);
        const Vec3 source = p.source + Vec3(0.0, 0.0, plane(p.source.x, p.source.y));
        cache.set_observer(ground.x, ground.y);
        const double dx = source.x - ground.x, dy = source.y - ground.y;
        const double distance = std::hypot(dx, dy);
        const GroundLine line = cache.line_to(dx, dy, distance);

        const double slope = (gx * dx + gy * dy) / distance;
        const double norm = std::sqrt(1.0 + slope * slope);
        CHECK(line.slope == Catch::Approx(slope).margin(2e-7));
        CHECK(line.intercept == Catch::Approx(ground.z).margin(1e-6));

        const RayPath point = terrain_ray(source, microphone, line, false);
        const double above = source.z - (ground.z + slope * distance);
        CHECK(point.source_height == Catch::Approx(above / norm).margin(1e-3));
        CHECK(point.receiver_height == Catch::Approx(p.microphone / norm).margin(1e-6));
        CHECK(point.ground_distance == Catch::Approx((distance + slope * (source.z - microphone.z)) / norm).margin(1e-3));
        // The rest of the record is the straight line's, through the real geometry.
        const RayPath straight = straight_ray(source, microphone);
        CHECK(point.range == straight.range);
        CHECK(point.direction.x == straight.direction.x);
        CHECK(point.separation.z == straight.separation.z);

        // A flush microphone's plate lies along the line, under the microphone.
        const RayPath flush = terrain_ray(source, ground, line, true);
        CHECK(flush.source_height == Catch::Approx(above / norm).margin(1e-3));
        CHECK(std::isnan(flush.receiver_height));
    }
}

TEST_CASE("plane_geometry is a rotation, and the identity on level ground", "[terrain]")
{
    SECTION("level ground at z = 0 changes nothing, to the bit")
    {
        const PlaneGeometry g = plane_geometry(1234.5678, 321.0987, 1.2, GroundLine {0.0, 0.0});
        CHECK(g.source_height == 321.0987);
        CHECK(g.receiver_height == 1.2);
        CHECK(g.ground_distance == 1234.5678);

        const Vec3 source(-812.25, 317.5, 455.125), microphone(12.5, -40.0, 1.2);
        const RayPath straight = straight_ray(source, microphone);
        const RayPath level = terrain_ray(source, microphone, GroundLine {0.0, 0.0}, false);
        CHECK(level.source_height == straight.source_height);
        CHECK(level.ground_distance == straight.ground_distance);
        CHECK(level.receiver_height == microphone.z);
        CHECK(level.range == straight.range);
    }

    SECTION("an upward slope toward the source lowers the source and stretches the base")
    {
        // Ground rising 1 in 10 toward a source 1000 m out at z = 210, the line
        // 10 m up at the observer: hand-computed.
        const PlaneGeometry g = plane_geometry(1000.0, 210.0, 11.2, GroundLine {10.0, 0.1});
        const double norm = std::sqrt(1.01);
        CHECK(g.source_height == Catch::Approx(100.0 / norm).epsilon(1e-14));
        CHECK(g.receiver_height == Catch::Approx(1.2 / norm).epsilon(1e-14));
        CHECK(g.ground_distance == Catch::Approx((1000.0 + 0.1 * 198.8) / norm).epsilon(1e-14));
        // And falling ground does the opposite.
        const PlaneGeometry down = plane_geometry(1000.0, 210.0, 11.2, GroundLine {10.0, -0.1});
        CHECK(down.source_height == Catch::Approx(300.0 / norm).epsilon(1e-14));
        CHECK(down.ground_distance == Catch::Approx((1000.0 - 0.1 * 198.8) / norm).epsilon(1e-14));
    }

    SECTION("the source-microphone distance is unchanged")
    {
        for (const double slope : {-0.3, -0.05, 0.02, 0.4})
        {
            const PlaneGeometry g = plane_geometry(750.0, 180.0, 3.0, GroundLine {-4.0, slope});
            CHECK(std::hypot(g.ground_distance, g.source_height - g.receiver_height)
                    == Catch::Approx(std::hypot(750.0, 177.0)).epsilon(1e-13));
        }
    }

    SECTION("heights below the line are raised to it, and say so")
    {
        // A receiver in a dip: its foot 0.5 m under the mean line, 1.2 m pole.
        bool below = false;
        const RayPath dip = terrain_ray(Vec3(500.0, 0.0, 100.0), Vec3(0.0, 0.0, -2.0), GroundLine {0.0, 0.0},
                false, &below);
        CHECK(below);
        CHECK(dip.receiver_height == 0.0);
        CHECK(dip.source_height == 100.0);
        // A flush microphone above the line keeps its height out of the plate's
        // geometry: the source is measured from the microphone.
        const RayPath plate = terrain_ray(Vec3(500.0, 0.0, 100.0), Vec3(0.0, 0.0, 0.7), GroundLine {0.0, 0.0},
                true, &below);
        CHECK(!below);
        CHECK(plate.source_height == Catch::Approx(99.3).epsilon(1e-14));
        // A source behind rising ground.
        terrain_ray(Vec3(500.0, 0.0, 20.0), Vec3(0.0, 0.0, 1.2), GroundLine {0.0, 0.05}, false, &below);
        CHECK(below);
    }
}

TEST_CASE("write the terrain files the --terrain CTests read", "[.terrain_files]")
{
    // Not a test: CMakeLists.txt's terrain_files runs it (by its hidden tag) to
    // write integrated/terrain_flat.nc and integrated/terrain_hills.nc around
    // the demo track's frame (origin 39.81043 N 74.42631 W, +x along 210
    // degrees, which reproduces the track's own x/y to 2 mm; its z = 0 is at
    // altitude 54.2651 ft).  The flat one is level everywhere at the height
    // --terrain_offset -1250 with --ground_elevation 54.2651 (ft) takes back to
    // z = 0, written as the program forms it, so it lands there to the bit.
    // The hills fall toward the approach's low end and stay under every state.
    const LocalFrame frame = LocalFrame::track_aligned(Geodetic {39.81043, -74.42631, 0.0}, 210.0);
    const auto lat = axis(39.75, 39.87, 241);
    const auto lon = axis(-74.51, -74.34, 341);
    std::filesystem::create_directories("integrated");
    const double flat = -((-1250.0 * 0.3048) - (54.2651 * 0.3048));
    write_terrain_file("integrated/terrain_flat.nc", lat, lon, [&](double, double) { return flat; });
    write_terrain_file("integrated/terrain_hills.nc", lat, lon, [&](double la, double lo) {
        const Vec3 p = frame.to_local(Geodetic {la, lo, 0.0});
        return -20.0 - 0.03 * p.x + 0.5 * rolling(p.x, p.y) - 6.0;
    });
    REQUIRE(Terrain("integrated/terrain_flat.nc").elevation(39.81043, -74.42631) == flat);

    // And arrays/coyle_approach.csv's observers (its first row is read as the
    // header, as the program reads it) standing on the hills, z in feet: the
    // no-terrain run over these is what a --terrain run would be if its ground
    // model ignored the terrain and only stood the observers on it.
    const Terrain hills("integrated/terrain_hills.nc");
    std::ifstream coyle("arrays/coyle_approach.csv");
    REQUIRE(coyle);
    std::ofstream out("integrated/terrain_hills_observers.csv");
    out << "id,x_ft,y_ft,z_ft\n" << std::setprecision(17);
    std::string line;
    std::getline(coyle, line);
    int rows = 0;
    while (std::getline(coyle, line))
    {
        std::stringstream fields(line);
        std::string id, x, y;
        std::getline(fields, id, ',');
        std::getline(fields, x, ',');
        std::getline(fields, y, ',');
        const double x_ft = std::stod(x), y_ft = std::stod(y);
        const Geodetic where = frame.to_geodetic(Vec3(x_ft * 0.3048, y_ft * 0.3048, 0.0));
        const double z = hills.elevation(where.latitude, where.longitude) + (0.0 - 54.2651 * 0.3048);
        out << std::stoi(id) << "," << x_ft << "," << y_ft << "," << z / 0.3048 << "\n";
        rows++;
    }
    REQUIRE(rows == 32);
}
