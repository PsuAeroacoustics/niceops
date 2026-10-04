#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <vector>

#include <netcdf.h>

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
