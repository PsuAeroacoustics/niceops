#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "buildings.hpp"
#include "geodetic.hpp"
#include "geojson_buildings.hpp"

using namespace niceops; // NOLINT

namespace {

std::string write_file(const std::string& name, const std::string& text)
{
    const auto path = (std::filesystem::temp_directory_path() / name).string();
    std::ofstream(path) << text;
    return path;
}

// "[lon, lat]" of a local point, in the frame.
std::string at(const LocalFrame& frame, double x, double y)
{
    const Geodetic g = frame.to_geodetic(Vec3(x, y, 0.0));
    char buffer[96];
    std::snprintf(buffer, sizeof buffer, "[%.12f, %.12f]", g.longitude, g.latitude);
    return buffer;
}

std::string ring(const LocalFrame& frame, std::initializer_list<std::pair<double, double>> points)
{
    std::string out = "[";
    for (const auto& [x, y] : points)
        out += at(frame, x, y) + ", ";
    out += at(frame, points.begin()->first, points.begin()->second) + "]"; // closed, as GeoJSON has it
    return out;
}

} // namespace

TEST_CASE("GeoJSON footprints land in the run's frame, courtyards and all", "[geojson]")
{
    // The demo case's frame: its origin, +x along 210 degrees.
    const LocalFrame frame = LocalFrame::track_aligned(Geodetic {39.8104276, -74.4263086, 0.0}, 210.0);
    const std::string block = ring(frame, {{0, 0}, {60, 0}, {60, 60}, {0, 60}});
    const std::string yard = ring(frame, {{20, 20}, {40, 20}, {40, 40}, {20, 40}});
    const std::string east = ring(frame, {{100, 0}, {120, 0}, {120, 20}, {100, 20}});
    const std::string west = ring(frame, {{-50, 0}, {-30, 0}, {-30, 20}, {-50, 20}});
    const std::string text = R"({"type": "FeatureCollection", "features": [
        {"type": "Feature", "id": "court", "properties": {"height": "15.5 m", "scattering": 0.2},
         "geometry": {"type": "Polygon", "coordinates": [)" + block + ", " + yard + R"(]}},
        {"type": "Feature", "properties": {"building:levels": "4"},
         "geometry": {"type": "MultiPolygon", "coordinates": [[)" + east + "], [" + west + R"(]]}},
        {"type": "Feature", "properties": {"height": 40, "min_height": 30},
         "geometry": {"type": "Polygon", "coordinates": [)" + east + R"(]}},
        {"type": "Feature", "properties": {"building": "yes"},
         "geometry": {"type": "Polygon", "coordinates": [)" + west + R"(]}},
        {"type": "Feature", "properties": {"height": "40 ft"},
         "geometry": {"type": "Point", "coordinates": [0, 0]}}
    ]})";
    const std::string path = write_file("niceops_geojson_test.geojson", text);
    GeojsonOptions options;
    GeojsonReport report;
    const auto buildings = read_geojson_buildings(path, frame, options, report);
    REQUIRE(report.features == 5);
    REQUIRE(report.footprints == 3);
    REQUIRE(report.raised == 1);
    REQUIRE(report.no_height == 1);
    REQUIRE(report.not_polygons == 1);
    REQUIRE(buildings.size() == 3);
    REQUIRE(buildings[0].id == "court");
    REQUIRE(buildings[0].footprint.top == 15.5);
    REQUIRE(buildings[0].scattering == 0.2);
    REQUIRE(std::isnan(buildings[0].absorption));
    REQUIRE(buildings[0].footprint.holes.size() == 1);
    REQUIRE(buildings[0].footprint.x.size() == 4); // the closing vertex dropped
    REQUIRE(buildings[1].footprint.top == 12.0); // four storeys of 3 m
    REQUIRE(buildings[1].id.find(".1") != std::string::npos);
    // The vertices come back where they were put, to well under a millimetre.
    for (std::size_t i = 0; i < 4; i++)
    {
        const double x[] {0, 60, 60, 0}, y[] {0, 0, 60, 60};
        REQUIRE(buildings[0].footprint.x[i] == Catch::Approx(x[i]).margin(1e-6));
        REQUIRE(buildings[0].footprint.y[i] == Catch::Approx(y[i]).margin(1e-6));
    }
    std::vector<Footprint> footprints;
    for (const auto& b : buildings)
        footprints.push_back(b.footprint);
    const Buildings city(footprints);
    REQUIRE(city.inside(30.0, 30.0) == city.size()); // the courtyard
    REQUIRE(city.inside(10.0, 30.0) == 0);

    // With a default height the unmeasured building comes in too.
    options.default_height = 9.0;
    GeojsonReport second;
    const auto all = read_geojson_buildings(path, frame, options, second);
    REQUIRE(second.defaulted == 1);
    REQUIRE(all.size() == 4);
    REQUIRE(all[3].footprint.top == 9.0);
}

TEST_CASE("heights are read as OpenStreetMap writes them", "[geojson]")
{
    REQUIRE(geojson_height("12") == 12.0);
    REQUIRE(geojson_height("12.5 m") == 12.5);
    REQUIRE(geojson_height("40 ft") == Catch::Approx(12.192));
    REQUIRE(geojson_height("40'") == Catch::Approx(12.192));
    REQUIRE(std::isnan(geojson_height("tall")));
    REQUIRE(std::isnan(geojson_height("12;14")));
    REQUIRE(std::isnan(geojson_height("3 storeys")));
}

TEST_CASE("GeoJSON in another reference system is refused", "[geojson]")
{
    const LocalFrame frame(Geodetic {39.8, -74.4, 0.0});
    const std::string path = write_file("niceops_geojson_crs.geojson",
            R"({"type": "FeatureCollection", "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::3857"}},
                "features": []})");
    GeojsonReport report;
    REQUIRE_THROWS_AS(read_geojson_buildings(path, frame, GeojsonOptions {}, report), std::runtime_error);
    const std::string not_geojson = write_file("niceops_geojson_bad.geojson", R"({"type": "Topology"})");
    REQUIRE_THROWS_AS(read_geojson_buildings(not_geojson, frame, GeojsonOptions {}, report), std::runtime_error);
}
