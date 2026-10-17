#include "geojson_buildings.hpp"

#include <cctype>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

namespace niceops {

namespace {

using json = nlohmann::json;

std::string lower(std::string s)
{
    for (char& c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// A property as metres, or NaN: a number, or text geojson_height reads.
double height_property(const json& properties, const char* key)
{
    const auto it = properties.find(key);
    if (it == properties.end() || it->is_null())
        return std::numeric_limits<double>::quiet_NaN();
    if (it->is_number())
        return it->get<double>();
    if (it->is_string())
        return geojson_height(it->get<std::string>());
    return std::numeric_limits<double>::quiet_NaN();
}

double number_property(const json& properties, const char* key)
{
    const auto it = properties.find(key);
    if (it == properties.end() || it->is_null())
        return std::numeric_limits<double>::quiet_NaN();
    if (it->is_number())
        return it->get<double>();
    if (it->is_string())
    {
        try
        {
            std::size_t used = 0;
            const double v = std::stod(it->get<std::string>(), &used);
            return used == it->get<std::string>().size() ? v : std::numeric_limits<double>::quiet_NaN();
        }
        catch (const std::exception&)
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

// One ring, projected: the closing vertex dropped, as is any vertex within a
// millimetre of the one before it.
Ring project_ring(const json& coordinates, const LocalFrame& frame, const std::string& where)
{
    if (!coordinates.is_array())
        throw std::runtime_error(where + ": a ring is not an array of positions");
    Ring ring;
    for (const json& position : coordinates)
    {
        if (!position.is_array() || position.size() < 2 || !position[0].is_number() || !position[1].is_number())
            throw std::runtime_error(where + ": a position is not [longitude, latitude]");
        const double longitude = position[0].get<double>(), latitude = position[1].get<double>();
        if (!(std::abs(latitude) <= 90.0 && std::abs(longitude) <= 180.0))
            throw std::runtime_error(where + ": a position is not a longitude and latitude in degrees");
        const Vec3 local = frame.to_local(Geodetic {latitude, longitude, 0.0});
        if (!ring.x.empty() && std::hypot(local.x - ring.x.back(), local.y - ring.y.back()) < 1e-3)
            continue;
        ring.x.push_back(local.x);
        ring.y.push_back(local.y);
    }
    while (ring.x.size() > 1 && std::hypot(ring.x.front() - ring.x.back(), ring.y.front() - ring.y.back()) < 1e-3)
    {
        ring.x.pop_back();
        ring.y.pop_back();
    }
    return ring;
}

} // namespace

double geojson_height(const std::string& text)
{
    std::size_t used = 0;
    double value;
    try
    {
        value = std::stod(text, &used);
    }
    catch (const std::exception&)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::string unit = lower(text.substr(used));
    unit.erase(0, unit.find_first_not_of(" \t"));
    unit.erase(unit.find_last_not_of(" \t") + 1);
    if (unit.empty() || unit == "m" || unit == "metre" || unit == "metres" || unit == "meter" || unit == "meters")
        return value;
    if (unit == "ft" || unit == "feet" || unit == "foot" || unit == "'")
        return value * 0.3048;
    return std::numeric_limits<double>::quiet_NaN();
}

std::vector<GeojsonBuilding> read_geojson_buildings(const std::string& path, const LocalFrame& frame,
        const GeojsonOptions& options, GeojsonReport& report)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("cannot open " + path);
    json document;
    try
    {
        in >> document;
    }
    catch (const json::exception& e)
    {
        throw std::runtime_error(path + " is not JSON: " + e.what());
    }
    if (const auto crs = document.find("crs"); crs != document.end())
    {
        // RFC 7946 has no crs member; the 2008 specification's must name
        // WGS84 longitude and latitude.
        const std::string name = lower(crs->value("properties", json::object()).value("name", std::string()));
        if (name.find("crs84") == std::string::npos && name.find("4326") == std::string::npos)
            throw std::runtime_error(path + ": its coordinates are in " + (name.empty() ? "an unnamed system" : name)
                    + ", not WGS84 longitude and latitude");
    }
    std::vector<json> features;
    const std::string type = document.value("type", std::string());
    if (type == "FeatureCollection" && document.contains("features") && document["features"].is_array())
        features.assign(document["features"].begin(), document["features"].end());
    else if (type == "Feature")
        features.push_back(document);
    else
        throw std::runtime_error(path + " is not a GeoJSON FeatureCollection or Feature");

    std::vector<GeojsonBuilding> out;
    for (std::size_t k = 0; k < features.size(); k++)
    {
        const json& feature = features[k];
        report.features++;
        std::string id = std::to_string(k + 1);
        if (const auto it = feature.find("id"); it != feature.end())
            id = it->is_string() ? it->get<std::string>() : it->dump();
        const std::string where = path + " feature " + id;
        const json properties = feature.value("properties", json::object());
        const json geometry = feature.contains("geometry") ? feature["geometry"] : json();
        const std::string kind = geometry.is_object() ? geometry.value("type", std::string()) : std::string();
        if (kind != "Polygon" && kind != "MultiPolygon")
        {
            report.not_polygons++;
            continue;
        }
        if (properties.is_object()
                && (height_property(properties, "min_height") > 0.0
                        || number_property(properties, "building:min_level") > 0.0))
        {
            report.raised++;
            continue;
        }
        double top = std::numeric_limits<double>::quiet_NaN();
        if (properties.is_object())
        {
            top = height_property(properties, "height");
            if (!(top > 0.0))
                top = height_property(properties, "building:height");
            if (!(top > 0.0))
            {
                const double levels = number_property(properties, "building:levels");
                if (levels > 0.0)
                    top = levels * options.level_height;
            }
        }
        if (!(top > 0.0))
        {
            if (!(options.default_height > 0.0))
            {
                report.no_height++;
                continue;
            }
            top = options.default_height;
            report.defaulted++;
        }
        std::vector<json> polygons;
        if (kind == "Polygon")
            polygons.push_back(geometry.at("coordinates"));
        else
            for (const json& polygon : geometry.at("coordinates"))
                polygons.push_back(polygon);
        for (std::size_t p = 0; p < polygons.size(); p++)
        {
            const json& rings = polygons[p];
            if (!rings.is_array() || rings.empty())
                throw std::runtime_error(where + ": a polygon has no rings");
            GeojsonBuilding building;
            building.id = polygons.size() > 1 ? id + "." + std::to_string(p + 1) : id;
            const Ring outer = project_ring(rings[0], frame, where);
            building.footprint.x = outer.x;
            building.footprint.y = outer.y;
            building.footprint.top = top;
            for (std::size_t r = 1; r < rings.size(); r++)
                building.footprint.holes.push_back(project_ring(rings[r], frame, where));
            if (properties.is_object())
            {
                building.absorption = number_property(properties, "absorption");
                building.scattering = number_property(properties, "scattering");
            }
            out.push_back(std::move(building));
            report.footprints++;
        }
    }
    return out;
}

} // namespace niceops
