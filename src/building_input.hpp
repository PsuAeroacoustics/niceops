#pragma once

// The district and microphone files the urban tools share
// (niceops_urban_table, niceops_point_field): a CSV of footprint vertices or
// a GeoJSON of WGS84 footprints, and the observer file's columns by name.

#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "buildings.hpp"
#include "geodetic.hpp"
#include "geojson_buildings.hpp"
#include "propagation_table.hpp"
#include "util/constants.hpp"
#include "util/csv.hpp"

namespace niceops {

// A column by name.  A length column may name its unit, as an observer
// file's do (x_ft, height_m, "width (ft)"): it must then be `units`.
inline std::size_t length_column(const LabeledCsv& csv, const std::string& file, const std::string& units,
        const std::string& name, bool required, bool length)
{
    const std::string ft = name + "_ft", m = name + "_m";
    const std::size_t c = length ? csv.column({name.c_str(), ft.c_str(), m.c_str()}) : csv.column({name.c_str()});
    if (c == LabeledCsv::npos)
    {
        if (required)
            throw std::runtime_error(file + " has no '" + name + "' column");
        return c;
    }
    std::string declared = csv.unit_token(c);
    if (declared.empty() && csv.column_matches(c, {ft.c_str()}))
        declared = "ft";
    if (declared.empty() && csv.column_matches(c, {m.c_str()}))
        declared = "m";
    if (length && !declared.empty() && declared != units)
        throw std::runtime_error(file + ": column '" + csv.header[c] + "' is in " + declared + ", and --units is "
                + units);
    return c;
}

// A CSV of numbers under a header, every row complete.
inline LabeledCsv read_numeric_rows(const std::string& file, const std::string& what)
{
    CsvReport report;
    LabeledCsv csv = read_labeled_csv(file, &report);
    if (report.dropped > 0)
        throw std::runtime_error(file + ": line " + std::to_string(report.first_dropped_line)
                + " is not a row of numbers");
    if (csv.rows.empty())
        throw std::runtime_error(file + " holds no " + what);
    for (std::size_t row = 0; row < csv.rows.size(); row++)
    {
        if (csv.rows[row].size() != csv.header.size())
            throw std::runtime_error(file + ": row " + std::to_string(row + 1) + " has "
                    + std::to_string(csv.rows[row].size()) + " fields where the header names "
                    + std::to_string(csv.header.size()));
        for (std::size_t column = 0; column < csv.rows[row].size(); column++)
            if (!std::isfinite(csv.rows[row][column]))
                throw std::runtime_error(file + ": row " + std::to_string(row + 1) + "'s " + csv.header[column]
                        + " is not a finite number");
    }
    return csv;
}

struct BuildingInputOptions
{
    std::string units = "ft";
    double to_meters = ft2m;
    // GeoJSON: the run's WGS84 origin (latitude, longitude) and the compass
    // bearing of its +x, and geojson_buildings.hpp's heights.
    std::vector<double> origin;
    double frame_bearing = 90.0;
    double level_height = 3.0;
    double default_height = 0.0;
    bool default_height_given = false;
    // Where a building gives none.
    double absorption = 0.0;
    double scattering = 0.0;
};

// The footprints, in metres in the run's frame, and each building's walls'
// absorption and scattering coefficients and its id.
struct BuildingInput
{
    std::vector<Footprint> footprints;
    std::vector<double> absorption;
    std::vector<double> scattering;
    std::vector<std::string> ids;
};

// A GeoJSON file (.geojson or .json) of WGS84 footprints, projected into the
// run's frame, or a CSV with one row per footprint vertex: building (an id,
// a building's rows consecutive, its vertices in order around it), x, y, top
// and optionally absorption and scattering.
inline BuildingInput read_buildings(const std::string& file, const BuildingInputOptions& o)
{
    BuildingInput in;
    auto& footprints = in.footprints;
    auto& building_absorption = in.absorption;
    auto& building_scattering = in.scattering;
    auto& ids = in.ids;
    const std::string extension = std::filesystem::path(file).extension().string();
    if (extension == ".geojson" || extension == ".json")
    {
        // WGS84 footprints, projected into the run's frame.
        if (o.origin.size() != 2)
            throw std::runtime_error(file + " is GeoJSON, in longitude and latitude: give the run's "
                    "--origin (and its --frame_bearing) to project it into the run's frame");
        const LocalFrame frame = LocalFrame::track_aligned(Geodetic {o.origin[0], o.origin[1], 0.0}, o.frame_bearing);
        GeojsonOptions geojson;
        geojson.level_height = o.level_height;
        if (o.default_height_given)
            geojson.default_height = o.default_height;
        GeojsonReport report;
        for (auto& building : read_geojson_buildings(file, frame, geojson, report))
        {
            const double alpha = std::isnan(building.absorption) ? o.absorption : building.absorption;
            const double scatter = std::isnan(building.scattering) ? o.scattering : building.scattering;
            if (!(alpha >= 0.0 && alpha <= 1.0) || !(scatter >= 0.0 && scatter <= 1.0))
                throw std::runtime_error(file + " feature " + building.id
                        + ": absorption and scattering coefficients must lie within [0, 1]");
            footprints.push_back(std::move(building.footprint));
            building_absorption.push_back(alpha);
            building_scattering.push_back(scatter);
            ids.push_back(building.id);
        }
        std::cerr << file << ": " << report.footprints << " footprints from " << report.features
                  << " features";
        if (report.defaulted > 0)
            std::cerr << "; " << report.defaulted << " given --default_height";
        if (report.raised > 0)
            std::cerr << "; " << report.raised << " raised off the ground left out";
        if (report.not_polygons > 0)
            std::cerr << "; " << report.not_polygons << " not polygons left out";
        std::cerr << "." << std::endl;
        if (report.no_height > 0)
            throw std::runtime_error(file + ": " + std::to_string(report.no_height)
                    + " buildings give no height (height, building:height or building:levels); give "
                      "--default_height to take them at one");
        if (footprints.empty())
            throw std::runtime_error(file + " holds no building footprints");
    }
    else
    {
        const LabeledCsv bcsv = read_numeric_rows(file, "buildings");
        const std::size_t cid = length_column(bcsv, file, o.units, "building", true, false);
        const std::size_t bx = length_column(bcsv, file, o.units, "x", true, true);
        const std::size_t by = length_column(bcsv, file, o.units, "y", true, true);
        const std::size_t bt = length_column(bcsv, file, o.units, "top", true, true);
        const std::size_t ba = length_column(bcsv, file, o.units, "absorption", false, false);
        const std::size_t bs = length_column(bcsv, file, o.units, "scattering", false, false);
        std::set<double> seen;
        double current = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t row = 0; row < bcsv.rows.size(); row++)
        {
            const auto& r = bcsv.rows[row];
            const std::string where = file + " row " + std::to_string(row + 1) + ": ";
            const double alpha = ba == LabeledCsv::npos ? o.absorption : r[ba];
            const double scatter = bs == LabeledCsv::npos ? o.scattering : r[bs];
            if (footprints.empty() || r[cid] != current)
            {
                if (!seen.insert(r[cid]).second)
                    throw std::runtime_error(where + "building " + exact_text(r[cid])
                            + "'s rows are not consecutive");
                if (!(alpha >= 0.0 && alpha <= 1.0))
                    throw std::runtime_error(where + "an absorption coefficient must lie within [0, 1]");
                if (!(scatter >= 0.0 && scatter <= 1.0))
                    throw std::runtime_error(where + "a scattering coefficient must lie within [0, 1]");
                footprints.emplace_back();
                footprints.back().top = r[bt] * o.to_meters;
                building_absorption.push_back(alpha);
                building_scattering.push_back(scatter);
                current = r[cid];
                ids.push_back(exact_text(r[cid]));
            }
            else if (r[bt] * o.to_meters != footprints.back().top || alpha != building_absorption.back()
                    || scatter != building_scattering.back())
                throw std::runtime_error(
                        where + "a building's rows must give one top, one absorption and one scattering");
            footprints.back().x.push_back(r[bx] * o.to_meters);
            footprints.back().y.push_back(r[by] * o.to_meters);
        }
    }
    return in;
}

} // namespace niceops
