#pragma once

// The district and microphone files the urban tools share
// (niceops_urban_table, niceops_point_field): a CSV of footprint vertices or
// a GeoJSON of WGS84 footprints, and the observer file's columns by name.

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
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

// A scattering coefficient band by band, as tools/profile_scattering.py
// writes it: a CSV of band_hz and s, returned sorted by frequency.
inline std::vector<std::pair<double, double>> read_scattering_spectrum(const std::string& file,
        const std::string& units)
{
    const LabeledCsv rows = read_numeric_rows(file, "scattering coefficients");
    const std::size_t fc = length_column(rows, file, units, "band_hz", true, false);
    const std::size_t sc = length_column(rows, file, units, "s", true, false);
    std::vector<std::pair<double, double>> spectrum;
    for (const auto& r : rows.rows)
    {
        if (!(r[fc] > 0.0) || !(r[sc] >= 0.0 && r[sc] <= 1.0))
            throw std::runtime_error(file + ": each band_hz must be positive and each s within [0, 1]");
        spectrum.emplace_back(r[fc], r[sc]);
    }
    std::sort(spectrum.begin(), spectrum.end());
    for (std::size_t i = 1; i < spectrum.size(); i++)
        if (spectrum[i].first == spectrum[i - 1].first)
            throw std::runtime_error(file + " gives " + exact_text(spectrum[i].first) + " Hz twice");
    return spectrum;
}

// That spectrum at each band, interpolated in log frequency.  A band more
// than 1 % beyond its ends is refused, `file` naming it.
inline std::vector<double> scattering_at(const std::vector<std::pair<double, double>>& spectrum,
        const std::vector<double>& bands, const std::string& file)
{
    std::vector<double> out(bands.size());
    for (std::size_t b = 0; b < bands.size(); b++)
    {
        const double f = bands[b];
        if (f < spectrum.front().first / 1.01 || f > spectrum.back().first * 1.01)
            throw std::runtime_error(file + " covers " + exact_text(spectrum.front().first) + " to "
                    + exact_text(spectrum.back().first) + " Hz, not the band at " + exact_text(f) + " Hz");
        std::size_t i = 1;
        while (i + 1 < spectrum.size() && spectrum[i].first < f)
            i++;
        if (spectrum.size() == 1)
            out[b] = spectrum[0].second;
        else
        {
            const double t = std::clamp(std::log(f / spectrum[i - 1].first)
                            / std::log(spectrum[i].first / spectrum[i - 1].first), 0.0, 1.0);
            out[b] = spectrum[i - 1].second + t * (spectrum[i].second - spectrum[i - 1].second);
        }
    }
    return out;
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
    // A CSV of roof planes (read_roofs), or empty for flat roofs.
    std::string roofs;
};

// The footprints, in metres in the run's frame, and each building's walls'
// absorption and scattering coefficients and its id.
struct BuildingInput
{
    std::vector<Footprint> footprints;
    std::vector<double> absorption;
    std::vector<double> scattering;
    std::vector<std::string> ids;
    // Which of the run's reliefs (facade_relief.hpp) a building's walls carry,
    // counting from 1, 0 for none; -1 where the file has no relief column.
    std::vector<int> relief;
};

// Which of `count` reliefs a building carries: its relief column's, or
// without one the first wherever it scatters.  Throws where the column names
// one past `count`.
inline int relief_of(const BuildingInput& in, std::size_t building, std::size_t count)
{
    const int given = in.relief[building];
    if (given < 0)
        return count > 0 && in.scattering[building] > 0.0 ? 1 : 0;
    if (static_cast<std::size_t>(given) > count)
        throw std::runtime_error("building " + in.ids[building] + " carries relief " + std::to_string(given)
                + ", and the run gives " + std::to_string(count));
    return given;
}

// What is wrong with where a source or microphone stands among the
// buildings, "" for nothing: below the floor, inside a building (under its
// roof, pitched or not), or on a building's edge, where no one surface's
// reflection stands for it (Buildings::support).  On a wall, a level roof or
// a slope is fine: the surface reflects at the point.
inline std::string misplaced_point(const Buildings& city, const std::vector<std::string>& ids, const Vec3& p)
{
    if (!(p.z >= 0.0))
        return "stands below the floor";
    const std::size_t b = city.inside(p.x, p.y);
    if (b < city.size() && p.z < roof_height(city.footprints()[b], p.x, p.y) - Buildings::contact)
        return "stands inside building " + ids[b];
    const Buildings::Support on = city.support(p);
    if (on.kind == Buildings::Support::Kind::edge)
        return "lies on a building's edge, or where two of its surfaces meet: move it off by a millimetre or more";
    return "";
}

// Pitched roofs (LoD2) from a CSV with one row per plane: building (an id of
// the district's), x, y and z, a point the plane passes through, and slope_x
// and slope_y, its rise per unit run along x and along y.  A building's roof
// is the least of its planes, capped at its top: a top below the planes'
// highest point is a flat there (a mansard's deck), and one within 1 mm of it
// or above changes nothing.  A building given no plane keeps its flat roof.
inline void read_roofs(const std::string& file, const BuildingInputOptions& o, BuildingInput& in)
{
    const LabeledCsv csv = read_numeric_rows(file, "roof planes");
    const std::size_t cid = length_column(csv, file, o.units, "building", true, false);
    const std::size_t cx = length_column(csv, file, o.units, "x", true, true);
    const std::size_t cy = length_column(csv, file, o.units, "y", true, true);
    const std::size_t cz = length_column(csv, file, o.units, "z", true, true);
    const std::size_t csx = length_column(csv, file, o.units, "slope_x", true, false);
    const std::size_t csy = length_column(csv, file, o.units, "slope_y", true, false);
    std::map<std::string, std::size_t> index;
    for (std::size_t b = 0; b < in.ids.size(); b++)
        index.emplace(in.ids[b], b);
    for (std::size_t row = 0; row < csv.rows.size(); row++)
    {
        const auto& r = csv.rows[row];
        const auto found = index.find(exact_text(r[cid]));
        if (found == index.end())
            throw std::runtime_error(file + " row " + std::to_string(row + 1) + ": building " + exact_text(r[cid])
                    + " is not in the district");
        const double x = r[cx] * o.to_meters, y = r[cy] * o.to_meters, z = r[cz] * o.to_meters;
        in.footprints[found->second].roof.push_back({r[csx], r[csy], z - r[csx] * x - r[csy] * y});
    }
    for (Footprint& f : in.footprints)
        if (!f.roof.empty() && f.top < roof_top(f) - 1e-3)
            f.roof.push_back({0.0, 0.0, f.top});
}

// A GeoJSON file (.geojson or .json) of WGS84 footprints, projected into the
// run's frame, or a CSV with one row per footprint vertex: building (an id,
// a building's rows consecutive, its vertices in order around it), x, y, top
// and optionally absorption, scattering and relief (a whole number); then
// `o.roofs`'s planes (read_roofs), which a GeoJSON district does not take.
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
            in.relief.push_back(-1);
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
        const std::size_t br = length_column(bcsv, file, o.units, "relief", false, false);
        std::set<double> seen;
        double current = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t row = 0; row < bcsv.rows.size(); row++)
        {
            const auto& r = bcsv.rows[row];
            const std::string where = file + " row " + std::to_string(row + 1) + ": ";
            const double alpha = ba == LabeledCsv::npos ? o.absorption : r[ba];
            const double scatter = bs == LabeledCsv::npos ? o.scattering : r[bs];
            const int relief = br == LabeledCsv::npos ? -1 : static_cast<int>(r[br]);
            if (br != LabeledCsv::npos && !(r[br] >= 0.0 && r[br] == std::round(r[br])))
                throw std::runtime_error(where + "a relief is a whole number, 0 for none");
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
                in.relief.push_back(relief);
                current = r[cid];
                ids.push_back(exact_text(r[cid]));
            }
            else if (r[bt] * o.to_meters != footprints.back().top || alpha != building_absorption.back()
                    || scatter != building_scattering.back() || relief != in.relief.back())
                throw std::runtime_error(
                        where + "a building's rows must give one top, one absorption, one scattering and one relief");
            footprints.back().x.push_back(r[bx] * o.to_meters);
            footprints.back().y.push_back(r[by] * o.to_meters);
        }
    }
    if (!o.roofs.empty())
    {
        if (extension == ".geojson" || extension == ".json")
            throw std::runtime_error(o.roofs + ": roof planes go with a CSV district, in its frame; " + file
                    + " is GeoJSON");
        read_roofs(o.roofs, o, in);
    }
    return in;
}

} // namespace niceops
