#include "ground_classes.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "util/csv.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

// Every column the file may carry, by its normalized name, and as the refusal
// lists them.
const std::vector<std::string>& known_columns()
{
    static const std::vector<std::string> names = {"class", "name", "impedancemodel", "flowresistance",
            "porosityrate", "layerdepth", "roughness"};
    return names;
}
constexpr const char* known_spelled = "class, name, impedance_model, flow_resistance, porosity_rate, layer_depth "
                                      "and roughness";

std::string trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t\r\n\f\v");
    if (first == std::string::npos)
        return "";
    const auto last = text.find_last_not_of(" \t\r\n\f\v");
    return text.substr(first, last - first + 1);
}

// Split on every comma, keeping empty fields, and trim each.
std::vector<std::string> split_fields(const std::string& line)
{
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true)
    {
        const auto comma = line.find(',', begin);
        fields.push_back(trim(line.substr(begin, comma == std::string::npos ? std::string::npos : comma - begin)));
        if (comma == std::string::npos)
            return fields;
        begin = comma + 1;
    }
}

// The whole field as a finite number, or a refusal naming the line and column.
double parse_number(const std::string& field, const std::string& path, std::size_t line, const std::string& column)
{
    double value = 0.0;
    std::size_t used = 0;
    try
    {
        value = std::stod(field, &used);
    }
    catch (const std::exception&)
    {
        used = 0;
    }
    if (used == 0 || used != field.size() || !std::isfinite(value))
        input_error(path, " line ", line, ": the ", column, " column holds \"", field,
                "\", which is not a finite number.");
    return value;
}

bool valid_name(const std::string& name)
{
    if (name.empty() || name.front() == '.')
        return false;
    return std::all_of(name.begin(), name.end(), [](unsigned char c) {
        return std::isalnum(c) || c == '_' || c == '.' || c == '-';
    });
}

std::optional<ImpedanceModel> impedance_model_named(const std::string& name)
{
    if (name == "delany_bazley") return ImpedanceModel::delany_bazley;
    if (name == "miki") return ImpedanceModel::miki;
    if (name == "variable_porosity") return ImpedanceModel::variable_porosity;
    if (name == "delany_bazley_layer") return ImpedanceModel::delany_bazley_layer;
    return std::nullopt;
}

} // namespace

GroundClassFile niceops::read_ground_classes(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
        input_error("Unable to open the ground classes file ", path, ".");
    GroundClassFile file;
    file.path = path;

    LabeledCsv header;
    std::map<std::string, std::size_t> index;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(in, line))
    {
        line_number++;
        detail::strip_carriage_return(line);
        const std::string content = trim(line);
        if (content.empty() || content.front() == '#')
            continue;
        if (content.find('"') != std::string::npos)
            input_error(path, " line ", line_number, " holds a quote. The ground classes file is plain "
                    "comma-separated text with no quoting.");
        const std::vector<std::string> fields = split_fields(content);
        if (header.header.empty())
        {
            header.header = fields;
            for (std::size_t i = 0; i < fields.size(); i++)
            {
                const std::string name = LabeledCsv::normalize(fields[i]);
                if (name == "turbulence" || name == "turbulencemodel" || name == "gammat")
                    input_error(path, " line ", line_number, ": column \"", fields[i], "\": turbulence is the air's, "
                            "not the ground's. Set it for the whole run with --turbulence, --turbulence_model or "
                            "--gamma_t.");
                if (std::find(known_columns().begin(), known_columns().end(), name) == known_columns().end())
                    input_error(path, " line ", line_number, ": unknown column \"", fields[i], "\". The columns are ",
                            known_spelled, ".");
                if (!index.emplace(name, i).second)
                    input_error(path, " line ", line_number, ": column \"", fields[i], "\" is given twice.");
                file.columns.insert(name);
            }
            if (index.count("class") == 0)
                input_error(path, " line ", line_number, ": the header has no class column; every row needs the "
                        "integer id an observer file's ground column names it by.");
            // Units: normalize drops them from the name, so each is checked.
            // The flow resistance may say kPa (s/m^2), the porosity rate 1/m,
            // and the two lengths are checked against --units by the caller.
            for (const auto& [column, value] : index)
            {
                const std::string token = header.unit_token(value);
                if (token.empty() || column == "layerdepth" || column == "roughness")
                    continue;
                if (column == "flowresistance" && token == "kpa")
                    continue;
                if (column == "porosityrate" && token == "1/m")
                    continue;
                const std::string what = column == "flowresistance" ? "it is in kPa s/m^2"
                        : column == "porosityrate"                  ? "it is in 1/m"
                                                                    : "it has no unit";
                input_error(path, ": column \"", header.header[value], "\" declares ", token, ", and ", what,
                        ". Nothing in the ground classes file is converted but its lengths, which are in --units.");
            }
            if (index.count("layerdepth") != 0)
                file.layer_depth_unit = header.unit_token(index.at("layerdepth"));
            if (index.count("roughness") != 0)
                file.roughness_unit = header.unit_token(index.at("roughness"));
            continue;
        }

        if (fields.size() != header.header.size())
            input_error(path, " line ", line_number, " has ", fields.size(), " fields, and the header has ",
                    header.header.size(), ".");
        const auto cell = [&](const std::string& column) -> std::string {
            const auto found = index.find(column);
            return found == index.end() ? std::string() : fields[found->second];
        };
        const auto number = [&](const std::string& column, const char* spelled) -> std::optional<double> {
            const std::string field = cell(column);
            if (field.empty())
                return std::nullopt;
            return parse_number(field, path, line_number, spelled);
        };

        GroundClassRow row;
        row.line = line_number;
        // The id: the digits of an integer from 0 to 2^31 - 1, the range a
        // land-cover code takes, and nothing else ("1.0", "-1", "1e3" are not).
        const std::string id = cell("class");
        constexpr std::int64_t largest_id = 2147483647;
        const bool digits = !id.empty() && id.size() <= 10
                && std::all_of(id.begin(), id.end(), [](unsigned char c) { return std::isdigit(c); });
        if (!digits || std::stoll(id) > largest_id)
            input_error(path, " line ", line_number, ": the class \"", id, "\" is not an integer from 0 to ",
                    largest_id, ".");
        row.id = std::stoll(id);
        row.name = cell("name").empty() ? id : cell("name");
        if (!valid_name(row.name))
            input_error(path, " line ", line_number, ": the name \"", row.name, "\" may hold only letters, digits, "
                    "_, . and -, and may not start with a dot.");
        if (!cell("impedancemodel").empty())
        {
            row.impedance_model = impedance_model_named(cell("impedancemodel"));
            if (!row.impedance_model)
                input_error(path, " line ", line_number, ": unknown impedance model \"", cell("impedancemodel"),
                        "\": use delany_bazley, miki, variable_porosity or delany_bazley_layer.");
        }
        row.flow_resistance = number("flowresistance", "flow_resistance");
        if (row.flow_resistance && !(*row.flow_resistance > 0.0))
            input_error(path, " line ", line_number, ": the flow resistance ", *row.flow_resistance,
                    " kPa s/m^2 is not positive.");
        row.porosity_rate = number("porosityrate", "porosity_rate");
        if (row.porosity_rate && !(*row.porosity_rate >= 0.0))
            input_error(path, " line ", line_number, ": the porosity rate ", *row.porosity_rate, " is negative.");
        row.layer_depth = number("layerdepth", "layer_depth");
        if (row.layer_depth && !(*row.layer_depth > 0.0))
            input_error(path, " line ", line_number, ": the layer depth ", *row.layer_depth, " is not positive.");
        row.roughness = number("roughness", "roughness");
        if (row.roughness && !(*row.roughness >= 0.0))
            input_error(path, " line ", line_number, ": the roughness ", *row.roughness, " is negative.");
        for (const auto& earlier : file.rows)
            if (earlier.id == row.id)
                input_error(path, ": the rows on lines ", earlier.line, " and ", line_number, " are both class ",
                        row.id, ".");
        file.rows.push_back(std::move(row));
    }
    if (header.header.empty())
        input_error("The ground classes file ", path, " has no header.");
    if (file.rows.empty())
        input_error("The ground classes file ", path, " has no classes.");
    return file;
}

std::vector<GroundClass> niceops::ground_classes_from(const GroundClassFile& file, const GroundModel& ground,
        LengthUnit unit)
{
    const std::set<std::string> spellings = unit.to_meters == 1.0
            ? std::set<std::string> {"m", "meter", "meters", "metre", "metres"}
            : std::set<std::string> {"ft", "feet", "foot"};
    for (const auto& [token, column] : {std::pair {file.layer_depth_unit, "layer_depth"},
                 std::pair {file.roughness_unit, "roughness"}})
        if (!token.empty() && spellings.count(token) == 0)
            input_error(file.path, ": the ", column, " column declares ", token, ", and this run reads lengths in ",
                    unit.name, " (--units). Nothing is converted: give the column in ", unit.name,
                    ", or change --units.");
    std::vector<GroundClass> classes;
    classes.reserve(file.rows.size());
    for (const auto& row : file.rows)
    {
        GroundClass surface;
        surface.id = row.id;
        surface.name = row.name;
        surface.line = row.line;
        surface.model = ground;
        if (row.impedance_model)
            surface.model.impedance_model = *row.impedance_model;
        if (row.flow_resistance)
            surface.model.flow_resistance = *row.flow_resistance;
        if (row.porosity_rate)
            surface.model.porosity_rate = *row.porosity_rate;
        if (row.layer_depth)
            surface.model.layer_depth = *row.layer_depth * unit.to_meters;
        if (row.roughness)
            surface.model.roughness = *row.roughness * unit.to_meters;
        classes.push_back(surface);
    }
    return classes;
}

std::vector<double> niceops::point_ground_gain_bounds(const std::vector<GroundModel>& models,
        const std::vector<std::vector<std::complex<double>>>& impedances, const std::vector<std::size_t>& classes,
        std::size_t first_band, std::size_t last_band)
{
    assert(!classes.empty());
    auto gains = ground_gain_bounds(models[classes[0]], impedances[classes[0]], first_band, last_band);
    for (std::size_t k = 1; k < classes.size(); k++)
    {
        const auto more = ground_gain_bounds(models[classes[k]], impedances[classes[k]], first_band, last_band);
        for (std::size_t band = 0; band < gains.size(); band++)
            gains[band] = std::max(gains[band], more[band]);
    }
    return gains;
}

const GroundClass* niceops::find_ground_class(const std::vector<GroundClass>& classes, std::int64_t id)
{
    const auto found = std::find_if(classes.begin(), classes.end(),
            [id](const GroundClass& surface) { return surface.id == id; });
    return found == classes.end() ? nullptr : &*found;
}
