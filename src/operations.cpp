#include "operations.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "util/csv.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

namespace {

// Every column the file may carry, by its normalized name, in the order the
// refusal lists them.
const std::vector<std::string>& known_columns()
{
    static const std::vector<std::string> names = {"name", "track", "count", "day", "evening", "night",
            "localtime", "database", "vehicle", "output", "atmosphere", "windspeed", "winddirection",
            "windheight", "humidity", "temperature", "pressure", "offset"};
    return names;
}

std::string trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t\r\n\f\v");
    if (first == std::string::npos)
        return "";
    const auto last = text.find_last_not_of(" \t\r\n\f\v");
    return text.substr(first, last - first + 1);
}

// Split on every comma, keeping empty fields (a trailing comma is one more,
// empty, field), and trim each.
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

std::string lower(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return text;
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

std::string resolve(const std::filesystem::path& folder, const std::string& given)
{
    const std::filesystem::path path(given);
    if (path.is_absolute())
        return path.lexically_normal().string();
    return (folder / path).lexically_normal().string();
}

bool valid_name(const std::string& name)
{
    if (name.empty() || name.front() == '.')
        return false;
    return std::all_of(name.begin(), name.end(), [](unsigned char c) {
        return std::isalnum(c) || c == '_' || c == '.' || c == '-';
    });
}

// HH:MM or HH:MM:SS as seconds after midnight, or nullopt.  One or two digits
// of hours, then two of minutes and optionally two of seconds: anything else
// ("7pm", "7:30 am", "24:00") is not a clock time.
std::optional<double> clock_seconds(const std::string& text)
{
    const auto digits = [&](std::size_t begin, std::size_t end) {
        if (end <= begin || end > text.size())
            return false;
        return std::all_of(text.begin() + static_cast<std::ptrdiff_t>(begin),
                text.begin() + static_cast<std::ptrdiff_t>(end), [](unsigned char c) { return std::isdigit(c); });
    };
    const auto colon = text.find(':');
    if (colon == std::string::npos || colon == 0 || colon > 2 || !digits(0, colon))
        return std::nullopt;
    const auto second_colon = text.find(':', colon + 1);
    const std::size_t minutes_end = second_colon == std::string::npos ? text.size() : second_colon;
    if (minutes_end != colon + 3 || !digits(colon + 1, minutes_end))
        return std::nullopt;
    int seconds = 0;
    if (second_colon != std::string::npos)
    {
        if (text.size() != second_colon + 3 || !digits(second_colon + 1, text.size()))
            return std::nullopt;
        seconds = std::stoi(text.substr(second_colon + 1));
    }
    const int hours = std::stoi(text.substr(0, colon));
    const int minutes = std::stoi(text.substr(colon + 1, 2));
    if (hours > 23 || minutes > 59 || seconds > 59)
        return std::nullopt;
    return hours * 3600.0 + minutes * 60.0 + seconds;
}

} // namespace

double niceops::parse_local_clock(const std::string& text)
{
    const auto seconds = clock_seconds(text);
    if (!seconds)
        input_error("local_time \"", text, "\" is not a clock time: give HH:MM or HH:MM:SS, 00:00 to 23:59:59.");
    return *seconds;
}

double niceops::operation_count(const Operation& row)
{
    if (row.periods)
        return (*row.periods)[0] + (*row.periods)[1] + (*row.periods)[2];
    return row.count.value_or(1.0);
}

OperationsFile niceops::read_operations(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
        input_error("Unable to open the operations file ", path, ".");
    OperationsFile file;
    file.path = path;
    const std::filesystem::path folder = std::filesystem::path(path).parent_path();

    // The header's columns: normalized name to index, and the raw header for
    // the declared units.
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
            input_error(path, " line ", line_number, " holds a quote. The operations file is plain comma-separated "
                    "text with no quoting: a field cannot contain a comma or a quote.");
        std::vector<std::string> fields = split_fields(content);
        if (header.header.empty())
        {
            header.header = fields;
            for (std::size_t i = 0; i < fields.size(); i++)
            {
                const std::string name = LabeledCsv::normalize(fields[i]);
                if (name == "startutc")
                    input_error(path, " line ", line_number, ": an operations file gives each row's clock as "
                            "local_time (HH:MM, local time by --utc_offset), not start_utc.");
                if (std::find(known_columns().begin(), known_columns().end(), name) == known_columns().end())
                {
                    std::string known;
                    for (const auto& k : {"name", "track", "count", "day", "evening", "night", "local_time",
                                 "database", "vehicle", "output", "atmosphere", "wind_speed", "wind_direction",
                                 "wind_height", "humidity", "temperature", "pressure", "offset"})
                        known += std::string(known.empty() ? "" : ", ") + k;
                    input_error(path, " line ", line_number, ": unknown column \"", fields[i], "\". The columns are ",
                            known, ".");
                }
                if (!index.emplace(name, i).second)
                    input_error(path, " line ", line_number, ": column \"", fields[i], "\" is given twice.");
                file.columns.insert(name);
            }
            if (!file.has("track"))
                input_error(path, " line ", line_number, ": the header has no track column; every operation "
                        "needs one.");
            const bool periods = file.has("day") || file.has("evening") || file.has("night");
            if (periods && file.has("count"))
                input_error(path, ": give each row's operations as count or as day, evening and night counts, "
                        "not both.");
            if (periods && file.has("localtime"))
                input_error(path, ": local_time dates an operation, and day, evening and night counts take the "
                        "clock's place; give one or the other.");
            if (file.has("windspeed") != file.has("winddirection"))
                input_error(path, ": wind_speed and wind_direction go together.");
            if (file.has("windheight") && !file.has("windspeed"))
                input_error(path, ": wind_height is the anemometer's height, and needs wind_speed and "
                        "wind_direction.");
            // Columns with a fixed unit: a header naming another unit means a
            // file in it, which nothing here converts.
            for (const auto& [column, units, expected] : {std::tuple{"humidity", std::set<std::string>{"%", "percent"}, "%"},
                         std::tuple{"temperature", std::set<std::string>{"k", "kelvin"}, "K"},
                         std::tuple{"pressure", std::set<std::string>{"kpa"}, "kPa"}})
            {
                const auto found = index.find(column);
                if (found == index.end())
                    continue;
                const std::string token = header.unit_token(found->second);
                if (!token.empty() && units.count(token) == 0)
                    input_error(path, ": column \"", header.header[found->second], "\" declares ", token, "; the ",
                            column, " column is in ", expected, ".");
            }
            // Columns with no unit to declare: a label on one would be dropped
            // by normalize, so `count (per year)` would read annual counts as
            // daily ones. The wind direction may say it is in degrees, and the
            // wind speed names its unit in each cell.
            for (const auto& [column, value] : index)
            {
                const std::string token = header.unit_token(value);
                if (token.empty() || column == "humidity" || column == "temperature" || column == "pressure"
                        || column == "offset" || column == "windheight")
                    continue;
                if (column == "winddirection" && (token == "deg" || token == "degree" || token == "degrees"))
                    continue;
                const std::string what = column == "windspeed" ? "the speed names its unit in each cell (8kt)"
                        : column == "winddirection"            ? "it is a compass bearing in degrees"
                        : column == "localtime"                ? "it is local time, by --utc_offset"
                        : column == "count" || column == "day" || column == "evening" || column == "night"
                        ? "it counts operations per average day"
                        : "it has no unit";
                input_error(path, ": column \"", header.header[value], "\" declares ", token, ", and ", what,
                        ". Nothing in the operations file is converted.");
            }
            if (file.has("offset"))
                file.offset_unit = header.unit_token(index.at("offset"));
            if (file.has("windheight"))
                file.wind_height_unit = header.unit_token(index.at("windheight"));
            continue;
        }

        if (fields.size() != header.header.size())
            input_error(path, " line ", line_number, " has ", fields.size(), " fields, and the header has ",
                    header.header.size(), ".");
        const auto cell = [&](const std::string& column) -> std::string {
            const auto found = index.find(column);
            return found == index.end() ? std::string() : fields[found->second];
        };
        const auto number = [&](const std::string& column) -> std::optional<double> {
            const std::string field = cell(column);
            if (field.empty())
                return std::nullopt;
            return parse_number(field, path, line_number, column);
        };

        Operation row;
        row.line = line_number;
        if (cell("track").empty())
            input_error(path, " line ", line_number, " names no track.");
        row.track = resolve(folder, cell("track"));
        for (auto [column, target] : {std::pair{"database", &row.database}, {"vehicle", &row.vehicle},
                     {"output", &row.output}, {"atmosphere", &row.atmosphere}})
            if (!cell(column).empty())
                *target = resolve(folder, cell(column));
        row.name = cell("name").empty() ? std::filesystem::path(row.track).stem().string() : cell("name");
        if (!valid_name(row.name))
            input_error(path, " line ", line_number, ": the name \"", row.name, "\" may hold only letters, digits, "
                    "_, . and -, and may not start with a dot.",
                    cell("name").empty() ? " It is the track's file name; add a name column." : "");

        row.count = number("count");
        if (row.count && !(*row.count > 0.0))
            input_error(path, " line ", line_number, ": count ", *row.count, " is not positive. Comment out a row "
                    "(#) to leave it out.");
        if (file.has("day") || file.has("evening") || file.has("night"))
        {
            std::array<double, 3> periods {};
            const char* names[] = {"day", "evening", "night"};
            for (std::size_t p = 0; p < 3; p++)
            {
                periods[p] = number(names[p]).value_or(0.0);
                if (periods[p] < 0.0)
                    input_error(path, " line ", line_number, ": the ", names[p], " count ", periods[p],
                            " is negative.");
            }
            if (!(periods[0] + periods[1] + periods[2] > 0.0))
                input_error(path, " line ", line_number, ": the day, evening and night counts add up to no "
                        "operations. Comment out a row (#) to leave it out.");
            row.periods = periods;
        }
        if (!cell("localtime").empty())
        {
            row.local_time = clock_seconds(cell("localtime"));
            if (!row.local_time)
                input_error(path, " line ", line_number, ": local_time \"", cell("localtime"), "\" is not a clock "
                        "time: give HH:MM or HH:MM:SS, 00:00 to 23:59:59.");
        }
        row.humidity = number("humidity");
        row.temperature = number("temperature");
        row.pressure = number("pressure");
        row.offset = number("offset");

        const std::string speed = cell("windspeed");
        const std::string direction = cell("winddirection");
        const std::string height = cell("windheight");
        if (speed.empty() != direction.empty() || (!height.empty() && speed.empty()))
            input_error(path, " line ", line_number, ": give wind_speed and wind_direction together (and "
                    "wind_height only with them), or leave all of them empty.");
        if (!speed.empty())
            row.wind_text = speed + "," + direction + (height.empty() ? "" : "," + height);

        // The files are checked here, before anything loads, with the path as
        // resolved: a relative path is against the operations file's folder,
        // which is not where a user typing -t would look.
        for (auto [what, file_path] : {std::pair{"track", &row.track}, {"database", &row.database},
                     {"vehicle", &row.vehicle}, {"atmosphere", &row.atmosphere}})
        {
            if (file_path->empty())
                continue;
            std::error_code ec;
            if (!std::filesystem::is_regular_file(*file_path, ec))
                input_error(path, " line ", line_number, " (", row.name, "): the ", what, " file ", *file_path,
                        " does not exist. Relative paths are taken from the operations file's folder.");
            info("Operation ", row.name, " (", path, " line ", line_number, "): ", what, " ", *file_path);
        }
        file.rows.push_back(std::move(row));
    }
    if (header.header.empty())
        input_error("The operations file ", path, " has no header.");
    if (file.rows.empty())
        input_error("The operations file ", path, " has no operations.");

    // Names are file names, and two rows writing one file is never meant.
    // Compared ignoring case, as some file systems do.
    std::map<std::string, const Operation*> seen;
    for (const auto& row : file.rows)
    {
        const auto [found, inserted] = seen.emplace(lower(row.name), &row);
        if (!inserted)
            input_error(path, ": the rows on lines ", found->second->line, " and ", row.line, " share the name ",
                    row.name, file.has("name") ? "; give them distinct names." : "; add a name column.");
    }
    return file;
}
