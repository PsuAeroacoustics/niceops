#pragma once

// Minimal numeric CSV reader replacing D's std.csv csvReader!(double).
// Inputs in this project are pure numeric (trajectory states, observer
// lists), so no quoting/escaping support is needed.  Lines that contain a
// non-numeric field (e.g. file headers) are skipped.

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace niceops {

namespace detail {

// Strip a trailing carriage return so CRLF files parse.  Without this every
// row's last field carries a '\r', fails to convert, and the whole line is
// silently dropped as non-numeric -- a Windows-exported CSV would read as zero
// records rather than failing loudly.
inline void strip_carriage_return(std::string& line)
{
    if (!line.empty() && line.back() == '\r')
        line.pop_back();
}

// Split one line on commas into doubles, returning false if any field is not a
// number.  Both readers below parse rows identically, so they share this.
//
// Surrounding whitespace is not "junk": std::stod already skips what precedes a
// number, and rejecting only what follows it made the two ends behave
// differently for no reason.  Every shipped AC track pads on the left and so
// parsed anyway, while arrays/coyle_approach.csv pads on the right and had
// every one of its rows dropped -- a scattered-observer run that read zero
// observers and still exited 0.
inline bool parse_numeric_row(const std::string& line, std::vector<double>& row)
{
    row.clear();
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ','))
    {
        const auto last = field.find_last_not_of(" \t\r\n\f\v");
        if (last == std::string::npos)
            return false;
        field.erase(last + 1);
        try
        {
            std::size_t pos = 0;
            const double value = std::stod(field, &pos);
            if (pos != field.size())
                return false;
            row.push_back(value);
        }
        catch (const std::exception&)
        {
            return false;
        }
    }
    return !row.empty();
}

} // namespace detail

// skip_first_record mirrors D's csvReader!(double)(contents, null): passing a
// null header makes std.csv consume the first record as a header row and
// discard it, unconditionally -- even when it holds numbers.  Callers that
// port such a csvReader must pass true, or they will read one extra observer
// that the D original silently dropped.
inline std::vector<std::vector<double>> read_numeric_csv(const std::string& path,
        bool skip_first_record = false)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Unable to open CSV file: " + path);
    std::vector<std::vector<double>> rows;
    std::string line;
    bool first = true;
    std::vector<double> row;
    while (std::getline(in, line))
    {
        detail::strip_carriage_return(line);
        if (line.empty())
            continue;
        if (first && skip_first_record)
        {
            first = false;
            continue;
        }
        first = false;
        if (detail::parse_numeric_row(line, row))
            rows.push_back(row);
    }
    return rows;
}

// A CSV whose first line names its columns, rather than fixing their order.
//
// The 27-column NASA/FAA AC state format is positional and stays that way, but
// the simpler trajectory formats are meant to be produced by whatever a mission
// planner already has, so they are read by column name instead -- any order,
// and only the columns that file actually carries.
struct LabelledCsv
{
    static constexpr std::size_t npos = static_cast<std::size_t>(-1);

    std::vector<std::string> header;
    std::vector<std::vector<double>> rows;

    // Index of the first column whose name matches any of `names`, or npos.
    // Both sides are normalised first (see normalise below), so "Latitude",
    // " lat (deg) " and "LAT" all find the same column.
    std::size_t column(std::initializer_list<const char*> names) const
    {
        for (std::size_t i = 0; i < header.size(); i++)
        {
            const std::string field = normalise(header[i]);
            for (const char* name : names)
            {
                if (field == normalise(name))
                    return i;
            }
        }
        return npos;
    }

    // Lower-case, drop any parenthesised or bracketed unit, and keep only
    // letters and digits.  That makes "Alt (ft)", "alt_ft" and "ALT" one name,
    // without the alias lists having to spell every punctuation variant.
    static std::string normalise(const std::string& field)
    {
        std::string out;
        for (char c : field)
        {
            if (c == '(' || c == '[')
                break;
            const unsigned char u = static_cast<unsigned char>(c);
            if (std::isalnum(u))
                out.push_back(static_cast<char>(std::tolower(u)));
        }
        return out;
    }
};

// Read `path` as a labelled CSV: the first non-empty line names the columns and
// every fully numeric line after it is a row.  Non-numeric lines are dropped as
// read_numeric_csv drops them, so a file with a units line under its header
// still reads.
inline LabelledCsv read_labelled_csv(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Unable to open CSV file: " + path);
    LabelledCsv csv;
    std::string line;
    bool first = true;
    std::vector<double> row;
    while (std::getline(in, line))
    {
        detail::strip_carriage_return(line);
        if (line.empty())
            continue;
        if (first)
        {
            first = false;
            std::stringstream ss(line);
            std::string field;
            while (std::getline(ss, field, ','))
                csv.header.push_back(field);
            continue;
        }
        if (detail::parse_numeric_row(line, row))
            csv.rows.push_back(row);
    }
    return csv;
}

} // namespace niceops
