#pragma once

// The operations file (--operations): one row per operation of an average day,
// in place of a -t per trajectory.  Each row names a track and how many times
// it is flown, and may name its own database, vehicle, clock, lateral offset
// and air.  The reader is pure: it checks the file's syntax and that the files
// it names exist, and converts nothing (lengths stay in the file's --units,
// paths are resolved against the operations file's own folder).

#include <array>
#include <cstddef>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace niceops {

struct Operation
{
    // 1-based line of the row in the file, for the messages.
    std::size_t line = 0;
    // The row's name: its output stem and its -m case name.  The track's stem
    // when the file gives none.
    std::string name;
    // Paths as resolved against the operations file's folder; empty when the
    // row (or the file) gives none.
    std::string track;
    std::string database;
    std::string vehicle;
    std::string output;
    std::string atmosphere;
    // --wind's text, "S,D" or "S,D,H", joined from the wind columns; empty when
    // the row gives no wind.
    std::string wind_text;
    // Operations per average day.  Unset means 1, or the period counts below.
    std::optional<double> count;
    // Day, evening and night counts, when the file has those columns.
    std::optional<std::array<double, 3>> periods;
    // Local seconds after midnight at the track's start (local_time).
    std::optional<double> local_time;
    std::optional<double> humidity;    // percent
    std::optional<double> temperature; // kelvin
    std::optional<double> pressure;    // kPa
    // Lateral offset in --units, + to the left of travel.
    std::optional<double> offset;
};

struct OperationsFile
{
    std::string path;
    std::vector<Operation> rows;
    // The normalized names of the columns the header carries (as
    // LabeledCsv::normalize gives them: "local_time" is "localtime").
    std::set<std::string> columns;
    // The length unit the offset and wind_height headers declare, lower-cased,
    // or empty: checked against --units by the caller.
    std::string offset_unit;
    std::string wind_height_unit;

    bool has(const std::string& column) const
    {
        return columns.count(column) != 0;
    }
};

// Reads and checks `path`.  Refuses, naming the line: an unknown, duplicate or
// missing track column, a row whose field count is not the header's, a quote, a
// count that is not finite and positive, period counts mixed with count or
// local_time, a malformed local_time, a name outside [A-Za-z0-9_.-] or shared
// with another row (ignoring case), and a track, database, vehicle or
// atmosphere file that does not exist.
OperationsFile read_operations(const std::string& path);

// HH:MM or HH:MM:SS, 00:00 to 23:59:59, as seconds after midnight; refuses
// anything else.
double parse_local_clock(const std::string& text);

// What a row's operations add up to: its count, its period counts' sum, or 1.
double operation_count(const Operation& row);

} // namespace niceops
