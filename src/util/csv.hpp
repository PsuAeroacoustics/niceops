#pragma once

// Minimal numeric CSV reader replacing D's std.csv csvReader!(double).
// Inputs in this project are pure numeric (trajectory states, observer
// lists), so no quoting/escaping support is needed.  Lines that contain a
// non-numeric field (e.g. file headers) are skipped.

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace niceops {

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
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        if (first && skip_first_record)
        {
            first = false;
            continue;
        }
        first = false;
        std::vector<double> row;
        std::stringstream ss(line);
        std::string field;
        bool numeric = true;
        while (std::getline(ss, field, ','))
        {
            try
            {
                size_t pos = 0;
                double value = std::stod(field, &pos);
                if (pos != field.size())
                    throw std::invalid_argument("trailing junk");
                row.push_back(value);
            }
            catch (const std::exception&)
            {
                numeric = false;
                break;
            }
        }
        if (numeric && !row.empty())
            rows.push_back(std::move(row));
    }
    return rows;
}

} // namespace niceops
