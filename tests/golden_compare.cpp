// Golden output comparison tool.  Compares two numeric CSV files (the C++
// port's output against a golden file captured from the D build) element by
// element with a relative tolerance.  Non-numeric header lines are skipped,
// so files may differ in their text headers.
//
// Usage: golden_compare <candidate.csv> <golden.csv>
// Exit code 0 if they agree within tolerance, 1 otherwise.

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Parse a CSV of numbers, skipping blank and non-numeric lines (headers).
std::vector<double> parse_csv(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file)
    {
        throw std::runtime_error("Could not open " + filename);
    }
    std::vector<double> values;
    std::string line;
    while (std::getline(file, line))
    {
        // Strip a trailing carriage return (CRLF files)
        if (!line.empty() && line.back() == '\r')
        {
            line.pop_back();
        }
        std::istringstream tokens(line);
        std::string token;
        bool numeric_line = true;
        while (std::getline(tokens, token, ','))
        {
            try
            {
                std::size_t consumed = 0;
                auto value = std::stod(token, &consumed);
                if (consumed != token.size())
                {
                    throw std::invalid_argument("trailing garbage");
                }
                values.push_back(value);
            }
            catch (const std::exception&)
            {
                numeric_line = false;
                break;
            }
        }
        if (!numeric_line)
        {
            continue; // header or comment line
        }
    }
    return values;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3 && argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <candidate.csv> <golden.csv> [absolute tolerance]\n";
        return 1;
    }
    // An absolute tolerance, for comparing two runs that are meant to agree
    // to within an approximation rather than to the bit.
    const double absolute = argc == 4 ? std::atof(argv[3]) : 0.0;

    std::vector<double> candidate, golden;
    try
    {
        candidate = parse_csv(argv[1]);
        golden = parse_csv(argv[2]);
    }
    catch (const std::exception& e)
    {
        std::cerr << "golden_compare: " << e.what() << "\n";
        return 1;
    }

    if (candidate.size() != golden.size())
    {
        std::cerr << "golden_compare: size mismatch, candidate has "
                  << candidate.size() << " values, golden has "
                  << golden.size() << " values.\n";
        return 1;
    }

    const double tolerance = 1e-9;
    for (std::size_t i = 0; i < golden.size(); i++)
    {
        auto scale = std::max(1.0, std::max(std::abs(candidate[i]), std::abs(golden[i])));
        if (std::abs(candidate[i] - golden[i]) > std::max(tolerance * scale, absolute))
        {
            std::cerr << "golden_compare: mismatch at value " << i << ": candidate "
                      << candidate[i] << " vs golden " << golden[i] << "\n";
            return 1;
        }
    }

    std::cout << "golden_compare: " << golden.size() << " values agree within "
              << (absolute > 0.0 ? absolute : tolerance) << (absolute > 0.0 ? " absolute" : " relative")
              << " tolerance.\n";
    return 0;
}
