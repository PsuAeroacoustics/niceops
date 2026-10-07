// Whether two propagation tables hold the same axes and the same values, bit
// for bit (NaN matching NaN).  The files themselves differ in their provenance
// (the command line is recorded), so they cannot be compared as bytes.
//
//   check_tables_identical a.nc b.nc
//
// Exits 0 when they agree, 1 naming the first difference otherwise.

#include <cstddef>
#include <cstring>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "propagation_table.hpp"

using namespace niceops; // NOLINT

namespace {

bool same_bits(double a, double b)
{
    return std::memcmp(&a, &b, sizeof a) == 0;
}

bool same_axis(const char* name, const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size())
    {
        std::cerr << name << ": " << a.size() << " values against " << b.size() << "\n";
        return false;
    }
    for (std::size_t k = 0; k < a.size(); k++)
        if (!same_bits(a[k], b[k]))
        {
            std::cerr << name << "[" << k << "]: " << a[k] << " against " << b[k] << "\n";
            return false;
        }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr << "usage: check_tables_identical a.nc b.nc\n";
        return 2;
    }
    try
    {
        const auto a = PropagationTable::load(argv[1]);
        const auto b = PropagationTable::load(argv[2]);
        const auto &x = a.axes(), &y = b.axes();
        if (!(same_axis("bearing", x.bearings, y.bearings) && same_axis("source_height", x.source_heights, y.source_heights)
                    && same_axis("receiver_height", x.receiver_heights, y.receiver_heights)
                    && same_axis("range", x.ranges, y.ranges) && same_axis("band", x.bands, y.bands)))
            return 1;
        if (a.method() != b.method() || a.includes_absorption() != b.includes_absorption())
        {
            std::cerr << "method or includes_absorption differs\n";
            return 1;
        }
        std::size_t compared = 0;
        for (std::size_t i = 0; i < x.bearings.size(); i++)
            for (std::size_t s = 0; s < x.source_heights.size(); s++)
                for (std::size_t h = 0; h < x.receiver_heights.size(); h++)
                    for (std::size_t r = 0; r < x.ranges.size(); r++)
                        for (std::size_t f = 0; f < x.bands.size(); f++, compared++)
                            if (!same_bits(a.at(i, s, h, r, f), b.at(i, s, h, r, f)))
                            {
                                std::cerr << "value differs at bearing " << i << ", source " << s << ", receiver "
                                          << h << ", range " << r << ", band " << f << ": " << a.at(i, s, h, r, f)
                                          << " against " << b.at(i, s, h, r, f) << "\n";
                                return 1;
                            }
        std::cout << "identical: " << compared << " values\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "check_tables_identical: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
