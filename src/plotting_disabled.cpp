// Stand-in for plotting.cpp in builds configured with -DNICEOPS_PLOTTING=OFF.
//
// Dropping plotting drops Python, pybind11 and matplotlib with it, which is the
// whole point: the acoustics need none of them, and a build that only computes
// footprints should not need a Python development install to link.
//
// Every entry point is present so the rest of the program still compiles and
// links unchanged; each one does nothing.  main() refuses the plotting flags up
// front via available(), so in practice none of these is ever reached -- they
// exist so that a future caller which forgets to check degrades quietly rather
// than failing to link.

#include "plotting.hpp"

#include "util/logger.hpp"

namespace plt {

bool available()
{
    return false;
}

namespace {

void unavailable()
{
    niceops::warn("This build has no plotting support; rebuild with "
                  "-DNICEOPS_PLOTTING=ON to draw figures.");
}

} // namespace

void figure()
{
    unavailable();
}

void contourf(const std::vector<double>&, const std::vector<double>&,
        const std::vector<std::vector<double>>&, const std::vector<double>&, const std::string&)
{
    unavailable();
}

void tricontourf(const std::vector<double>&, const std::vector<double>&,
        const std::vector<double>&, const std::vector<double>&, const std::string&)
{
    unavailable();
}

void plot(const std::vector<double>&, const std::vector<double>&)
{
    unavailable();
}

void plot(const std::vector<double>&, const std::vector<double>&, const std::string&)
{
    unavailable();
}

void plot(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
        const std::vector<double>&)
{
    unavailable();
}

void xlabel(const std::string&, int)
{
    unavailable();
}

void ylabel(const std::string&, int)
{
    unavailable();
}

void axis(const std::string&)
{
    unavailable();
}

void colorbar()
{
    unavailable();
}

void tight_layout()
{
    unavailable();
}

void savefig(const std::string&)
{
    unavailable();
}

void show()
{
    unavailable();
}

void clear()
{
    unavailable();
}

} // namespace plt
