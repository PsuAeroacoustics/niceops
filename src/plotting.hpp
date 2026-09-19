#pragma once

// Replacement for the vendored matplotlib-d binding.  Instead of shelling out
// to a D-to-Python bridge, this embeds a Python interpreter (pybind11) and
// calls matplotlib directly.  The interpreter is created lazily on the first
// plotting call, so non-plotting runs never pay the import cost.

#include <string>
#include <vector>

namespace plt {

// Whether this build has plotting compiled in.  Building with
// -DNICEOPS_PLOTTING=OFF drops Python, pybind11 and matplotlib from the
// dependency list entirely, and every call below becomes a no-op that says so.
bool available();

void figure();

void contourf(const std::vector<double>& x, const std::vector<double>& y,
        const std::vector<std::vector<double>>& z, const std::vector<double>& levels,
        const std::string& cmap);

void tricontourf(const std::vector<double>& x, const std::vector<double>& y,
        const std::vector<double>& z, const std::vector<double>& levels,
        const std::string& cmap);

void plot(const std::vector<double>& x, const std::vector<double>& y);

void plot(const std::vector<double>& x, const std::vector<double>& y,
        const std::string& format);

void plot(const std::vector<double>& x1, const std::vector<double>& y1,
        const std::vector<double>& x2, const std::vector<double>& y2);

void xlabel(const std::string& label, int fontsize = 0);

void ylabel(const std::string& label, int fontsize = 0);

void axis(const std::string& what);

void colorbar();

void tight_layout();

void savefig(const std::string& filename);

void show();

// Close all figures (matplotlib-d's clear())
void clear();

} // namespace plt
