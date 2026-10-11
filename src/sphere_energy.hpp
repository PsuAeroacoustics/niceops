#pragma once

// --provenance_spheres: each observer's SEL energy shared among the database's
// spheres, the whole matrix where --provenance_file writes only the top
// sphere, so that energy can be summed over any grouping of spheres (by speed,
// by flight day, measured or synthetic) after the run.

#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace niceops {

struct SphereEnergyFile
{
    std::vector<double> x, y;     // each observer's position, local meters
    std::vector<double> exposure; // each observer's energy shared out (the SEL's, linear)
    // Each sphere's stored condition: advance ratio, flight path angle (deg),
    // thrust coefficient, in the database's group order.
    std::vector<std::array<double, 3>> conditions;
    std::string database;         // the database's path
    bool grid = false;            // observers in -o's grid order (x fastest), else -b's
};

// Writes a netCDF file: per observer its position and exposure_db, per sphere
// its condition, and energy_share(observer, sphere), the share of the
// observer's energy each sphere supplied (summing to one; NaN for an observer
// that heard nothing).  `energy[i]` is observer i's sparse list of (sphere,
// energy).  Throws (input_error) if the file cannot be written.
void write_sphere_energy(const std::string& path, const SphereEnergyFile& meta,
        const std::vector<std::vector<std::pair<std::uint32_t, double>>>& energy);

} // namespace niceops
