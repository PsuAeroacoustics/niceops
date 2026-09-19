#pragma once

// Create a path's parent directory if it does not exist.
//
// Four places needed exactly this -- the log file, the CSV writers, and both
// netCDF writers -- and each had its own copy.

#include <filesystem>
#include <string>

namespace niceops {

inline void ensure_parent_directory(const std::string& path)
{
    const auto parent = std::filesystem::path(path).parent_path();
    if (!parent.empty())
        std::filesystem::create_directories(parent);
}

} // namespace niceops
