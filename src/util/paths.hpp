#pragma once

// Create a path's parent directory if it does not exist.
//
// Four places needed exactly this -- the log file, the CSV writers, and both
// netCDF writers -- and each had its own copy.  It also records the path for
// the run manifest.

#include <filesystem>
#include <mutex>
#include <string>
#include <unordered_set>
#include <vector>

namespace niceops {

namespace detail {

// Every writer passes its file through ensure_parent_directory, so that is
// where the run manifest (--manifest) learns what a run wrote.  Off unless
// the manifest asks, so a normal run records nothing.
struct WrittenFiles
{
    std::mutex mutex;
    bool enabled = false;
    std::vector<std::string> paths; // in the order first written
    std::unordered_set<std::string> seen;
};

inline WrittenFiles& written_files()
{
    static WrittenFiles files;
    return files;
}

} // namespace detail

inline void ensure_parent_directory(const std::string& path)
{
    const auto parent = std::filesystem::path(path).parent_path();
    if (!parent.empty())
        std::filesystem::create_directories(parent);
    auto& written = detail::written_files();
    std::lock_guard<std::mutex> lock(written.mutex);
    if (written.enabled && written.seen.insert(path).second)
        written.paths.push_back(path);
}

} // namespace niceops
