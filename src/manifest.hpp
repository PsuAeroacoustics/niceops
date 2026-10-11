#pragma once

// The run manifest (--manifest): one JSON file saying what produced a run's
// outputs, so a result can be traced to the binary, the options and the exact
// input files behind it.  A harness that ran a stale default binary for two
// days, fitting corrections to the wrong physics, had nothing that would have
// shown it.

#include <future>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace niceops {

// Whether two paths name one file: the same file when both exist, else the same
// path once made absolute and normalized.
bool same_file(const std::string& a, const std::string& b);

// Every root attribute and scalar variable of an open netCDF file, text as it
// is and numbers as text: how a database was built, for the manifest.  The
// loader reads them while it has the file open (Database::root_attributes).
std::map<std::string, std::string> netcdf_root_attributes(int ncid);

// What the build recorded about itself (src/build_info.cpp, from a header the
// build regenerates whenever the git state changes).
struct BuildInfo
{
    std::string version;   // `git describe --always --dirty --tags`, or "unknown"
    std::string compiler;  // e.g. "AppleClang 17.0.0"
    std::string build_type;
    bool asserts = false;
    bool plotting = false;
};
const BuildInfo& build_info();

class RunManifest
{
public:
    void set_command_line(int argc, const char* const* argv);
    // The options given, on the command line or in a --config file, as CLI11
    // writes them back (config_to_str).  Defaults are the build's, which the
    // version names.
    void set_options_given(std::string text)
    {
        options_given_ = std::move(text);
    }
    // The log file is appended to across runs, so it is named, not hashed.
    void set_log_file(std::string path)
    {
        log_file_ = std::move(path);
    }
    // A file the run read.  The same path under the same role is listed once,
    // and a file that does not exist when the manifest is written is listed
    // without a digest.  What the run wrote needs no listing: every writer's
    // file passes through ensure_parent_directory (util/paths.hpp), which
    // records it from the moment this is constructed.
    RunManifest();
    void add_input(const std::string& role, const std::string& path);
    // A database: listed as an input, with its root attributes and scalar
    // variables (how it was built: speed_reference, the build atmosphere,
    // hover_source, database_version, fixed_load_factor, ...), which the run
    // passes on once it has loaded it.  One that never loaded is read at write.
    void add_database(const std::string& path);
    void set_database_attributes(const std::string& path, std::map<std::string, std::string> attributes);

    // Starts hashing the inputs, one after another on one background thread,
    // so the digests are ready by the time the run ends.  Called once every
    // input is added.
    void start_hashing();

    // Refuses a manifest path that is one of the inputs, which the write
    // would truncate.  Called before the run, so it costs nothing; the
    // inputs do not change after it.
    void refuse_overwriting_an_input(const std::string& path) const;

    // Hashes every listed file and writes the manifest, the outputs the first
    // max_listed_outputs in the order written (a run writing per-observer
    // histories can write thousands; the rest are counted).  Throws
    // (input_error) if the file cannot be written.
    static constexpr std::size_t max_listed_outputs = 1000;
    void write(const std::string& path, int exit_status, double seconds) const;

private:
    struct Entry
    {
        std::string role;
        std::string path;
    };

    std::vector<std::string> argv_;
    std::string options_given_;
    std::string log_file_;
    std::vector<Entry> inputs_;
    std::map<std::string, std::map<std::string, std::string>> databases_;
    // The background hashing's result: each input's digest, by path.
    using Digests = std::map<std::string, std::optional<std::string>>;
    std::shared_future<Digests> digests_;
};

} // namespace niceops
