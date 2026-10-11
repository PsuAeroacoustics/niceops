#include "manifest.hpp"

#include <algorithm>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <system_error>

#include <netcdf.h>
#include <nlohmann/json.hpp>

#include "util/logger.hpp"
#include "util/paths.hpp"
#include "util/sha256.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

RunManifest::RunManifest()
{
    auto& written = detail::written_files();
    std::lock_guard<std::mutex> lock(written.mutex);
    written.enabled = true;
}

void RunManifest::set_command_line(int argc, const char* const* argv)
{
    argv_.assign(argv, argv + argc);
}

void RunManifest::add_input(const std::string& role, const std::string& path)
{
    if (path.empty())
        return;
    const bool listed = std::any_of(inputs_.begin(), inputs_.end(),
            [&](const Entry& e) { return e.role == role && e.path == path; });
    if (!listed)
        inputs_.push_back({role, path});
}

std::map<std::string, std::string> niceops::netcdf_root_attributes(int file)
{
    std::map<std::string, std::string> attributes;
    int count = 0;
    if (nc_inq_natts(file, &count) == NC_NOERR)
    {
        for (int k = 0; k < count; k++)
        {
            char name[NC_MAX_NAME + 1] = {};
            nc_type type = NC_NAT;
            std::size_t length = 0;
            if (nc_inq_attname(file, NC_GLOBAL, k, name) != NC_NOERR
                    || nc_inq_att(file, NC_GLOBAL, name, &type, &length) != NC_NOERR)
                continue;
            if (type == NC_CHAR)
            {
                std::string text(length, '\0');
                if (nc_get_att_text(file, NC_GLOBAL, name, text.data()) == NC_NOERR)
                    attributes[name] = text.substr(0, text.find('\0'));
            }
            else if (length >= 1 && type != NC_STRING)
            {
                std::vector<double> values(length);
                if (nc_get_att_double(file, NC_GLOBAL, name, values.data()) == NC_NOERR)
                {
                    std::string text;
                    for (std::size_t v = 0; v < length; v++)
                        text += (v ? " " : "") + nlohmann::json(values[v]).dump();
                    attributes[name] = text;
                }
            }
        }
    }
    // Root scalar variables: database_version, fixed_load_factor, the rotor.
    int variables = 0;
    if (nc_inq_nvars(file, &variables) == NC_NOERR)
    {
        for (int v = 0; v < variables; v++)
        {
            char name[NC_MAX_NAME + 1] = {};
            int dims = -1;
            nc_type type = NC_NAT;
            double value = 0.0;
            if (nc_inq_var(file, v, name, &type, &dims, nullptr, nullptr) == NC_NOERR && dims == 0
                    && type != NC_CHAR && type != NC_STRING && nc_get_var_double(file, v, &value) == NC_NOERR
                    && attributes.count(name) == 0)
                attributes[name] = nlohmann::json(value).dump();
        }
    }
    return attributes;
}

void RunManifest::add_database(const std::string& path)
{
    add_input("database", path);
}

void RunManifest::set_database_attributes(const std::string& path, std::map<std::string, std::string> attributes)
{
    databases_.emplace(path, std::move(attributes));
}

void RunManifest::start_hashing()
{
    for (const auto& input : inputs_)
    {
        if (digests_.count(input.path) == 0)
            digests_[input.path] = std::async(std::launch::async, sha256_of_file, input.path).share();
    }
}

namespace {

// Each file's size and digest: from `digests` when it was hashed ahead,
// otherwise now.
template <typename Entries>
nlohmann::ordered_json describe_files(const Entries& files,
        const std::map<std::string, std::shared_future<std::optional<std::string>>>& digests)
{
    auto list = nlohmann::ordered_json::array();
    for (const auto& file : files)
    {
        nlohmann::ordered_json entry;
        if (!file.role.empty())
            entry["role"] = file.role;
        entry["path"] = file.path;
        std::error_code ec;
        const auto size = std::filesystem::file_size(file.path, ec);
        if (ec)
        {
            entry["exists"] = false;
        }
        else
        {
            entry["bytes"] = size;
            const auto ahead = digests.find(file.path);
            const auto digest = ahead != digests.end() ? ahead->second.get() : sha256_of_file(file.path);
            if (digest)
                entry["sha256"] = *digest;
            else
                entry["sha256"] = nullptr;
        }
        list.push_back(entry);
    }
    return list;
}

std::string utc_now()
{
    const auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm tm {};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char text[32];
    std::strftime(text, sizeof(text), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return text;
}

} // namespace

bool niceops::same_file(const std::string& a, const std::string& b)
{
    if (a.empty() || b.empty())
        return false;
    std::error_code ec;
    if (std::filesystem::exists(a, ec) && std::filesystem::exists(b, ec))
        return std::filesystem::equivalent(a, b, ec);
    const auto ca = std::filesystem::weakly_canonical(a, ec);
    const auto cb = std::filesystem::weakly_canonical(b, ec);
    return !ec && ca == cb;
}

void RunManifest::refuse_overwriting_an_input(const std::string& path) const
{
    for (const auto& input : inputs_)
        if (same_file(path, input.path))
            input_error("--manifest would overwrite the ", input.role, " ", input.path,
                    ". Give another file name.");
}

void RunManifest::write(const std::string& path, int exit_status, double seconds) const
{
    const auto& build = build_info();
    nlohmann::ordered_json manifest;
    manifest["niceops"] = {
            {"version", build.version},
            {"compiler", build.compiler},
            {"build_type", build.build_type},
            {"asserts", build.asserts},
            {"plotting", build.plotting},
    };
    manifest["written_utc"] = utc_now();
    std::error_code ec;
    const auto cwd = std::filesystem::current_path(ec);
    manifest["working_directory"] = ec ? std::string() : cwd.string();
    manifest["argv"] = argv_;
    manifest["options_given"] = options_given_;
    manifest["log_file"] = log_file_;

    manifest["inputs"] = describe_files(inputs_, digests_);
    // Every database, in path order: as the loader read it, or read now if it
    // never loaded (empty if it cannot be opened; the loader said why).
    std::map<std::string, std::map<std::string, std::string>> databases = databases_;
    for (const auto& input : inputs_)
    {
        if (input.role != "database" || databases.count(input.path) != 0)
            continue;
        auto& attributes = databases[input.path];
        int file = 0;
        if (nc_open(input.path.c_str(), NC_NOWRITE, &file) == NC_NOERR)
        {
            attributes = netcdf_root_attributes(file);
            nc_close(file);
        }
    }
    auto listed_databases = nlohmann::ordered_json::object();
    for (const auto& [database, attributes] : databases)
        listed_databases[database] = attributes;
    manifest["databases"] = listed_databases;
    std::vector<Entry> outputs;
    std::size_t unlisted = 0;
    {
        auto& written = detail::written_files();
        std::lock_guard<std::mutex> lock(written.mutex);
        for (const auto& output : written.paths)
        {
            if (output == log_file_ || output == path)
                continue;
            if (outputs.size() == max_listed_outputs)
            {
                unlisted++;
                continue;
            }
            outputs.push_back({"", output});
        }
    }
    manifest["outputs"] = describe_files(outputs, {});
    manifest["outputs_not_listed"] = unlisted;

    const auto& logger = Logger::instance();
    manifest["warnings"] = logger.warning_count();
    manifest["errors"] = logger.error_count();
    manifest["first_problems"] = logger.first_problems();
    manifest["exit_status"] = exit_status;
    manifest["seconds"] = seconds;

    // An output is known only once written; the run's results are kept, and
    // the manifest refused.
    {
        auto& written = detail::written_files();
        std::lock_guard<std::mutex> lock(written.mutex);
        for (const auto& output : written.paths)
            if (output != log_file_ && same_file(path, output))
                input_error("--manifest would overwrite the output ", output, ". Give another file name.");
    }
    ensure_parent_directory(path);
    std::ofstream out(path);
    if (!out)
        input_error("Cannot write --manifest file ", path, ".");
    out << manifest.dump(2) << '\n';
    if (!out.flush())
        input_error("Could not finish writing --manifest file ", path, ".");
}
