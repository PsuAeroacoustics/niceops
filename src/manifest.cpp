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

void RunManifest::add(std::vector<Entry>& list, const std::string& role, const std::string& path)
{
    if (path.empty())
        return;
    const bool listed = std::any_of(list.begin(), list.end(),
            [&](const Entry& e) { return e.role == role && e.path == path; });
    if (!listed)
        list.push_back({role, path});
}

void RunManifest::add_input(const std::string& role, const std::string& path)
{
    add(inputs_, role, path);
}


namespace {

// Every root attribute and scalar variable of a netCDF file, text as it is and
// numbers as text;
// empty if the file cannot be opened (the loader will have said why).
std::map<std::string, std::string> root_attributes(const std::string& path)
{
    std::map<std::string, std::string> attributes;
    int file = 0;
    if (nc_open(path.c_str(), NC_NOWRITE, &file) != NC_NOERR)
        return attributes;
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
    nc_close(file);
    return attributes;
}

} // namespace

void RunManifest::add_database(const std::string& path)
{
    add_input("database", path);
    if (!path.empty() && databases_.count(path) == 0)
        databases_[path] = root_attributes(path);
}

namespace {

nlohmann::ordered_json describe_files(const std::vector<std::string>& roles,
        const std::vector<std::string>& paths)
{
    auto list = nlohmann::ordered_json::array();
    for (std::size_t k = 0; k < paths.size(); k++)
    {
        nlohmann::ordered_json entry;
        if (!roles[k].empty())
            entry["role"] = roles[k];
        entry["path"] = paths[k];
        std::error_code ec;
        const auto size = std::filesystem::file_size(paths[k], ec);
        if (ec)
        {
            entry["exists"] = false;
        }
        else
        {
            entry["bytes"] = size;
            const auto digest = sha256_of_file(paths[k]);
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

void RunManifest::refuse_overwriting_an_input(const std::string& path) const
{
    for (const auto& input : inputs_)
    {
        std::error_code ignored;
        if (std::filesystem::exists(path, ignored) && std::filesystem::equivalent(path, input.path, ignored))
            input_error("--manifest would overwrite the ", input.role, " ", input.path,
                    ". Give another file name.");
    }
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

    const auto split = [](const std::vector<Entry>& list, std::vector<std::string>& roles,
                               std::vector<std::string>& paths) {
        for (const auto& e : list)
        {
            roles.push_back(e.role);
            paths.push_back(e.path);
        }
    };
    std::vector<std::string> roles, paths;
    split(inputs_, roles, paths);
    manifest["inputs"] = describe_files(roles, paths);
    auto databases = nlohmann::ordered_json::object();
    for (const auto& [database, attributes] : databases_)
        databases[database] = attributes;
    manifest["databases"] = databases;
    roles.clear();
    paths.clear();
    std::size_t unlisted = 0;
    {
        auto& written = detail::written_files();
        std::lock_guard<std::mutex> lock(written.mutex);
        for (const auto& output : written.paths)
        {
            if (output == path || output == log_file_)
                continue;
            if (paths.size() == max_listed_outputs)
            {
                unlisted++;
                continue;
            }
            roles.emplace_back();
            paths.push_back(output);
        }
    }
    manifest["outputs"] = describe_files(roles, paths);
    manifest["outputs_not_listed"] = unlisted;

    const auto& logger = Logger::instance();
    manifest["warnings"] = logger.warning_count();
    manifest["errors"] = logger.error_count();
    manifest["first_problems"] = logger.first_problems();
    manifest["exit_status"] = exit_status;
    manifest["seconds"] = seconds;

    refuse_overwriting_an_input(path);
    ensure_parent_directory(path);
    std::ofstream out(path);
    if (!out)
        input_error("Cannot write --manifest file ", path, ".");
    out << manifest.dump(2) << '\n';
    if (!out.flush())
        input_error("Could not finish writing --manifest file ", path, ".");
}
