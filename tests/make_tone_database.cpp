// Writes a format 2 database for the CTest cases that run tones (CMakeLists.txt,
// tones_*): a copy of a spectral database with pure-tone sources added to every
// condition group, as panam's PT build lays them out.
//
//   make_tone_database IN.nod OUT.nod [--sources N] [--sync_group NAME|none]
//                      [--phase_reference rotor_azimuth|relative|none]
//                      [--band_content broadband|total]
//
// Each source is a rotor at 43.2 rad/s (shaft orders 16 and 32: 110 and 220 Hz),
// its hub 1.5 m to either side, with a smooth directivity (95 dB less 6 cos^2 of
// the polar angle, the phase turning with azimuth) and alternate senses of
// rotation.  Exits nonzero with the first failure on stderr.

#include <netcdf.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr double pi = 3.14159265358979323846;

void check(int status, const char* what)
{
    if (status != NC_NOERR)
    {
        std::cerr << "make_tone_database: " << what << ": " << nc_strerror(status) << "\n";
        std::exit(1);
    }
}

void put_text(int ncid, int varid, const char* name, const std::string& value)
{
    check(nc_put_att_text(ncid, varid, name, value.size(), value.c_str()), name);
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "usage: make_tone_database IN.nod OUT.nod [--sources N] [--sync_group NAME|none] "
                     "[--phase_reference KIND] [--band_content broadband|total]\n";
        return 2;
    }
    const std::string in = argv[1], out = argv[2];
    int n_sources = 2;
    std::string sync_group = "main", phase_reference = "rotor_azimuth", band_content = "broadband";
    for (int i = 3; i + 1 < argc; i += 2)
    {
        const std::string key = argv[i], value = argv[i + 1];
        if (key == "--sources")
            n_sources = std::atoi(value.c_str());
        else if (key == "--sync_group")
            sync_group = value == "none" ? "" : value;
        else if (key == "--phase_reference")
            phase_reference = value;
        else if (key == "--band_content")
            band_content = value;
        else
        {
            std::cerr << "make_tone_database: unknown option " << key << "\n";
            return 2;
        }
    }
    std::filesystem::copy_file(in, out, std::filesystem::copy_options::overwrite_existing);

    int ncid = 0;
    check(nc_open(out.c_str(), NC_WRITE, &ncid), "open");
    check(nc_redef(ncid), "redef");
    std::vector<std::string> tags;
    std::string list;
    for (int s = 0; s < n_sources; s++)
    {
        tags.push_back(s % 2 == 0 ? "left" + std::to_string(s / 2) : "right" + std::to_string(s / 2));
        list += (list.empty() ? "" : ",") + tags.back();
    }
    put_text(ncid, NC_GLOBAL, "sources", list);
    put_text(ncid, NC_GLOBAL, "band_content", band_content);
    int n_groups = 0;
    check(nc_inq_grps(ncid, &n_groups, nullptr), "groups");
    std::vector<int> groups(static_cast<std::size_t>(n_groups));
    check(nc_inq_grps(ncid, &n_groups, groups.data()), "groups");

    struct Vars
    {
        int group, center, rate, harmonic, frequency, amplitude, phase;
        int s;
        std::size_t n_phi, n_theta;
    };
    std::vector<Vars> vars;
    const std::vector<int> harmonics = {16, 32};
    for (const int g : groups)
    {
        int phi_dim = 0, theta_dim = 0;
        check(nc_inq_dimid(g, "PHI", &phi_dim), "PHI");
        check(nc_inq_dimid(g, "THETA", &theta_dim), "THETA");
        std::size_t n_phi = 0, n_theta = 0;
        check(nc_inq_dimlen(g, phi_dim, &n_phi), "PHI length");
        check(nc_inq_dimlen(g, theta_dim, &n_theta), "THETA length");
        for (int s = 0; s < n_sources; s++)
        {
            Vars v {};
            v.s = s;
            v.n_phi = n_phi;
            v.n_theta = n_theta;
            check(nc_def_grp(g, tags[static_cast<std::size_t>(s)].c_str(), &v.group), "source group");
            int tone_dim = 0, xyz_dim = 0;
            check(nc_def_dim(v.group, "tone", harmonics.size(), &tone_dim), "tone");
            check(nc_def_dim(v.group, "xyz", 3, &xyz_dim), "xyz");
            check(nc_def_var(v.group, "center", NC_DOUBLE, 1, &xyz_dim, &v.center), "center");
            check(nc_def_var(v.group, "shaft_rate", NC_DOUBLE, 0, nullptr, &v.rate), "shaft_rate");
            check(nc_def_var(v.group, "harmonic", NC_INT, 1, &tone_dim, &v.harmonic), "harmonic");
            check(nc_def_var(v.group, "frequency", NC_DOUBLE, 1, &tone_dim, &v.frequency), "frequency");
            const int dims[3] = {phi_dim, theta_dim, tone_dim};
            check(nc_def_var(v.group, "amplitude", NC_DOUBLE, 3, dims, &v.amplitude), "amplitude");
            check(nc_def_var(v.group, "phase", NC_DOUBLE, 3, dims, &v.phase), "phase");
            put_text(v.group, NC_GLOBAL, "phase_reference", phase_reference);
            if (!sync_group.empty())
                put_text(v.group, NC_GLOBAL, "sync_group", sync_group);
            vars.push_back(v);
        }
    }
    check(nc_enddef(ncid), "enddef");
    int version_var = 0;
    check(nc_inq_varid(ncid, "database_version", &version_var), "database_version");
    const int version = 2;
    check(nc_put_var_int(ncid, version_var, &version), "database_version");
    for (const auto& v : vars)
    {
        const double rate = (v.s % 2 == 0 ? 1.0 : -1.0) * 43.2;
        const double center[3] = {0.0, v.s % 2 == 0 ? -1.5 : 1.5, 0.0};
        check(nc_put_var_double(v.group, v.center, center), "center");
        check(nc_put_var_double(v.group, v.rate, &rate), "shaft_rate");
        check(nc_put_var_int(v.group, v.harmonic, harmonics.data()), "harmonic");
        std::vector<double> frequency;
        for (const int m : harmonics)
            frequency.push_back(m * std::fabs(rate) / (2.0 * pi));
        check(nc_put_var_double(v.group, v.frequency, frequency.data()), "frequency");
        std::vector<double> amplitude, phase;
        for (std::size_t i = 0; i < v.n_phi; i++)
            for (std::size_t j = 0; j < v.n_theta; j++)
                for (std::size_t k = 0; k < harmonics.size(); k++)
                {
                    const double polar = pi * static_cast<double>(j) / static_cast<double>(v.n_theta - 1);
                    const double azimuth = 2.0 * pi * static_cast<double>(i) / static_cast<double>(v.n_phi);
                    amplitude.push_back(95.0 - 6.0 * std::cos(polar) * std::cos(polar) - 4.0 * static_cast<double>(k));
                    phase.push_back(20.0 * std::sin(azimuth) + 10.0 * static_cast<double>(k));
                }
        check(nc_put_var_double(v.group, v.amplitude, amplitude.data()), "amplitude");
        check(nc_put_var_double(v.group, v.phase, phase.data()), "phase");
    }
    check(nc_close(ncid), "close");
    return 0;
}
