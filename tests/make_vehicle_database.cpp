// Writes a database that records its whole vehicle at the root, for the CTest
// cases that run without -v (CMakeLists.txt, root_vehicle_*): a copy of IN.nod
// with main_rotor_radius_meters, main_rotor_tip_speed_meters_per_sec,
// vehicle_weight_newtons and air_density_kg_m3 set, as panam writes them.
//
//   make_vehicle_database IN.nod OUT.nod RADIUS_M TIP_SPEED_MPS WEIGHT_N DENSITY_KG_M3
//
// Exits nonzero with the first failure on stderr.

#include <netcdf.h>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <utility>

namespace {

void check(int status, const char* what)
{
    if (status != NC_NOERR)
    {
        std::cerr << "make_vehicle_database: " << what << ": " << nc_strerror(status) << "\n";
        std::exit(1);
    }
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc != 7)
    {
        std::cerr << "usage: make_vehicle_database IN.nod OUT.nod RADIUS_M TIP_SPEED_MPS WEIGHT_N DENSITY_KG_M3\n";
        return 2;
    }
    const std::string in = argv[1], out = argv[2];
    std::error_code error;
    std::filesystem::create_directories(std::filesystem::path(out).parent_path(), error);
    std::filesystem::copy_file(in, out, std::filesystem::copy_options::overwrite_existing, error);
    if (error)
    {
        std::cerr << "make_vehicle_database: copying " << in << " to " << out << ": " << error.message() << "\n";
        return 1;
    }
    const std::pair<const char*, double> values[] = {{"main_rotor_radius_meters", std::atof(argv[3])},
            {"main_rotor_tip_speed_meters_per_sec", std::atof(argv[4])}, {"vehicle_weight_newtons", std::atof(argv[5])},
            {"air_density_kg_m3", std::atof(argv[6])}};
    int ncid = 0;
    check(nc_open(out.c_str(), NC_WRITE, &ncid), "opening the copy");
    check(nc_redef(ncid), "entering define mode");
    int ids[4] = {};
    for (int i = 0; i < 4; i++)
    {
        if (nc_inq_varid(ncid, values[i].first, &ids[i]) != NC_NOERR)
            check(nc_def_var(ncid, values[i].first, NC_DOUBLE, 0, nullptr, &ids[i]), values[i].first);
    }
    check(nc_enddef(ncid), "leaving define mode");
    for (int i = 0; i < 4; i++)
        check(nc_put_var_double(ncid, ids[i], &values[i].second), values[i].first);
    check(nc_close(ncid), "closing the copy");
    return 0;
}
