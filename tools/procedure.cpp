// Fly a procedure into a track niceops reads (src/procedure.hpp).
//
//   niceops_procedure approach.json -o tracks/approach.csv
//   niceops -t tracks/approach.csv -d databases/B407_fixed_load_factor.nod ...
//
// The procedure is JSON: a start (lat, lon, alt, speed, course) and a list of
// segments, every length and speed naming its unit in its key.  The manual's
// Companion tools section has the format.  The track is WGS84 with velocity and
// heading columns, lengths in --units, its alt on the start's datum.

#include <CLI/CLI.hpp>

#include <cstdio>
#include <exception>
#include <string>

#include "procedure.hpp"
#include "util/logger.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

int main(int argc, char* argv[])
{
    CLI::App app {"Fly a procedure into a WGS84 track for niceops."};
    std::string input, output;
    std::string unit_name = "ft";
    app.add_option("procedure", input, "Procedure file (JSON).")->required();
    app.add_option("-o,--output", output, "Track to write (CSV).")->required();
    app.add_option("--units", unit_name, "Length unit of the track: ft (default) or m.")
            ->check(CLI::IsMember({"ft", "m"}));
    CLI11_PARSE(app, argc, argv);
    try
    {
        const auto procedure = read_procedure(input);
        const auto states = fly(procedure);
        write_procedure_track(output, states, unit_name == "m" ? meters : feet);
        std::printf("Wrote %s: %zu states over %.1f s.\n", output.c_str(), states.size(),
                states.empty() ? 0.0 : states.back().t);
    }
    catch (const std::exception& problem)
    {
        error("Stopping: ", problem.what());
        return 2;
    }
    return 0;
}
