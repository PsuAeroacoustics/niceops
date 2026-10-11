// The operations file reader (--operations): its columns, defaults, path
// resolution, clock parsing and every refusal, on files written to a scratch
// folder.  What a run does with the rows is checked by the operations_* CTests.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <unistd.h>

#include "operations.hpp"

using namespace niceops;
using Catch::Matchers::ContainsSubstring;

namespace {

// A scratch folder holding an operations file and the (empty) files it names,
// removed with everything in it.
struct Folder
{
    std::filesystem::path root;
    Folder()
    {
        char pattern[] = "/tmp/niceops_operations_XXXXXX";
        REQUIRE(::mkdtemp(pattern) != nullptr);
        root = pattern;
        std::filesystem::create_directories(root / "tracks");
        for (const char* name : {"tracks/a.csv", "tracks/b.csv", "a.nod", "b.nod", "a.json", "air.csv"})
            std::ofstream(root / name) << "";
    }
    ~Folder()
    {
        std::error_code ec;
        std::filesystem::remove_all(root, ec);
    }
    std::string write(const std::string& contents) const
    {
        const auto path = root / "ops.csv";
        std::ofstream(path) << contents;
        return path.string();
    }
};

} // namespace

TEST_CASE("an operations file reads its columns by their normalized names", "[operations]")
{
    Folder folder;
    const auto path = folder.write("# comment above the header\n"
                                   "\n"
                                   "Name, Track, COUNT, Local_Time, Database, Vehicle, Offset (ft), Humidity (%)\r\n"
                                   "dep, tracks/a.csv, 6.5, 22:30, a.nod, a.json, -500, 40\n"
                                   "# a row left out\n"
                                   "arr, tracks/b.csv, 2, 07:05:09, b.nod, a.json, 0, 60\n");
    const OperationsFile file = read_operations(path);
    REQUIRE(file.rows.size() == 2);
    REQUIRE(file.has("localtime"));
    REQUIRE(file.has("offset"));
    REQUIRE(file.offset_unit == "ft");
    const Operation& dep = file.rows[0];
    REQUIRE(dep.line == 4);
    REQUIRE(dep.name == "dep");
    REQUIRE(dep.track == (folder.root / "tracks/a.csv").string());
    REQUIRE(dep.database == (folder.root / "a.nod").string());
    REQUIRE(dep.count == 6.5);
    REQUIRE(dep.local_time == 22.0 * 3600.0 + 30.0 * 60.0);
    REQUIRE(dep.offset == -500.0);
    REQUIRE(dep.humidity == 40.0);
    REQUIRE(file.rows[1].line == 6);
    REQUIRE(file.rows[1].local_time == 7.0 * 3600.0 + 5.0 * 60.0 + 9.0);
}

TEST_CASE("an operations file's empty cells and absent columns take the defaults", "[operations]")
{
    Folder folder;
    const auto path = folder.write("track,count,name,database\ntracks/a.csv,,,\n");
    const Operation row = read_operations(path).rows.front();
    REQUIRE(row.name == "a");
    REQUIRE_FALSE(row.count);
    REQUIRE(operation_count(row) == 1.0);
    REQUIRE(row.database.empty());
    REQUIRE(row.vehicle.empty());
    REQUIRE_FALSE(row.local_time);
    REQUIRE_FALSE(row.periods);
    REQUIRE_FALSE(row.offset);
    REQUIRE(row.wind_text.empty());
}

TEST_CASE("an operations file's paths are its own folder's, and absolute ones stand", "[operations]")
{
    Folder folder;
    std::filesystem::create_directories(folder.root / "sub");
    const auto absolute = (folder.root / "tracks/b.csv").string();
    std::ofstream(folder.root / "sub/ops.csv") << "track\n../tracks/a.csv\n" << absolute << "\n";
    const auto file = read_operations((folder.root / "sub/ops.csv").string());
    REQUIRE(file.rows[0].track == (folder.root / "tracks/a.csv").string());
    REQUIRE(file.rows[1].track == absolute);
}

TEST_CASE("an operations file joins its wind columns into --wind's text", "[operations]")
{
    Folder folder;
    const auto path = folder.write("track,wind_speed,wind_direction,wind_height (m)\n"
                                   "tracks/a.csv,8kt,200,10\ntracks/b.csv,3mps,180,\n");
    const auto file = read_operations(path);
    REQUIRE(file.rows[0].wind_text == "8kt,200,10");
    REQUIRE(file.rows[1].wind_text == "3mps,180");
    REQUIRE(file.wind_height_unit == "m");
}

TEST_CASE("period counts default to zero and add up to the row's operations", "[operations]")
{
    Folder folder;
    const auto path = folder.write("track,day,evening,night\ntracks/a.csv,2,,0.5\n");
    const Operation row = read_operations(path).rows.front();
    REQUIRE(row.periods);
    REQUIRE((*row.periods)[0] == 2.0);
    REQUIRE((*row.periods)[1] == 0.0);
    REQUIRE((*row.periods)[2] == 0.5);
    REQUIRE(operation_count(row) == 2.5);
}

TEST_CASE("local_time is HH:MM or HH:MM:SS within one day", "[operations]")
{
    REQUIRE(parse_local_clock("00:00") == 0.0);
    REQUIRE(parse_local_clock("7:30") == 7.0 * 3600.0 + 30.0 * 60.0);
    REQUIRE(parse_local_clock("23:59:59") == 86399.0);
    for (const char* bad : {"24:00", "25:00", "7pm", "7:30 am", "12:60", "12:00:60", "12", "12:0", ":30", "123:00",
                 "12:00:", "-1:00", ""})
        REQUIRE_THROWS_WITH(parse_local_clock(bad), ContainsSubstring("is not a clock time"));
}

TEST_CASE("an operations file refuses what it cannot read, naming the line", "[operations]")
{
    Folder folder;
    const auto refuses = [&](const std::string& contents, const std::string& message) {
        INFO(contents);
        REQUIRE_THROWS_WITH(read_operations(folder.write(contents)), ContainsSubstring(message));
    };
    refuses("track,cnt\ntracks/a.csv,1\n", "unknown column \"cnt\"");
    refuses("track,Count,count\ntracks/a.csv,1,1\n", "column \"count\" is given twice");
    refuses("name,count\na,1\n", "the header has no track column");
    refuses("track,count\ntracks/a.csv\n", "line 2 has 1 fields, and the header has 2");
    refuses("track,count\ntracks/a.csv,1,\n", "line 2 has 3 fields");
    refuses("track\n\"tracks/a.csv\"\n", "line 2 holds a quote");
    refuses("track,count\ntracks/a.csv,0\n", "count 0 is not positive");
    refuses("track,count\ntracks/a.csv,-2\n", "count -2 is not positive");
    refuses("track,count\ntracks/a.csv,nan\n", "\"nan\", which is not a finite number");
    refuses("track,count\ntracks/a.csv,inf\n", "not a finite number");
    refuses("track,count\ntracks/a.csv,6.2 ops\n", "\"6.2 ops\", which is not a finite number");
    refuses("track,count,day\ntracks/a.csv,1,1\n", "as count or as day, evening and night counts, not both");
    refuses("track,local_time,night\ntracks/a.csv,12:00,1\n", "give one or the other");
    refuses("track,day,evening,night\ntracks/a.csv,0,0,0\n", "add up to no operations");
    refuses("track,day\ntracks/a.csv,-1\n", "the day count -1 is negative");
    refuses("track,local_time\ntracks/a.csv,25:00\n", "line 2: local_time \"25:00\" is not a clock time");
    refuses("track,local_time\ntracks/a.csv,7pm\n", "local_time \"7pm\" is not a clock time");
    refuses("track,start_utc\ntracks/a.csv,3600\n", "not start_utc");
    refuses("track\ntracks/a.csv\ntracks/../tracks/a.csv\n", "lines 2 and 3 share the name a; add a name column");
    refuses("name,track\ndep,tracks/a.csv\nDEP,tracks/b.csv\n", "share the name DEP; give them distinct names");
    refuses("name,track\n.hidden,tracks/a.csv\n", "may not start with a dot");
    refuses("name,track\na b,tracks/a.csv\n", "may hold only letters, digits");
    refuses("track\ntracks/missing.csv\n", "line 2 (missing): the track file " + (folder.root / "tracks/missing.csv").string()
                    + " does not exist");
    refuses("track,database\ntracks/a.csv,c.nod\n", "the database file");
    refuses("track,atmosphere\ntracks/a.csv,nowhere.csv\n", "the atmosphere file");
    refuses("track,wind_speed\ntracks/a.csv,8kt\n", "wind_speed and wind_direction go together");
    refuses("track,wind_height\ntracks/a.csv,10\n", "needs wind_speed");
    refuses("track,wind_speed,wind_direction\ntracks/a.csv,8kt,\n", "give wind_speed and wind_direction together");
    refuses("track,temperature (C)\ntracks/a.csv,15\n", "the temperature column is in K");
    refuses("track,count\n", "has no operations");
    refuses("# only a comment\n", "has no header");
    REQUIRE_THROWS_WITH(read_operations((folder.root / "absent.csv").string()),
            ContainsSubstring("Unable to open the operations file"));
}
