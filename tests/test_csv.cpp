// Pins the header-consumption semantics of read_numeric_csv.
//
// D reads both the observer file and the trajectory CSV with
// csvReader!(double)(contents, null).  A null header argument makes std.csv
// consume the first record as a header and discard it UNCONDITIONALLY --
// including when that record holds perfectly good numbers.  The port must do
// the same, or a headerless observer file yields one more observer than the D
// original produced.  Verified against the D build: a 40-row headerless file
// gives 39 observers in both implementations.

#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <unistd.h>

#include "util/csv.hpp"

using namespace niceops;

namespace {

// Minimal scoped temp file so the test leaves nothing behind.
struct TempCsv
{
    std::string path;
    explicit TempCsv(const std::string& contents)
    {
        // mkstemp rather than tmpnam: tmpnam is deprecated and races.
        char pattern[] = "/tmp/niceops_test_XXXXXX";
        const int fd = ::mkstemp(pattern);
        REQUIRE(fd != -1);
        ::close(fd);
        path = pattern;
        std::ofstream out(path);
        out << contents;
    }
    ~TempCsv()
    {
        std::remove(path.c_str());
    }
};

} // namespace

TEST_CASE("read_numeric_csv skips a text header either way", "[csv]")
{
    TempCsv f("x,y,z\n1.0,2.0,3.0\n4.0,5.0,6.0\n");
    // Non-numeric lines are dropped regardless of the flag, so both agree here.
    REQUIRE(read_numeric_csv(f.path, false).size() == 2);
    REQUIRE(read_numeric_csv(f.path, true).size() == 2);
}

TEST_CASE("read_numeric_csv drops a numeric first record only when asked", "[csv]")
{
    TempCsv f("1.0,2.0,3.0\n4.0,5.0,6.0\n7.0,8.0,9.0\n");

    // Without the flag every numeric row is data.
    auto all = read_numeric_csv(f.path, false);
    REQUIRE(all.size() == 3);
    REQUIRE(all[0][0] == 1.0);

    // With it, the first record is eaten as a header -- this is the D
    // behaviour the observer and trajectory readers must reproduce.
    auto skipped = read_numeric_csv(f.path, true);
    REQUIRE(skipped.size() == 2);
    REQUIRE(skipped[0][0] == 4.0);
}

TEST_CASE("read_numeric_csv defaults to keeping every numeric record", "[csv]")
{
    TempCsv f("1.0,2.0\n3.0,4.0\n");
    REQUIRE(read_numeric_csv(f.path).size() == 2);
}
