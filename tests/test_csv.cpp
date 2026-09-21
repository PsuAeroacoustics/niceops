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
#include <catch2/catch_approx.hpp>

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

TEST_CASE("read_numeric_csv tolerates CRLF line endings", "[csv]")
{
    // A Windows-exported CSV used to parse as zero records: the trailing '\r'
    // made each line's last field non-numeric, so every row was silently
    // dropped and the crash surfaced much later, in the acoustics.
    TempCsv f("x,y,z\r\n1.0,2.0,3.0\r\n4.0,5.0,6.0\r\n");
    auto rows = read_numeric_csv(f.path, /*skip_first_record=*/true);
    REQUIRE(rows.size() == 2);
    REQUIRE(rows[0][2] == 3.0);
    REQUIRE(rows[1][0] == 4.0);
}

TEST_CASE("Labelled CSV keeps its header and matches names loosely", "[csv]")
{
    // The simpler trajectory formats are read by column name, so whatever a
    // planner already exports can be used without being reordered -- but only
    // when the name is recognised, and the matching has to survive the ways a
    // spreadsheet decorates a header.
    const TempCsv file("Time (s), Latitude ,LONGITUDE,alt_ft,junk\n"
                       "0,39.8,-74.4,1000,x\n"
                       "1,39.81,-74.41,1100,y\n");
    const auto csv = read_labelled_csv(file.path);
    REQUIRE(csv.header.size() == 5);
    // The rows carry a non-numeric field, so they are dropped exactly as
    // read_numeric_csv drops them -- the header still parses.
    REQUIRE(csv.rows.empty());

    REQUIRE(csv.column({"time"}) == 0);
    REQUIRE(csv.column({"lat", "latitude"}) == 1);
    REQUIRE(csv.column({"lon", "longitude"}) == 2);
    REQUIRE(csv.column({"alt ft"}) == 3);
    // Not every spelling is a match, and that is the point: an unrecognised
    // altitude column fails the file rather than being guessed at.
    REQUIRE(csv.column({"alt m"}) == LabelledCsv::npos);
    REQUIRE(csv.column({"utcsec"}) == LabelledCsv::npos);

    // Normalisation drops case, surrounding space, punctuation and a
    // parenthesised or bracketed unit.
    REQUIRE(LabelledCsv::normalise("  Alt (ft) ") == "alt");
    REQUIRE(LabelledCsv::normalise("alt_ft") == "altft");
    REQUIRE(LabelledCsv::normalise("V East [m/s]") == "veast");

    const TempCsv numeric("t,lat,lon,alt\n0,39.8,-74.4,1000\n1,39.81,-74.41,1100\n");
    const auto rows = read_labelled_csv(numeric.path);
    REQUIRE(rows.header.size() == 4);
    REQUIRE(rows.rows.size() == 2);
    REQUIRE(rows.rows[1][3] == 1100.0);
    // The same rows read_numeric_csv sees: one parser, one set of rows, so the
    // positional AC path is unchanged by being read through this.
    REQUIRE(read_numeric_csv(numeric.path, /*skip_first_record=*/true).size() == rows.rows.size());
}

TEST_CASE("fields padded with whitespace still parse", "[csv]")
{
    // std::stod already skips what precedes a number, so rejecting only what
    // follows it made the two ends behave differently.  Every shipped AC track
    // pads on the left and so parsed anyway; arrays/coyle_approach.csv pads on
    // the right and had every row dropped, leaving a scattered-observer run
    // that read zero observers and still exited 0.
    const TempCsv padded("id,x,y,z\n1 ,-794.12,-1928.64 ,10.69\n2 , -788.33 ,-1085.09,-1.08\n");
    const auto rows = read_numeric_csv(padded.path, /*skip_first_record=*/true);
    REQUIRE(rows.size() == 2);
    REQUIRE(rows[0].size() == 4);
    REQUIRE(rows[0][1] == Catch::Approx(-794.12));
    REQUIRE(rows[0][2] == Catch::Approx(-1928.64));
    REQUIRE(rows[1][1] == Catch::Approx(-788.33));

    // Genuinely non-numeric fields are still rejected, whitespace or not.
    const TempCsv junk("a,b\n1,2\n3, four \n");
    REQUIRE(read_numeric_csv(junk.path, /*skip_first_record=*/true).size() == 1);
}

TEST_CASE("the shipped scattered observer file reads", "[csv]")
{
    // It is the only scattered-observer set in the repository and was entirely
    // unreadable: every row is padded before its commas.
    const auto rows = read_numeric_csv("arrays/coyle_approach.csv", /*skip_first_record=*/true);
    REQUIRE(rows.size() > 24);
    for (const auto& row : rows)
        REQUIRE(row.size() == 4);
}
