#pragma once

// Ground classes (--ground_classes): a few surfaces, each with its own ground
// model, that an observer file's sixth column (id,x,y,z,height,ground) picks
// per microphone -- a concrete pad among grass, say.  The file is a CSV keyed
// by an integer class id; each row gives the columns where its surface differs
// from the command line's ground, and inherits the rest.  Turbulence is the
// air's and stays run-wide; roughness is the surface's and goes with the class.
//
// read_ground_classes is a pure reader, as read_operations is: it checks the
// file's syntax and converts nothing.  ground_classes_from applies --units and
// the inheritance.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "ground_effect.hpp"
#include "util/constants.hpp"

namespace niceops {

// One row as the file gives it: empty optionals are cells (or columns) left to
// the command line.
struct GroundClassRow
{
    // 1-based line of the row in the file, for the messages.
    std::size_t line = 0;
    std::int64_t id = 0;
    // For the log; the id as text when the file gives none.
    std::string name;
    std::optional<ImpedanceModel> impedance_model;
    std::optional<double> flow_resistance; // kPa s/m^2
    std::optional<double> porosity_rate;   // 1/m, whatever --units says
    std::optional<double> layer_depth;     // in --units
    std::optional<double> roughness;       // in --units
};

struct GroundClassFile
{
    std::string path;
    std::vector<GroundClassRow> rows;
    // The normalized names of the header's columns.
    std::set<std::string> columns;
    // The length units the layer_depth and roughness headers declare,
    // lower-cased, or empty: checked against --units by ground_classes_from.
    std::string layer_depth_unit;
    std::string roughness_unit;
};

// Reads and checks `path`.  Refuses, naming the line: an unknown or repeated
// column (a turbulence column with its own reason), a missing class column, a
// header unit on a column that has none, a quote, a row whose field count is
// not the header's, an id that is not an integer from 0 to 2^31 - 1 or that
// repeats, a name outside [A-Za-z0-9_.-], an unknown impedance model, a flow
// resistance that is not positive, a negative porosity rate or roughness, a
// layer depth that is not positive, and a file with no rows.
GroundClassFile read_ground_classes(const std::string& path);

// A surface the run can evaluate: the command line's ground with the row's
// overrides applied, lengths in meters.
struct GroundClass
{
    std::int64_t id = 0;
    std::string name;
    std::size_t line = 0;
    GroundModel model;
};

// The file's classes in file order.  `ground` is the run's ground, already in
// meters; the file's layer_depth and roughness are in `unit`, and a header
// declaring another unit is refused, since nothing converts it.
std::vector<GroundClass> ground_classes_from(const GroundClassFile& file, const GroundModel& ground,
        LengthUnit unit);

// The class with `id`, or nullptr.
const GroundClass* find_ground_class(const std::vector<GroundClass>& classes, std::int64_t id);

} // namespace niceops
