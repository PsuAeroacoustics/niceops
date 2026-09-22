#pragma once

// Mathematical constants shared across src/.
//
// M_PI is a POSIX extension, not part of the C++ standard: glibc's math.h hides
// it under a strict -std=c++17, which is what CXX_EXTENSIONS OFF selects.  Four
// files had each defined their own PI for that reason while two others used
// M_PI anyway, so the codebase both documented the hazard and walked into it.

namespace niceops {

inline constexpr double pi = 3.14159265358979323846;

// Degrees and radians, spelled once.  The port converts in both directions all
// over, and `180.0 / PI` written out is easy to misread at a glance.
inline constexpr double deg_to_rad = pi / 180.0;
inline constexpr double rad_to_deg = 180.0 / pi;

// Feet to metres, exactly.  The computation is metric throughout; feet exist
// only at the edges, where files and the command line are read and written.
inline constexpr double ft2m = 0.3048;

// The length unit the interface speaks, selected by --units.  `to_metres`
// multiplies a user-supplied number to get metres, and divides an internal one
// on the way back out; `name` is what the output files and figure axes are
// labelled with.
//
// It covers the command line and every file that carries a length -- trajectory
// input, observer lists, level output, time histories.  The one exception is
// the state file written by -S, which is metric whatever this says: it is a
// dump of the internal state and golden/state_197421.csv, the D-derived oracle,
// encodes it that way.
struct LengthUnit
{
    double to_metres = ft2m;
    const char* name = "ft";
};

inline constexpr LengthUnit feet {ft2m, "ft"};
inline constexpr LengthUnit metres {1.0, "m"};

} // namespace niceops
