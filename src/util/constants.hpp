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

} // namespace niceops
