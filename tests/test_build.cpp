#include <catch2/catch_test_macros.hpp>

TEST_CASE("asserts are live when NICEOPS_ASSERTS asks for them", "[build]")
{
    // The test build (NICEOPS_ASSERTS=ON, as CI and CLAUDE.md configure it)
    // keeps assert() live, and the sphere database's isfinite guards are part of
    // what the suite checks.  This file is compiled with the project's flags, so
    // it fails if an explicit -DCMAKE_CXX_FLAGS=-DNDEBUG, or CMake's own flags,
    // compile the asserts out of a build that asked for them.  A shipped build
    // (NICEOPS_ASSERTS off) has NDEBUG and no requirement here.
#if defined(NICEOPS_ASSERTS_REQUIRED) && defined(NDEBUG)
    FAIL("NDEBUG is defined in a build that asked for live asserts (NICEOPS_ASSERTS=ON)");
#else
    SUCCEED();
#endif
}
