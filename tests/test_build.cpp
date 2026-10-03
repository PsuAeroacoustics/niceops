#include <catch2/catch_test_macros.hpp>

TEST_CASE("assert() is live in this build", "[build]")
{
    // CMakeLists.txt and CLAUDE.md say the asserts stay on (the port contract),
    // and sphere_database.cpp's isfinite guards rely on it.  CMake's default
    // Release flags define NDEBUG, which compiled all 17 of them out; the project
    // now strips it.  This file is compiled with the project's flags, so it fails
    // if that stops being true -- an explicit -DCMAKE_CXX_FLAGS=-DNDEBUG included.
#ifdef NDEBUG
    FAIL("NDEBUG is defined: assert() is compiled out, contrary to CMakeLists.txt and CLAUDE.md");
#else
    SUCCEED();
#endif
}
