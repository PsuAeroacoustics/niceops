#pragma once

// The demo database, loaded once per test process: a load costs a large share of a case that only
// reads it, and the cases that read it only call its const members.

#include "sphere_database.hpp"

namespace niceops::test {

inline const Database& demo_database()
{
    static const Database database("databases/S-76D_M3.nod");
    return database;
}

} // namespace niceops::test
