#include <catch2/catch_test_macros.hpp>

#include <string>

#include "util/sha256.hpp"

using namespace niceops; // NOLINT(google-build-using-namespace) -- ported code reads like the D original

TEST_CASE("SHA-256 matches the FIPS 180-4 examples", "[manifest]")
{
    const auto digest = [](const std::string& text) {
        Sha256 hash;
        hash.update(text);
        return hash.hex();
    };
    CHECK(digest("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    CHECK(digest("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    // Two blocks: the padding spills into a second one.
    CHECK(digest("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")
            == "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
    // A million bytes, fed in pieces that straddle the block boundary.
    Sha256 million;
    const std::string piece(1000, 'a');
    for (int k = 0; k < 1000; k++)
        million.update(piece);
    CHECK(million.hex() == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}
