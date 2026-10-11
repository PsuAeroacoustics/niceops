#pragma once

// SHA-256 (FIPS 180-4), for the run manifest's input and output digests.  Small
// and dependency-free: the manifest needs a digest a user can check with
// `shasum -a 256`, not speed, and one hash of the largest shipped database
// takes about 0.3 s (72 MB).

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <optional>
#include <string>
#include <vector>

namespace niceops {

class Sha256
{
public:
    Sha256()
    {
        state_ = {0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au, 0x510e527fu, 0x9b05688cu,
                0x1f83d9abu, 0x5be0cd19u};
    }

    void update(const unsigned char* data, std::size_t size)
    {
        bits_ += static_cast<std::uint64_t>(size) * 8u;
        while (size > 0)
        {
            const std::size_t take = std::min(size, block_.size() - filled_);
            for (std::size_t k = 0; k < take; k++)
                block_[filled_ + k] = data[k];
            filled_ += take;
            data += take;
            size -= take;
            if (filled_ == block_.size())
            {
                compress();
                filled_ = 0;
            }
        }
    }

    void update(const std::string& text)
    {
        update(reinterpret_cast<const unsigned char*>(text.data()), text.size());
    }

    // The digest as 64 lowercase hex digits.  Ends the computation.
    std::string hex()
    {
        const std::uint64_t bits = bits_;
        const unsigned char one = 0x80;
        update(&one, 1);
        const unsigned char zero = 0;
        while (filled_ != 56)
            update(&zero, 1);
        unsigned char length[8];
        for (int k = 0; k < 8; k++)
            length[k] = static_cast<unsigned char>(bits >> (56 - 8 * k));
        update(length, 8);
        static const char* digits = "0123456789abcdef";
        std::string out;
        for (const std::uint32_t word : state_)
            for (int shift = 28; shift >= 0; shift -= 4)
                out += digits[(word >> shift) & 0xfu];
        return out;
    }

private:
    static std::uint32_t rotr(std::uint32_t x, int n)
    {
        return (x >> n) | (x << (32 - n));
    }

    void compress()
    {
        static constexpr std::array<std::uint32_t, 64> k = {0x428a2f98u, 0x71374491u, 0xb5c0fbcfu,
                0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u,
                0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u,
                0xc19bf174u, 0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu,
                0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau, 0x983e5152u, 0xa831c66du, 0xb00327c8u,
                0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u, 0x27b70a85u,
                0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu,
                0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u,
                0xd6990624u, 0xf40e3585u, 0x106aa070u, 0x19a4c116u, 0x1e376c08u, 0x2748774cu,
                0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu,
                0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u,
                0xc67178f2u};
        std::array<std::uint32_t, 64> w {};
        for (std::size_t t = 0; t < 16; t++)
            w[t] = (static_cast<std::uint32_t>(block_[4 * t]) << 24)
                    | (static_cast<std::uint32_t>(block_[4 * t + 1]) << 16)
                    | (static_cast<std::uint32_t>(block_[4 * t + 2]) << 8)
                    | static_cast<std::uint32_t>(block_[4 * t + 3]);
        for (std::size_t t = 16; t < 64; t++)
        {
            const std::uint32_t s0 = rotr(w[t - 15], 7) ^ rotr(w[t - 15], 18) ^ (w[t - 15] >> 3);
            const std::uint32_t s1 = rotr(w[t - 2], 17) ^ rotr(w[t - 2], 19) ^ (w[t - 2] >> 10);
            w[t] = w[t - 16] + s0 + w[t - 7] + s1;
        }
        std::uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        std::uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
        for (std::size_t t = 0; t < 64; t++)
        {
            const std::uint32_t t1 = h + (rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25)) + ((e & f) ^ (~e & g))
                    + k[t] + w[t];
            const std::uint32_t t2 = (rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
            h = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
        }
        state_[0] += a;
        state_[1] += b;
        state_[2] += c;
        state_[3] += d;
        state_[4] += e;
        state_[5] += f;
        state_[6] += g;
        state_[7] += h;
    }

    std::array<std::uint32_t, 8> state_ {};
    std::array<unsigned char, 64> block_ {};
    std::size_t filled_ = 0;
    std::uint64_t bits_ = 0;
};

// The SHA-256 of a file's bytes, or nullopt if it cannot be read.
inline std::optional<std::string> sha256_of_file(const std::string& path)
{
    std::FILE* file = std::fopen(path.c_str(), "rb");
    if (!file)
        return std::nullopt;
    Sha256 hash;
    std::vector<unsigned char> buffer(1 << 20);
    std::size_t got = 0;
    while ((got = std::fread(buffer.data(), 1, buffer.size(), file)) > 0)
        hash.update(buffer.data(), got);
    const bool failed = std::ferror(file) != 0;
    std::fclose(file);
    if (failed)
        return std::nullopt;
    return hash.hex();
}

} // namespace niceops
