#pragma once

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

namespace dr3_embed {

inline std::vector<float> seeded_row(std::size_t d, std::uint32_t seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(-1.f, 1.f);
    std::vector<float> row(d);
    for (auto& v : row) {
        v = dist(rng);
    }
    return row;
}

// Deterministic row-major corpus. Seed 7 is part of the fixture contract.
inline void make_corpus(std::size_t n, std::size_t d, std::vector<float>& query, std::vector<float>& corpus)
{
    query = seeded_row(d, 7u);
    corpus.resize(n * d);
    for (std::size_t i = 0; i < n; ++i) {
        auto row = seeded_row(d, static_cast<std::uint32_t>(100 + i));
        std::copy(row.begin(), row.end(), corpus.begin() + static_cast<std::ptrdiff_t>(i * d));
    }
}

} // namespace dr3_embed
