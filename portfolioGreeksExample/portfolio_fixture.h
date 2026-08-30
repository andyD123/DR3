#pragma once

#include "black_scholes_reference.h"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace dr3_example {

enum class OptionType { Call, Put };

struct Position {
    double spot = 0.0;
    double strike = 0.0;
    double maturity = 0.0;
    double vol = 0.0;
    double rate = 0.0;
    double quantity = 0.0;
    OptionType type = OptionType::Call;
};

inline bool is_valid_position(const Position& p)
{
    return p.spot > 0.0 && p.strike > 0.0 && p.maturity > 0.0 && p.vol > 0.0
        && std::isfinite(p.rate) && std::isfinite(p.quantity);
}

// Deterministic LCG so fixtures do not depend on <random> engines.
inline std::uint32_t lcg(std::uint32_t& state)
{
    state = state * 1664525u + 1013904223u;
    return state;
}

inline double unit(std::uint32_t& state)
{
    return (lcg(state) >> 8) * (1.0 / 16777216.0);
}

inline std::vector<Position> make_portfolio(int n, std::uint32_t seed = 1)
{
    if (n <= 0) {
        throw std::invalid_argument("portfolio size must be positive");
    }
    std::vector<Position> book;
    book.reserve(static_cast<size_t>(n));
    std::uint32_t state = seed;
    for (int i = 0; i < n; ++i) {
        Position p;
        p.spot = 80.0 + 40.0 * unit(state);
        p.strike = 70.0 + 50.0 * unit(state);
        p.maturity = 0.25 + 1.75 * unit(state);
        p.vol = 0.10 + 0.40 * unit(state);
        p.rate = 0.01 + 0.09 * unit(state);
        p.quantity = (unit(state) < 0.5 ? -1.0 : 1.0) * (1.0 + 9.0 * unit(state));
        p.type = (unit(state) < 0.5) ? OptionType::Call : OptionType::Put;
        if (!is_valid_position(p)) {
            throw std::logic_error("fixture produced an invalid position");
        }
        book.push_back(p);
    }
    return book;
}

inline std::vector<Position> exclude_invalid(const std::vector<Position>& in, int& rejected)
{
    std::vector<Position> out;
    rejected = 0;
    for (const auto& p : in) {
        if (is_valid_position(p)) {
            out.push_back(p);
        } else {
            rejected += 1;
        }
    }
    return out;
}

}  // namespace dr3_example
