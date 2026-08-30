#pragma once

#include <cstddef>

namespace dr3::lattice
{

enum class OptionType
{
    Call,
    Put
};

struct VanillaOption
{
    double spot;
    double strike;
    double volatility;
    double rate;
    double dividendYield;
    double maturity;
    OptionType type;
};

struct TreeConfig
{
    std::size_t steps;
};

} // namespace dr3::lattice
