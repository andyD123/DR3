#include "european_pde.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace
{

using dr3::lattice::OptionType;
using dr3::lattice::PdeConfig;
using dr3::lattice::PdeScheme;
using dr3::lattice::VanillaOption;

VanillaOption standardCall()
{
    return {100.0, 100.0, 0.20, 0.05, 0.0, 1.0, OptionType::Call};
}

PdeConfig config(PdeScheme scheme,
                 std::size_t spaceSteps = 401,
                 std::size_t timeSteps = 400)
{
    return {scheme, spaceSteps, timeSteps, 0.0, 400.0};
}

double normalCdf(double value)
{
    return 0.5 * std::erfc(-value / std::sqrt(2.0));
}

double blackScholes(const VanillaOption& option)
{
    const double rootTime = std::sqrt(option.maturity);
    const double d1 = (std::log(option.spot / option.strike)
                       + (option.rate - option.dividendYield
                          + 0.5 * option.volatility * option.volatility) * option.maturity)
        / (option.volatility * rootTime);
    const double d2 = d1 - option.volatility * rootTime;
    const double discountedSpot = option.spot
        * std::exp(-option.dividendYield * option.maturity);
    const double discountedStrike = option.strike
        * std::exp(-option.rate * option.maturity);
    return option.type == OptionType::Call
        ? discountedSpot * normalCdf(d1) - discountedStrike * normalCdf(d2)
        : discountedStrike * normalCdf(-d2) - discountedSpot * normalCdf(-d1);
}

} // namespace

TEST(EuropeanPde, TerminalCallPayoff)
{
    const auto option = standardCall();
    const dr3::lattice::Grid1D grid{0.0, 200.0, 5};
    const auto payoff = dr3::lattice::europeanTerminalPayoff(option, grid);
    EXPECT_EQ(payoff, (std::vector<double>{0.0, 0.0, 0.0, 50.0, 100.0}));
}

TEST(EuropeanPde, TerminalPutPayoff)
{
    auto option = standardCall();
    option.type = OptionType::Put;
    const dr3::lattice::Grid1D grid{0.0, 200.0, 5};
    const auto payoff = dr3::lattice::europeanTerminalPayoff(option, grid);
    EXPECT_EQ(payoff, (std::vector<double>{100.0, 50.0, 0.0, 0.0, 0.0}));
}

TEST(EuropeanPde, CallLowerBoundary)
{
    const auto option = standardCall();
    const auto boundary = dr3::lattice::europeanBoundaryValues(option, 0.0, 400.0, 0.75);
    EXPECT_DOUBLE_EQ(boundary.lower, 0.0);
}

TEST(EuropeanPde, CallUpperBoundary)
{
    const auto option = standardCall();
    const auto boundary = dr3::lattice::europeanBoundaryValues(option, 0.0, 400.0, 0.75);
    EXPECT_NEAR(boundary.upper, 400.0 - 100.0 * std::exp(-0.05 * 0.75), 1.0e-13);
}

TEST(EuropeanPde, PutLowerBoundary)
{
    auto option = standardCall();
    option.type = OptionType::Put;
    const auto boundary = dr3::lattice::europeanBoundaryValues(option, 0.0, 400.0, 0.75);
    EXPECT_NEAR(boundary.lower, 100.0 * std::exp(-0.05 * 0.75), 1.0e-13);
}

TEST(EuropeanPde, PutUpperBoundary)
{
    auto option = standardCall();
    option.type = OptionType::Put;
    const auto boundary = dr3::lattice::europeanBoundaryValues(option, 0.0, 400.0, 0.75);
    EXPECT_DOUBLE_EQ(boundary.upper, 0.0);
}

TEST(EuropeanPde, ExplicitRejectsUnstableConfiguration)
{
    EXPECT_THROW(dr3::lattice::europeanPdePrice(
                     standardCall(), config(PdeScheme::Explicit, 401, 10)),
                 std::invalid_argument);
}

TEST(EuropeanPde, ExplicitStableConfigurationMatchesBlackScholes)
{
    auto option = standardCall();
    option.dividendYield = 0.02;
    const auto result = dr3::lattice::europeanPdePrice(
        option, config(PdeScheme::Explicit, 101, 800));
    EXPECT_NEAR(result.price, blackScholes(option), 0.10);
}

TEST(EuropeanPde, BackwardEulerMatchesBlackScholes)
{
    const auto option = standardCall();
    const auto result = dr3::lattice::europeanPdePrice(
        option, config(PdeScheme::BackwardEuler));
    EXPECT_NEAR(result.price, blackScholes(option), 0.02);
}

TEST(EuropeanPde, CrankNicolsonMatchesBlackScholes)
{
    const auto option = standardCall();
    const auto result = dr3::lattice::europeanPdePrice(
        option, config(PdeScheme::CrankNicolson));
    EXPECT_NEAR(result.price, blackScholes(option), 0.02);
}

TEST(EuropeanPde, RannacherMatchesBlackScholes)
{
    const auto option = standardCall();
    const auto result = dr3::lattice::europeanPdePrice(
        option, config(PdeScheme::CrankNicolsonRannacher));
    EXPECT_NEAR(result.price, blackScholes(option), 0.02);
}

TEST(EuropeanPde, PutCallParity)
{
    auto call = standardCall();
    call.dividendYield = 0.03;
    auto put = call;
    put.type = OptionType::Put;
    const auto pdeConfig = config(PdeScheme::CrankNicolsonRannacher);
    const double difference = dr3::lattice::europeanPdePrice(call, pdeConfig).price
        - dr3::lattice::europeanPdePrice(put, pdeConfig).price;
    const double parity = call.spot * std::exp(-call.dividendYield * call.maturity)
        - call.strike * std::exp(-call.rate * call.maturity);
    EXPECT_NEAR(difference, parity, 0.02);
}

TEST(EuropeanPde, CallSurfaceIsIncreasingInSpot)
{
    const auto result = dr3::lattice::europeanPdePrice(
        standardCall(), config(PdeScheme::CrankNicolsonRannacher));
    for (std::size_t index = 1; index < result.values.size(); ++index)
    {
        EXPECT_GE(result.values[index] + 1.0e-11, result.values[index - 1]);
    }
}

TEST(EuropeanPde, CallSurfaceIsConvexInSpot)
{
    const auto result = dr3::lattice::europeanPdePrice(
        standardCall(), config(PdeScheme::CrankNicolsonRannacher));
    for (std::size_t index = 1; index + 1 < result.values.size(); ++index)
    {
        EXPECT_GE(result.values[index - 1] - 2.0 * result.values[index]
                      + result.values[index + 1],
                  -1.0e-9);
    }
}

TEST(EuropeanPde, PutSurfaceIsDecreasingInSpot)
{
    auto put = standardCall();
    put.type = OptionType::Put;
    const auto result = dr3::lattice::europeanPdePrice(
        put, config(PdeScheme::CrankNicolsonRannacher));
    for (std::size_t index = 1; index < result.values.size(); ++index)
    {
        EXPECT_LE(result.values[index], result.values[index - 1] + 1.0e-11);
    }
}

TEST(EuropeanPde, RefinementReducesPricingError)
{
    const auto option = standardCall();
    const double reference = blackScholes(option);
    const double coarseError = std::abs(dr3::lattice::europeanPdePrice(
        option, config(PdeScheme::CrankNicolsonRannacher, 51, 50)).price - reference);
    const double fineError = std::abs(dr3::lattice::europeanPdePrice(
        option, config(PdeScheme::CrankNicolsonRannacher, 401, 400)).price - reference);
    EXPECT_LT(fineError, coarseError * 0.20);
}

TEST(EuropeanPde, NonMultipleOfSimdWidth)
{
    const auto option = standardCall();
    const auto result = dr3::lattice::europeanPdePrice(
        option, config(PdeScheme::CrankNicolsonRannacher, 402, 400));
    EXPECT_EQ(result.values.size(), 402u);
    EXPECT_NE(402u % 4u, 0u);
    EXPECT_NEAR(result.price, blackScholes(option), 0.02);
}

TEST(EuropeanPde, RejectsInvalidConfiguration)
{
    auto option = standardCall();
    EXPECT_THROW(dr3::lattice::europeanPdePrice(
                     option, config(PdeScheme::BackwardEuler, 2, 10)),
                 std::invalid_argument);
    auto invalidTime = config(PdeScheme::BackwardEuler);
    invalidTime.timeSteps = 0;
    EXPECT_THROW(dr3::lattice::europeanPdePrice(option, invalidTime), std::invalid_argument);
    auto outside = config(PdeScheme::BackwardEuler);
    outside.maximumSpot = 90.0;
    EXPECT_THROW(dr3::lattice::europeanPdePrice(option, outside), std::invalid_argument);
}

TEST(EuropeanPde, RepresentativeParameterGridProducesFiniteNonnegativePrices)
{
    constexpr std::array<double, 5> spotRatios{0.6, 0.8, 1.0, 1.2, 1.5};
    constexpr std::array<double, 3> volatilities{0.10, 0.20, 0.60};
    constexpr std::array<double, 3> maturities{0.10, 1.00, 5.00};
    constexpr std::array<double, 3> rates{-0.01, 0.00, 0.05};
    constexpr std::array<double, 2> dividends{0.00, 0.03};
    const auto pdeConfig = config(PdeScheme::CrankNicolsonRannacher, 101, 100);

    for (const double spotRatio : spotRatios)
    for (const double volatility : volatilities)
    for (const double maturity : maturities)
    for (const double rate : rates)
    for (const double dividend : dividends)
    {
        SCOPED_TRACE(::testing::Message()
                     << "S/K=" << spotRatio << ", vol=" << volatility
                     << ", T=" << maturity << ", r=" << rate
                     << ", q=" << dividend);
        const VanillaOption option{100.0 * spotRatio, 100.0, volatility,
                                   rate, dividend, maturity, OptionType::Call};
        const auto result = dr3::lattice::europeanPdePrice(option, pdeConfig);
        EXPECT_TRUE(std::isfinite(result.price));
        EXPECT_GE(result.price, 0.0);
        EXPECT_TRUE(std::all_of(result.values.begin(), result.values.end(),
                                [](double value)
                                {
                                    return std::isfinite(value) && value >= 0.0;
                                }));
    }
}
