#include "tree_pricers.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{

using dr3::lattice::OptionType;
using dr3::lattice::TreeConfig;
using dr3::lattice::VanillaOption;

VanillaOption callOption()
{
    return {100.0, 100.0, 0.20, 0.05, 0.0, 1.0, OptionType::Call};
}

double normalCdf(double value)
{
    return 0.5 * std::erfc(-value / std::sqrt(2.0));
}

double blackScholes(const VanillaOption& option)
{
    if (option.maturity == 0.0)
    {
        return option.type == OptionType::Call
            ? std::max(option.spot - option.strike, 0.0)
            : std::max(option.strike - option.spot, 0.0);
    }

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
    if (option.type == OptionType::Call)
    {
        return discountedSpot * normalCdf(d1) - discountedStrike * normalCdf(d2);
    }
    return discountedStrike * normalCdf(-d2) - discountedSpot * normalCdf(-d1);
}

} // namespace

TEST(TreePricer, RejectsZeroSteps)
{
    EXPECT_THROW(dr3::lattice::europeanBinomial(callOption(), {0}), std::invalid_argument);
}

TEST(TreePricer, RejectsInvalidStrike)
{
    auto option = callOption();
    option.strike = 0.0;
    EXPECT_THROW(dr3::lattice::americanTrinomial(option, {10}), std::invalid_argument);
}

TEST(TreePricer, RejectsInvalidVolatility)
{
    auto option = callOption();
    option.volatility = 0.0;
    EXPECT_THROW(dr3::lattice::europeanBinomial(option, {10}), std::invalid_argument);
}

TEST(TreePricer, RejectsOtherInvalidOptionInputs)
{
    auto option = callOption();
    option.spot = -1.0;
    EXPECT_THROW(dr3::lattice::europeanTrinomial(option, {10}), std::invalid_argument);
    option = callOption();
    option.maturity = -1.0;
    EXPECT_THROW(dr3::lattice::europeanTrinomial(option, {10}), std::invalid_argument);
    option = callOption();
    option.rate = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(dr3::lattice::europeanBinomial(option, {10}), std::invalid_argument);
}

TEST(TreePricer, RejectsInvalidProbabilitiesAndBarrierInputs)
{
    auto option = callOption();
    option.rate = 1.0;
    option.volatility = 0.01;
    EXPECT_THROW(dr3::lattice::europeanBinomial(option, {1}), std::domain_error);

    option = callOption();
    EXPECT_THROW(dr3::lattice::upAndOutTrinomial(option, {10}, 0.0), std::invalid_argument);
    EXPECT_THROW(dr3::lattice::upAndOutTrinomial(option, {10}, 120.0, -1.0),
                 std::invalid_argument);
}

TEST(TreePricer, MaturityZeroReturnsIntrinsic)
{
    auto call = callOption();
    call.spot = 110.0;
    call.maturity = 0.0;
    auto put = call;
    put.spot = 90.0;
    put.type = OptionType::Put;

    EXPECT_DOUBLE_EQ(dr3::lattice::europeanBinomial(call, {1}), 10.0);
    EXPECT_DOUBLE_EQ(dr3::lattice::europeanTrinomial(call, {2}), 10.0);
    EXPECT_DOUBLE_EQ(dr3::lattice::americanTrinomial(put, {3}), 10.0);
    EXPECT_DOUBLE_EQ(dr3::lattice::upAndOutTrinomial(call, {3}, 120.0), 10.0);
}

TEST(TreePricer, SupportsOddStepCounts)
{
    constexpr std::array<std::size_t, 8> stepCounts{
        1, 3, 63, 65, 255, 257, 1023, 1025};
    const auto option = callOption();

    for (const auto steps : stepCounts)
    {
        SCOPED_TRACE(steps);
        const double binomial = dr3::lattice::europeanBinomial(option, {steps});
        const double trinomial = dr3::lattice::europeanTrinomial(option, {steps});
        EXPECT_TRUE(std::isfinite(binomial));
        EXPECT_TRUE(std::isfinite(trinomial));
        EXPECT_GE(binomial, 0.0);
        EXPECT_GE(trinomial, 0.0);
    }
}

TEST(TreePricer, SupportsEvenStepCounts)
{
    constexpr std::array<std::size_t, 4> stepCounts{2, 64, 256, 1024};
    const auto option = callOption();

    for (const auto steps : stepCounts)
    {
        SCOPED_TRACE(steps);
        const double binomial = dr3::lattice::europeanBinomial(option, {steps});
        const double trinomial = dr3::lattice::europeanTrinomial(option, {steps});
        EXPECT_TRUE(std::isfinite(binomial));
        EXPECT_TRUE(std::isfinite(trinomial));
        EXPECT_GE(binomial, 0.0);
        EXPECT_GE(trinomial, 0.0);
    }
}

TEST(TreePricer, EuropeanBinomialConvergesToBlackScholes)
{
    const auto option = callOption();
    const double reference = blackScholes(option);
    const double binomialCoarseError = std::abs(
        dr3::lattice::europeanBinomial(option, {16}) - reference);
    const double binomialFineError = std::abs(
        dr3::lattice::europeanBinomial(option, {1024}) - reference);
    EXPECT_LT(binomialFineError, 0.02);
    EXPECT_LT(binomialFineError, binomialCoarseError * 0.20);
}

TEST(TreePricer, EuropeanTrinomialConvergesToBlackScholes)
{
    const auto option = callOption();
    const double reference = blackScholes(option);
    const double trinomialCoarseError = std::abs(
        dr3::lattice::europeanTrinomial(option, {16}) - reference);
    const double trinomialFineError = std::abs(
        dr3::lattice::europeanTrinomial(option, {1024}) - reference);

    EXPECT_LT(trinomialFineError, 0.02);
    EXPECT_LT(trinomialFineError, trinomialCoarseError * 0.20);
}

TEST(TreePricer, PutCallParity)
{
    auto call = callOption();
    call.dividendYield = 0.02;
    auto put = call;
    put.type = OptionType::Put;
    const double parity = call.spot * std::exp(-call.dividendYield * call.maturity)
        - call.strike * std::exp(-call.rate * call.maturity);

    const double binomialDifference = dr3::lattice::europeanBinomial(call, {512})
        - dr3::lattice::europeanBinomial(put, {512});
    const double trinomialDifference = dr3::lattice::europeanTrinomial(call, {512})
        - dr3::lattice::europeanTrinomial(put, {512});
    EXPECT_NEAR(binomialDifference, parity, 1.0e-8);
    EXPECT_NEAR(trinomialDifference, parity, 0.01);
}

TEST(TreePricer, CallPriceWithinNoArbitrageBounds)
{
    auto call = callOption();
    const double callPrice = dr3::lattice::europeanTrinomial(call, {512});

    EXPECT_GE(callPrice, std::max(call.spot * std::exp(-call.dividendYield * call.maturity)
                                  - call.strike * std::exp(-call.rate * call.maturity),
                                  0.0));
    EXPECT_LE(callPrice, call.spot * std::exp(-call.dividendYield * call.maturity));
}

TEST(TreePricer, PutPriceWithinNoArbitrageBounds)
{
    auto put = callOption();
    put.type = OptionType::Put;
    const double putPrice = dr3::lattice::europeanTrinomial(put, {512});

    EXPECT_GE(putPrice, std::max(put.strike * std::exp(-put.rate * put.maturity)
                                 - put.spot * std::exp(-put.dividendYield * put.maturity),
                                 0.0));
    EXPECT_LE(putPrice, put.strike * std::exp(-put.rate * put.maturity));
}

TEST(TreePricer, CallIncreasesWithSpot)
{
    auto lowCall = callOption();
    lowCall.spot = 90.0;
    auto highCall = lowCall;
    highCall.spot = 110.0;

    EXPECT_LT(dr3::lattice::europeanTrinomial(lowCall, {256}),
              dr3::lattice::europeanTrinomial(highCall, {256}));
}

TEST(TreePricer, PutDecreasesWithSpot)
{
    auto lowPut = callOption();
    lowPut.spot = 90.0;
    lowPut.type = OptionType::Put;
    auto highPut = lowPut;
    highPut.spot = 110.0;

    EXPECT_GT(dr3::lattice::europeanTrinomial(lowPut, {256}),
              dr3::lattice::europeanTrinomial(highPut, {256}));
}

TEST(TreePricer, AmericanPriceNotBelowEuropeanPrice)
{
    auto put = callOption();
    put.type = OptionType::Put;
    put.spot = 90.0;
    const double european = dr3::lattice::europeanTrinomial(put, {512});
    const double american = dr3::lattice::americanTrinomial(put, {512});
    EXPECT_GE(american, european);
    EXPECT_GT(american, european + 0.01);
}

TEST(TreePricer, AmericanCallWithoutDividendMatchesEuropeanCall)
{
    const auto option = callOption();
    const double european = dr3::lattice::europeanTrinomial(option, {512});
    const double american = dr3::lattice::americanTrinomial(option, {512});
    EXPECT_NEAR(american, european, 1.0e-10);
}

TEST(TreePricer, UpAndOutPriceNotAboveVanillaPrice)
{
    const auto option = callOption();
    const double vanilla = dr3::lattice::europeanTrinomial(option, {512});
    const double barrier = dr3::lattice::upAndOutTrinomial(option, {512}, 120.0);
    EXPECT_GE(barrier, 0.0);
    EXPECT_LE(barrier, vanilla);
    EXPECT_LT(barrier, vanilla - 0.01);
}
