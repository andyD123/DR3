#include "pch.h"

#include "../Vectorisation/VecX/dr3.h"
#include "testNamespace.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

using ForwardADValueAccess = decltype(std::declval<VecxD&>().value());
using ForwardADDerivativeAccess = decltype(std::declval<VecxD&>().derivative());

static_assert(std::is_same<ForwardADValueAccess, const VecXX&>::value,
    "VecD value access must remain read-only");
static_assert(std::is_same<ForwardADDerivativeAccess, const VecXX&>::value,
    "VecD derivative access must remain read-only");
static_assert(!std::is_assignable<ForwardADValueAccess, VecXX>::value,
    "VecD callers must not replace only the value through its accessor");
static_assert(!std::is_assignable<ForwardADDerivativeAccess, VecXX>::value,
    "VecD callers must not replace only the derivative through its accessor");

namespace {

constexpr double kAlgebraTolerance = 1.0e-12;
constexpr double kAnalyticGreekTolerance = 1.0e-7;
constexpr double kFiniteDifferenceTolerance = 1.0e-5;

const std::size_t kLengths[] = {1, 3, 4, 5, 9};

std::vector<double> values(std::size_t n, double offset = 0.0)
{
    std::vector<double> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = offset + 0.4 + 0.13 * static_cast<double>(i);
    }
    return out;
}

template <typename Function>
double central_difference(Function&& f, double x, double h = 1.0e-6)
{
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

double norm_pdf(double x)
{
    constexpr double inv_sqrt_2pi = 0.39894228040143267794;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

double norm_cdf(double x)
{
    return 0.5 * (1.0 + std::erf(x * 0.70710678118654752440));
}

struct BlackScholesReference {
    double price;
    double delta;
    double vega;
    double rho;
    double dPrice_dMaturity;
};

BlackScholesReference black_scholes_reference(
    double spot, double strike, double maturity, double rate, double volatility)
{
    const double root_t = std::sqrt(maturity);
    const double sigma_root_t = volatility * root_t;
    const double d1 = (std::log(spot / strike) +
        (rate + 0.5 * volatility * volatility) * maturity) / sigma_root_t;
    const double d2 = d1 - sigma_root_t;
    const double discount = std::exp(-rate * maturity);
    const double nd1 = norm_cdf(d1);
    const double nd2 = norm_cdf(d2);
    const double pdf_d1 = norm_pdf(d1);
    return {
        spot * nd1 - strike * discount * nd2,
        nd1,
        spot * pdf_d1 * root_t,
        strike * maturity * discount * nd2,
        spot * pdf_d1 * volatility / (2.0 * root_t) +
            rate * strike * discount * nd2,
    };
}

void require_positive_finite(const VecxD& x, const char* name)
{
    auto check = [name](double value) {
        if (!(value > 0.0) || !std::isfinite(value)) {
            throw std::invalid_argument(std::string("Black-Scholes invalid ") + name);
        }
    };
    if (x.isScalar()) {
        check(x.getScalarValue());
        return;
    }
    for (std::size_t i = 0; i < x.size(); ++i) {
        check(x.value()[i]);
    }
}

void require_finite(const VecxD& x, const char* name)
{
    auto check = [name](double value) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(std::string("Black-Scholes invalid ") + name);
        }
    };
    if (x.isScalar()) {
        check(x.getScalarValue());
        return;
    }
    for (std::size_t i = 0; i < x.size(); ++i) {
        check(x.value()[i]);
    }
}

VecxD black_scholes_call_forward_ad(
    const VecxD& spot, const VecxD& strike, const VecxD& maturity,
    const VecxD& rate, const VecxD& volatility)
{
    require_positive_finite(spot, "spot");
    require_positive_finite(strike, "strike");
    require_positive_finite(maturity, "maturity");
    require_finite(rate, "rate");
    require_positive_finite(volatility, "volatility");

    const auto root_t = sqrt(maturity);
    const auto sigma_root_t = volatility * root_t;
    const auto d1 = (log(spot / strike) +
        (rate + 0.5 * volatility * volatility) * maturity) / sigma_root_t;
    const auto d2 = d1 - sigma_root_t;
    return spot * cdfnorm(d1) - strike * exp(-rate * maturity) * cdfnorm(d2);
}

void expect_all(const VecXX& actual, const std::vector<double>& expected, double tolerance)
{
    ASSERT_FALSE(actual.isScalar());
    ASSERT_EQ(static_cast<std::size_t>(actual.size()), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(actual[i], expected[i], tolerance) << "i=" << i;
    }
}

std::vector<double> constant(std::size_t n, double value)
{
    return std::vector<double>(n, value);
}

} // namespace

TEST(ForwardAD, SeedVariableWithOnes)
{
    for (std::size_t n : kLengths) {
        VecXX primal(values(n));
        const auto variable = D(primal);
        expect_all(variable.value(), values(n), 0.0);
        expect_all(variable.derivative(), constant(n, 1.0), 0.0);
    }
}

TEST(ForwardAD, SeedConstantWithZeros)
{
    for (std::size_t n : kLengths) {
        const auto input = values(n);
        const auto fixed = C(VecXX(input));
        expect_all(fixed.value(), input, 0.0);
        expect_all(fixed.derivative(), constant(n, 0.0), 0.0);
    }
}

TEST(ForwardAD, FactoriesAndRvalueConstructorInstantiate)
{
    constexpr int n = 5;
    const auto ones = VecxD::makeDVecOnes(2.5, n);
    const auto zeros = VecxD::makeDVecZero(2.5, n);
    const auto ones_v = VecxD::makeDVecOnesV(2.5, n);
    const auto zeros_v = VecxD::makeDVecZeroV(2.5, n);
    expect_all(ones.derivative(), constant(n, 1.0), 0.0);
    expect_all(zeros.derivative(), constant(n, 0.0), 0.0);
    expect_all(ones_v.derivative(), constant(n, 1.0), 0.0);
    expect_all(zeros_v.derivative(), constant(n, 0.0), 0.0);

    const auto input = values(n);
    const VecxD moved{VecXX(input)};
    expect_all(moved.value(), input, 0.0);
    expect_all(moved.derivative(), constant(n, 0.0), 0.0);
}

TEST(ForwardAD, ScalarMinusVariableHasNegativeDerivative)
{
    const auto input = values(5);
    const auto result = 3.0 - D(VecXX(input));
    std::vector<double> expected_value(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        expected_value[i] = 3.0 - input[i];
    }
    expect_all(result.value(), expected_value, kAlgebraTolerance);
    expect_all(result.derivative(), constant(input.size(), -1.0), kAlgebraTolerance);
}

TEST(ForwardAD, AdditionRule)
{
    const auto x_values = values(9, 0.2);
    const auto y_values = values(9, 1.1);
    const auto result = D(VecXX(x_values)) + 2.0 * D(VecXX(y_values));
    expect_all(result.derivative(), constant(x_values.size(), 3.0), kAlgebraTolerance);
}

TEST(ForwardAD, SubtractionRule)
{
    const auto input = values(5);
    const auto result = 4.0 * D(VecXX(input)) - D(VecXX(input));
    expect_all(result.derivative(), constant(input.size(), 3.0), kAlgebraTolerance);
}

TEST(ForwardAD, ProductRule)
{
    const auto x_values = values(5, 0.5);
    const auto y_values = values(5, 1.3);
    const auto result = D(VecXX(x_values)) * D(VecXX(y_values));
    std::vector<double> expected(x_values.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        expected[i] = x_values[i] + y_values[i];
    }
    expect_all(result.derivative(), expected, kAlgebraTolerance);
}

TEST(ForwardAD, QuotientRule)
{
    const auto x_values = values(5, 1.0);
    const auto y_values = values(5, 2.0);
    const auto result = D(VecXX(x_values)) / D(VecXX(y_values));
    std::vector<double> expected(x_values.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        expected[i] = (y_values[i] - x_values[i]) / (y_values[i] * y_values[i]);
    }
    expect_all(result.derivative(), expected, kAlgebraTolerance);
}

TEST(ForwardAD, ExpLogChainRule)
{
    const auto input = values(9, 0.7);
    const auto result = exp(log(D(VecXX(input))));
    expect_all(result.value(), input, kAlgebraTolerance);
    expect_all(result.derivative(), constant(input.size(), 1.0), kAlgebraTolerance);
}

TEST(ForwardAD, SqrtRuleForPositiveInputs)
{
    const auto input = values(5, 0.8);
    const auto result = sqrt(D(VecXX(input)));
    std::vector<double> expected(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        expected[i] = 0.5 / std::sqrt(input[i]);
    }
    expect_all(result.derivative(), expected, kAlgebraTolerance);
}

TEST(ForwardAD, CompositeExpressionMatchesCentralDifference)
{
    const auto input = values(9, 0.9);
    const auto x = D(VecXX(input));
    const auto result = exp(log(x * x + 1.7)) / sqrt(x + 0.8);
    for (std::size_t i = 0; i < input.size(); ++i) {
        const auto scalar = [](double v) {
            return std::exp(std::log(v * v + 1.7)) / std::sqrt(v + 0.8);
        };
        EXPECT_NEAR(result.derivative()[i], central_difference(scalar, input[i]),
            kFiniteDifferenceTolerance) << "i=" << i;
    }
}

TEST(ForwardAD, ScalarBroadcast)
{
    const auto input = values(5);
    const VecxD scalar(2.0, 3.0);
    const auto result = D(VecXX(input)) + scalar;
    expect_all(result.derivative(), constant(input.size(), 4.0), kAlgebraTolerance);
}

TEST(ForwardAD, NonMultipleOfSimdWidth)
{
    const auto input = values(5, 0.6);
    const auto result = exp(D(VecXX(input)) * 1.5 - 0.2);
    ASSERT_EQ(result.size(), input.size());
    ASSERT_EQ(result.derivative().size(), static_cast<int>(input.size()));
    for (std::size_t i = 0; i < input.size(); ++i) {
        const double expected_value = std::exp(1.5 * input[i] - 0.2);
        EXPECT_NEAR(result.value()[i], expected_value, kAlgebraTolerance);
        EXPECT_NEAR(result.derivative()[i], 1.5 * expected_value, kAlgebraTolerance);
    }
}

TEST(ForwardAD, RejectsMismatchedLogicalSizes)
{
    const VecXX values3(std::vector<double>{1.0, 2.0, 3.0});
    const VecXX values2(std::vector<double>{1.0, 2.0});
    EXPECT_THROW((VecxD(values3, values2)), std::invalid_argument);
    EXPECT_THROW((VecxD(values3, VecXX(0.0))), std::invalid_argument);
    const auto x = D(values3);
    const auto y = D(values2);
    EXPECT_THROW((void)(x + y), std::invalid_argument);
    EXPECT_THROW((void)(x - y), std::invalid_argument);
    EXPECT_THROW((void)(x * y), std::invalid_argument);
    EXPECT_THROW((void)(x / y), std::invalid_argument);
}

TEST(ForwardAD, ReplacementAPIsPreserveShapeInvariant)
{
    const std::vector<double> original_values = {1.0, 2.0, 3.0};
    const std::vector<double> original_derivatives = {0.1, 0.2, 0.3};
    VecxD value{VecXX(original_values), VecXX(original_derivatives)};

    EXPECT_THROW(value.replaceValue(VecXX(std::vector<double>{4.0, 5.0})),
        std::invalid_argument);
    EXPECT_THROW(value.replaceDerivative(VecXX(std::vector<double>{0.4, 0.5})),
        std::invalid_argument);
    EXPECT_THROW(value.replaceDerivative(VecXX(1.0)), std::invalid_argument);
    EXPECT_THROW(value.replaceValueAndDerivative(
        VecXX(std::vector<double>{4.0, 5.0}),
        VecXX(std::vector<double>{0.4, 0.5, 0.6})), std::invalid_argument);

    expect_all(value.value(), original_values, 0.0);
    expect_all(value.derivative(), original_derivatives, 0.0);

    const std::vector<double> replacement_values = {4.0, 5.0, 6.0};
    const std::vector<double> replacement_derivatives = {0.4, 0.5, 0.6};
    value.replaceValue(VecXX(replacement_values));
    value.replaceDerivative(VecXX(replacement_derivatives));
    expect_all(value.value(), replacement_values, 0.0);
    expect_all(value.derivative(), replacement_derivatives, 0.0);

    const std::vector<double> resized_values = {7.0, 8.0};
    const std::vector<double> resized_derivatives = {0.7, 0.8};
    value.replaceValueAndDerivative(
        VecXX(resized_values), VecXX(resized_derivatives));
    expect_all(value.value(), resized_values, 0.0);
    expect_all(value.derivative(), resized_derivatives, 0.0);

    VecxD assigned{VecXX(original_values), VecXX(original_derivatives)};
    assigned = value;
    expect_all(assigned.value(), resized_values, 0.0);
    expect_all(assigned.derivative(), resized_derivatives, 0.0);

    VecxD scalar(2.0, 3.0);
    EXPECT_THROW(scalar.replaceValue(VecXX(resized_values)), std::invalid_argument);
    EXPECT_THROW(scalar.replaceDerivative(VecXX(resized_derivatives)),
        std::invalid_argument);
    EXPECT_DOUBLE_EQ(scalar.getScalarValue(), 2.0);
    EXPECT_DOUBLE_EQ(scalar.getScalarDeriv(), 3.0);
}

TEST(ForwardAD, BlackScholesDeltaMatchesAnalytic)
{
    const std::vector<double> spots = {80.0, 95.0, 100.0, 110.0, 125.0};
    const auto result = black_scholes_call_forward_ad(
        D(VecXX(spots)), C(VecXX(constant(spots.size(), 100.0))),
        C(VecXX(constant(spots.size(), 1.25))), C(VecXX(constant(spots.size(), 0.04))),
        C(VecXX(constant(spots.size(), 0.22))));
    for (std::size_t i = 0; i < spots.size(); ++i) {
        const auto ref = black_scholes_reference(spots[i], 100.0, 1.25, 0.04, 0.22);
        EXPECT_NEAR(result.derivative()[i], ref.delta, kAnalyticGreekTolerance);
    }
}

TEST(ForwardAD, BlackScholesVegaMatchesAnalytic)
{
    const std::vector<double> vols = {0.10, 0.15, 0.20, 0.30, 0.60};
    const auto result = black_scholes_call_forward_ad(
        C(VecXX(constant(vols.size(), 102.0))), C(VecXX(constant(vols.size(), 100.0))),
        C(VecXX(constant(vols.size(), 0.8))), C(VecXX(constant(vols.size(), 0.03))),
        D(VecXX(vols)));
    for (std::size_t i = 0; i < vols.size(); ++i) {
        const auto ref = black_scholes_reference(102.0, 100.0, 0.8, 0.03, vols[i]);
        EXPECT_NEAR(result.derivative()[i], ref.vega, kAnalyticGreekTolerance);
    }
}

TEST(ForwardAD, BlackScholesRateSensitivityMatchesAnalytic)
{
    const std::vector<double> rates = {-0.01, 0.0, 0.02, 0.05, 0.10};
    const auto result = black_scholes_call_forward_ad(
        C(VecXX(constant(rates.size(), 100.0))), C(VecXX(constant(rates.size(), 105.0))),
        C(VecXX(constant(rates.size(), 1.5))), D(VecXX(rates)),
        C(VecXX(constant(rates.size(), 0.25))));
    for (std::size_t i = 0; i < rates.size(); ++i) {
        const auto ref = black_scholes_reference(100.0, 105.0, 1.5, rates[i], 0.25);
        EXPECT_NEAR(result.derivative()[i], ref.rho, kAnalyticGreekTolerance);
    }
}

TEST(ForwardAD, BlackScholesMaturitySensitivityMatchesAnalytic)
{
    const std::vector<double> maturities = {0.1, 0.25, 0.75, 1.0, 3.0};
    const auto result = black_scholes_call_forward_ad(
        C(VecXX(constant(maturities.size(), 98.0))),
        C(VecXX(constant(maturities.size(), 100.0))), D(VecXX(maturities)),
        C(VecXX(constant(maturities.size(), 0.05))),
        C(VecXX(constant(maturities.size(), 0.20))));
    for (std::size_t i = 0; i < maturities.size(); ++i) {
        const auto ref = black_scholes_reference(98.0, 100.0, maturities[i], 0.05, 0.20);
        EXPECT_NEAR(result.derivative()[i], ref.dPrice_dMaturity,
            kAnalyticGreekTolerance);
    }
}

TEST(ForwardAD, BlackScholesSensitivitiesMatchCentralDifference)
{
    constexpr double s = 101.0;
    constexpr double k = 100.0;
    constexpr double t = 1.2;
    constexpr double r = 0.04;
    constexpr double v = 0.23;
    const auto scalar_price = [](double spot, double strike, double maturity,
                                  double rate, double volatility) {
        return black_scholes_reference(spot, strike, maturity, rate, volatility).price;
    };
    auto evaluate = [](const VecxD& spot, const VecxD& strike, const VecxD& maturity,
                       const VecxD& rate, const VecxD& volatility) {
        return black_scholes_call_forward_ad(spot, strike, maturity, rate, volatility)
            .getScalarDeriv();
    };
    EXPECT_NEAR(evaluate(VecxD(s, 1.0), VecxD(k), VecxD(t), VecxD(r), VecxD(v)),
        central_difference([&](double x) { return scalar_price(x, k, t, r, v); }, s),
        kFiniteDifferenceTolerance);
    EXPECT_NEAR(evaluate(VecxD(s), VecxD(k), VecxD(t), VecxD(r), VecxD(v, 1.0)),
        central_difference([&](double x) { return scalar_price(s, k, t, r, x); }, v),
        kFiniteDifferenceTolerance);
    EXPECT_NEAR(evaluate(VecxD(s), VecxD(k), VecxD(t), VecxD(r, 1.0), VecxD(v)),
        central_difference([&](double x) { return scalar_price(s, k, t, x, v); }, r),
        kFiniteDifferenceTolerance);
    EXPECT_NEAR(evaluate(VecxD(s), VecxD(k), VecxD(t, 1.0), VecxD(r), VecxD(v)),
        central_difference([&](double x) { return scalar_price(s, k, x, r, v); }, t),
        kFiniteDifferenceTolerance);
}

TEST(ForwardAD, BlackScholesRejectsInvalidInputs)
{
    EXPECT_THROW(black_scholes_call_forward_ad(
        VecxD(100.0), VecxD(0.0), VecxD(1.0), VecxD(0.05), VecxD(0.2)),
        std::invalid_argument);
    EXPECT_THROW(black_scholes_call_forward_ad(
        VecxD(100.0), VecxD(100.0), VecxD(0.0), VecxD(0.05), VecxD(0.2)),
        std::invalid_argument);
    EXPECT_THROW(black_scholes_call_forward_ad(
        VecxD(100.0), VecxD(100.0), VecxD(1.0), VecxD(0.05), VecxD(-0.2)),
        std::invalid_argument);
    EXPECT_THROW(black_scholes_call_forward_ad(
        VecxD(100.0), VecxD(100.0), VecxD(1.0),
        VecxD(std::numeric_limits<double>::infinity()), VecxD(0.2)),
        std::invalid_argument);
}
