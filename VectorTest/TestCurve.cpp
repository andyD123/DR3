#include "pch.h"

#include "../Vectorisation/Curves/curve.h"
#include "testNamespace.h"

#include <cmath>
#include <limits>
#include <vector>

namespace {

using DRC::Curves::Curve;
using DRC::Curves::Extrapolation;
using DRC::Curves::Interpolation;

Curve<double> scalar_curve(
    Interpolation interpolation = Interpolation::Linear,
    Extrapolation left = Extrapolation::Flat,
    Extrapolation right = Extrapolation::Flat)
{
    return Curve<double>({0.5, 1.0, 2.0}, {0.01, 0.02, 0.04},
        interpolation, left, right);
}

} // namespace

TEST(Curve, RejectsEmptyPillars)
{
    EXPECT_THROW((Curve<double>({}, {}, Interpolation::Linear,
        Extrapolation::Flat, Extrapolation::Flat)), std::invalid_argument);
}

TEST(Curve, RejectsMismatchedPillarAndValueCounts)
{
    EXPECT_THROW((Curve<double>({1.0, 2.0}, {0.01}, Interpolation::Linear,
        Extrapolation::Flat, Extrapolation::Flat)), std::invalid_argument);
}

TEST(Curve, RejectsUnsortedAndDuplicatePillars)
{
    EXPECT_THROW((Curve<double>({1.0, 0.5}, {0.01, 0.02}, Interpolation::Linear,
        Extrapolation::Flat, Extrapolation::Flat)), std::invalid_argument);
    EXPECT_THROW((Curve<double>({1.0, 1.0}, {0.01, 0.02}, Interpolation::Linear,
        Extrapolation::Flat, Extrapolation::Flat)), std::invalid_argument);
}

TEST(Curve, RejectsNonFinitePillarsAndValues)
{
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW((Curve<double>({1.0, inf}, {0.01, 0.02}, Interpolation::Linear,
        Extrapolation::Flat, Extrapolation::Flat)), std::invalid_argument);
    EXPECT_THROW((Curve<double>({1.0, 2.0}, {0.01, inf}, Interpolation::Linear,
        Extrapolation::Flat, Extrapolation::Flat)), std::invalid_argument);
}

TEST(Curve, ExactFirstInteriorAndFinalPillars)
{
    const auto curve = scalar_curve();
    EXPECT_DOUBLE_EQ(curve.value_at(0.5), 0.01);
    EXPECT_DOUBLE_EQ(curve.value_at(1.0), 0.02);
    EXPECT_DOUBLE_EQ(curve.value_at(2.0), 0.04);
}

TEST(Curve, LinearMidpointAndKnownFraction)
{
    const auto curve = scalar_curve();
    EXPECT_DOUBLE_EQ(curve.value_at(1.5), 0.03);
    EXPECT_DOUBLE_EQ(curve.value_at(1.25), 0.025);
}

TEST(Curve, FlatInterpolation)
{
    const auto curve = scalar_curve(Interpolation::Flat);
    EXPECT_DOUBLE_EQ(curve.value_at(0.75), 0.01);
    EXPECT_DOUBLE_EQ(curve.value_at(1.75), 0.02);
}

TEST(Curve, ExplicitFlatAndThrowExtrapolation)
{
    const auto flat = scalar_curve();
    EXPECT_DOUBLE_EQ(flat.value_at(-1.0), 0.01);
    EXPECT_DOUBLE_EQ(flat.value_at(5.0), 0.04);
    const auto throwing = scalar_curve(
        Interpolation::Linear, Extrapolation::Throw, Extrapolation::Throw);
    EXPECT_THROW(throwing.value_at(0.0), std::out_of_range);
    EXPECT_THROW(throwing.value_at(3.0), std::out_of_range);
}

TEST(Curve, ResetReplacesAllDerivedState)
{
    auto curve = scalar_curve();
    curve.reset({2.0, 4.0}, {0.10, 0.30});
    EXPECT_DOUBLE_EQ(curve.value_at(3.0), 0.20);
    EXPECT_DOUBLE_EQ(curve.value_at(1.0), 0.10);
}

TEST(Curve, BulkSortedEvaluationMatchesScalarEvaluation)
{
    const auto curve = scalar_curve();
    const double queries[] = {0.0, 0.5, 0.75, 1.25, 2.0, 3.0};
    double output[6] = {};
    curve.evaluate_sorted(queries, 6, output, 6);
    for (std::size_t i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(output[i], curve.value_at(queries[i]));
    }
    EXPECT_THROW(curve.evaluate_sorted(queries, 6, output, 5), std::invalid_argument);
    const double unsorted[] = {1.0, 0.5};
    EXPECT_THROW(curve.evaluate_sorted(unsorted, 2, output, 6), std::invalid_argument);
    EXPECT_NO_THROW(curve.evaluate_sorted(nullptr, 0, nullptr, 0));
}

TEST(Curve, VectorValuesMatchScalarLaneByLaneAndTail)
{
    const std::vector<double> first = {0.01, 0.02, 0.03, 0.04, 0.05};
    const std::vector<double> second = {0.03, 0.04, 0.05, 0.06, 0.07};
    Curve<VecXX> curve({1.0, 2.0}, {VecXX(first), VecXX(second)},
        Interpolation::Linear, Extrapolation::Flat, Extrapolation::Flat);
    const VecXX result = curve.value_at(1.25);
    ASSERT_EQ(result.size(), 5);
    for (std::size_t lane = 0; lane < first.size(); ++lane) {
        EXPECT_NEAR(result[lane], first[lane] + 0.25 * (second[lane] - first[lane]), 1e-15);
    }
}

TEST(Curve, RejectsMismatchedVectorValueShapes)
{
    EXPECT_THROW((Curve<VecXX>({1.0, 2.0},
        {VecXX(std::vector<double>{1.0, 2.0}), VecXX(std::vector<double>{1.0})},
        Interpolation::Linear, Extrapolation::Flat, Extrapolation::Flat)),
        std::invalid_argument);
}

TEST(CurveSensitivity, LinearNodeWeights)
{
    Curve<VecxD> lower_active({1.0, 2.0, 3.0},
        {VecxD(0.01, 1.0), VecxD(0.02), VecxD(0.04)},
        Interpolation::Linear, Extrapolation::Flat, Extrapolation::Flat);
    Curve<VecxD> upper_active({1.0, 2.0, 3.0},
        {VecxD(0.01), VecxD(0.02, 1.0), VecxD(0.04)},
        Interpolation::Linear, Extrapolation::Flat, Extrapolation::Flat);
    Curve<VecxD> non_adjacent({1.0, 2.0, 3.0},
        {VecxD(0.01), VecxD(0.02), VecxD(0.04, 1.0)},
        Interpolation::Linear, Extrapolation::Flat, Extrapolation::Flat);
    EXPECT_NEAR(lower_active.value_at(1.25).getScalarDeriv(), 0.75, 1e-15);
    EXPECT_NEAR(upper_active.value_at(1.25).getScalarDeriv(), 0.25, 1e-15);
    EXPECT_NEAR(non_adjacent.value_at(1.25).getScalarDeriv(), 0.0, 1e-15);
    EXPECT_NEAR(lower_active.value_at(1.25).getScalarDeriv() +
        upper_active.value_at(1.25).getScalarDeriv(), 1.0, 1e-15);
}

TEST(CurveSensitivity, DiscountFactorSensitivityMatchesFormulaAndCentralBump)
{
    constexpr double query = 1.25;
    Curve<VecxD> curve({1.0, 2.0},
        {VecxD(0.02, 1.0), VecxD(0.04)},
        Interpolation::Linear, Extrapolation::Flat, Extrapolation::Flat);
    const auto discount = exp(-curve.value_at(query) * query);
    const double zero_rate = 0.02 + 0.25 * (0.04 - 0.02);
    const double expected = -query * std::exp(-zero_rate * query) * 0.75;
    EXPECT_NEAR(discount.getScalarDeriv(), expected, 1e-14);
    const auto bumped = [query](double lower) {
        const double z = lower + 0.25 * (0.04 - lower);
        return std::exp(-z * query);
    };
    const double h = 1e-6;
    const double central = (bumped(0.02 + h) - bumped(0.02 - h)) / (2.0 * h);
    EXPECT_NEAR(discount.getScalarDeriv(), central, 1e-9);
}

TEST(CurveSensitivity, SimdBlockSeedsSeveralPillarsAtOnce)
{
    const std::vector<double> primal0(5, 0.01);
    const std::vector<double> primal1(5, 0.02);
    const std::vector<double> primal2(5, 0.04);
    const std::vector<double> d0 = {1.0, 0.0, 0.0, 0.0, 0.0};
    const std::vector<double> d1 = {0.0, 1.0, 0.0, 0.0, 0.0};
    const std::vector<double> d2 = {0.0, 0.0, 1.0, 0.0, 0.0};
    Curve<VecxD> curve({1.0, 2.0, 3.0},
        {VecxD(VecXX(primal0), VecXX(d0)), VecxD(VecXX(primal1), VecXX(d1)),
         VecxD(VecXX(primal2), VecXX(d2))},
        Interpolation::Linear, Extrapolation::Flat, Extrapolation::Flat);
    const auto result = curve.value_at(1.25);
    EXPECT_NEAR(result.derivative()[0], 0.75, 1e-15);
    EXPECT_NEAR(result.derivative()[1], 0.25, 1e-15);
    EXPECT_NEAR(result.derivative()[2], 0.0, 1e-15);
    EXPECT_NEAR(result.derivative()[3], 0.0, 1e-15);
    EXPECT_NEAR(result.derivative()[4], 0.0, 1e-15);
}
