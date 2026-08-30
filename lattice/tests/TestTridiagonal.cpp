#include "tridiagonal.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace
{

std::vector<double> multiply(const dr3::lattice::TridiagonalSystem& system,
                             const std::vector<double>& values)
{
    std::vector<double> result(values.size());
    for (std::size_t row = 0; row < values.size(); ++row)
    {
        result[row] = system.main[row] * values[row];
        if (row > 0)
        {
            result[row] += system.lower[row - 1] * values[row - 1];
        }
        if (row + 1 < values.size())
        {
            result[row] += system.upper[row] * values[row + 1];
        }
    }
    return result;
}

std::vector<double> solve(const dr3::lattice::TridiagonalSystem& system,
                          const std::vector<double>& rightHandSide)
{
    const dr3::lattice::ThomasFactorization factorization{system};
    std::vector<double> result(system.size());
    std::vector<double> workspace(system.size());
    factorization.solve(rightHandSide, result, workspace);
    return result;
}

void expectVectorNear(const std::vector<double>& actual,
                      const std::vector<double>& expected,
                      double tolerance = 1.0e-12)
{
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t index = 0; index < actual.size(); ++index)
    {
        EXPECT_NEAR(actual[index], expected[index], tolerance) << "index " << index;
    }
}

} // namespace

TEST(Tridiagonal, SolvesKnownThreeByThreeSystem)
{
    const dr3::lattice::TridiagonalSystem system{{-1.0, -1.0},
                                                  {2.0, 2.0, 2.0},
                                                  {-1.0, -1.0}};
    const std::vector<double> expected{1.0, 2.0, 3.0};
    expectVectorNear(solve(system, multiply(system, expected)), expected);
}

TEST(Tridiagonal, SolvesKnownFiveByFiveSystem)
{
    const dr3::lattice::TridiagonalSystem system{{-1.0, -1.0, -1.0, -1.0},
                                                  {4.0, 4.0, 4.0, 4.0, 4.0},
                                                  {-1.0, -1.0, -1.0, -1.0}};
    const std::vector<double> expected{1.0, -2.0, 3.0, -4.0, 5.0};
    expectVectorNear(solve(system, multiply(system, expected)), expected);
}

TEST(Tridiagonal, SingleElementSystem)
{
    const dr3::lattice::TridiagonalSystem system{{}, {4.0}, {}};
    expectVectorNear(solve(system, {10.0}), {2.5});
}

TEST(Tridiagonal, RejectsDimensionMismatch)
{
    EXPECT_THROW((dr3::lattice::ThomasFactorization{{{}, {1.0, 2.0}, {1.0}}}),
                 std::invalid_argument);
    EXPECT_THROW((dr3::lattice::ThomasFactorization{{{1.0}, {1.0}, {}}}),
                 std::invalid_argument);

    const dr3::lattice::ThomasFactorization factorization{{{}, {2.0}, {}}};
    std::vector<double> rhs(2), output(1), workspace(1);
    EXPECT_THROW(factorization.solve(rhs, output, workspace), std::invalid_argument);
}

TEST(Tridiagonal, RejectsNearZeroPivot)
{
    EXPECT_THROW((dr3::lattice::ThomasFactorization{{{}, {1.0e-16}, {}}}),
                 std::domain_error);
    EXPECT_THROW((dr3::lattice::ThomasFactorization{{{1.0}, {1.0, 1.0}, {1.0}}}),
                 std::domain_error);
}

TEST(Tridiagonal, RandomDiagonallyDominantResidualIsSmall)
{
    constexpr std::size_t size = 97;
    std::mt19937_64 generator(0xD3C0FFEEu);
    std::uniform_real_distribution<double> distribution(-0.75, 0.75);
    dr3::lattice::TridiagonalSystem system;
    system.lower.resize(size - 1);
    system.main.resize(size);
    system.upper.resize(size - 1);
    std::vector<double> rightHandSide(size);
    for (double& value : system.lower) value = distribution(generator);
    for (double& value : system.upper) value = distribution(generator);
    for (double& value : rightHandSide) value = distribution(generator);
    for (std::size_t row = 0; row < size; ++row)
    {
        const double offDiagonal = (row > 0 ? std::abs(system.lower[row - 1]) : 0.0)
            + (row + 1 < size ? std::abs(system.upper[row]) : 0.0);
        system.main[row] = offDiagonal + 1.0;
    }

    const auto result = solve(system, rightHandSide);
    const auto reconstructed = multiply(system, result);
    double maximumResidual = 0.0;
    for (std::size_t row = 0; row < size; ++row)
    {
        maximumResidual = std::max(maximumResidual,
                                   std::abs(reconstructed[row] - rightHandSide[row]));
    }
    EXPECT_LE(maximumResidual, 1.0e-11);
}

TEST(Tridiagonal, FactorizationCanBeReused)
{
    const dr3::lattice::TridiagonalSystem system{{-1.0, -1.0},
                                                  {3.0, 3.0, 3.0},
                                                  {-1.0, -1.0}};
    const dr3::lattice::ThomasFactorization factorization{system};
    std::vector<double> output(3), workspace(3);
    const std::vector<double> firstExpected{1.0, 2.0, 3.0};
    const std::vector<double> secondExpected{-3.0, 0.5, 4.0};

    factorization.solve(multiply(system, firstExpected), output, workspace);
    expectVectorNear(output, firstExpected);
    factorization.solve(multiply(system, secondExpected), output, workspace);
    expectVectorNear(output, secondExpected);
}

TEST(Tridiagonal, MultipleRightHandSidesMatchIndependentSolves)
{
    const dr3::lattice::TridiagonalSystem system{{0.5, -0.25, 0.75},
                                                  {3.0, 4.0, 5.0, 6.0},
                                                  {-0.5, 0.25, -0.75}};
    const dr3::lattice::ThomasFactorization reused{system};
    std::vector<double> output(4), workspace(4);
    const std::vector<std::vector<double>> rightHandSides{
        {1.0, 2.0, 3.0, 4.0}, {-4.0, 3.0, -2.0, 1.0}, {0.25, 0.5, 0.75, 1.0}};

    for (const auto& rightHandSide : rightHandSides)
    {
        reused.solve(rightHandSide, output, workspace);
        expectVectorNear(output, solve(system, rightHandSide));
    }
}

TEST(Tridiagonal, InputMatrixRemainsUnchanged)
{
    const dr3::lattice::TridiagonalSystem original{{-1.0, -2.0},
                                                    {4.0, 5.0, 6.0},
                                                    {0.5, 1.5}};
    auto system = original;
    const dr3::lattice::ThomasFactorization factorization{system};

    EXPECT_EQ(system.lower, original.lower);
    EXPECT_EQ(system.main, original.main);
    EXPECT_EQ(system.upper, original.upper);
}
