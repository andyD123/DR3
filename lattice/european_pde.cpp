#include "european_pde.h"

#include "pde_workspace.h"
#include "tree_pricer_detail.h"
#include "tridiagonal.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace dr3::lattice
{
namespace
{

struct SpatialOperator
{
    std::vector<double> lower;
    std::vector<double> main;
    std::vector<double> upper;
};

void validateConfig(const VanillaOption& option, const PdeConfig& config)
{
    detail::validate(option, {1});
    if (config.spaceSteps < 3)
    {
        throw std::invalid_argument("PDE spaceSteps must be at least three");
    }
    if (config.timeSteps < 1)
    {
        throw std::invalid_argument("PDE timeSteps must be at least one");
    }
    if (!std::isfinite(config.minimumSpot) || !std::isfinite(config.maximumSpot)
        || config.maximumSpot <= config.minimumSpot)
    {
        throw std::invalid_argument("PDE spot bounds must be finite and increasing");
    }
    if (option.spot < config.minimumSpot || option.spot > config.maximumSpot)
    {
        throw std::invalid_argument("pricing spot must lie inside the PDE grid");
    }
    switch (config.scheme)
    {
    case PdeScheme::Explicit:
    case PdeScheme::BackwardEuler:
    case PdeScheme::CrankNicolson:
    case PdeScheme::CrankNicolsonRannacher:
        break;
    default:
        throw std::invalid_argument("PDE scheme is invalid");
    }
}

SpatialOperator makeSpatialOperator(const VanillaOption& option, const Grid1D& grid)
{
    SpatialOperator result{
        std::vector<double>(grid.nodeCount()),
        std::vector<double>(grid.nodeCount()),
        std::vector<double>(grid.nodeCount())};
    const double inverseSpacing = 1.0 / grid.spacing();
    const double inverseSpacingSquared = inverseSpacing * inverseSpacing;
    const double variance = option.volatility * option.volatility;
    const double carry = option.rate - option.dividendYield;
    for (std::size_t index = grid.interiorBegin(); index < grid.interiorEnd(); ++index)
    {
        const double spot = grid[index];
        const double diffusion = 0.5 * variance * spot * spot * inverseSpacingSquared;
        const double convection = 0.5 * carry * spot * inverseSpacing;
        result.lower[index] = diffusion - convection;
        result.main[index] = -2.0 * diffusion - option.rate;
        result.upper[index] = diffusion + convection;
    }
    return result;
}

TridiagonalSystem makeImplicitSystem(const SpatialOperator& op,
                                     std::size_t nodeCount,
                                     double theta,
                                     double timeStep)
{
    const std::size_t interiorCount = nodeCount - 2;
    TridiagonalSystem system;
    system.lower.resize(interiorCount - 1);
    system.main.resize(interiorCount);
    system.upper.resize(interiorCount - 1);
    for (std::size_t row = 0; row < interiorCount; ++row)
    {
        const std::size_t index = row + 1;
        system.main[row] = 1.0 - theta * timeStep * op.main[index];
        if (row > 0)
        {
            system.lower[row - 1] = -theta * timeStep * op.lower[index];
        }
        if (row + 1 < interiorCount)
        {
            system.upper[row] = -theta * timeStep * op.upper[index];
        }
    }
    return system;
}

void applyBoundaries(std::vector<double>& values,
                     const VanillaOption& option,
                     const Grid1D& grid,
                     double timeToMaturity)
{
    const auto boundary = europeanBoundaryValues(
        option, grid.minimum(), grid.maximum(), timeToMaturity);
    values.front() = boundary.lower;
    values.back() = boundary.upper;
}

void ensureFiniteAndNonnegative(std::vector<double>& values)
{
    for (double& value : values)
    {
        if (!std::isfinite(value))
        {
            throw std::runtime_error("PDE produced a non-finite value");
        }
        value = std::max(value, 0.0);
    }
}

void validateExplicitStability(const SpatialOperator& op,
                               const Grid1D& grid,
                               double timeStep)
{
    constexpr double tolerance = 1.0e-14;
    for (std::size_t index = grid.interiorBegin(); index < grid.interiorEnd(); ++index)
    {
        const double lowerWeight = timeStep * op.lower[index];
        const double mainWeight = 1.0 + timeStep * op.main[index];
        const double upperWeight = timeStep * op.upper[index];
        if (!std::isfinite(lowerWeight) || !std::isfinite(mainWeight)
            || !std::isfinite(upperWeight) || lowerWeight < -tolerance
            || mainWeight < -tolerance || upperWeight < -tolerance)
        {
            throw std::invalid_argument("explicit PDE configuration is unstable");
        }
    }
}

void explicitStep(const VanillaOption& option,
                  const Grid1D& grid,
                  const SpatialOperator& op,
                  double timeStep,
                  double newTime,
                  PdeWorkspace& workspace)
{
    applyBoundaries(workspace.next, option, grid, newTime);
    for (std::size_t index = grid.interiorBegin(); index < grid.interiorEnd(); ++index)
    {
        workspace.next[index] = workspace.previous[index] + timeStep
            * (op.lower[index] * workspace.previous[index - 1]
               + op.main[index] * workspace.previous[index]
               + op.upper[index] * workspace.previous[index + 1]);
    }
    ensureFiniteAndNonnegative(workspace.next);
    workspace.previous.swap(workspace.next);
}

void thetaStep(const VanillaOption& option,
               const Grid1D& grid,
               const SpatialOperator& op,
               const ThomasFactorization& factorization,
               double theta,
               double timeStep,
               double newTime,
               PdeWorkspace& workspace)
{
    applyBoundaries(workspace.next, option, grid, newTime);
    const double explicitWeight = (1.0 - theta) * timeStep;
    for (std::size_t index = grid.interiorBegin(); index < grid.interiorEnd(); ++index)
    {
        workspace.rightHandSide[index] = workspace.previous[index] + explicitWeight
            * (op.lower[index] * workspace.previous[index - 1]
               + op.main[index] * workspace.previous[index]
               + op.upper[index] * workspace.previous[index + 1]);
    }
    workspace.rightHandSide[1] += theta * timeStep * op.lower[1] * workspace.next.front();
    const std::size_t finalInterior = grid.nodeCount() - 2;
    workspace.rightHandSide[finalInterior] += theta * timeStep
        * op.upper[finalInterior] * workspace.next.back();

    const std::size_t interiorCount = grid.nodeCount() - 2;
    factorization.solve(
        ConstDoubleSpan(workspace.rightHandSide.data() + 1, interiorCount),
        DoubleSpan(workspace.next.data() + 1, interiorCount),
        DoubleSpan(workspace.solver.data(), interiorCount));
    ensureFiniteAndNonnegative(workspace.next);
    workspace.previous.swap(workspace.next);
}

} // namespace

std::vector<double> europeanTerminalPayoff(const VanillaOption& option, const Grid1D& grid)
{
    detail::validate(option, {1});
    std::vector<double> payoff(grid.nodeCount());
    for (std::size_t index = 0; index < grid.nodeCount(); ++index)
    {
        payoff[index] = detail::intrinsic(option, grid[index]);
    }
    return payoff;
}

PdeBoundaryValues europeanBoundaryValues(const VanillaOption& option,
                                         double minimumSpot,
                                         double maximumSpot,
                                         double timeToMaturity)
{
    detail::validate(option, {1});
    if (!std::isfinite(minimumSpot) || !std::isfinite(maximumSpot)
        || maximumSpot <= minimumSpot || !std::isfinite(timeToMaturity)
        || timeToMaturity < 0.0)
    {
        throw std::invalid_argument("PDE boundary inputs are invalid");
    }
    const double discountedStrike = option.strike * std::exp(-option.rate * timeToMaturity);
    const double discountedMinimum = minimumSpot
        * std::exp(-option.dividendYield * timeToMaturity);
    const double discountedMaximum = maximumSpot
        * std::exp(-option.dividendYield * timeToMaturity);
    PdeBoundaryValues values;
    if (option.type == OptionType::Call)
    {
        values = {std::max(discountedMinimum - discountedStrike, 0.0),
                  std::max(discountedMaximum - discountedStrike, 0.0)};
    }
    else
    {
        values = {std::max(discountedStrike - discountedMinimum, 0.0),
                  std::max(discountedStrike - discountedMaximum, 0.0)};
    }
    values.validate();
    return values;
}

PdeResult europeanPdePrice(const VanillaOption& option, const PdeConfig& config)
{
    validateConfig(option, config);
    Grid1D grid(config.minimumSpot, config.maximumSpot, config.spaceSteps);
    PdeWorkspace workspace(grid.nodeCount());
    workspace.previous = europeanTerminalPayoff(option, grid);
    if (option.maturity == 0.0)
    {
        return {grid.interpolate(workspace.previous, option.spot),
                std::move(grid), std::move(workspace.previous)};
    }

    const SpatialOperator op = makeSpatialOperator(option, grid);
    const double timeStep = option.maturity / static_cast<double>(config.timeSteps);
    if (config.scheme == PdeScheme::Explicit)
    {
        validateExplicitStability(op, grid, timeStep);
        for (std::size_t step = 1; step <= config.timeSteps; ++step)
        {
            explicitStep(option, grid, op, timeStep,
                         static_cast<double>(step) * timeStep, workspace);
        }
    }
    else if (config.scheme == PdeScheme::BackwardEuler)
    {
        const ThomasFactorization factorization{
            makeImplicitSystem(op, grid.nodeCount(), 1.0, timeStep)};
        for (std::size_t step = 1; step <= config.timeSteps; ++step)
        {
            thetaStep(option, grid, op, factorization, 1.0, timeStep,
                      static_cast<double>(step) * timeStep, workspace);
        }
    }
    else
    {
        const ThomasFactorization crankNicolson{
            makeImplicitSystem(op, grid.nodeCount(), 0.5, timeStep)};
        std::size_t firstCrankNicolsonStep = 1;
        if (config.scheme == PdeScheme::CrankNicolsonRannacher)
        {
            const double halfStep = 0.5 * timeStep;
            const ThomasFactorization backwardEulerHalfStep{
                makeImplicitSystem(op, grid.nodeCount(), 1.0, halfStep)};
            thetaStep(option, grid, op, backwardEulerHalfStep, 1.0, halfStep,
                      halfStep, workspace);
            thetaStep(option, grid, op, backwardEulerHalfStep, 1.0, halfStep,
                      timeStep, workspace);
            firstCrankNicolsonStep = 2;
        }
        for (std::size_t step = firstCrankNicolsonStep;
             step <= config.timeSteps; ++step)
        {
            thetaStep(option, grid, op, crankNicolson, 0.5, timeStep,
                      static_cast<double>(step) * timeStep, workspace);
        }
    }

    const double price = grid.interpolate(workspace.previous, option.spot);
    return {price, std::move(grid), std::move(workspace.previous)};
}

} // namespace dr3::lattice
