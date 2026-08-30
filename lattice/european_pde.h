#pragma once

#include "pde_boundary.h"
#include "pde_grid.h"
#include "pricing_types.h"

#include <cstddef>
#include <vector>

namespace dr3::lattice
{

enum class PdeScheme
{
    Explicit,
    BackwardEuler,
    CrankNicolson,
    CrankNicolsonRannacher
};

struct PdeConfig
{
    PdeScheme scheme{PdeScheme::CrankNicolsonRannacher};
    // Logical spot-grid node count, including both boundary nodes.
    std::size_t spaceSteps{401};
    std::size_t timeSteps{400};
    double minimumSpot{0.0};
    double maximumSpot{400.0};
};

struct PdeResult
{
    double price;
    Grid1D grid;
    std::vector<double> values;
};

std::vector<double> europeanTerminalPayoff(const VanillaOption& option,
                                           const Grid1D& grid);

PdeBoundaryValues europeanBoundaryValues(const VanillaOption& option,
                                         double minimumSpot,
                                         double maximumSpot,
                                         double timeToMaturity);

PdeResult europeanPdePrice(const VanillaOption& option, const PdeConfig& config);

} // namespace dr3::lattice
