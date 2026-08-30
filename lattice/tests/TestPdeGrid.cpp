#include "pde_grid.h"
#include "pde_workspace.h"

#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>
#include <vector>

TEST(Grid1D, RejectsInvalidNodeCount)
{
    EXPECT_THROW((dr3::lattice::Grid1D{0.0, 1.0, 2}), std::invalid_argument);
    EXPECT_THROW((dr3::lattice::Grid1D{1.0, 1.0, 3}), std::invalid_argument);
    EXPECT_THROW((dr3::lattice::Grid1D{2.0, 1.0, 3}), std::invalid_argument);
    EXPECT_THROW((dr3::lattice::Grid1D{0.0,
                                       std::numeric_limits<double>::infinity(),
                                       3}),
                 std::invalid_argument);
}

TEST(Grid1D, ContainsExpectedEndpoints)
{
    const dr3::lattice::Grid1D grid{-2.5, 7.5, 11};
    EXPECT_DOUBLE_EQ(grid.minimum(), -2.5);
    EXPECT_DOUBLE_EQ(grid.maximum(), 7.5);
    EXPECT_DOUBLE_EQ(grid[grid.lowerBoundary()], -2.5);
    EXPECT_DOUBLE_EQ(grid[grid.upperBoundary()], 7.5);
    EXPECT_EQ(grid.interiorBegin(), 1u);
    EXPECT_EQ(grid.interiorEnd(), 10u);
}

TEST(Grid1D, HasUniformSpacing)
{
    const dr3::lattice::Grid1D grid{0.0, 2.0, 9};
    EXPECT_DOUBLE_EQ(grid.spacing(), 0.25);
    for (std::size_t index = 1; index < grid.nodeCount(); ++index)
    {
        EXPECT_NEAR(grid[index] - grid[index - 1], grid.spacing(), 1.0e-15);
    }
}

TEST(Grid1D, LocatesExactNode)
{
    const dr3::lattice::Grid1D grid{0.0, 10.0, 6};
    const auto middle = grid.locate(4.0);
    EXPECT_EQ(middle.lowerIndex, 2u);
    EXPECT_EQ(middle.upperIndex, 2u);
    EXPECT_DOUBLE_EQ(middle.upperWeight, 0.0);

    const auto final = grid.locate(10.0);
    EXPECT_EQ(final.lowerIndex, 5u);
    EXPECT_EQ(final.upperIndex, 5u);
}

TEST(Grid1D, InterpolatesBetweenNodes)
{
    const dr3::lattice::Grid1D grid{0.0, 4.0, 5};
    const std::vector<double> values{0.0, 10.0, 20.0, 30.0, 40.0};
    const auto location = grid.locate(2.25);

    EXPECT_EQ(location.lowerIndex, 2u);
    EXPECT_EQ(location.upperIndex, 3u);
    EXPECT_DOUBLE_EQ(location.upperWeight, 0.25);
    EXPECT_DOUBLE_EQ(grid.interpolate(values, 2.25), 22.5);
    EXPECT_DOUBLE_EQ(grid.interpolate(values, 4.0), 40.0);
    EXPECT_THROW(grid.interpolate(std::vector<double>(4), 2.0), std::invalid_argument);
    EXPECT_THROW(grid.locate(4.1), std::invalid_argument);
}

TEST(PdeWorkspace, OwnsReusableLogicalBuffers)
{
    const dr3::lattice::PdeWorkspace workspace{17};
    EXPECT_EQ(workspace.size(), 17u);
    EXPECT_EQ(workspace.previous.size(), 17u);
    EXPECT_EQ(workspace.next.size(), 17u);
    EXPECT_EQ(workspace.rightHandSide.size(), 17u);
    EXPECT_EQ(workspace.solution.size(), 17u);
    EXPECT_EQ(workspace.solver.size(), 17u);
    EXPECT_THROW((dr3::lattice::PdeWorkspace{0}), std::invalid_argument);
}
