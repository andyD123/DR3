#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace dr3::lattice
{

struct PdeWorkspace
{
    explicit PdeWorkspace(std::size_t nodeCount)
        : previous(nodeCount), next(nodeCount), rightHandSide(nodeCount),
          solution(nodeCount), solver(nodeCount)
    {
        if (nodeCount == 0)
        {
            throw std::invalid_argument("PDE workspace must contain at least one node");
        }
    }

    std::size_t size() const noexcept { return previous.size(); }

    std::vector<double> previous;
    std::vector<double> next;
    std::vector<double> rightHandSide;
    std::vector<double> solution;
    std::vector<double> solver;
};

} // namespace dr3::lattice
