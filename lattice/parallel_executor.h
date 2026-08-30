#pragma once

#include "european_pde.h"

#include <cstddef>
#include <functional>
#include <vector>

namespace dr3::lattice
{

using IndexedJob = std::function<void(std::size_t)>;

// Executes every index exactly once using deterministic contiguous partitions.
// The call returns only after all workers have joined. If workers throw, the
// exception from the lowest worker partition is rethrown on the calling thread.
void parallelFor(std::size_t jobCount,
                 std::size_t threadCount,
                 const IndexedJob& job);

std::vector<double> parallelEuropeanPdePrices(
    const std::vector<VanillaOption>& options,
    const PdeConfig& config,
    std::size_t threadCount);

// Deterministic Neumaier compensated reduction in input order.
double deterministicCompensatedSum(const std::vector<double>& values);

} // namespace dr3::lattice
