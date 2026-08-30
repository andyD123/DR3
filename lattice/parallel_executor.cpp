#include "parallel_executor.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <thread>
#include <vector>

namespace dr3::lattice
{

void parallelFor(std::size_t jobCount,
                 std::size_t threadCount,
                 const IndexedJob& job)
{
    if (threadCount < 1)
    {
        throw std::invalid_argument("parallel executor threadCount must be at least one");
    }
    if (!job)
    {
        throw std::invalid_argument("parallel executor job must be callable");
    }
    if (jobCount == 0)
    {
        return;
    }

    const std::size_t workerCount = std::min(threadCount, jobCount);
    std::vector<std::thread> workers;
    std::vector<std::exception_ptr> exceptions(workerCount);
    workers.reserve(workerCount);

    try
    {
        for (std::size_t worker = 0; worker < workerCount; ++worker)
        {
            const std::size_t begin = worker * jobCount / workerCount;
            const std::size_t end = (worker + 1) * jobCount / workerCount;
            workers.emplace_back([&, worker, begin, end]
            {
                try
                {
                    for (std::size_t index = begin; index < end; ++index)
                    {
                        job(index);
                    }
                }
                catch (...)
                {
                    exceptions[worker] = std::current_exception();
                }
            });
        }
    }
    catch (...)
    {
        for (auto& worker : workers)
        {
            if (worker.joinable()) worker.join();
        }
        throw;
    }

    for (auto& worker : workers) worker.join();
    for (const auto& exception : exceptions)
    {
        if (exception) std::rethrow_exception(exception);
    }
}

std::vector<double> parallelEuropeanPdePrices(
    const std::vector<VanillaOption>& options,
    const PdeConfig& config,
    std::size_t threadCount)
{
    std::vector<double> prices(options.size());
    parallelFor(options.size(), threadCount, [&](std::size_t index)
    {
        prices[index] = europeanPdePrice(options[index], config).price;
    });
    return prices;
}

double deterministicCompensatedSum(const std::vector<double>& values)
{
    double sum = 0.0;
    double correction = 0.0;
    for (const double value : values)
    {
        if (!std::isfinite(value))
        {
            throw std::invalid_argument("compensated sum inputs must be finite");
        }
        const double tentative = sum + value;
        if (std::abs(sum) >= std::abs(value))
        {
            correction += (sum - tentative) + value;
        }
        else
        {
            correction += (value - tentative) + sum;
        }
        sum = tentative;
    }
    return sum + correction;
}

} // namespace dr3::lattice
