#include "parallel_executor.h"
#include "tree_pricers.h"

#include <gtest/gtest.h>

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace
{

using dr3::lattice::OptionType;
using dr3::lattice::PdeConfig;
using dr3::lattice::PdeScheme;
using dr3::lattice::VanillaOption;

class StartGate
{
public:
    explicit StartGate(std::size_t participants) : remaining_(participants) {}

    void wait()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        if (--remaining_ == 0)
        {
            open_ = true;
            condition_.notify_all();
            return;
        }
        condition_.wait(lock, [this] { return open_; });
    }

private:
    std::mutex mutex_;
    std::condition_variable condition_;
    std::size_t remaining_;
    bool open_{false};
};

std::vector<VanillaOption> options(std::size_t count)
{
    std::vector<VanillaOption> result;
    result.reserve(count);
    for (std::size_t index = 0; index < count; ++index)
    {
        result.push_back({80.0 + static_cast<double>(index % 9) * 5.0,
                          100.0,
                          0.15 + static_cast<double>(index % 4) * 0.05,
                          0.01 + static_cast<double>(index % 3) * 0.02,
                          static_cast<double>(index % 2) * 0.01,
                          0.5 + static_cast<double>(index % 5) * 0.25,
                          index % 2 == 0 ? OptionType::Call : OptionType::Put});
    }
    return result;
}

PdeConfig pdeConfig()
{
    return {PdeScheme::CrankNicolsonRannacher, 101, 100, 0.0, 400.0};
}

std::vector<double> serialPdePrices(const std::vector<VanillaOption>& instruments)
{
    std::vector<double> prices;
    prices.reserve(instruments.size());
    const auto config = pdeConfig();
    for (const auto& instrument : instruments)
    {
        prices.push_back(dr3::lattice::europeanPdePrice(instrument, config).price);
    }
    return prices;
}

} // namespace

TEST(ParallelExecutor, OneThreadMatchesSerial)
{
    const auto instruments = options(9);
    EXPECT_EQ(dr3::lattice::parallelEuropeanPdePrices(instruments, pdeConfig(), 1),
              serialPdePrices(instruments));
}

TEST(ParallelExecutor, TwoThreadsMatchSerial)
{
    const auto instruments = options(10);
    EXPECT_EQ(dr3::lattice::parallelEuropeanPdePrices(instruments, pdeConfig(), 2),
              serialPdePrices(instruments));
}

TEST(ParallelExecutor, MoreThreadsThanJobs)
{
    const auto instruments = options(3);
    EXPECT_EQ(dr3::lattice::parallelEuropeanPdePrices(instruments, pdeConfig(), 16),
              serialPdePrices(instruments));
}

TEST(ParallelExecutor, EmptyBatch)
{
    EXPECT_TRUE(dr3::lattice::parallelEuropeanPdePrices({}, pdeConfig(), 4).empty());
}

TEST(ParallelExecutor, OddNumberOfBatches)
{
    const auto instruments = options(17);
    EXPECT_EQ(dr3::lattice::parallelEuropeanPdePrices(instruments, pdeConfig(), 4),
              serialPdePrices(instruments));
}

TEST(ParallelExecutor, NonMultipleOfSimdWidth)
{
    const auto instruments = options(5);
    ASSERT_NE(instruments.size() % 4, 0u);
    EXPECT_EQ(dr3::lattice::parallelEuropeanPdePrices(instruments, pdeConfig(), 3),
              serialPdePrices(instruments));
}

TEST(ParallelExecutor, EachJobProcessedExactlyOnce)
{
    constexpr std::size_t jobCount = 37;
    std::vector<std::atomic<unsigned int>> counts(jobCount);
    for (auto& count : counts) count.store(0, std::memory_order_relaxed);
    dr3::lattice::parallelFor(jobCount, 8, [&](std::size_t index)
    {
        counts[index].fetch_add(1, std::memory_order_relaxed);
    });
    for (const auto& count : counts) EXPECT_EQ(count.load(), 1u);
}

TEST(ParallelExecutor, PropagatesWorkerException)
{
    EXPECT_THROW(dr3::lattice::parallelFor(8, 4, [](std::size_t index)
                 {
                     if (index == 5) throw std::runtime_error("worker failure");
                 }),
                 std::runtime_error);
}

TEST(ParallelExecutor, RepeatedExecutionIsDeterministic)
{
    const auto instruments = options(21);
    const auto expected = dr3::lattice::parallelEuropeanPdePrices(
        instruments, pdeConfig(), 4);
    for (int repetition = 0; repetition < 20; ++repetition)
    {
        EXPECT_EQ(dr3::lattice::parallelEuropeanPdePrices(
                      instruments, pdeConfig(), 4),
                  expected);
    }
}

TEST(ParallelExecutor, RejectsZeroThreads)
{
    EXPECT_THROW(dr3::lattice::parallelFor(1, 0, [](std::size_t) {}),
                 std::invalid_argument);
    EXPECT_THROW(dr3::lattice::parallelEuropeanPdePrices({}, pdeConfig(), 0),
                 std::invalid_argument);
}

TEST(ParallelExecutor, CompensatedAggregateIsDeterministic)
{
    const std::vector<double> values{1.0e16, 1.0, -1.0e16, 3.0};
    EXPECT_DOUBLE_EQ(dr3::lattice::deterministicCompensatedSum(values), 4.0);
    EXPECT_DOUBLE_EQ(dr3::lattice::deterministicCompensatedSum(values), 4.0);
}

TEST(ConcurrencySafety, IndependentTreePricesMatchSerial)
{
    constexpr std::size_t threadCount = 8;
    const auto instruments = options(threadCount);
    std::vector<double> serial(threadCount), concurrent(threadCount);
    for (std::size_t index = 0; index < threadCount; ++index)
    {
        serial[index] = dr3::lattice::europeanTrinomial(instruments[index], {257});
    }

    StartGate start(threadCount);
    std::vector<std::thread> workers;
    for (std::size_t index = 0; index < threadCount; ++index)
    {
        workers.emplace_back([&, index]
        {
            start.wait();
            concurrent[index] = dr3::lattice::europeanTrinomial(
                instruments[index], {257});
        });
    }
    for (auto& worker : workers) worker.join();
    EXPECT_EQ(concurrent, serial);
}

TEST(ConcurrencySafety, IndependentEuropeanPdeMatchesSerial)
{
    constexpr std::size_t threadCount = 8;
    const auto instruments = options(threadCount);
    const auto serial = serialPdePrices(instruments);
    std::vector<double> concurrent(threadCount);
    StartGate start(threadCount);
    std::vector<std::thread> workers;
    for (std::size_t index = 0; index < threadCount; ++index)
    {
        workers.emplace_back([&, index]
        {
            start.wait();
            concurrent[index] = dr3::lattice::europeanPdePrice(
                instruments[index], pdeConfig()).price;
        });
    }
    for (auto& worker : workers) worker.join();
    EXPECT_EQ(concurrent, serial);
}

TEST(ConcurrencySafety, NonMultipleOfSimdWidth)
{
    const auto instruments = options(13);
    ASSERT_NE(instruments.size() % 4, 0u);
    EXPECT_EQ(dr3::lattice::parallelEuropeanPdePrices(instruments, pdeConfig(), 4),
              serialPdePrices(instruments));
}

TEST(ConcurrencySafety, MoreThreadsThanJobs)
{
    const auto instruments = options(2);
    EXPECT_EQ(dr3::lattice::parallelEuropeanPdePrices(instruments, pdeConfig(), 32),
              serialPdePrices(instruments));
}

TEST(ConcurrencySafety, RepeatedRunsAreDeterministic)
{
    const auto instruments = options(11);
    const auto expected = serialPdePrices(instruments);
    for (int repetition = 0; repetition < 10; ++repetition)
    {
        EXPECT_EQ(dr3::lattice::parallelEuropeanPdePrices(
                      instruments, pdeConfig(), 4),
                  expected);
    }
}

TEST(ConcurrencySafety, ExceptionsReachCallingThread)
{
    auto instruments = options(8);
    instruments[6].volatility = 0.0;
    EXPECT_THROW(dr3::lattice::parallelEuropeanPdePrices(
                     instruments, pdeConfig(), 4),
                 std::invalid_argument);
}

TEST(ConcurrencySafety, EveryJobCompletesExactlyOnce)
{
    constexpr std::size_t jobCount = 65;
    std::vector<std::atomic<unsigned int>> counts(jobCount);
    for (auto& count : counts) count.store(0, std::memory_order_relaxed);
    dr3::lattice::parallelFor(jobCount, 8, [&](std::size_t index)
    {
        counts[index].fetch_add(1, std::memory_order_relaxed);
    });
    for (const auto& count : counts) EXPECT_EQ(count.load(), 1u);
}
