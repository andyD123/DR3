#include "pch.h"

#include "../Vectorisation/VecX/alloc_policy.h"
#include "../Vectorisation/VecX/alloc_policy_imp.h"
#include "../Vectorisation/VecX/dr3.h"
#include "testNamespace.h"

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <set>
#include <thread>
#include <vector>

namespace
{

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

bool allocateFillVerifyFree(std::size_t logicalSize, double marker)
{
    std::size_t allocatedSize = logicalSize;
    double* values = nullptr;
    allocPool(allocatedSize, values);
    if (values == nullptr || reinterpret_cast<std::uintptr_t>(values) % ByteAllignment != 0)
    {
        return false;
    }
    for (std::size_t index = 0; index < allocatedSize; ++index)
    {
        values[index] = marker + static_cast<double>(index);
    }
    bool valid = true;
    for (std::size_t index = 0; index < allocatedSize; ++index)
    {
        valid = valid && values[index] == marker + static_cast<double>(index);
    }
    freePool(allocatedSize, values);
    return valid;
}

} // namespace

TEST(ParallelAllocator, ConcurrentFirstUse)
{
    freeAllAllocators(double{});
    constexpr std::size_t threadCount = 8;
    StartGate start(threadCount);
    std::atomic<bool> valid{true};
    std::vector<std::thread> workers;
    for (std::size_t worker = 0; worker < threadCount; ++worker)
    {
        workers.emplace_back([&, worker]
        {
            start.wait();
            if (!allocateFillVerifyFree(65, 1000.0 * static_cast<double>(worker)))
            {
                valid.store(false, std::memory_order_relaxed);
            }
        });
    }
    for (auto& worker : workers) worker.join();
    EXPECT_TRUE(valid.load());
}

TEST(ParallelAllocator, ConcurrentAllocateAndFree)
{
    constexpr std::size_t threadCount = 8;
    constexpr std::size_t iterations = 1000;
    StartGate start(threadCount);
    std::atomic<bool> valid{true};
    std::vector<std::thread> workers;
    for (std::size_t worker = 0; worker < threadCount; ++worker)
    {
        workers.emplace_back([&, worker]
        {
            start.wait();
            for (std::size_t iteration = 0; iteration < iterations; ++iteration)
            {
                if (!allocateFillVerifyFree(63 + iteration % 3,
                                             static_cast<double>(worker * iterations + iteration)))
                {
                    valid.store(false, std::memory_order_relaxed);
                    return;
                }
            }
        });
    }
    for (auto& worker : workers) worker.join();
    EXPECT_TRUE(valid.load());
}

TEST(ParallelAllocator, ConcurrentMixedLogicalSizes)
{
    const std::vector<std::size_t> sizes{1, 3, 4, 5, 63, 64, 65, 255, 256, 257};
    constexpr std::size_t threadCount = 8;
    StartGate start(threadCount);
    std::atomic<bool> valid{true};
    std::vector<std::thread> workers;
    for (std::size_t worker = 0; worker < threadCount; ++worker)
    {
        workers.emplace_back([&, worker]
        {
            start.wait();
            for (std::size_t iteration = 0; iteration < 1000; ++iteration)
            {
                const auto size = sizes[(worker + iteration) % sizes.size()];
                if (!allocateFillVerifyFree(size, static_cast<double>(worker + iteration)))
                {
                    valid.store(false, std::memory_order_relaxed);
                    return;
                }
            }
        });
    }
    for (auto& worker : workers) worker.join();
    EXPECT_TRUE(valid.load());
}

TEST(ParallelAllocator, NoDuplicateLiveBlocks)
{
    constexpr std::size_t threadCount = 8;
    StartGate start(threadCount);
    StartGate allocated(threadCount + 1);
    StartGate release(threadCount + 1);
    std::vector<double*> addresses(threadCount, nullptr);
    std::vector<std::size_t> allocatedSizes(threadCount);
    std::vector<std::thread> workers;
    for (std::size_t worker = 0; worker < threadCount; ++worker)
    {
        workers.emplace_back([&, worker]
        {
            start.wait();
            allocatedSizes[worker] = 65;
            allocPool(allocatedSizes[worker], addresses[worker]);
            allocated.wait();
            release.wait();
            freePool(allocatedSizes[worker], addresses[worker]);
        });
    }

    allocated.wait();
    const std::set<double*> unique(addresses.begin(), addresses.end());
    EXPECT_EQ(unique.size(), threadCount);
    EXPECT_EQ(unique.count(nullptr), 0u);
    release.wait();
    for (auto& worker : workers) worker.join();
}

TEST(ParallelAllocator, AlignmentIsPreserved)
{
    for (const std::size_t logicalSize : {1u, 3u, 5u, 63u, 64u, 65u, 257u})
    {
        std::size_t allocatedSize = logicalSize;
        double* values = nullptr;
        allocPool(allocatedSize, values);
        EXPECT_EQ(reinterpret_cast<std::uintptr_t>(values) % ByteAllignment, 0u);
        freePool(allocatedSize, values);
    }
}

TEST(ParallelAllocator, SameThreadAllocateFreeCycles)
{
    for (std::size_t iteration = 0; iteration < 5000; ++iteration)
    {
        ASSERT_TRUE(allocateFillVerifyFree(1 + iteration % 257,
                                           static_cast<double>(iteration)));
    }
}

TEST(ParallelAllocator, CrossThreadDestructionFollowsContract)
{
    std::size_t allocatedSize = 65;
    double* values = nullptr;
    std::thread allocateThread([&] { allocPool(allocatedSize, values); });
    allocateThread.join();
    ASSERT_NE(values, nullptr);
    std::thread freeThread([&] { freePool(allocatedSize, values); });
    freeThread.join();
    SUCCEED();
}

TEST(ParallelAllocator, WorkerTeardownReleasesWorkspace)
{
    constexpr std::size_t threadCount = 8;
    std::atomic<std::size_t> completed{0};
    std::vector<std::thread> workers;
    for (std::size_t worker = 0; worker < threadCount; ++worker)
    {
        workers.emplace_back([&]
        {
            VecXX workspace(1.0, 257);
            workspace[256] = 42.0;
            if (workspace[256] == 42.0)
            {
                completed.fetch_add(1, std::memory_order_relaxed);
            }
        });
    }
    for (auto& worker : workers) worker.join();
    EXPECT_EQ(completed.load(), threadCount);
}

TEST(ParallelAllocator, CleanupAfterWorkersJoin)
{
    std::atomic<bool> valid{false};
    std::thread worker([&valid]
    {
        VecXX value(2.0, 65);
        valid.store(value[64] == 2.0, std::memory_order_relaxed);
    });
    worker.join();
    EXPECT_TRUE(valid.load());
    freeAllAllocators(double{});
    EXPECT_TRUE(allocateFillVerifyFree(65, 7.0));
}

TEST(ParallelAllocator, RepeatedStartupAndShutdown)
{
    for (std::size_t repetition = 0; repetition < 10; ++repetition)
    {
        std::atomic<bool> valid{true};
        std::vector<std::thread> workers;
        for (std::size_t worker = 0; worker < 4; ++worker)
        {
            workers.emplace_back([&, worker]
            {
                if (!allocateFillVerifyFree(63 + worker, static_cast<double>(repetition)))
                {
                    valid.store(false, std::memory_order_relaxed);
                }
            });
        }
        for (auto& worker : workers) worker.join();
        ASSERT_TRUE(valid.load());
        freeAllAllocators(double{});
    }
}

TEST(ConcurrencySafety, IndependentVecOperationsMatchSerial)
{
    constexpr std::size_t threadCount = 8;
    StartGate start(threadCount);
    std::vector<double> results(threadCount);
    std::vector<std::thread> workers;
    for (std::size_t worker = 0; worker < threadCount; ++worker)
    {
        workers.emplace_back([&, worker]
        {
            start.wait();
            std::vector<double> input(65);
            for (std::size_t index = 0; index < input.size(); ++index)
            {
                input[index] = static_cast<double>(worker + index);
            }
            VecXX values(input);
            auto transformed = (values + 2.0) * 3.0;
            results[worker] = transformed[64];
        });
    }
    for (auto& worker : workers) worker.join();
    for (std::size_t worker = 0; worker < threadCount; ++worker)
    {
        EXPECT_DOUBLE_EQ(results[worker], 3.0 * (static_cast<double>(worker + 64) + 2.0));
    }
}

TEST(ConcurrencySafety, IndependentForwardAdMatchesSerial)
{
    constexpr std::size_t threadCount = 8;
    StartGate start(threadCount);
    std::vector<double> primal(threadCount), derivative(threadCount);
    std::vector<std::thread> workers;
    for (std::size_t worker = 0; worker < threadCount; ++worker)
    {
        workers.emplace_back([&, worker]
        {
            start.wait();
            VecXX inputs(1.0 + static_cast<double>(worker), 65);
            auto variables = D(inputs);
            auto result = variables * variables + 2.0 * variables;
            primal[worker] = result.value()[64];
            derivative[worker] = result.derivative()[64];
        });
    }
    for (auto& worker : workers) worker.join();
    for (std::size_t worker = 0; worker < threadCount; ++worker)
    {
        const double input = 1.0 + static_cast<double>(worker);
        EXPECT_DOUBLE_EQ(primal[worker], input * input + 2.0 * input);
        EXPECT_DOUBLE_EQ(derivative[worker], 2.0 * input + 2.0);
    }
}
