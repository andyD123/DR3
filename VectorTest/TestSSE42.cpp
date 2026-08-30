// Compile this translation unit with SSE4.2 so VCL INSTRSET >= 6 and
// the existing 128-bit DRC::VecD2D path uses SSE4.2 Vec2d operations.
// The AVX2 VectorTest suite is unchanged.

#include "pch.h"

#include "../Vectorisation/VCL/instrset.h"
#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/alloc_policy.h"
#include "../Vectorisation/VecX/target_name_space.h"
#include "../Vectorisation/VecX/dr3.h"

#include <cmath>
#include <numeric>
#include <vector>

using namespace DRC::VecD2D;

AllAllocatorsGuard<typename VecXX::SCALA_TYPE> allocGuard;

TEST(TestSSE42, CompiledInstructionSetIsAtLeastSse42)
{
    EXPECT_GE(INSTRSET, 6);
}

TEST(TestSSE42, VecD2DWidthIsTwo)
{
    EXPECT_EQ(static_cast<int>(VecXX::INS::size()), 2);
    EXPECT_EQ(InstructionTraits<VecXX::INS>::width, 2);
}

TEST(TestSSE42, ArithmeticAndTail)
{
    std::vector<double> mix{1.0, 2.0, 3.0};
    VecXX vec(mix);
    EXPECT_EQ(vec.size(), 3);
    EXPECT_NE(vec.size() % 2, 0);

    auto added = vec + vec;
    EXPECT_DOUBLE_EQ(added[0], 2.0);
    EXPECT_DOUBLE_EQ(added[1], 4.0);
    EXPECT_DOUBLE_EQ(added[2], 6.0);

    auto scaled = vec * 0.5;
    EXPECT_DOUBLE_EQ(scaled[0], 0.5);
    EXPECT_DOUBLE_EQ(scaled[2], 1.5);
}

TEST(TestSSE42, TransformNamedLambda)
{
    std::vector<double> input(5);
    std::iota(input.begin(), input.end(), 1.0);
    VecXX vec(input);
    auto square = [](auto x) { return x * x; };
    auto out = transform(square, vec);
    EXPECT_EQ(out.size(), 5);
    EXPECT_DOUBLE_EQ(out[0], 1.0);
    EXPECT_DOUBLE_EQ(out[4], 25.0);
}

TEST(TestSSE42, FloorAndMax)
{
    std::vector<double> input{1.2, -1.8, 2.5};
    VecXX vec(input);
    auto floored = floor(vec);
    EXPECT_EQ(floored.size(), 3);
    EXPECT_DOUBLE_EQ(floored[0], std::floor(1.2));
    EXPECT_DOUBLE_EQ(floored[1], std::floor(-1.8));
    EXPECT_DOUBLE_EQ(floored[2], std::floor(2.5));

    VecXX ones(std::vector<double>{1.0, 1.0, 1.0});
    auto hi = max(vec, ones);
    EXPECT_DOUBLE_EQ(hi[0], 1.2);
    EXPECT_DOUBLE_EQ(hi[1], 1.0);
    EXPECT_DOUBLE_EQ(hi[2], 2.5);
}
