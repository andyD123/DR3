#include "pch.h"

#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/alloc_policy.h"
#include "../Vectorisation/VecX/target_name_space.h"
#include "../Vectorisation/VecX/dr3.h"
#include "testNamespace.h"
#include "dr3TestUtil.h"

#include <numeric>
#include <vector>

TEST(TestBasicPortable, MakeAndIndex)
{
    std::vector<Numeric> three(3, asNumber(42.0));
    std::vector<Numeric> mix{ asNumber(1.0), asNumber(2.0), asNumber(3.0) };
    VecXX vec2(three);
    VecXX vec1(mix);

    EXPECT_EQ(vec2.size(), 3);
    EXPECT_NUMERIC_EQ(vec1[0], asNumber(1.0));
    EXPECT_NUMERIC_EQ(vec1[1], asNumber(2.0));
    EXPECT_NUMERIC_EQ(vec1[2], asNumber(3.0));
}

TEST(TestBasicPortable, AddMultiply)
{
    std::vector<Numeric> mix{ asNumber(1.0), asNumber(2.0), asNumber(3.0) };
    VecXX vec1(mix);
    VecXX vec2(std::vector<Numeric>(3, asNumber(42.0)));

    auto added = vec1 + vec2;
    EXPECT_EQ(added.size(), 3);
    EXPECT_NUMERIC_EQ(added[0], asNumber(43.0));
    EXPECT_NUMERIC_EQ(added[1], asNumber(44.0));
    EXPECT_NUMERIC_EQ(added[2], asNumber(45.0));

    auto scaled = vec1 * asNumber(2.0);
    EXPECT_NUMERIC_EQ(scaled[0], asNumber(2.0));
    EXPECT_NUMERIC_EQ(scaled[1], asNumber(4.0));
    EXPECT_NUMERIC_EQ(scaled[2], asNumber(6.0));
}

TEST(TestBasicPortable, LengthNotMultipleOfSimdWidth)
{
    const int width = static_cast<int>(VecXX::INS::size());
    const int size = width * 2 + 1;
    std::vector<Numeric> input(static_cast<size_t>(size));
    std::iota(input.begin(), input.end(), asNumber(0.0));
    VecXX vec(input);

    EXPECT_EQ(vec.size(), size);
    EXPECT_NE(size % width, 0);

    auto doubled = vec + vec;
    EXPECT_EQ(doubled.size(), size);
    EXPECT_NUMERIC_EQ(doubled[0], asNumber(0.0));
    EXPECT_NUMERIC_EQ(doubled[size - 1], asNumber(2.0 * (size - 1)));
}

TEST(TestBasicPortable, TransformLambda)
{
    std::vector<Numeric> input{ asNumber(1.0), asNumber(2.0), asNumber(3.0),
        asNumber(4.0), asNumber(5.0) };
    VecXX vec(input);
    auto doubleIt = [](auto x) { return x + x; };
    auto doubled = transform(doubleIt, vec);
    EXPECT_EQ(doubled.size(), 5);
    EXPECT_NUMERIC_EQ(doubled[0], asNumber(2.0));
    EXPECT_NUMERIC_EQ(doubled[4], asNumber(10.0));
}
