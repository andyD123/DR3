#include "pch.h"

#include "../Vectorisation/VecX/alloc_policy.h"
#include "../Vectorisation/VecX/similarity.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

using DRC::AI::SearchHit;
using DRC::AI::kZeroNormCosine;

namespace {

AllAllocatorsGuard<float> allocGuard;

bool close(double a, double b, double abs_tol, double rel_tol)
{
    const double diff = std::abs(a - b);
    return diff <= abs_tol || diff <= rel_tol * std::max({1.0, std::abs(a), std::abs(b)});
}

std::vector<float> filled(std::size_t n, float v)
{
    return std::vector<float>(n, v);
}

std::vector<float> iota_f(std::size_t n, float start = 1.f)
{
    std::vector<float> x(n);
    std::iota(x.begin(), x.end(), start);
    return x;
}

std::vector<float> seeded(std::size_t n, std::uint32_t seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(-1.f, 1.f);
    std::vector<float> x(n);
    for (auto& v : x) {
        v = dist(rng);
    }
    return x;
}

} // namespace

TEST(TestSimilarity, MismatchedLengthsAreRejected)
{
    const float a[] = {1.f, 2.f};
    const float b[] = {1.f};
    EXPECT_THROW(DRC::AI::dot_product(a, 2, b, 1), std::invalid_argument);
    EXPECT_THROW(DRC::AI::squared_l2_distance(a, 2, b, 1), std::invalid_argument);
    EXPECT_THROW(DRC::AI::cosine_similarity(a, 2, b, 1), std::invalid_argument);
}

TEST(TestSimilarity, NonOwningSpansMatchPointerKernels)
{
    std::vector<float> a = {1.f, -2.f, 3.f, 0.5f, 7.f, -4.f, 2.f, 9.f, -1.f};
    std::vector<float> b = {-3.f, 1.f, 2.f, 4.f, -1.f, 5.f, 0.25f, 2.f, 8.f};
    DRC::AI::Simd::SpanXX as(a.data(), a.size());
    DRC::AI::Simd::SpanXX bs(b.data(), b.size());
    EXPECT_NEAR(DRC::AI::dot_product(as, bs),
        DRC::AI::dot_product(a.data(), a.size(), b.data(), b.size()), 1e-5f);
    EXPECT_NEAR(DRC::AI::squared_l2_distance(as, bs),
        DRC::AI::squared_l2_distance(a.data(), a.size(), b.data(), b.size()), 1e-5f);
    EXPECT_NEAR(DRC::AI::cosine_similarity(as, bs),
        DRC::AI::cosine_similarity(a.data(), a.size(), b.data(), b.size()), 1e-5f);
}

TEST(TestSimilarity, NonOwningSpanTailsAreAlignmentSafe)
{
    constexpr std::size_t dims[] = {1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15};
    alignas(32) std::array<float, 20> a_storage{};
    alignas(32) std::array<float, 20> b_storage{};
    for (std::size_t offset = 0; offset < 4; ++offset) {
        float* a = a_storage.data() + offset;
        float* b = b_storage.data() + offset;
        for (std::size_t i = 0; i < 16; ++i) {
            a[i] = static_cast<float>(i + 1) * 0.125f;
            b[i] = static_cast<float>(17 - i) * -0.0625f;
        }
        for (std::size_t d : dims) {
            DRC::AI::Simd::SpanXX as(a, d);
            DRC::AI::Simd::SpanXX bs(b, d);
            EXPECT_TRUE(close(
                DRC::AI::dot_product(as, bs),
                DRC::AI::ref::dot_product(a, d, b, d),
                1e-5, 1e-5))
                << "dot offset=" << offset << " d=" << d;
            EXPECT_TRUE(close(
                DRC::AI::squared_l2_distance(as, bs),
                DRC::AI::ref::squared_l2_distance(a, d, b, d),
                1e-5, 1e-5))
                << "l2 offset=" << offset << " d=" << d;
        }
    }
}

TEST(TestSimilarity, IdenticalVectors)
{
    const std::size_t dims[] = {1, 7, 8, 9, 15, 16, 17, 127, 128, 129};
    for (std::size_t d : dims) {
        auto x = iota_f(d);
        const float dot = DRC::AI::dot_product(x.data(), d, x.data(), d);
        const float l2 = DRC::AI::squared_l2_distance(x.data(), d, x.data(), d);
        const float cos = DRC::AI::cosine_similarity(x.data(), d, x.data(), d);
        const double ref_dot = DRC::AI::ref::dot_product(x.data(), d, x.data(), d);
        EXPECT_TRUE(close(dot, ref_dot, 1e-5, 1e-5)) << "d=" << d;
        EXPECT_NEAR(l2, 0.f, 1e-5f) << "d=" << d;
        EXPECT_TRUE(close(cos, 1.0, 1e-5, 1e-5)) << "d=" << d;
        EXPECT_TRUE(close(dot, ref_dot, 1e-5, 1e-5));
    }
}

TEST(TestSimilarity, OrthogonalAndOpposite)
{
    const float e0[] = {1.f, 0.f, 0.f};
    const float e1[] = {0.f, 1.f, 0.f};
    const float neg[] = {-1.f, 0.f, 0.f};
    EXPECT_NEAR(DRC::AI::cosine_similarity(e0, 3, e1, 3), 0.f, 1e-6f);
    EXPECT_NEAR(DRC::AI::cosine_similarity(e0, 3, neg, 3), -1.f, 1e-6f);
    EXPECT_NEAR(DRC::AI::ref::cosine_similarity(e0, 3, e1, 3), 0.0, 1e-12);
}

TEST(TestSimilarity, ZeroVectorsNoNaN)
{
    auto z = filled(8, 0.f);
    auto x = iota_f(8);
    const float cos_zz = DRC::AI::cosine_similarity(z.data(), 8, z.data(), 8);
    const float cos_zx = DRC::AI::cosine_similarity(z.data(), 8, x.data(), 8);
    EXPECT_FALSE(std::isnan(cos_zz));
    EXPECT_FALSE(std::isnan(cos_zx));
    EXPECT_EQ(cos_zz, kZeroNormCosine);
    EXPECT_EQ(cos_zx, kZeroNormCosine);
    EXPECT_NEAR(DRC::AI::dot_product(z.data(), 8, x.data(), 8), 0.f, 1e-7f);
    EXPECT_FALSE(std::isnan(DRC::AI::squared_l2_distance(z.data(), 8, x.data(), 8)));
}

TEST(TestSimilarity, SmallAndLargeMagnitudes)
{
    auto tiny = filled(9, 1.0e-20f);
    auto huge = filled(9, 1.0e10f);
    auto ones = filled(9, 1.f);
    EXPECT_FALSE(std::isnan(DRC::AI::cosine_similarity(tiny.data(), 9, tiny.data(), 9)));
    EXPECT_NEAR(DRC::AI::cosine_similarity(huge.data(), 9, huge.data(), 9), 1.f, 1e-5f);
    const float d = DRC::AI::dot_product(ones.data(), 9, huge.data(), 9);
    EXPECT_TRUE(close(d, DRC::AI::ref::dot_product(ones.data(), 9, huge.data(), 9), 1.0, 1e-5));
}

TEST(TestSimilarity, TailsAcrossWidths)
{
    const std::size_t dims[] = {1, 7, 8, 9, 15, 16, 17, 127, 128, 129};
    for (std::size_t d : dims) {
        auto a = seeded(d, 11u);
        auto b = seeded(d, 23u);
        EXPECT_TRUE(close(
            DRC::AI::dot_product(a.data(), d, b.data(), d),
            DRC::AI::ref::dot_product(a.data(), d, b.data(), d),
            1e-4, 1e-5))
            << "dot d=" << d;
        EXPECT_TRUE(close(
            DRC::AI::squared_l2_distance(a.data(), d, b.data(), d),
            DRC::AI::ref::squared_l2_distance(a.data(), d, b.data(), d),
            1e-4, 1e-5))
            << "l2 d=" << d;
        EXPECT_TRUE(close(
            DRC::AI::cosine_similarity(a.data(), d, b.data(), d),
            DRC::AI::ref::cosine_similarity(a.data(), d, b.data(), d),
            1e-5, 1e-5))
            << "cos d=" << d;
    }
}

TEST(TestSimilarity, NormalizeL2)
{
    auto x = iota_f(17);
    DRC::AI::normalize_l2_inplace(x.data(), x.size());
    const float n2 = DRC::AI::l2_norm2(x.data(), x.size());
    EXPECT_NEAR(n2, 1.f, 1e-5f);
    auto z = filled(4, 0.f);
    DRC::AI::normalize_l2_inplace(z.data(), z.size());
    EXPECT_EQ(z[0], 0.f);
}

TEST(TestSimilarity, TieBreakBySmallerIndex)
{
    const float query[] = {1.f, 0.f};
    const float corpus[] = {
        1.f, 0.f,
        1.f, 0.f,
        0.f, 1.f,
    };
    auto hits = DRC::AI::top_k_inner_product(query, 2, corpus, 3, 2);
    ASSERT_EQ(hits.size(), 2u);
    EXPECT_EQ(hits[0].index, 0);
    EXPECT_EQ(hits[1].index, 1);
    EXPECT_NEAR(hits[0].score, hits[1].score, 1e-6f);
}

TEST(TestSimilarity, CorpusSmallerThanKReturnsAllSorted)
{
    const float query[] = {1.f, 0.f};
    const float corpus[] = {
        0.f, 1.f,
        2.f, 0.f,
        1.f, 0.f,
    };
    auto hits = DRC::AI::top_k_inner_product(query, 2, corpus, 3, 10);
    ASSERT_EQ(hits.size(), 3u);
    EXPECT_EQ(hits[0].index, 1);
    EXPECT_EQ(hits[1].index, 2);
    EXPECT_EQ(hits[2].index, 0);
}

TEST(TestSimilarity, PaddedRowStrideDoesNotCopyOrScorePadding)
{
    const float query[] = {1.f, 0.f};
    const float corpus[] = {
        1.f, 0.f, 1000.f, 1000.f,
        3.f, 0.f, -1000.f, -1000.f,
        2.f, 0.f, 500.f, 500.f,
    };
    auto hits = DRC::AI::top_k_inner_product(query, 2, corpus, 3, 3, nullptr, 4);
    ASSERT_EQ(hits.size(), 3u);
    EXPECT_EQ(hits[0].index, 1);
    EXPECT_EQ(hits[1].index, 2);
    EXPECT_EQ(hits[2].index, 0);
    EXPECT_FLOAT_EQ(hits[0].score, 3.f);
}

TEST(TestSimilarity, ThresholdOmitsLowScores)
{
    const float query[] = {1.f, 0.f};
    const float corpus[] = {
        1.f, 0.f,
        0.1f, 0.f,
        0.5f, 0.f,
    };
    const float cutoff = 0.4f;
    auto hits = DRC::AI::top_k_inner_product(query, 2, corpus, 3, 3, &cutoff);
    ASSERT_EQ(hits.size(), 2u);
    EXPECT_EQ(hits[0].index, 0);
    EXPECT_EQ(hits[1].index, 2);
}

TEST(TestSimilarity, L2TopKPrefersNearest)
{
    const float query[] = {0.f, 0.f};
    const float corpus[] = {
        3.f, 0.f,
        1.f, 0.f,
        2.f, 0.f,
    };
    auto hits = DRC::AI::top_k_l2(query, 2, corpus, 3, 2);
    ASSERT_EQ(hits.size(), 2u);
    EXPECT_EQ(hits[0].index, 1);
    EXPECT_EQ(hits[1].index, 2);
}

TEST(TestSimilarity, DeterministicRandomFixtures)
{
    constexpr std::size_t n = 32;
    constexpr std::size_t d = 16;
    auto query = seeded(d, 99u);
    std::vector<float> corpus(n * d);
    for (std::size_t i = 0; i < n; ++i) {
        auto row = seeded(d, static_cast<std::uint32_t>(1000 + i));
        std::copy(row.begin(), row.end(), corpus.begin() + static_cast<std::ptrdiff_t>(i * d));
    }
    auto a = DRC::AI::top_k_cosine_similarity(query.data(), d, corpus.data(), n, 5);
    auto b = DRC::AI::top_k_cosine_similarity(query.data(), d, corpus.data(), n, 5);
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        EXPECT_EQ(a[i].index, b[i].index);
        EXPECT_EQ(a[i].score, b[i].score);
    }
    // SIMD top-1 matches scalar argmax of cosine.
    double best = -2.0;
    int best_i = -1;
    for (std::size_t i = 0; i < n; ++i) {
        const double c = DRC::AI::ref::cosine_similarity(
            query.data(), d, corpus.data() + i * d, d);
        if (c > best || (c == best && static_cast<int>(i) < best_i)) {
            best = c;
            best_i = static_cast<int>(i);
        }
    }
    EXPECT_EQ(a[0].index, best_i);
}

TEST(TestSimilarity, EmptyInputs)
{
    EXPECT_EQ(DRC::AI::dot_product(nullptr, 0, nullptr, 0), 0.f);
    auto hits = DRC::AI::top_k_inner_product(nullptr, 0, nullptr, 0, 3);
    EXPECT_TRUE(hits.empty());
}
