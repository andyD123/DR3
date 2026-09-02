/****************************  similarity.h   *******************************
 * Exact float32 embedding similarity and top-k retrieval for DR3.
 * AVX2 (VecF8F) is the tested SIMD path. Double is the scalar reference.
 *****************************************************************************/
#pragma once

#include "dr3.h"
#include "target_name_space.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace DRC {
namespace AI {

namespace Simd = VecF8F;

struct SearchHit {
    int index;
    float score;
};

namespace detail {

inline Vec8f load_zero_padded_tail(const float* values, std::size_t n)
{
    alignas(32) float padded[8] = {};
    std::copy_n(values, n, padded);
    Vec8f result;
    result.load_a(padded);
    return result;
}

} // namespace detail

// Cosine of a zero-norm vector is this sentinel (never NaN).
constexpr float kZeroNormCosine = 0.f;
constexpr float kZeroNormEps = 1.0e-30f;

inline void require_same_length(std::size_t a, std::size_t b, const char* what)
{
    if (a != b) {
        throw std::invalid_argument(std::string(what) + ": mismatched lengths");
    }
}

inline Simd::VecXX load_vec(const float* p, std::size_t n)
{
    if (p == nullptr && n != 0) {
        throw std::invalid_argument("null embedding pointer");
    }
    if (n == 0) {
        return Simd::VecXX();
    }
    return Simd::VecXX(std::vector<float>(p, p + n));
}

inline float dot_product(const float* a, std::size_t na, const float* b, std::size_t nb)
{
    require_same_length(na, nb, "dot_product");
    if (na == 0) {
        return 0.f;
    }
    auto va = load_vec(a, na);
    auto vb = load_vec(b, nb);
    auto mul = [](auto x, auto y) { return x * y; };
    auto add = [](auto x, auto y) { return x + y; };
    return transformReduce(va, vb, mul, add);
}

inline float dot_product(Simd::SpanXX a, Simd::SpanXX b)
{
    require_same_length(a.size(), b.size(), "dot_product");
    if (a.empty()) {
        return 0.f;
    }
    constexpr std::size_t width = 8;
    Vec8f sum(0.f);
    std::size_t i = 0;
    for (; i + width <= a.size(); i += width) {
        Vec8f av, bv;
        av.load(a.start() + i);
        bv.load(b.start() + i);
        sum += av * bv;
    }
    if (i < a.size()) {
        const Vec8f av = detail::load_zero_padded_tail(a.start() + i, a.size() - i);
        const Vec8f bv = detail::load_zero_padded_tail(b.start() + i, b.size() - i);
        sum += av * bv;
    }
    return horizontal_add(sum);
}

inline float squared_l2_distance(const float* a, std::size_t na, const float* b, std::size_t nb)
{
    require_same_length(na, nb, "squared_l2_distance");
    if (na == 0) {
        return 0.f;
    }
    auto va = load_vec(a, na);
    auto vb = load_vec(b, nb);
    auto sqdiff = [](auto x, auto y) {
        auto d = x - y;
        return d * d;
    };
    auto add = [](auto x, auto y) { return x + y; };
    return transformReduce(va, vb, sqdiff, add);
}

inline float squared_l2_distance(Simd::SpanXX a, Simd::SpanXX b)
{
    require_same_length(a.size(), b.size(), "squared_l2_distance");
    if (a.empty()) {
        return 0.f;
    }
    constexpr std::size_t width = 8;
    Vec8f sum(0.f);
    std::size_t i = 0;
    for (; i + width <= a.size(); i += width) {
        Vec8f av, bv;
        av.load(a.start() + i);
        bv.load(b.start() + i);
        const Vec8f d = av - bv;
        sum += d * d;
    }
    if (i < a.size()) {
        const Vec8f av = detail::load_zero_padded_tail(a.start() + i, a.size() - i);
        const Vec8f bv = detail::load_zero_padded_tail(b.start() + i, b.size() - i);
        const Vec8f d = av - bv;
        sum += d * d;
    }
    return horizontal_add(sum);
}

inline float l2_norm2(const float* a, std::size_t n)
{
    return dot_product(a, n, a, n);
}

inline float cosine_similarity(const float* a, std::size_t na, const float* b, std::size_t nb)
{
    require_same_length(na, nb, "cosine_similarity");
    if (na == 0) {
        return kZeroNormCosine;
    }
    const float n2a = l2_norm2(a, na);
    const float n2b = l2_norm2(b, nb);
    if (!(n2a > kZeroNormEps) || !(n2b > kZeroNormEps)) {
        return kZeroNormCosine;
    }
    const float dot = dot_product(a, na, b, nb);
    const double denom = std::sqrt(static_cast<double>(n2a)) * std::sqrt(static_cast<double>(n2b));
    return static_cast<float>(static_cast<double>(dot) / denom);
}

inline float cosine_similarity(Simd::SpanXX a, Simd::SpanXX b)
{
    require_same_length(a.size(), b.size(), "cosine_similarity");
    if (a.empty()) {
        return kZeroNormCosine;
    }
    const float n2a = dot_product(a, a);
    const float n2b = dot_product(b, b);
    if (!(n2a > kZeroNormEps) || !(n2b > kZeroNormEps)) {
        return kZeroNormCosine;
    }
    const double denom = std::sqrt(static_cast<double>(n2a)) * std::sqrt(static_cast<double>(n2b));
    return static_cast<float>(static_cast<double>(dot_product(a, b)) / denom);
}

inline void normalize_l2_inplace(float* x, std::size_t n)
{
    if (n == 0) {
        return;
    }
    if (x == nullptr) {
        throw std::invalid_argument("normalize_l2_inplace: null pointer");
    }
    const float n2 = l2_norm2(x, n);
    if (!(n2 > kZeroNormEps)) {
        return;
    }
    const float inv = static_cast<float>(1.0 / std::sqrt(static_cast<double>(n2)));
    auto v = load_vec(x, n);
    auto scaled = v * inv;
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = scaled[static_cast<int>(i)];
    }
}

enum class TopKMetric { InnerProduct, Cosine, SquaredL2 };

inline bool better_hit(TopKMetric metric, const SearchHit& neu, const SearchHit& other)
{
    if (neu.score != other.score) {
        if (metric == TopKMetric::SquaredL2) {
            return neu.score < other.score;
        }
        return neu.score > other.score;
    }
    return neu.index < other.index;
}

struct HitWorse {
    TopKMetric metric;
    bool operator()(const SearchHit& a, const SearchHit& b) const
    {
        // Worst of the current best sits on top of the bounded heap.
        return better_hit(metric, a, b);
    }
};

inline bool passes_threshold(TopKMetric metric, float score, const float* threshold)
{
    if (threshold == nullptr) {
        return true;
    }
    if (metric == TopKMetric::SquaredL2) {
        return score <= *threshold;
    }
    return score >= *threshold;
}

inline std::vector<SearchHit> top_k(
    TopKMetric metric,
    const float* query,
    std::size_t d,
    const float* corpus,
    std::size_t n,
    std::size_t stride,
    std::size_t k,
    const float* threshold = nullptr)
{
    if (query == nullptr && d != 0) {
        throw std::invalid_argument("top_k: null query");
    }
    if (n > 0 && corpus == nullptr) {
        throw std::invalid_argument("top_k: null corpus");
    }
    if (stride < d) {
        throw std::invalid_argument("top_k: stride smaller than dimension");
    }
    if (k == 0 || n == 0 || d == 0) {
        return {};
    }

    const std::size_t keep = std::min(k, n);
    std::priority_queue<SearchHit, std::vector<SearchHit>, HitWorse> best{HitWorse{metric}};

    Simd::SpanXX query_span(const_cast<float*>(query), d);
    for (std::size_t i = 0; i < n; ++i) {
        const float* row = corpus + i * stride;
        Simd::SpanXX row_span(const_cast<float*>(row), d);
        float s = 0.f;
        switch (metric) {
        case TopKMetric::InnerProduct:
            s = dot_product(query_span, row_span);
            break;
        case TopKMetric::Cosine:
            s = cosine_similarity(query_span, row_span);
            break;
        case TopKMetric::SquaredL2:
            s = squared_l2_distance(query_span, row_span);
            break;
        }
        if (!passes_threshold(metric, s, threshold)) {
            continue;
        }
        SearchHit hit{static_cast<int>(i), s};
        if (best.size() < keep) {
            best.push(hit);
            continue;
        }
        if (better_hit(metric, hit, best.top())) {
            best.pop();
            best.push(hit);
        }
    }

    std::vector<SearchHit> out;
    out.reserve(best.size());
    while (!best.empty()) {
        out.push_back(best.top());
        best.pop();
    }
    std::sort(out.begin(), out.end(), [metric](const SearchHit& a, const SearchHit& b) {
        return better_hit(metric, a, b);
    });
    return out;
}

inline std::vector<SearchHit> top_k_inner_product(
    const float* query, std::size_t d,
    const float* corpus, std::size_t n,
    std::size_t k, const float* threshold = nullptr, std::size_t stride = 0)
{
    return top_k(TopKMetric::InnerProduct, query, d, corpus, n, stride == 0 ? d : stride, k, threshold);
}

inline std::vector<SearchHit> top_k_cosine_similarity(
    const float* query, std::size_t d,
    const float* corpus, std::size_t n,
    std::size_t k, const float* threshold = nullptr, std::size_t stride = 0)
{
    return top_k(TopKMetric::Cosine, query, d, corpus, n, stride == 0 ? d : stride, k, threshold);
}

inline std::vector<SearchHit> top_k_l2(
    const float* query, std::size_t d,
    const float* corpus, std::size_t n,
    std::size_t k, const float* threshold = nullptr, std::size_t stride = 0)
{
    return top_k(TopKMetric::SquaredL2, query, d, corpus, n, stride == 0 ? d : stride, k, threshold);
}

namespace ref {

inline double dot_product(const float* a, std::size_t na, const float* b, std::size_t nb)
{
    require_same_length(na, nb, "ref::dot_product");
    double s = 0.0;
    for (std::size_t i = 0; i < na; ++i) {
        s += static_cast<double>(a[i]) * static_cast<double>(b[i]);
    }
    return s;
}

inline double squared_l2_distance(const float* a, std::size_t na, const float* b, std::size_t nb)
{
    require_same_length(na, nb, "ref::squared_l2_distance");
    double s = 0.0;
    for (std::size_t i = 0; i < na; ++i) {
        const double d = static_cast<double>(a[i]) - static_cast<double>(b[i]);
        s += d * d;
    }
    return s;
}

inline double cosine_similarity(const float* a, std::size_t na, const float* b, std::size_t nb)
{
    require_same_length(na, nb, "ref::cosine_similarity");
    const double n2a = dot_product(a, na, a, na);
    const double n2b = dot_product(b, nb, b, nb);
    if (!(n2a > static_cast<double>(kZeroNormEps)) || !(n2b > static_cast<double>(kZeroNormEps))) {
        return static_cast<double>(kZeroNormCosine);
    }
    return dot_product(a, na, b, nb) / (std::sqrt(n2a) * std::sqrt(n2b));
}

} // namespace ref

} // namespace AI
} // namespace DRC
