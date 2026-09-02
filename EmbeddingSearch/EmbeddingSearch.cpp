// Exact embedding similarity and top-k search using DR3 AVX2 float kernels.
// Printed timings are measurements on this machine, not a claimed speedup.

#include "../Vectorisation/VecX/alloc_policy.h"
#include "../Vectorisation/VecX/similarity.h"
#include "embedding_fixture.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

AllAllocatorsGuard<float> allocGuard;

namespace {

int fail(const std::string& msg)
{
    std::cerr << "FAIL: " << msg << '\n';
    return 1;
}

bool close(double a, double b, double abs_tol, double rel_tol)
{
    const double diff = std::abs(a - b);
    return diff <= abs_tol || diff <= rel_tol * std::max({1.0, std::abs(a), std::abs(b)});
}

int run_self_test()
{
    {
        const float a[] = {1.f, 2.f};
        const float b[] = {1.f};
        try {
            DRC::AI::dot_product(a, 2, b, 1);
            return fail("mismatched lengths were not rejected");
        } catch (const std::invalid_argument&) {
        }
    }

    const float ident[] = {3.f, 4.f};
    const float cos = DRC::AI::cosine_similarity(ident, 2, ident, 2);
    const float l2 = DRC::AI::squared_l2_distance(ident, 2, ident, 2);
    const float dot = DRC::AI::dot_product(ident, 2, ident, 2);
    if (!close(cos, 1.0, 1e-5, 1e-5) || l2 > 1e-6f || !close(dot, 25.0, 1e-5, 1e-5)) {
        return fail("identical vector identities");
    }

    const float z[] = {0.f, 0.f};
    const float cz = DRC::AI::cosine_similarity(z, 2, ident, 2);
    if (std::isnan(cz) || cz != DRC::AI::kZeroNormCosine) {
        return fail("zero-norm cosine must be the documented sentinel");
    }

    const std::size_t dims[] = {1, 7, 8, 9, 15, 16, 17, 127, 128, 129};
    for (std::size_t d : dims) {
        std::vector<float> query;
        std::vector<float> corpus;
        dr3_embed::make_corpus(4, d, query, corpus);
        const float* row0 = corpus.data();
        if (!close(
                DRC::AI::dot_product(query.data(), d, row0, d),
                DRC::AI::ref::dot_product(query.data(), d, row0, d),
                1e-4, 1e-5)) {
            return fail("SIMD vs scalar dot at d=" + std::to_string(d));
        }
    }

    {
        std::vector<float> query;
        std::vector<float> corpus;
        dr3_embed::make_corpus(16, 8, query, corpus);
        auto hits = DRC::AI::top_k_inner_product(query.data(), 8, corpus.data(), 16, 5);
        auto again = DRC::AI::top_k_inner_product(query.data(), 8, corpus.data(), 16, 5);
        if (hits.size() != 5 || again.size() != 5) {
            return fail("top-k size");
        }
        for (std::size_t i = 0; i < hits.size(); ++i) {
            if (hits[i].index != again[i].index || hits[i].score != again[i].score) {
                return fail("top-k was not deterministic");
            }
        }
        const float cutoff = hits.back().score;
        auto filtered = DRC::AI::top_k_inner_product(
            query.data(), 8, corpus.data(), 16, 16, &cutoff);
        for (const auto& h : filtered) {
            if (h.score < cutoff) {
                return fail("threshold leaked a low score");
            }
        }
    }

    {
        const float q[] = {1.f, 0.f};
        const float c[] = {1.f, 0.f, 1.f, 0.f};
        auto hits = DRC::AI::top_k_inner_product(q, 2, c, 2, 2);
        if (hits.size() != 2 || hits[0].index != 0 || hits[1].index != 1) {
            return fail("tie-break is not smaller index first");
        }
    }

    std::cout << "EmbeddingSearch self-test passed\n";
    return 0;
}

void run_bench(std::size_t n, std::size_t d, std::size_t k)
{
    std::vector<float> query;
    std::vector<float> corpus;
    dr3_embed::make_corpus(n, d, query, corpus);

    auto t0 = std::chrono::steady_clock::now();
    auto hits = DRC::AI::top_k_cosine_similarity(query.data(), d, corpus.data(), n, k);
    auto t1 = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    double max_abs = 0.0;
    for (const auto& h : hits) {
        const double ref = DRC::AI::ref::cosine_similarity(
            query.data(), d, corpus.data() + static_cast<std::size_t>(h.index) * d, d);
        max_abs = std::max(max_abs, std::abs(ref - static_cast<double>(h.score)));
    }

    std::cout << "n=" << n << " d=" << d << " k=" << k
              << " topk_ms=" << ms
              << " hits=" << hits.size()
              << " max_abs_vs_scalar=" << max_abs << '\n';
    std::cout << "Times depend on CPU, compiler, instruction set, and build type.\n";
}

} // namespace

int main(int argc, char** argv)
{
    std::size_t n = 1024;
    std::size_t d = 128;
    std::size_t k = 10;
    bool self_test = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--self-test") {
            self_test = true;
        } else if (arg == "--size" && i + 1 < argc) {
            n = static_cast<std::size_t>(std::strtoull(argv[++i], nullptr, 10));
        } else if (arg == "--dim" && i + 1 < argc) {
            d = static_cast<std::size_t>(std::strtoull(argv[++i], nullptr, 10));
        } else if (arg == "--k" && i + 1 < argc) {
            k = static_cast<std::size_t>(std::strtoull(argv[++i], nullptr, 10));
        }
    }
    if (self_test) {
        return run_self_test();
    }
    run_bench(n, d, k);
    return 0;
}
