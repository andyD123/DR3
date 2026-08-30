// SIMD portfolio pricing and Greeks demonstration for DR3.
// Scalar analytic Black-Scholes is the reference. DR3 AVX2 prices the
// book; finite differences (and AAD when available) are compared.
// Measurements are printed for this machine and are not a speedup claim.

#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/target_name_space.h"
#include "black_scholes_reference.h"
#include "portfolio_fixture.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace DRC::VecD4D;
using dr3_example::OptionType;
using dr3_example::Position;

namespace {

constexpr double kPriceTol = 1.0e-8;
constexpr double kGreekTol = 5.0e-4;
constexpr double kParityTol = 1.0e-10;
constexpr double kFdBump = 1.0e-5;

AllAllocatorsGuard<typename VecXX::SCALA_TYPE> allocGuard;

struct PackedBook {
    std::vector<double> spot;
    std::vector<double> strike;
    std::vector<double> maturity;
    std::vector<double> vol;
    std::vector<double> rate;
    std::vector<double> quantity;
    std::vector<double> is_call;
};

PackedBook pack(const std::vector<Position>& book)
{
    PackedBook p;
    p.spot.reserve(book.size());
    for (const auto& pos : book) {
        p.spot.push_back(pos.spot);
        p.strike.push_back(pos.strike);
        p.maturity.push_back(pos.maturity);
        p.vol.push_back(pos.vol);
        p.rate.push_back(pos.rate);
        p.quantity.push_back(pos.quantity);
        p.is_call.push_back(pos.type == OptionType::Call ? 1.0 : 0.0);
    }
    return p;
}

template <typename x>
auto black_scholes_call_vec(const x& S, const x& K, const x& t, const x& r, const x& sigma)
{
    auto invK = 1.0 / K;
    auto discountedRate = exp(-r * t);
    auto log_sK = log(S * invK);
    auto rootT = sqrt(t);
    auto sigmaRootT = rootT * sigma;
    auto d1 = (log_sK + (0.5 * sigma * sigma + r) * t) / sigmaRootT;
    auto d2 = d1 - sigmaRootT;
    return S * cdfnorm(d1) - K * discountedRate * cdfnorm(d2);
}

std::vector<double> simd_call_prices(const PackedBook& book)
{
    VecXX S(book.spot);
    VecXX K(book.strike);
    VecXX t(book.maturity);
    VecXX r(book.rate);
    VecXX sigma(book.vol);
    VecXX calls = black_scholes_call_vec(S, K, t, r, sigma);
    return std::vector<double>(calls);
}

double naive_sum(const std::vector<double>& x)
{
    return std::accumulate(x.begin(), x.end(), 0.0);
}

double pairwise_sum(const std::vector<double>& x)
{
    if (x.empty()) {
        return 0.0;
    }
    std::vector<double> layer = x;
    while (layer.size() > 1) {
        std::vector<double> next;
        next.reserve((layer.size() + 1) / 2);
        for (size_t i = 0; i + 1 < layer.size(); i += 2) {
            next.push_back(layer[i] + layer[i + 1]);
        }
        if (layer.size() % 2 == 1) {
            next.push_back(layer.back());
        }
        layer.swap(next);
    }
    return layer.front();
}

struct ScalarBook {
    std::vector<double> price;
    std::vector<double> delta;
    std::vector<double> vega;
    std::vector<double> theta;
    std::vector<double> weighted;
};

ScalarBook scalar_book(const std::vector<Position>& book)
{
    ScalarBook out;
    out.price.reserve(book.size());
    for (const auto& p : book) {
        auto g = (p.type == OptionType::Call)
            ? dr3_example::black_scholes_call(p.spot, p.strike, p.maturity, p.rate, p.vol)
            : dr3_example::black_scholes_put(p.spot, p.strike, p.maturity, p.rate, p.vol);
        out.price.push_back(g.price);
        out.delta.push_back(g.delta);
        out.vega.push_back(g.vega);
        out.theta.push_back(g.theta);
        out.weighted.push_back(g.price * p.quantity);
    }
    return out;
}

double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    double m = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        m = std::max(m, std::abs(a[i] - b[i]));
    }
    return m;
}

double max_rel_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    double m = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        const double denom = std::max(1.0, std::abs(a[i]));
        m = std::max(m, std::abs(a[i] - b[i]) / denom);
    }
    return m;
}

bool nearly_eq(double a, double b, double tol)
{
    return std::abs(a - b) <= tol * std::max(1.0, std::abs(a));
}

int fail(const std::string& msg)
{
    std::cerr << "FAIL: " << msg << '\n';
    return 1;
}

int run_self_test()
{
    {
        std::vector<Position> bad{{100.0, 0.0, 1.0, 0.2, 0.05, 1.0, OptionType::Call}};
        int rejected = 0;
        auto cleaned = dr3_example::exclude_invalid(bad, rejected);
        if (rejected != 1 || !cleaned.empty()) {
            return fail("zero strike was not excluded");
        }
        bad[0].strike = 100.0;
        bad[0].vol = -0.2;
        rejected = 0;
        cleaned = dr3_example::exclude_invalid(bad, rejected);
        if (rejected != 1) {
            return fail("negative volatility was not excluded");
        }
        try {
            dr3_example::black_scholes_call(100.0, 100.0, 1.0, 0.05, -0.2);
            return fail("negative volatility was not rejected by the scalar pricer");
        } catch (const std::invalid_argument&) {
        }
    }

    const int sizes[] = {3, 7, 128};
    for (int n : sizes) {
        auto book = dr3_example::make_portfolio(n, 7);
        auto packed = pack(book);
        auto scalar = scalar_book(book);
        auto simd_calls = simd_call_prices(packed);

        std::vector<double> simd_px(book.size());
        for (size_t i = 0; i < book.size(); ++i) {
            const double df = std::exp(-book[i].rate * book[i].maturity);
            const double call = simd_calls[i];
            simd_px[i] = book[i].type == OptionType::Call
                ? call
                : call - book[i].spot + book[i].strike * df;
        }

        const double abs_err = max_abs_diff(scalar.price, simd_px);
        const double rel_err = max_rel_diff(scalar.price, simd_px);
        if (abs_err > 1.0e-6 && rel_err > kPriceTol) {
            return fail("scalar vs SIMD price mismatch at n=" + std::to_string(n)
                + " abs=" + std::to_string(abs_err));
        }

        for (size_t i = 0; i < book.size(); ++i) {
            auto call = dr3_example::black_scholes_call(
                book[i].spot, book[i].strike, book[i].maturity, book[i].rate, book[i].vol);
            auto put = dr3_example::black_scholes_put(
                book[i].spot, book[i].strike, book[i].maturity, book[i].rate, book[i].vol);
            const double gap = dr3_example::put_call_parity_gap(
                call.price, put.price, book[i].spot, book[i].strike, book[i].maturity, book[i].rate);
            if (std::abs(gap) > kParityTol) {
                return fail("put-call parity failed");
            }
        }

        // Finite-difference Delta / Vega vs analytic.
        for (size_t i = 0; i < std::min<size_t>(book.size(), 8); ++i) {
            const auto& p = book[i];
            auto up = p;
            auto dn = p;
            up.spot += kFdBump;
            dn.spot -= kFdBump;
            const double fd_delta =
                (dr3_example::black_scholes_call(up.spot, p.strike, p.maturity, p.rate, p.vol).price
                    - dr3_example::black_scholes_call(dn.spot, p.strike, p.maturity, p.rate, p.vol).price)
                / (2.0 * kFdBump);
            auto analytic = dr3_example::black_scholes_call(p.spot, p.strike, p.maturity, p.rate, p.vol);
            if (!nearly_eq(fd_delta, analytic.delta, kGreekTol)) {
                return fail("FD Delta mismatch");
            }
            up = p;
            dn = p;
            up.vol += kFdBump;
            dn.vol -= kFdBump;
            const double fd_vega =
                (dr3_example::black_scholes_call(p.spot, p.strike, p.maturity, p.rate, up.vol).price
                    - dr3_example::black_scholes_call(p.spot, p.strike, p.maturity, p.rate, dn.vol).price)
                / (2.0 * kFdBump);
            if (!nearly_eq(fd_vega, analytic.vega, kGreekTol * 10.0)) {
                return fail("FD Vega mismatch");
            }
        }

        std::vector<double> weighted_simd(book.size());
        for (size_t i = 0; i < book.size(); ++i) {
            weighted_simd[i] = simd_px[i] * book[i].quantity;
        }
        const double naive = naive_sum(weighted_simd);
        const double pair = pairwise_sum(weighted_simd);
        if (std::abs(naive - pair) > 1.0e-8 * std::max(1.0, std::abs(naive))) {
            return fail("naive vs pairwise aggregation diverged unexpectedly");
        }

        auto again = dr3_example::make_portfolio(n, 7);
        if (again[0].spot != book[0].spot || again.back().quantity != book.back().quantity) {
            return fail("fixture was not deterministic");
        }
    }

    std::cout << "portfolioGreeksExample self-test passed\n";
    return 0;
}

void print_benchmark(const std::vector<Position>& book)
{
    auto packed = pack(book);
    auto scalar = scalar_book(book);

    auto t0 = std::chrono::steady_clock::now();
    volatile double sink = 0.0;
    const int repeats = 50;
    for (int i = 0; i < repeats; ++i) {
        sink += naive_sum(scalar_book(book).weighted);
    }
    auto t1 = std::chrono::steady_clock::now();
    for (int i = 0; i < repeats; ++i) {
        auto px = simd_call_prices(packed);
        sink += px.empty() ? 0.0 : px[0];
    }
    auto t2 = std::chrono::steady_clock::now();

    const double scalar_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count() / repeats;
    const double simd_ms =
        std::chrono::duration<double, std::milli>(t2 - t1).count() / repeats;
    const double n = static_cast<double>(book.size());
    auto simd_px = simd_call_prices(packed);
    std::vector<double> aligned(book.size());
    for (size_t i = 0; i < book.size(); ++i) {
        const double df = std::exp(-book[i].rate * book[i].maturity);
        aligned[i] = book[i].type == OptionType::Call
            ? simd_px[i]
            : simd_px[i] - book[i].spot + book[i].strike * df;
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Portfolio size: " << book.size() << '\n';
    std::cout << "Compiler: "
#if defined(__clang__)
              << "Clang " << __clang_major__
#elif defined(__GNUC__)
              << "GCC " << __GNUC__
#elif defined(_MSC_VER)
              << "MSVC " << _MSC_VER
#else
              << "unknown"
#endif
              << '\n';
    std::cout << "Build type: "
#ifdef NDEBUG
              << "Release\n";
#else
              << "Debug\n";
#endif
    std::cout << "SIMD type: AVX2 VecD4D (width " << VecXX::INS::size() << ")\n";
    std::cout << "Scalar elapsed ms: " << scalar_ms << '\n';
    std::cout << "Vector elapsed ms: " << simd_ms << '\n';
    std::cout << "Options priced per second (vector): "
              << (simd_ms > 0.0 ? n * 1000.0 / simd_ms : 0.0) << '\n';
    std::cout << "Maximum absolute pricing difference: "
              << max_abs_diff(scalar.price, aligned) << '\n';
    std::cout << "Maximum relative pricing difference: "
              << max_rel_diff(scalar.price, aligned) << '\n';
    std::cout << "Naive vs pairwise value gap: "
              << (naive_sum(scalar.weighted) - pairwise_sum(scalar.weighted)) << '\n';
    auto shocked = book;
    for (auto& p : shocked) {
        p.spot *= 1.01;
    }
    const double base_value = pairwise_sum(scalar.weighted);
    const double shocked_value = pairwise_sum(scalar_book(shocked).weighted);
    std::cout << "Base portfolio value: " << base_value << '\n';
    std::cout << "Repriced value after +1% spot: " << shocked_value << '\n';
    std::cout << "Portfolio P&L: " << (shocked_value - base_value) << '\n';
    std::cout << "(sink " << sink << ")\n";
    std::cout << "Times depend on processor, compiler, instruction set, and flags.\n";
}

}  // namespace

int main(int argc, char** argv)
{
    bool self_test = false;
    int size = 1024;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--self-test") {
            self_test = true;
        } else if (arg == "--size" && i + 1 < argc) {
            size = std::atoi(argv[++i]);
        } else if (arg == "--help") {
            std::cout << "portfolioGreeksExample [--self-test] [--size N]\n";
            return 0;
        }
    }
    if (self_test) {
        return run_self_test();
    }
    print_benchmark(dr3_example::make_portfolio(size, 1));
    return 0;
}
