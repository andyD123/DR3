// Vectorised Monte Carlo European call with inverse-normal uniforms,
// antithetic variates, and a 95% confidence interval. Demonstration only.

#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/target_name_space.h"
#include "../portfolioGreeksExample/black_scholes_reference.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace DRC::VecD4D;

namespace {

AllAllocatorsGuard<typename VecXX::SCALA_TYPE> allocGuard;

constexpr double kPi = 3.14159265358979323846;

struct McResult {
    double mean = 0.0;
    double standard_error = 0.0;
    double ci_low = 0.0;
    double ci_high = 0.0;
    int paths = 0;
};

std::uint32_t lcg(std::uint32_t& state)
{
    state = state * 1664525u + 1013904223u;
    return state;
}

double open_unit(std::uint32_t& state)
{
    // (0, 1) exclusive so cdfnorminv stays finite.
    const double u = ((lcg(state) >> 8) + 0.5) * (1.0 / 16777216.0);
    return std::min(1.0 - 1.0e-12, std::max(1.0e-12, u));
}

std::vector<double> uniforms(int n, std::uint32_t seed)
{
    std::vector<double> u(static_cast<size_t>(n));
    std::uint32_t state = seed;
    for (int i = 0; i < n; ++i) {
        u[static_cast<size_t>(i)] = open_unit(state);
    }
    return u;
}

double payoff(double S, double K, double t, double r, double sigma, double z)
{
    const double drift = (r - 0.5 * sigma * sigma) * t;
    const double vol = sigma * std::sqrt(t);
    const double ST = S * std::exp(drift + vol * z);
    const double disc = std::exp(-r * t);
    return disc * std::max(ST - K, 0.0);
}

McResult summarise(const std::vector<double>& samples)
{
    McResult r;
    r.paths = static_cast<int>(samples.size());
    const double n = static_cast<double>(samples.size());
    r.mean = std::accumulate(samples.begin(), samples.end(), 0.0) / n;
    double var = 0.0;
    for (double x : samples) {
        const double d = x - r.mean;
        var += d * d;
    }
    var /= std::max(1.0, n - 1.0);
    r.standard_error = std::sqrt(var / n);
    r.ci_low = r.mean - 1.96 * r.standard_error;
    r.ci_high = r.mean + 1.96 * r.standard_error;
    return r;
}

McResult mc_scalar(double S, double K, double t, double r, double sigma, int paths,
    std::uint32_t seed, bool antithetic)
{
    auto u = uniforms(paths, seed);
    std::vector<double> samples;
    samples.reserve(static_cast<size_t>(antithetic ? paths : paths));
    for (int i = 0; i < paths; ++i) {
        const double z = cdfnorminv(u[static_cast<size_t>(i)]);
        const double p = payoff(S, K, t, r, sigma, z);
        if (antithetic) {
            samples.push_back(0.5 * (p + payoff(S, K, t, r, sigma, -z)));
        } else {
            samples.push_back(p);
        }
    }
    return summarise(samples);
}

McResult mc_simd(double S, double K, double t, double r, double sigma, int paths,
    std::uint32_t seed, bool antithetic)
{
    auto u = uniforms(paths, seed);
    VecXX uv(u);
    VecXX z = cdfnorminv(uv);
    std::vector<double> zv(z);
    std::vector<double> samples;
    samples.reserve(static_cast<size_t>(paths));
    for (int i = 0; i < paths; ++i) {
        const double zi = zv[static_cast<size_t>(i)];
        const double p = payoff(S, K, t, r, sigma, zi);
        if (antithetic) {
            samples.push_back(0.5 * (p + payoff(S, K, t, r, sigma, -zi)));
        } else {
            samples.push_back(p);
        }
    }
    return summarise(samples);
}

int fail(const std::string& msg)
{
    std::cerr << "FAIL: " << msg << '\n';
    return 1;
}

int run_self_test()
{
    const double S = 100.0, K = 100.0, t = 1.0, r = 0.05, sigma = 0.20;
    const auto analytic = dr3_example::black_scholes_call(S, K, t, r, sigma);

    auto a = mc_scalar(S, K, t, r, sigma, 17, 42, true);
    auto b = mc_scalar(S, K, t, r, sigma, 17, 42, true);
    if (a.mean != b.mean) {
        return fail("scalar MC was not deterministic");
    }

    auto simd = mc_simd(S, K, t, r, sigma, 17, 42, true);
    if (std::abs(simd.mean - a.mean) > 1.0e-10) {
        return fail("SIMD and scalar means diverged on identical uniforms");
    }
    if (17 % static_cast<int>(VecXX::INS::size()) == 0) {
        return fail("expected a path count that is not a SIMD multiple");
    }

    auto small = mc_simd(S, K, t, r, sigma, 512, 7, true);
    auto large = mc_simd(S, K, t, r, sigma, 4096, 7, true);
    if (!(large.standard_error < small.standard_error)) {
        return fail("standard error did not narrow with more paths");
    }

    auto fat = mc_simd(S, K, t, r, sigma, 65536, 1, true);
    const bool inside = analytic.price >= fat.ci_low && analytic.price <= fat.ci_high;
    const bool close = std::abs(fat.mean - analytic.price) < 0.15;
    if (!inside && !close) {
        return fail("analytic price missed the CI and the absolute error bound");
    }

    std::cout << "monteCarloExample self-test passed\n";
    return 0;
}

void print_report(int paths)
{
    const double S = 100.0, K = 100.0, t = 1.0, r = 0.05, sigma = 0.20;
    const auto analytic = dr3_example::black_scholes_call(S, K, t, r, sigma);
    auto scalar = mc_scalar(S, K, t, r, sigma, paths, 1, true);
    auto simd = mc_simd(S, K, t, r, sigma, paths, 1, true);
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Paths: " << paths << '\n';
    std::cout << "SIMD type: AVX2 VecD4D (width " << VecXX::INS::size() << ")\n";
    std::cout << "Analytic BS: " << analytic.price << '\n';
    std::cout << "Scalar mean: " << scalar.mean << "  se=" << scalar.standard_error
              << "  CI=[" << scalar.ci_low << ", " << scalar.ci_high << "]\n";
    std::cout << "SIMD mean:   " << simd.mean << "  se=" << simd.standard_error
              << "  CI=[" << simd.ci_low << ", " << simd.ci_high << "]\n";
    std::cout << "Times and errors depend on processor, compiler, and flags.\n";
}

}  // namespace

int main(int argc, char** argv)
{
    bool self_test = false;
    int paths = 65536;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--self-test") {
            self_test = true;
        } else if (arg == "--paths" && i + 1 < argc) {
            paths = std::atoi(argv[++i]);
        }
    }
    if (self_test) {
        return run_self_test();
    }
    print_report(paths);
    return 0;
}
