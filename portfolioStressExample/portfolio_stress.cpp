// Vectorised portfolio stress and P&L aggregation for DR3.
// Filters, reductions, and naive / pairwise / Kahan summation on a
// deterministic book. Demonstration only.

#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/target_name_space.h"
#include "../portfolioGreeksExample/black_scholes_reference.h"
#include "../portfolioGreeksExample/portfolio_fixture.h"

#include <algorithm>
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

AllAllocatorsGuard<typename VecXX::SCALA_TYPE> allocGuard;

constexpr int kBuckets = 4;
constexpr double kLossThreshold = 5.0;

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

double kahan_sum(const std::vector<double>& x)
{
    double sum = 0.0;
    double c = 0.0;
    for (double v : x) {
        double y = v - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

double unit_price(const Position& p)
{
    return (p.type == OptionType::Call)
        ? dr3_example::black_scholes_call(p.spot, p.strike, p.maturity, p.rate, p.vol).price
        : dr3_example::black_scholes_put(p.spot, p.strike, p.maturity, p.rate, p.vol).price;
}

std::vector<double> position_values(const std::vector<Position>& book)
{
    std::vector<double> v(book.size());
    for (size_t i = 0; i < book.size(); ++i) {
        v[i] = unit_price(book[i]) * book[i].quantity;
    }
    return v;
}

std::vector<Position> shock_spot(const std::vector<Position>& book, double factor)
{
    auto out = book;
    for (auto& p : out) {
        p.spot *= factor;
    }
    return out;
}

std::vector<Position> shock_vol(const std::vector<Position>& book, double factor)
{
    auto out = book;
    for (auto& p : out) {
        p.vol *= factor;
    }
    return out;
}

std::vector<double> pnl(const std::vector<double>& base, const std::vector<double>& shocked)
{
    std::vector<double> out(base.size());
    for (size_t i = 0; i < base.size(); ++i) {
        out[i] = shocked[i] - base[i];
    }
    return out;
}

std::vector<double> filter_losses_scalar(const std::vector<double>& x, double threshold)
{
    std::vector<double> out;
    for (double v : x) {
        if (v < -threshold) {
            out.push_back(v);
        }
    }
    return out;
}

std::vector<double> filter_losses_dr3(const std::vector<double>& x, double threshold)
{
    VecXX data(x);
    auto below = [=](auto v) { return v < VecXX::scalar(-threshold); };
    auto view = ApplyFilter(below, data);
    std::vector<double> out(static_cast<size_t>(view.size()));
    for (int i = 0; i < view.size(); ++i) {
        out[static_cast<size_t>(i)] = view[i];
    }
    return out;
}

std::vector<double> bucket_totals(const std::vector<double>& x, int buckets)
{
    std::vector<double> tot(static_cast<size_t>(buckets), 0.0);
    for (size_t i = 0; i < x.size(); ++i) {
        tot[i % static_cast<size_t>(buckets)] += x[i];
    }
    return tot;
}

int fail(const std::string& msg)
{
    std::cerr << "FAIL: " << msg << '\n';
    return 1;
}

bool vec_close(const std::vector<double>& a, const std::vector<double>& b, double tol)
{
    if (a.size() != b.size()) {
        return false;
    }
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::abs(a[i] - b[i]) > tol * std::max(1.0, std::abs(a[i]))) {
            return false;
        }
    }
    return true;
}

int run_self_test()
{
    {
        std::vector<Position> bad{{100.0, 0.0, 1.0, 0.2, 0.01, 1.0, OptionType::Call}};
        int rejected = 0;
        auto cleaned = dr3_example::exclude_invalid(bad, rejected);
        if (rejected != 1 || !cleaned.empty()) {
            return fail("invalid strike not excluded");
        }
    }

    auto book = dr3_example::make_portfolio(17, 3);
    auto base = position_values(book);
    auto down10 = position_values(shock_spot(book, 0.90));
    auto p = pnl(base, down10);

    if (naive_sum(p) == 0.0) {
        return fail("spot -10% produced no P&L");
    }

    auto filtered_s = filter_losses_scalar(p, kLossThreshold);
    auto filtered_v = filter_losses_dr3(p, kLossThreshold);
    if (!vec_close(filtered_s, filtered_v, 1.0e-12)) {
        return fail("DR3 loss filter did not match scalar predicate");
    }

    auto buckets = bucket_totals(p, kBuckets);
    double bucket_sum = naive_sum(buckets);
    if (std::abs(bucket_sum - naive_sum(p)) > 1.0e-10) {
        return fail("bucket totals do not reconstruct book P&L");
    }

    const double n = naive_sum(p);
    const double pair = pairwise_sum(p);
    const double kahan = kahan_sum(p);
    if (std::abs(n - pair) > 1.0e-10 || std::abs(n - kahan) > 1.0e-10) {
        return fail("naive/pairwise/Kahan disagree on the fixture");
    }

    auto vol_up = position_values(shock_vol(book, 1.20));
    if (naive_sum(pnl(base, vol_up)) == 0.0) {
        return fail("vol shock produced no P&L");
    }

    auto again = dr3_example::make_portfolio(17, 3);
    if (again[0].spot != book[0].spot) {
        return fail("fixture was not deterministic");
    }

    std::cout << "portfolioStressExample self-test passed\n";
    return 0;
}

void print_report(const std::vector<Position>& book)
{
    auto base = position_values(book);
    const double shocks[] = {0.90, 0.95, 1.05, 1.10};
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Portfolio size: " << book.size() << '\n';
    std::cout << "Loss threshold: " << kLossThreshold << '\n';
    for (double f : shocks) {
        auto p = pnl(base, position_values(shock_spot(book, f)));
        auto losses = filter_losses_dr3(p, kLossThreshold);
        std::cout << "spot x" << f
                  << "  PnL=" << pairwise_sum(p)
                  << "  loss_count=" << losses.size()
                  << "  loss_sum=" << kahan_sum(losses)
                  << "  naive-pairwise=" << (naive_sum(p) - pairwise_sum(p))
                  << "  naive-kahan=" << (naive_sum(p) - kahan_sum(p))
                  << '\n';
    }
    auto pvol = pnl(base, position_values(shock_vol(book, 1.25)));
    std::cout << "vol x1.25 PnL=" << pairwise_sum(pvol) << '\n';
    std::cout << "Times and residuals depend on processor, compiler, and flags.\n";
}

}  // namespace

int main(int argc, char** argv)
{
    bool self_test = false;
    int size = 256;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--self-test") {
            self_test = true;
        } else if (arg == "--size" && i + 1 < argc) {
            size = std::atoi(argv[++i]);
        }
    }
    if (self_test) {
        return run_self_test();
    }
    print_report(dr3_example::make_portfolio(size, 1));
    return 0;
}
