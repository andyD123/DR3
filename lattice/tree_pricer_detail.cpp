#include "tree_pricer_detail.h"

#include "../Vectorisation/VecX/alloc_policy_imp.h"
#include "../Vectorisation/VecX/dr3.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{

using namespace DRC::VecD4D;

AllAllocatorsGuard<double> latticeAllocators;

void requireProbability(double value, const char* name)
{
    if (!std::isfinite(value) || value < 0.0 || value > 1.0)
    {
        throw std::domain_error(std::string("invalid tree probability: ") + name);
    }
}

} // namespace

namespace dr3::lattice::detail
{

void validate(const VanillaOption& option, TreeConfig config)
{
    if (!std::isfinite(option.spot) || !std::isfinite(option.strike)
        || !std::isfinite(option.volatility) || !std::isfinite(option.rate)
        || !std::isfinite(option.dividendYield) || !std::isfinite(option.maturity))
    {
        throw std::invalid_argument("option inputs must be finite");
    }
    if (option.spot < 0.0)
    {
        throw std::invalid_argument("spot must be non-negative");
    }
    if (option.strike <= 0.0)
    {
        throw std::invalid_argument("strike must be positive");
    }
    if (option.volatility <= 0.0)
    {
        throw std::invalid_argument("volatility must be positive");
    }
    if (option.maturity < 0.0)
    {
        throw std::invalid_argument("maturity must be non-negative");
    }
    if (config.steps < 1)
    {
        throw std::invalid_argument("tree steps must be at least one");
    }
    if (config.steps > static_cast<std::size_t>(std::numeric_limits<int>::max() / 2 - 1))
    {
        throw std::invalid_argument("tree step count is too large");
    }
    if (option.type != OptionType::Call && option.type != OptionType::Put)
    {
        throw std::invalid_argument("option type is invalid");
    }
}

double intrinsic(const VanillaOption& option, double spot)
{
    return option.type == OptionType::Call
        ? std::max(spot - option.strike, 0.0)
        : std::max(option.strike - spot, 0.0);
}

double priceBinomial(const VanillaOption& option, TreeConfig config)
{
    using namespace DRC::VecD4D;

    const int steps = static_cast<int>(config.steps);
    const double dt = option.maturity / static_cast<double>(steps);
    const double up = std::exp(option.volatility * std::sqrt(dt));
    const double down = 1.0 / up;
    const double growth = std::exp((option.rate - option.dividendYield) * dt);
    const double probabilityUp = (growth - down) / (up - down);
    const double probabilityDown = 1.0 - probabilityUp;
    const double discount = std::exp(-option.rate * dt);

    requireProbability(probabilityUp, "up");
    requireProbability(probabilityDown, "down");

    VecXX prices(0.0, steps + 1);
    double stock = option.spot * std::pow(down, steps);
    const double ratio = up / down;
    for (int index = 0; index <= steps; ++index)
    {
        prices[index] = intrinsic(option, stock);
        stock *= ratio;
    }

    VecXX scratch = prices;
    BinomialSampler<VecXX::INS> sampler;
    const VecXX::INS pu(probabilityUp);
    const VecXX::INS pd(probabilityDown);
    const VecXX::INS disc(discount);
    auto rollback = [=](BinomialSampler<VecXX::INS>& values)
    {
        return disc * (values.X_1.value * pu + values.X_0.value * pd);
    };

    for (int remaining = steps; remaining > 0; --remaining)
    {
        assert(remaining + 1 <= prices.size());
        transform(prices, scratch, rollback, sampler, 0, remaining + 1);
        std::swap(prices, scratch);
    }
    return prices[0];
}

double priceTrinomial(const VanillaOption& option,
                      TreeConfig config,
                      bool american,
                      double barrier,
                      double rebate,
                      bool barrierEnabled)
{
    using namespace DRC::VecD4D;

    const int steps = static_cast<int>(config.steps);
    const int nodeCount = 2 * steps + 1;
    const double dt = option.maturity / static_cast<double>(steps);
    const double dx = option.volatility * std::sqrt(3.0 * dt);
    const double drift = option.rate - option.dividendYield
        - 0.5 * option.volatility * option.volatility;
    const double varianceTerm = (dt * option.volatility * option.volatility
                                 + drift * drift * dt * dt) / (dx * dx);
    const double driftTerm = drift * dt / dx;
    const double probabilityUp = 0.5 * (varianceTerm + driftTerm);
    const double probabilityDown = 0.5 * (varianceTerm - driftTerm);
    const double probabilityMiddle = 1.0 - varianceTerm;
    const double discount = std::exp(-option.rate * dt);

    requireProbability(probabilityUp, "up");
    requireProbability(probabilityMiddle, "middle");
    requireProbability(probabilityDown, "down");
    if (std::abs(probabilityUp + probabilityMiddle + probabilityDown - 1.0) > 1.0e-12)
    {
        throw std::domain_error("trinomial probabilities must sum to one");
    }

    VecXX stockPrices(0.0, nodeCount);
    VecXX prices(0.0, nodeCount);
    double stock = option.spot * std::exp(-steps * dx);
    const double stockStep = std::exp(dx);
    for (int index = 0; index < nodeCount; ++index)
    {
        stockPrices[index] = stock;
        prices[index] = barrierEnabled && stock >= barrier
            ? rebate
            : intrinsic(option, stock);
        stock *= stockStep;
    }

    VecXX scratch = prices;
    TrinomialSampler<VecXX::INS> sampler;
    const VecXX::INS pu(probabilityUp);
    const VecXX::INS pm(probabilityMiddle);
    const VecXX::INS pd(probabilityDown);
    const VecXX::INS disc(discount);
    auto rollback = [=](TrinomialSampler<VecXX::INS>& values)
    {
        return disc * (values.X_1.value * pu
                       + values.X_0.value * pm
                       + values.X_Minus_1.value * pd);
    };

    for (int level = 0; level < steps; ++level)
    {
        const int begin = level;
        const int end = nodeCount - level;
        assert(begin >= 0 && begin < end && end <= prices.size());
        transform(prices, scratch, rollback, sampler, begin, end);
        std::swap(prices, scratch);

        const int activeBegin = level + 1;
        const int activeEnd = nodeCount - level - 1;
        for (int index = activeBegin; index < activeEnd; ++index)
        {
            if (barrierEnabled && stockPrices[index] >= barrier)
            {
                prices[index] = rebate;
            }
            else if (american)
            {
                prices[index] = std::max(prices[index], intrinsic(option, stockPrices[index]));
            }
        }
    }
    return prices[steps];
}

} // namespace dr3::lattice::detail
