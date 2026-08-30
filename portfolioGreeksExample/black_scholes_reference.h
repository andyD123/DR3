#pragma once

#include <cmath>
#include <stdexcept>

namespace dr3_example {

inline double norm_pdf(double x)
{
    static constexpr double inv_sqrt_2pi = 0.39894228040143267794;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

inline double norm_cdf(double x)
{
    return 0.5 * (1.0 + std::erf(x * 0.70710678118654752440));
}

struct BlackScholesPoint {
    double price = 0.0;
    double delta = 0.0;
    double vega = 0.0;
    double theta = 0.0;
};

inline void reject_invalid(double S, double K, double t, double r, double sigma)
{
    if (!(S > 0.0) || !(K > 0.0) || !(t > 0.0) || !(sigma > 0.0) || !std::isfinite(r)) {
        throw std::invalid_argument(
            "Black-Scholes inputs require positive spot, strike, maturity and volatility");
    }
}

inline BlackScholesPoint black_scholes_call(double S, double K, double t, double r, double sigma)
{
    reject_invalid(S, K, t, r, sigma);
    const double sqrt_t = std::sqrt(t);
    const double sigma_sqrt_t = sigma * sqrt_t;
    const double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * t) / sigma_sqrt_t;
    const double d2 = d1 - sigma_sqrt_t;
    const double df = std::exp(-r * t);
    const double nd1 = norm_cdf(d1);
    const double nd2 = norm_cdf(d2);
    const double npd1 = norm_pdf(d1);
    BlackScholesPoint g;
    g.price = S * nd1 - K * df * nd2;
    g.delta = nd1;
    g.vega = S * npd1 * sqrt_t;
    g.theta = -S * npd1 * sigma / (2.0 * sqrt_t) - r * K * df * nd2;
    return g;
}

inline BlackScholesPoint black_scholes_put(double S, double K, double t, double r, double sigma)
{
    auto call = black_scholes_call(S, K, t, r, sigma);
    const double df = std::exp(-r * t);
    BlackScholesPoint g;
    g.price = call.price - S + K * df;
    g.delta = call.delta - 1.0;
    g.vega = call.vega;
    g.theta = call.theta + r * K * df;
    return g;
}

inline double put_call_parity_gap(double call, double put, double S, double K, double t, double r)
{
    return (call - put) - (S - K * std::exp(-r * t));
}

}  // namespace dr3_example
