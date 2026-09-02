/************************ neural_kernels.h **********************************
 * Stable float32 preprocessing kernels for small CPU inference workloads.
 * AVX2 (VecF8F) is the current implementation and tested path; callers own
 * all storage.
 *
 * Public contract:
 * - Buffers are contiguous one-dimensional float32 arrays; caller alignment
 *   is not required. Only the first n logical values participate, including
 *   when n is not a multiple of the VecF8F width.
 * - Inputs and optional affine arrays contain finite values. Epsilon is finite
 *   and positive; LeakyReLU alpha is finite and non-negative; clip bounds are
 *   finite and ordered.
 * - LayerNorm uses population variance and RMSNorm uses mean square. Both add
 *   epsilon before the square root.
 * - Exact input/output aliasing is supported. Other overlapping ranges,
 *   including overlap with affine arrays, are not supported.
 * - Empty operations are no-ops and may use null buffers; empty log-sum-exp
 *   returns negative infinity. Parameters are still validated.
 * - Sum-like reductions visit logical values in index order and accumulate in
 *   double. The operation order is deterministic, but floating results should
 *   be compared with a tolerance rather than for cross-platform bit identity.
 *****************************************************************************/
#pragma once

#include "dr3.h"
#include "target_name_space.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace DRC {
namespace AI {

namespace NeuralSimd = VecF8F;

enum class Activation {
    Relu,
    LeakyRelu,
    Sigmoid,
    Tanh,
    Gelu,
    Silu
};

inline void require_buffer(const float* p, std::size_t n, const char* what)
{
    if (p == nullptr && n != 0) {
        throw std::invalid_argument(what);
    }
}

namespace detail {

inline void require_finite_buffer(const float* p, std::size_t n, const char* what)
{
    require_buffer(p, n, what);
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(p[i])) {
            throw std::invalid_argument(what);
        }
    }
}

inline void require_finite_optional_buffer(const float* p, std::size_t n, const char* what)
{
    if (p != nullptr) {
        require_finite_buffer(p, n, what);
    }
}

inline void require_activation_parameters(Activation kind, float leaky_alpha)
{
    switch (kind) {
    case Activation::LeakyRelu:
        if (!std::isfinite(leaky_alpha) || leaky_alpha < 0.f) {
            throw std::invalid_argument("LeakyReLU alpha must be finite and non-negative");
        }
        return;
    case Activation::Relu:
    case Activation::Sigmoid:
    case Activation::Tanh:
    case Activation::Gelu:
    case Activation::Silu:
        return;
    }
    throw std::invalid_argument("unknown activation");
}

} // namespace detail

inline NeuralSimd::VecXX load_neural_vec(const float* p, std::size_t n)
{
    require_buffer(p, n, "null input buffer");
    if (n == 0) {
        return NeuralSimd::VecXX();
    }
    return NeuralSimd::VecXX(std::vector<float>(p, p + n));
}

inline void store_neural_vec(const NeuralSimd::VecXX& v, float* out, std::size_t n)
{
    require_buffer(out, n, "null output buffer");
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = v[static_cast<int>(i)];
    }
}

inline NeuralSimd::VecXX activate_vec(
    const NeuralSimd::VecXX& x, Activation activation, float leaky_alpha = 0.01f)
{
    detail::require_activation_parameters(activation, leaky_alpha);
    switch (activation) {
    case Activation::Relu:
        return max(x, 0.f);
    case Activation::LeakyRelu:
        return max(x, 0.f) + leaky_alpha * min(x, 0.f);
    case Activation::Sigmoid:
        return 1.f / (1.f + exp(-x));
    case Activation::Tanh: {
        auto op = [](auto z) { return tanh(z); };
        return transform(op, x);
    }
    case Activation::Gelu: {
        // Hendrycks/Gimpel tanh approximation used by many inference runtimes.
        auto op = [](auto z) {
            using Value = decltype(z);
            const Value inner = Value(0.7978845608028654f) *
                (z + Value(0.044715f) * z * z * z);
            return Value(0.5f) * z * (Value(1.f) + tanh(inner));
        };
        return transform(op, x);
    }
    case Activation::Silu:
        return x / (1.f + exp(-x));
    }
    throw std::invalid_argument("unknown activation");
}

inline void activation(
    const float* input, std::size_t n, float* output,
    Activation kind, float leaky_alpha = 0.01f)
{
    detail::require_finite_buffer(input, n, "activation: input must contain only finite values");
    require_buffer(output, n, "activation: null output");
    detail::require_activation_parameters(kind, leaky_alpha);
    if (n == 0) {
        return;
    }
    const auto x = load_neural_vec(input, n);
    store_neural_vec(activate_vec(x, kind, leaky_alpha), output, n);
}

inline void activation_inplace(
    float* values, std::size_t n, Activation kind, float leaky_alpha = 0.01f)
{
    activation(values, n, values, kind, leaky_alpha);
}

inline void relu(const float* x, std::size_t n, float* out)
{
    activation(x, n, out, Activation::Relu);
}

inline void leaky_relu(const float* x, std::size_t n, float* out, float alpha = 0.01f)
{
    activation(x, n, out, Activation::LeakyRelu, alpha);
}

inline void sigmoid(const float* x, std::size_t n, float* out)
{
    activation(x, n, out, Activation::Sigmoid);
}

inline void tanh_activation(const float* x, std::size_t n, float* out)
{
    activation(x, n, out, Activation::Tanh);
}

inline void gelu(const float* x, std::size_t n, float* out)
{
    activation(x, n, out, Activation::Gelu);
}

inline void silu(const float* x, std::size_t n, float* out)
{
    activation(x, n, out, Activation::Silu);
}

inline void clip(const float* x, std::size_t n, float lo, float hi, float* out)
{
    if (!std::isfinite(lo) || !std::isfinite(hi) || lo > hi) {
        throw std::invalid_argument("clip: bounds must be finite and ordered");
    }
    detail::require_finite_buffer(x, n, "clip: input must contain only finite values");
    require_buffer(out, n, "clip: null output");
    if (n == 0) {
        return;
    }
    store_neural_vec(max(min(load_neural_vec(x, n), hi), lo), out, n);
}

inline void clip_inplace(float* x, std::size_t n, float lo, float hi)
{
    clip(x, n, lo, hi, x);
}

inline void bias_plus_activation(
    const float* input, const float* bias, std::size_t n, float* output,
    Activation kind, float leaky_alpha = 0.01f)
{
    detail::require_finite_buffer(
        input, n, "bias_plus_activation: input must contain only finite values");
    detail::require_finite_buffer(
        bias, n, "bias_plus_activation: bias must contain only finite values");
    require_buffer(output, n, "bias_plus_activation: null output");
    detail::require_activation_parameters(kind, leaky_alpha);
    if (n == 0) {
        return;
    }
    const auto z = load_neural_vec(input, n) + load_neural_vec(bias, n);
    store_neural_vec(activate_vec(z, kind, leaky_alpha), output, n);
}

inline float log_sum_exp(const float* input, std::size_t n)
{
    detail::require_finite_buffer(input, n, "log_sum_exp: input must contain only finite values");
    if (n == 0) {
        return -std::numeric_limits<float>::infinity();
    }
    const float xmax = *std::max_element(input, input + n);
    const auto shifted_exp = exp(load_neural_vec(input, n) - xmax);
    // Deliberately exclude VecXX padding and preserve logical index order.
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sum += static_cast<double>(shifted_exp[static_cast<int>(i)]);
    }
    return xmax + static_cast<float>(std::log(sum));
}

inline void softmax(const float* input, std::size_t n, float* output)
{
    detail::require_finite_buffer(input, n, "softmax: input must contain only finite values");
    require_buffer(output, n, "softmax: null output");
    if (n == 0) {
        return;
    }
    const float xmax = *std::max_element(input, input + n);
    const auto numerators = exp(load_neural_vec(input, n) - xmax);
    // Deliberately exclude VecXX padding and preserve logical index order.
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sum += static_cast<double>(numerators[static_cast<int>(i)]);
    }
    store_neural_vec(numerators / static_cast<float>(sum), output, n);
}

inline void softmax_inplace(float* values, std::size_t n)
{
    softmax(values, n, values);
}

inline void log_softmax(const float* input, std::size_t n, float* output)
{
    detail::require_finite_buffer(input, n, "log_softmax: input must contain only finite values");
    require_buffer(output, n, "log_softmax: null output");
    if (n == 0) {
        return;
    }
    const float lse = log_sum_exp(input, n);
    store_neural_vec(load_neural_vec(input, n) - lse, output, n);
}

inline void log_softmax_inplace(float* values, std::size_t n)
{
    log_softmax(values, n, values);
}

inline void layer_norm(
    const float* input, std::size_t n, float epsilon, float* output,
    const float* scale = nullptr, const float* bias = nullptr)
{
    detail::require_finite_buffer(input, n, "layer_norm: input must contain only finite values");
    detail::require_finite_optional_buffer(
        scale, n, "layer_norm: scale must contain only finite values");
    detail::require_finite_optional_buffer(
        bias, n, "layer_norm: bias must contain only finite values");
    require_buffer(output, n, "layer_norm: null output");
    if (!std::isfinite(epsilon) || !(epsilon > 0.f)) {
        throw std::invalid_argument("layer_norm: epsilon must be finite and positive");
    }
    if (n == 0) {
        return;
    }
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sum += input[i];
    }
    const double mean = sum / static_cast<double>(n);
    double squared = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double d = static_cast<double>(input[i]) - mean;
        squared += d * d;
    }
    const float inv_std = static_cast<float>(1.0 /
        std::sqrt(squared / static_cast<double>(n) + epsilon));
    const float mean_hi = static_cast<float>(mean);
    const float mean_lo = static_cast<float>(mean - static_cast<double>(mean_hi));
    auto normalized = ((load_neural_vec(input, n) - mean_hi) - mean_lo) * inv_std;
    if (scale != nullptr) {
        normalized *= load_neural_vec(scale, n);
    }
    if (bias != nullptr) {
        normalized += load_neural_vec(bias, n);
    }
    store_neural_vec(normalized, output, n);
}

inline void layer_norm_inplace(
    float* values, std::size_t n, float epsilon,
    const float* scale = nullptr, const float* bias = nullptr)
{
    layer_norm(values, n, epsilon, values, scale, bias);
}

inline void rms_norm(
    const float* input, std::size_t n, float epsilon, float* output,
    const float* scale = nullptr)
{
    detail::require_finite_buffer(input, n, "rms_norm: input must contain only finite values");
    detail::require_finite_optional_buffer(
        scale, n, "rms_norm: scale must contain only finite values");
    require_buffer(output, n, "rms_norm: null output");
    if (!std::isfinite(epsilon) || !(epsilon > 0.f)) {
        throw std::invalid_argument("rms_norm: epsilon must be finite and positive");
    }
    if (n == 0) {
        return;
    }
    double squared = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double x = input[i];
        squared += x * x;
    }
    const float inv_rms = static_cast<float>(1.0 /
        std::sqrt(squared / static_cast<double>(n) + epsilon));
    auto normalized = load_neural_vec(input, n) * inv_rms;
    if (scale != nullptr) {
        normalized *= load_neural_vec(scale, n);
    }
    store_neural_vec(normalized, output, n);
}

inline void rms_norm_inplace(
    float* values, std::size_t n, float epsilon, const float* scale = nullptr)
{
    rms_norm(values, n, epsilon, values, scale);
}

} // namespace AI
} // namespace DRC
