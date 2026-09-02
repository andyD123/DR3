#include "pch.h"

#include "../Vectorisation/VecX/neural_kernels.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

namespace {

constexpr std::array<std::size_t, 7> kBoundarySizes = {1, 7, 8, 9, 15, 16, 17};
constexpr float kCanary = -12345.5f;

struct MisalignedBuffer {
    alignas(32) std::array<float, 32> storage;

    MisalignedBuffer()
    {
        storage.fill(kCanary);
    }

    float* data() { return storage.data() + 1; }
    const float* data() const { return storage.data() + 1; }
};

void copy_values(const std::vector<float>& source, MisalignedBuffer& destination)
{
    std::copy(source.begin(), source.end(), destination.data());
}

bool close(double actual, double expected, double abs_tol = 2e-5, double rel_tol = 2e-5)
{
    if (!std::isfinite(actual) || !std::isfinite(expected)) {
        return false;
    }
    const double diff = std::abs(actual - expected);
    return diff <= abs_tol ||
        diff <= rel_tol * std::max({1.0, std::abs(actual), std::abs(expected)});
}

double relative_error(double actual, double expected)
{
    return std::abs(actual - expected) / std::max(1e-30, std::abs(expected));
}

std::vector<float> deterministic_values(std::size_t n)
{
    static constexpr float pattern[] = {
        -8.f, -2.f, -0.5f, 0.f, 0.25f, 1.f, 3.f, 8.f, -1.f,
        4.5f, -6.25f, 2.75f, -3.5f, 0.125f, 6.f, -7.f, 1.5f};
    std::vector<float> values(n);
    for (std::size_t i = 0; i < n; ++i) {
        values[i] = pattern[i % (sizeof(pattern) / sizeof(pattern[0]))];
    }
    return values;
}

double log_sum_exp_ref(const float* x, std::size_t n)
{
    const double xmax = static_cast<double>(*std::max_element(x, x + n));
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sum += std::exp(static_cast<double>(x[i]) - xmax);
    }
    return xmax + std::log(sum);
}

std::vector<double> softmax_ref(const float* x, std::size_t n)
{
    const double lse = log_sum_exp_ref(x, n);
    std::vector<double> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = std::exp(static_cast<double>(x[i]) - lse);
    }
    return out;
}

std::vector<double> log_softmax_ref(const float* x, std::size_t n)
{
    const double lse = log_sum_exp_ref(x, n);
    std::vector<double> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = static_cast<double>(x[i]) - lse;
    }
    return out;
}

std::vector<double> layer_norm_ref(
    const float* x, std::size_t n, double epsilon,
    const float* scale = nullptr, const float* bias = nullptr)
{
    const double mean = std::accumulate(x, x + n, 0.0) / static_cast<double>(n);
    double variance = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double centered = static_cast<double>(x[i]) - mean;
        variance += centered * centered;
    }
    variance /= static_cast<double>(n);
    const double inv_std = 1.0 / std::sqrt(variance + epsilon);
    std::vector<double> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = (static_cast<double>(x[i]) - mean) * inv_std;
        if (scale != nullptr) {
            out[i] *= static_cast<double>(scale[i]);
        }
        if (bias != nullptr) {
            out[i] += static_cast<double>(bias[i]);
        }
    }
    return out;
}

std::vector<double> rms_norm_ref(
    const float* x, std::size_t n, double epsilon, const float* scale = nullptr)
{
    double mean_square = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double value = static_cast<double>(x[i]);
        mean_square += value * value;
    }
    mean_square /= static_cast<double>(n);
    const double inv_rms = 1.0 / std::sqrt(mean_square + epsilon);
    std::vector<double> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = static_cast<double>(x[i]) * inv_rms;
        if (scale != nullptr) {
            out[i] *= static_cast<double>(scale[i]);
        }
    }
    return out;
}

double activation_ref(float value, DRC::AI::Activation kind, double leaky_alpha = 0.01)
{
    const double x = static_cast<double>(value);
    switch (kind) {
    case DRC::AI::Activation::Relu:
        return std::max(x, 0.0);
    case DRC::AI::Activation::LeakyRelu:
        return x >= 0.0 ? x : leaky_alpha * x;
    case DRC::AI::Activation::Sigmoid:
        return 1.0 / (1.0 + std::exp(-x));
    case DRC::AI::Activation::Tanh:
        return std::tanh(x);
    case DRC::AI::Activation::Gelu:
        return 0.5 * x *
            (1.0 + std::tanh(0.7978845608028654 * (x + 0.044715 * x * x * x)));
    case DRC::AI::Activation::Silu:
        return x / (1.0 + std::exp(-x));
    }
    return std::numeric_limits<double>::quiet_NaN();
}

} // namespace

TEST(Softmax, UniformSingleAndTailsExcludePadding)
{
    for (std::size_t n : kBoundarySizes) {
        std::vector<float> x(n, -10000.f);
        std::vector<float> storage(n + 2, kCanary);
        DRC::AI::softmax(x.data(), n, storage.data() + 1);
        EXPECT_FLOAT_EQ(storage.front(), kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(storage.back(), kCanary) << "n=" << n;
        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_NEAR(storage[i + 1], 1.f / static_cast<float>(n), 2e-6f)
                << "n=" << n << " i=" << i;
        }
    }
}

TEST(Softmax, VecF8FMatchesIndependentDoubleAcrossTails)
{
    double max_softmax_abs_error = 0.0;
    double max_softmax_rel_error = 0.0;
    double max_log_softmax_abs_error = 0.0;
    double max_lse_abs_error = 0.0;

    for (std::size_t n : kBoundarySizes) {
        const auto logical = deterministic_values(n);
        MisalignedBuffer input_storage;
        MisalignedBuffer soft_storage;
        MisalignedBuffer log_storage;
        MisalignedBuffer soft_inplace;
        MisalignedBuffer log_inplace;
        copy_values(logical, input_storage);
        copy_values(logical, soft_inplace);
        copy_values(logical, log_inplace);
        const auto original_storage = input_storage.storage;
        const float* input = input_storage.data();
        float* soft = soft_storage.data();
        float* logs = log_storage.data();

        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(input) % 32, 4u);
        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(soft) % 32, 4u);
        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(soft_inplace.data()) % 32, 4u);

        DRC::AI::softmax(input, n, soft);
        DRC::AI::log_softmax(input, n, logs);
        DRC::AI::softmax_inplace(soft_inplace.data(), n);
        DRC::AI::log_softmax_inplace(log_inplace.data(), n);
        const float lse = DRC::AI::log_sum_exp(input, n);

        const auto soft_ref = softmax_ref(input, n);
        const auto logs_ref = log_softmax_ref(input, n);
        const double lse_ref = log_sum_exp_ref(input, n);
        double probability_sum = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            probability_sum += static_cast<double>(soft[i]);
            EXPECT_TRUE(close(soft[i], soft_ref[i])) << "softmax n=" << n << " i=" << i;
            EXPECT_TRUE(close(logs[i], logs_ref[i], 3e-5, 3e-5))
                << "log_softmax n=" << n << " i=" << i;
            EXPECT_NEAR(logs[i], std::log(soft[i]), 3e-5f)
                << "log(softmax) n=" << n << " i=" << i;
            EXPECT_FLOAT_EQ(soft[i], soft_inplace.data()[i]) << "n=" << n << " i=" << i;
            EXPECT_FLOAT_EQ(logs[i], log_inplace.data()[i]) << "n=" << n << " i=" << i;
            max_softmax_abs_error = std::max(
                max_softmax_abs_error, std::abs(static_cast<double>(soft[i]) - soft_ref[i]));
            max_softmax_rel_error = std::max(
                max_softmax_rel_error, relative_error(soft[i], soft_ref[i]));
            max_log_softmax_abs_error = std::max(
                max_log_softmax_abs_error, std::abs(static_cast<double>(logs[i]) - logs_ref[i]));
        }
        max_lse_abs_error = std::max(
            max_lse_abs_error, std::abs(static_cast<double>(lse) - lse_ref));
        EXPECT_TRUE(close(lse, lse_ref, 3e-5, 3e-5)) << "n=" << n;
        EXPECT_NEAR(probability_sum, 1.0, 2e-6) << "n=" << n;
        EXPECT_EQ(input_storage.storage, original_storage) << "n=" << n;
        EXPECT_FLOAT_EQ(soft_storage.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(soft_storage.storage[n + 1], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(log_storage.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(log_storage.storage[n + 1], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(soft_inplace.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(soft_inplace.storage[n + 1], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(log_inplace.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(log_inplace.storage[n + 1], kCanary) << "n=" << n;
    }

    RecordProperty("softmax_max_abs_error", max_softmax_abs_error);
    RecordProperty("softmax_max_rel_error", max_softmax_rel_error);
    RecordProperty("log_softmax_max_abs_error", max_log_softmax_abs_error);
    RecordProperty("log_sum_exp_max_abs_error", max_lse_abs_error);
    EXPECT_TRUE(std::isfinite(max_softmax_abs_error));
    EXPECT_TRUE(std::isfinite(max_softmax_rel_error));
    EXPECT_TRUE(std::isfinite(max_log_softmax_abs_error));
    EXPECT_TRUE(std::isfinite(max_lse_abs_error));
    EXPECT_LT(max_softmax_rel_error, 3e-5);
}

TEST(Softmax, ExtremeFiniteLogitsRemainStable)
{
    const std::vector<float> x = {
        10000.f, 9999.f, -10000.f, -9999.f, 0.f, 7.f, -4.f, 10001.f, 2.f};
    std::vector<float> soft(x.size()), logs(x.size());
    DRC::AI::softmax(x.data(), x.size(), soft.data());
    DRC::AI::log_softmax(x.data(), x.size(), logs.data());
    const auto soft_ref = softmax_ref(x.data(), x.size());
    const auto logs_ref = log_softmax_ref(x.data(), x.size());
    const double lse_ref = log_sum_exp_ref(x.data(), x.size());
    const float lse = DRC::AI::log_sum_exp(x.data(), x.size());
    double sum = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        EXPECT_TRUE(std::isfinite(soft[i]));
        EXPECT_TRUE(std::isfinite(logs[i]));
        EXPECT_TRUE(close(soft[i], soft_ref[i], 2e-5, 2e-5));
        EXPECT_NEAR(logs[i], logs_ref[i], 2e-3);
        sum += static_cast<double>(soft[i]);
    }
    EXPECT_TRUE(std::isfinite(lse));
    EXPECT_NEAR(lse, lse_ref, 2e-3);
    EXPECT_NEAR(sum, 1.0, 2e-6);
}

TEST(Softmax, EmptyAndNonFiniteInputsHaveExplicitBehavior)
{
    EXPECT_EQ(DRC::AI::log_sum_exp(nullptr, 0), -std::numeric_limits<float>::infinity());
    EXPECT_NO_THROW(DRC::AI::softmax(nullptr, 0, nullptr));
    EXPECT_NO_THROW(DRC::AI::log_softmax(nullptr, 0, nullptr));

    const std::array<float, 3> exceptional = {
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::infinity(),
        -std::numeric_limits<float>::infinity()};
    for (float bad : exceptional) {
        const float x[] = {0.f, bad, 1.f};
        float out[] = {kCanary, kCanary, kCanary};
        EXPECT_THROW(DRC::AI::softmax(x, 3, out), std::invalid_argument);
        EXPECT_THROW(DRC::AI::log_softmax(x, 3, out), std::invalid_argument);
        EXPECT_THROW(DRC::AI::log_sum_exp(x, 3), std::invalid_argument);
        EXPECT_FLOAT_EQ(out[0], kCanary);
        EXPECT_FLOAT_EQ(out[1], kCanary);
        EXPECT_FLOAT_EQ(out[2], kCanary);
    }
}

TEST(LayerNorm, VecF8FMatchesDoubleAcrossTailsAndAffineInPlace)
{
    double max_abs_error = 0.0;
    double max_rel_error = 0.0;
    for (std::size_t n : kBoundarySizes) {
        const auto logical = deterministic_values(n);
        std::vector<float> scale_values(n), bias_values(n);
        for (std::size_t i = 0; i < n; ++i) {
            scale_values[i] = 0.8f + static_cast<float>(i) * 0.03f;
            bias_values[i] = -0.2f + static_cast<float>(i) * 0.01f;
        }
        MisalignedBuffer input_storage;
        MisalignedBuffer scale;
        MisalignedBuffer bias;
        MisalignedBuffer output_storage;
        MisalignedBuffer inplace;
        copy_values(logical, input_storage);
        copy_values(scale_values, scale);
        copy_values(bias_values, bias);
        copy_values(logical, inplace);
        const auto original_storage = input_storage.storage;
        const float* input = input_storage.data();
        float* output = output_storage.data();
        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(scale.data()) % 32, 4u);
        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(bias.data()) % 32, 4u);
        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(inplace.data()) % 32, 4u);
        DRC::AI::layer_norm(input, n, 1e-5f, output, scale.data(), bias.data());
        DRC::AI::layer_norm_inplace(
            inplace.data(), n, 1e-5f, scale.data(), bias.data());
        const auto ref = layer_norm_ref(input, n, 1e-5, scale.data(), bias.data());
        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_TRUE(close(output[i], ref[i], 4e-5, 4e-5)) << "n=" << n << " i=" << i;
            EXPECT_FLOAT_EQ(output[i], inplace.data()[i]) << "n=" << n << " i=" << i;
            max_abs_error = std::max(
                max_abs_error, std::abs(static_cast<double>(output[i]) - ref[i]));
            max_rel_error = std::max(max_rel_error, relative_error(output[i], ref[i]));
        }
        EXPECT_EQ(input_storage.storage, original_storage) << "n=" << n;
        EXPECT_FLOAT_EQ(output_storage.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(output_storage.storage[n + 1], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(inplace.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(inplace.storage[n + 1], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(scale.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(scale.storage[n + 1], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(bias.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(bias.storage[n + 1], kCanary) << "n=" << n;
    }
    RecordProperty("layer_norm_max_abs_error", max_abs_error);
    RecordProperty("layer_norm_max_rel_error", max_rel_error);
    EXPECT_TRUE(std::isfinite(max_abs_error));
    EXPECT_TRUE(std::isfinite(max_rel_error));
}

TEST(LayerNorm, ZeroVarianceUsesEpsilonAndAffineBias)
{
    std::vector<float> x(9, 3.f), scale(x.size()), bias(x.size()), out(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        scale[i] = 1.f + static_cast<float>(i) * 0.1f;
        bias[i] = -0.4f + static_cast<float>(i) * 0.05f;
    }
    std::vector<float> inplace = x;
    DRC::AI::layer_norm(
        x.data(), x.size(), 1e-5f, out.data(), scale.data(), bias.data());
    DRC::AI::layer_norm_inplace(
        inplace.data(), inplace.size(), 1e-5f, scale.data(), bias.data());
    for (std::size_t i = 0; i < x.size(); ++i) {
        EXPECT_TRUE(std::isfinite(out[i]));
        EXPECT_FLOAT_EQ(out[i], bias[i]);
        EXPECT_FLOAT_EQ(out[i], inplace[i]);
    }
}

TEST(LayerNorm, PreservesMeanResidualAndAppliesEpsilonInsideSquareRoot)
{
    {
        const float x[] = {100000000.f, 100000008.f};
        float out[2] = {};
        DRC::AI::layer_norm(x, 2, 1e-5f, out);
        const auto ref = layer_norm_ref(x, 2, 1e-5);
        EXPECT_NEAR(out[0], ref[0], 2e-6);
        EXPECT_NEAR(out[1], ref[1], 2e-6);
        EXPECT_LT(out[0], -0.99f);
        EXPECT_GT(out[1], 0.99f);
    }
    {
        const float x[] = {1.f, 1.0002f};
        float out[2] = {};
        DRC::AI::layer_norm(x, 2, 1e-2f, out);
        const auto ref = layer_norm_ref(x, 2, 1e-2);
        EXPECT_NEAR(out[0], ref[0], 2e-6);
        EXPECT_NEAR(out[1], ref[1], 2e-6);
    }
}

TEST(LayerNorm, RejectsInvalidEpsilonAndNonFiniteArrays)
{
    const float x[] = {-1.f, 0.f, 1.f};
    float out[3] = {};
    const std::array<float, 4> invalid_epsilon = {
        0.f, -1.f, std::numeric_limits<float>::infinity(),
        std::numeric_limits<float>::quiet_NaN()};
    for (float epsilon : invalid_epsilon) {
        EXPECT_THROW(DRC::AI::layer_norm(x, 3, epsilon, out), std::invalid_argument);
    }
    const float bad_input[] = {0.f, std::numeric_limits<float>::quiet_NaN(), 1.f};
    const float bad_scale[] = {1.f, std::numeric_limits<float>::infinity(), 1.f};
    const float bad_bias[] = {0.f, -std::numeric_limits<float>::infinity(), 0.f};
    EXPECT_THROW(DRC::AI::layer_norm(bad_input, 3, 1e-5f, out), std::invalid_argument);
    EXPECT_THROW(DRC::AI::layer_norm(x, 3, 1e-5f, out, bad_scale), std::invalid_argument);
    EXPECT_THROW(DRC::AI::layer_norm(x, 3, 1e-5f, out, nullptr, bad_bias), std::invalid_argument);
    EXPECT_NO_THROW(DRC::AI::layer_norm(nullptr, 0, 1e-5f, nullptr));
}

TEST(RMSNorm, VecF8FMatchesDoubleAcrossTailsScaleAndInPlace)
{
    double max_abs_error = 0.0;
    double max_rel_error = 0.0;
    for (std::size_t n : kBoundarySizes) {
        const auto logical = deterministic_values(n);
        std::vector<float> scale_values(n);
        for (std::size_t i = 0; i < n; ++i) {
            scale_values[i] = 0.75f + static_cast<float>(i) * 0.025f;
        }
        MisalignedBuffer input_storage;
        MisalignedBuffer scale;
        MisalignedBuffer output_storage;
        MisalignedBuffer inplace;
        copy_values(logical, input_storage);
        copy_values(scale_values, scale);
        copy_values(logical, inplace);
        const auto original_storage = input_storage.storage;
        const float* input = input_storage.data();
        float* output = output_storage.data();
        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(input) % 32, 4u);
        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(output) % 32, 4u);
        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(scale.data()) % 32, 4u);
        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(inplace.data()) % 32, 4u);
        DRC::AI::rms_norm(input, n, 1e-5f, output, scale.data());
        DRC::AI::rms_norm_inplace(inplace.data(), n, 1e-5f, scale.data());
        const auto ref = rms_norm_ref(input, n, 1e-5, scale.data());
        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_TRUE(close(output[i], ref[i], 4e-5, 4e-5)) << "n=" << n << " i=" << i;
            EXPECT_FLOAT_EQ(output[i], inplace.data()[i]) << "n=" << n << " i=" << i;
            max_abs_error = std::max(
                max_abs_error, std::abs(static_cast<double>(output[i]) - ref[i]));
            max_rel_error = std::max(max_rel_error, relative_error(output[i], ref[i]));
        }
        EXPECT_EQ(input_storage.storage, original_storage) << "n=" << n;
        EXPECT_FLOAT_EQ(output_storage.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(output_storage.storage[n + 1], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(inplace.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(inplace.storage[n + 1], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(scale.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(scale.storage[n + 1], kCanary) << "n=" << n;
    }
    RecordProperty("rms_norm_max_abs_error", max_abs_error);
    RecordProperty("rms_norm_max_rel_error", max_rel_error);
    EXPECT_TRUE(std::isfinite(max_abs_error));
    EXPECT_TRUE(std::isfinite(max_rel_error));
}

TEST(RMSNorm, ZeroVectorUsesEpsilon)
{
    const std::vector<float> x(9, 0.f);
    const std::vector<float> scale(9, 2.f);
    std::vector<float> out(9, kCanary);
    DRC::AI::rms_norm(x.data(), x.size(), 1e-5f, out.data(), scale.data());
    for (float value : out) {
        EXPECT_FLOAT_EQ(value, 0.f);
    }
}

TEST(RMSNorm, AppliesEpsilonInsideSquareRoot)
{
    const float x[] = {1e-4f, -1e-4f};
    float out[2] = {};
    DRC::AI::rms_norm(x, 2, 1e-2f, out);
    const auto ref = rms_norm_ref(x, 2, 1e-2);
    EXPECT_NEAR(out[0], ref[0], 2e-7);
    EXPECT_NEAR(out[1], ref[1], 2e-7);
}

TEST(RMSNorm, RejectsInvalidEpsilonAndNonFiniteArrays)
{
    const float x[] = {-1.f, 0.f, 1.f};
    float out[3] = {};
    const std::array<float, 4> invalid_epsilon = {
        0.f, -1.f, std::numeric_limits<float>::infinity(),
        std::numeric_limits<float>::quiet_NaN()};
    for (float epsilon : invalid_epsilon) {
        EXPECT_THROW(DRC::AI::rms_norm(x, 3, epsilon, out), std::invalid_argument);
    }
    const float bad_input[] = {0.f, std::numeric_limits<float>::infinity(), 1.f};
    const float bad_scale[] = {1.f, std::numeric_limits<float>::quiet_NaN(), 1.f};
    EXPECT_THROW(DRC::AI::rms_norm(bad_input, 3, 1e-5f, out), std::invalid_argument);
    EXPECT_THROW(DRC::AI::rms_norm(x, 3, 1e-5f, out, bad_scale), std::invalid_argument);
    EXPECT_NO_THROW(DRC::AI::rms_norm(nullptr, 0, 1e-5f, nullptr));
}

TEST(Activation, VecF8FMatchesScalarReferenceAcrossTailsAndInPlace)
{
    const std::array<DRC::AI::Activation, 6> kinds = {
        DRC::AI::Activation::Relu, DRC::AI::Activation::LeakyRelu,
        DRC::AI::Activation::Sigmoid, DRC::AI::Activation::Tanh,
        DRC::AI::Activation::Gelu, DRC::AI::Activation::Silu};
    constexpr float alpha = 0.125f;
    double max_abs_error = 0.0;
    double max_rel_error = 0.0;
    for (std::size_t n : kBoundarySizes) {
        const auto logical = deterministic_values(n);
        MisalignedBuffer input_storage;
        copy_values(logical, input_storage);
        const auto original_storage = input_storage.storage;
        for (DRC::AI::Activation kind : kinds) {
            MisalignedBuffer output_storage;
            MisalignedBuffer inplace;
            copy_values(logical, inplace);
            ASSERT_EQ(reinterpret_cast<std::uintptr_t>(input_storage.data()) % 32, 4u);
            ASSERT_EQ(reinterpret_cast<std::uintptr_t>(output_storage.data()) % 32, 4u);
            ASSERT_EQ(reinterpret_cast<std::uintptr_t>(inplace.data()) % 32, 4u);
            DRC::AI::activation(
                input_storage.data(), n, output_storage.data(), kind, alpha);
            DRC::AI::activation_inplace(inplace.data(), n, kind, alpha);
            for (std::size_t i = 0; i < n; ++i) {
                const double ref = activation_ref(logical[i], kind, alpha);
                const double actual = output_storage.data()[i];
                EXPECT_TRUE(close(actual, ref, 4e-5, 4e-5))
                    << "n=" << n << " i=" << i << " kind=" << static_cast<int>(kind);
                EXPECT_FLOAT_EQ(output_storage.data()[i], inplace.data()[i])
                    << "n=" << n << " i=" << i << " kind=" << static_cast<int>(kind);
                max_abs_error = std::max(max_abs_error, std::abs(actual - ref));
                max_rel_error = std::max(max_rel_error, relative_error(actual, ref));
            }
            EXPECT_FLOAT_EQ(output_storage.storage[0], kCanary);
            EXPECT_FLOAT_EQ(output_storage.storage[n + 1], kCanary);
            EXPECT_FLOAT_EQ(inplace.storage[0], kCanary);
            EXPECT_FLOAT_EQ(inplace.storage[n + 1], kCanary);
        }
        EXPECT_EQ(input_storage.storage, original_storage) << "n=" << n;
    }
    RecordProperty("activation_max_abs_error", max_abs_error);
    RecordProperty("activation_max_rel_error", max_rel_error);
    EXPECT_TRUE(std::isfinite(max_abs_error));
    EXPECT_TRUE(std::isfinite(max_rel_error));
}

TEST(Activation, ClipAndBiasHelperMatchScalarReference)
{
    for (std::size_t n : kBoundarySizes) {
        const auto logical = deterministic_values(n);
        std::vector<float> bias_values(n);
        for (std::size_t i = 0; i < n; ++i) {
            bias_values[i] = 0.5f - static_cast<float>(i) * 0.015625f;
        }
        MisalignedBuffer input;
        MisalignedBuffer bias;
        MisalignedBuffer output;
        MisalignedBuffer clip_inplace;
        MisalignedBuffer bias_alias;
        copy_values(logical, input);
        copy_values(bias_values, bias);
        copy_values(logical, clip_inplace);
        copy_values(logical, bias_alias);
        const auto original_input = input.storage;

        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(bias.data()) % 32, 4u);
        ASSERT_EQ(reinterpret_cast<std::uintptr_t>(bias_alias.data()) % 32, 4u);
        DRC::AI::clip(input.data(), n, -1.f, 2.f, output.data());
        DRC::AI::clip_inplace(clip_inplace.data(), n, -1.f, 2.f);
        for (std::size_t i = 0; i < n; ++i) {
            const float ref = std::max(-1.f, std::min(2.f, logical[i]));
            EXPECT_FLOAT_EQ(output.data()[i], ref) << "clip n=" << n << " i=" << i;
            EXPECT_FLOAT_EQ(clip_inplace.data()[i], ref) << "clip inplace n=" << n << " i=" << i;
        }

        DRC::AI::bias_plus_activation(
            input.data(), bias.data(), n, output.data(), DRC::AI::Activation::Relu);
        DRC::AI::bias_plus_activation(
            bias_alias.data(), bias.data(), n, bias_alias.data(), DRC::AI::Activation::Relu);
        for (std::size_t i = 0; i < n; ++i) {
            const float ref = std::max(0.f, logical[i] + bias_values[i]);
            EXPECT_FLOAT_EQ(output.data()[i], ref) << "bias n=" << n << " i=" << i;
            EXPECT_FLOAT_EQ(bias_alias.data()[i], ref) << "bias alias n=" << n << " i=" << i;
        }
        EXPECT_EQ(input.storage, original_input) << "n=" << n;
        EXPECT_FLOAT_EQ(output.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(output.storage[n + 1], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(bias_alias.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(bias_alias.storage[n + 1], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(bias.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(bias.storage[n + 1], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(clip_inplace.storage[0], kCanary) << "n=" << n;
        EXPECT_FLOAT_EQ(clip_inplace.storage[n + 1], kCanary) << "n=" << n;
    }
}

TEST(Activation, RejectsNonFiniteInputsAndInvalidParameters)
{
    const float x[] = {-1.f, 0.f, 1.f};
    const float bias[] = {0.5f, 0.5f, 0.5f};
    float out[3] = {};
    const std::array<float, 3> invalid_alpha = {
        -0.1f, std::numeric_limits<float>::infinity(),
        std::numeric_limits<float>::quiet_NaN()};
    for (float alpha : invalid_alpha) {
        EXPECT_THROW(
            DRC::AI::activation(x, 3, out, DRC::AI::Activation::LeakyRelu, alpha),
            std::invalid_argument);
        EXPECT_THROW(DRC::AI::leaky_relu(x, 3, out, alpha), std::invalid_argument);
        EXPECT_THROW(
            DRC::AI::bias_plus_activation(
                x, bias, 3, out, DRC::AI::Activation::LeakyRelu, alpha),
            std::invalid_argument);
    }

    const float bad_input[] = {0.f, std::numeric_limits<float>::quiet_NaN(), 1.f};
    const float bad_bias[] = {0.f, std::numeric_limits<float>::infinity(), 0.f};
    EXPECT_THROW(DRC::AI::relu(bad_input, 3, out), std::invalid_argument);
    EXPECT_THROW(
        DRC::AI::bias_plus_activation(
            x, bad_bias, 3, out, DRC::AI::Activation::Relu),
        std::invalid_argument);
    EXPECT_THROW(DRC::AI::clip(bad_input, 3, -1.f, 1.f, out), std::invalid_argument);
    EXPECT_THROW(DRC::AI::clip(x, 3, 2.f, -1.f, out), std::invalid_argument);
    EXPECT_THROW(
        DRC::AI::clip(x, 3, -std::numeric_limits<float>::infinity(), 1.f, out),
        std::invalid_argument);
    EXPECT_THROW(
        DRC::AI::clip(
            x, 3, -1.f, std::numeric_limits<float>::quiet_NaN(), out),
        std::invalid_argument);

    EXPECT_NO_THROW(
        DRC::AI::activation(nullptr, 0, nullptr, DRC::AI::Activation::Relu));
    EXPECT_THROW(
        DRC::AI::activation(nullptr, 0, nullptr, DRC::AI::Activation::LeakyRelu, -1.f),
        std::invalid_argument);
    const auto unknown = static_cast<DRC::AI::Activation>(999);
    EXPECT_THROW(DRC::AI::activation(x, 3, out, unknown), std::invalid_argument);
    EXPECT_THROW(DRC::AI::activation(nullptr, 0, nullptr, unknown), std::invalid_argument);
}
