// Deterministic small-batch MLP inference using DR3 float32 SIMD kernels.

#include "../Vectorisation/VecX/alloc_policy.h"
#include "../Vectorisation/VecX/neural_kernels.h"
#include "DenseLayer.h"
#include "expected_outputs.h"
#include "model_weights.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

AllAllocatorsGuard<float> allocGuard;

namespace {

constexpr float kEpsilon = 1.0e-5f;
constexpr double kSelfTestTolerance = 5.0e-5;
constexpr std::size_t kMaximumBatch = 16;
constexpr std::size_t kBenchmarkIterations = 200000;

using namespace dr3_mlp::fixture;

struct FixedModel {
    dr3_mlp::DenseLayerView first;
    dr3_mlp::DenseLayerView second;
};

enum class Path { Scalar, Avx2 };

struct ErrorSummary {
    double max_abs = 0.0;
    double max_rel = 0.0;
};

FixedModel load_fixed_model()
{
    return {
        {kWeights1.data(), kBias1.data(), kInput, kHidden},
        {kWeights2.data(), kBias2.data(), kHidden, kOutput},
    };
}

double gelu_scalar(double x)
{
    return 0.5 * x * (1.0 + std::tanh(
        0.7978845608028654 * (x + 0.044715 * x * x * x)));
}

void layer_norm_scalar(const float* input, std::size_t n, float* output)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sum += input[i];
    }
    const double mean = sum / static_cast<double>(n);
    double squared = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double centered = static_cast<double>(input[i]) - mean;
        squared += centered * centered;
    }
    const double inv_std = 1.0 /
        std::sqrt(squared / static_cast<double>(n) + kEpsilon);
    for (std::size_t i = 0; i < n; ++i) {
        output[i] = static_cast<float>(
            (static_cast<double>(input[i]) - mean) * inv_std);
    }
}

void predict_scalar(const FixedModel& model, const float* input, float* output)
{
    std::array<float, kInput> normalized{};
    std::array<float, kHidden> hidden{};
    std::array<float, kOutput> logits{};

    layer_norm_scalar(input, normalized.size(), normalized.data());
    dr3_mlp::dense_forward_scalar(model.first, normalized.data(), normalized.size(),
        hidden.data(), hidden.size());
    for (float& value : hidden) {
        value = static_cast<float>(gelu_scalar(value));
    }
    dr3_mlp::dense_forward_scalar(model.second, hidden.data(), hidden.size(),
        logits.data(), logits.size());

    const float max_logit = *std::max_element(logits.begin(), logits.end());
    double sum = 0.0;
    for (float value : logits) {
        sum += std::exp(static_cast<double>(value) - max_logit);
    }
    for (std::size_t i = 0; i < logits.size(); ++i) {
        output[i] = static_cast<float>(
            std::exp(static_cast<double>(logits[i]) - max_logit) / sum);
    }
}

void predict_avx2(const FixedModel& model, const float* input, float* output)
{
    std::array<float, kInput> normalized{};
    std::array<float, kHidden> hidden{};
    std::array<float, kOutput> logits{};

    DRC::AI::layer_norm(input, normalized.size(), kEpsilon, normalized.data());
    dr3_mlp::dense_forward_simd(model.first, normalized.data(), normalized.size(),
        hidden.data(), hidden.size());
    DRC::AI::activation_inplace(
        hidden.data(), hidden.size(), DRC::AI::Activation::Gelu);
    dr3_mlp::dense_forward_simd(model.second, hidden.data(), hidden.size(),
        logits.data(), logits.size());
    DRC::AI::softmax(logits.data(), logits.size(), output);
}

void predict_batch(
    Path path, const FixedModel& model, const float* inputs,
    std::size_t batch, float* outputs)
{
    if (batch == 0 || batch > kMaximumBatch) {
        throw std::invalid_argument("predict_batch: batch must be in [1, 16]");
    }
    for (std::size_t sample = 0; sample < batch; ++sample) {
        const float* input = inputs + sample * kInput;
        float* output = outputs + sample * kOutput;
        if (path == Path::Avx2) {
            predict_avx2(model, input, output);
        } else {
            predict_scalar(model, input, output);
        }
    }
}

ErrorSummary compare_outputs(const float* expected, const float* actual, std::size_t n)
{
    ErrorSummary errors;
    for (std::size_t i = 0; i < n; ++i) {
        const double abs_error = std::abs(
            static_cast<double>(actual[i]) - expected[i]);
        const double rel_error = abs_error /
            std::max(1.0e-12, std::abs(static_cast<double>(expected[i])));
        errors.max_abs = std::max(errors.max_abs, abs_error);
        errors.max_rel = std::max(errors.max_rel, rel_error);
    }
    return errors;
}

double checksum(const float* values, std::size_t n)
{
    double result = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        result += static_cast<double>(i + 1) * values[i];
    }
    return result;
}

int fail(const std::string& message)
{
    std::cerr << "FAIL: " << message << '\n';
    return 1;
}

int run_self_test()
{
    static_assert(kInput % 8 != 0 && kHidden % 8 != 0 && kOutput % 8 != 0,
        "fixture dimensions must exercise AVX2 tails");
    const FixedModel model = load_fixed_model();
    std::array<float, kBatch * kOutput> scalar{};
    std::array<float, kBatch * kOutput> avx2{};
    std::array<float, kBatch * kOutput> repeated{};

    predict_batch(Path::Scalar, model, kInputs.data(), kBatch, scalar.data());
    predict_batch(Path::Avx2, model, kInputs.data(), kBatch, avx2.data());
    predict_batch(Path::Avx2, model, kInputs.data(), kBatch, repeated.data());
    const ErrorSummary errors = compare_outputs(scalar.data(), avx2.data(), avx2.size());
    if (errors.max_abs > kSelfTestTolerance) {
        return fail("AVX2 output differs from scalar-double reference");
    }
    if (compare_outputs(kExpectedOutputs.data(), avx2.data(), avx2.size()).max_abs >
        kSelfTestTolerance) {
        return fail("AVX2 output differs from committed fixture");
    }
    if (avx2 != repeated) {
        return fail("AVX2 output is not deterministic");
    }

    for (std::size_t sample = 0; sample < kBatch; ++sample) {
        double probability_sum = 0.0;
        for (std::size_t i = 0; i < kOutput; ++i) {
            probability_sum += avx2[sample * kOutput + i];
        }
        if (std::abs(probability_sum - 1.0) > 2.0e-6) {
            return fail("softmax probabilities do not sum to one");
        }
    }

    try {
        std::array<float, kHidden> output{};
        dr3_mlp::dense_forward_simd(model.first, kInputs.data(), kInput - 1,
            output.data(), output.size());
        return fail("dense layer accepted a mismatched input length");
    } catch (const std::invalid_argument&) {
    }
    try {
        predict_batch(Path::Avx2, model, kInputs.data(), 0, avx2.data());
        return fail("zero batch was accepted");
    } catch (const std::invalid_argument&) {
    }
    try {
        predict_batch(Path::Avx2, model, kInputs.data(), kMaximumBatch + 1, avx2.data());
        return fail("batch above 16 was accepted");
    } catch (const std::invalid_argument&) {
    }

    std::array<float, kMaximumBatch * kInput> batch_inputs{};
    std::array<float, kMaximumBatch * kOutput> batch_scalar{};
    std::array<float, kMaximumBatch * kOutput> batch_avx2{};
    for (std::size_t sample = 0; sample < kMaximumBatch; ++sample) {
        std::copy_n(kInputs.data() + (sample % kBatch) * kInput, kInput,
            batch_inputs.data() + sample * kInput);
    }
    for (std::size_t batch = 1; batch <= kMaximumBatch; ++batch) {
        predict_batch(Path::Scalar, model, batch_inputs.data(), batch, batch_scalar.data());
        predict_batch(Path::Avx2, model, batch_inputs.data(), batch, batch_avx2.data());
        if (compare_outputs(batch_scalar.data(), batch_avx2.data(), batch * kOutput).max_abs >
            kSelfTestTolerance) {
            return fail("scalar/AVX2 agreement failed for batch " + std::to_string(batch));
        }
    }

    std::cout << "SmallBatchInference self-test passed: max_abs=" << errors.max_abs
              << " max_rel=" << errors.max_rel
              << " tolerance=" << kSelfTestTolerance << '\n';
    return 0;
}

void emit_expected()
{
    const FixedModel model = load_fixed_model();
    std::array<float, kBatch * kOutput> output{};
    predict_batch(Path::Scalar, model, kInputs.data(), kBatch, output.data());
    std::cout << std::setprecision(9);
    for (float value : output) {
        std::cout << value << "f,\n";
    }
}

const char* compiler_name()
{
#if defined(_MSC_VER)
    return "MSVC";
#elif defined(__clang__)
    return "Clang";
#elif defined(__GNUC__)
    return "GCC";
#else
    return "unknown";
#endif
}

const char* build_type()
{
#ifdef NDEBUG
    return "Release";
#else
    return "Debug";
#endif
}

struct Measurement {
    double microseconds_per_batch;
    double output_checksum;
};

Measurement measure(Path path, const FixedModel& model)
{
    std::array<float, kBatch * kOutput> output{};
    for (std::size_t i = 0; i < 1000; ++i) {
        predict_batch(path, model, kInputs.data(), kBatch, output.data());
    }

    volatile float guard = 0.f;
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < kBenchmarkIterations; ++i) {
        predict_batch(path, model, kInputs.data(), kBatch, output.data());
        guard += output[i % output.size()];
    }
    const auto elapsed = std::chrono::duration<double, std::micro>(
        std::chrono::steady_clock::now() - start).count();
    if (!std::isfinite(guard)) {
        throw std::runtime_error("benchmark guard is not finite");
    }
    return {elapsed / static_cast<double>(kBenchmarkIterations),
        checksum(output.data(), output.size())};
}

void run_benchmark()
{
    const FixedModel model = load_fixed_model();
    std::array<float, kBatch * kOutput> scalar_output{};
    std::array<float, kBatch * kOutput> avx2_output{};
    predict_batch(Path::Scalar, model, kInputs.data(), kBatch, scalar_output.data());
    predict_batch(Path::Avx2, model, kInputs.data(), kBatch, avx2_output.data());
    const ErrorSummary errors = compare_outputs(
        scalar_output.data(), avx2_output.data(), avx2_output.size());
    const Measurement scalar = measure(Path::Scalar, model);
    const Measurement avx2 = measure(Path::Avx2, model);

    std::cout << std::fixed << std::setprecision(6)
              << "model=" << kInput << "->" << kHidden << "->" << kOutput
              << " batch=" << kBatch
              << " compiler=" << compiler_name()
              << " INSTRSET=" << INSTRSET
              << " build=" << build_type()
              << " iterations=" << kBenchmarkIterations << '\n'
              << "path=scalar_double latency_us_per_batch="
              << scalar.microseconds_per_batch
              << " checksum=" << scalar.output_checksum << '\n'
              << "path=avx2_vecf8f latency_us_per_batch="
              << avx2.microseconds_per_batch
              << " checksum=" << avx2.output_checksum << '\n'
              << "agreement max_abs=" << errors.max_abs
              << " max_rel=" << errors.max_rel << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    bool self_test = false;
    bool emit_fixture = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--self-test") {
            self_test = true;
        } else if (arg == "--emit-expected") {
            emit_fixture = true;
        } else {
            std::cerr << "Unknown option: " << arg << '\n';
            return 2;
        }
    }
    if (emit_fixture) {
        emit_expected();
        return 0;
    }
    if (self_test) {
        return run_self_test();
    }
    run_benchmark();
    return 0;
}
