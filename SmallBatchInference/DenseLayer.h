#pragma once

#include "../Vectorisation/VecX/similarity.h"

#include <cstddef>
#include <stdexcept>

namespace dr3_mlp {

struct DenseLayerView {
    const float* weights;
    const float* bias;
    std::size_t input_size;
    std::size_t output_size;
};

inline void validate_dense(
    const DenseLayerView& layer, const float* input, std::size_t input_size,
    float* output, std::size_t output_size)
{
    if (input_size != layer.input_size || output_size != layer.output_size) {
        throw std::invalid_argument("dense_forward: mismatched dimensions");
    }
    if ((input == nullptr && input_size != 0) ||
        (output == nullptr && output_size != 0) ||
        (layer.weights == nullptr && input_size * output_size != 0) ||
        (layer.bias == nullptr && output_size != 0)) {
        throw std::invalid_argument("dense_forward: null buffer");
    }
}

inline void dense_forward_simd(
    const DenseLayerView& layer, const float* input, std::size_t input_size,
    float* output, std::size_t output_size)
{
    validate_dense(layer, input, input_size, output, output_size);
    DRC::AI::Simd::SpanXX input_span(const_cast<float*>(input), input_size);
    for (std::size_t row = 0; row < output_size; ++row) {
        float* weights = const_cast<float*>(layer.weights + row * input_size);
        DRC::AI::Simd::SpanXX weight_span(weights, input_size);
        output[row] = DRC::AI::dot_product(input_span, weight_span) + layer.bias[row];
    }
}

inline void dense_forward_scalar(
    const DenseLayerView& layer, const float* input, std::size_t input_size,
    float* output, std::size_t output_size)
{
    validate_dense(layer, input, input_size, output, output_size);
    for (std::size_t row = 0; row < output_size; ++row) {
        double sum = layer.bias[row];
        const float* weights = layer.weights + row * input_size;
        for (std::size_t col = 0; col < input_size; ++col) {
            sum += static_cast<double>(weights[col]) * input[col];
        }
        output[row] = static_cast<float>(sum);
    }
}

} // namespace dr3_mlp
