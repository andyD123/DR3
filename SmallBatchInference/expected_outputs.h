#pragma once

#include "model_weights.h"

#include <array>

namespace dr3_mlp {
namespace fixture {

// Generated once from the scalar-double reference in SmallBatchInference.cpp.
// Values are probabilities in row-major batch order.
constexpr std::array<float, kBatch * kOutput> kExpectedOutputs = {{
    0.280332625f, 0.339604199f, 0.380063176f,
    0.362993747f, 0.327188492f, 0.309817761f,
}};

} // namespace fixture
} // namespace dr3_mlp
