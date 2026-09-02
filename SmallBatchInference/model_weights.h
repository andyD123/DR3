#pragma once

#include <array>

namespace dr3_mlp {
namespace fixture {

constexpr std::size_t kInput = 5;
constexpr std::size_t kHidden = 7;
constexpr std::size_t kOutput = 3;
constexpr std::size_t kBatch = 2;

// Row-major weights for a deterministic 5 -> 7 -> 3 LayerNorm/GELU MLP.
constexpr std::array<float, kHidden * kInput> kWeights1 = {{
    0.12f, -0.31f, 0.08f, 0.24f, -0.17f,
    -0.22f, 0.15f, 0.41f, -0.09f, 0.05f,
    0.33f, 0.07f, -0.28f, 0.19f, 0.11f,
    -0.14f, 0.26f, 0.09f, -0.37f, 0.21f,
    0.05f, -0.18f, 0.32f, 0.14f, -0.29f,
    0.27f, 0.03f, -0.11f, 0.35f, 0.16f,
    -0.08f, 0.29f, -0.24f, 0.06f, 0.31f,
}};

constexpr std::array<float, kHidden> kBias1 = {{
    0.03f, -0.02f, 0.07f, -0.05f, 0.01f, 0.04f, -0.06f,
}};

constexpr std::array<float, kOutput * kHidden> kWeights2 = {{
    0.21f, -0.13f, 0.08f, 0.32f, -0.17f, 0.11f, 0.05f,
    -0.16f, 0.27f, 0.19f, -0.04f, 0.23f, -0.09f, 0.14f,
    0.07f, 0.18f, -0.25f, 0.12f, 0.09f, 0.28f, -0.21f,
}};

constexpr std::array<float, kOutput> kBias2 = {{0.02f, -0.03f, 0.01f}};

constexpr std::array<float, kBatch * kInput> kInputs = {{
    0.25f, -0.50f, 1.25f, 0.75f, -1.00f,
    -0.80f, 0.40f, 0.10f, -0.30f, 1.10f,
}};

} // namespace fixture
} // namespace dr3_mlp
