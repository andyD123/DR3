#pragma once

namespace DRC {
namespace Curves {

enum class Interpolation {
    Linear,
    Flat
};

enum class Extrapolation {
    Flat,
    Throw
};

} // namespace Curves
} // namespace DRC
