#pragma once

#include "interpolation.h"
#include "../VecX/dr3.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace DRC {
namespace Curves {
namespace detail {

inline bool finite_value(double value)
{
    return std::isfinite(value);
}

template <typename INS_VEC>
bool finite_value(const Vec<INS_VEC>& value)
{
    if (value.isScalar()) {
        return std::isfinite(static_cast<double>(value.getScalarValue()));
    }
    for (int i = 0; i < value.size(); ++i) {
        if (!std::isfinite(static_cast<double>(value[static_cast<std::size_t>(i)]))) {
            return false;
        }
    }
    return true;
}

template <typename INS_VEC>
bool finite_value(const VecD<INS_VEC>& value)
{
    return finite_value(value.value()) && finite_value(value.derivative());
}

inline bool same_shape(double, double)
{
    return true;
}

template <typename INS_VEC>
bool same_shape(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
    return lhs.isScalar() == rhs.isScalar() &&
        (lhs.isScalar() || lhs.size() == rhs.size());
}

template <typename INS_VEC>
bool same_shape(const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
    return same_shape(lhs.value(), rhs.value()) &&
        same_shape(lhs.derivative(), rhs.derivative());
}

} // namespace detail

template <typename Value>
class Curve {
public:
    Curve(
        std::vector<double> pillars,
        std::vector<Value> values,
        Interpolation interpolation,
        Extrapolation left_extrapolation,
        Extrapolation right_extrapolation)
        : interpolation_(interpolation),
          left_extrapolation_(left_extrapolation),
          right_extrapolation_(right_extrapolation)
    {
        reset(std::move(pillars), std::move(values));
    }

    void reset(std::vector<double> pillars, std::vector<Value> values)
    {
        validate(pillars, values);
        pillars_ = std::move(pillars);
        values_ = std::move(values);
    }

    Value value_at(double query) const
    {
        if (!std::isfinite(query)) {
            throw std::invalid_argument("Curve: query must be finite");
        }
        if (query < pillars_.front()) {
            return extrapolate(query, true);
        }
        if (query > pillars_.back()) {
            return extrapolate(query, false);
        }

        const auto upper = std::lower_bound(pillars_.begin(), pillars_.end(), query);
        const std::size_t upper_index = static_cast<std::size_t>(
            std::distance(pillars_.begin(), upper));
        if (upper != pillars_.end() && *upper == query) {
            return values_[upper_index];
        }
        // query is inside the pillar range and not exact, so both neighbours exist.
        const std::size_t lower_index = upper_index - 1;
        if (interpolation_ == Interpolation::Flat) {
            return values_[lower_index];
        }
        const double fraction = (query - pillars_[lower_index]) /
            (pillars_[upper_index] - pillars_[lower_index]);
        return values_[lower_index] +
            (values_[upper_index] - values_[lower_index]) * fraction;
    }

    void evaluate_sorted(
        const double* queries, std::size_t query_count,
        Value* output, std::size_t output_capacity) const
    {
        if ((queries == nullptr && query_count != 0) ||
            (output == nullptr && query_count != 0)) {
            throw std::invalid_argument("Curve: null bulk buffer");
        }
        if (output_capacity < query_count) {
            throw std::invalid_argument("Curve: bulk output capacity is too small");
        }
        if (query_count == 0) {
            return;
        }
        if (!std::is_sorted(queries, queries + query_count)) {
            throw std::invalid_argument("Curve: bulk queries must be sorted");
        }
        for (std::size_t i = 0; i < query_count; ++i) {
            output[i] = value_at(queries[i]);
        }
    }

    const std::vector<double>& pillars() const noexcept { return pillars_; }
    const std::vector<Value>& values() const noexcept { return values_; }

private:
    static void validate(
        const std::vector<double>& pillars, const std::vector<Value>& values)
    {
        if (pillars.empty()) {
            throw std::invalid_argument("Curve: at least one pillar is required");
        }
        if (pillars.size() != values.size()) {
            throw std::invalid_argument("Curve: pillar and value counts must match");
        }
        for (std::size_t i = 0; i < pillars.size(); ++i) {
            if (!std::isfinite(pillars[i])) {
                throw std::invalid_argument("Curve: pillars must be finite");
            }
            if (i != 0 && !(pillars[i] > pillars[i - 1])) {
                throw std::invalid_argument("Curve: pillars must be strictly increasing");
            }
            if (!detail::finite_value(values[i])) {
                throw std::invalid_argument("Curve: values must be finite");
            }
            if (i != 0 && !detail::same_shape(values[0], values[i])) {
                throw std::invalid_argument("Curve: value shapes must match");
            }
        }
    }

    Value extrapolate(double, bool left) const
    {
        const Extrapolation policy = left ? left_extrapolation_ : right_extrapolation_;
        if (policy == Extrapolation::Throw) {
            throw std::out_of_range(left
                ? "Curve: query precedes first pillar"
                : "Curve: query exceeds final pillar");
        }
        return left ? values_.front() : values_.back();
    }

    std::vector<double> pillars_;
    std::vector<Value> values_;
    Interpolation interpolation_;
    Extrapolation left_extrapolation_;
    Extrapolation right_extrapolation_;
};

} // namespace Curves
} // namespace DRC
