#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace dr3::lattice
{

struct GridLocation
{
    std::size_t lowerIndex;
    std::size_t upperIndex;
    double upperWeight;
};

class Grid1D
{
public:
    Grid1D(double minimum, double maximum, std::size_t nodeCount)
        : coordinates_(nodeCount), spacing_((maximum - minimum) / (nodeCount - 1))
    {
        if (nodeCount < 3)
        {
            throw std::invalid_argument("a one-dimensional grid requires at least three nodes");
        }
        if (!std::isfinite(minimum) || !std::isfinite(maximum) || maximum <= minimum)
        {
            throw std::invalid_argument("grid bounds must be finite and increasing");
        }
        if (!std::isfinite(spacing_) || spacing_ <= 0.0)
        {
            throw std::invalid_argument("grid spacing must be finite and positive");
        }
        for (std::size_t index = 0; index < nodeCount; ++index)
        {
            coordinates_[index] = minimum + static_cast<double>(index) * spacing_;
        }
        coordinates_.back() = maximum;
    }

    std::size_t nodeCount() const noexcept { return coordinates_.size(); }
    double spacing() const noexcept { return spacing_; }
    double minimum() const noexcept { return coordinates_.front(); }
    double maximum() const noexcept { return coordinates_.back(); }
    std::size_t lowerBoundary() const noexcept { return 0; }
    std::size_t upperBoundary() const noexcept { return nodeCount() - 1; }
    std::size_t interiorBegin() const noexcept { return 1; }
    std::size_t interiorEnd() const noexcept { return nodeCount() - 1; }
    const std::vector<double>& coordinates() const noexcept { return coordinates_; }
    double operator[](std::size_t index) const { return coordinates_.at(index); }

    GridLocation locate(double coordinate) const
    {
        if (!std::isfinite(coordinate) || coordinate < minimum() || coordinate > maximum())
        {
            throw std::invalid_argument("coordinate must be finite and inside the grid");
        }
        const auto upper = std::lower_bound(coordinates_.begin(), coordinates_.end(), coordinate);
        if (upper == coordinates_.end())
        {
            const auto last = nodeCount() - 1;
            return {last, last, 0.0};
        }
        const auto upperIndex = static_cast<std::size_t>(upper - coordinates_.begin());
        if (*upper == coordinate || upperIndex == 0)
        {
            return {upperIndex, upperIndex, 0.0};
        }
        const auto lowerIndex = upperIndex - 1;
        const double weight = (coordinate - coordinates_[lowerIndex])
            / (coordinates_[upperIndex] - coordinates_[lowerIndex]);
        return {lowerIndex, upperIndex, weight};
    }

    double interpolate(const std::vector<double>& values, double coordinate) const
    {
        if (values.size() != nodeCount())
        {
            throw std::invalid_argument("grid values must match the logical node count");
        }
        const auto location = locate(coordinate);
        if (location.lowerIndex == location.upperIndex)
        {
            return values[location.lowerIndex];
        }
        return values[location.lowerIndex] * (1.0 - location.upperWeight)
            + values[location.upperIndex] * location.upperWeight;
    }

private:
    std::vector<double> coordinates_;
    double spacing_;
};

} // namespace dr3::lattice
