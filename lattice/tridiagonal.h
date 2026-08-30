#pragma once

#include <cstddef>
#include <vector>

namespace dr3::lattice
{

struct ConstDoubleSpan
{
    ConstDoubleSpan(const double* values, std::size_t count) : data(values), size(count) {}
    explicit ConstDoubleSpan(const std::vector<double>& values)
        : data(values.data()), size(values.size()) {}

    const double& operator[](std::size_t index) const { return data[index]; }

    const double* data;
    std::size_t size;
};

struct DoubleSpan
{
    DoubleSpan(double* values, std::size_t count) : data(values), size(count) {}
    explicit DoubleSpan(std::vector<double>& values) : data(values.data()), size(values.size()) {}

    double& operator[](std::size_t index) const { return data[index]; }

    double* data;
    std::size_t size;
};

struct TridiagonalSystem
{
    std::vector<double> lower;
    std::vector<double> main;
    std::vector<double> upper;

    std::size_t size() const noexcept { return main.size(); }
    void validate() const;
};

class ThomasFactorization
{
public:
    explicit ThomasFactorization(const TridiagonalSystem& system,
                                 double pivotTolerance = 1.0e-14);

    std::size_t size() const noexcept { return pivots_.size(); }

    // No allocation occurs here: output and forwardWorkspace are caller-owned
    // and must both have exactly size() elements.
    void solve(ConstDoubleSpan rightHandSide,
               DoubleSpan output,
               DoubleSpan forwardWorkspace) const;

    void solve(const std::vector<double>& rightHandSide,
               std::vector<double>& output,
               std::vector<double>& forwardWorkspace) const
    {
        solve(ConstDoubleSpan(rightHandSide), DoubleSpan(output), DoubleSpan(forwardWorkspace));
    }

private:
    std::vector<double> lowerMultipliers_;
    std::vector<double> pivots_;
    std::vector<double> upper_;
};

} // namespace dr3::lattice
