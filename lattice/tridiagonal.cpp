#include "tridiagonal.h"

#include <cmath>
#include <stdexcept>

namespace dr3::lattice
{

void TridiagonalSystem::validate() const
{
    if (main.empty() || lower.size() + 1 != main.size() || upper.size() + 1 != main.size())
    {
        throw std::invalid_argument("tridiagonal dimensions are inconsistent");
    }
    const auto finite = [](const std::vector<double>& values)
    {
        for (const double value : values)
        {
            if (!std::isfinite(value))
            {
                return false;
            }
        }
        return true;
    };
    if (!finite(lower) || !finite(main) || !finite(upper))
    {
        throw std::invalid_argument("tridiagonal coefficients must be finite");
    }
}

ThomasFactorization::ThomasFactorization(const TridiagonalSystem& system,
                                         double pivotTolerance)
    : lowerMultipliers_(system.lower.size()), pivots_(system.main), upper_(system.upper)
{
    system.validate();
    if (!std::isfinite(pivotTolerance) || pivotTolerance <= 0.0)
    {
        throw std::invalid_argument("pivot tolerance must be finite and positive");
    }
    const auto checkPivot = [pivotTolerance](double pivot)
    {
        if (!std::isfinite(pivot) || std::abs(pivot) <= pivotTolerance)
        {
            throw std::domain_error("tridiagonal pivot is too close to zero");
        }
    };

    checkPivot(pivots_[0]);
    for (std::size_t row = 1; row < pivots_.size(); ++row)
    {
        lowerMultipliers_[row - 1] = system.lower[row - 1] / pivots_[row - 1];
        pivots_[row] -= lowerMultipliers_[row - 1] * upper_[row - 1];
        checkPivot(pivots_[row]);
    }
}

void ThomasFactorization::solve(ConstDoubleSpan rightHandSide,
                                DoubleSpan output,
                                DoubleSpan forwardWorkspace) const
{
    if (rightHandSide.size != size() || output.size != size()
        || forwardWorkspace.size != size())
    {
        throw std::invalid_argument("tridiagonal solve buffers must exactly match the system size");
    }
    if ((rightHandSide.data == nullptr || output.data == nullptr
         || forwardWorkspace.data == nullptr) && size() != 0)
    {
        throw std::invalid_argument("tridiagonal solve buffers must not be null");
    }
    for (std::size_t row = 0; row < size(); ++row)
    {
        if (!std::isfinite(rightHandSide[row]))
        {
            throw std::invalid_argument("tridiagonal right-hand side must be finite");
        }
    }

    forwardWorkspace[0] = rightHandSide[0];
    for (std::size_t row = 1; row < size(); ++row)
    {
        forwardWorkspace[row] = rightHandSide[row]
            - lowerMultipliers_[row - 1] * forwardWorkspace[row - 1];
    }

    output[size() - 1] = forwardWorkspace[size() - 1] / pivots_[size() - 1];
    for (std::size_t row = size() - 1; row-- > 0;)
    {
        output[row] = (forwardWorkspace[row] - upper_[row] * output[row + 1]) / pivots_[row];
    }
}

} // namespace dr3::lattice
