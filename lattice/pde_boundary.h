#pragma once

#include <cmath>
#include <stdexcept>

namespace dr3::lattice
{

struct PdeBoundaryValues
{
    double lower;
    double upper;

    void validate() const
    {
        if (!std::isfinite(lower) || !std::isfinite(upper))
        {
            throw std::invalid_argument("PDE boundary values must be finite");
        }
    }
};

} // namespace dr3::lattice
