#pragma once

#include "pricing_types.h"

namespace dr3::lattice::detail
{

void validate(const VanillaOption& option, TreeConfig config);
double intrinsic(const VanillaOption& option, double spot);
double priceBinomial(const VanillaOption& option, TreeConfig config);
double priceTrinomial(const VanillaOption& option,
                      TreeConfig config,
                      bool american,
                      double barrier,
                      double rebate,
                      bool barrierEnabled);

} // namespace dr3::lattice::detail
