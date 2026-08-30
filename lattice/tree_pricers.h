#pragma once

#include "pricing_types.h"

namespace dr3::lattice
{

double europeanBinomial(const VanillaOption& option, TreeConfig config);
double europeanTrinomial(const VanillaOption& option, TreeConfig config);
double americanTrinomial(const VanillaOption& option, TreeConfig config);

// Discretely monitored at every tree date. If the spot is at or above
// barrier, the option is knocked out and is worth rebate.
double upAndOutTrinomial(const VanillaOption& option,
                        TreeConfig config,
                        double barrier,
                        double rebate = 0.0);

} // namespace dr3::lattice
