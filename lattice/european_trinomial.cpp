#include "tree_pricers.h"

#include "tree_pricer_detail.h"

namespace dr3::lattice
{

double europeanTrinomial(const VanillaOption& option, TreeConfig config)
{
    detail::validate(option, config);
    if (option.maturity == 0.0)
    {
        return detail::intrinsic(option, option.spot);
    }
    return detail::priceTrinomial(option, config, false, 0.0, 0.0, false);
}

} // namespace dr3::lattice
