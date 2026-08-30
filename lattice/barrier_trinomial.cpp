#include "tree_pricers.h"

#include "tree_pricer_detail.h"

#include <cmath>
#include <stdexcept>

namespace dr3::lattice
{

double upAndOutTrinomial(const VanillaOption& option,
                        TreeConfig config,
                        double barrier,
                        double rebate)
{
    detail::validate(option, config);
    if (!std::isfinite(barrier) || barrier <= 0.0)
    {
        throw std::invalid_argument("barrier must be finite and positive");
    }
    if (!std::isfinite(rebate) || rebate < 0.0)
    {
        throw std::invalid_argument("rebate must be finite and non-negative");
    }
    if (option.spot >= barrier)
    {
        return rebate;
    }
    if (option.maturity == 0.0)
    {
        return detail::intrinsic(option, option.spot);
    }
    return detail::priceTrinomial(option, config, false, barrier, rebate, true);
}

} // namespace dr3::lattice
