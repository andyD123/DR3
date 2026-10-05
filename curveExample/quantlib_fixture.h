#pragma once
// Shared deterministic example data, not part of the adapter API.
#include "quantlib_adapter.h"
#include <ql/pricingengines/bond/discountingbondengine.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/termstructures/yield/zerospreadedtermstructure.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/daycounters/actualactual.hpp>
#include <ql/time/schedule.hpp>

namespace curve_demo::ql_fixture {
namespace ql = QuantLib;
using namespace dr3_curve::quantlib;
inline std::vector<CurveHandle> scenarioCurves(ql::Date reference, std::size_t count) {
    ql::TARGET calendar;
    std::vector<ql::Date> pillars{reference};
    for(int year:{1,3,7,15,31}) pillars.push_back(calendar.advance(reference,ql::Period(year,ql::Years)));
    const std::vector<ql::Rate> baseRates{.02,.022,.025,.027,.029,.03};
    const CurveHandle base(ql::ext::make_shared<ql::ZeroCurve>(pillars,baseRates,ql::Actual365Fixed(),calendar));
    std::vector<CurveHandle> curves;
    for(std::size_t s=0;s<count;++s) {
        const double shift=s==0?0.:(s%2?1.:-1.)*.0001*(1.+s/2);
        if(s%2) {
            const ql::Handle<ql::Quote> quote(ql::ext::make_shared<ql::SimpleQuote>(shift));
            curves.emplace_back(ql::ext::make_shared<ql::ZeroSpreadedTermStructure>(base,quote));
        } else {
            auto rates=baseRates;
            for(std::size_t j=0;j<rates.size();++j) rates[j]+=shift*(.5+.1*j);
            auto dates=pillars;
            // Different interior grids are intentional: this adapter does not
            // assume a shared grid or replace QuantLib's interpolation.
            if(s) dates[2]+=static_cast<ql::Integer>(s%11);
            curves.emplace_back(ql::ext::make_shared<ql::ZeroCurve>(dates,rates,ql::Actual365Fixed(),calendar));
        }
    }
    return curves;
}
inline std::vector<BondPtr> portfolio(ql::Date reference,std::size_t count) {
    ql::TARGET calendar;
    std::vector<BondPtr> bonds; bonds.reserve(count);
    for(std::size_t i=0;i<count;++i) {
        const auto issue=calendar.advance(reference,ql::Period(-static_cast<int>(1+i%15),ql::Months));
        const auto maturity=calendar.advance(reference,ql::Period(static_cast<int>(2+i%19),ql::Years));
        const ql::Schedule schedule(issue,maturity,ql::Period(ql::Semiannual),calendar,
            ql::ModifiedFollowing,ql::ModifiedFollowing,ql::DateGeneration::Backward,false);
        bonds.push_back(ql::ext::make_shared<ql::FixedRateBond>(2,1000.+5.*(i%11),schedule,
            std::vector<ql::Rate>{.02+.0007*(i%40)},ql::ActualActual(ql::ActualActual::ISMA),
            ql::ModifiedFollowing,100.,issue));
    }
    return bonds;
}
inline BondPtr boundaryBond(ql::Date issue,ql::Date maturity,bool exCoupon=false) {
    ql::TARGET calendar;
    const ql::Schedule schedule(issue,maturity,ql::Period(ql::Annual),calendar,
        ql::ModifiedFollowing,ql::ModifiedFollowing,ql::DateGeneration::Backward,false);
    return ql::ext::make_shared<ql::FixedRateBond>(2,1000.,schedule,std::vector<ql::Rate>{.04},
        ql::ActualActual(ql::ActualActual::ISMA),ql::ModifiedFollowing,100.,issue,
        calendar,exCoupon?ql::Period(7,ql::Days):ql::Period(),calendar,ql::Unadjusted,false);
}
}
