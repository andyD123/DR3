#include "quantlib_benchmark.h"
#include <iostream>

namespace ql=QuantLib;
namespace adapter=dr3_curve::quantlib;
namespace fixture=curve_demo::ql_fixture;
int checks=0;
void check(bool condition,const char* what) { ++checks; if(!condition) throw std::runtime_error(what); }
template<class F> void rejects(F f) {bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,"Expected rejection");}

void reconcilePortfolio(const std::vector<adapter::CurveHandle>& curves,
                        const std::vector<adapter::BondPtr>& bonds) {
    const auto snapshot=adapter::snapshotCurves(curves,bonds);
    for(std::size_t s=0;s<curves.size();++s) {
        const auto engine=ql::ext::make_shared<ql::DiscountingBondEngine>(curves[s],false);
        for(const auto& bond:bonds) {
            bond->setPricingEngine(engine);
            bond->recalculate();
            const auto adapted=adapter::adaptBond(*bond,snapshot).prices();
            double error=0;
            curve_demo::near(adapted.npv[s],bond->NPV(),error); ++checks;
            curve_demo::near(adapted.dirty[s],bond->dirtyPrice(),error); ++checks;
            curve_demo::near(adapted.clean[s],bond->cleanPrice(),error); ++checks;
        }
    }
}

int main() {
    try {
        ql::SavedSettings saved;
        const ql::Date today(5,ql::October,2026);
        ql::Settings::instance().evaluationDate()=today;
        ql::Settings::instance().includeReferenceDateEvents()=false;
        ql::Settings::instance().includeTodaysCashFlows()=false;
        auto bonds=fixture::portfolio(today,9);
        for(std::size_t n:{1u,3u,9u,17u,65u})
            reconcilePortfolio(fixture::scenarioCurves(today,n),bonds);

        // Actual QuantLib convention boundaries, not hand-created cash amounts.
        std::vector<adapter::BondPtr> boundaries{
            fixture::boundaryBond(ql::Date(5,ql::October,2024),ql::Date(5,ql::October,2030)),
            fixture::boundaryBond(ql::Date(7,ql::October,2024),ql::Date(7,ql::October,2030)),
            fixture::boundaryBond(ql::Date(10,ql::October,2024),ql::Date(10,ql::October,2030),true),
            fixture::boundaryBond(ql::Date(30,ql::September,2024),ql::Date(30,ql::September,2026))};
        bool referenceFlow=false,settlementFlow=false,exCoupon=false;
        for(const auto& cf:boundaries[0]->cashflows()) referenceFlow|=cf->date()==today;
        for(const auto& cf:boundaries[1]->cashflows()) settlementFlow|=cf->date()==boundaries[1]->settlementDate();
        for(const auto& cf:boundaries[2]->cashflows()) exCoupon|=cf->tradingExCoupon(today);
        check(referenceFlow,"Fixture did not cover reference-date flow");
        check(settlementFlow,"Fixture did not cover settlement-date flow");
        check(exCoupon,"Fixture did not cover ex-coupon period");
        check(boundaries[3]->isExpired(),"Fixture did not cover expired bond");
        reconcilePortfolio(fixture::scenarioCurves(today,9),boundaries);

        // Frozen snapshots intentionally survive source quote updates. A fresh
        // snapshot must change and reconcile with fresh QuantLib prices.
        auto quote=ql::ext::make_shared<ql::SimpleQuote>(-.005);
        std::vector<adapter::CurveHandle> live{
            adapter::CurveHandle(ql::ext::make_shared<ql::FlatForward>(today,
                ql::Handle<ql::Quote>(quote),ql::Actual365Fixed()))};
        const auto oldSnapshot=adapter::snapshotCurves(live,bonds);
        const auto oldPricer=adapter::adaptBond(*bonds.front(),oldSnapshot);
        const double oldPrice=oldPricer.prices().npv[0];
        reconcilePortfolio(live,bonds);
        quote->setValue(.07);
        check(oldPricer.prices().npv[0]==oldPrice,"Frozen snapshot changed silently");
        const auto freshSnapshot=adapter::snapshotCurves(live,bonds);
        const double freshPrice=adapter::adaptBond(*bonds.front(),freshSnapshot).prices().npv[0];
        check(freshPrice<oldPrice,"Refreshed snapshot failed to change");
        reconcilePortfolio(live,bonds);

        rejects([&]{adapter::snapshotCurves({},bonds);});
        rejects([&]{adapter::snapshotCurves({adapter::CurveHandle()},bonds);});
        auto differentDate=live;
        differentDate.emplace_back(ql::ext::make_shared<ql::FlatForward>(today+1,.03,ql::Actual365Fixed()));
        rejects([&]{adapter::snapshotCurves(differentDate,bonds);});
        rejects([&]{adapter::adaptBond(*bonds.front(),adapter::Snapshot{});});
        rejects([&]{adapter::Prices wrong(2);oldPricer.priceInto(wrong);});
        ql::Settings::instance().evaluationDate()=today+1;
        rejects([&]{adapter::adaptBond(*bonds.front(),oldSnapshot);});
        check(oldPricer.prices().npv[0]==oldPrice,"Explicit frozen valuation changed with global date");
        std::cout<<"quantlib_adapter_tests PASS QuantLib="<<QL_VERSION
                 <<" backend="<<dr3_curve::backend_name<<" checks="<<checks<<'\n';
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"FAIL after "<<checks<<" checks: "<<e.what()<<'\n'; return 1;
    }
}
