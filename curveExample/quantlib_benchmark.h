#pragma once
#include "quantlib_fixture.h"
#include "demo_support.h"
#include <ql/version.hpp>

namespace curve_demo {
namespace ql = QuantLib;
using namespace dr3_curve::quantlib;

struct Comparison {
    double snapshotMs=0, prepareMs=0, quantlibMs=0, scalarCachedMs=0, dr3Ms=0;
    double maxNpvError=0,maxDirtyError=0,maxCleanError=0;
    std::size_t dates=0;
    std::vector<double> quantlibValues, dr3Values, scalarValues;
};
inline std::size_t priceIndex(std::size_t bond,std::size_t scenario,std::size_t scenarios) {
    return (bond*scenarios+scenario)*3;
}
inline void reconcile(Comparison& report) {
    for(std::size_t i=0;i<report.dr3Values.size();i+=3) {
        near(report.dr3Values[i],report.quantlibValues[i],report.maxNpvError);
        near(report.dr3Values[i+1],report.quantlibValues[i+1],report.maxDirtyError);
        near(report.dr3Values[i+2],report.quantlibValues[i+2],report.maxCleanError);
        double ignored=0;
        for(std::size_t k=0;k<3;++k) near(report.scalarValues[i+k],report.quantlibValues[i+k],ignored);
    }
}
inline Comparison comparePricing(const std::vector<CurveHandle>& curves,
                                 const std::vector<BondPtr>& bonds,std::size_t repeats) {
    const auto scenarios=curves.size();
    if(curves.empty() || bonds.empty() || !repeats) throw std::invalid_argument("Empty benchmark");
    Comparison report;
    const auto size=3*bonds.size()*scenarios;
    report.quantlibValues.resize(size);report.dr3Values.resize(size);report.scalarValues.resize(size);

    Snapshot snapshot;
    report.snapshotMs=timeMs([&]{snapshot=snapshotCurves(curves,bonds);});
    report.dates=snapshot.discounts->dates();
    std::vector<AdaptedBond> prepared;
    report.prepareMs=timeMs([&]{
        prepared.reserve(bonds.size());
        for(const auto& bond:bonds) prepared.push_back(adaptBond(*bond,snapshot));
    });

    ql::RelinkableHandle<ql::YieldTermStructure> discount;
    const auto engine=ql::ext::make_shared<ql::DiscountingBondEngine>(discount,false);
    for(const auto& bond:bonds) bond->setPricingEngine(engine);
    Prices scratch(scenarios);
    auto normalQuantLib=[&]{
        for(std::size_t s=0;s<scenarios;++s) {
            discount.linkTo(curves[s].currentLink());
            for(std::size_t b=0;b<bonds.size();++b) {
                // Required even for a one-scenario repeat: never time LazyObject cache hits.
                bonds[b]->recalculate();
                const auto i=priceIndex(b,s,scenarios);
                report.quantlibValues[i]=bonds[b]->NPV();
                report.quantlibValues[i+1]=bonds[b]->dirtyPrice();
                report.quantlibValues[i+2]=bonds[b]->cleanPrice();
            }
        }
    };
    auto scalarCached=[&]{
        for(std::size_t b=0;b<prepared.size();++b) for(std::size_t s=0;s<scenarios;++s) {
            const auto i=priceIndex(b,s,scenarios);
            prepared[b].scalarScenario(s,report.scalarValues[i],report.scalarValues[i+1],report.scalarValues[i+2]);
        }
    };
    auto dr3Cached=[&]{
        for(std::size_t b=0;b<prepared.size();++b) {
            prepared[b].priceInto(scratch);
            for(std::size_t s=0;s<scenarios;++s) {
                const auto i=priceIndex(b,s,scenarios);
                report.dr3Values[i]=scratch.npv[s];
                report.dr3Values[i+1]=scratch.dirty[s];
                report.dr3Values[i+2]=scratch.clean[s];
            }
        }
    };
    normalQuantLib();scalarCached();dr3Cached();reconcile(report); // correctness before timing
    report.quantlibMs=medianMs(repeats,normalQuantLib);
    report.scalarCachedMs=medianMs(repeats,scalarCached);
    report.dr3Ms=medianMs(repeats,dr3Cached);
    reconcile(report); // observe/check the timed outputs too
    return report;
}
inline void printComparison(const Comparison& r,const Options& opt) {
    std::cout<<"QUANTLIB -> DR3: SAME CURVES, SAME BONDS, SCENARIO PRICING\n"
             <<"QuantLib="<<QL_VERSION<<" backend="<<dr3_curve::backend_name
             <<" bonds="<<opt.bonds<<" scenarios="<<opt.scenarios<<" unique_dates="<<r.dates<<" repeats="<<opt.repeats<<'\n'
             <<std::setprecision(10)
             <<"snapshot_ms="<<r.snapshotMs<<" adapt_pricers_ms="<<r.prepareMs<<'\n'
             <<"quantlib_engine_ms="<<r.quantlibMs<<" scalar_cached_ms="<<r.scalarCachedMs<<" dr3_cached_ms="<<r.dr3Ms<<'\n'
             <<"adapted_cold_components_ms="<<r.snapshotMs+r.prepareMs+r.dr3Ms
             <<" warm_speedup="<<r.quantlibMs/r.dr3Ms
             <<" cold_components_speedup="<<r.quantlibMs/(r.snapshotMs+r.prepareMs+r.dr3Ms)
             <<" scalar_cached_over_dr3="<<r.scalarCachedMs/r.dr3Ms<<'\n'
             <<"reconciliation=PASS max_abs_npv="<<r.maxNpvError
             <<" max_abs_dirty="<<r.maxDirtyError<<" max_abs_clean="<<r.maxCleanError<<'\n'
             <<"bond 0 (NPV currency; dirty/clean per 100 face):\n"
             <<"scenario     QuantLib NPV        DR3 NPV        dirty        clean\n";
    for(std::size_t s=0;s<std::min<std::size_t>(opt.scenarios,5);++s) {
        const auto i=priceIndex(0,s,opt.scenarios);
        std::cout<<s<<"            "<<r.quantlibValues[i]<<"     "<<r.dr3Values[i]<<"     "
                 <<r.dr3Values[i+1]<<"     "<<r.dr3Values[i+2]<<'\n';
    }
    printHash(r.dr3Values);
    std::cout<<"Timing: pricing includes writing all NPV/dirty/clean results; input construction is shared and excluded.\n"
             <<"Cold components sum one measured snapshot + preparation + median pricing; not a direct cold-run latency.\n"
             <<"Scalar-cached control separates prepared-data reuse from SIMD. No universal speedup claim.\n";
}
}
