#include "quantlib_benchmark.h"

int main(int argc,char** argv) {
    try {
        const auto opt=curve_demo::options(argc,argv);
        QuantLib::SavedSettings saved;
        const QuantLib::Date today(5,QuantLib::October,2026);
        QuantLib::Settings::instance().evaluationDate()=today;
        QuantLib::Settings::instance().includeReferenceDateEvents()=false;
        QuantLib::Settings::instance().includeTodaysCashFlows()=false;

        // These are ordinary QuantLib curves and FixedRateBond instruments.
        const auto curves=curve_demo::ql_fixture::scenarioCurves(today,opt.scenarios);
        const auto bonds=curve_demo::ql_fixture::portfolio(today,opt.bonds);

        // The complete migration interface:
        const auto snapshot=dr3_curve::quantlib::snapshotCurves(curves,bonds);
        const auto pricer=dr3_curve::quantlib::adaptBond(*bonds.front(),snapshot);
        const auto prices=pricer.prices();
        std::cout<<"first_bond_base_npv="<<std::setprecision(12)<<prices.npv[0]<<'\n';

        // The harness checks every bond/scenario against DiscountingBondEngine
        // and reports preparation separately from actual repricing.
        const auto result=curve_demo::comparePricing(curves,bonds,opt.repeats);
        curve_demo::printComparison(result,opt);
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"FAIL "<<e.what()<<'\n'; return 1;
    }
}
