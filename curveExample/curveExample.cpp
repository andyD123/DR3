#include "dr3_pricing.h"
#include "demo_support.h"
int main(int argc,char** argv) {
 try {
    const auto opt=curve_demo::options(argc,argv);
    std::vector<double> dates; for(int i=1;i<=40;++i) dates.push_back(i*.5);
    const auto discount=[](double t,std::size_t s){return std::exp(-(.025+.0001*s+.0002*t)*t);};
    std::shared_ptr<const dr3_curve::DiscountSnapshot<double>> table;
    const double preparation=curve_demo::timeMs([&]{table=dr3_curve::DiscountSnapshot<double>::build(dates,opt.scenarios,discount);});
    std::vector<dr3_curve::PreparedLeg<double>> bonds;
    for(std::size_t b=0;b<opt.bonds;++b) {
        std::vector<dr3_curve::Payment> flows; const int periods=4+2*(b%19); const double coupon=1.+.05*(b%40);
        for(int i=1;i<=periods;++i) flows.push_back({table->index(.5*i),coupon+(i==periods?100.:0.)});
        bonds.emplace_back(table,std::move(flows));
    }
    std::vector<double> scalar(opt.bonds*opt.scenarios),vector(scalar.size());
    dr3_curve::Vector scratch(0.0,static_cast<int>(opt.scenarios));
    auto scalarRun=[&]{for(std::size_t b=0;b<opt.bonds;++b)for(std::size_t s=0;s<opt.scenarios;++s)scalar[b*opt.scenarios+s]=bonds[b].scalarScenario(s);};
    auto vectorRun=[&]{for(std::size_t b=0;b<opt.bonds;++b){bonds[b].priceInto(scratch);std::copy(scratch.begin(),scratch.end(),vector.begin()+b*opt.scenarios);}};
    scalarRun();vectorRun();double error=0;
    for(std::size_t i=0;i<scalar.size();++i)curve_demo::near(vector[i],scalar[i],error);
    const double scalarMs=curve_demo::medianMs(opt.repeats,scalarRun),vectorMs=curve_demo::medianMs(opt.repeats,vectorRun);
    std::cout<<"DR3 VECTOR CURVE / FIXED CASH-FLOW DEMONSTRATOR\nbackend="<<dr3_curve::backend_name<<" bonds="<<opt.bonds<<" scenarios="<<opt.scenarios<<" dates="<<table->dates()<<'\n'
             <<std::setprecision(10)<<"discount_snapshot_ms="<<preparation<<" scalar_cached_ms="<<scalarMs<<" dr3_cached_ms="<<vectorMs<<" scalar_cached_over_dr3="<<scalarMs/vectorMs<<'\n'
             <<"reconciliation=PASS max_abs_error="<<error<<" base_bond_price="<<vector[0]<<" last_scenario_price="<<vector[opt.scenarios-1]<<'\n';
    curve_demo::printHash(vector); return 0;
 }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
