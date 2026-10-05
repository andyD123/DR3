#include "dr3_pricing.h"
#include "demo_support.h"
#include <iostream>
#include <limits>
using dr3_curve::Vector;
int checks=0;
void check(bool x) {++checks; if(!x) throw std::runtime_error("Vector test assertion");}
template<class F> void throws(F f) {bool caught=false; try{f();}catch(const std::exception&){caught=true;}check(caught);}
int main() {
 try {
    double error=0;
    for(int n:{1,3,7,8,9,17,65,200}) {
        const std::vector<double> times{0,1,5};
        std::vector<Vector> rates;
        for(double t: times) {
            Vector r(0.0,n); for(int s=0;s<n;++s) r[s]=-.005+.0002*s+.001*t;
            rates.push_back(std::move(r));
        }
        Curve2<double,Vector,ZeroInterpCached<double,Vector>> c(3);
        c.setValues(times.begin(),times.end(),rates.begin(),rates.end());
        for(double t:{0.,.25,1.,2.5,5.,7.}) {
            const auto& actual=c.valueAtRef(t);
            for(int s=0;s<n;++s) {
                curve_demo::near(actual[s],std::exp(-(-.005+.0002*s+.001*std::min(t,5.))*t),error); ++checks;
            }
        }
        auto cCopy=c;
        for(auto& r:rates) for(int s=0;s<n;++s) r[s]=.03;
        c.setValues(times.begin(),times.end(),rates.begin(),rates.end());
        curve_demo::near(c.valueAt(2)[0],std::exp(-.06),error); ++checks;
        curve_demo::near(cCopy.valueAt(2)[0],std::exp(.006),error); ++checks;
        const auto table=dr3_curve::DiscountSnapshot<double>::build({0.,.5,1.,2.,2.},n,
            [](double t,std::size_t s){return std::exp(-(.01+.0003*s+.001*t)*t);});
        check(table->dates()==4); const auto* row=&table->row(2); check(row==&table->row(2));
        dr3_curve::PreparedLeg<double> bond(table,{{table->index(.5),2.},{table->index(1.),2.},{table->index(2.),102.}});
        auto pv=bond.prices();
        for(int s=0;s<n;++s) {
            const double r=.01+.0003*s;
            const double expected=2*std::exp(-(r+.0005)*.5)+2*std::exp(-(r+.001))+102*std::exp(-(r+.002)*2);
            curve_demo::near(pv[s],expected,error); ++checks;
        }
        // Scenario-dependent floating coupon: projection and discount curves differ.
        const auto fwd=dr3_curve::DiscountSnapshot<double>::build({.5,1.},n,
            [](double t,std::size_t s){return std::exp(-(.025+.0004*s)*t);});
        Vector floating=(fwd->row(0)/fwd->row(1)-1.0+.001*.5)*1000000.0*table->row(table->index(1.));
        for(int s=0;s<n;++s) {
            const double expected=(std::exp((.025+.0004*s)*.5)-1.+.0005)*1000000.*std::exp(-(.011+.0003*s));
            curve_demo::near(floating[s],expected,error); ++checks;
        }
        dr3_curve::PreparedLeg<double> empty(table,{}); auto zero=empty.prices(); for(int s=0;s<n;++s) check(zero[s]==0);
        throws([&]{table->index(3);});
        throws([&]{Vector wrong(0.0,n+1);bond.priceInto(wrong);});
        throws([&]{dr3_curve::PreparedLeg<double> bad(table,{{99,1}});});
        const std::vector<double> two{0,1};
        std::vector<Vector> mismatch{Vector(.02,n),Vector(.02,n+1)};
        throws([&]{c.setValues(two.begin(),two.end(),mismatch.begin(),mismatch.end());});
    }
    throws([]{dr3_curve::DiscountSnapshot<double>::build({1.},0,[](double,std::size_t){return 1.;});});
    throws([]{dr3_curve::DiscountSnapshot<double>::build({1.},3,[](double,std::size_t){return -1.;});});
    throws([]{dr3_curve::DiscountSnapshot<double>::build({1.},3,[](double,std::size_t){return std::numeric_limits<double>::quiet_NaN();});});
    std::cout<<"curve_vector_tests PASS backend="<<dr3_curve::backend_name<<" checks="<<checks<<" max_abs_error="<<std::setprecision(12)<<error<<'\n';
    return 0;
 } catch(const std::exception& e) {std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
