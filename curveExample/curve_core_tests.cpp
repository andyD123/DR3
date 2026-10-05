#include "curve.h"
#include <functional>
#include <iostream>
#include <limits>
#include <string>

int checks = 0;
void check(bool ok, const char* what) { ++checks; if (!ok) throw std::runtime_error(what); }
void near(double a, double b) { check(std::isfinite(a) && std::isfinite(b) && std::abs(a-b) <= 2e-13*(1+std::abs(b)), "numerical mismatch"); }
template<class E, class F> void throws(F f) { bool caught=false; try { f(); } catch(const E&) { caught=true; } check(caught,"expected exception"); }
template<class C, class T, class V> void set(C& c, std::vector<T> x, std::vector<V> y) { c.setValues(x.begin(),x.end(),y.begin(),y.end()); }
int main() {
 try {
    Curve<long,double> c;
    throws<std::logic_error>([&]{c.valueAt(0);});
    set(c,std::vector<long>{0,10,20},std::vector<double>{1,3,7});
    near(c.valueAt(0),1); near(c.valueAt(5),2); near(c.valueAt(10),3); near(c.valueAt(20),7);
    near(c.valueAt(-1),1); near(c.valueAt(21),7);
    for (int i=0;i<=20;++i) near(c.valueAt(i),i<=10?1+.2*i:3+.4*(i-10));
    throws<std::invalid_argument>([&]{set(c,std::vector<long>{},std::vector<double>{});});
    throws<std::invalid_argument>([&]{set(c,std::vector<long>{1,2},std::vector<double>{1});});
    throws<std::invalid_argument>([&]{set(c,std::vector<long>{2,1},std::vector<double>{1,2});});
    throws<std::invalid_argument>([&]{set(c,std::vector<long>{1,1},std::vector<double>{1,2});});
    near(c.valueAt(5),2); // failed update is transactional
    set(c,std::vector<long>{5},std::vector<double>{2});
    near(c.valueAt(4),2); near(c.valueAt(5),2); near(c.valueAt(6),2);
    Curve<double,double,FlatInterp<double,double>> flat;
    set(flat,std::vector<double>{0,1,2},std::vector<double>{3,4,5});
    near(flat.valueAt(.5),3); near(flat.valueAt(1),4); near(flat.valueAt(2),5);
    const double nan=std::numeric_limits<double>::quiet_NaN(), inf=std::numeric_limits<double>::infinity();
    throws<std::invalid_argument>([&]{flat.valueAt(nan);});
    throws<std::invalid_argument>([&]{set(flat,std::vector<double>{0,inf},std::vector<double>{1,2});});
    throws<std::invalid_argument>([&]{set(flat,std::vector<double>{0,1},std::vector<double>{nan,2});});
    Curve2<double,double,ZeroInterpCached<double,double>> cached(2);
    Curve<double,double,ZeroInterp<double,double>> zero;
    for(double r: {-.02,0.,.04}) {
        set(zero,std::vector<double>{0,1,10},std::vector<double>{r,r+.001,r+.01});
        set(cached,std::vector<double>{0,1,10},std::vector<double>{r,r+.001,r+.01});
        for(double t:{0.,.5,1.,5.,10.,12.}) near(cached.valueAt(t),zero.valueAt(t));
        near(zero.valueAt(12),std::exp(-(r+.01)*12));
    }
    throws<std::out_of_range>([&]{zero.valueAt(-.1);});
    set(cached,std::vector<double>{0,10},std::vector<double>{.03,.03});
    auto* first=&cached.valueAtRef(2); check(first==&cached.valueAtRef(2),"cache hit copied");
    cached.valueAt(3); check(first==&cached.valueAtRef(2),"reference changed on unrelated insertion");
    auto copy=cached; near(copy.valueAt(2),std::exp(-.06));
    set(cached,std::vector<double>{0,10},std::vector<double>{.04,.04});
    check(cached.cachedDates()==0,"stale cache"); near(cached.valueAt(2),std::exp(-.08)); near(copy.valueAt(2),std::exp(-.06));
    cached.valueAt(11); set(cached,std::vector<double>{0,10},std::vector<double>{.05,.05}); near(cached.valueAt(11),std::exp(-.55));
    throws<std::invalid_argument>([]{lru_cache<int,double> bad(0);});
    lru_cache<int,double> lru(2); lru.put(1,10); lru.put(2,20); lru.get(1); lru.put(3,30);
    check(lru.exists(1)&&!lru.exists(2)&&lru.exists(3),"LRU eviction incorrect");
    auto lruCopy=lru; lru.clear(); near(lruCopy.get(1),10);
    lru=lruCopy; lruCopy.clear(); near(lru.get(3),30);
    auto moved=std::move(lru); near(moved.get(1),10); moved.put(1,11); near(moved.get(1),11);
    throws<std::out_of_range>([&]{moved.get(8);});
    std::cout<<"curve_core_tests PASS checks="<<checks<<'\n';
    return 0;
 } catch(const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<" after "<<checks<<" checks\n"; return 1; }
}
