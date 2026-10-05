#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace curve_demo {
struct Options { std::size_t bonds=2000, scenarios=65, repeats=3; };
inline Options options(int argc, char** argv) {
    Options out;
    for (int i=1; i<argc; ++i) {
        std::string name=argv[i];
        if (i+1==argc) throw std::invalid_argument("Option requires a value: "+name);
        const std::string text=argv[++i];
        if (text.empty() || text.find_first_not_of("0123456789")!=std::string::npos)
            throw std::invalid_argument("Options require positive integers");
        const auto n=std::stoull(text);
        if (name=="--bonds") out.bonds=n;
        else if(name=="--scenarios") out.scenarios=n;
        else if(name=="--repeats") out.repeats=n;
        else throw std::invalid_argument("Unknown option: "+name);
    }
    if (!out.bonds || out.bonds>100000 || !out.scenarios || out.scenarios>4096 ||
        !out.repeats || out.repeats>50 || out.bonds>10000000/out.scenarios)
        throw std::invalid_argument("Bounds: bonds 1..100000, scenarios 1..4096, repeats 1..50, pairs <= 10000000");
    return out;
}
template<class F> double timeMs(F f) {
    const auto start=std::chrono::steady_clock::now(); f();
    return std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
}
template<class F> double medianMs(std::size_t repeats, F f) {
    std::vector<double> times; times.reserve(repeats);
    for(std::size_t i=0;i<repeats;++i) times.push_back(timeMs(f));
    std::sort(times.begin(),times.end());
    const auto middle=times.size()/2;
    return times.size()%2 ? times[middle] : (times[middle-1]+times[middle])*.5;
}
inline void near(double actual, double expected, double& maxError) {
    if (!std::isfinite(actual) || !std::isfinite(expected)) throw std::runtime_error("Nonfinite price");
    const double error=std::abs(actual-expected);
    maxError=std::max(maxError,error);
    if(error>2e-11*(1+std::abs(expected))) throw std::runtime_error("Scalar-reference reconciliation failed");
}
// Canonical little-endian IEEE-754 VALUE bytes only. No file names or metadata.
inline std::uint64_t numericHash(const std::vector<double>& values) {
    static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,"IEEE binary64 required");
    std::uint64_t h=14695981039346656037ull;
    for(const double x: values) {
        std::uint64_t bits; std::memcpy(&bits,&x,sizeof bits);
        for(int i=0;i<8;++i) { h^=(bits>>(8*i))&255u; h*=1099511628211ull; }
    }
    return h;
}
inline void printHash(const std::vector<double>& values) {
    std::cout<<"computed_values_fnv64=0x"<<std::hex<<numericHash(values)<<std::dec<<'\n';
}
}
