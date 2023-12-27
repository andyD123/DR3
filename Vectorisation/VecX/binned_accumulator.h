#pragma once
#include "dr3.h"
#include "instruction_traits.h"
#include <type_traits>


template<typename INS_T>
struct BinsT
{

    inline static constexpr bool isDbl = std::is_same<double, typename InstructionTraits<INS_T>::FloatType >::value;

    inline static const INS_T TINY_C{ isDbl ? 1e-32  : 1.0e-14} ;
    inline static const INS_T VerySmall{isDbl ?  1e-16  : 1.0e-7};
    inline static const INS_T SMALL{ isDbl ?  1.0  : 1.0 };
    inline static const INS_T BIG{ isDbl ? 1.0e16  : 1.0e7 };
  


    static inline auto roundIt(INS_T X, INS_T LEVEL)
    {
        auto INV_LEVEL = 1.0l / LEVEL;
        auto big = (LEVEL * round(X * INV_LEVEL));
        auto small = X - big;
        return std::pair(big, small);

    };



    INS_T veryBigSummV;// = 0.0l;
    INS_T bigSummV;// = 0.0l;
    INS_T smallSumV; //= 0.0l;
    INS_T tinyV;// = 0.0l;


    BinsT(){}


    BinsT(INS_T x)
    {

        auto resRoundVeryBig = roundIt(x, BIG);
        auto resRoundBig = roundIt(resRoundVeryBig.second, SMALL);
        auto resRoundSmall = roundIt(resRoundBig.second, VerySmall);

        veryBigSummV = resRoundVeryBig.first;
        bigSummV = resRoundBig.first;
        smallSumV = resRoundSmall.first;
        tinyV = resRoundSmall.second;

    }

    BinsT& operator *(INS_T rhs)
    {

        veryBigSummV *= rhs;
        bigSummV *= rhs;
        smallSumV *= rhs;
        tinyV *= rhs;

        return *this;
    }

    BinsT(double x)
    {

        auto resRoundVeryBig = roundIt(x, BIG);
        auto resRoundBig = roundIt(resRoundVeryBig.second, SMALL);
        auto resRoundSmall = roundIt(resRoundBig.second, VerySmall);

        veryBigSummV = resRoundVeryBig.first;
        bigSummV = resRoundBig.first;
        smallSumV = resRoundSmall.first;
        tinyV = resRoundSmall.second;

    }

    BinsT(BinsT&& x) noexcept
    {
        veryBigSummV = x.veryBigSummV;
        bigSummV = x.bigSummV;
        smallSumV = x.smallSumV;
        tinyV = x.tinyV;

    };


    BinsT& operator =(const BinsT& x) //noexcept
    {
        veryBigSummV = x.veryBigSummV;
        bigSummV = x.bigSummV;
        smallSumV = x.smallSumV;
        tinyV = x.tinyV;
        return *this;
    };


    BinsT& operator += (const BinsT& rhs)
    {
        auto resRoundTiny = roundIt(tinyV + rhs.tinyV, TINY_C);
        tinyV = resRoundTiny.second;
        auto smallRound = roundIt(smallSumV + rhs.smallSumV + resRoundTiny.first, SMALL);
        smallSumV = smallRound.second;
        auto bigRound = roundIt(smallRound.first + bigSummV + rhs.bigSummV, BIG);
        bigSummV = bigRound.second;
        veryBigSummV = bigRound.first + veryBigSummV + rhs.veryBigSummV;

        return *this;
    }



    double hsum()
    {
        auto lambdaBinSum = [this]() {return (((horizontal_add(tinyV)) + horizontal_add(smallSumV)) + horizontal_add(bigSummV)) + horizontal_add(veryBigSummV); };

        return lambdaBinSum();
    }

    using  INS =  INS_T;
};




static auto BinnedAdd = [](auto& bin, auto x) mutable
{

    using  INS_T =  decltype(x);

    auto roundIt = BinsT< INS_T>::roundIt;

  
    const auto  VerySmall = BinsT<INS_T>::VerySmall;
    const auto  SMALL = BinsT<INS_T>::SMALL;
    const auto  BIG = BinsT<INS_T>::BIG;
    const auto  TINY_C = BinsT<INS_T>::TINY_C;



    auto resRoundVeryBig = roundIt(x, BIG);
    auto resRoundBig = roundIt(resRoundVeryBig.second, SMALL);
    auto resRoundSmall = roundIt(resRoundBig.second, VerySmall);
    auto tinyRound = roundIt(resRoundSmall.second, TINY_C);
    bin.tinyV += (tinyRound.first + tinyRound.second);

    bin.smallSumV += (resRoundSmall.first + tinyRound.first);
    auto smallRound = roundIt(bin.smallSumV, SMALL);
    bin.smallSumV = smallRound.second;
    bin.bigSummV += resRoundBig.first + smallRound.first;
    auto bigRound = roundIt(bin.bigSummV, BIG);
    bin.bigSummV = bigRound.second;
    bin.veryBigSummV += bigRound.first;

    auto NULL_Vec = INS_T(InstructionTraits<INS_T>::nullValue);

    return  NULL_Vec;

};