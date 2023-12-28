#pragma once
#include "dr3.h"
#include "instruction_traits.h"
#include <type_traits>


template<typename INS_T>
struct BinsT
{
    using  INS = INS_T;
    //using FLOAT = InstructionTraits<INS_T>::FloatType;

    inline static constexpr bool isDbl = std::is_same<double, typename InstructionTraits<INS_T>::FloatType >::value;

    inline static const INS_T TINY_C{ isDbl ? 1e-32  : 1.0e-14f} ;
    inline static const INS_T VerySmall{isDbl ?  1e-16  : 1.0e-7f};
    inline static const INS_T SMALL{ isDbl ?  1.0  : 1.0f };
    inline static const INS_T BIG{ isDbl ? 1.0e16  : 1.0e7f };
  


    static inline auto roundIt(INS_T X, INS_T LEVEL)
    {
        auto INV_LEVEL = 1.0l / LEVEL;
        auto big = (LEVEL * round(X * INV_LEVEL));
        auto small = X - big;
        return std::pair(big, small);

    };



    INS_T veryBigSummV{ 0.0 };
    INS_T bigSummV{ 0.0 };
    INS_T smallSumV{ 0.0 };
    INS_T tinyV{ 0.0 };


    BinsT(){}


    BinsT(typename InstructionTraits<INS_T>::FloatType x)
    {
        alignas(512) std::array< InstructionTraits<INS_T>::FloatType, InstructionTraits<INS_T>::width> msk;

        //static alligned array for each type not quick fix below
        for(size_t i=1; i < InstructionTraits<INS_T>::width;++i)
        {
            msk[i] = 0;
        }
        msk[0] = 1;
        INS_T MASK;
        MASK.load_a(&msk[0]);
        set(MASK * x);

    }


    void set(INS_T x)
    {
        auto resRoundVeryBig = roundIt(x, BIG);
        auto resRoundBig = roundIt(resRoundVeryBig.second, SMALL);
        auto resRoundSmall = roundIt(resRoundBig.second, VerySmall);

        veryBigSummV = resRoundVeryBig.first;
        bigSummV = resRoundBig.first;
        smallSumV = resRoundSmall.first;
        tinyV = resRoundSmall.second;
    }

    BinsT(INS_T x)
    {
        set(x);
    }

    BinsT& operator *(INS_T rhs)
    {

        veryBigSummV *= rhs;
        bigSummV *= rhs;
        smallSumV *= rhs;
        tinyV *= rhs;

        return *this;
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



    auto hsum()
    {
        auto lambdaBinSum = [this]() {return (((horizontal_add(tinyV)) + horizontal_add(smallSumV)) + horizontal_add(bigSummV)) + horizontal_add(veryBigSummV); };
        return lambdaBinSum();
    }

   
};




static auto BinnedAdd = [](auto& bin, auto x) mutable
{
    bin += x;
    using  INS_T = decltype(x);
    auto NULL_Vec = INS_T(InstructionTraits<INS_T>::nullValue);
    return  NULL_Vec;

};