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
/*
    inline static const INS_T TINY_C{ isDbl ? 1e-30 : 1.0e-14f };
    inline static const INS_T VERY_SMALL_C{ isDbl ? 1e-15 : 1.0e-7f };
    inline static const INS_T SMALL_C{ isDbl ? 1.0 : 1.0f };
    inline static const INS_T BIG_C{ isDbl ? 1.0e15 : 1.0e7f };
*/

    inline static const INS_T TINY_C{ isDbl ? pow(1024.0 , -10.0) : 1.0e-14f };
    inline static const INS_T VERY_SMALL_C{ isDbl ?   pow(1024.0,-5.0) : 1.0e-7f };
    inline static const INS_T SMALL_C{ isDbl ? 1.0 : 1.0f };
    inline static const INS_T BIG_C{ isDbl ? pow(1024.0, 5.0) : 1.0e7f };


    static inline auto roundIt(INS_T X, INS_T LEVEL)
    {
        INS_T ten = 10.0;
        auto INV_LEVEL = 1.0l / LEVEL;
        auto big = (LEVEL * truncate(X * INV_LEVEL));
        auto small = X - big;
       
        /*
        auto smallOrig = small;
        auto bigLog = iff(big >0.0, truncate(log10(big)),big);
        auto places = 16.0 - bigLog;
        auto multPow = pow(ten, places);

        small *= multPow;
        small = iff(big > 0.0, truncate(small) / multPow, smallOrig);
        */

       //auto truncatedSmall = truncate(small * INV_LEVEL*10.) * LEVEL;
       // small = select(abs(big) > 0.0, truncatedSmall, small);


        return std::pair(big, small);

    };


    INS_T m_scaleFactor{ InstructionTraits<INS_T>::oneValue };
    INS_T veryBigSummV{ InstructionTraits<INS_T>::nullValue };
    INS_T bigSummV{ InstructionTraits<INS_T>::nullValue };
    INS_T smallSumV{ InstructionTraits<INS_T>::nullValue };
    INS_T tinyV{ InstructionTraits<INS_T>::nullValue };


    INS_T TINY{ TINY_C };
    INS_T VERY_SMALL{ VERY_SMALL_C };
    INS_T SMALL{ SMALL_C };
    INS_T BIG{ BIG_C };





    BinsT() :
        m_scaleFactor{ InstructionTraits<INS_T>::oneValue },
        TINY{ m_scaleFactor * TINY_C },
        VERY_SMALL{ m_scaleFactor * VERY_SMALL_C },
        SMALL{ m_scaleFactor * SMALL_C },
        BIG{ m_scaleFactor * BIG_C }
    {}





    BinsT(typename InstructionTraits<INS_T>::FloatType x, typename InstructionTraits<INS_T>::FloatType scaleFactor = InstructionTraits<INS_T>::oneValue) :
        m_scaleFactor{ scaleFactor },
        TINY{ m_scaleFactor * TINY_C },
        VERY_SMALL{ m_scaleFactor * VERY_SMALL_C },
        SMALL{ m_scaleFactor * SMALL_C },
        BIG{ m_scaleFactor * BIG_C }
    {
    
        INS_T MASK(InstructionTraits<INS_T>::nullValue);
        MASK.insert(0, InstructionTraits<INS_T>::oneValue);

        set(MASK * x);

    }


    void set(INS_T x)
    {
        auto resRoundVeryBig = roundIt(x, BIG);
        auto resRoundBig = roundIt(resRoundVeryBig.second, SMALL);
        auto resRoundSmall = roundIt(resRoundBig.second, VERY_SMALL);

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

        m_scaleFactor = x.m_scaleFactor;
        TINY = x.TINY;
        VERY_SMALL = x.VERY_SMALL;
        SMALL = x.SMALL;
        BIG = x.BIG;


    };


    BinsT& operator =(const BinsT& x) 
    {
        veryBigSummV = x.veryBigSummV;
        bigSummV = x.bigSummV;
        smallSumV = x.smallSumV;
        tinyV = x.tinyV;

        m_scaleFactor = x.m_scaleFactor;
        TINY = x.TINY;
        VERY_SMALL = x.VERY_SMALL;
        SMALL = x.SMALL;
        BIG = x.BIG;


        return *this;
    };


    BinsT& operator += (const BinsT& rhs)
    {
        auto resRoundTiny = roundIt(tinyV + rhs.tinyV, VERY_SMALL);
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