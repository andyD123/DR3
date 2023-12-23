// GettingStarted.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <algorithm>
#include <random>
#include <numeric>
#include <iterator>
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>  
#include <functional>

#include <utility>

#include "../Vectorisation/VecX/dr3.h"

#include "../Vectorisation/intrinsic_utils.h"



//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
//using namespace DRC::VecF4F;
//using namespace DRC::VecD4D;
using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;

//using namespace DRC::VecLDb;

//#include <iostream>

#include <iostream>
#include <algorithm>
//#include <bits/stdc++.h>




template<typename INS_T>
struct BinsT
{

    auto roundIt(VecXX::INS X, VecXX::INS LEVEL)
        {
            auto INV_LEVEL = 1.0l / LEVEL;

         //   auto correct = INV_LEVEL * LEVEL;
           
            auto big = (LEVEL * round(X * INV_LEVEL));//
      //      big *= correct;
            auto small = X - big;

            return std::pair(big, small);

        };
  


    INS_T veryBigSummV;// = 0.0l;
    INS_T bigSummV;// = 0.0l;
    INS_T smallSumV; //= 0.0l;
    INS_T tinyV;// = 0.0l;


    const VecXX::INS VerySmall = 1e-16;
    const VecXX::INS SMALL = 1.0;
    const VecXX::INS BIG = 1.0e16;
    const VecXX::INS TINY_C = 1e-32;

    BinsT() {}


     BinsT(INS_T x)
    {


        const INS_T VerySmall = 1e-15;
        const INS_T SMALL = 1.0;
        const INS_T BIG = 1.0e15;
    //    const INS_T TINY_C = 1e-30; // ?

        auto roundIt = [&](auto X, auto LEVEL)
        {
            auto INV_LEVEL = 1.0l / LEVEL;
           
            auto big = (LEVEL * round(X * INV_LEVEL));
            auto small = X - big;

            return std::pair(big, small);

        };



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


        const VecXX::INS VerySmall = 1e-15;
        const VecXX::INS SMALL = 1.0;
        const VecXX::INS BIG = 1.0e15;
      //  const VecXX::INS TINY_C = 1e-30; // ?




        auto roundIt = [&](auto X, auto LEVEL)
        {
            auto INV_LEVEL = 1.0l / LEVEL;
            auto big = (LEVEL * round(X * INV_LEVEL));
            auto small = X - big;
            return std::pair(big, small);

        };



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
        auto resRoundTiny = roundIt(tinyV +rhs.tinyV, TINY_C);
        tinyV = resRoundTiny.second;
         
        auto smallRound = roundIt(smallSumV+ rhs.smallSumV+ resRoundTiny.first, SMALL);
        smallSumV = smallRound.second;
       
        auto bigRound = roundIt(smallRound.first +bigSummV + rhs.bigSummV, BIG);

        bigSummV = bigRound.second;

        veryBigSummV = bigRound.first + veryBigSummV + rhs.veryBigSummV;

        return *this;
    }



    double hsum()
    {
        auto lambdaBinSum = [this]() {return horizontal_add(tinyV) + horizontal_add(smallSumV) + horizontal_add(bigSummV) + horizontal_add(veryBigSummV); };

        return lambdaBinSum();
    }
};



int main()
{

    using FLOAT = InstructionTraits<VecXX::INS>::FloatType;
    //  int i = 32;
    long SZ = 10 *1024 * 1024 +2;// +3;// +7;// +7;// +1;// //  *100;// *100;// 100000000;///  512*512;// 63;// 900000;//512 * 512;// 1024 * 1024; // ///800;// 

   // long SZ = 100000;///  512*512;// 63;// 900000;//512 * 512;// 1024 * 1024; // ///800;// 

   //SZ = 8;

  //  FLOAT a = 10.;
    FLOAT b = 1.0 / 100.;
    FLOAT  c = 1.0 / (b * b);
 //   FLOAT d = c * c;

   // int n = 16;
   // int m = 1;
    /**/
    //__m256d Number= { a,b,c,d };
    //auto res = getExponent(Number);

    //auto res2 = getMantissa(Number);//, n, m);


    //Vec4d  number_t(Number);

    //Vec4d resVCL = getExponent(number_t);


    //Vec4d LVL = 32;
    //Vec4d zero = 0.;

    auto roundIt2 = [&](auto X, auto LEVEL)
    {
     
        auto integ = roundi(X / LEVEL);
         auto exp_level = getExponent(LEVEL);
        auto exp_X = getExponent(X);

        auto zero = X;
        zero = 0.0;

        auto trunc = exp_X - exp_level;
        auto bigger = max(zero, trunc);


        auto big = LEVEL * to_double(integ);
        auto small = X - big;

        return std::pair(big, small);
    };


  
    VecXX data(10000.0/3.0, SZ);

   // bool USE_BIG_CANCELLATION =  false;
    bool USE_BIG_CANCELLATION = true;

    int i = 0;

    {


        auto mixed = data;


        long count = 0;
        for (auto& x : mixed)
        {

            count++;
            auto a =  100000000000000000.l;
            auto b = -100000000000000000.l;

            
            if (!USE_BIG_CANCELLATION)
            {
                a = 0.0;
                b = 0.0;
            }
          
            if ((count > 15) && (count % 15 == 0))
            {
                auto c = mixed[count] + mixed[count - 1] + mixed[count - 2];
                mixed[count] = c;
                mixed[count - 1] = b;
                mixed[count - 2] = a;
            }
           

        }


        for(int kkk = 0; kkk < 100; kkk++)
        {

        std::random_device rd;
        std::mt19937 g(rd());

 
        std::shuffle(mixed.begin(), mixed.end(), g);

  
        auto std_acc =  std::accumulate(mixed.begin(), mixed.end(), 0.0l);


        auto sumIt = [](auto x, auto y) {return x + y; };
        auto sumPairwiseDr3 =  ApplyAccumulate2UR_X_pairwise(mixed, sumIt);

        
        auto DRCubedAccum =  ApplyAccumulate2UR_X(mixed, sumIt);

        auto trad = 0.0l;
        for (auto x : mixed)
        {
            trad += x;
        }



        auto NULL_Vec = VecXX::INS(0.0);
      
        auto KhanAddV = [c = NULL_Vec](auto sum, auto rhs) mutable
        {
            auto y = rhs - c;
            auto t = sum + y;
            c = (t - sum);
            c = c - y;
            sum = t;
            return t;
        };


        auto sumKhan =  reduce(mixed, KhanAddV);



        auto roundIt = [&](VecXX::INS X, VecXX::INS LEVEL)
            {
                auto INV_LEVEL = 1.0l / LEVEL;
                auto big = (LEVEL * round(X * INV_LEVEL));
                auto small = X - big;
                return std::pair(big, small);
            };


        auto BinnedAdd = [&](auto& bin, auto x) mutable
         {

 
            const VecXX::INS VerySmall = 1e-16;
            const VecXX::INS SMALL = 1.0;
            const VecXX::INS BIG = 1.0e16;
            const VecXX::INS TINY_C = 1e-32; // ?

  

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

            return  NULL_Vec;

        };

      



       
        BinsT< VecXX::INS> binsT;
      

        auto lambdaBinSum = reduceII(binsT, mixed, BinnedAdd);

     
        std::cout << ++i << "\n" << std::setprecision(16) << trad << "\t for loop sum  trad = \n"
                 << std_acc << "\t std::accumulate sum  trad = \n"
                 << DRCubedAccum << "\t accumulate DR3 = \n" 
                 << sumPairwiseDr3 << "\t  sum pairwise = \n"
                 << sumKhan << " ,\t Kahan acc ,\n"
                 << lambdaBinSum << "\t \t lambda bin sum acc, \n \n \n \n";
           
    }
}
}

