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

    auto roundIt(VecXX::INS X, VecXX::INS LEVEL)//, auto INV_LEVEL)
        {
            auto INV_LEVEL = 1.0l / LEVEL;

            auto correct = INV_LEVEL * LEVEL;
            //auto rd = long( X / LEVEL);
            auto big = (LEVEL * round(X * INV_LEVEL));//
            big *= correct;
            auto small = X - big;

            return std::pair(big, small);

        };
   // VecXX::INS veryBigSummV = 0.0l;
   // VecXX::INS bigSummV = 0.0l;
   // VecXX::INS smallSumV = 0.0l;
   // VecXX::INS tinyV = 0.0l;


    INS_T veryBigSummV;// = 0.0l;
    INS_T bigSummV;// = 0.0l;
    INS_T smallSumV; //= 0.0l;
    INS_T tinyV;// = 0.0l;


  //  const VecXX::INS VerySmall = 1e-15;
  //  const VecXX::INS SMALL = 1.0;
  //  const VecXX::INS BIG = 1.0e15;
  //  const VecXX::INS TINY_C = 1e-30; // ?

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
        const INS_T TINY_C = 1e-30; // ?

        auto roundIt = [&](auto X, auto LEVEL)
        {
            auto INV_LEVEL = 1.0l / LEVEL;

            auto correct = INV_LEVEL * LEVEL;
            //auto rd = long( X / LEVEL);
            auto big = (LEVEL * round(X * INV_LEVEL)) * correct;
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
        const VecXX::INS TINY_C = 1e-30; // ?




        //const VecXX::INS VerySmall = 1e-8;
        //const VecXX::INS SMALL = 1.0;
        //const VecXX::INS BIG = 1.0e8;
        //const VecXX::INS TINY_C = 1e-15; // ?


        auto roundIt = [&](auto X, auto LEVEL)
        {
            auto INV_LEVEL = 1.0l / LEVEL;

            auto correct = INV_LEVEL * LEVEL;
            //auto rd = long( X / LEVEL);
            auto big = (LEVEL * round(X * INV_LEVEL)) * correct;
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
    long SZ = 1024 * 1024+5;// +3;// +7;// +7;// +1;// //  *100;// *100;// 100000000;///  512*512;// 63;// 900000;//512 * 512;// 1024 * 1024; // ///800;// 

   // long SZ = 100000;///  512*512;// 63;// 900000;//512 * 512;// 1024 * 1024; // ///800;// 

   //SZ = 8;

    FLOAT a = 10.;
    FLOAT b = 1.0 / 100.;
    FLOAT  c = 1.0 / (b * b);
    FLOAT d = c * c;

    int n = 16;
    int m = 1;
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
        //auto rd = long( X / LEVEL);
        auto integ = roundi(X / LEVEL);
        //auto rmdr = x % LEVEL; //roundi(X / LEVEL);
        //auto integ = (x - rmdr) / LEVEL;


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


    //auto rnd_res = roundIt2(number_t, LVL);



    VecXX data(0., SZ);// 64 * 1024);



    FLOAT count = 1.;// *1 / 3.;

    long idx = 0;

    for (auto& x : data)
    {


      //  x = FLOAT( 1.0);// count;
       /*
        if (idx % 4 == 0)
        {
            x += 100000000.111111;
        }
        else if ( idx % 2 == 0)
        {
           x -= 100000000.111111;
        }
        */

        if /* (idx > 500000)*/ (true) ///
        {
           //  x = 8e-17;
           //   x *= 8e2l;
             // x = 0.0;
            //x = count * 2.0 / 3.0;
            x = 10000* 1.0 / 3.0;//
          //  x += 100.*count * 1./100000000;;
            //  x += 1.l / 3.0l;// 0.1 * 2.0 / 3.0;
        }
        else
        {
            x = 10.11l;
        }

       // idx++;
        count++;
    }



    int i = 0;

    {


        auto mixed = data;// 10 * data / 9. * 0.00000001;

        /*   */    

      
        //   mixed[259] -=5055010.0;// 0;

/* */
        long count = 0;
        for (auto& x : mixed)
        {

            count++;
            auto a =  1000000000000.l;
            auto b = -1000000000000.l;

            

         //  a = 0.0;
         //  b = 0.0;

            if ((count > 15) && (count % 15 == 0))
            {
                auto c = mixed[count] + mixed[count - 1] + mixed[count - 2];
                mixed[count] = a;
                mixed[count - 1] = b;
                mixed[count - 2] = c;
            }


        }
  

        for(int kkk = 0; kkk < 100; kkk++)
        {

        std::random_device rd;
        std::mt19937 g(rd());

      // std::shuffle(v.begin(), v.end(), g);

       // mixed[800252] = -1;// a;
       // mixed[800253] = 1;// c;// b;

        shuffle(mixed.begin(), mixed.end(), g);

        /**/
        /* 
        auto c = mixed[800252] + mixed[800253] + mixed[800257];


        auto a = mixed[800252] + 10000000000000.0;
        auto b = mixed[800253] - 10000000000000.0;

          a =  100000000000000000000.0l;
          b = -100000000000000000000.0l;

        //a = 0.;
       //b = 0.;
          
        long count = 0;
        for (auto& x : mixed)
        {

            count++;
           auto a =  1000000000.l;
           auto b = -1000000000.l;

           if ( (count >97) && (count % 97 == 0) )
           {
               mixed[count] = a;
               mixed[count - 1] = b;
           }
           

        }
       
  //      mixed[800252] = a;// a;
  //      mixed[800253] =  b;// c;// b;
  //      mixed[800257] = c;// -5055010.0;





        shuffle(mixed.begin(), mixed.end(), g);

   */   
//   mixed[258] += 1. / 60000000.;

      //  std::vector<FLOAT> dbg = mixed;

        auto std_acc =  std::accumulate(mixed.begin(), mixed.end(), 0.0l);


        auto sumIt = [](auto x, auto y) {return x + y; };
   //     auto sumPairwiseDr3 =  ApplyAccumulate2UR_X_pairwise(mixed, sumIt);

        auto DRCubedAccum =  ApplyAccumulate2UR_X(mixed, sumIt);

        auto trad = 0.0l;
        for (auto x : mixed)
        {
            trad += x;
        }



        auto NULL_Vec = VecXX::INS(0.0);
        const auto zero = InstructionTraits<VecXX::INS>::nullValue;
        //lambdas
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




        SpanXX startspan(mixed.begin(), 32);



        double startSum = 0;// reduce(startspan, sumIt);


       // template<typename T> 
        auto roundIt = [&](VecXX::INS X, VecXX::INS LEVEL)//, auto INV_LEVEL)
            {
                auto INV_LEVEL = 1.0l / LEVEL;

                auto correct = INV_LEVEL * LEVEL;
                //auto rd = long( X / LEVEL);
                auto big = (LEVEL * round(X * INV_LEVEL));//
                big *= correct;
                auto small = X - big;

                return std::pair(big, small);

            };


        auto BinnedAdd = [&](auto& bin, auto x) mutable
        //auto BinnedAdd = [bin=Bins(0.0)](auto y, auto x) mutable
        {

 
            const VecXX::INS VerySmall = 1e-16;
            const VecXX::INS SMALL = 1.0;
            const VecXX::INS BIG = 1.0e16;
            const VecXX::INS TINY_C = 1e-32; // ?

        /*

            auto roundIt = [&](auto X, auto LEVEL)//, auto INV_LEVEL)
            {
                auto INV_LEVEL = 1.0l / LEVEL;

                auto correct = INV_LEVEL * LEVEL;
                //auto rd = long( X / LEVEL);
                auto big = (LEVEL * round(X * INV_LEVEL));//
                big *= correct;
                auto small = X - big;

                return std::pair(big, small);

            };

            */
            auto resRoundVeryBig = roundIt(x, BIG);
            auto resRoundBig = roundIt(resRoundVeryBig.second, SMALL);
            auto resRoundSmall = roundIt(resRoundBig.second, VerySmall);
            //bin.tinyV += resRoundSmall.second;
            auto tinyRound = roundIt(resRoundSmall.second, TINY_C);
            bin.tinyV += (tinyRound.first + tinyRound.second);

            // do truncation here 
            //auto resRoundTiny = roundIt(bin.tinyV, TINY_C);
            //bin.tinyV = resRoundTiny.first;

            bin.smallSumV += (resRoundSmall.first + tinyRound.first);
            auto smallRound = roundIt(bin.smallSumV, SMALL);
            bin.smallSumV = smallRound.second;
            bin.bigSummV += resRoundBig.first + smallRound.first;
            auto bigRound = roundIt(bin.bigSummV, BIG);

            bin.bigSummV = bigRound.second;

            bin.veryBigSummV += bigRound.first;

            // return (bin.tinyV + bin.smallSumV) + bin.bigSummV;
             //reduce across these registers
         //   auto lambdaBinSum = ((horizontal_add(bin.tinyV) + horizontal_add(bin.smallSumV)) + horizontal_add(bin.bigSummV)) + horizontal_add(bin.veryBigSummV);

            return  NULL_Vec;

        };

      


        VecXX a(1.0, 31);
        auto rK = reduceI(a, KhanAddV);


        // not working yet
       // auto rrr = reduceI(mixed, BinnedAdd);
       
        BinsT< VecXX::INS> binsT;
        auto nullVec = VecXX::INS(0.);

        auto lambdaBinSum = reduceII(binsT, mixed, BinnedAdd);// , nullVec, 0);



        // std::cout << i << "\t" <<  std::setprecision(9)       << "standard sum "         << stdaccum  << ",   sum pairwise ="   << sum << "\n";
        std::cout << i << "\n" << std::setprecision(16) << trad << "\t for loop sum  trad = \n"
            << std_acc << "\t std::accumulate sum  trad = \n"
            << DRCubedAccum << "\t accumulate DR3 = \n"  /* << ",\n split aggregate =" << bigSumm + smallSum + verySmallSum*/ /* + sumS*/ /* << ", \n pairwise split aggregate =" << sumB_pairWise + sumS_pairWise_rnd + sumS_pairWise_cry*/
         //   << sumPairwiseDr3 << "\t  sum pairwise = \n"
            << sumKhan << " ,\t Kahan acc ,\n"
            << lambdaBinSum << "\t lambda bin sum acc, \n \n \n \n";
    }
}
}

