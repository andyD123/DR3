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


#include "../Vectorisation/VecX/dr3.h"

#include "../Vectorisation/intrinsic_utils.h"



//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
//using namespace DRC::VecF4F;
using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;

//#include <iostream>

#include <iostream>

int main()
{
  //  int i = 32;
    long SZ = 1024 * 1024; // 1000000;//

    double a = 10.;
    double b = 1.0/100.;
    double  c =1.0/( b * b);
    double d = c * c;

/*
    __m256d Number= { a,b,c,d };
    auto res = getExponent(Number);

    Vec4d  number_t( Number);

    Vec4d resVCL = getExponent(number_t);

    */
    VecXX data(1., SZ);// 64 * 1024);

   

    double count = 1. *1/6.;

    long idx = 0;
   
    for (auto& x : data)
    {
          
        
        x = 1.0;// count;
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

        if (idx > 500000)
        {
            x *= 0.7e-11;
           // x = 0.0;
        }
        else
        {
            x = 10.11;
        }

        idx++;
    }

  

    int i = 0;
 
    {

       
        auto mixed = data;// 10 * data / 9. * 0.00000001;
        std::vector<double> dbg = mixed;

        auto c = mixed[252] + mixed[256];

        /*  */
     //   auto a = mixed[252] +  10000000000000.0;
     //   auto b = mixed[100256]  -10000000000000.0;
       
     a =   100000000000000000000.0;
     b =  -100000000000000000000.0;

    a = 0.;
    b = 0.;
    
     mixed[259] -=5055010.0;// 0;

  
  
       mixed[252] = a;
       mixed[256] = b;

       mixed[257] += c;// -5055010.0;

       mixed[258] += 1. / 60000000.;

   
        auto sumIt = [](auto x, auto y) {return x + y; };
        auto sum = ApplyAccumulate2UR_X_pairwise(mixed, sumIt);

        auto stdaccum = ApplyAccumulate2UR_X(mixed, sumIt);

        auto trad = 0.0;
        for (auto x : mixed)
        {
            trad += x;
        }

  
        const  double level = -10;

        const auto Mult = 100000000.0;
        const auto Div = 1. / Mult;


        auto sumItBig = [&](auto x, auto y)
        {
                  return select( log(abs(y)) > level,  x + round(y),x);
        };


        auto sumItSmall = [&](auto x, auto y)
        {
            /*
                if (log(abs(y)) <= level)
                {
                    y *= Mult;
                    y = y- round(y);
                   y*=Div;
                }
                */
           // auto x_dash = x;// -round(x);

                return select(log(abs(y)) <= level, x +( y-round(y)), x);
         //   };

            
            //return select(log(abs(y)) < level, (x - round(x)) + (y - round(y)), x );
        };


        auto sumItSmallCarry = [&](auto x, auto y)
        {
            /*
                if (log(abs(y)) <= level)
                {
                    y *= Mult;
                    y = y- round(y);
                   y*=Div;
                }
                */
            auto x_dash = x - round(x);
            auto y_dash = y - round(y);

           // auto x_dash = x - round(x);
            //auto unrounded_x_ret = x_dash + (y - round(y));

            auto sum = x_dash + y_dash;

            auto rt = sum - round(sum);


            return select(log(abs(y)) <= level, rt, x);
            //   };


               //return select(log(abs(y)) < level, (x - round(x)) + (y - round(y)), x );
        };



        auto sumItSmallRounded = [&](auto x, auto y)
        {
            auto x_dash = x - round(x);
            auto y_dash = y - round(y);

            return select(log(abs(y)) <= level, round(x_dash+y_dash), x);

        };



         /*     
        auto sumMiddle = [](auto x, auto y)
        {
            return  select(log(abs(y)) > 0.0, (x-round(x)) + (y-round(y)), x);
            //return select(log(abs(y)) < 0.0, x + y, x);
        };


        auto sumMiddleRnd = [](auto x, auto y)
        {
            return  select(log(abs(y)) > 0.0,  round(x) +  round(y), x);
            //return select(log(abs(y)) < 0.0, x + y, x);
        };
*/
        

        /*   
        auto sumItBig = [](auto x, auto y)
        {
            return select(getExponent(y) > 0.0, x + y, x);
        };
        auto sumItSmall = [](auto x, auto y)
        {
            return select(getExponent(y) <= 0.0, x + y, x);
        };
        */
   

        auto sumB =  reduceI(mixed, sumItBig);
        auto sumS =  reduceI(mixed, sumItSmall);
        auto sumM = 0.;
        //   reduceI(mixed, sumMiddle);
        auto sumMR = 0.;// reduceI(mixed, sumMiddleRnd);


        auto sumB_pairWise = ApplyAccumulate2UR_X_pairwise(mixed, sumItBig);
        auto sumS_pairWise_rnd = ApplyAccumulate2UR_X_pairwise(mixed, sumItSmallRounded);

        auto sumS_pairWise_cry = ApplyAccumulate2UR_X_pairwise(mixed, sumItSmallCarry);


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

    
        auto sumKhan = reduce(mixed, KhanAddV);

     


        SpanXX startspan(mixed.begin(), 32);



       double startSum = reduce(startspan, sumIt);


       auto roundIt = [&](auto X, auto LEVEL)
       {
           //auto rd = long( X / LEVEL);
           auto big = LEVEL * std::round(X / LEVEL);
           auto small = X - big;

           return std::pair( big , small );
       };

       
       double bigSumm = 0;
       double smallSum = 0;
       double verySmallSum = 0.0;
       
       double VerySmall = 1e-10;
       double SMALL = 1.0;
       double tiny = 0.0;
       double  BIG = 1.0e10;

       for (double x : mixed)
       {
           auto resRoundBig = roundIt(x, SMALL);
           auto resRoundSmall = roundIt(resRoundBig.second, VerySmall);

           tiny += resRoundSmall.second;
           auto tinyRound = roundIt(tiny, VerySmall);

           tiny = tinyRound.second;

           smallSum +=( resRoundSmall.first + tinyRound.first);

           auto smallRound = roundIt(smallSum, SMALL);

           smallSum = smallRound.second  ;

           bigSumm += resRoundBig.first + smallRound.first;
          



       }



      // const auto zero = InstructionTraits<VecXX::INS>::nullValue;
       struct Bins
       {


           VecXX::INS veryBigSummV = 0.0;
           VecXX::INS bigSummV = 0.0;
           VecXX::INS smallSumV = 0.0;
           VecXX::INS tinyV = 0.0;

          // struct LEVELS
           //{  };
            const VecXX::INS VerySmall = 1e-10;
            const VecXX::INS SMALL = 1.0;
            const VecXX::INS BIG = 1.0e10;
            //const  VecXX::INS tiny = 0.0; // ?
         



           Bins(double x)
           {
               auto roundIt = [&](auto X, auto LEVEL)
               {
                   //auto rd = long( X / LEVEL);
                   auto big = LEVEL * round(X / LEVEL);
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

           Bins( Bins&& x) noexcept
           {
               veryBigSummV = x.veryBigSummV;
               bigSummV = x.bigSummV;
               smallSumV = x.smallSumV;
               tinyV = x.tinyV;

           };
       };


       //doing accum in Bins
       Bins bin(0.);

      
       auto BinnedAdd = [ &](auto y,  auto x) mutable
       {
           auto roundIt = [&](auto X, auto LEVEL)
           {
               //auto rd = long( X / LEVEL);
               auto big = LEVEL * round(X / LEVEL);
               auto small = X - big;

               return std::pair(big, small);
           };

           auto resRoundVeryBig = roundIt(x, BIG);
           auto resRoundBig = roundIt(resRoundVeryBig.second, SMALL);
           auto resRoundSmall = roundIt(resRoundBig.second, VerySmall);
           bin.tinyV += resRoundSmall.second;
           auto tinyRound = roundIt(bin.tinyV, VerySmall);
           bin.tinyV = tinyRound.second;
           bin.smallSumV += (resRoundSmall.first + tinyRound.first);
           auto smallRound = roundIt(bin.smallSumV, SMALL);
           bin.smallSumV = smallRound.second;
           bin.bigSummV += resRoundBig.first + smallRound.first;
           auto bigRound = roundIt(bin.bigSummV, BIG);

           bin.bigSummV = bigRound.second;

           bin.veryBigSummV += bigRound.first;

          // return (bin.tinyV + bin.smallSumV) + bin.bigSummV;
           //reduce across these registers
           auto lambdaBinSum = ((horizontal_add(bin.tinyV) + horizontal_add(bin.smallSumV)) + horizontal_add(bin.bigSummV)) + horizontal_add(bin.veryBigSummV);

           return NULL_Vec;

       };

       // not working yet
       auto rrr=  reduceI(mixed, BinnedAdd);
       
       //need to do sorted add across registers or use binned add lambda
       auto lambdaBinSum =(horizontal_add(bin.tinyV) + horizontal_add(bin.smallSumV))+ horizontal_add(bin.bigSummV);

       // need to sort add for bbig





       // std::cout << i << "\t" <<  std::setprecision(9)       << "standard sum "         << stdaccum  << ",   sum pairwise ="   << sum << "\n";
        std::cout << i << "\t" << std::setprecision(16) << "standard sum  trad" << trad <<",accumulate " <<stdaccum << ",   sum pairwise =" << sum << ", split aggregate =" << bigSumm + smallSum + verySmallSum /* + sumS*/ << ", \n pairwise split aggregate =" << sumB_pairWise + sumS_pairWise_rnd + sumS_pairWise_cry << " , Kahan acc ," << sumKhan << " , lambda bin sum acc ," << lambdaBinSum << "\n";
    }
}

