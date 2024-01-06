// scratch.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "../Vectorisation/VecX/dr3.h"

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


//using namespace DRC::VecD2D; 
//using namespace DRC::VecF4F;
using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;



void doAdd()
{

    

    using BINNED_ACCUMULATOR = BinsT<VecXX::INS>;
    using FLOAT = InstructionTraits<VecXX::INS>::FloatType;
    FLOAT oneThird = static_cast<FLOAT>(1.0 / 3.0);
    

    BINNED_ACCUMULATOR bin;

  

    bin += oneThird;
    std::cout << std::setprecision(8) << bin.hsum();

    bin += 100000.0f;

    auto t = bin.hsum();

    std::cout << std::setprecision(8)<< "\n" << t << "\n" << std::endl;

    bin += -100000.0f;
    t=bin.hsum();

    std::cout << std::setprecision(8) << "\n" << t << "\n" << std::endl;

    //////////////////////////
    FLOAT bigThird = oneThird * 1.e20f;

    BINNED_ACCUMULATOR bigbin;

 
    bigbin += bigThird;
    std::cout << std::setprecision(8) << bigbin.hsum();

    bigbin += 100000.0f  *1.e20f;

    t = bigbin.hsum();

    std::cout << std::setprecision(8) << "\n" << t << "\n" << std::endl;

    bigbin += -100000.0f* 1.e20f;
    t = bigbin.hsum();

    std::cout << std::setprecision(8) << "\n" << t << "\n" << std::endl;


}

/*
 example summation using std reduction, and for loop
 then pairwise_reduce  and  reduce with Kahan summation 
 and finally using binned summation.

 generally for large sums we get the same rounding error for  the for loop and
 std accumulation

 with pairwise _reduce and  kahan, we tend to last digit level precision.

 However when we add pairs of large  +ve and -ve numbers which cancel each other
 they destroy accuracy of intermediate sums.  This is achieved by running 
 using BIG_CANCELLATION  = true.

 the cancelation flag does not change the theoretical value of the sum.
 and the sums should return the same result  as when CANCELLATION =false.

 Binned arithmetic can achieve this sort of stability under the right circumstances.

 The input is randomly permuted  and re calculated. The  ideal result is that we get the same
 answer irrespective of the ordering.

 This can happen with pairwise and kahan summation if we have no significan cancellation.

*/

void setCancelInput(float& flt)
{
    flt = 100000.f;
}

void setCancelInput(double& dbl)
{
    dbl = 100000000000.0;
}


int main()
{
   
    doAdd();

    long SZ = 10000 * 1024 ;
  
    using FLOAT = InstructionTraits<VecXX::INS>::FloatType;
    FLOAT initVal = static_cast<FLOAT>(1.0 / 3.0);
    
    VecXX data(initVal, SZ);

    double scale = 1.0;// us power of 2   eg 1.0 / 1024.0 * 1.0 / 1024.0 * 1.0 / 1024.0;
    data *= scale;

 
     bool USE_BIG_CANCELLATION =  false;
    //bool USE_BIG_CANCELLATION = true;

    int i = 0;
    {

        auto mixed = data;
        long count = 0;
        for (auto& x : mixed)
        {

            count++;
            x += count * 0.0001f;
            FLOAT a;
            setCancelInput(a);
            FLOAT b = -a;

            if (!USE_BIG_CANCELLATION)
            {
                a = 0.0;
                b = 0.0;
            }
            
            if ((count > 17) && (count % 17 == 0))
            {
                auto c = mixed[count] + mixed[count - 1] + mixed[count - 2];
                mixed[count] = c;
                mixed[count - 1] = b;
                mixed[count - 2] = a;
            }
        
            ignore(x);
        }

        std::vector<FLOAT> scaledVec = mixed;

        for(int kkk = 0; kkk < 100; kkk++)
        {

        std::random_device rd;
        std::mt19937 g(rd());

 
        std::shuffle(mixed.begin(), mixed.end(), g);
        std::vector<FLOAT> obs= mixed;

  
        auto std_acc =  std::accumulate(mixed.begin(), mixed.end(), static_cast<FLOAT>(0.0));


        auto sumIt = [](auto x, auto y) {return x + y; };
        auto sumPairwiseDr3 =  ApplyAccumulate2UR_X_pairwise(mixed, sumIt);

        
        auto DRCubedAccum =  ApplyAccumulate2UR_X(mixed, sumIt);

        FLOAT trad_for_loop = 0.0f;
        for (auto x : mixed)
        {
            trad_for_loop += x;
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


        auto sumKahan =  reduce(mixed, KhanAddV);

        NULL_Vec = VecXX::INS(0.0);

        auto sumPairwiseWithKahan = ApplyAccumulate2UR_X_pairwise(mixed, KhanAddV);

    
        // reduce with binned accumulator
        using BINNED_ACCUMULATOR =  BinsT<VecXX::INS>;
        BINNED_ACCUMULATOR Bin(0.0);     

        auto binned_Sum = reduceWithAccumulator(Bin, mixed, BinnedAdd);

     
        std::cout << ++i << "\n" << std::setprecision(16) 
                 << trad_for_loop << "\t for loop sum  trad = \n"
                 << std_acc << "\t std::accumulate sum  trad = \n"
                 << DRCubedAccum << "\t accumulate DR3 = \n" 
                 << sumPairwiseDr3 << "\t  sum pairwise = \n"
                 << sumKahan << " ,\t Kahan acc ,\n"
                 << sumPairwiseWithKahan << " ,\t pairwise_sum  using Kahan acc ,\n"
                 << binned_Sum << "\t \t binned sum acc, \n \n \n \n";
           
    }
        
    }
}

