// scratch.cpp : This file contains the 'main' function. Program execution begins and ends there.
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



//using namespace DRC::VecDb;  // broken
//using namespace DRC::VecD2D; 
//using namespace DRC::VecF4F;
//using namespace DRC::VecD4D;
using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;

//using namespace DRC::VecLDb; //broken

#include <iostream>
#include <algorithm>





int main()
{


    long SZ = 10 *1024 * 1024 +3;
    
  
    VecXX data(10000.0/3.0, SZ);

    //bool USE_BIG_CANCELLATION =  false;
    bool USE_BIG_CANCELLATION = true;

    int i = 0;

    {

        auto mixed = data;
        long count = 0;
        for (auto& x : mixed)
        {

            count++;
            x += count;
            auto a =  100000000000000000.l;
            auto b = -100000000000000000.l;

            
            if (!USE_BIG_CANCELLATION)
            {
                a = 0.0;
                b = 0.0;
            }
 
            if ((count > 5) && (count % 5 == 0))
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

        auto trad_for_loop = 0.0l;
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
    
        // reduce with user defined accumulator
        using BINNED_ACCUMULATOR =  BinsT<VecXX::INS>;
        auto binned_Sum = reduce< BINNED_ACCUMULATOR >( mixed, BinnedAdd);

     
        std::cout << ++i << "\n" << std::setprecision(16) 
                 << trad_for_loop << "\t for loop sum  trad = \n"
                 << std_acc << "\t std::accumulate sum  trad = \n"
                 << DRCubedAccum << "\t accumulate DR3 = \n" 
                 << sumPairwiseDr3 << "\t  sum pairwise = \n"
                 << sumKahan << " ,\t Kahan acc ,\n"
                 << binned_Sum << "\t \t binned sum acc, \n \n \n \n";
           
    }
}
}

