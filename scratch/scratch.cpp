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


void doTrucate()
{

 VecXX::INS values = 1.0 / 7.0;

 auto expon = truncate(log2(values));

 auto lg2Vals = log2(values);

 auto mant = lg2Vals- expon  ;

 auto shft1 = exp2((mant+ expon));

 auto x = truncate(lg2Vals + 10) - 10;
 auto y = exp2(x);

 auto z = values - x;




 //trucate  1024 , 2^10  



 auto expon1 = floor(log2(values)) -10;

 auto  mant1 = lg2Vals + expon1;

 auto shft = exp2(mant1);





}



void doAdd()
{

    

    using BINNED_ACCUMULATOR = BinsT<VecXX::INS>;
    using FLOAT = InstructionTraits<VecXX::INS>::FloatType;
    FLOAT oneThird = static_cast<FLOAT>(1.0 / 3.0);
    

    BINNED_ACCUMULATOR bin;

    auto initVal = bin.hsum();

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

    initVal = bigbin.hsum();

    bigbin += bigThird;
    std::cout << std::setprecision(8) << bigbin.hsum();

    bigbin += 100000.0f  *1.e20f;

    t = bigbin.hsum();

    std::cout << std::setprecision(8) << "\n" << t << "\n" << std::endl;

    bigbin += -100000.0f* 1.e20f;
    t = bigbin.hsum();

    std::cout << std::setprecision(8) << "\n" << t << "\n" << std::endl;









}


int main()
{
    doTrucate();
    doAdd();

    long SZ = 10000 * 1024 -1;// *1024 + 5;
  //  SZ = 1000 * 1024 + 4;
    using FLOAT = InstructionTraits<VecXX::INS>::FloatType;
    FLOAT initVal = static_cast<FLOAT>(1.0 / 3.0);
    
    VecXX data(initVal, SZ);

    double scale = 1.0/1024.0 * 1.0 / 1024.0 * 1.0 / 1024.0;// 1 / (1024.0 * 1024.0 * 1024.);// 1.0e3;// 1.e16;// 1.0e-3;// 200;// 0.005;// 1.0e-18; //problem for 1.0e-3 -ve
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
           // x += count * 0.0000001f;// *0.0;
            FLOAT a = static_cast < FLOAT> (1000000000000.f);
            FLOAT b = static_cast<FLOAT>( - 1000000000000.f);

  /*  */
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
        

        }

        std::vector<double> scaledVec = mixed;

        for(int kkk = 0; kkk < 100; kkk++)
        {

        std::random_device rd;
        std::mt19937 g(rd());

 
        std::shuffle(mixed.begin(), mixed.end(), g);
        std::vector<FLOAT> obs= mixed;

  
        auto std_acc =  std::accumulate(mixed.begin(), mixed.end(), 0.0);


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
    
        // reduce with user defined accumulator
        using BINNED_ACCUMULATOR =  BinsT<VecXX::INS>;

        BINNED_ACCUMULATOR Bin(0.0);// , 8.);
     

        auto binned_Sum = reduceWithAccumulator(Bin, mixed, BinnedAdd);

        //auto binned_Sum = reduce< BINNED_ACCUMULATOR >( mixed, BinnedAdd,1.0e30);

     
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

