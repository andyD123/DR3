// scratch.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "../Vectorisation/VecX/dr3.h"


#include <random>
#include <numeric>
#include <iostream>
#include <vector>
#include <iomanip>
#include <functional>



//using namespace DRC::VecD2D; 
//using namespace DRC::VecF4F;
using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;



void doAddWithoutCancellation()
{

    using BINNED_ACCUMULATOR = BinsT<VecXX::INS>;
    using FLOAT = InstructionTraits<VecXX::INS>::FloatType;
    FLOAT oneThird = static_cast<FLOAT>(1.0 / 3.0);
    
    std::cout << "create empty bin, value 0.0 \n";
    BINNED_ACCUMULATOR bin;
    std::cout << "add one third to it , value =";

    bin += oneThird;
    std::cout << std::setprecision(8) << bin.hsum() << "\n";

    std::cout << "add one hundred thousand to it , value =";
    bin += 100000.0f;

    auto t = bin.hsum();

    std::cout << std::setprecision(8)<< t << std::endl;

    std::cout << "dd minus one hundred thousand to it , value =";
    bin += -100000.0f;
    t=bin.hsum();

    std::cout << std::setprecision(8) << t <<  std::endl;

    std::cout << "no cancellation! \n \n \n";


}

/*
 Example summation using std reduction, and for loop
 then pairwise_reduce  and  reduce with Kahan summation 
 and finally using binned summation.

 Generally for large sums, we get the same rounding error for  the for loop and
 std accumulation/ reduce

 However,  with pairwise reduce and  kahan accumulation, we tend to last digit level precision.

 Unfortunately, when we add pairs of large  +ve and -ve numbers which cancel each other
 they destroy accuracy of intermediate sums.  In the code we turn this on by
 setting  BIG_CANCELLATION  = true.

 The cancellation flag does not change the theoretical value of the sum,
 the sums should return the same result  as when CANCELLATION =false.
 However, we find that only the binned arithmetic scheme can achieve this
 sort of stability with this example.

 The input data set is randomly permuted  and the sum re calculated. The  ideal result is that we get the same
 answer irrespective of the ordering. The actual results differ to varying degrees.

 Both pairwise and kahan summation can significantly reduce rounding errors, however
 binned summation tends to work much better if we  have significant cancellation.
 For loops and std::accumulate /reduce are generally less accurate.

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
   //simple binned sum example
    doAddWithoutCancellation();

    long SZ = 10000 * 1024 ; // size of data set to be summed
    using FLOAT = InstructionTraits<VecXX::INS>::FloatType;
    FLOAT initVal = static_cast<FLOAT>(1.0 / 3.0);
    
    VecXX data(initVal, SZ);
    double scale = 1.0;// us power of 2   eg 1.0 / 1024.0 * 1.0 / 1024.0 * 1.0 / 1024.0;
    data *= scale;

    bool USE_BIG_CANCELLATION =  false;

    for( int C = 0; C < 2;C++) // iterate using cancellation data set
    {
        if (C >0) { USE_BIG_CANCELLATION =  true;}

        int i = 0;


        auto mixed = data;
        long count = 0;

        //make data members slightly different and add
        //cancelling values if required
        for (auto &x: mixed) {
            count++;
            x += count * 0.0001f;
            FLOAT a;
            setCancelInput(a);
            FLOAT b = -a;

            if (!USE_BIG_CANCELLATION) {
                a = 0.0;
                b = 0.0;
            }

            if ((count > 17) && (count % 17 == 0)) //every 17'th element set up cancellation
            {
                auto c = mixed[count] + mixed[count - 1] + mixed[count - 2];
                mixed[count] = c;
                mixed[count - 1] = b;
                mixed[count - 2] = a;
            }
            ignore(x);
        }

        //std::vector<FLOAT> scaledVec = mixed; //for debug observation
        // run ten permutations of  data set and do summation
        for (int kkk = 0; kkk < 10; kkk++)
        {
            std::random_device rd;
            std::mt19937 g(rd());

            std::shuffle(mixed.begin(), mixed.end(), g);
            // std::vector<FLOAT> obs= mixed;

            auto std_acc = std::accumulate(mixed.begin(), mixed.end(), static_cast<FLOAT>(0.0));
            auto std_reduce = std::reduce(mixed.begin(), mixed.end(), static_cast<FLOAT>(0.0));

            auto sumIt = [](auto x, auto y) { return x + y; };
            auto sumPairwiseDr3 = ApplyAccumulate2UR_X_pairwise(mixed, sumIt);

            auto DRCubedAccum = ApplyAccumulate2UR_X(mixed, sumIt);

            FLOAT trad_for_loop = 0.0f;
            for (auto x: mixed) {
                trad_for_loop += x;
            }

            //correcting summation  lambda
            auto NULL_Vec = VecXX::INS(0.0);
            auto KhanAddV = [c = NULL_Vec](auto sum, auto rhs) mutable {
                auto y = rhs - c;
                auto t = sum + y;
                c = (t - sum);
                c = c - y;
                sum = t;
                return t;
            };


            auto sumKahan = reduce(mixed, KhanAddV);

            NULL_Vec = VecXX::INS(0.0);
            auto sumPairwiseWithKahan = ApplyAccumulate2UR_X_pairwise(mixed, KhanAddV);

            // reduce with binned accumulator
            using BINNED_ACCUMULATOR = BinsT<VecXX::INS>;
            BINNED_ACCUMULATOR Bin(0.0);

            auto binned_Sum = reduceWithAccumulator(Bin, mixed, BinnedAdd);

            std::cout << "\nUsing Significant Cancellation Data = " << std::boolalpha << USE_BIG_CANCELLATION << "  \n";
            std::cout << "shuffled version " << ++i << "\n" << std::setprecision(16)
                      << trad_for_loop << "\t for loop sum   \n"
                      << std_acc << "\t std::accumulate sum  \n"
                      << std_reduce << "\t std::reduce  \n"
                      << DRCubedAccum << "\t accumulate DR3 \n"
                      << sumPairwiseDr3 << "\t sum pairwise  \n"
                      << sumKahan << "\t sum Kahan acc \n"
                      << sumPairwiseWithKahan << " \t pairwise_sum  using Kahan acc \n"
                      << binned_Sum << "\t binned sum acc \n \n \n \n";


        }
    }
}

