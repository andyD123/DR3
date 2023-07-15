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


#include "../Vectorisation/VecX/dr3.h"



//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
//using namespace DRC::VecF4F;
//using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
using namespace DRC::VecF8F;

//#include <iostream>

#include <iostream>

int main()
{
  //  int i = 32;
   // VecXX data(10., 64 *1024);
    for (int i = 1; i < 300000; i++)
    {
        VecXX data(1., i);
        data /= 9.;

        auto sumIt = [](auto x, auto y) {return x + y; };
        auto sum = ApplyAccumulate2UR_X_pairwise(data, sumIt);

        auto stdaccum = ApplyAccumulate2UR_X(data, sumIt);

       // std::cout << i << "\t" <<  std::setprecision(9)       << "standard sum "         << stdaccum  << ",   sum pairwise ="   << sum << "\n";
        std::cout << i << "\t" << std::setprecision(16) << "standard sum " << stdaccum << ",   sum pairwise =" << sum << "\n";
    }
}

