// LeaveOneOutRegression.cpp : This file contains the 'main' function. Program execution begins and ends there.
//



#include "../Vectorisation/VecX/dr3.h"


#include <random>
#include <numeric>
#include <iostream>
#include <vector>
#include <iomanip>
#include <functional>
#include <chrono>



//using namespace DRC::VecD2D; 
//using namespace DRC::VecF4F;
//using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
using namespace DRC::VecF8F;




int main()
{


   // std::vector<double> vc(300, 1.1);

  //  VecXX myvec = vc;

    int N = 1000000;

    VecXX data_X(N);
    VecXX data_Y(N);

    for (int i = 0; i < N; i++)
    {
        data_X[i] = i;
        data_Y[i] = 0.5 * i + 0.25 +rand() / double(RAND_MAX);
    }

  //  std::vector<double>  debg;

  const auto  startTme = std::chrono::high_resolution_clock::now();

    auto MULT = [](auto x, auto y) { return x * y; };
    auto SUM = [](auto x, auto y) { return x + y; };
    auto SQR = [](auto x) {return x * x; };

    auto S_x = reduce(data_X,SUM);
    auto S_y = reduce(data_Y, SUM);

    auto S_xx = transformReduce(data_X, SQR, SUM);

    auto S_xy = transformReduce(data_X, data_Y, MULT, SUM);

    auto SX_loo = S_x - data_X;

 //   debg = SX_loo;

    auto SY_loo = S_y - data_Y;

//    debg = SY_loo;

    auto data_X_squared  = data_X * data_X;
    auto SXxx_loo = S_xx - data_X_squared;

 //   debg = SXxx_loo;

    auto data_X_Y = data_X * data_Y;
    auto SXY_loo = S_xy - data_X_Y;

  //  debg = SXY_loo;


    ///
    double lambda = 0.0;// 0.1;
    double Sz = data_X.size()-1.0;

    auto denominator = (Sz * (SXxx_loo + lambda)) - (SX_loo * SX_loo);

    auto Beta_0_numerator = (SXxx_loo + lambda) * SY_loo - SX_loo * SXY_loo;

    auto Beta_0 = Beta_0_numerator / denominator;

    auto Beta_1_numerator = Sz * SXY_loo - (SX_loo * SY_loo);

    auto Beta_1 = Beta_1_numerator / denominator;

 //   debg = denominator;
 //   debg = Beta_0_numerator;
 //  debg = Beta_1_numerator;

 //   std::vector<double> offset = Beta_0;
 //   std::vector<double> slope = Beta_1;
 //   std::vector<double> x = data_X;
 //   std::vector<double> y = data_Y;


   const auto endTme = std::chrono::high_resolution_clock::now();
   
  //  const std::chrono::duration<double, std::milli> runtime = endTme - startTme;
  
   // auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(runtime);

    const std::chrono::duration<double, std::milli> fp_ms = endTme - startTme;

    // integral duration: requires duration_cast
   // const auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    std::cout << fp_ms.count() <<  " milli seconds";
}

