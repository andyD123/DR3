// LeaveOneOutRegression.cpp : This file contains the 'main' function. Program execution begins and ends there.
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
//using namespace DRC::VecD4D;
using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;




int main()
{


    std::vector<double> vc(300, 1.1);

    VecXX myvec = vc;

    int N = 1000000;

    VecXX data_X(N);
    VecXX data_Y(N);

    for (int i = 0; i < N; i++)
    {
        data_X[i] = i;
        data_Y[i] = 0.5 * i + 0.25;// +rand() / double(RAND_MAX);
    }



    auto MULT = [](auto x, auto y) { return x * y; };
    auto SUM = [](auto x, auto y) { return x + y; };
    auto SQR = [](auto x) {return x * x; };

    auto S_x = reduce(data_X,SUM);
    auto S_y = reduce(data_Y, SUM);

    auto S_xx = transformReduce(data_X, SQR, SUM);

    auto S_xy = transformReduce(data_X, data_Y, MULT, SUM);

    auto SX_loo = S_x - data_X;
    auto SY_loo = S_y - data_Y;

    auto data_X_squared  = data_X * data_X;
    auto SXxx_loo = S_xx - data_X_squared;

    auto data_X_Y = data_X * data_Y;
    auto SXY_loo = S_xy - data_X_Y;


    ///
    double lambda = 0.0;// 0.1;
    double Sz = data_X.size();

    auto denominator = Sz * (SXxx_loo + lambda) + SX_loo * SX_loo;

    auto Beta_0_numerator = (SXxx_loo + lambda) * SY_loo - SX_loo * SXY_loo;

    auto Beta_0 = Beta_0_numerator / denominator;

    auto Beta_1_numerator = Sz * SXY_loo - (SXxx_loo * SY_loo);

    auto Beta_1 = Beta_1_numerator / denominator;


    std::vector<double> offset = Beta_0;
    std::vector<double> slope = Beta_1;
    std::vector<double> x = data_X;
    std::vector<double> y = data_Y;





    std::cout << "Hello World!\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
