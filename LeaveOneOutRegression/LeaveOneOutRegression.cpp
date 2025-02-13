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




#include <iostream>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include <tuple>

// Function to perform ridge regression on the provided data.
// It computes the intercept and slope by solving the regularized normal equations:
//   (X^T X + lambda * I) * w = X^T y,
// where X is the design matrix with a column of ones and x values.
void ridgeRegression(const std::vector<double>& x,
    const std::vector<double>& y,
    double lambda,
    double& slope,
    double& intercept)
{
    if (x.size() != y.size() || x.size() < 2) {
        throw std::invalid_argument("Vectors must have the same size and contain at least two elements.");
    }

    size_t n = x.size();
    double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumXX = 0.0;

    // Compute the necessary sums.
    for (size_t i = 0; i < n; ++i) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumXX += x[i] * x[i];
    }

    // For our model:
    //   y = intercept + slope * x,
    // the design matrix X is:
    //   [1  x1]
    //   [1  x2]
    //   ...
    //   [1  xn]
    //
    // Therefore:
    //   X^T X = [ [n,      sumX],
    //             [sumX,   sumXX] ]
    //
    // With ridge regularization (adding lambda to the diagonal):
    //   A = [ [n + lambda,    sumX      ],
    //         [sumX,          sumXX + lambda] ]
    //
    // and X^T y = [ sumY, sumXY ]^T.
    double A00 = n + lambda;
    double A01 = sumX;
    double A10 = sumX;
    double A11 = sumXX + lambda;

    double det = A00 * A11 - A01 * A10;
    if (det == 0) {
        throw std::runtime_error("The regularized matrix is singular; cannot solve the regression.");
    }

    // Compute the inverse of A:
    // A^{-1} = (1/det) * [ [ A11, -A01 ],
    //                      [ -A10, A00 ] ]
    double invA00 = A11 / det;
    double invA01 = -A01 / det;
    double invA10 = -A10 / det;
    double invA11 = A00 / det;

    // Multiply A^{-1} by X^T y to get w = [intercept, slope].
    intercept = invA00 * sumY + invA01 * sumXY;
    slope = invA10 * sumY + invA11 * sumXY;
}



// Function to perform leave-one-out cross-validation for ridge regression.
// For each observation, we remove it from the data, compute the ridge regression model
// on the remaining points, and then compute the squared error for the left-out point.
// Returns the sum of squared errors and also sets avgError to the average squared error.
auto leaveOneOutRidgeCV(const std::vector<double>& x,
    const std::vector<double>& y,
    double lambda,
    double& avgError)
{
    size_t n = x.size();
    if (n < 3) {
        throw std::invalid_argument("Leave-one-out regression requires at least three data points.");
    }

    double totalSquaredError = 0.0;
   


    std::vector<double> B_0;
    std::vector<double> B_1;

    // Loop over each observation as the left-out sample.
    for (size_t i = 0; i < n; ++i) {
        // Create training sets excluding the i-th data point.
        std::vector<double> xTrain;
        std::vector<double> yTrain;
        xTrain.reserve(n - 1);
        yTrain.reserve(n - 1);

        for (size_t j = 0; j < n; ++j) {
            if (j == i)
                continue;
            xTrain.push_back(x[j]);
            yTrain.push_back(y[j]);
        }

        double slope = 0.0, intercept = 0.0;
        ridgeRegression(xTrain, yTrain, lambda, slope, intercept);
        B_0.push_back(intercept);
        B_1.push_back(slope);


        // Predict the left-out observation.
        double prediction = intercept + slope * x[i];
        double error = y[i] - prediction;
        totalSquaredError += error * error;
    }

    avgError = totalSquaredError / n;
    //return totalSquaredError;
    auto ret = std::tuple{ B_0, B_1 };
    return   ret;
}




//switch the use of different instruction sets

//using namespace DRC::VecD2D;   // sse2 double
//using namespace DRC::VecF4F;   // sse2 float
//using namespace DRC::VecD4D;   // AVX2 double
using namespace DRC::VecD8D;     // AVX512 double 
//using namespace DRC::VecF16F;  // AVX512 float 
//using namespace DRC::VecF8F;   // AVX2    float



void fastLOO()
{


    int N = 1000000;

    VecXX data_X(N);
    VecXX data_Y(N);


    std::cout << "generating data set size " << N << std::endl;

    for (int i = 0; i < N; i++)
    {
        data_X[i] = i;
        data_Y[i] = 0.5 * i + 0.25 + rand() / double(RAND_MAX);
    }

    std::cout << "fitting data" << std::endl;

    double lambda = 0.0;
    const int LOOP_MAX =  200;

    //  std::vector<double>  debg;  

    const auto  startTme = std::chrono::high_resolution_clock::now();



    for (int LOOP = 0; LOOP < LOOP_MAX; LOOP++)
    {

        auto MULT = [](auto x, auto y) { return x * y; };
        auto SUM = [](auto x, auto y) { return x + y; };
        auto SQR = [](auto x) {return x * x; };

        //compute reductions for Sx,Sy,Sxx,Sxy
        auto S_x = reduce(data_X, SUM);
        auto S_y = reduce(data_Y, SUM);
        auto S_xx = transformReduce(data_X, SQR, SUM);
        auto S_xy = transformReduce(data_X, data_Y, MULT, SUM);

        // compute leave one out vectors

        auto SX_loo = S_x - data_X;   //leave one out SX
        auto SY_loo = S_y - data_Y;   //leave one out SY

        auto data_X_squared = data_X * data_X;
        auto SXX_loo = S_xx - data_X_squared; //leave one out SXX

        auto data_X_Y = data_X * data_Y;
        auto SXY_loo = S_xy - data_X_Y;  //leave one out SXY

        double lambda = 0.0;// 0.1;  //regularisation parameter
        double Sz = data_X.size() - 1.0;

        // Compute the fit parameters
        auto denominator = (Sz * (SXX_loo + lambda)) - (SX_loo * SX_loo);

        auto Beta_0_numerator = (SXX_loo + lambda) * SY_loo - SX_loo * SXY_loo;
        auto Beta_0 = Beta_0_numerator / denominator; //vector of fits for Beta 0 offsets

        auto Beta_1_numerator = Sz * SXY_loo - (SX_loo * SY_loo);
        auto Beta_1 = Beta_1_numerator / denominator;  //vector of Beta_1   slopes

    }

    
//    std::cout << "setting data" << std::endl;

//    std::vector<double> x = data_X;
 //   std::vector<double> y = data_Y;

   

      //   debg = denominator;
      //   debg = Beta_0_numerator;
      //   debg = Beta_1_numerator;

      //   std::vector<double> offset = Beta_0;
      //   std::vector<double> slope = Beta_1;
      //   std::vector<double> x = data_X;
      //   std::vector<double> y = data_Y;


    const auto endTme = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double, std::milli> fp_ms = endTme - startTme;

    std::cout << fp_ms.count() / LOOP_MAX << " milli seconds per fit";

}

void stdLOO()
{


    int N = 1000000; // takes about 4 hours 

    VecXX data_X(N);
    VecXX data_Y(N);


    std::cout << "generating data set size " << N << std::endl;

    for (int i = 0; i < N; i++)
    {
        data_X[i] = i;
        data_Y[i] = 0.5 * i + 0.25 + rand() / double(RAND_MAX);
    }


    double lambda = 0.0;
    const int LOOP_MAX = 1;// 200;

     std::cout << "setting data" << std::endl;

    std::vector<double> x = data_X;
    std::vector<double> y = data_Y;

    std::cout << "fitting data" << std::endl;

    const auto  startTme = std::chrono::high_resolution_clock::now();

    double sum_er = 0.0;
    auto X = leaveOneOutRidgeCV(x, y, lambda, sum_er);

    std::vector<double> B_0 = get<0>(X);

    std::vector<double> B_1 = get<1>(X);
 
    const auto endTme = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double, std::milli> fp_ms = endTme - startTme;

    std::cout << fp_ms.count() / LOOP_MAX << " milli seconds per fit";

}

int main()
{


    fastLOO();

  // stdLOO();

    return 0;


}

