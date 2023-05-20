#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/apply_operation.h"
#include "../Vectorisation/VecX/vec_d.h"
#include "../Vectorisation/VecX/vec_bool.h"
#include "../Vectorisation/VecX/vec_view.h"
#include "../Vectorisation/VecX/target_name_space.h"


#include "cumNormal.h"

#include <algorithm>
#include <random>
#include <numeric>
#include <iterator>
#include <iostream>
#include <vector>
#include <chrono>
#include <functional>
#include <iomanip>
#include <cmath>
#include <immintrin.h>




#include "../Vectorisation/VecX/vcl_latest.h"




auto numOps = [](int TEST_LOOP_SZ, int SZ) { return  static_cast<int>(double(TEST_LOOP_SZ) * double(SZ)); };

double phi(double x);


void compareNumerics();
void doPerformanceComparison();
void profileExample();



template<typename T>
bool vectorsEqual(const std::vector<T>& C1, const std::vector<T>& C2, const std::vector<T>& C3)
{
    bool  testOK = true;
    const double ERR = 1e-15; //for examples
    if (C1.size() != C2.size())
    {
        return false;
    }

    if (C3.size() != C2.size())
    {
        return false;
    }

    for (size_t i = 0; i < C3.size(); i++)
    {
        auto err1 = fabs((C1[i] - C2[i]) / (C2[i] + C1[i]));
        auto err2 = fabs((C1[i] - C3[i]) / (C1[i] + C3[i]));

        if ((err1 > ERR) || (err2 > ERR))
        {
            testOK = false;
            std::cout << "\n err diff@ " << i << " err1 =" << err1 << ", err2 = " << err2 << "\n";
            break;
        }
    }

    return testOK;

}


double normalCDF(double value)
{
    return 0.5 * erfc(-value * std::sqrt(0.5));
}


double phi(double x);


//compare different implementations
void compareNumerics()
{

    for (double x = 0.001; x < 1.0; x += 0.003)
    {
        VecXX X(x, 16);
        auto phi_res = phi(x);
        auto res = calcCDFNormalFMA(X)[0]; 

        std::cout<< std::setprecision (16)  << x << "," << res << "expected" << normalCDF(x) <<  "phi impl ="<< phi_res << ",err" << res - normalCDF(x) << std::endl;
    }

}



//set instruction set via namespace in cdfNormalInverse.h
int main()
{
   doPerformanceComparison();
 //   compareNumerics();
 //   profileExample();
    return 0;
}


void doPerformanceComparison()
{

    long loop = 10000;
  
    using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
    std::random_device rd;
    std::mt19937 g(rd());
   

    for (int SZ = 100; SZ < 10000; SZ += 100)
    {
        std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
        std::vector<FloatType>  result(SZ, VecXX::SCALA_TYPE(6.66));

        for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ)-0.5; }
        std::shuffle(&v[0], &v[SZ - 1], g);


        VecXX test(v);
        auto test2 = test;
        {
            //warm up
            for (long l = 0; l < loop; l++)
            {
                auto 	res = calcCDFNormalFMA(test);
            }

            auto startTme = std::chrono::high_resolution_clock::now();
            for (long l = 0; l < loop; l++)
            {
                auto 	res = calcCDFNormalFMA(test);
            }
            auto endTime = std::chrono::high_resolution_clock::now();
            auto runtime = endTime - startTme;
            std::cout << SZ << "," << " calcCDFNormalFMA ," << runtime.count() / 1000000000.0;
            std::cout << " , evluations run , " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << ",  evals per sec |, ";
        }

        // timing using STD LIB Function
        {
            //warmup 
            for (long l = 0; l < loop; l++)
            {
                std::transform(begin(v), end(v), begin(result), [](FloatType x) { return normalCDF(x); });
            }

            auto startTme = std::chrono::high_resolution_clock::now();
            for (long l = 0; l < loop; l++)
            {
                  std::transform(begin(v), end(v), begin(result), [](FloatType x) { return normalCDF(x); });
              
            }
            auto endTime = std::chrono::high_resolution_clock::now();
            auto runtime = endTime - startTme;
            std::cout << SZ << "," << " standard lib erfc approach ," << runtime.count() / 1000000000.0;
            std::cout << " , evluations run , " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << ",  evals per sec |, \n";
        }


    }



}



void profileExample()
{
   
    long loop = 1000000;// 100;
    int SZ = 2000;
    {

        using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
        std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));


        for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ) -0.5; }

        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(&v[0], &v[SZ - 1], g);

        VecXX test(v);
        auto test2 = test;


        {
            auto startTme = std::chrono::high_resolution_clock::now();
            for (long l = 0; l < loop; l++)
            {
                auto 	res = calcCDFNormalFMA(test);      
            }

            auto endTime = std::chrono::high_resolution_clock::now();
            auto runtime = endTime - startTme;
            std::cout << " , evluations run , " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << ",  evals per sec |,";// << std::endl;

        }
    }
}






double phi(double x)
{
    // https://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
    // references a Wests's implementation in Willmot.

    constexpr double inv_RT2PI = 0.39894228040143267794;

 //   static const double SPLIT = 7.07106781186547; //orig
    constexpr double SPLIT = 7.42;// 7106781186547; //play  appears to give less error

    constexpr double N0 = 220.206867912376;
    constexpr double N1 = 221.213596169931;
    constexpr double N2 = 112.079291497871;
    constexpr double N3 = 33.912866078383;
    constexpr double N4 = 6.37396220353165;
    constexpr double N5 = 0.700383064443688;
    constexpr double N6 = 3.52624965998911e-02;
    constexpr double M0 = 440.413735824752;
    constexpr double M1 = 793.826512519948;
    constexpr double M2 = 637.333633378831;
    constexpr double M3 = 296.564248779674;
    constexpr double M4 = 86.7807322029461;
    constexpr double M5 = 16.064177579207;
    constexpr double M6 = 1.75566716318264;
    constexpr double M7 = 8.83883476483184e-02;

    double z = fabs(x);
    double c = 0.0;

    if (z <= 37.0)
    {
         double e = exp(-z * z / 2.0);
        if (z < SPLIT)
        {
             double n = (((((N6 * z + N5) * z + N4) * z + N3) * z + N2) * z + N1) * z + N0;
             double d = ((((((M7 * z + M6) * z + M5) * z + M4) * z + M3) * z + M2) * z + M1) * z + M0;
          c = e * n / d;
        }
        else
        {
          //  const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
          //  this expands to rational polynomial 

            /*
            
            39 +300Z + 78Z^2 +  200Z^3  +  13 Z^4   + 20 Z^5
            ------------------------------------------------
            160 + 65Z  + 180 Z^2   + 13 Z^3  + 20 z^4
            
            */
            //we could extract first common term from both polynomials

            double n = (((((20. * z) * z + 13.) * z + 200.) * z + 78.) * z + 300.) * z + 39.;
            double d = ((((20. * z) * z + 13.) * z + 180.)* z + 65.)* z + 160.;

           // const double f = n / d;
           double inv_f = d / n;
            c = e * inv_f * inv_RT2PI;
           
        }
    }
    return x <= 0.0 ? c : 1.0 - c;
}


