#include "cumNormal.h"

#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/apply_operation.h"
#include "../Vectorisation/VecX/vec_d.h"
#include "../Vectorisation/VecX/vec_bool.h"
#include "../Vectorisation/VecX/vec_view.h"

#include "../Vectorisation/VecX/target_name_space.h"


#include <algorithm>
#include <random>
#include <numeric>
#include <iterator>
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>  
#include <chrono>
#include <iostream>
#include <functional>


#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <format>


#include <immintrin.h>




#include "../Vectorisation/VecX/vcl_latest.h"




//#include "norm.h"

//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
//using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;





auto numOps = [](int TEST_LOOP_SZ, int SZ) { return  static_cast<int>(double(TEST_LOOP_SZ) * double(SZ)); };



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
        //auto res =calcCDFNormsSparseFMAOnePass(X)[0]; //OK
        //auto res = calcCDFNormsSparseFMA(X)[0]; // good
        //auto res = calcCDFNormWichuraViewsAndFMA(X)[0]; // superb accuract
        //auto res = calcCDFNormWithViewsAndFMA(X)[0]; // superb accuract
        //auto res = calcCDFNormWichuraViewsAndFMA(X)[0]; // superb accuract 
        auto phi_res = phi(x);
        auto res = calcCDFNormal(X)[0]; // superb accuract




        std::cout<< std::setprecision (16)  << x << "," << res << "expected" << normalCDF(x) <<  "phi impl ="<< phi_res << ",err" << res - normalCDF(x) << std::endl;
    }

}



//set instruction set via namespace in cdfNormalInverse.h
int main()
{
    //doPerformanceComparison();
    	compareNumerics();
    //profileExample();
    return 0;
}


void doPerformanceComparison()
{

    long loop = 10000;
    //for (int SZ = 100; SZ < 10000; SZ += 100)
    //{

    using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
    //	std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));


        //for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }

    std::random_device rd;
    std::mt19937 g(rd());
    //	std::shuffle(&v[0], &v[SZ - 1], g);


        //VecXX test(v);
        //auto test2 = test;

        //VecXX resultsWichura;
        //std::vector<double> resT(test.size(), .0);

        //
        ///*
    for (int SZ = 100; SZ < 10000; SZ += 100)
    {
        std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
        for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }
        std::shuffle(&v[0], &v[SZ - 1], g);


        VecXX test(v);
        auto test2 = test;

        //warm
        for (long l = 0; l < loop; l++)
        {
            auto 	res = calcCDFNormal(test);
        }

        auto startTme = std::chrono::high_resolution_clock::now();
        for (long l = 0; l < loop; l++)
        {
            auto 	res = calcCDFNormal(test);
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto runtime = endTime - startTme;
        std::cout << SZ << "," << " calcCDFNormsSparseFMA ," << runtime.count() / 1000000000.0;//<< std::endl;
        std::cout << " , evluations run , " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << ",  evals per sec |, \n";// << std::endl;
    }

}



///////////////////////////////////



double phi(double x)
{
    // https://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
    // references a Wests's implementation in Willmot.

    static  double RT2PI = std::sqrt(4.0 * acos(0.0));

    static  double inv_RT2PI = 1.0 / RT2PI;

 //   static const double SPLIT = 7.07106781186547; //orig
    static const double SPLIT = 7.42;// 7106781186547; //play  appears to give less error

    static const double N0 = 220.206867912376;
    static const double N1 = 221.213596169931;
    static const double N2 = 112.079291497871;
    static const double N3 = 33.912866078383;
    static const double N4 = 6.37396220353165;
    static const double N5 = 0.700383064443688;
    static const double N6 = 3.52624965998911e-02;
    static const double M0 = 440.413735824752;
    static const double M1 = 793.826512519948;
    static const double M2 = 637.333633378831;
    static const double M3 = 296.564248779674;
    static const double M4 = 86.7807322029461;
    static const double M5 = 16.064177579207;
    static const double M6 = 1.75566716318264;
    static const double M7 = 8.83883476483184e-02;

    const double z = fabs(x);
    double c = 0.0;

    if (z <= 37.0)
    {
        const double e = exp(-z * z / 2.0);
        if (z < SPLIT)
        {
            const double n = (((((N6 * z + N5) * z + N4) * z + N3) * z + N2) * z + N1) * z + N0;
            const double d = ((((((M7 * z + M6) * z + M5) * z + M4) * z + M3) * z + M2) * z + M1) * z + M0;
          c = e * n / d;
        }
        else
        {
          //  const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
          // this expand to rational polynomial 

            /*
            
            39 +300Z + 78Z^2 +  200Z^3  +  13 Z^4   + 20 Z^5
            ------------------------------------------------
            160 + 65Z  + 180 Z^2   + 13 Z^3  + 20 z^4
            
            */

            const double n = (((((20. * z) * z + 13.) * z + 200.) * z + 78.) * z + 300.) * z + 39.;

            const double d = ((((20. * z) * z + 13.) * z + 180.)* z + 65.)* z + 160.;

           // const double f = n / d;
            const double inv_f = d / n;

           // c = e / (RT2PI * f);

            c = e * inv_f * inv_RT2PI;
           
        }
    }
    return x <= 0.0 ? c : 1.0 - c;
}




 /*



template<typename VEC>
 VEC phi(VEC x)
{
    // https://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
    // references a Wests's implementation in Willmot.

   //static  const double rtpi = std::sqrt(4.0 * std::acos(0.0));
   // = 2.506628274631000502415765284811

   // const VEC RT2PI(rtpi);// = std::sqrt(4.0 * acos(0.0));

  //  const VEC inv_RT2PI = 1.0 / RT2PI;
    // 0.39894228040143267793994605993438

   

   

    const VEC z = abs(x);

    //auto centralLambda = [&](VEC z)
    //{
        constexpr double N0 = 220.206867912376 ;
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


        VEC n_c = (((((N6 * z + N5) * z + N4) * z + N3) * z + N2) * z + N1) * z + N0;
        VEC d_c = ((((((M7 * z + M6) * z + M5) * z + M4) * z + M3) * z + M2) * z + M1) * z + M0;
        VEC  central = n_c / d_c;
 //       return c;

   // };

  //  auto outerLambda = [](VEC z)
  //  {
        constexpr double inv_RT2PI(0.39894228040143267793994605993438);
        VEC d_outer = (((((20. * z) * z + 13.) * z + 200.) * z + 78.) * z + 300.) * z + 39.;
        VEC n_outer = ((((20. * z) * z + 13.) * z + 180.) * z + 65.) * z + 160.;        
        VEC outer = inv_RT2PI * n_outer / d_outer ;
  //      return c;

    //};

   
   
   //const VEC half(1.0 / 2.0);

  //  const VEC one(1.0 );
    VEC e = exp(-z * z * 0.5);

   // VEC central = centralLambda(z);
   // VEC outer = outerLambda(z);

    //   static const double SPLIT = 7.07106781186547; //orig
    const VEC SPLIT(7.42);// 7106781186547; //play  appears to give less error

    VEC RES = select( (z < SPLIT), central, outer);
    RES *= e; 

    return select(x <= VEC(0.0), RES, 1.0- RES);


}



 int main(int c, char** argv)
 {


     return 0;
 }





int main(int c, char** argv)
{

    Vec4d test = { 1.0, 2.0,3.0,3.142 };

    auto  norm = phi(test);



    std::cout << std::setprecision(16);
    for (double x = -6.; x < 9.0; x += 0.001)
    {

        Vec4d tst = x;
        double vres = phi(tst)[0];

        std::cout << " blah \n" << x << "\n" << phi(x) << "\n" << normalCDF(x) << "\n"  << vres <<"\n  err"  << phi(x) - normalCDF(x) << "\n \n";
    }



	return 0;
}

*/