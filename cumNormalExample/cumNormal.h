#pragma once






#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/apply_operation.h"
#include "../Vectorisation/VecX/vec_d.h"
#include "../Vectorisation/VecX/vec_bool.h"
#include "../Vectorisation/VecX/vec_view.h"

#include "../Vectorisation/VecX/target_name_space.h"

#include <immintrin.h>



//#include "norm.h"

//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
//using namespace DRC::VecD4D;
using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;


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



double getnull(double);




template<typename VEC>
VEC phi(VEC x)
{
	// https://stackoverflow.com/questions/2328258/cumulative-normal-distribution-function-in-c-c
	// references a Wests's implementation in Willmot.

 	const VEC z = abs(x);


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


	VEC n_c = (((((N6 * z + N5) * z + N4) * z + N3) * z + N2) * z + N1) * z + N0;
	VEC d_c = ((((((M7 * z + M6) * z + M5) * z + M4) * z + M3) * z + M2) * z + M1) * z + M0;
	VEC  central = n_c / d_c;

	constexpr double inv_RT2PI(0.39894228040143267793994605993438);
	VEC d_outer = (((((20. * z) * z + 13.) * z + 200.) * z + 78.) * z + 300.) * z + 39.;
	VEC n_outer = ((((20. * z) * z + 13.) * z + 180.) * z + 65.) * z + 160.;
	VEC outer = inv_RT2PI * n_outer / d_outer;

	
	VEC e = exp(-z * z * 0.5);

	 //   static const double SPLIT = 7.07106781186547; //orig
	const VEC SPLIT(7.42);// 7106781186547; //play  appears to give less error

	VEC RES = select((z < SPLIT), central, outer);
	RES *= e;

	return select(x <= VEC(0.0), RES, 1.0 - RES);

}




template <typename VecXX>
VecXX calcCDFNormal(const VecXX& X)
{
	   //TO DO FMA
	   auto centralLambda = [&](auto z)
	   {

			const static double N[] = { 3.52624965998911e-02 , 0.700383064443688,   6.37396220353165, 33.912866078383,  112.079291497871,  221.213596169931, 220.206867912376 };
			const static double M[] = { 8.83883476483184e-02, 1.75566716318264, 16.064177579207, 86.7807322029461 , 296.564248779674,  637.333633378831, 793.826512519948,440.413735824752 };

			auto n_c = (((((N[0] * z + N[1]) * z + N[2]) * z + N[3]) * z + N[4]) * z + N[5]) * z + N[6];
			auto d_c = ((((((M[0] * z + M[1]) * z + M[2]) * z + M[3]) * z + M[4]) * z + M[5]) * z + M[6]) * z + M[7];
			return n_c / d_c;
		};


	   auto outerLambda = [](auto z)
       {
			constexpr double inv_RT2PI(0.39894228040143267793994605993438);
			auto d_outer = (((((20. * z) * z + 13.) * z + 200.) * z + 78.) * z + 300.) * z + 39.;
			auto n_outer = ((((20. * z) * z + 13.) * z + 180.) * z + 65.) * z + 160.;
			return  inv_RT2PI * n_outer / d_outer;
	   };

	
	auto onePass = [=](auto x)
	{

		auto z = abs(x);

		auto e = exp(-z * z * 0.5);

		auto outer = outerLambda(z);

		auto central = centralLambda(z);

		auto  SPLIT =7.42;// 7106781186547; //play  appears to give less error

		auto RES = select((z < SPLIT), central, outer);
		RES *= e;

		return select(x <= 0.0, RES, 1.0 - RES);

	};



	auto res = ApplyTransformUR_XX(X, onePass);
	return res;

}



