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
using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
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



long double qnorm8(long double p);



template <typename VecXX>
VecXX calcCDFNormWichuraViewsAndFMA(VecXX& X)
{
	//wichura
	const static  double c[] = { 7.7454501427834140764e-4 , .0227238449892691845833 ,.24178072517745061177, 1.27045825245236838258 ,  3.64784832476320460504, 5.7694972214606914055, 4.6303378461565452959, 1.42343711074968357734 };
	const static  double d[] = { 1.05075007164441684324e-9 , 5.475938084995344946e-4, .0151986665636164571966, .14810397642748007459, .68976733498510000455,  1.6763848301838038494,  2.05319162663775882187,1. };
	const  static  double e[] = { 2.01033439929228813265e-7 ,   2.71155556874348757815e-5,   .0012426609473880784386, .026532189526576123093, .29656057182850489123,   1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772 };
	const static  double f[] = { 2.04426310338993978564e-15 , 1.4215117583164458887e-7, 1.8463183175100546818e-5,  7.868691311456132591e-4, .0148753612908506148525,.13692988092273580531, .59983220655588793769, 1. };



	/////////FMA Lambda ///////////////////////////
	auto lambdaMain = [=](const auto& p)
	{
		auto q = p - 0.5;
		auto r = .180625 - q * q;
		auto denom = 1. / (mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(5226.495278852854561, r, 28729.085735721942674), r, 39307.89580009271061), r, 21213.794301586595867), r, 5394.1960214247511077), r, 687.1870074920579083), r, 42.313330701600911252), r, 1));
		auto num = (mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(2509.0809287301226727, r, 33430.575583588128105), r, 67265.770927008700853), r, 45921.953931549871457), r, 13731.693765509461125), r, 1971.5909503065514427), r, 133.14166789178437745), r, 3.387132872796366608) * q);
		return denom * num;
	};




	auto isOuterRangelambda = [=](auto p)
	{
		const auto p_low = 0.5 - 0.425;
		const auto p_high = 1.0 - p_low;
		return (p_high < p) || (p < p_low);
	};


	auto dolambdaLow = [=](auto r, auto p)
	{
		auto q = p - 0.5;
		r += -1.6;
		auto	val = mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(c[0], r, c[1]), r, c[2]), r, c[3]), r, c[4]), r, c[5]), r, c[6]), r, c[7]);
		auto denom = (decltype(p))(1.0) / mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(d[0], r, d[1]), r, d[2]), r, d[3]), r, d[4]), r, d[5]), r, d[6]), r, (decltype(p))(1.0));
		val = val * denom;
		auto valMult = iff(q < (decltype(p))(0.0), (decltype(p))(-1.0), (decltype(p))(1.0));
		val *= valMult;
		return val;
	};




	auto dolambdaHi = [=](auto r, auto p)
	{

		auto	val = (((((((r * e[0] + e[1]) * r + e[2]) * r + e[3]) * r + e[4]) * r + e[5]) * r + e[6]) * r + e[7])
			/ (((((((r * f[0] + f[1]) * r + f[2]) * r + f[3]) * r + f[4]) * r + f[5]) * r + f[6]) * r + 1.);

		auto valMult = iff((p - 0.5) < (decltype(p))(0.0), (decltype(p))(-1.0), (decltype(p))(1.0));
		val *= valMult;
		return val;
	};



	//fast new version
	auto tple = ApplyOperationAndFilter(lambdaMain, isOuterRangelambda, X);

	auto& res = std::get<0>(tple);
	auto& outside = std::get<1>(tple);
	auto& r = outside;
	auto sign = outside;

	auto sqrtMinusLog = [](auto p)
	{
		auto  r = min(p, 1 - p);
		r = sqrt(-log(r));
		return r;
	};

	ApplyUnitaryOperation(sqrtMinusLog, r);

	auto calcSignMultiplier = [](auto p) { const decltype(p) half = 0.5; const decltype(p) zero = 0.0; const decltype(p) one = 1.0;
	auto   q = p - half; return iff(q < zero, -one, one); };

	ApplyUnitaryOperation(calcSignMultiplier, X);// sign);
	auto lowerUpperFilter = [](auto r) { return r <= 5.0; };
	auto lowerUpperTpl = ApplyBinaryFilter(lowerUpperFilter, r);
	auto& lower = std::get<0>(lowerUpperTpl);
	VecVW lower_p(lower, X);

	lower = ApplyBinaryOperation(dolambdaLow, lower, lower_p);
	lower.write(res);
	auto& upper = std::get<1>(lowerUpperTpl);
	if (upper.size() < 1)
		return res;

	VecVW upper_p(upper, X);
	upper = ApplyBinaryOperation(dolambdaHi, upper, upper_p);
	upper.write(res);
	return res;
}


template <typename VecXX>
VecXX calcCDFNormWichuraViewsAndFMA2splits(const VecXX& inputVecX)
{
	//Extrema boundary constants
	constexpr auto ExtrmMin =0.00000000001388;// exp(-25.0);
	constexpr auto ExtrmMax = 1. - ExtrmMin;

	// region test lambdas
	auto isExtremeLambda = [&](auto p)
	{
		return (p < ExtrmMin) || (p > ExtrmMax);
	};

	auto isOuterRangelambda = [&](auto p)
	{
		constexpr auto p_low = 0.5 - 0.425;
		constexpr auto p_high = 1.0 - p_low;
		return ((p_high < p) && (ExtrmMax > p)) || ((p > ExtrmMin) && (p < p_low));

	};

	//region evaluation lambdas
	auto centralRegionLambda = [](auto p)
	{
		auto q = p - 0.5;
		auto r = .180625 - q * q;
		auto denom = 1. / (mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(5226.495278852854561, r, 28729.085735721942674), r, 39307.89580009271061), r, 21213.794301586595867), r, 5394.1960214247511077), r, 687.1870074920579083), r, 42.313330701600911252), r, 1));
		auto num = (mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(2509.0809287301226727, r, 33430.575583588128105), r, 67265.770927008700853), r, 45921.953931549871457), r, 13731.693765509461125), r, 1971.5909503065514427), r, 133.14166789178437745), r, 3.387132872796366608) * q);
		return denom * num;
	};

	auto outerRegionlambda = [](auto p)
	{
		//wichura polynomial coefficients
		constexpr static  double c[] = { 7.7454501427834140764e-4 , .0227238449892691845833 ,.24178072517745061177, 1.27045825245236838258 ,  3.64784832476320460504, 5.7694972214606914055, 4.6303378461565452959, 1.42343711074968357734 };
		constexpr static  double d[] = { 1.05075007164441684324e-9 , 5.475938084995344946e-4, .0151986665636164571966, .14810397642748007459, .68976733498510000455,  1.6763848301838038494,  2.05319162663775882187,1. };
		auto  r = min(p, 1 - p);
		r = sqrt(-log(r));

		auto q = p - 0.5;
		r += -1.6;
		auto denom = (decltype(p))(1.0) / mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(d[0], r, d[1]), r, d[2]), r, d[3]), r, d[4]), r, d[5]), r, d[6]), r, (decltype(p))(1.0));
		auto	val = mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(c[0], r, c[1]), r, c[2]), r, c[3]), r, c[4]), r, c[5]), r, c[6]), r, c[7]);	
		val = val * denom;
		auto valMult = iff(q < (decltype(p))(0.0), (decltype(p))(-1.0), (decltype(p))(1.0));
		val *= valMult;
		return val;
	};


	auto extremaRegionLambda = []( auto p)
	{
		//wichura polynomial coefficients
		constexpr static  double e[] = { 2.01033439929228813265e-7 ,   2.71155556874348757815e-5,   .0012426609473880784386, .026532189526576123093, .29656057182850489123,   1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772 };
		constexpr static  double f[] = { 2.04426310338993978564e-15 , 1.4215117583164458887e-7, 1.8463183175100546818e-5,  7.868691311456132591e-4, .0148753612908506148525,.13692988092273580531, .59983220655588793769, 1. };

		auto  r = min(p, 1 - p);
		r = sqrt(-log(r));

		//auto q = p - 0.5;
		// very close to  0 or 1 
		r += -5.;

		auto	val = (((((((r * e[0] + e[1]) * r + e[2]) * r + e[3]) * r + e[4]) * r + e[5]) * r + e[6]) * r + e[7])
			/ (((((((r * f[0] + f[1]) * r + f[2]) * r + f[3]) * r + f[4]) * r + f[5]) * r + f[6]) * r + 1.);

		auto valMult = iff((p - 0.5) < (decltype(p))(0.0), (decltype(p))(-1.0), (decltype(p))(1.0));
		val *= valMult;

		return val;

	};


    // Apply central region lambda to all elements and filter outer range and extrema range elements to views
	// return all inside a tupple. 
	auto tple = ApplyOperationAndFilter(centralRegionLambda, isOuterRangelambda, isExtremeLambda, inputVecX);

	auto& res = std::get<0>(tple);  // main result initially filed with values from applying central Region Lambda
	auto& outside = std::get<1>(tple); // view containing outer range elements
	auto& extreme = std::get<2>(tple); // view containing extreme range elements ( usually empty)

	// use outerRegionLambda to transform  filtered  outer range values in  the view (outside) and write the results directly to result object res
	ApplyUnitaryOperationWrite(outerRegionlambda, outside, res);
    if (extreme.size() < 1)
	{
		return res;
	}

	//  transform values filtered to the extrema view using  lambda  extremaRegionLambda  and write transformed values to res
 	ApplyUnitaryOperationWrite(extremaRegionLambda, extreme, res);
	return res;
}










template <typename VecXX>
VecXX calcCDFNormsSparse(VecXX& X)
{

	static double a[] = { 0.0,  -3.969683028665376e+01, 2.209460984245205e+02,-2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00 };
	static double b[] = { 0.0, -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,  6.680131188771972e+01, -1.328068155288572e+01 };
	static double c[] = { 0.0,-7.784894002430293e-03,-3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
	static double d[] = { 0.0,  7.784695709041462e-03, 3.224671290700398e-01,  2.445134137142996e+00, 3.754408661907416e+00 };

	auto aclambdaMain = [=](auto p)
	{
		auto X = p;
		auto q = p - 0.5;
		auto r = q * q;
		X = (((((a[1] * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * r + a[6]) * q /
			(((((b[1] * r + b[2]) * r + b[3]) * r + b[4]) * r + b[5]) * r + 1);

		return X;
	};




	auto aclambdaLow = [=](auto initVal, auto p)
	{
		const auto p_low = 0.02425;
		auto condLo = (0.0 < p) && (p < p_low);

		if (!horizontal_or(condLo))
			return initVal;

		auto q = sqrt(-2.0 * log(p));
		auto X = (((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
			((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1.0);

		return select(condLo, X, initVal);

	};


	auto aclambdaHi = [=](auto initVal, auto p)
	{
		const auto p_low = 0.02425;
		const auto p_high = 1 - p_low;
		auto condHi = (p_high < p) && (p < 1);
		if (!horizontal_or(condHi))
			return initVal;

		auto q = sqrt(-2.0 * log(1 - p));
		const auto X = -(((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6])
			/ ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1.0);
		return select(condHi, X, initVal);
	};


	auto res = ApplyUnitaryOperation1(X, aclambdaMain);
	SparseUpdateWithLambda1(res, X, aclambdaLow);
	SparseUpdateWithLambda1(res, X, aclambdaHi);
	return res;

}


template <typename VecXX>
VecXX calcCDFNormsSparseFMA(VecXX& X)
{
	static double a[] = { 0.0,  -3.969683028665376e+01, 2.209460984245205e+02,-2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00 };
	static double b[] = { 0.0, -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,  6.680131188771972e+01, -1.328068155288572e+01 };
	static double c[] = { 0.0,-7.784894002430293e-03,-3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
	static double d[] = { 0.0,  7.784695709041462e-03, 3.224671290700398e-01,  2.445134137142996e+00, 3.754408661907416e+00 };

	auto aclambdaMain = [=](auto p)
	{
		auto X = p;
		auto q = p - 0.5;
		auto r = q * q;

		X = (mul_add(mul_add(mul_add(mul_add(mul_add(a[1], r, a[2]), r, a[3]), r, a[4]), r, a[5]), r, a[6]) * q) /
			mul_add(mul_add(mul_add(mul_add(mul_add(b[1], r, b[2]), r, b[3]), r, b[4]), r, b[5]), r, 1);

		return X;

	};


	auto aclambdaLow = [=](auto initVal, auto p)
	{
		const auto p_low = 0.02425;
		auto condLo = (0.0 < p) && (p < p_low);

		if (!horizontal_or(condLo))
			return initVal;

		auto q = sqrt(-2.0 * log(p));
		auto X = mul_add(mul_add(mul_add(mul_add(mul_add(c[1], q, c[2]), q, c[3]), q, c[4]), q, c[5]), q, c[6]) 
			/ mul_add(mul_add(mul_add(mul_add(d[1], q, d[2]), q, d[3]), q, d[4]), q, 1.0);
		return select(condLo, X, initVal);

	};


	auto aclambdaHi = [=](auto initVal, auto p)
	{
		const auto p_low = 0.02425;
		const auto p_high = 1 - p_low;
		auto condHi = (p_high < p) && (p < 1);
		if (!horizontal_or(condHi))
			return initVal;

		auto q = sqrt(-2.0 * log(1 - p));
		const auto X = -(((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6])
			/ ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1.0);
		return select(condHi, X, initVal);
	};


	auto res = ApplyUnitaryOperation1(X, aclambdaMain);
	SparseUpdateWithLambda1(res, X, aclambdaLow);
	SparseUpdateWithLambda1(res, X, aclambdaHi);
	return res;

};


template <typename VecXX>
VecXX calcCDFNormsSparseFMA2(VecXX& X)
{

	/// acklams inverse cdf normal
	static double a[] = { 0.0,  -3.969683028665376e+01, 2.209460984245205e+02,-2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00 };
	static double b[] = { 0.0, -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,  6.680131188771972e+01, -1.328068155288572e+01 };
	static double c[] = { 0.0,-7.784894002430293e-03,-3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
	static double d[] = { 0.0,  7.784695709041462e-03, 3.224671290700398e-01,  2.445134137142996e+00, 3.754408661907416e+00 };

	auto aclambdaMain = [=](auto p)
	{
		auto X = p;
		auto q = p - 0.5;
		auto r = q * q;

		X = (mul_add(mul_add(mul_add(mul_add(mul_add(a[1], r, a[2]), r, a[3]), r, a[4]), r, a[5]), r, a[6]) * q) /
			mul_add(mul_add(mul_add(mul_add(mul_add(b[1], r, b[2]), r, b[3]), r, b[4]), r, b[5]), r, 1);

		return X;

	};


	auto isLow = [](auto p) { const auto p_low = 0.02425; return  (0.0 < p) && (p < p_low); };
	auto doLow = [&](auto p)
	{
		auto q = sqrt(-2.0 * log(p));
		auto X = mul_add(mul_add(mul_add(mul_add(mul_add(c[1], q, c[2]), q, c[3]), q, c[4]), q, c[5]), q, c[6])
			/ mul_add(mul_add(mul_add(mul_add(d[1], q, d[2]), q, d[3]), q, d[4]), q, 1.0);
		return X; };



	auto isHi = [](auto x) { const auto p_low =  1.0 -0.02425; return  (x >= p_low) && (x  < 1.0); };

	auto doHi = [&](auto p)
	{
		auto q = sqrt(-2.0 * log(1 - p));
		const auto X = -(((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6])
			/ ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1.0);
		return X;
	};

	auto res = ApplyUnitaryOperation1(X, aclambdaMain);
	ApplySparseUnitaryOperationU(X, res, doLow, isLow);
	ApplySparseUnitaryOperationU(X, res, doHi, isHi);
	return res;

};



template <typename VecXX>
VecXX calcCDFNormsSparseFMAOnePass(VecXX& X)
{

	/// acklams inverse cdf normal
	const static double a[] = { 0.0,  -3.969683028665376e+01, 2.209460984245205e+02,-2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00 };
	const static double b[] = { 0.0, -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,  6.680131188771972e+01, -1.328068155288572e+01 };
	const static double c[] = { 0.0,-7.784894002430293e-03,-3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
	const static double d[] = { 0.0,  7.784695709041462e-03, 3.224671290700398e-01,  2.445134137142996e+00, 3.754408661907416e+00 };


	auto aclambdaLow = [=](auto condLo, auto initVal, auto p)
	{
		auto q = sqrt(-2.0 * log(p));
		auto X = mul_add(mul_add(mul_add(mul_add(mul_add(c[1], q, c[2]), q, c[3]), q, c[4]), q, c[5]), q, c[6]) / mul_add(mul_add(mul_add(mul_add(d[1], q, d[2]), q, d[3]), q, d[4]), q, 1.0);
		return select(condLo, X, initVal);

	};


	auto aclambdaHi = [=](auto condHi, auto initVal, auto p)
	{
		auto q = sqrt(-2.0 * log(1 - p));
		const auto X = -(mul_add(mul_add(mul_add(mul_add(mul_add(c[1], q, c[2]), q, c[3]), q, c[4]), q, c[5]), q, c[6]))
			/ mul_add(mul_add(mul_add(mul_add(d[1], q, d[2]), q, d[3]), q, d[4]), q, 1.0);
		return select(condHi, X, initVal);
	};


	auto aclam_sparseOnePass = [=](auto p)
	{
		auto X = p;
		auto q = p - 0.5;
		auto r = q * q;

		const auto p_low = 0.02425;
		const auto p_high = 1 - p_low;

		auto aboveLow = (p_low < p);
		auto belowHi = (p < p_high);
		auto condAllDone = aboveLow && belowHi;


		X = (mul_add(mul_add(mul_add(mul_add(mul_add(a[1], r, a[2]), r, a[3]), r, a[4]), r, a[5]), r, a[6]) * q) /
			mul_add(mul_add(mul_add(mul_add(mul_add(b[1], r, b[2]), r, b[3]), r, b[4]), r, b[5]), r, 1);

		if (horizontal_and(condAllDone))
		{
			return X;
		}

		if (!horizontal_and(belowHi))
		{
			X = aclambdaHi(belowHi, X, p);
		}
		if (!horizontal_and(aboveLow))
		{
			X = aclambdaLow(aboveLow, X, p);
		}
		return X;
	};

	auto res =ApplyTransformUR_XX(X, aclam_sparseOnePass);
	return res;

}


template <typename VecXX>
VecXX calcCDFNormWithViews(VecXX& X)
{

	/// acklams inverse cdf normal
	static double a[] = { 0.0,  -3.969683028665376e+01, 2.209460984245205e+02,-2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00 };
	static double b[] = { 0.0, -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,  6.680131188771972e+01, -1.328068155288572e+01 };
	static double c[] = { 0.0,-7.784894002430293e-03,-3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
	static double d[] = { 0.0,  7.784695709041462e-03, 3.224671290700398e-01,  2.445134137142996e+00, 3.754408661907416e+00 };

	auto aclambdaMain = [=](auto p)
	{
		auto X = p;
		auto q = p - 0.5;
		auto r = q * q;
		X = (((((a[1] * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * r + a[6]) * q /
			(((((b[1] * r + b[2]) * r + b[3]) * r + b[4]) * r + b[5]) * r + 1);

		return X;
	};

	auto islambdaLow = [=](auto p)
	{
		const auto p_low = 0.02425;
		return (0.0 < p) && (p < p_low);
	};

	auto islambdaHi = [=](auto p)
	{
		const auto p_low = 0.02425;
		const auto p_high = 1.0 - p_low;
		return (p_high < p) && (p < 1.0);
	};


	auto dolambdaLow = [=](auto p)
	{

		auto q = sqrt(-2.0 * log(p));
		return (((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
			((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1.0);

	};

	auto dolambdaHi = [=](auto p)
	{
		auto q = sqrt(-2.0 * log(1 - p));
		return -(((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
			((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1.0);
	};


	VecXX res = ApplyUnitaryOperation1(X, aclambdaMain);
	auto lambada = getLambdaBool(islambdaLow);
	auto lowVw = ApplyFilter(lambada, X); // a view 
	ApplyUnitaryOperation(lowVw,dolambdaLow);
	lowVw.write(res);

	auto lambadaHi = getLambdaBool(islambdaHi);
	auto hiVw = ApplyFilter(lambadaHi, X); // a view 
	ApplyUnitaryOperation(hiVw,dolambdaHi);
	hiVw.write(res);


	return res;

}


template <typename VecXX>
VecXX calcCDFNormWithViewsAndFMA(VecXX& X)
{

	/// acklams inverse cdf normal
	const static double a[] = { 0.0,  -3.969683028665376e+01, 2.209460984245205e+02,-2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00 };
	const static double b[] = { 0.0, -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,  6.680131188771972e+01, -1.328068155288572e+01 };
	const static double c[] = { 0.0,-7.784894002430293e-03,-3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
	const static double d[] = { 0.0,  7.784695709041462e-03, 3.224671290700398e-01,  2.445134137142996e+00, 3.754408661907416e+00 };

	/////////FMA Lambda ///////////////////////////
	auto aclambdaMain = [=](auto p)
	{
		auto X = p;
		auto q = p - 0.5;
		auto r = q * q;
		X = (mul_add(mul_add(mul_add(mul_add(mul_add(a[1], r, a[2]), r, a[3]), r, a[4]), r, a[5]), r, a[6]) * q) /
			mul_add(mul_add(mul_add(mul_add(mul_add(b[1], r, b[2]), r, b[3]), r, b[4]), r, b[5]), r, 1);

		return X;

	};

	auto islambdaLow = [=](auto p)
	{
		const auto p_low = 0.02425;
		return (0.0 < p) && (p < p_low);
	};


	auto isOuterRangelambda = [=](auto p)
	{
		const auto p_low = 0.02425;
		const auto p_high = 1.0 - p_low;
		return (p_high < p) || (p < p_low);
	};


	auto dolambdaLow = [=](auto p)
	{
		auto q = sqrt(-2.0 * log(p));
		return mul_add(mul_add(mul_add(mul_add(mul_add(c[1], q, c[2]), q, c[3]), q, c[4]), q, c[5]), q, c[6])
			/ mul_add(mul_add(mul_add(mul_add(d[1], q, d[2]), q, d[3]), q, d[4]), q, 1.0);


	};

	auto dolambdaHi = [=](auto p)
	{
		auto q = sqrt(-2.0 * log(1 - p));
		return  -(mul_add(mul_add(mul_add(mul_add(mul_add(c[1], q, c[2]), q, c[3]), q, c[4]), q, c[5]), q, c[6]))
			/ mul_add(mul_add(mul_add(mul_add(d[1], q, d[2]), q, d[3]), q, d[4]), q, 1.0);
	};

	//fast new version
	auto tple = ApplyOperationAndFilter(aclambdaMain, isOuterRangelambda, X);

	auto& res = std::get<0>(tple);
	auto& outside = std::get<1>(tple);

	auto loHiVwTple = ApplyBinaryFilter(islambdaLow, outside); // a view 

	auto& lowVw = std::get<0>(loHiVwTple);
	ApplyUnitaryOperation(lowVw,dolambdaLow );
	lowVw.write(res);
	auto& hiVw = std::get<1>(loHiVwTple);
	ApplyUnitaryOperation(hiVw,dolambdaHi);
	hiVw.write(res);
	return res;

}


template <typename VecXX>
VecXX calcCDFNormWithViewsAndFMAWriteFromViewCalc(VecXX& X)
{

	/// acklams inverse cdf normal
	const static double a[] = { 0.0,  -3.969683028665376e+01, 2.209460984245205e+02,-2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00 };
	const static double b[] = { 0.0, -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,  6.680131188771972e+01, -1.328068155288572e+01 };
	const static double c[] = { 0.0,-7.784894002430293e-03,-3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
	const static double d[] = { 0.0,  7.784695709041462e-03, 3.224671290700398e-01,  2.445134137142996e+00, 3.754408661907416e+00 };

	/////////FMA Lambda ///////////////////////////
	auto aclambdaMain = [=](auto p)
	{
		auto X = p;
		auto q = p - VecXX::scalar(0.5);
		auto r = q * q;

		X = (mul_add(mul_add(mul_add(mul_add(mul_add(a[1], r, a[2]), r, a[3]), r, a[4]), r, a[5]), r, a[6]) * q) /
			mul_add(mul_add(mul_add(mul_add(mul_add(b[1], r, b[2]), r, b[3]), r, b[4]), r, b[5]), r, 1);

		return X;

	};


	auto islambdaLow = [=](auto p)
	{
		const auto p_low = VecXX::scalar(0.02425);
		return (VecXX::scalar(0.0) < p) && (p < p_low);
	};



	auto isOuterRangelambda = [=](auto p)
	{
		const auto p_low =  VecXX::scalar(0.02425);
		const auto p_high = VecXX::scalar(1.0) - p_low;

		return (p_high < p) || (p < p_low);
	};



	auto dolambdaLow = [=](auto p)
	{
		auto q = sqrt(VecXX::scalar(-2.0) * log(p));
		return mul_add(mul_add(mul_add(mul_add(mul_add(c[1], q, c[2]), q, c[3]), q, c[4]), q, c[5]), q, c[6]) 
			/ mul_add(mul_add(mul_add(mul_add(d[1], q, d[2]), q, d[3]), q, d[4]), q, 1.0);
	};

	auto dolambdaHi = [=](auto p)
	{
		auto q = sqrt(-2.0 * log(1 - p));
		return  -(mul_add(mul_add(mul_add(mul_add(mul_add(c[1], q, c[2]), q, c[3]), q, c[4]), q, c[5]), q, c[6]))
			/ mul_add(mul_add(mul_add(mul_add(d[1], q, d[2]), q, d[3]), q, d[4]), q, 1.0);
	};


	auto tple = ApplyOperationAndFilter(aclambdaMain, isOuterRangelambda, X);

	auto& res = std::get<0>(tple);
	auto& outside = std::get<1>(tple);
	auto loHiVwTple = ApplyBinaryFilter(islambdaLow, outside); // a view 
	auto& lowVw = std::get<0>(loHiVwTple);
	ApplyUnitaryOperationWrite(dolambdaLow, lowVw, res);
	auto& hiVw = std::get<1>(loHiVwTple);
	ApplyUnitaryOperationWrite(dolambdaHi, hiVw, res);
	return res;

}

//float version 
template <typename VecXX>
VecXX calcCDFNormBranchlessFloat(VecXX& X)
{


	auto lambdaMain = [](auto u)
	{

		auto ushift = u - 0.5f;
		iff(ushift > 0.0f, 1.0f - u, u);

		typename VecXX::INS v = -log(u + u);
		typename  VecXX::INS p = 1.68267776058639e-6f;
		p = p * v + 0.0007404314351202936f;
		p = p * v + 0.03602364419560667f;
		p = p * v + 0.4500443083534446f;
		p = p * v + 1.861100468283588f;
		p = p * v + 2.748475794390544f;
		p = p * v + 1.253314132218524f;
		typename  VecXX::INS q = 0.00003709787159774307f;
		q = q * v + 0.004513659269519104f;
		q = q * v + 0.1101701640048184f;
		q = q * v + 0.8410203004476538f;
		q = q * v + 2.402969434512837f;
		q = q * v + 2.692965915568952f;
		q = q * v + 1.0f;

		static const typename  VecXX::INS one(1.0);
		static const typename  VecXX::INS minus_one(-1.0);
		static const typename  VecXX::INS zero(0.0);

		typename VecXX::INS dirnMult = select(v < zero, v, -v);
		typename VecXX::INS dirnMult_ushift = select(u < 0.5, one, minus_one);

		return dirnMult_ushift *dirnMult * p / q;
	};


	auto res = ApplyUnitaryOperation(X, lambdaMain);
	return res;

}




