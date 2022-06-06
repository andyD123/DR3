/****************************  apply_operation.h  *******************************
* Author:        Andrew Drakeford
* Date created:  2021-04-10
* Last modified: 2021-04-10
* Version:       1.0
* Project:       DR Cubed
* Description:
*
* (c) Copyright 2019 Andrew Drakeford
* Apache License version 2.0 or later.
*****************************************************************************/
#pragma once
#include "vec.h"
#include "vec_double.h"
#include "instruction_traits.h"
#include "boolean_operations.h"
#include "accumulate_transform.h"
#include "binary_unitary_operations.h"
#include "math_ops.h"
#include "filter_select.h"
#include "conditional_select_eval.h"
#include "vec_view.h"
#include "vcl_latest.h"



template<typename INS_VEC>
static INS_VEC cdfnormD(INS_VEC x)
{

	//   https://mathworld.wolfram.com/Erfc.html
	constexpr typename  InstructionTraits<INS_VEC>::FloatType invRootPi = 0.564189583547756;
	constexpr typename InstructionTraits<INS_VEC>::FloatType invRootTwo = 0.707106781186548;
	return invRootTwo * invRootPi*exp(-0.5*x*x);
}




template<typename INS_VEC>
static INS_VEC cdfnorm(const INS_VEC& z)  
{

	//   https://mathworld.wolfram.com/Erfc.html
	INS_VEC b1 = 0.31938153;
	INS_VEC b2 = -0.356563782;
	INS_VEC b3 = 1.781477937;
	INS_VEC b4 = -1.821255978;
	INS_VEC b5 = 1.330274429;
	INS_VEC p = 0.2316419;
	INS_VEC c2 = 0.3989423;
	auto x = select(z > (const INS_VEC)(6.0), (INS_VEC)(1.0), z); // this guards against overflow
	x = x;
	auto y = select(z < (INS_VEC )-6.0, (INS_VEC)0.0, z);
	y = y;
	INS_VEC a = abs(z);
	INS_VEC t = 1.0 / (1.0 + a*p);
	INS_VEC b = c2*exp((-z)*(z / 2.0));
	INS_VEC n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b*n;
	n = select(z < 0.0, 1.0 - n,n);
	return n;
}






template<typename INS_VEC>
VecD<INS_VEC> cdfnorm(const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(cdfnorm(rhs.value()), rhs.derivative()*cdfnormD(rhs.value()));
}

//to do replace with WS 16 digit impl
template<typename INS_VEC>
Vec<INS_VEC> cdfnorminv(const Vec<INS_VEC>& X)
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
		const auto X = -(((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
			((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1.0);
		return select(condHi, X, initVal);
	};



	auto res = ApplyUnitaryOperation1(X, aclambdaMain);
	SparseUpdateWithLambda1(res, X, aclambdaLow);
	SparseUpdateWithLambda1(res, X, aclambdaHi);

	return res;
}
//

/*
template <typename VecXX>
VecXX cdfnorminv(VecXX& X)
{
	//wichura
	const static  double a[] = { 2509.0809287301226727 , 33430.575583588128105, 67265.770927008700853, 45921.953931549871457, 13731.693765509461125,  1971.5909503065514427, 133.14166789178437745,3.387132872796366608 };
	const static  double b[] = { 5226.495278852854561, 28729.085735721942674,    39307.89580009271061, 21213.794301586595867, 5394.1960214247511077,   687.1870074920579083, 42.313330701600911252 };
	const static  double c[] = { 7.7454501427834140764e-4 , .0227238449892691845833 ,.24178072517745061177, 1.27045825245236838258 ,  3.64784832476320460504, 5.7694972214606914055, 4.6303378461565452959, 1.42343711074968357734 };
	const static  double d[] = { 1.05075007164441684324e-9 , 5.475938084995344946e-4, .0151986665636164571966, .14810397642748007459, .68976733498510000455,  1.6763848301838038494,  2.05319162663775882187,1. };
	const  static double e[] = { 2.01033439929228813265e-7 ,   2.71155556874348757815e-5,   .0012426609473880784386, .026532189526576123093, .29656057182850489123,   1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772 };
	const static  double f[] = { 2.04426310338993978564e-15 , 1.4215117583164458887e-7, 1.8463183175100546818e-5,  7.868691311456132591e-4, .0148753612908506148525,.13692988092273580531, .59983220655588793769, 1. };


	constexpr auto ExtrmMin = 0.00000000001388;// exp(-25.0);
	constexpr auto ExtrmMax = 1 - ExtrmMin;

	/////////FMA Lambda ///////////////////////////
	auto aclambdaMain = [=](const auto& p)
	{

		auto q = p - 0.5;
		auto r = .180625 - q * q;

		auto denom = 1. / (mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(5226.495278852854561, r, 28729.085735721942674), r, 39307.89580009271061), r, 21213.794301586595867), r, 5394.1960214247511077), r, 687.1870074920579083), r, 42.313330701600911252), r, 1));
		auto num = (mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(2509.0809287301226727, r, 33430.575583588128105), r, 67265.770927008700853), r, 45921.953931549871457), r, 13731.693765509461125), r, 1971.5909503065514427), r, 133.14166789178437745), r, 3.387132872796366608) * q);

		return denom * num;

	};


	auto doExtrema = [=](auto p)
	{

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

	auto isExtremeLambda = [=](auto p)
	{
		return (p < ExtrmMin) || (p > ExtrmMax);
	};

	auto isOuterRangelambda = [=](auto p)
	{
		constexpr auto p_low = 0.5 - 0.425;
		constexpr auto p_high = 1.0 - p_low;

		return ((p_high < p) && (ExtrmMax > p)) || ((p > ExtrmMin) && (p < p_low));

	};

	auto dolambdaLowNew = [=](auto p)
	{
		auto  r = min(p, 1 - p);
		r = sqrt(-log(r));

		auto q = p - 0.5;
		r += -1.6;
		auto	val = mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(c[0], r, c[1]), r, c[2]), r, c[3]), r, c[4]), r, c[5]), r, c[6]), r, c[7]);
		auto denom = (decltype(p))(1.0) / mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(d[0], r, d[1]), r, d[2]), r, d[3]), r, d[4]), r, d[5]), r, d[6]), r, (decltype(p))(1.0));
		val = val * denom;

		auto valMult = iff(q < (decltype(p))(0.0), (decltype(p))(-1.0), (decltype(p))(1.0));
		val *= valMult;
		return val;
	};

	auto tple = ApplyOperationAndFilter2W(aclambdaMain, isOuterRangelambda, isExtremeLambda, X);

	auto& res = std::get<0>(tple);
	auto& outside = std::get<1>(tple);
	auto& extreme = std::get<2>(tple);

	ApplyUnitaryOperationXWrite(dolambdaLowNew, outside, res);
	if (extreme.size() < 1)
	{
		return res;
	}

	ApplyUnitaryOperationXWrite(doExtrema, extreme, res);
	return res;

}
//
*/


/*

// @WichuraQuantile   
// see R implementation KCL implementation
double cdfnorminv( double p)
{

	static long double a[] = { 2509.0809287301226727 , 33430.575583588128105, 67265.770927008700853, 45921.953931549871457, 13731.693765509461125,  1971.5909503065514427, 133.14166789178437745,3.387132872796366608 };
	static long double b[] = { 5226.495278852854561, 28729.085735721942674,    39307.89580009271061, 21213.794301586595867, 5394.1960214247511077,   687.1870074920579083, 42.313330701600911252 };
	static long double c[] = { 7.7454501427834140764e-4 , .0227238449892691845833 ,.24178072517745061177, 1.27045825245236838258 ,  3.64784832476320460504, 5.7694972214606914055, 4.6303378461565452959, 1.42343711074968357734 };
	static long double d[] = { 1.05075007164441684324e-9 , 5.475938084995344946e-4, .0151986665636164571966, .14810397642748007459, .68976733498510000455,  1.6763848301838038494,  2.05319162663775882187,1. };
	static long double e[] = { 2.01033439929228813265e-7 ,   2.71155556874348757815e-5,   .0012426609473880784386, .026532189526576123093, .29656057182850489123,   1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772 };
	static long double f[] = { 2.04426310338993978564e-15 , 1.4215117583164458887e-7, 1.8463183175100546818e-5,  7.868691311456132591e-4, .0148753612908506148525,.13692988092273580531, .59983220655588793769, 1. };


	long double val = 0.0;
	long double   q = p - 0.5;


	if (fabs(q) <= .425)
	{
		long double r = .180625 - q * q;
		val =

			q * (((((((r * a[0] +
				a[1]) * r + a[2]) * r +
				a[3]) * r + a[4]) * r +
				a[5]) * r + a[6]) * r +
				a[7])
			/ (((((((r * b[0] +
				b[1]) * r + b[2]) * r +
				b[3]) * r + b[4]) * r +
				b[5]) * r + b[6]) * r + 1.);

	}
	else
	{ // closer than 0.075 from {0,1} boundary 

	   // r = min(p, 1-p) < 0.075 
		long double  r = std::min(p, 1 - p);
		r = sqrt(-log(r));



		// <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 
		if (r <= 5.)
		{
			r += -1.6;
			val = (((((((r * c[0] + c[1]) * r + c[2]) * r + c[3]) * r + c[4]) * r + c[5]) * r + c[6]) * r + c[7])
				/ (((((((r * d[0] + d[1]) * r + d[2]) * r + d[3]) * r + d[4]) * r + d[5]) * r + d[6]) * r + 1.);


		}
		else
		{ // very close to  0 or 1 
			r += -5.;

			val = (((((((r * e[0] + e[1]) * r + e[2]) * r + e[3]) * r + e[4]) * r + e[5]) * r + e[6]) * r + e[7])
				/ (((((((r * f[0] + f[1]) * r + f[2]) * r + f[3]) * r + f[4]) * r + f[5]) * r + f[6]) * r + 1.);

		}

		long double valMult = (q < 0.0) ? -1.0 : 1.0;
		val *= valMult;
	}

	return val;
}

*/