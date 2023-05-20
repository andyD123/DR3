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

#include <type_traits>



template<typename INS_VEC>
static INS_VEC cdfnormD(INS_VEC x)
{

	auto asNumber = [](auto x)
	{
		return static_cast<typename InstructionTraits<INS_VEC>::FloatType>(x);
	};

	//   https://mathworld.wolfram.com/Erfc.html
	constexpr typename  InstructionTraits<INS_VEC>::FloatType invRootPi = asNumber(0.564189583547756);
	constexpr typename InstructionTraits<INS_VEC>::FloatType invRootTwo =asNumber( 0.707106781186548);
	return invRootTwo * invRootPi*exp(-0.5*x*x);
}



template<typename INS_VEC>
static INS_VEC cdfnorm(const INS_VEC& z)  
{

	auto asNumber = [](auto x)
	{
		return static_cast<typename InstructionTraits<INS_VEC>::FloatType>(x);
	};

	auto asInsVec = [&](auto x){  return  INS_VEC(asNumber(x) ); };


	//   https://mathworld.wolfram.com/Erfc.html
	INS_VEC b1 = asInsVec(0.31938153);
	INS_VEC b2 = asInsVec(-0.356563782);
	INS_VEC b3 = asInsVec(1.781477937);
	INS_VEC b4 = asInsVec(-1.821255978);
	INS_VEC b5 = asInsVec(1.330274429);
	INS_VEC p = asInsVec(0.2316419);	
	INS_VEC c2 = asInsVec(0.3989423);

	const auto  cond1 = (z > asInsVec(6.0));
	INS_VEC x = select(cond1, asInsVec(1.0), z);
	x = x;
	
	INS_VEC y = select( (z < asInsVec(-6.0)),asInsVec(0.0), z);
	y = y;
	INS_VEC a = abs(z);
	INS_VEC t = asInsVec(1.0) / (asInsVec(1.0) + a*p);
	INS_VEC b = c2*exp((-z)*(z / asInsVec(2.0)));
	INS_VEC n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = asInsVec(1.0) - b*n;
	n = select( (z < asInsVec(0.0) ), asInsVec(1.0) - n,n);
	return n;
}



template<typename INS_VEC>
Vec<INS_VEC> cdfnorm(const Vec<INS_VEC>& X)

{

	auto centralLambda = [&](auto z)
	{

		constexpr double N[] = { 3.52624965998911e-02 , 0.700383064443688,   6.37396220353165, 33.912866078383,  112.079291497871,  221.213596169931, 220.206867912376 };
		constexpr double M[] = { 8.83883476483184e-02, 1.75566716318264, 16.064177579207, 86.7807322029461 , 296.564248779674,  637.333633378831, 793.826512519948,440.413735824752 };

		auto inv_dc = 1.0 / mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(M[0], z, M[1]), z, M[2]), z, M[3]), z, M[4]), z, M[5]), z, M[6]), z, M[7]);
		auto n_c = mul_add(mul_add(mul_add(mul_add(mul_add(mul_add(N[0], z, N[1]), z, N[2]), z, N[3]), z, N[4]), z, N[5]), z, N[6]);

		return n_c * inv_dc;
	};


	auto outerLambda = [](auto z)
	{
		constexpr double inv_RT2PI(0.39894228040143267793994605993438);
		constexpr double  d[] = { 20. , 13., 200., 78., 300., 39. };
		constexpr double  n[] = { 20., 13., 180., 65., 160. };

		auto d_outer = mul_add(mul_add(mul_add(mul_add(mul_add((d[0] * z), z, d[1]), z, d[2]), z, d[3]), z, d[4]), z, d[5]);
		auto inv_d_outer = inv_RT2PI / d_outer;

		auto n_outer = mul_add(mul_add(mul_add(mul_add((n[0] * z), z, n[1]), z, n[2]), z, n[3]), z, n[4]);
		return   n_outer * inv_d_outer;
	};



	auto onePass = [=](auto x)
	{
		auto z = abs(x);
		auto e = exp(-z * z * 0.5);
		auto central = centralLambda(z);
		auto  SPLIT = 7.42;// 7106781186547; // appears to give less error
		auto condAllDone = (x * x < SPLIT* SPLIT);

		if (horizontal_and(condAllDone))
		{
			central *= e;
			return select(x <= 0.0, central, 1.0 - central);
		}

		auto outer = outerLambda(z);
		auto RES = select((z < SPLIT), central, outer);
		RES *= e;
		return select(x <= 0.0, RES, 1.0 - RES);

	};

	return ApplyTransformUR_X(X, onePass);

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

	auto asNumber = [](auto x)
	{
		return static_cast<typename InstructionTraits<INS_VEC>::FloatType>(x);
	};


	/// acklams inverse cdf normal
	static typename InstructionTraits<INS_VEC>::FloatType a[] = { asNumber(0.0), asNumber( -3.969683028665376e+01), asNumber(2.209460984245205e+02), asNumber(-2.759285104469687e+02), asNumber(1.383577518672690e+02), asNumber(-3.066479806614716e+01) ,  asNumber(2.506628277459239e+00)};
	static typename InstructionTraits<INS_VEC>::FloatType b[] = { asNumber(0.0), asNumber(-5.447609879822406e+01),  asNumber(1.615858368580409e+02), asNumber(-1.556989798598866e+02), asNumber(6.680131188771972e+01), asNumber(-1.328068155288572e+01) };
	static typename InstructionTraits<INS_VEC>::FloatType c[] = { asNumber(0.0), asNumber(-7.784894002430293e-03), asNumber(-3.223964580411365e-01), asNumber(-2.400758277161838e+00), asNumber(-2.549732539343734e+00), asNumber(4.374664141464968e+00), asNumber(2.938163982698783e+00) };
	static typename InstructionTraits<INS_VEC>::FloatType d[] = { asNumber(0.0), asNumber(7.784695709041462e-03), asNumber(3.224671290700398e-01),  asNumber(2.445134137142996e+00), asNumber(3.754408661907416e+00) };

	auto aclambdaMain = [=](auto p)
	{
		auto X = p;
		auto q = p - asNumber(0.5);
		auto r = q * q;
		X = (((((a[1] * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * r + a[6]) * q /
			(((((b[1] * r + b[2]) * r + b[3]) * r + b[4]) * r + b[5]) * r + 1);

		return X;
	};


	auto aclambdaLow = [=](auto initVal, auto p)
	{
		const auto p_low = asNumber(0.02425);
		auto condLo = (asNumber(0.0) < p) && (p < p_low);

		if (!horizontal_or(condLo))
			return initVal;

		auto q = sqrt(-2.0 * log(p));
		auto X = (((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
			((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + asNumber(1.0));

		return select(condLo, X, initVal);

	};


	auto aclambdaHi = [=](auto initVal, auto p)
	{
		const auto p_low = asNumber(0.02425);
		const auto p_high = 1 - p_low;
		auto condHi = (p_high < p) && (p < 1);
		if (!horizontal_or(condHi))
			return initVal;

		auto q = sqrt(asNumber(-2.0) * log(asNumber(1.) - p));
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
