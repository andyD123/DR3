/****************************  sampler.h   *******************************
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
#include "instruction_traits.h"
#include "unroll_operators.h"

#include <stdexcept>
#include <tuple>


/*
Samplers introduce multiple variables corresponding to offset  values  into the lambda
 eg  X[i-1] x[0] X[i+1]  they use unaligned loads.

 would be nice for variadic form

 the template integers are the offsets 

*/

template<typename INS_VEC, int X_Minus1 = -1, int X0 = 0, int X1 = 1 >
struct TrinomialSampler :public std::tuple< RegisterElement< INS_VEC, X_Minus1,false>,
											RegisterElement< INS_VEC, X0,false>,
											RegisterElement< INS_VEC, X1,false> >
{
	int stride() const { return 1; }

	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	using TypeXMinus1 = RegisterElement< INS_VEC, X_Minus1,false>;
	using TypeX0 = RegisterElement< INS_VEC, X0,false>;
	using TypeX1 = RegisterElement< INS_VEC, X1,false>;


	TypeX1& X_1;
	TypeX0& X_0;
	TypeXMinus1& X_Minus_1;


	TrinomialSampler() :
		X_1(std::get<2>(*this)),
		X_0(std::get<1>(*this)),
		X_Minus_1(std::get<0>(*this))
	{

	}

	enum class  DIR { DOWN = 0, MID = 1, UP = 2 };

	inline void load(Float* pData)
	{
		X_1.load_u(pData);
		X_0.load_u(pData);
		X_Minus_1.load_u(pData);
	}

	static constexpr int max() { return std::max(X_Minus1, std::max(X0, X1)); }
	static constexpr int min() { return std::min(X_Minus1, std::min(X0, X1)); }

};

template<typename INS_VEC, int ITEM, int X_Minus1 = -1, int X0 = 0, int X1 = 1 >
INS_VEC get(TrinomialSampler< INS_VEC, X_Minus1, X0, X1>& sampler)
{
	return get<ITEM>(sampler);
}




template<typename INS_VEC, int X0 = 0, int X1 = 1>
struct BinomialSampler :public std::tuple< 
	RegisterElement< INS_VEC, X0, false>,
	RegisterElement< INS_VEC, X1, false> >
{
	int stride() const { return 1; }

	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	using TypeX0 = RegisterElement< INS_VEC, X0, false>;
	using TypeX1 = RegisterElement< INS_VEC, X1, false>;

	TypeX1& X_1;
	TypeX0& X_0;
	
	BinomialSampler() :
		X_1(std::get<1>(*this)),
		X_0(std::get<0>(*this))
	{
	}

	inline void load(Float* pData)
	{
		X_1.load_u(pData);
		X_0.load_u(pData);
	}

	static constexpr int max() { return std::max(X0, X1); }
	static constexpr int min() { return std::min(X0, X1); }

};


template<typename INS_VEC, int X0 = 0>
struct UnitarySampler : public std::tuple< RegisterElement< INS_VEC, X0, false> >
{
	int stride() const { return 1; }

	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	using TypeX0 = RegisterElement< INS_VEC, X0, false>;
	TypeX0& X_0;

	UnitarySampler() :
		X_0(std::get<0>(*this))
	{

	}


	inline void load(Float* pData)
	{
		X_0.load_u(pData);
	}

	static constexpr int max() { return X0; }
	static constexpr int min() { return X0; }

};


/*
used for wrapping and sampler and making its loaded register available via 
conversion
*/
template<typename INS_VEC>
struct Convertable : public  UnitarySampler<INS_VEC, 0>
{

	//convert to ins_vec
	operator INS_VEC ()
	{
		return UnitarySampler<INS_VEC, 0>::X_0.value;
	}

};


template <typename CALL_1, typename CALL_2>
struct Overloaded :public CALL_1,  CALL_2
{
	Overloaded(CALL_1& c1, CALL_2& c2) :CALL_1(c1), CALL_2(c2) {}

	using CALL_1::operator();
	using CALL_2::operator();
};



template <typename CALL_1, typename CALL_2>
Overloaded< CALL_1,  CALL_2> makeOverloaded(CALL_1& c1,  CALL_2& c2)
{
	return Overloaded< CALL_1, CALL_2>(c1, c2);
};



/*
 fills a SIMD register with values spread at a fixed stride length from each other
 as in the case with a matrix

*/

template<typename INS_VEC, int X0 = 0>
struct StridedSampler : public std::tuple< RegisterElement< INS_VEC, X0, false> >
{

	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	using TypeX0 = RegisterElement< INS_VEC, X0, false>;
	
	TypeX0& X_0;
	static constexpr int  width = InstructionTraits< INS_VEC>::width;

	const size_t  m_stride;
	const int  step;

	size_t stride() const
	{
		return m_stride;
	}


	explicit StridedSampler(size_t stride) : X_0(std::get<0>(*this)),m_stride(stride), step(static_cast<int>(m_stride))
	{
		if ((width > stride) || (stride % width > 0)) throw std::runtime_error("bad stride size for width");

	}


	template <typename T>
	void load2(T* pbase)
	{
		INS_VEC loaded = { *pbase, *(pbase + step) };
		X_0.value = loaded;
	}

	template <typename T>
	void load4(T* pbase)
	{
		INS_VEC loaded = { *pbase, *(pbase + step), *(pbase + 2 * step),*(pbase + 3 * step) };
		X_0.value = loaded;
	}

	template <typename T>
	void load8(T* pbase)
	{
		INS_VEC loaded = { *pbase, *(pbase + step), *(pbase + 2 * step),*(pbase + 3 * step) ,
			*(pbase + 4 * step),*(pbase + 5 * step),*(pbase + 6 * step),*(pbase + 7 * step) };

		X_0.value = loaded;
	}

	template <typename T>
	void load16(T* pbase)
	{
		INS_VEC loaded = { *pbase, *(pbase + step), *(pbase + 2 * step),*(pbase + 3 * step) ,
			*(pbase + 4 * step),*(pbase + 5 * step),*(pbase + 6 * step),*(pbase + 7 * step),
			*(pbase + 8 * step), *(pbase + 9 * step), *(pbase + 10 * step), *(pbase + 11 * step),
			*(pbase + 12 * step), *(pbase + 13 * step), *(pbase + 14 * step), *(pbase + 15 * step) };

		X_0.value = loaded;
	}

	inline void load(Float* pData)
	{

		if constexpr (width == 2)
		{
			load2(pData);
		}
		if constexpr (width == 4)
		{
			load4(pData);
		}
		if constexpr (width == 8)
		{
			load8(pData);
		}
		if constexpr (width == 16)
		{
			load16(pData);
		}
	}

	constexpr int min() { return X0; }
	constexpr int max() { return X0; };

};


//experimental 
//strided span
template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP, typename SAMPLER>
typename InstructionTraits<INS_VEC>::FloatType ApplyAccumulate2UR_X_EX_STRD(const VEC_TYPE<INS_VEC>& rhs1, OP& oper, SAMPLER& sampler, long i = 0, int impSZ = -1)
{

	auto pRhs1 = rhs1.start();

	long sz  = (impSZ < 0) ? static_cast<long>(rhs1.paddedSize()) : impSZ;

	const int width = InstructionTraits<INS_VEC>::width;
	int stride = static_cast<int>(sampler.stride());
	int step = 4 * width * stride;


	SAMPLER RHS1(sampler);
	INS_VEC RES = 0.;

	SAMPLER RHS2(sampler);
	INS_VEC RES1 = 0.;

	SAMPLER RHS3(sampler);
	INS_VEC RES2 = 0.;

	SAMPLER RHS4(sampler);
	INS_VEC RES3 = 0.;

	//we can only get a starting position bigger than  zero when we access points in the 
	// data preceeding the starting point, so we advance to a popint where we sample valid /existing data
	i = i + std::max(0, -sampler.min());


	//similarly if we are sampling  points beyond current index, we need to reduce maximum value iterated to so
	// that we stay in a valid range 
	impSZ = impSZ - std::max(0, sampler.max());



	if (sz >= step * 2)
	{
		//initialise first set of registers
		{
			RHS1.load(pRhs1 + i);
			RES = RHS1.X_0.value;
			RHS2.load(pRhs1 + i + width*stride);
			RES1 = RHS2.X_0.value;
			RHS3.load(pRhs1 + i + width * 2 *stride);
			RES2 = RHS3.X_0.value;
			RHS4.load(pRhs1 + i + width * 3* stride);
			RES3 = RHS4.X_0.value;
		}

		i += step;
		long rhsSZ = static_cast<long>(rhs1.size());
		for (; i <= (rhsSZ - step); i += step)
		{
			RHS1.load(pRhs1 + i);
			RES = oper(RES, RHS1);

			RHS2.load(pRhs1 + i + width* stride);
			RES1 = oper(RES1, RHS2);

			RHS3.load(pRhs1 + i + width * stride * 2);
			RES2 = oper(RES2, RHS3);

			RHS4.load(pRhs1 + i + width * stride * 3);
			RES3 = oper(RES3, RHS4);

		}

		// odd bits whole register loads
		for (; i <= rhsSZ - stride * width; i += stride* width)
		{
			RHS1.load(pRhs1 + i);
			RES = oper(RES, RHS1);
		}

		RES = oper(RES, RES1);
		RES2 = oper(RES2, RES3);
		RES = oper(RES, RES2);

	}
	else if (sz >= stride * width * 2)
	{
		RHS2.load(pRhs1);
		RES = RHS2.X_0.value;

		i += width * stride;
		// odd bits
		for (; i <= sz - stride* width; i += stride* width)
		{
			RHS1.load(pRhs1 + i);
			RES = oper(RES, RHS1);
		}

	}

	typename InstructionTraits<INS_VEC>::FloatType result = RES[0];
	long min_wdth = std::min(sz,(long) width);
	//across vectors lanes  // not assuming horizontal versoion exist
	for (long j = 1; j < min_wdth; ++j)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(result, RES[j], oper);
	}

	//end bits for vecs not filling padding
	for (; i < sz; i+=stride)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(pRhs1[i], result, oper);
	}

	return result;
}


//experimental 
//strided span


template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP, typename TRANSFORM, typename SAMPLER>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate2UR_X_EX_STRD(const VEC_TYPE<INS_VEC>& rhs1, OP& oper, TRANSFORM& trans, SAMPLER& sampler, long i = 0, int impSZ = -1)
{

	auto pRhs1 = rhs1.start();

	long sz = (impSZ < 0) ? static_cast<long>(rhs1.paddedSize()) : impSZ;

	const int width = InstructionTraits<INS_VEC>::width;
	int stride = static_cast<int>(sampler.stride());
	int step = 4 * width * stride;


	SAMPLER RHS1(sampler);
	INS_VEC RES = 0.;

	SAMPLER RHS2(sampler);
	INS_VEC RES1 = 0.;

	SAMPLER RHS3(sampler);
	INS_VEC RES2 = 0.;

	SAMPLER RHS4(sampler);
	INS_VEC RES3 = 0.;

	//we can only get a starting position bigger than  zero when we access points in the 
	// data preceeding the starting point, so we advance to a popint where we sample valid /existing data
	i = i + std::max(0, -sampler.min());


	//similarly if we are sampling  points beyond current index, we need to reduce maximum value iterated to so
	// that we stay in a valid range 
	impSZ = impSZ - std::max(0, sampler.max());



	if (sz >= step * 2)
	{
		//initialise first set of registers
		{
			RHS1.load(pRhs1 + i);
			RES = trans(RHS1);
			RHS2.load(pRhs1 + i + width * stride);
			RES1 = trans(RHS2);
			RHS3.load(pRhs1 + i + width * 2 * stride);
			RES2 = trans(RHS3);
			RHS4.load(pRhs1 + i + width * 3 * stride);
			RES3 = trans(RHS4);
		}

		i += step;
		long rhsSZ = static_cast<long>(rhs1.size());
		for (; i <= (rhsSZ - step); i += step)
		{
			RHS1.load(pRhs1 + i);
			RES = oper(RES, trans(RHS1));

			RHS2.load(pRhs1 + i + width * stride);
			RES1 = oper(RES1, trans(RHS2));

			RHS3.load(pRhs1 + i + width * stride * 2);
			RES2 = oper(RES2, trans(RHS3));

			RHS4.load(pRhs1 + i + width * stride * 3);
			RES3 = oper(RES3, trans(RHS4));

		}

		// odd bits whole register loads
		for (; i <= rhsSZ - stride * width; i += stride * width)
		{
			RHS1.load(pRhs1 + i);
			RES = oper(RES, trans(RHS1));
		}

		RES = oper(RES, RES1);
		RES2 = oper(RES2, RES3);
		RES = oper(RES, RES2);

	}
	else if (sz >= stride * width * 2)
	{
		RHS2.load(pRhs1);
		RES = trans(RHS2);

		i += width * stride;
		// odd bits
		for (; i <= sz - stride * width; i += stride * width)
		{
			RHS1.load(pRhs1 + i);
			RES = oper(RES, trans(RHS1));
		}

	}

	typename InstructionTraits<INS_VEC>::FloatType result = RES[0];
	long min_wdth = std::min(sz, (long)width);
	//across vectors lanes  // not assuming horizontal versoion exist
	for (long j = 1; j < min_wdth; ++j)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(result, RES[j], oper);
	}

	//end bits for vecs not filling padding
	for (; i < sz; i += stride)
	{
		RHS1.X_0.value = pRhs1[i];
		typename InstructionTraits<INS_VEC>::FloatType transformed = trans(RHS1)[0];
		result = ApplyBinaryOperationVec<INS_VEC, OP>(transformed, result, oper);
	}

	return result;
}

