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
		if ((width > stride) || (stride % width > 0)) throw std::exception("bad stride size for width");

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
