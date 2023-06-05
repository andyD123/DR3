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





template<typename INS_VEC, int OFFSET>
struct SampleElement :public RegisterElement< INS_VEC, OFFSET, false>
{

	template <int VAL>
	INS_VEC  get() {};// cant instantiate


	template<>
	INS_VEC  get<OFFSET>() { return ::m_offsetData1; };


};

template<typename INS_VEC, int X_Minus1 = -1, int X0 = 0, int X1 = 1 >
struct TrinomialSampler
{
	int stride() const { return 1; }

	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	template <int VAL>
	INS_VEC  get() {};// cant instantiate

	template<>
	INS_VEC  get<X1>() { return X_1.value; };

	template<>
	INS_VEC  get<X0>() { return X_0.value; };

	template<>
	INS_VEC  get<X_Minus1>() { return X_Minus_1.value; };


	inline void load(Float* pData)
	{
		X_1.load_u(pData);
		X_0.load_u(pData);
		X_Minus_1.load_u(pData);
	}

	static constexpr int max() { return std::max(X_Minus1, std::max(X0, X1)); }
	static constexpr int min() { return std::min(X_Minus1, std::min(X0, X1)); }


	RegisterElement< INS_VEC, X1, false> X_1;
	RegisterElement< INS_VEC, X0, false> X_0;
	RegisterElement< INS_VEC, X_Minus1, false> X_Minus_1;


};



template<typename INS_VEC, int X0 = 0, int X1 = 1>
struct BinomialSampler
{
	int stride() const { return 1; }

	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	template <int VAL>
	INS_VEC  get() {  };// cant instantiate

	template<>
	INS_VEC  get<X1>() { return X_1.value; };

	template<>
	INS_VEC  get<X0>() { return X_0.value; };



	inline void load(Float* pData)
	{
		X_1.load_u(pData);
		X_0.load_u(pData);
	}

	static constexpr int max() { return std::max(X0, X1); }
	static constexpr int min() { return std::min(X0, X1); }


	RegisterElement< INS_VEC, 1, false> X_1;
	RegisterElement< INS_VEC, 0, false> X_0;

};


template<typename INS_VEC, int X0 = 0>
struct UnitarySampler
{

	int stride() const { return 1; }

	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	template <int VAL>
	INS_VEC  get() {  };// cant instantiate


	template<>
	INS_VEC  get<X0>() { return X_0.value; };


	inline void load(Float* pData)
	{
		X_0.load_u(pData);
	}

	static constexpr int max() { return X0; }
	static constexpr int min() { return X0; }

	RegisterElement< INS_VEC, 0, false> X_0;

};


template<typename INS_VEC>
struct Convertable : public  UnitarySampler<INS_VEC, 0>
{

	//convert to ins_vec
	operator INS_VEC ()
	{
		return UnitarySampler<INS_VEC, 0>::X_0.value;
	}

};


//would be nice for variadic


template<typename INS_VEC, int X0 = 0>
struct StridedSampler
{

	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	static constexpr int  width = InstructionTraits< INS_VEC>::width;
	const size_t  m_stride;
	const int  step;

	size_t stride() const
	{
		return m_stride;
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


	explicit StridedSampler(size_t stride) :m_stride(stride), step(static_cast<int>(m_stride))
	{
		if ((width > stride) || (stride % width > 0)) throw std::exception("bad stride size for width");

	}

	template <int VAL>
	INS_VEC  get() {  };// cant instantiate


	template<>
	INS_VEC  get<X0>() { return X_0.value; };


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

	//static 
	constexpr int min() { return X0; }
	//static
	constexpr int max() { return X0; };// +(width - 1) * step;}


	RegisterElement< INS_VEC, 0, false> X_0;

};

