/****************************  instruction_traits.h   *******************************
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
#pragma warning(suppress:4984)

#include "vec.h"
#include "vec_double.h"
#include "../VCL/vectormath_common.h"


template<typename INS_VEC> class  Vec;
template<typename INS_VEC> class  VecBool;
template<typename INS_VEC> class  VecD;
template< typename INS_VEC> class VecView;


template< typename INS_VEC>
struct InstructionTraits
{
	using BoolType = VecBoolD;
	using FloatType = double;
	static const int width = 2;
	static const double nullValue;
	static const double oneValue;
	static constexpr  bool alignedLoadStore = false;
	static constexpr  bool boolTypeIsAlignedLoadStore = false;
	static constexpr  bool useScatter = false;

	static constexpr  bool isCompact = false;
	using RegBoolType = VecBoolD;
	using MemBoolType = VecDouble;

};




template<>
struct InstructionTraits<VecDouble>
{
	using IdxType = Vec2q;
	using BoolType = VecBoolD;
	using FloatType = double;
	static const int width = 2;
	static const double nullValue;
	static const double oneValue;
	static constexpr  bool alignedLoadStore = true;
	static constexpr  bool useScatter = false;
	static constexpr  bool boolTypeIsAlignedLoadStore = true;

	static constexpr  bool isCompact = false;
	using RegBoolType = VecBoolD;
	using MemBoolType = VecDouble;
};


template<>
struct InstructionTraits<Vec2d>
{
	using IdxType = Vec2q;//Vec2q;
	using BoolType = Vec2db;// VecBoolD;
	using FloatType = double;
	static const int width = 2;
	static const double nullValue;
	static const double oneValue;
	static constexpr  bool alignedLoadStore = true;
	static constexpr  bool boolTypeIsAlignedLoadStore = false;
	static constexpr  bool useScatter = false;

	static constexpr  bool isCompact = false;
	using RegBoolType = Vec2db;
	using MemBoolType = Vec2d;


};


template<>
struct InstructionTraits<Vec4f>
{
	using IdxType = Vec4i;
	using BoolType = Vec4fb;
	using FloatType = float;
	static const int width = 4;
	static const float nullValue;
	static constexpr  bool alignedLoadStore = true;
	static constexpr  bool boolTypeIsAlignedLoadStore = false;
	static const float oneValue;
	static constexpr  bool useScatter = false;

	static constexpr  bool isCompact = false;
	using RegBoolType = Vec4fb;
	using MemBoolType = Vec4f;
};





template<>
struct InstructionTraits<Vec4d>
{
	using IdxType = Vec4i;
	using BoolType = Vec4db;
	using FloatType = double;
	static const int width = 4;
	static const double nullValue;
	static constexpr  bool alignedLoadStore = true;
	static constexpr  bool boolTypeIsAlignedLoadStore = false;
	static const double oneValue;
	static constexpr  bool useScatter = true;

	static constexpr  bool isCompact = false;
	using RegBoolType = Vec4db;
	using MemBoolType = Vec4d;
};


template<>
struct InstructionTraits<Vec8f>
{
	using IdxType = Vec8i;
	using BoolType = Vec8fb;
	using FloatType = float;

	static const int width = 8;
	static const float nullValue;
	static const float oneValue;
	static constexpr  bool alignedLoadStore = true;
	static constexpr  bool boolTypeIsAlignedLoadStore = false;
	static constexpr  bool useScatter = true;

	static constexpr  bool isCompact = false;
	using RegBoolType = Vec8fb;
	using MemBoolType = Vec8f;

};



template<>
struct InstructionTraits<Vec8d>
{
	using IdxType = Vec8i;

	using BoolType = Vec8db;

	using FloatType = double;
	static const int width = 8;
	static const double nullValue;
	static const double oneValue;
	static constexpr  bool alignedLoadStore = false;
	static constexpr  bool boolTypeIsAlignedLoadStore = false;
	static constexpr  bool useScatter = true;

	static constexpr  bool isCompact =  true;// false;// true;
	using RegBoolType = Vec8db;// Vec8db;
	using MemBoolType = Vec8d;

	//using BoolType = NOT_KNOWN_HERE;

};



template<>
struct InstructionTraits<Vec16f>
{
	using IdxType = Vec16i;
	using BoolType = Vec16fb;
	using FloatType = float;
	static const int width = 16;
	static const float nullValue;
	static const float oneValue;
	static constexpr  bool alignedLoadStore = false;
	static constexpr  bool boolTypeIsAlignedLoadStore = false;
	static constexpr  bool useScatter = true;

	static constexpr  bool isCompact = true; //to do
	using RegBoolType = Vec16b;
	using MemBoolType = Vec16b;
};




template<typename TRAIT>
inline  typename InstructionTraits<TRAIT>::MemBoolType boolCompactSave(typename InstructionTraits<TRAIT>::RegBoolType regVal )
{
	return  static_cast<typename InstructionTraits<TRAIT>::MemBoolType>(regVal);
}


template<typename TRAIT  >
inline  auto boolConvert(typename  InstructionTraits<TRAIT>::RegBoolType regVal)
{
	if constexpr (! InstructionTraits<TRAIT>::isCompact )
	{
		return regVal;
	}
	else
	{
		return  boolCompactSave< TRAIT>(regVal);
	}
}



static inline Vec8d boolCompactSave(Vec8db regVal)
{
	//auto in = to_bits(regVal);
	Vec8db ret = regVal;
	//ret.load_bits(in);
	Vec8d ones = 11111111111.;
	Vec8d zeros = 0;
	auto result = select(ret, ones, zeros);
	return result;
}


/*
static inline Vec8d boolCompactSave( Vec8b regVal)
{
	auto in = to_bits(regVal);
	Vec8db ret;
	ret.load_bits(in);
	Vec8d ones = 11111111111.;
	Vec8d zeros = 0;
	auto result = select(ret, ones, zeros);
	return result;

}
*/