/****************************  boolean_operations.h   *******************************
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

#include "vec_bool.h"
#include "instruction_traits.h"
#include "unroll_operators.h"
#include "error_utils.h"

/****************************************************************************
*   *OpElement  is a primative  which loads data to register, applies an operator  and stores result back to memory
*  This can be binary, unitary  for numeric or boolean vector types 
*  The *OpElement types instantiate unroll  blocks 
*  for applying unrolled/ "De-Rolled" code blocks for applying the operator over a larger vector

**********************************************************************************/

template< typename INS_VEC>
typename InstructionTraits<INS_VEC>::BoolType loadBool(typename InstructionTraits<INS_VEC>::FloatType* pdata)
{
	INS_VEC data;
	return InstructionTraits<INS_VEC>::BoolType(data.load_a(pdata));
}

template< typename INS_VEC>
void storeBool_a(typename InstructionTraits<INS_VEC>::BoolType& toStore, typename InstructionTraits<INS_VEC>::FloatType* pdata)
{
	INS_VEC data(toStore);
	data.store_a(pdata);
}

//probably same as above 
template< typename INS_VEC>
void storeBool2_a(typename InstructionTraits<INS_VEC>::RegBoolType& toStore, typename InstructionTraits<INS_VEC>::FloatType* pdata)
{

	auto cnvrt =boolConvert<INS_VEC>(toStore);

	INS_VEC(cnvrt).store_a(pdata);

	//perhaps
	//boolConvert<INS_VEC>(toStore).store_a(pdata);

	//INS_VEC data(tostore);
	//data.store_a(pdata);
}


/*
template< typename INS_VEC>
void storeBool(typename InstructionTraits<INS_VEC>::BoolType& tostore, typename InstructionTraits<INS_VEC>::FloatType* pdata)
{
	INS_VEC data(tostore);
	data.store_a(pdata);
}
*/

template<typename INS_VEC, typename OP, int OFFSET>
struct  UnitaryBoolOpElement
{
	typename InstructionTraits<INS_VEC>::BoolType RHS;
	typename InstructionTraits<INS_VEC>::BoolType RES;
	using Float = typename InstructionTraits<INS_VEC>::FloatType;
	const int width = InstructionTraits< INS_VEC>::width;
	const int relativeOffset = InstructionTraits< INS_VEC>::width * OFFSET;

	inline void apply(Float* pRhs, int i, Float* pRes, OP oper)
	{
		RHS = loadBool<INS_VEC>(pRhs  + relativeOffset +i);
		RES = oper.apply(RHS);
		storeBool_a<INS_VEC>(RES, pRes + relativeOffset +i);
	}
};


template< typename INS_VEC, typename OP>
struct Unroll_UnitaryBoolean
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	static inline void apply_4(int sz, Float* pLhs, Float* pOut, OP oper)
	{
		UnrollUnitaryBody_4< UnitaryBoolOpElement, INS_VEC, OP> unitaryUnrolled;
		unitaryUnrolled.apply(sz, pLhs, pOut, oper);
	}
};


template<typename INS_VEC, typename OP, int OFFSET>
struct  BinaryBoolOpElement
{
	typename InstructionTraits<INS_VEC>::BoolType LHS;
	typename InstructionTraits<INS_VEC>::BoolType RHS;
	typename InstructionTraits<INS_VEC>::BoolType RES;
	using Float = typename InstructionTraits<INS_VEC>::FloatType;
	const int width = InstructionTraits< INS_VEC>::width;
	const int relativeOffset = InstructionTraits< INS_VEC>::width * OFFSET;

	inline void apply(Float* pLhs, Float* pRhs, int i, Float* pRes, OP oper)
	{
		LHS = loadBool<INS_VEC>(pLhs + relativeOffset +i);
		RHS = loadBool<INS_VEC>(pRhs + relativeOffset +i);
		RES = oper.apply(LHS, RHS);
		storeBool_a<INS_VEC>(RES, pRes + relativeOffset +i);
	}

	inline void apply(typename InstructionTraits<INS_VEC>::BoolType LHS , Float* pRhs, int i ,Float* pRes, OP oper)
	{
		RHS = loadBool<INS_VEC>(pRhs + relativeOffset +i);
		RES = oper.apply(LHS, RHS);
		storeBool_a<INS_VEC>(RES, pRes + relativeOffset +i);
	}

	inline void apply(Float* pLhs, typename InstructionTraits<INS_VEC>::BoolType RHS, int i, Float* pRes, OP oper)
	{
		LHS = loadBool<INS_VEC>(pLhs + relativeOffset +i);
		RES = oper.apply(LHS, RHS);
		storeBool_a<INS_VEC>(RES, pRes + relativeOffset +i);
	}

};


template<typename INS_VEC, typename OP, int OFFSET>
struct  BinaryBoolNumericOpElement
{
	INS_VEC LHS;
	INS_VEC RHS;
	typename InstructionTraits<INS_VEC>::BoolType RES;
	using Float = typename InstructionTraits<INS_VEC>::FloatType;
	const int width = InstructionTraits< INS_VEC>::width;
	const int relativeOffset = InstructionTraits< INS_VEC>::width * OFFSET;

	inline void apply(Float* pLhs, Float* pRhs, int  i, Float* pRes, OP oper)
	{
		LHS.load_a(pLhs + relativeOffset +i);
		RHS.load_a(pRhs + relativeOffset +i);
		RES = oper.apply(LHS, RHS);
		//INS_VEC(RES).store_a(pRes + relativeOffset +i);
		storeBool2_a<INS_VEC>(RES, pRes + relativeOffset + i);

		//storeBool2_a(RES, pRes + relativeOffset + i);

		
	}

	inline void apply(typename InstructionTraits<INS_VEC>::BoolType LHS, Float* pRhs, int i, Float* pRes, OP oper)
	{
		RHS.load_a(pRhs + relativeOffset + i);
		RES = oper.apply(LHS, RHS);
		//INS_VEC(RES).store_a(pRes + relativeOffset +i);
		storeBool2_a<INS_VEC>(RES, pRes + relativeOffset + i);
	}

	inline void apply(Float* pLhs, typename InstructionTraits<INS_VEC>::BoolType RHS, int i, Float* pRes,  OP oper)
	{
		LHS.load_a(pLhs + relativeOffset + i);
		RES = oper.apply(LHS, RHS);
		//INS_VEC(RES).store_a(pRes + relativeOffset + i);
		storeBool2_a<INS_VEC>(RES, pRes + relativeOffset + i);
	}

};


template<typename INS_VEC, typename OP>
struct BinaryBoolUnroll
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	static inline void apply_4(int sz, typename InstructionTraits<INS_VEC>::BoolType LHS, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< BinaryBoolOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, LHS, pRhs, pOut, oper);
	}

	static inline void apply_4(int sz, Float* pLhs, typename InstructionTraits<INS_VEC>::BoolType RHS, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< BinaryBoolOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, pLhs, RHS, pOut, oper);
	}

	static inline void apply_4(int sz, Float* pLhs, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< BinaryBoolOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, pLhs, pRhs, pOut, oper);
	}
};

template<typename INS_VEC, typename OP>
struct BinaryBoolNumericUnroll
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	static inline void apply_4(int sz, typename InstructionTraits<INS_VEC>::BoolType LHS, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< BinaryBoolNumericOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, LHS, pRhs, pOut, oper);
	}

	static inline void apply_4(int sz, Float* pLhs, typename InstructionTraits<INS_VEC>::BoolType RHS, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< BinaryBoolNumericOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, pLhs, RHS, pOut, oper);
	}

	static inline void apply_4(int sz, Float* pLhs, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< BinaryBoolNumericOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, pLhs, pRhs, pOut, oper);
	}
};

//for operators such as  !   not
template< typename INS_VEC, typename OP>
VecBool<INS_VEC> ApplyBooleanUnitaryOperation(const VecBool<INS_VEC>& rhs)
{
	check_vector(rhs);
	VecBool<INS_VEC> result(static_cast<int>(rhs.size()));
	auto pOut = result.start();
	auto pLhs = rhs.start();
	int sz = static_cast<int>(rhs.paddedSize());
	OP oper;
	Unroll_UnitaryBoolean<INS_VEC, OP>::apply_4(sz, pLhs, pOut, oper);
	return result;
}

//for binary boolean logical operators &&, OR  
template< typename INS_VEC, typename OP>
VecBool<INS_VEC> ApplyBooleanBinaryOperation(const VecBool<INS_VEC>& lhs, const VecBool<INS_VEC>& rhs)
{
	check_pair(lhs, rhs);
	VecBool<INS_VEC> result(rhs.size());
	auto pRes = result.start();
	auto pLhs = lhs.start();
	auto pRhs = rhs.start();
	int sz = lhs.paddedSize();
	OP oper;
	BinaryBoolUnroll<INS_VEC, OP>::apply_4( sz, pLhs, pRhs, pRes, oper);
	return result;
}

// for binary numeric conditional operators   eg A < b
//aligned load
template< typename INS_VEC, typename OP>
VecBool<INS_VEC> ApplyBooleanBinaryOperation(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	check_pair(lhs, rhs);
	VecBool<INS_VEC> result(static_cast<int>(lhs.size()));
	auto pRes = result.start();
	auto pLhs = lhs.start();
	auto pRhs = rhs.start();
	int sz = lhs.paddedSize();
	OP oper;
	BinaryBoolNumericUnroll<INS_VEC, OP>::apply_4(sz, pLhs, pRhs, pRes, oper);
	return result;
}

// for binary numeric conditional operators   eg A < b  overload 
//aligned load
template< typename INS_VEC, typename OP>
VecBool<INS_VEC> ApplyBooleanBinaryOperation(typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	check_vector(rhs);
	VecBool<INS_VEC> result(rhs.size());
	OP oper;
	auto pRes = result.start();
	auto pRhs = rhs.start();
	int sz = rhs.paddedSize();
	INS_VEC LHS(lhs);
	BinaryBoolNumericUnroll<INS_VEC, OP>::apply_4(sz, LHS, pRhs, pRes, oper);
	return result;
}

// for binary numeric conditional operators   eg A < b  overload 
template< typename INS_VEC, typename OP>
VecBool<INS_VEC> ApplyBooleanBinaryOperation(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	check_vector(lhs);

	VecBool< INS_VEC> result(static_cast<int>(lhs.size()) );

	auto pRes = result.start();
	auto pLhs = lhs.start();
	const int width = InstructionTraits<INS_VEC>::width; 
	int step = 4 * width;
	INS_VEC RHS(rhs);

	INS_VEC LHS;
	typename InstructionTraits<INS_VEC>::BoolType RES;

	INS_VEC LHS1;
	typename InstructionTraits<INS_VEC>::BoolType RES1;

	INS_VEC LHS2;
	typename InstructionTraits<INS_VEC>::BoolType RES2;

	INS_VEC LHS3;
	typename InstructionTraits<INS_VEC>::BoolType RES3;

	OP oper;

	int sz = lhs.paddedSize();


	if constexpr (InstructionTraits<INS_VEC>::alignedLoadStore)
	{
		if constexpr (InstructionTraits<INS_VEC>::boolTypeIsAlignedLoadStore)
		{
			int i = 0;
			if (sz > (4 * width))
			{
				for (; i < sz - (4 * width); i += step)
				{
					LHS.load_a(pLhs + i);
					RES = oper.apply(LHS, RHS);
					RES.store_a(pRes + i);

					LHS1.load_a(pLhs + i + width);
					RES1 = oper.apply(LHS1, RHS);
					RES1.store_a(pRes + i + width);

					LHS2.load_a(pLhs + i + 2 * width);
					RES2 = oper.apply(LHS2, RHS);
					RES2.store_a(pRes + i + 2 * width);

					LHS3.load_a(pLhs + i + 3 * width);
					RES3 = oper.apply(LHS3, RHS);
					RES3.store_a(pRes + i + 3 * width);
				}
			}
			for (; i <= sz - width; i += width)
			{
				LHS.load_a(pLhs + i);
				RES = oper.apply(LHS, RHS);
				RES.store_a(pRes + i);
			}

		}
		else
		{
			int i = 0;
			if (sz > (4 * width))
			{
				for (; i < sz - (4 * width); i += step)
				{
					LHS.load_a(pLhs + i);
					RES = oper.apply(LHS, RHS);
					storeBool_a<INS_VEC>(RES, pRes + i);

					LHS1.load_a(pLhs + i + width);
					RES1 = oper.apply(LHS1, RHS);
					storeBool_a<INS_VEC>(RES1, pRes + i + width);

					LHS2.load_a(pLhs + i + 2 * width);
					RES2 = oper.apply(LHS2, RHS);
					storeBool_a<INS_VEC>(RES2, pRes + i + 2 * width);

					LHS3.load_a(pLhs + i + 3 * width);
					RES3 = oper.apply(LHS3, RHS);
					storeBool_a<INS_VEC>(RES3, pRes + i + 3 * width);
				}
			}
			for (; i <= sz - width; i += width)
			{
				LHS.load_a(pLhs + i);
				RES = oper.apply(LHS, RHS);
				storeBool_a<INS_VEC>(RES, pRes + i);
			}

		}
	}
	else
	{
		int i = 0;
		if (sz > (4 * width))
		{
			for (; i < sz - (4 * width); i += step)
			{
				LHS.load(pLhs + i);
				RES = oper.apply(LHS, RHS);
				//INS_VEC(RES).store(pRes + i);
				storeBool2_a<INS_VEC>(RES, pRes + i);

				LHS1.load(pLhs + i + width);
				RES1 = oper.apply(LHS1, RHS);
				//INS_VEC(boolConvert<INS_VEC> (RES1)  ).store(pRes + i + width);
				storeBool2_a<INS_VEC>(RES1, pRes + i + width);

				LHS2.load(pLhs + i + 2 * width);
				RES2 = oper.apply(LHS2, RHS);
				//INS_VEC(RES2).store(pRes + i + 2 * width);
				storeBool2_a<INS_VEC>(RES2, pRes + i + 2* width);

				LHS3.load(pLhs + i + 3 * width);
				RES3 = oper.apply(LHS3, RHS);
				//INS_VEC(RES3).store(pRes + i + 3 * width);
				storeBool2_a<INS_VEC>(RES3, pRes + i + 3 * width);
			}
		}
		for (; i <= sz - width; i += width)
		{
			LHS.load(pLhs + i);
			RES = oper.apply(LHS, RHS);
			//INS_VEC(RES).store(pRes + i);
			storeBool2_a<INS_VEC>(RES, pRes + i);
		}
	}

	return result;

}



