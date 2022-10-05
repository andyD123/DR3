/****************************  vec.h   *******************************
* Author:        Andrew Drakeford
* Date created:  2021-04-01
* Last modified: 2021-04-01
* Version:       1.0
* Project:       DR Cubed
* Description:
*
* (c) Copyright 2019 Andrew Drakeford
* Apache License version 2.0 or later.
*****************************************************************************/
#pragma once
#include "vec.h"
#include "vec_bool.h"
#include "instruction_traits.h"
#include "binary_unitary_operations.h"
#include "unroll_operators.h"
#include "error_utils.h"

#include <stdexcept>

// split the accumulate and transform


/////////////////////////////////////////////////////////////////
//  we have two cases   example find max   or sum elements
// find max, we can introduce the max value as first element value to all compute lanes
// where as if singular value to be incorporated once  such a summation, null the start values and add at the end.
//
////////////////////////////////////////////////////////////////////
template< typename INS_VEC, typename OP, typename OP_SCALAR>
typename InstructionTraits<INS_VEC>::FloatType ApplyAccumulate2(const Vec<INS_VEC>& rhs1, OP& oper, OP_SCALAR& operS, typename InstructionTraits<INS_VEC>::FloatType initVal, bool singularInit = true)
{
	check_vector(rhs1);

	if (isScalar(rhs1))
	{
		return ApplyBinaryOperation1<INS_VEC, OP_SCALAR>(typename InstructionTraits<INS_VEC>::FloatType(rhs1.getScalarValue()), typename InstructionTraits<INS_VEC>::FloatType(initVal), operS);
	}

	auto pRhs1 = rhs1.start();
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 1 * width;


	INS_VEC RHS1;
	INS_VEC RES = singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue;
	int i = 0;
	int rhsSZ = static_cast<int>(rhs1.size() - step);
	for (; i < rhsSZ; i += step)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RES, RHS1);

	}

	typename InstructionTraits<INS_VEC>::FloatType result = RES[0];
	//across vectors
	for (int j = 1; j < width; ++j)
	{
		result = operS(result, RES[j]);
	}

	//odd bit at end
	for (; i < static_cast<int>(rhs1.size()); ++i)
	{
		result = operS(result, pRhs1[i]);
	}


	if (singularInit)
	{
		return result;
	}

	return result = operS(result, initVal);

}

//unroll version
template< typename INS_VEC, typename OP, typename OP_SCALAR>
typename InstructionTraits<INS_VEC>::FloatType ApplyAccumulate2UR(const Vec<INS_VEC>& rhs1, OP& oper, OP_SCALAR& operS, typename InstructionTraits<INS_VEC>::FloatType initVal)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		return ApplyBinaryOperation1<INS_VEC, OP_SCALAR>(typename InstructionTraits<INS_VEC>::FloatType(rhs1.getScalarValue()), typename InstructionTraits<INS_VEC>::FloatType(initVal), operS);
	}

	auto pRhs1 = rhs1.start();
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;


	INS_VEC RHS1;
	INS_VEC RES;

	INS_VEC RHS2;
	INS_VEC RES1;

	INS_VEC RHS3;
	INS_VEC RES2;

	INS_VEC RHS4;
	INS_VEC RES3;

	//int sz = static_cast<int>( rhs1.paddedSize());

	typename InstructionTraits<INS_VEC>::FloatType res = initVal;
	RES = InstructionTraits<INS_VEC>::nullValue;

	RES1 = RES;
	RES2 = RES;
	RES3 = RES;

	int i = 0;
	int rhsSZ = static_cast<int>(rhs1.size() - step);
	for (; i < rhsSZ; i += step)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RES, RHS1);

		RHS2.load_a(pRhs1 + i + width);
		RES1 = oper(RES1, RHS2);

		RHS3.load_a(pRhs1 + i + width * 2);
		RES2 = oper(RES2, RHS3);

		RHS4.load_a(pRhs1 + i + width * 3);
		RES3 = oper(RES3, RHS4);

	}

	RES = oper(RES, RES1);
	RES2 = oper(RES2, RES3);
	RES = oper(RES, RES2);


	typename InstructionTraits<INS_VEC>::FloatType result = InstructionTraits<INS_VEC>::nullValue;
	//across vectors
	for (int j = 0; j < width; ++j)
	{
		result = operS(result, RES[j]);
	}

	//odd bit at end
	for (; i < static_cast<int>(rhs1.size()); ++i)
	{
		result = operS(result, pRhs1[i]);
	}

	return result = operS(result, res);


}


//accumulate needs a scalar operator  since is a reduction
template< typename INS_VEC, typename OP, typename OPT_SCALAR, typename OPT, typename OP_ACC_SCALAR>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate(const Vec<INS_VEC>& rhs1, OPT& operTransform, OPT_SCALAR& operTransformScalar, OP& operAcc, OP_ACC_SCALAR& operAccScalar, typename InstructionTraits<INS_VEC>::FloatType initVal)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		auto trform = ApplyUnitaryOperation(rhs1, operTransform);
		return ApplyBinaryOperation1<INS_VEC, OP>(trform.getScalarValue(), initVal, operAcc);
	}

	auto pRhs1 = rhs1.start();
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 1 * width;


	INS_VEC RHS1;
	INS_VEC RES;

	//int sz = rhs1.paddedSize();
	int sz = static_cast<int>(rhs1.size());


	RES = InstructionTraits<INS_VEC>::nullValue;
	int i = 0;
	int rhsSZ = sz - step;
	for (; i < rhsSZ; i += step)
	{
		RHS1.load_a(pRhs1 + i);
		RES = operAcc(RES, operTransform(RHS1));
	}

	typename InstructionTraits<INS_VEC>::FloatType result = initVal;
	//sum across reg elements  scalar// could use hadd  ?
	for (int j = 0; j < width; ++j)
	{
		result = ApplyBinaryOperation1<INS_VEC, OP>(RES[j], result, operAcc);
	}


	//odd end bit ?
	for (; i < sz; ++i)
	{
		result = operAccScalar(result, operTransformScalar(pRhs1[i]));

	}


	return result;

}



template< typename INS_VEC, typename OP, typename OPT_SCALAR, typename OPT, typename OP_ACC_SCALAR>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulateUR(const Vec<INS_VEC>& rhs1, OPT& operTransform, OPT_SCALAR& operTransformScalar, OP& operAcc, OP_ACC_SCALAR& operAccScalar, typename InstructionTraits<INS_VEC>::FloatType initVal)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		auto trform = ApplyUnitaryOperation(rhs1, operTransform);
		return ApplyBinaryOperation1<INS_VEC, OP>(trform.getScalarValue(), initVal, operAcc);
	}

	auto pRhs1 = rhs1.start();
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	INS_VEC RHS1;
	INS_VEC RES;

	INS_VEC RHS2;
	INS_VEC RES1;

	INS_VEC RHS3;
	INS_VEC RES2;

	INS_VEC RHS4;
	INS_VEC RES3;

	RES = InstructionTraits<INS_VEC>::nullValue;
	RES1 = RES;
	RES2 = RES;
	RES3 = RES;

	int i = 0;
	int sz = static_cast<int>(rhs1.size());
	int rhsSZ = static_cast<int>(rhs1.size()) - step;
	for (; i < rhsSZ; i += step)
	{

		RHS1.load_a(pRhs1 + i);
		RES = operAcc(RES, operTransform(RHS1));

		RHS2.load_a(pRhs1 + i + width);
		RES1 = operAcc(RES1, operTransform(RHS2));

		RHS3.load_a(pRhs1 + i + width * 2);
		RES2 = operAcc(RES2, operTransform(RHS3));

		RHS4.load_a(pRhs1 + i + width * 3);
		RES3 = operAcc(RES3, operTransform(RHS4));

	}
	//accumulate unrolled 
	RES = operAcc(RES, RES1);
	RES2 = operAcc(RES2, RES3);
	RES = operAcc(RES, RES2);


	typename InstructionTraits<INS_VEC>::FloatType result = initVal;
	for (int j = 0; j < width; ++j)
	{
		result = operAccScalar(RES[j], result);
	}

	//odd end bit
	for (; i < sz; ++i)
	{
		result = operAccScalar(operTransformScalar(pRhs1[i]), result);
	}

	return result;

}



//scalar operations will be implemented using vector types operating on single lane, saves switching instructions sets and giving multiple lambda
//deprecate the above lambda
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyAccumulate(const Vec<INS_VEC>& rhs1, OP& oper, typename InstructionTraits<INS_VEC>::FloatType initVal)
{
	check_vector(rhs1);
	auto scalarCopyOfOp = oper;
	auto scalarCopyOfOp1 = oper;
	auto scalarCopyOfOp2 = oper;


	if (isScalar(rhs1))
	{
		return ApplyBinaryOperationVec<INS_VEC, OP>(rhs1.getScalarValue(), initVal, oper);
	}

	auto pRhs1 = rhs1.start();
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 1 * width;


	INS_VEC RHS1;
	INS_VEC RES;


	RES = InstructionTraits<INS_VEC>::nullValue;

	int i = 0;
	int range = std::max(0, rhs1.size() - step);
	for (; i < range; i += step)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RES, RHS1);
	}
	typename InstructionTraits<INS_VEC>::FloatType result = InstructionTraits<INS_VEC>::nullValue; // initVal;
	for (int j = 0; j < width; ++j)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(RES[j], result, scalarCopyOfOp);
	}

	//end bits for vecs not filling padding
	for (; i < rhs1.size(); ++i)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(pRhs1[i], result, scalarCopyOfOp1);
	}

	result = ApplyBinaryOperationVec<INS_VEC, OP>(result, initVal, scalarCopyOfOp2);
	return result;
}



template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyAccumulate2(const Vec<INS_VEC>& rhs1, OP& oper, typename InstructionTraits<INS_VEC>::FloatType initVal, bool singularInit = true)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		return ApplyBinaryOperationVec<INS_VEC, OP>(typename InstructionTraits<INS_VEC>::FloatType(rhs1.getScalarValue()), typename InstructionTraits<INS_VEC>::FloatType(initVal), oper);
	}

	auto pRhs = rhs1.start();
	int sz = rhs1.size();
	return Unroll_Accumulate<INS_VEC, OP>::apply_1(sz, pRhs, oper, initVal, singularInit);
}


//unroll version
template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyAccumulate2UR(const Vec<INS_VEC>& rhs1, OP& oper, typename InstructionTraits<INS_VEC>::FloatType initVal, bool singularInit = true)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		return ApplyBinaryOperationVec<INS_VEC, OP>(typename InstructionTraits<INS_VEC>::FloatType(rhs1.getScalarValue()), typename InstructionTraits<INS_VEC>::FloatType(initVal), oper);
	}

	auto pRhs = rhs1.start();
	int sz = rhs1.size();
	return Unroll_Accumulate<INS_VEC, OP>::apply_4(sz, pRhs, oper, initVal, singularInit);
}

//////////////////////

//unrolled version helps greatly with VC2019
template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyAccumulate2UR_X(const Vec<INS_VEC>& rhs1, OP& oper )
{
	check_vector(rhs1);
	if (isScalar(rhs1)) // nothing to accumulate with so just return  value
	{
		return rhs1.getScalarValue();
	}

	int sz = rhs1.size();
	auto pRhs1 = rhs1.start();
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	INS_VEC RHS1;
	INS_VEC RES;

	INS_VEC RHS2;
	INS_VEC RES1;

	INS_VEC RHS3;
	INS_VEC RES2;

	INS_VEC RHS4;
	INS_VEC RES3;

	int i = 0;

	if (sz >= step * 2)
	{
		//initialise first set of registers
		{
			RES.load_a(pRhs1 + i);
			RES1.load_a(pRhs1 + i + width);
			RES2.load_a(pRhs1 + i + width * 2);
			RES3.load_a(pRhs1 + i + width * 3);
		}

		i += step;
		int rhsSZ = rhs1.size();
		for (; i <= (rhsSZ - step); i += step)
		{
			RHS1.load_a(pRhs1 + i);
			RES = oper(RES, RHS1);

			RHS2.load_a(pRhs1 + i + width);
			RES1 = oper(RES1, RHS2);

			RHS3.load_a(pRhs1 + i + width * 2);
			RES2 = oper(RES2, RHS3);

			RHS4.load_a(pRhs1 + i + width * 3);
			RES3 = oper(RES3, RHS4);

		}

		// odd bits
		for (; i <= rhsSZ-width; i += width)
		{
			RHS1.load_a(pRhs1 + i);
			RES = oper(RES, RHS1);
		}

		RES = oper(RES, RES1);
		RES2 = oper(RES2, RES3);
		RES = oper(RES, RES2);

	}
	else
	{
		RES.load_a(pRhs1);

		i += width; 
		// odd bits
		for (; i <= sz - width; i += width)
		{
			RHS1.load_a(pRhs1 + i);
			RES = oper(RES, RHS1);
		}

	}

	typename InstructionTraits<INS_VEC>::FloatType result = RES[0];
	int min_wdth = std::min(sz, width);
	//across vectors lanes  // not assuming horizontal versoion exist
	for (int j = 1; j < min_wdth; ++j)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(result, RES[j], oper);
	}

	//end bits for vecs not filling padding
	for (; i < rhs1.size(); ++i)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(pRhs1[i], result, oper);
	}

	return result;
}




template< typename INS_VEC,  typename OPT, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate(const Vec<INS_VEC>& rhs1, OPT& operTransform,  OP& operAcc,  typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		auto trform = ApplyUnitaryOperation(rhs1, operTransform);
		return ApplyBinaryOperationVec<INS_VEC, OP>(trform.getScalarValue(),  initVal, operAcc);
	}

	auto pRhs = rhs1.start();
	int sz = static_cast<int>(rhs1.size());
	return  Unroll_TransformAccumulate<INS_VEC, OP, OPT>::apply_1(sz, pRhs, operAcc, operTransform, initVal, singularInit);
}




template< typename INS_VEC, typename OP, typename OPT>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulateUR(const Vec<INS_VEC>& rhs1, OPT& operTransform,  OP& operAcc,  typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		auto trform = ApplyUnitaryOperation(rhs1, operTransform);
		return ApplyBinaryOperation1<INS_VEC, OP>(trform.getScalarValue(), initVal, operAcc);
	}

	auto pRhs = rhs1.start();
	int sz = rhs1.size();
	return  Unroll_TransformAccumulate<INS_VEC, OP, OPT>::apply_4(sz, pRhs, operAcc, operTransform, initVal, singularInit);
}



//////////////////// replacement version //////////////
//unitary 
//unrolled version helps greatly with VC2019
template< typename INS_VEC, typename OPT, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate2UR_X(const Vec<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc)
{
	check_vector(rhs1);
	if (isScalar(rhs1)) // nothing to accumulate with so just transform  and  return  value
	{
		auto trform = ApplyUnitaryOperation(rhs1, operTransform);
		return trform.getScalarValue();

	}

	int sz = rhs1.size();
	auto pRhs1 = rhs1.start();
	const int width = InstructionTraits<INS_VEC>::width;

	auto zero = InstructionTraits<INS_VEC>::nullValue;

	int step = 4 * width;

	INS_VEC RHS1 = zero;
	INS_VEC RES = zero;

	INS_VEC RHS2 = zero;
	INS_VEC RES1 = zero;

	INS_VEC RHS3 = zero;
	INS_VEC RES2 = zero;

	INS_VEC RHS4 = zero;
	INS_VEC RES3 = zero;

	int i = 0;

	if (sz >= step * 2)
	{
		//initialise first set of registers
		{
			RHS1.load_a(pRhs1 + i);
			RHS2.load_a(pRhs1 + i + width);
			RHS3.load_a(pRhs1 + i + width * 2);
			RHS4.load_a(pRhs1 + i + width * 3);

			RES =  operTransform(RHS1);
			RES1 = operTransform(RHS2);
			RES2 = operTransform(RHS3);
			RES3 = operTransform(RHS4);
		}

		i += step;
		int rhsSZ = rhs1.size();
		for (; i <= (rhsSZ - step); i += step)
		{
			RHS1.load_a(pRhs1 + i);
			RHS2.load_a(pRhs1 + i + width);
			RHS3.load_a(pRhs1 + i + width * 2);
			RHS4.load_a(pRhs1 + i + width * 3);

			RES = operAcc(RES, operTransform(RHS1) );
			RES1 = operAcc(RES1, operTransform(RHS2) );
			RES2 = operAcc(RES2, operTransform(RHS3) );
			RES3 = operAcc(RES3, operTransform(RHS4));

		}

		// odd bits
		for (; i <= rhsSZ - width; i += width)
		{
			RHS1.load_a(pRhs1 + i);
			RES = operAcc(RES, operTransform(RHS1));
		}

		RES = operAcc(RES, RES1);
		RES2 = operAcc(RES2, RES3);
		RES = operAcc(RES, RES2);

	}
	else
	{
		RHS1.load_a(pRhs1);
		RES = operTransform(RHS1);

		i += width;
		// odd bits
		for (; i <= sz - width; i += width)
		{
			RHS1.load_a(pRhs1 + i);
			RES = operAcc(RES, operTransform(RHS1) );
		}

	}

	typename InstructionTraits<INS_VEC>::FloatType result = RES[0];
	int min_wdth = std::min(sz, width);
	//across vectors lanes  // not assuming horizontal versoion exist
	for (int j = 1; j < min_wdth; ++j)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(result, RES[j], operAcc);
	}

	//end bits for vecs not filling padding
	for (; i < rhs1.size(); ++i)
	{
		//need to transform
		typename InstructionTraits<INS_VEC>::FloatType trfmResult = operTransform(INS_VEC(pRhs1[i]))[0];
		result = ApplyBinaryOperationVec<INS_VEC, OP>(result, trfmResult, operAcc);
	}

	return result;
}




/////////////////////////////////////////
template< typename INS_VEC, typename OP, typename OPT>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulateUR(const Vec<INS_VEC>& lhs1, const Vec<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{
	check_pair(lhs1, rhs1);
	//to do cover case of scalars
	if (isScalar(rhs1) && isScalar(lhs1))
	{
		return ApplyBinaryOperation1<INS_VEC, OP>(operTransform(lhs1.getScalarValue(), rhs1.getScalarValue()), initVal, operAcc);
	}

	auto pLhs = lhs1.start();
	auto pRhs = rhs1.start();
	int sz = rhs1.size();
	return  Unroll_TransformAccumulate<INS_VEC, OP, OPT>::apply_4(sz,pLhs, pRhs, operAcc, operTransform, initVal, singularInit);
}

///////////////////////////////
template< typename INS_VEC, typename OPT, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate2UR_X(const Vec<INS_VEC>& lhs1, const Vec<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc)
{

	check_pair(lhs1, rhs1);
	auto zero = InstructionTraits<INS_VEC>::nullValue;
	//to do cover case of scalars
	if (isScalar(rhs1) && isScalar(lhs1))
	{
		return ApplyBinaryOperation1<INS_VEC, OP>(operTransform(lhs1.getScalarValue(), rhs1.getScalarValue()) , zero, operAcc);
	}


	auto pRhs1 = rhs1.start();
	auto pLhs1 = lhs1.start();
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	int sz = rhs1.size();
	Vec<INS_VEC> ret(sz);

	INS_VEC LHS1 = zero;
	INS_VEC RHS1 = zero;
	INS_VEC RES = zero;

	INS_VEC LHS2 = zero;
	INS_VEC RHS2 = zero;
	INS_VEC RES1 = zero;

	INS_VEC LHS3 = zero;
	INS_VEC RHS3 = zero;
	INS_VEC RES2 = zero;

	INS_VEC LHS4 = zero;
	INS_VEC RHS4 = zero;
	INS_VEC RES3 = zero;

	int i = 0;

	if (sz >= step * 2)
	{
		//initialise first set of registers
		{
			RHS1.load_a(pRhs1 + i);
			RHS2.load_a(pRhs1 + i + width);
			RHS3.load_a(pRhs1 + i + width * 2);
			RHS4.load_a(pRhs1 + i + width * 3);

			LHS1.load_a(pLhs1 + i);
			LHS2.load_a(pLhs1 + i + width);
			LHS3.load_a(pLhs1 + i + width * 2);
			LHS4.load_a(pLhs1 + i + width * 3);

			RES = operTransform(LHS1,RHS1);
			RES1 = operTransform(LHS2,RHS2);
			RES2 = operTransform(LHS3,RHS3);
			RES3 = operTransform(LHS4,RHS4);
		}

		i += step;
		//int rhsSZ = rhs1.size();
		int impSZ = rhs1.paddedSize();
		//int rhsSZ = sz - step;
		int rhsSZ = impSZ - step;

		for (; i <= (rhsSZ - step); i += step)
		{
			LHS1.load_a(pLhs1 + i);
			RHS1.load_a(pRhs1 + i);
			RES = operAcc(RES, operTransform(LHS1,RHS1));

			LHS2.load_a(pLhs1 + i + width);
			RHS2.load_a(pRhs1 + i + width);
			RES1 = operAcc(RES1, operTransform(LHS2,RHS2));

			LHS3.load_a(pLhs1 + i + width * 2); 
			RHS3.load_a(pRhs1 + i + width * 2);
			RES2 = operAcc(RES2, operTransform(LHS3,RHS3));

			LHS4.load_a(pLhs1 + i + width * 3);
			RHS4.load_a(pRhs1 + i + width * 3);
			RES3 = operAcc(RES3, operTransform(LHS4,RHS4));

		}

		// odd bits
		for (; i <= rhsSZ - width; i += width)
		{
			LHS1.load_a(pLhs1 + i);
			RHS1.load_a(pRhs1 + i);
			RES = operAcc(RES, operTransform(LHS1,RHS1));
		}

		RES = operAcc(RES, RES1);
		RES2 = operAcc(RES2, RES3);
		RES = operAcc(RES, RES2);

	}
	else
	{
		LHS1.load_a(pLhs1 );
		RHS1.load_a(pRhs1);
		RES = operTransform(LHS1,RHS1);

		i += width;
		// odd bits
		for (; i <= sz - width; i += width)
		{
			LHS1.load_a(pLhs1 + i);
			RHS1.load_a(pRhs1 + i);
			RES = operAcc(RES, operTransform(LHS1,RHS1));
		}

	}

	typename InstructionTraits<INS_VEC>::FloatType result = RES[0];
	int min_wdth = std::min(sz, width);
	//across vectors lanes  // not assuming horizontal versoion exist
	for (int j = 1; j < min_wdth; ++j)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(result, RES[j], operAcc);
	}

	//end bits for vecs not filling padding
	for (; i < rhs1.size(); ++i)
	{
		//need to transform before aggregate
		typename InstructionTraits<INS_VEC>::FloatType trfmResult = operTransform(INS_VEC(pLhs1[i]), INS_VEC(pRhs1[i]))[0];
		result = ApplyBinaryOperationVec<INS_VEC, OP>(result, trfmResult, operAcc);
	}
	return result;

}




///////////////////////////////



template< typename INS_VEC, typename OP >
Vec<INS_VEC>  ApplyTransformUR_X(const Vec<INS_VEC>& rhs1, OP& oper)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		return ApplyUnitaryOperation<INS_VEC, OP>(rhs1.getScalarValue(), oper);
	}

	int sz = rhs1.size();
	auto pRhs1 = rhs1.start();
	Vec<INS_VEC> ret(sz);// ret(rhs1); //???
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	INS_VEC RHS1;
	INS_VEC RES;

	INS_VEC RHS2;
	INS_VEC RES1;

	INS_VEC RHS3;
	INS_VEC RES2;

	INS_VEC RHS4;
	INS_VEC RES3;

	int i = 0;
	
	//int rhsSZ = sz - step;
	int impSZ = rhs1.paddedSize();
	//int rhsSZ = sz - step;
	int rhsSZ = impSZ - step;
	for (; i < rhsSZ ; i += step)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RHS1);
		RES.store_a(pRet + i);

		RHS2.load_a(pRhs1 + i + width);
		RES1 = oper(RHS2);
		RES1.store_a(pRet + i + width);

		RHS3.load_a(pRhs1 + i + width * 2);
		RES2 =  oper(RHS3);
		RES2.store_a(pRet + i + width * 2);

		RHS4.load_a(pRhs1 + i + width * 3);
		RES3 =  oper(RHS4);
		RES3.store_a(pRet + i + width * 3);
	}

	for (; i <= impSZ -width; i += width)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RHS1);
		RES.store_a(pRet + i);
	}

	//to do odd end bit
	//no odd bit
	return ret;
}




template< typename INS_VEC, typename OP >
void ApplyTransformUR_X( VecView<INS_VEC>& rhs1, OP& oper)
{

	if (!rhs1.isScalar())
	{

	
	check_vector(rhs1); //calls overload with a view
	//views are not scalar

	//int sz = rhs1.size();
	auto pRhs1 = rhs1.start();
	auto pRet = pRhs1;

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	INS_VEC RHS1;
	INS_VEC RES;

	INS_VEC RHS2;
	INS_VEC RES1;

	INS_VEC RHS3;
	INS_VEC RES2;

	INS_VEC RHS4;
	INS_VEC RES3;

	int i = 0;

	//int rhsSZ = sz - step;
	int impSZ = rhs1.paddedSize();
	//int impSZ = rhs1.fillSize();
	//int rhsSZ = sz - step;
	int rhsSZ = impSZ - step;

	for (; i < rhsSZ; i += step)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RHS1);
		RES.store_a(pRet + i);

		RHS2.load_a(pRhs1 + i + width);
		RES1 = oper(RHS2);
		RES1.store_a(pRet + i + width);

		RHS3.load_a(pRhs1 + i + width * 2);
		RES2 = oper(RHS3);
		RES2.store_a(pRet + i + width * 2);

		RHS4.load_a(pRhs1 + i + width * 3);
		RES3 = oper(RHS4);
		RES3.store_a(pRet + i + width * 3);
	}

	for (; i <= impSZ - width; i += width)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RHS1);
		RES.store_a(pRet + i);
	}

 //views are padded and filled to width of register 
 // so no end bits

	}
	else
	{
		auto val = rhs1.getScalarValue();
		auto scalarRes = oper(INS_VEC(val))[0];
		VecView<INS_VEC> result;
		result = scalarRes;
		rhs1 = result;
	}
 
	
}




template< typename INS_VEC, typename OP >
VecView<INS_VEC> ApplyTransformUR_X(const VecView<INS_VEC>& rhs1, OP& oper)
{
	check_vector(rhs1); //calls overload with a view
	//views are not scalar



	int sz = rhs1.size();
	auto pRhs1 = rhs1.start();
	VecView<INS_VEC> ret(rhs1);// static_cast<size_t>(sz));
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	INS_VEC RHS1;
	INS_VEC RES;

	INS_VEC RHS2;
	INS_VEC RES1;

	INS_VEC RHS3;
	INS_VEC RES2;

	INS_VEC RHS4;
	INS_VEC RES3;

	int i = 0;

	//int rhsSZ = sz - step;
	int impSZ = rhs1.paddedSize();
	//int impSZ = rhs1.fillSize();
	//int rhsSZ = sz - step;
	int rhsSZ = impSZ - step;

	for (; i < rhsSZ; i += step)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RHS1);
		RES.store_a(pRet + i);

		RHS2.load_a(pRhs1 + i + width);
		RES1 = oper(RHS2);
		RES1.store_a(pRet + i + width);

		RHS3.load_a(pRhs1 + i + width * 2);
		RES2 = oper(RHS3);
		RES2.store_a(pRet + i + width * 2);

		RHS4.load_a(pRhs1 + i + width * 3);
		RES3 = oper(RHS4);
		RES3.store_a(pRet + i + width * 3);
	}

	for (; i <= impSZ - width; i += width)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RHS1);
		RES.store_a(pRet + i);
	}

	//views are padded and filled to width of register 
	// so no end bits
	return ret;

}




template< typename INS_VEC, typename OP >
Vec<INS_VEC>  ApplyTransformUR_XX(const Vec<INS_VEC>& rhs1, OP& oper)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		return ApplyUnitaryOperation<INS_VEC, OP>(rhs1.getScalarValue(), oper);
	}

	int sz = rhs1.size();
	auto pRhs1 = rhs1.start();
	Vec<INS_VEC> ret(sz);
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 8 * width;

	INS_VEC RHS1;
	INS_VEC RES;

	INS_VEC RHS2;
	INS_VEC RES1;

	INS_VEC RHS3;
	INS_VEC RES2;

	INS_VEC RHS4;
	INS_VEC RES3;


	INS_VEC RHS5;
	INS_VEC RES4;

	INS_VEC RHS6;
	INS_VEC RES5;

	INS_VEC RHS7;
	INS_VEC RES6;

	INS_VEC RHS8;
	INS_VEC RES7;

	int i = 0;

	//int rhsSZ = sz - step;
	int impSZ = rhs1.paddedSize();
	//int rhsSZ = sz - step;
	int rhsSZ = impSZ - step;

	for (; i < rhsSZ; i += step)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RHS1);
		RES.store_a(pRet + i);

		RHS2.load_a(pRhs1 + i + width);
		RES1 = oper(RHS2);
		RES1.store_a(pRet + i + width);

		RHS3.load_a(pRhs1 + i + width * 2);
		RES2 = oper(RHS3);
		RES2.store_a(pRet + i + width * 2);

		RHS4.load_a(pRhs1 + i + width * 3);
		RES3 = oper(RHS4);
		RES3.store_a(pRet + i + width * 3);

		RHS5.load_a(pRhs1 + i + width * 4);
		RES4 = oper(RHS5);
		RES4.store_a(pRet + i + width * 4);

		RHS6.load_a(pRhs1 + i + width*5);
		RES5 = oper(RHS6);
		RES5.store_a(pRet + i + width*5);

		RHS7.load_a(pRhs1 + i + width * 6);
		RES6 = oper(RHS7);
		RES6.store_a(pRet + i + width * 6);

		RHS8.load_a(pRhs1 + i + width * 7);
		RES7 = oper(RHS8);
		RES7.store_a(pRet + i + width * 7);
	}

	for (; i <= impSZ -width; i += width)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RHS1);
		RES.store_a(pRet + i);
	}

	return ret;
}

////////////////////

//TO DO BINARY


template< typename INS_VEC, typename OPER >
Vec<INS_VEC>  ApplyBinaryTransformUR_X(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs, OPER& oper)
{
	check_pair(lhs, rhs);
	if (isScalar(lhs))
	{
		return ApplyBinaryTransformUR_X(lhs.getScalarValue(), rhs, oper);
	}
	if (isScalar(rhs))
	{
		return ApplyBinaryTransformUR_X(lhs, rhs.getScalarValue(), oper);
	}

	int sz = rhs.size();
	auto pRhs = rhs.start();
	auto pLhs = lhs.start();
	Vec<INS_VEC> ret(sz);
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	INS_VEC RHS;
	INS_VEC LHS;
	INS_VEC RES;
	
	INS_VEC RHS1;
	INS_VEC LHS1;
	INS_VEC RES1;
	
	INS_VEC RHS2;
	INS_VEC LHS2;
	INS_VEC RES2;

	INS_VEC RHS3;
	INS_VEC LHS3;
	INS_VEC RES3;


	int i = 0;

	//int rhsSZ = sz - step;
	int impSZ = lhs.paddedSize();
	//int rhsSZ = sz - step;
	int rhsSZ = impSZ - step;

	for (; i < rhsSZ; i += step)
	{
		LHS.load_a(pLhs + i);
		RHS.load_a(pRhs + i);
		RES = oper(LHS, RHS);
		RES.store_a(pRet + i);

		LHS1.load_a(pLhs + i+ width);
		RHS1.load_a(pRhs + i+ width);
		RES1 = oper(LHS1, RHS1);
		RES1.store_a(pRet + i+ width);

		LHS2.load_a(pLhs + i + 2* width);
		RHS2.load_a(pRhs + i + 2*width);
		RES2 = oper(LHS2, RHS2);
		RES2.store_a(pRet + i +2*width);

		LHS3.load_a(pLhs + i + 3*width);
		RHS3.load_a(pRhs + i + 3*width);
		RES3 = oper(LHS3, RHS3);
		RES3.store_a(pRet + i +3*width);
	}

	for (; i <= impSZ - width; i += width)
	{
		LHS.load_a(pLhs + i);
		RHS.load_a(pRhs + i);
		RES = oper(LHS, RHS);
		RES.store_a(pRet + i);
	}

	//Since vector is padded no odd end bits, just unused end bits 
	return ret;
}



template< typename INS_VEC, typename OPER >
Vec<INS_VEC>  ApplyBinaryTransformUR_X(typename InstructionTraits<INS_VEC>::FloatType  lhs, const Vec<INS_VEC>& rhs, OPER& oper)
{
	check_vector(rhs);

	INS_VEC RHS;
	INS_VEC LHS;
	INS_VEC RES;
	LHS = lhs;

	if (isScalar(rhs))
	{
		// both are scalar
		RHS = rhs.getScalarValue();
		RES = oper(LHS, RHS);
		return Vec<INS_VEC>(RES[0]);
	}


	int sz = rhs.size();
	auto pRhs = rhs.start();
	Vec<INS_VEC> ret(sz);
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	INS_VEC RHS1;
	INS_VEC RES1;

	INS_VEC RHS2;
	INS_VEC RES2;

	INS_VEC RHS3;
	INS_VEC RES3;

	int i = 0;

	int impSZ = rhs.paddedSize();
	//int rhsSZ = sz - step;
	int rhsSZ = impSZ - step;
	for (; i < rhsSZ; i += step)
	{
		RHS.load_a(pRhs + i);
		RES = oper(LHS, RHS);
		RES.store_a(pRet + i);

		RHS1.load_a(pRhs + i + width);
		RES1 = oper(LHS, RHS1);
		RES1.store_a(pRet + i + width);

		RHS2.load_a(pRhs + i + 2 * width);
		RES2 = oper(LHS, RHS2);
		RES2.store_a(pRet + i +2 * width);

		RHS3.load_a(pRhs + i + 3 * width);
		RES3 = oper(LHS, RHS3);
		RES3.store_a(pRet + i + 3 * width);
	}

	for (; i <= impSZ - width; i += width)
	{
		RHS.load_a(pRhs + i);
		RES = oper(LHS, RHS);
		RES.store_a(pRet + i);
	}

	//Since vector is padded no odd end bits, just unused end bits 
	return ret;
}



template< typename INS_VEC, typename OPER >
Vec<INS_VEC>  ApplyBinaryTransformUR_X(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType  rhs,  OPER& oper)
{
	check_vector(lhs);

	INS_VEC LHS;
	INS_VEC RHS;
	INS_VEC RES;
	RHS = rhs;

	if (isScalar(lhs))
	{
		// both are scalar
		LHS = lhs.getScalarValue();
		RES = oper(LHS, RHS);
		return Vec<INS_VEC>(RES[0]);
	}


	int sz = lhs.size();
	auto pLhs = lhs.start();
	Vec<INS_VEC> ret(sz);
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	INS_VEC LHS1;
	INS_VEC RES1;

	INS_VEC LHS2;
	INS_VEC RES2;

	INS_VEC LHS3;
	INS_VEC RES3;

	int i = 0;

	int impSZ = lhs.paddedSize();
	//int rhsSZ = sz - step;
	int rhsSZ = impSZ - step;
	for (; i < rhsSZ; i += step)
	{
		LHS.load_a(pLhs + i);
		RES = oper(LHS, RHS);
		RES.store_a(pRet + i);

		LHS1.load_a(pLhs + i + width);
		RES1 = oper(LHS1, RHS);
		RES1.store_a(pRet + i +width);

		LHS2.load_a(pLhs + i + 2 * width);
		RES2 = oper(LHS2, RHS);
		RES2.store_a(pRet + i +2 * width);

		LHS3.load_a(pLhs + i + 3 * width);
		RES3 = oper(LHS3, RHS);
		RES3.store_a(pRet + i + 3 * width);
	}

	for (; i <= impSZ - width; i += width)
	{
		LHS.load_a(pLhs + i);
		RES = oper(LHS, RHS);
		RES.store_a(pRet + i);
	}

	//Since vector is padded no odd end bits, just unused end bits 
	return ret;
}




//////////////////////////////


template< typename INS_VEC, typename OP, typename OPER_TRUE, typename OPER_FALSE >
Vec<INS_VEC>  ApplySelectTransformUR_X(const Vec<INS_VEC>& rhs1, OP& cond, OPER_TRUE& trueOper, OPER_FALSE& falseOper  )
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		//apply same operation across whole register with same values in all lanes
		//and return value from only one lane
		INS_VEC RHS = rhs1.getScalarValue(); //set with scalar value
		INS_VEC TRU = trueOper(RHS);
		INS_VEC FLS = falseOper(RHS);
		INS_VEC RES = select(cond(RHS), TRU, FLS);
		return Vec<INS_VEC>(RES[0]);
	}

	int sz = rhs1.size();
	auto pRhs1 = rhs1.start();
	Vec<INS_VEC> ret(sz);
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	INS_VEC RHS;
	INS_VEC RES;
	INS_VEC TRU;
	INS_VEC FLS;

	INS_VEC RHS1;
	INS_VEC RES1;
	INS_VEC TRU1;
	INS_VEC FLS1;

	INS_VEC RHS2;
	INS_VEC RES2;
	INS_VEC TRU2;
	INS_VEC FLS2;

	INS_VEC RHS3;
	INS_VEC RES3;
	INS_VEC TRU3;
	INS_VEC FLS3;

	int i = 0;

	//int rhsSZ = sz - step;
	int impSZ = rhs1.paddedSize();
	//int rhsSZ = sz - step;
	int rhsSZ = impSZ - step;


	for (; i < rhsSZ; i += step)
	{
		RHS.load_a(pRhs1 + i);
		TRU = trueOper(RHS);
		FLS = falseOper(RHS);
		RES = select(cond(RHS), TRU, FLS);
		RES.store_a(pRet + i);

		RHS1.load_a(pRhs1 + i + width);
		TRU1 = trueOper(RHS1);
		FLS1 = falseOper(RHS1);
		RES1 = select(cond(RHS1), TRU1, FLS1);
		RES1.store_a(pRet + i + width);

		RHS2.load_a(pRhs1 + i + width * 2);
		TRU2 = trueOper(RHS2);
		FLS2 = falseOper(RHS2);
		RES2 = select(cond(RHS2), TRU2, FLS2);
		RES2.store_a(pRet + i + width * 2);

		RHS3.load_a(pRhs1 + i + width * 3);
		TRU3 = trueOper(RHS3);
		FLS3 = falseOper(RHS3);
		RES3 = select(cond(RHS3), TRU3, FLS3);
		RES3.store_a(pRet + i + width * 3);
	}

	for (; i <= impSZ - width; i += width)
	{
		RHS.load_a(pRhs1 + i);
		TRU = trueOper(RHS);
		FLS = falseOper(RHS);
		RES = select(cond(RHS), TRU, FLS);
		RES.store_a(pRet + i);
	}

	//Since vector is padded no odd end bits, just unused end bits 
	return ret;
}


template< typename INS_VEC, typename BOOL_OPER, typename TRUE_OPER, typename FALSE_OPER>
Vec<INS_VEC> ApplySelectionOperationFuncUR_X(BOOL_OPER& COND, const Vec<INS_VEC>& testData, TRUE_OPER& trueOper, FALSE_OPER& falseOper)
{
	return ApplySelectTransformUR_X(testData, COND, trueOper, falseOper);

	/*
	check_vector(testData);
	Vec<INS_VEC> result(testData.size());
	auto pRes = result.start();
	auto pX = testData.start();
	int sz = testData.paddedSize();
	Unroll_Select< INS_VEC, BOOL_OPER, TRUE_OPER, FALSE_OPER>::apply_4(sz, pX, trueOper, falseOper, pRes, COND);
	return result;
	*/
}

template< typename INS_VEC, typename OP >
Vec<INS_VEC>  ApplySelectTransformUR_XC(const Vec<INS_VEC>& rhs1, OP& cond, typename InstructionTraits<INS_VEC>::FloatType trueVal, typename InstructionTraits<INS_VEC>::FloatType& falseVal)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		INS_VEC RHS = rhs1.getScalarValue();
		INS_VEC TRU = trueVal;
		INS_VEC FLS = falseVal;
		INS_VEC RES = select(cond(RHS), TRU, FLS);
		return Vec<INS_VEC>(RES[0]);
	}

	int sz = rhs1.size();
	auto pRhs1 = rhs1.start();
	Vec<INS_VEC> ret(sz);
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	INS_VEC RHS;
	INS_VEC RES;
	INS_VEC TRU = trueVal;
	INS_VEC FLS = falseVal;

	INS_VEC RHS1;
	INS_VEC RES1;

	INS_VEC RHS2;
	INS_VEC RES2;

	INS_VEC RHS3;
	INS_VEC RES3;

	int i = 0;

	//int rhsSZ = sz - step;
	int impSZ = rhs1.paddedSize();
	//int rhsSZ = sz - step;
	int rhsSZ = impSZ - step;

	for (; i < rhsSZ; i += step)
	{
		RHS.load_a(pRhs1 + i);
		RES = select(cond(RHS), TRU, FLS);
		RES.store_a(pRet + i);

		RHS1.load_a(pRhs1 + i + width);
		RES1 = select(cond(RHS1), TRU, FLS);
		RES1.store_a(pRet + i + width);

		RHS2.load_a(pRhs1 + i + width * 2);
		RES2 = select(cond(RHS2), TRU, FLS);
		RES2.store_a(pRet + i + width * 2);

		RHS3.load_a(pRhs1 + i + width * 3);
		RES3 = select(cond(RHS3), TRU, FLS);
		RES3.store_a(pRet + i + width * 3);
	}

	for (; i <= impSZ - width; i += width)
	{
		RHS.load_a(pRhs1 + i);
		RES = select(cond(RHS), TRU, FLS);
		RES.store_a(pRet + i);
	}

	//no odd end bit
	return ret;
}


template< typename INS_VEC, typename BOOL_OPER>
Vec<INS_VEC> ApplySelectionOperationCUR_X(BOOL_OPER& COND, const Vec<INS_VEC>& testData, typename InstructionTraits<INS_VEC>::FloatType trueVal, typename InstructionTraits<INS_VEC>::FloatType& falseVal)
{
	return ApplySelectTransformUR_XC(testData, COND, trueVal, falseVal);

}





///////////////////////////

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyAdjacentDiff(const Vec<INS_VEC>& rhs1, OP& oper)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		throw std::range_error("ApplyAdjacentDiff called on scalar");
	}


	Vec<INS_VEC> result(rhs1.size());

	auto pRes = result.start();
	auto pRhs1 = rhs1.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 1 * width;


	INS_VEC RHS1;
	INS_VEC RHS2;
	INS_VEC RES;

	int sz = rhs1.paddedSize();

	int i = 0;
	for (; i < (sz - step); i += step)
	{
		RHS1.load_a(pRhs1 + i);
		RHS2.load(pRhs1 + i + 1); // load un-alligned
		RES = oper(RHS1, RHS2);
		RES.store_a(pRes + i);
	}
	for (int j = i; j < (sz - 1); j++)
	{
		result[j] = ApplyBinaryOperation1<INS_VEC, OP>(rhs1[j], rhs1[j + 1], oper);
	}

	return result;

};

