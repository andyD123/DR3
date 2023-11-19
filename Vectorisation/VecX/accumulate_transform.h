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
#include "scan.h"
#include "span.h"
#include "sampler.h"

#include <stdexcept>
#include <tuple>


/*
when we init a vec for transform we  just allocate a  new one using size
when we init a view for transform. if this is a transform of the existrig
view we need to copy the index so we do a copy construct.

*/
template<typename INS_VEC>
Vec<INS_VEC> init(const Vec<INS_VEC>& rhs)
{
	return Vec<INS_VEC>(rhs.size());
}

template<typename INS_VEC>
Vec<INS_VEC> initTransformer(const Vec<INS_VEC>& rhs)
{
	return Vec<INS_VEC>(rhs.size());
}

template<typename INS_VEC>
VecView<INS_VEC> init(const VecView<INS_VEC>& rhs)
{
	return VecView<INS_VEC>(rhs.size());
}

//we need to copy the indexing
template<typename INS_VEC>
VecView<INS_VEC> initTransformer(const VecView<INS_VEC>& rhs)
{
	return VecView<INS_VEC>( rhs);
}



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


//unrolled version helps greatly with VC2019
template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyAccumulate2UR_X(const VEC_TYPE<INS_VEC>& rhs1, OP& oper, int  ZeroInitTag = 1)
{
	check_vector(rhs1);
	if (isScalar(rhs1)) // nothing to accumulate with so just return  value
	{
		return rhs1.getScalarValue();
	}

	long sz = static_cast<long>(rhs1.size());
	auto pRhs1 = rhs1.start();
	const long width = InstructionTraits<INS_VEC>::width;
	long step = 4 * width;

	INS_VEC RHS1;
	INS_VEC RES;

	INS_VEC RHS2;
	INS_VEC RES1;

	INS_VEC RHS3;
	INS_VEC RES2;

	INS_VEC RHS4;
	INS_VEC RES3;

	long i = 0;

	if (sz >= step * 2)
	{
	
		if  (ZeroInitTag==0)
		{
			//initialise first set of registers
			{
				RES = InstructionTraits<INS_VEC>::nullValue;
				RES1 = InstructionTraits<INS_VEC>::nullValue;
				RES2 = InstructionTraits<INS_VEC>::nullValue;
				RES3 = InstructionTraits<INS_VEC>::nullValue;
			}
		}
		else
		{
			//initialise first set of registers
			{
				RES.load_a(pRhs1 + i);
				RES1.load_a(pRhs1 + i + width);
				RES2.load_a(pRhs1 + i + width * 2);
				RES3.load_a(pRhs1 + i + width * 3);
			}
			i += step;
		}

		

		long rhsSZ = static_cast<long>(rhs1.size());
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
	long min_wdth = std::min(sz, width);
	//across vectors lanes  // not assuming horizontal versoion exist
	for (long j = 1; j < min_wdth; ++j)
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



//unrolled version helps greatly with VC2019
template<  template <class> typename VEC_TYPE, typename INS_VEC, template <class> typename ACCUMULATOR_TYPE, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyAccumulate2UR_X_Accum(const ACCUMULATOR_TYPE< INS_VEC>& acc, const VEC_TYPE<INS_VEC>& rhs1, OP& oper, int  ZeroInitTag = 1)
{
	check_vector(rhs1);
	if (isScalar(rhs1)) // nothing to accumulate with so just return  value
	{
		return rhs1.getScalarValue();
	}

	long sz = static_cast<long>(rhs1.size());
	auto pRhs1 = rhs1.start();
	const long width = InstructionTraits<INS_VEC>::width;
	long step = 4 * width;

	using ACC = ACCUMULATOR_TYPE< INS_VEC>;
	
	INS_VEC RHS1;
	ACC RES = InstructionTraits<INS_VEC>::nullValue;

	INS_VEC RHS2;
	ACC RES1;

	INS_VEC RHS3;
	ACC RES2;

	INS_VEC RHS4;
	ACC RES3;

	long i = 0;

	if (sz >= step * 2)
	{

		if (ZeroInitTag == 0)
		{
			//initialise first set of registers
			{
				RES = InstructionTraits<INS_VEC>::nullValue;
				RES1 = InstructionTraits<INS_VEC>::nullValue;
				RES2 = InstructionTraits<INS_VEC>::nullValue;
				RES3 = InstructionTraits<INS_VEC>::nullValue;
			}
		}
		else
		{
			//initialise first set of registers
			{
				INS_VEC temp;
				temp.load_a(pRhs1 + i);
				RES =temp;
				temp.load_a(pRhs1 + i + width);
				RES1 =temp;
				temp.load_a(pRhs1 + i + width * 2);
				RES2 =temp;
				temp.load_a(pRhs1 + i + width * 3);
				RES3 =temp;
			}
			i += step;
		}



		long rhsSZ = static_cast<long>(rhs1.size());
		for (; i <= (rhsSZ - step); i += step)
		{
			RHS1.load_a(pRhs1 + i);
			oper(RES, RHS1);

			RHS2.load_a(pRhs1 + i + width);
			oper(RES1, RHS2);

			RHS3.load_a(pRhs1 + i + width * 2);
			oper(RES2, RHS3);

			RHS4.load_a(pRhs1 + i + width * 3);
			oper(RES3, RHS4);

		}

		// odd bits
		for (; i <= rhsSZ - width; i += width)
		{
			RHS1.load_a(pRhs1 + i);
			RES = oper(RES, RHS1);
		}

		//RES = oper(RES, RES1);
		//RES2 = oper(RES2, RES3);
		//RES = oper(RES, RES2);

		RES += RES1;
		RES2 += RES3;
		RES += RES2;
	}
	else
	{
		//TO DO

		
		ACC RES= InstructionTraits<INS_VEC>::nullValue;
		RHS1.load_a(pRhs1);
		RES = oper(RES, RHS1);
		i += width;
		// odd bits
		for (; i <= sz - width; i += width)
		{
			RHS1.load_a(pRhs1 + i);
			RES = oper(RES, RHS1);
		}

	}

	//TO DO need masked aggregation

	//should used masked operator id below reg width.
	//also consideration of remainder bits
	return RES.hsum();

	//horizontal add across lanes

	// then some masked add across lanes


	/*
	typename InstructionTraits<INS_VEC>::FloatType result = RES[0];
	long min_wdth = std::min(sz, width);
	//across vectors lanes  // not assuming horizontal versoion exist
	for (long j = 1; j < min_wdth; ++j)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(result, RES[j], oper);
	}

	//end bits for vecs not filling padding
	for (; i < rhs1.size(); ++i)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(pRhs1[i], result, oper);
	}
	*/

	//ACC

	//return result;
}











//experimental pairwise unrolled version  avx512
template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyAccumulate2UR_X_pairwise(const VEC_TYPE<INS_VEC>& rhs1, OP& oper)
{
	check_vector(rhs1);
	if (isScalar(rhs1)) // nothing to accumulate with so just return  value
	{
		return rhs1.getScalarValue();
	}

	long sz = static_cast<long>(rhs1.size());
	auto pRhs = rhs1.start();
	const long width = InstructionTraits<INS_VEC>::width;

	using Float = typename InstructionTraits<INS_VEC>::FloatType;

	
	INS_VEC RHS0;
	INS_VEC RES0;

	INS_VEC RHS1;
	INS_VEC RES1;

	INS_VEC RHS2;
	INS_VEC RES2;

	INS_VEC RHS3;


	auto accum_8 = [&](Float* pRhs1)
	{

		RHS0.load_a(pRhs1);
		RHS1.load_a(pRhs1 + width);
		RHS2.load_a(pRhs1 + width * 2);
		RHS3.load_a(pRhs1 + width * 3);

		RES0 = oper(RHS0, RHS1);
		RES1 = oper(RHS2, RHS3);

		INS_VEC RES = oper(RES0, RES1);

		RHS0.load_a(pRhs1 + width * 4);
		RHS1.load_a(pRhs1 + width *5);
		RHS2.load_a(pRhs1 + width * 6);
		RHS3.load_a(pRhs1 + width * 7);

		RES0 = oper(RHS0, RHS1);
		RES1 = oper(RHS2, RHS3);

		return oper(RES , oper(RES0,RES1));
	
	};



	auto accum_16 = [&](Float* pRhs1)
	{

		RHS0.load_a(pRhs1);
		RHS1.load_a(pRhs1 + width);
		RHS2.load_a(pRhs1 + width * 2);
		RHS3.load_a(pRhs1 + width * 3);

		RES0 = oper(RHS0, RHS1);
		RES1 = oper(RHS2, RHS3);
		INS_VEC RES = oper(RES0, RES1);

		RHS0.load_a(pRhs1 + width * 4);
		RHS1.load_a(pRhs1 + width * 5);
		RHS2.load_a(pRhs1 + width * 6);
		RHS3.load_a(pRhs1 + width * 7);

		RES0 = oper(RHS0, RHS1);
		RES1 = oper(RHS2, RHS3);
		RES =  oper(RES, oper(RES0, RES1));


		RHS0.load_a(pRhs1 + width * 8);
		RHS1.load_a(pRhs1 + width * 9);
		RHS2.load_a(pRhs1 + width * 10);
		RHS3.load_a(pRhs1 + width * 11);

		RES0 = oper(RHS0, RHS1);
		RES1 = oper(RHS2, RHS3);
		INS_VEC RES_T =  oper(RES0, RES1);

		RHS0.load_a(pRhs1 + width * 12);
		RHS1.load_a(pRhs1 + width * 13);
		RHS2.load_a(pRhs1 + width * 14);
		RHS3.load_a(pRhs1 + width * 15);

		RES0 = oper(RHS0, RHS1);
		RES1 = oper(RHS2, RHS3);
		RES_T = oper(RES_T, oper(RES0, RES1));

		return oper(RES, RES_T);
	};


	
	
	auto accum_4 = [&](Float* pRhs1)
	{
		
		RHS0.load_a(pRhs1 );
		RHS1.load_a(pRhs1 + width);
		RHS2.load_a(pRhs1 + width * 2);
		RHS3.load_a(pRhs1 + width * 3);

		RES0 = oper(RHS0, RHS1);
		RES1 = oper(RHS2, RHS3);

		RES2 = oper(RES0, RES1);
		return RES2;

	};


	auto accum_2 = [&](Float* pRhs1)
	{
		RHS0.load_a(pRhs1);
		RHS1.load_a(pRhs1 + width);
		
		RES0 = oper(RHS0, RHS1);	
		return RES0;
	};

    
	INS_VEC ZERO = 0.0;
	INS_VEC val = 0.0;
	INS_VEC RES = 0.0;

	size_t working_Size = sz;

	size_t RemainSZ = 31;
	int remainder = static_cast<int>(working_Size & RemainSZ);

	while (remainder >= InstructionTraits< INS_VEC>::width)
	{ 	
		val.load(pRhs);
		RES = oper(RES, val);

		remainder -= InstructionTraits< INS_VEC>::width;
		pRhs+= InstructionTraits< INS_VEC>::width;
	}

	auto whole_reg_sum = scanN(RES, ZERO, oper)[InstructionTraits< INS_VEC>::width - 1];

	
	//RES += val;
	if (remainder > 0)
	{	
		val = 0.0;
		val.load_partial(remainder, pRhs);
		RES = scanN(val, ZERO, oper);
		whole_reg_sum += RES[remainder];
		pRhs += remainder;
	}


	
	auto accum_recurse = [&](Float* pRhs1, size_t extent, auto&& self)
	{
		if (extent == 16 * width)
		{
			return  accum_16(pRhs1);
		}

		if (extent == 8 * width)
		{
			return  accum_8(pRhs1);
		}

		if (extent == 4 * width)
		{
			return  accum_4(pRhs1);
		}

		if (extent == 2 * width)
		{
			return  accum_2(pRhs1);
		}

		// split  maybe use << for division and bitwise  for size test
		INS_VEC LHS2 = self(pRhs1, extent / 2, self);
		INS_VEC RHS2 = self(pRhs1 + extent / 2, extent / 2, self);
		return oper(LHS2, RHS2);

	};


	//working_Size
	/*
	* keep doubling block size untill reach input extent and accumulate if is power of 2 
	* cf integer power  in math
	*/
	INS_VEC  accumulations_of_all_blocks = 0.0;
	size_t current_block_SZ = RemainSZ + 1;
	while (working_Size >= current_block_SZ)
	{

		if (working_Size & current_block_SZ)
		{
			accumulations_of_all_blocks += accum_recurse(pRhs, current_block_SZ, accum_recurse);
			pRhs += current_block_SZ;
		}
		//pRhs += current_block_SZ;
		current_block_SZ *= 2;
	}

	auto accum_over_all_blocks = scanN(accumulations_of_all_blocks, ZERO, oper)[InstructionTraits< INS_VEC>::width - 1];

	return accum_over_all_blocks + whole_reg_sum;

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











template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OPT, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate2UR_X_Impl(const VEC_TYPE<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc)
{
	check_vector(rhs1);
	if (isScalar(rhs1)) // nothing to accumulate with so just transform  and  return  value
	{
		auto trform = ApplyUnitaryOperation(rhs1, operTransform);
		return trform.getScalarValue();
	}

	long sz = static_cast<long>(rhs1.size());
	auto pRhs1 = rhs1.start();
	const long width = InstructionTraits<INS_VEC>::width;

	auto zero = InstructionTraits<INS_VEC>::nullValue;

	long step = 4 * width;

	INS_VEC RHS1 = zero;
	INS_VEC RES = zero;

	INS_VEC RHS2 = zero;
	INS_VEC RES1 = zero;

	INS_VEC RHS3 = zero;
	INS_VEC RES2 = zero;

	INS_VEC RHS4 = zero;
	INS_VEC RES3 = zero;

	long i = 0;

	if (sz >= step * 2)
	{
		//initialise first set of registers
		{
			RHS1.load_a(pRhs1 + i);
			RHS2.load_a(pRhs1 + i + width);
			RHS3.load_a(pRhs1 + i + width * 2);
			RHS4.load_a(pRhs1 + i + width * 3);

			RES = operTransform(RHS1);
			RES1 = operTransform(RHS2);
			RES2 = operTransform(RHS3);
			RES3 = operTransform(RHS4);
		}

		i += step;
		long rhsSZ = static_cast<long>(rhs1.size());
		for (; i <= (rhsSZ - step); i += step)
		{
			RHS1.load_a(pRhs1 + i);
			RHS2.load_a(pRhs1 + i + width);
			RHS3.load_a(pRhs1 + i + width * 2);
			RHS4.load_a(pRhs1 + i + width * 3);

			RES = operAcc(RES, operTransform(RHS1));
			RES1 = operAcc(RES1, operTransform(RHS2));
			RES2 = operAcc(RES2, operTransform(RHS3));
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
			RES = operAcc(RES, operTransform(RHS1));
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

 
//unitary 
//unrolled version helps greatly with VC2019
template< typename INS_VEC, typename OPT, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate2UR_X(const Vec<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc)
{
	return ApplyTransformAccumulate2UR_X_Impl(rhs1, operTransform, operAcc);
}


template< typename INS_VEC, typename OPT, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate2UR_X(const VecView<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc)
{
	return ApplyTransformAccumulate2UR_X_Impl(rhs1, operTransform, operAcc);
}


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


template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OPT, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate2UR_X_ImplBin(const VEC_TYPE<INS_VEC>& lhs1, const VEC_TYPE<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc)
{

	check_pair(lhs1, rhs1);
	auto zero = InstructionTraits<INS_VEC>::nullValue;
	//to do cover case of scalars
	if (isScalar(rhs1) && isScalar(lhs1))
	{
		return ApplyBinaryOperation1<INS_VEC, OP>(operTransform(lhs1.getScalarValue(), rhs1.getScalarValue()), zero, operAcc);
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

			RES = operTransform(LHS1, RHS1);
			RES1 = operTransform(LHS2, RHS2);
			RES2 = operTransform(LHS3, RHS3);
			RES3 = operTransform(LHS4, RHS4);
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
			RES = operAcc(RES, operTransform(LHS1, RHS1));

			LHS2.load_a(pLhs1 + i + width);
			RHS2.load_a(pRhs1 + i + width);
			RES1 = operAcc(RES1, operTransform(LHS2, RHS2));

			LHS3.load_a(pLhs1 + i + width * 2);
			RHS3.load_a(pRhs1 + i + width * 2);
			RES2 = operAcc(RES2, operTransform(LHS3, RHS3));

			LHS4.load_a(pLhs1 + i + width * 3);
			RHS4.load_a(pRhs1 + i + width * 3);
			RES3 = operAcc(RES3, operTransform(LHS4, RHS4));

		}

		// odd bits
		for (; i <= rhsSZ - width; i += width)
		{
			LHS1.load_a(pLhs1 + i);
			RHS1.load_a(pRhs1 + i);
			RES = operAcc(RES, operTransform(LHS1, RHS1));
		}

		RES = operAcc(RES, RES1);
		RES2 = operAcc(RES2, RES3);
		RES = operAcc(RES, RES2);

	}
	else
	{
		LHS1.load_a(pLhs1);
		RHS1.load_a(pRhs1);
		RES = operTransform(LHS1, RHS1);

		i += width;
		// odd bits
		for (; i <= sz - width; i += width)
		{
			LHS1.load_a(pLhs1 + i);
			RHS1.load_a(pRhs1 + i);
			RES = operAcc(RES, operTransform(LHS1, RHS1));
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




template< typename INS_VEC, typename OPT, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate2UR_XBin(const Vec<INS_VEC>& lhs1, const Vec<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc)
{

	return ApplyTransformAccumulate2UR_X_ImplBin(lhs1, rhs1, operTransform, operAcc);

}

template< typename INS_VEC, typename OPT, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate2UR_XBin(const VecView<INS_VEC>& lhs1, const VecView<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc)
{
	return ApplyTransformAccumulate2UR_X_ImplBin(lhs1, rhs1, operTransform, operAcc);
}


template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP >
VEC_TYPE<INS_VEC>  ApplyTransformUR_X_Impl(const VEC_TYPE<INS_VEC>& rhs1, OP& oper)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		return ApplyUnitaryOperation<INS_VEC, OP>(rhs1.getScalarValue(), oper);
	}

	auto pRhs1 = rhs1.start();
	VEC_TYPE<INS_VEC> ret(initTransformer(rhs1)); //ret(sz,1,1);// ret(rhs1); //???
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

	//to do odd end bit
	//no odd bit
	return ret;
}


template< typename INS_VEC, typename OP >
Vec<INS_VEC>  ApplyTransformUR_X(const Vec<INS_VEC>& rhs1, OP& oper)
{
	return ApplyTransformUR_X_Impl(rhs1, oper);
}


template< typename INS_VEC, typename OP >
VecView<INS_VEC> ApplyTransformUR_X(const VecView<INS_VEC>& rhs1, OP& oper)
{
	return ApplyTransformUR_X_Impl(rhs1, oper);
}


template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP >
void ApplyTransformUR_X_Impl(VEC_TYPE<INS_VEC>& rhs1, OP& oper)
{

	if (!rhs1.isScalar())
	{

		check_vector(rhs1); //calls overload with a view
		//views are not scalar

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


		int impSZ = rhs1.paddedSize();
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
		VEC_TYPE<INS_VEC> result;
		result = scalarRes;
		rhs1 = result;
	}

}

// experimental 
template<  template <class> typename VEC_TYPE, template <class> typename VEC_TYPE_RET, typename INS_VEC, typename OP, typename SAMPLER >
void ApplyTransformUR_X_Impl_EX(VEC_TYPE<INS_VEC>& rhs1, VEC_TYPE_RET<INS_VEC>& ret, OP& oper, SAMPLER& sampler, int i = 0, int impSZ = -1)
{

	impSZ = (impSZ < 0) ? static_cast<long>(rhs1.paddedSize()) : impSZ;


	auto pRhs1 = rhs1.start();
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int stride = static_cast<int>(sampler.stride());
	stride = (stride == 1) ? 1 : stride * width;
	int step = 4 * width * stride;




	SAMPLER RHS1(sampler);
	INS_VEC RES;


	SAMPLER RHS2(sampler);
	INS_VEC RES1;

	SAMPLER RHS3(sampler);
	INS_VEC RES2;

	SAMPLER RHS4(sampler);
	INS_VEC RES3;

	//we can only get a starting position bigger than  zero when we access points in the 
	// data preceeding the starting point, so we advance to a popint where we sample valid /existing data
	i = i + std::max(0, -sampler.min());



	//similarly if we are sampling  points beyond current index, we need to reduce maximum value iterated to so
	// that we stay in a valid range 
	impSZ = impSZ - std::max(0, sampler.max());

	int ld_offset = 0;
	int SV_offset = 0;

	int rhsSZ = impSZ - step;
	for (; i < rhsSZ; i += step)
	{
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(RHS1);
		RES.store(pRet + i + SV_offset);

		RHS2.load(pRhs1 + i + ld_offset + width * stride);
		RES1 = oper(RHS2);
		RES1.store(pRet + i + SV_offset + width);

		RHS3.load(pRhs1 + i + ld_offset + width * stride * 2);
		RES2 = oper(RHS3);
		RES2.store(pRet + i + SV_offset + width * 2);

		RHS4.load(pRhs1 + i + ld_offset + width * stride * 3);
		RES3 = oper(RHS4);
		RES3.store(pRet + i + SV_offset + width * 3);
	}
	//
	for (; i <= impSZ - width * stride; i += width * stride)
	{
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(RHS1);
		RES.store(pRet + i + SV_offset);
	}

	//one register case or do it scalar ops ?

	if (i < (impSZ - width * stride))
	{
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(RHS1);
		RES.store(pRet + i + SV_offset);
	}


	//iterate over remaining values and store
	RHS1.load(pRhs1 + i + ld_offset);
	RES = oper(RHS1);

	int j = 0;
	for (; i < impSZ; ++i, ++j)
	{
		(pRet + SV_offset)[i] = RES[j];
	}

}



// experimental  binary version for use with spans
template<  template <class> typename VEC_TYPE, template <class> typename VEC_TYPE_RET, typename INS_VEC, typename OP, typename SAMPLER >
void ApplyTransformUR_X_Impl_EX_BN(VEC_TYPE<INS_VEC>& lhs1, VEC_TYPE<INS_VEC>& rhs1, VEC_TYPE_RET<INS_VEC>& ret, OP& oper, SAMPLER& sampler_lhs, SAMPLER& sampler, int i = 0, int impSZ = -1)
{

	//TO DO SOME CHECKING on sizes etc
	//TAKING rhs as corect size in this initial go

	//NB strides can be different onleft and right too  2 DO

	impSZ = (impSZ < 0) ? static_cast<long>(rhs1.paddedSize()) : impSZ;

	auto pLhs1 = lhs1.start();
	auto pRhs1 = rhs1.start();
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int stride = static_cast<int>(sampler.stride());
	stride = (stride == 1) ? 1 : stride * width;
	int step = 4 * width * stride;



	SAMPLER LHS1(sampler_lhs);
	SAMPLER RHS1(sampler);
	INS_VEC RES;

	SAMPLER LHS2(sampler_lhs);
	SAMPLER RHS2(sampler);
	INS_VEC RES1;

	SAMPLER LHS3(sampler_lhs);
	SAMPLER RHS3(sampler);
	INS_VEC RES2;

	SAMPLER LHS4(sampler_lhs);
	SAMPLER RHS4(sampler);
	INS_VEC RES3;

	//we can only get a starting position bigger than  zero when we access points in the 
	// data preceeding the starting point, so we advance to a popint where we sample valid /existing data
	i = i + std::max(0, -sampler.min());



	//similarly if we are sampling  points beyond current index, we need to reduce maximum value iterated to so
	// that we stay in a valid range 
	impSZ = impSZ - std::max(0, sampler.max());

	int ld_offset = 0;
	int SV_offset = 0;

	int rhsSZ = impSZ - step;
	for (; i < rhsSZ; i += step)
	{
		LHS1.load(pLhs1 + i + ld_offset);
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(LHS1,RHS1);
		RES.store(pRet + i + SV_offset);

		LHS2.load(pLhs1 + i + ld_offset + width * stride);
		RHS2.load(pRhs1 + i + ld_offset + width * stride);
		RES1 = oper(LHS2,RHS2);
		RES1.store(pRet + i + SV_offset + width);

		LHS3.load(pLhs1 + i + ld_offset + width * stride * 2);
		RHS3.load(pRhs1 + i + ld_offset + width * stride * 2);
		RES2 = oper(LHS3,RHS3);
		RES2.store(pRet + i + SV_offset + width * 2);

		LHS4.load(pLhs1 + i + ld_offset + width * stride * 3);
		RHS4.load(pRhs1 + i + ld_offset + width * stride * 3);
		RES3 = oper(LHS4,RHS4);
		RES3.store(pRet + i + SV_offset + width * 3);
	}
	//
	for (; i <= impSZ - width * stride; i += width * stride)
	{
		LHS1.load(pLhs1 + i + ld_offset);
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(LHS1,RHS1);
		RES.store(pRet + i + SV_offset);
	}

	//one register case or do it scalar ops ?
	//

	if (i < (impSZ - width * stride))
	{
		LHS1.load(pLhs1 + i + ld_offset);
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(LHS1,RHS1);
		RES.store(pRet + i + SV_offset);
	}

	//in the case that all input ranges are al.igned  and padded and that 
	//spans cone from inside these ranges and the sampler is prevented 
	// from going out of these ranges then this is fine.


	// These mighty not be padded
	// 
	/*
	//iterate over remaining values and store
	RHS1.load(pRhs1 + i + ld_offset);
	RES = oper(RHS1);

	int j = 0;
	for (; i < impSZ; ++i, ++j)
	{
		(pRet + SV_offset)[i] = RES[j];
	}

	*/

	//iterate over remaining values and store
	/*
	if (stride == 1) // and on last element
	{
		LHS1.load(pLhs1 + i + ld_offset);
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(LHS1,RHS1);

		int j = 0;
		for (; i < impSZ; ++i, ++j)
		{
			(pRet + SV_offset)[i] = RES[j];
		}
	}
	else //load each element at a  time
	*/
	{
		//wrong 
		/*
		int K = 0;
		int j = 0;
		for (; i < impSZ; ++j)
		{
			LHS1.X_0.value = pLhs1[i + ld_offset];
			RHS1.X_0.value = pRhs1[i + ld_offset];
			RES = oper(LHS,RHS1);
			(pRet + i+ SV_offset)[K] = RES[0];
			i += stride;
			K++;
		}
		*/


		//int j = 0;
		for (; i < impSZ; ++i)
		{
			LHS1.X_0.value = pLhs1[i + ld_offset];
			RHS1.X_0.value = pRhs1[i + ld_offset];
			RES = oper(LHS1, RHS1);
			(pRet + i + SV_offset)[0] = RES[0];
			//i++;
			
			
		}


	}

}





// experimental 
template<  template <class> typename VEC_TYPE, template <class> typename VEC_TYPE_RET,typename INS_VEC, typename OP, typename SAMPLER >
void ApplyTransformUR_X_Impl_EX_STRD(VEC_TYPE<INS_VEC>& rhs1, VEC_TYPE_RET<INS_VEC>& ret, OP& oper, SAMPLER& sampler, int i = 0, int impSZ = -1)
{

	impSZ = (impSZ < 0) ? static_cast<long>(rhs1.paddedSize()) : impSZ;


	auto pRhs1 = rhs1.start();
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int stride = static_cast<int>(sampler.stride());


	int step = 4 * width * stride;


	SAMPLER RHS1(sampler);
	INS_VEC RES;


	SAMPLER RHS2(sampler);
	INS_VEC RES1;

	SAMPLER RHS3(sampler);
	INS_VEC RES2;

	SAMPLER RHS4(sampler);
	INS_VEC RES3;

	//we can only get a starting position bigger than  zero when we access points in the 
	// data preceeding the starting point, so we advance to a popint where we sample valid /existing data
	i = i+ std::max(0,-sampler.min());

	
	
	//similarly if we are sampling  points beyond current index, we need to reduce maximum value iterated to so
	// that we stay in a valid range 
	impSZ = impSZ -std::max(0,sampler.max());

	int ld_offset = 0;
	int SV_offset = 0;

	int rhsSZ = impSZ - step;
	int K = 0;
	for (; i < rhsSZ; i += step, K+=4*width)
	{
		RHS1.load(pRhs1 + i+ ld_offset);
		RES = oper(RHS1);
		RES.store(pRet + K + SV_offset);

		RHS2.load(pRhs1 + i + ld_offset + width * stride);
		RES1 = oper(RHS2);
		RES1.store(pRet + K + SV_offset + width);

		RHS3.load(pRhs1 + i + ld_offset + width * stride * 2);
		RES2 = oper(RHS3);
		RES2.store(pRet + K + SV_offset + width * 2);

		RHS4.load(pRhs1 + i + ld_offset + width * stride * 3);
		RES3 = oper(RHS4);
		RES3.store(pRet + K + SV_offset + width * 3);
	}
	

	for (; i <= impSZ - width * stride; i += width * stride, K+=width)
	{
		RHS1.load(pRhs1 + i + ld_offset );
		RES = oper(RHS1);
		RES.store(pRet + K + SV_offset);
	}

	//one register case or do it scalar ops ?

	if ( i < (impSZ- width * stride) )
	{
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(RHS1);
		RES.store(pRet + K + SV_offset);
	}

	

	//iterate over remaining values and store
	if (stride == 1)
	{
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(RHS1);

		int j = 0;
		for (; i < impSZ; ++i, ++j)
		{
			(pRet + SV_offset)[i] = RES[j];
		}
	}
	else //load each element at a  time
	{
	
		int j = 0;
		for (; i < impSZ;  ++j)
		{

			RHS1.X_0.value = pRhs1[i + ld_offset];
			RES = oper(RHS1);
			(pRet + SV_offset)[K] = RES[0];
			i += stride;
			K++;
		}

	}

}



// experimental binary transform  // need to check handling remainders ie not whole registers left at end of loop 
// since might get unpadded cases
//assumed to be unaligned
template<  template <class> typename VEC_TYPE, template <class> typename VEC_TYPE_AUX, template <class> typename VEC_TYPE_RET,typename INS_VEC, typename OP, typename SAMPLER >
void ApplyTransformUR_X_Impl_EX(const VEC_TYPE<INS_VEC>& rhs1,const VEC_TYPE_AUX<INS_VEC>& aux, VEC_TYPE_RET<INS_VEC>& ret, OP& oper, SAMPLER& sampler, int i = 0, int impSZ = -1)
{

	impSZ = (impSZ < 0) ? rhs1.paddedSize() : impSZ;

	auto pAux1 = aux.start();
	auto pRhs1 = rhs1.start();
	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	INS_VEC AUX1;
	SAMPLER RHS1(sampler);
	INS_VEC RES;

	INS_VEC AUX2;
	SAMPLER RHS2(sampler);
	INS_VEC RES1;

	INS_VEC AUX3;
	SAMPLER RHS3(sampler);
	INS_VEC RES2;

	INS_VEC AUX4;
	SAMPLER RHS4(sampler);
	INS_VEC RES3;

	//int i = 0;

	//we can only get a starting position bigger than  zero when we access points in the 
	// data preceeding the starting point, so we advance to a popint where we sample valid /existing data
	i = i + std::max(0, -sampler.min());

	//similarly if we are sampling  points beyond current index, we need to reduce maximum value iterated to so
	// that we stay in a valid range 
	impSZ = impSZ - std::max(0, sampler.max());

	int ld_offset = 0;
	int SV_offset = 0;

	int rhsSZ = impSZ - step;
	for (; i < rhsSZ; i += step)
	{
		AUX1.load(pAux1 + i + ld_offset);
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(RHS1,AUX1);
		RES.store(pRet + i + SV_offset);

		AUX2.load(pAux1 + i + ld_offset + width);
		RHS2.load(pRhs1 + i + ld_offset + width);
		RES1 = oper(RHS2, AUX2);
		RES1.store(pRet + i + SV_offset + width);

		AUX3.load(pAux1 + i + ld_offset + width * 2);
		RHS3.load(pRhs1 + i + ld_offset + width * 2);
		RES2 = oper(RHS3, AUX3);
		RES2.store(pRet + i + SV_offset + width * 2);

		AUX4.load(pAux1 + i + ld_offset + width * 3);
		RHS4.load(pRhs1 + i + ld_offset + width * 3);
		RES3 = oper(RHS4, AUX4);
		RES3.store(pRet + i + SV_offset + width * 3);
	}

	for (; i <= impSZ - width; i += width)
	{
		AUX1.load(pAux1 + i + ld_offset);
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(RHS1, AUX1);
		RES.store(pRet + i + SV_offset);
	}

	//one register case or do it scalar ops ?

	if (i < (impSZ - width))
	{
		AUX1.load(pAux1 + i + ld_offset);
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(RHS1,AUX1);
		RES.store(pRet + i + SV_offset);
	}

	//move to one register width from last valid
	//point to calculate // fine if both are padded and using aligned so probably an error
	i = impSZ - width;
	AUX1.load(pAux1 + i + ld_offset);
	RHS1.load(pRhs1 + i + ld_offset);
	RES = oper(RHS1, AUX1);
	RES.store(pRet + i + SV_offset);

}


template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OPT, typename OP, typename SAMPLER>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformAccumulate2UR_X_ImplBin_EX(const VEC_TYPE<INS_VEC>& lhs1, const VEC_TYPE<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc,
	SAMPLER sampler_lhs, SAMPLER sampler_rhs, int i = 0, int impSZ = -1)
{

	check_pair(lhs1, rhs1);
	auto zero = InstructionTraits<INS_VEC>::nullValue;
	//to do cover case of scalars
//	if (isScalar(rhs1) && isScalar(lhs1))
//	{
//		return ApplyBinaryOperation1<INS_VEC, OP>(operTransform(lhs1.getScalarValue(), rhs1.getScalarValue()), zero, operAcc);
//	}


	auto pRhs1 = rhs1.start();
	auto pLhs1 = lhs1.start();
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;

	int sz = static_cast<int>(rhs1.size());
	Vec<INS_VEC> ret(sz);

	SAMPLER LHS1;
	LHS1.X_0.value = zero;
	SAMPLER RHS1;
	RHS1.X_0.value = zero;
	INS_VEC RES = zero;

	SAMPLER LHS2;
	LHS2.X_0.value = zero;
	SAMPLER RHS2;
	RHS2.X_0.value = zero;
	INS_VEC RES1 = zero;

	SAMPLER LHS3;
	LHS3.X_0.value = zero;
	SAMPLER RHS3;
	RHS3.X_0.value = zero;
	INS_VEC RES2 = zero;

	SAMPLER LHS4;
	LHS4.X_0.value = zero;
	SAMPLER RHS4;
	RHS4.X_0.value = zero;
	INS_VEC RES3 = zero;



	if (sz >= step * 2)
	{
		//initialise first set of registers
		{
			RHS1.load(pRhs1 + i);
			RHS2.load(pRhs1 + i + width);
			RHS3.load(pRhs1 + i + width * 2);
			RHS4.load(pRhs1 + i + width * 3);

			LHS1.load(pLhs1 + i);
			LHS2.load(pLhs1 + i + width);
			LHS3.load(pLhs1 + i + width * 2);
			LHS4.load(pLhs1 + i + width * 3);

			RES = operTransform(LHS1, RHS1);
			RES1 = operTransform(LHS2, RHS2);
			RES2 = operTransform(LHS3, RHS3);
			RES3 = operTransform(LHS4, RHS4);
		}

		i += step;
		
		int impSZ = static_cast<int>(rhs1.paddedSize());
		
		int rhsSZ = impSZ - step;

		
		for (; i <= (impSZ - step); i += step)
		{
			LHS1.load(pLhs1 + i);
			RHS1.load(pRhs1 + i);
			RES = operAcc(RES, operTransform(LHS1, RHS1));

			LHS2.load(pLhs1 + i + width);
			RHS2.load(pRhs1 + i + width);
			RES1 = operAcc(RES1, operTransform(LHS2, RHS2));

			LHS3.load(pLhs1 + i + width * 2);
			RHS3.load(pRhs1 + i + width * 2);
			RES2 = operAcc(RES2, operTransform(LHS3, RHS3));

			LHS4.load(pLhs1 + i + width * 3);
			RHS4.load(pRhs1 + i + width * 3);
			RES3 = operAcc(RES3, operTransform(LHS4, RHS4));

		}

		// odd bits
		for (; i <= impSZ - width; i += width)
		{
			LHS1.load(pLhs1 + i);
			RHS1.load(pRhs1 + i);
			RES = operAcc(RES, operTransform(LHS1, RHS1));
		}

		RES = operAcc(RES, RES1);
		RES2 = operAcc(RES2, RES3);
		RES = operAcc(RES, RES2);

	}
	else
	{
		LHS1.load(pLhs1);
		RHS1.load(pRhs1);
		RES = operTransform(LHS1, RHS1);

		i += width;
		// odd bits
		for (; i <= sz - width; i += width)
		{
			LHS1.load(pLhs1 + i);
			RHS1.load(pRhs1 + i);
			RES = operAcc(RES, operTransform(LHS1, RHS1));
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





template< typename INS_VEC, typename OP >
void ApplyTransformUR_X(VecView<INS_VEC>& rhs1, OP& oper)
{
	ApplyTransformUR_X_Impl(rhs1, oper);

}

template< typename INS_VEC, typename OP >
void ApplyTransformUR_X(Vec<INS_VEC>& rhs1, OP& oper)
{
	ApplyTransformUR_X_Impl(rhs1, oper);
}


template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP >
VEC_TYPE<INS_VEC>  ApplyTransformUR_XX_Impl(const VEC_TYPE<INS_VEC>& rhs1, OP& oper)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		return ApplyUnitaryOperation<INS_VEC, OP>(rhs1.getScalarValue(), oper);
	}

	
	auto pRhs1 = rhs1.start();
	VEC_TYPE<INS_VEC> ret(initTransformer(rhs1));
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


	int impSZ = rhs1.paddedSize();

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

		RHS6.load_a(pRhs1 + i + width * 5);
		RES5 = oper(RHS6);
		RES5.store_a(pRet + i + width * 5);

		RHS7.load_a(pRhs1 + i + width * 6);
		RES6 = oper(RHS7);
		RES6.store_a(pRet + i + width * 6);

		RHS8.load_a(pRhs1 + i + width * 7);
		RES7 = oper(RHS8);
		RES7.store_a(pRet + i + width * 7);
	}

	for (; i <= impSZ - width; i += width)
	{
		RHS1.load_a(pRhs1 + i);
		RES = oper(RHS1);
		RES.store_a(pRet + i);
	}

	return ret;
}


template< typename INS_VEC, typename OP >
Vec<INS_VEC>  ApplyTransformUR_XX(const Vec<INS_VEC>& rhs1, OP& oper)
{
	return ApplyTransformUR_XX_Impl(rhs1, oper);
}


template< typename INS_VEC, typename OP >
VecView<INS_VEC>  ApplyTransformUR_XX( const VecView<INS_VEC>& rhs1, OP& oper)
{
	return ApplyTransformUR_XX_Impl(rhs1, oper);

}


template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP >
VEC_TYPE<INS_VEC>  ApplyBinaryTransformUR_X_Impl(const VEC_TYPE<INS_VEC>& lhs, const VEC_TYPE<INS_VEC>& rhs, OP& oper)
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

	//int sz = rhs.size();
	auto pRhs = rhs.start();
	auto pLhs = lhs.start();
	VEC_TYPE<INS_VEC> ret(initTransformer(lhs));
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

	int impSZ = lhs.paddedSize();
	int rhsSZ = impSZ - step;

	for (; i < rhsSZ; i += step)
	{
		LHS.load_a(pLhs + i);
		RHS.load_a(pRhs + i);
		RES = oper(LHS, RHS);
		RES.store_a(pRet + i);

		LHS1.load_a(pLhs + i + width);
		RHS1.load_a(pRhs + i + width);
		RES1 = oper(LHS1, RHS1);
		RES1.store_a(pRet + i + width);

		LHS2.load_a(pLhs + i + 2 * width);
		RHS2.load_a(pRhs + i + 2 * width);
		RES2 = oper(LHS2, RHS2);
		RES2.store_a(pRet + i + 2 * width);

		LHS3.load_a(pLhs + i + 3 * width);
		RHS3.load_a(pRhs + i + 3 * width);
		RES3 = oper(LHS3, RHS3);
		RES3.store_a(pRet + i + 3 * width);
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
Vec<INS_VEC>  ApplyBinaryTransformUR_X(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs, OPER& oper)
{
	return ApplyBinaryTransformUR_X_Impl(lhs, rhs, oper);
}


template< typename INS_VEC, typename OPER >
VecView<INS_VEC>  ApplyBinaryTransformUR_X(const VecView<INS_VEC>& lhs, const VecView<INS_VEC>& rhs, OPER& oper)
{
	return ApplyBinaryTransformUR_X_Impl(lhs, rhs, oper);
}


template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP >
VEC_TYPE<INS_VEC>  ApplyBinaryTransformUR_X_Impl(typename InstructionTraits<INS_VEC>::FloatType  lhs, const VEC_TYPE<INS_VEC>& rhs, OP& oper)
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


	//int sz = rhs.size();
	auto pRhs = rhs.start();
	VEC_TYPE<INS_VEC> ret(initTransformer(rhs));
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
		RES2.store_a(pRet + i + 2 * width);

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
VecView<INS_VEC>  ApplyBinaryTransformUR_X(typename InstructionTraits<INS_VEC>::FloatType  lhs, const VecView<INS_VEC>& rhs, OPER& oper)
{
	return ApplyBinaryTransformUR_X_Impl(lhs,rhs, oper);
}


template< typename INS_VEC, typename OPER >
Vec<INS_VEC>  ApplyBinaryTransformUR_X(typename InstructionTraits<INS_VEC>::FloatType  lhs, const Vec<INS_VEC>& rhs, OPER& oper)
{
	return ApplyBinaryTransformUR_X_Impl(lhs,rhs, oper);
}


template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP >
VEC_TYPE<INS_VEC>  ApplyBinaryTransformUR_X_Impl(const VEC_TYPE<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType  rhs, OP& oper)
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
		return VEC_TYPE<INS_VEC>(RES[0]);
	}


	//int sz = lhs.size();
	auto pLhs = lhs.start();
	VEC_TYPE<INS_VEC> ret(initTransformer(lhs));// ret(sz);
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
		RES1.store_a(pRet + i + width);

		LHS2.load_a(pLhs + i + 2 * width);
		RES2 = oper(LHS2, RHS);
		RES2.store_a(pRet + i + 2 * width);

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


template< typename INS_VEC, typename OPER >
VecView<INS_VEC>  ApplyBinaryTransformUR_X(const VecView<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType  rhs, OPER& oper)
{
	return ApplyBinaryTransformUR_X_Impl(lhs, rhs, oper);
}


template< typename INS_VEC, typename OPER >
Vec<INS_VEC>  ApplyBinaryTransformUR_X(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType  rhs,  OPER& oper)
{
	return ApplyBinaryTransformUR_X_Impl(lhs, rhs, oper);

}


template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP, typename OPER_TRUE, typename OPER_FALSE >
VEC_TYPE<INS_VEC>  ApplySelectTransformUR_X_Impl(const VEC_TYPE<INS_VEC>& rhs1, OP& cond, OPER_TRUE& trueOper, OPER_FALSE& falseOper)
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
		return VEC_TYPE<INS_VEC>(RES[0]);
	}

	auto pRhs1 = rhs1.start();
	VEC_TYPE<INS_VEC> ret(initTransformer(rhs1));
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
	int impSZ = rhs1.paddedSize();
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


template< typename INS_VEC, typename OP, typename OPER_TRUE, typename OPER_FALSE >
Vec<INS_VEC>  ApplySelectTransformUR_X(const Vec<INS_VEC>& rhs1, OP& cond, OPER_TRUE& trueOper, OPER_FALSE& falseOper  )
{
	return ApplySelectTransformUR_X_Impl(rhs1, cond, trueOper, falseOper);
}


template< typename INS_VEC, typename OP, typename OPER_TRUE, typename OPER_FALSE >
VecView<INS_VEC>  ApplySelectTransformUR_X(const VecView<INS_VEC>& rhs1, OP& cond, OPER_TRUE& trueOper, OPER_FALSE& falseOper)
{
	return ApplySelectTransformUR_X_Impl(rhs1, cond, trueOper, falseOper);
}


template< typename INS_VEC, typename BOOL_OPER, typename TRUE_OPER, typename FALSE_OPER>
Vec<INS_VEC> ApplySelectionOperationFuncUR_X(BOOL_OPER& COND, const Vec<INS_VEC>& testData, TRUE_OPER& trueOper, FALSE_OPER& falseOper)
{
	return ApplySelectTransformUR_X(testData, COND, trueOper, falseOper);
}


template<  template <class> typename VEC_TYPE , typename INS_VEC, typename OP >
VEC_TYPE<INS_VEC>  ApplySelectTransformUR_XC_Impl(const VEC_TYPE<INS_VEC >& rhs1, OP& cond, typename InstructionTraits<INS_VEC>::FloatType trueVal, typename InstructionTraits<INS_VEC>::FloatType& falseVal)
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

	auto pRhs1 = rhs1.start();
	VEC_TYPE<INS_VEC> ret(initTransformer(rhs1));
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
	int impSZ = rhs1.paddedSize();
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


template< typename INS_VEC, typename OP >
Vec<INS_VEC>  ApplySelectTransformUR_XC(const Vec<INS_VEC>& rhs1, OP& cond, typename InstructionTraits<INS_VEC>::FloatType trueVal, typename InstructionTraits<INS_VEC>::FloatType& falseVal)
{
	return ApplySelectTransformUR_XC_Impl(rhs1, cond, trueVal, falseVal);
}


template< typename INS_VEC, typename OP >
VecView<INS_VEC>  ApplySelectTransformUR_XC(const VecView<INS_VEC>& rhs1, OP& cond, typename InstructionTraits<INS_VEC>::FloatType trueVal, typename InstructionTraits<INS_VEC>::FloatType& falseVal)
{
	return ApplySelectTransformUR_XC_Impl(rhs1, cond, trueVal, falseVal);
}





//////experimental unrolled  double accumulation
template< typename INS_VEC, typename OP, typename OP2>
typename std::tuple<typename InstructionTraits<INS_VEC>::FloatType, typename InstructionTraits<INS_VEC>::FloatType>
 ApplyAccumulate2UR_X2(const Vec<INS_VEC>& rhs1, OP& oper, OP2& oper2)
{
	check_vector(rhs1);
	if (isScalar(rhs1)) // nothing to accumulate with so just return  value
	{
		return { rhs1.getScalarValue(),rhs1.getScalarValue() };
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

	//
	INS_VEC RHS1_A;
	INS_VEC RES_A;

	INS_VEC RHS2_A;
	INS_VEC RES1_A;

	INS_VEC RHS3_A;
	INS_VEC RES2_A;

	INS_VEC RHS4_A;
	INS_VEC RES3_A;

	//

	int i = 0;

	if (sz >= step * 2)
	{
		//initialise first set of registers
		{
			RES.load_a(pRhs1 + i);
			RES1.load_a(pRhs1 + i + width);
			RES2.load_a(pRhs1 + i + width * 2);
			RES3.load_a(pRhs1 + i + width * 3);


			RES_A = RES;
			RES1_A = RES1;
			RES2_A = RES2;
			RES3_A = RES3;
		}

		i += step;
		int rhsSZ = rhs1.size();
		for (; i <= (rhsSZ - step); i += step)
		{
			RHS1.load_a(pRhs1 + i);
			RES = oper(RES, RHS1);
			RES_A = oper2(RES_A, RHS1);

			RHS2.load_a(pRhs1 + i + width);
			RES1 = oper(RES1, RHS2);
			RES1_A = oper2(RES1_A, RHS2);

			RHS3.load_a(pRhs1 + i + width * 2);
			RES2 = oper(RES2, RHS3);
			RES2_A = oper2(RES2_A, RHS3);

			RHS4.load_a(pRhs1 + i + width * 3);
			RES3 = oper(RES3, RHS4);
			RES3_A = oper2(RES3_A, RHS4);

		}

		// odd bits
		for (; i <= rhsSZ - width; i += width)
		{
			RHS1.load_a(pRhs1 + i);
			RES = oper(RES, RHS1);
			RES_A = oper2(RES_A, RHS1);
		}

		RES = oper(RES, RES1);
		RES2 = oper(RES2, RES3);
		RES = oper(RES, RES2);

		RES_A = oper2(RES_A, RES1_A);
		RES2_A = oper2(RES2_A, RES3_A);
		RES_A = oper2(RES_A, RES2_A);

	}
	else
	{
		RES.load_a(pRhs1);
		RES_A = RES;

		i += width;
		// odd bits
		for (; i <= sz - width; i += width)
		{
			RHS1.load_a(pRhs1 + i);
			RES = oper(RES, RHS1);
			RES_A = oper2(RES_A, RHS1);
		}

	}

	typename InstructionTraits<INS_VEC>::FloatType result = RES[0];
	typename InstructionTraits<INS_VEC>::FloatType result_A = RES_A[0];
	int min_wdth = std::min(sz, width);
	//across vectors lanes  // not assuming horizontal version exist
	for (int j = 1; j < min_wdth; ++j)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(result, RES[j], oper);
		result_A = ApplyBinaryOperationVec<INS_VEC, OP2>(result_A, RES_A[j], oper2);
	}

	//end bits for vecs not filling padding
	for (; i < rhs1.size(); ++i)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(pRhs1[i], result, oper);
		result_A = ApplyBinaryOperationVec<INS_VEC, OP2>(pRhs1[i], result_A, oper2);
	}

	return { result, result_A };
}




//////////////

//////experimental unrolled  double transform accumulation
template< typename INS_VEC, typename TF, typename TF2, typename OP, typename OP2>
typename std::tuple<typename InstructionTraits<INS_VEC>::FloatType, typename InstructionTraits<INS_VEC>::FloatType>
ApplyTransformAccumulate2UR_X2(const Vec<INS_VEC>& rhs1, TF& opTrans, OP& oper, TF2& opTrans2, OP2& oper2)
{
	check_vector(rhs1);
	if (isScalar(rhs1)) // nothing to accumulate with so just return  value
	{
		return { rhs1.getScalarValue(),rhs1.getScalarValue() };
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

	//
	INS_VEC RHS1_A;
	INS_VEC RES_A;

	INS_VEC RHS2_A;
	INS_VEC RES1_A;

	INS_VEC RHS3_A;
	INS_VEC RES2_A;

	INS_VEC RHS4_A;
	INS_VEC RES3_A;

	//

	int i = 0;

	if (sz >= step * 2)
	{
		//initialise first set of registers
		{
			RES.load_a(pRhs1 + i);
			RES1.load_a(pRhs1 + i + width);
			RES2.load_a(pRhs1 + i + width * 2);
			RES3.load_a(pRhs1 + i + width * 3);


			RES_A = opTrans2(RES);
			RES1_A = opTrans2(RES1);
			RES2_A = opTrans2(RES2);
			RES3_A = opTrans2(RES3);


			RES = opTrans(RES);
			RES1 = opTrans(RES1);
			RES2 = opTrans(RES2);
			RES3 = opTrans(RES3);
		}

		i += step;
		int rhsSZ = rhs1.size();
		for (; i <= (rhsSZ - step); i += step)
		{
			RHS1.load_a(pRhs1 + i);
			RES = oper(RES, opTrans(RHS1));
			RES_A = oper2(RES_A, opTrans2(RHS1));

			RHS2.load_a(pRhs1 + i + width);
			RES1 = oper(RES1, opTrans(RHS2));
			RES1_A = oper2(RES1_A, opTrans2(RHS2));

			RHS3.load_a(pRhs1 + i + width * 2);
			RES2 = oper(RES2, opTrans(RHS3));
			RES2_A = oper2(RES2_A, opTrans2(RHS3));

			RHS4.load_a(pRhs1 + i + width * 3);
			RES3 = oper(RES3, opTrans(RHS4));
			RES3_A = oper2(RES3_A, opTrans2(RHS4));

		}

		// odd bits
		for (; i <= rhsSZ - width; i += width)
		{
			RHS1.load_a(pRhs1 + i);
			RES = oper(RES, opTrans(RHS1));
			RES_A = oper2(RES_A, opTrans2(RHS1));
		}

		RES = oper(RES, RES1);
		RES2 = oper(RES2, RES3);
		RES = oper(RES, RES2);

		RES_A = oper2(RES_A, RES1_A);
		RES2_A = oper2(RES2_A, RES3_A);
		RES_A = oper2(RES_A, RES2_A);

	}
	else
	{
		RES.load_a(pRhs1);
		RES_A = opTrans2(RES);
		RES= opTrans(RES);

		i += width;
		// odd bits
		for (; i <= sz - width; i += width)
		{
			RHS1.load_a(pRhs1 + i);
			RES = oper(RES, opTrans(RHS1));
			RES_A = oper2(RES_A, opTrans2(RHS1));
		}

	}

	typename InstructionTraits<INS_VEC>::FloatType result = RES[0];
	typename InstructionTraits<INS_VEC>::FloatType result_A = RES_A[0];
	int min_wdth = std::min(sz, width);
	//across vectors lanes  // not assuming horizontal version exist
	for (int j = 1; j < min_wdth; ++j)
	{
		result = ApplyBinaryOperationVec<INS_VEC, OP>(result, RES[j], oper);
		result_A = ApplyBinaryOperationVec<INS_VEC, OP2>(result_A, RES_A[j], oper2);
	}

	//end bits for vecs not filling padding
	for (; i < rhs1.size(); ++i)
	{
		//this is a bug need to transform
		result = ApplyBinaryOperationVec<INS_VEC, OP>(pRhs1[i], result, oper);
		result_A = ApplyBinaryOperationVec<INS_VEC, OP2>(pRhs1[i], result_A, oper2);
	}

	return { result, result_A };
}


