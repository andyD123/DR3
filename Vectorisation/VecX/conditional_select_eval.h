/****************************  conditional_select_eval.h  *******************************
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
#include  "filter_select.h"
#include "accumulate_transform.h"
#include "error_utils.h"
#include <tuple>


/*
 applies the tesFunc to the vector val, if the func returns true  it applies the trueLambda to the value otherwise it applies the falseLambda
 can be slow if true/false Lambda functions are not heavyweight

  // rename filter Transform
*/
template< typename INS_VEC, typename  BOOL_TEST_OP, typename  TRUE_LAMBDA, typename  FALSE_LAMBDA>
Vec<INS_VEC> splitConditionalCalculate(const Vec<INS_VEC>& val, BOOL_TEST_OP& testFunc, TRUE_LAMBDA& trueLambda, FALSE_LAMBDA& falseLambda)
{

	check_vector(val);

	if (!val.isScalar())
	{
		
		auto  vwTupple = ApplyBinaryFilter(testFunc, val);
		ApplyUnitaryOperation(trueLambda, std::get<0>(vwTupple));
		ApplyUnitaryOperation(falseLambda, std::get<1>(vwTupple));
		return  merge(vwTupple);
	}
	else
	{
		return selectTransform(testFunc, val, trueLambda, falseLambda);
	}
}




/*
 applies the tesFunc to the vector val, if the func returns true  it applies the trueLambda to the value otherwise it applies the falseLambda
 can be slow if true/false Lambda functions are not heavyweight
 //using unrolled transform for VS2019
   // rename filter Transform
*/
template< typename INS_VEC, typename  BOOL_TEST_OP, typename  TRUE_LAMBDA, typename  FALSE_LAMBDA>
Vec<INS_VEC> splitConditionalCalculate_X(const Vec<INS_VEC>& val, BOOL_TEST_OP& testFunc, TRUE_LAMBDA& trueLambda, FALSE_LAMBDA& falseLambda)
{
	check_vector(val);

	if (!val.isScalar())
	{
		auto  vwTupple = ApplyBinaryFilter(testFunc, val);
		ApplyTransformUR_X(std::get<0>(vwTupple), trueLambda);
		ApplyTransformUR_X(std::get<1>(vwTupple), falseLambda);
		return  merge(vwTupple);
	}
	else
	{
		return selectTransform(testFunc, val, trueLambda, falseLambda);
	}
}




/*
 applies the tesFunc to the view val, if the func returns true  it applies the trueLambda to the value otherwise it applies the falseLambda
 can be slow if true/false Lambda functions are not heavyweight
 // rename filter Transform

*/
template< typename INS_VEC, typename  BOOL_TEST_OP, typename  TRUE_LAMBDA, typename  FALSE_LAMBDA>
VecView<INS_VEC> splitConditionalCalculate(const VecView<INS_VEC>& val, BOOL_TEST_OP& testFunc, TRUE_LAMBDA& trueLambda, FALSE_LAMBDA& falseLambda)
{

	check_vector(val);

	//if (!val.isScalar())
	//{
		auto vwTupple = ApplyBinaryFilter(testFunc, val);
		ApplyUnitaryOperation(trueLambda, std::get<0>(vwTupple));
		ApplyUnitaryOperation(falseLambda, std::get<1>(vwTupple));
		return  merge(vwTupple);

/*
	}
	else
	{
		
		//return selectTransform(testFunc, val, trueLambda, falseLambda);
	}
	*/
}

/*
if the condition at point [i] in the bool condition vector is true,  sets the ith element of result value to the  trueVals[i] otherwise falseVals[i]
*/
template< typename INS_VEC>
Vec<INS_VEC> ApplySelectionOperation(const VecBool<INS_VEC>& condition, const Vec<INS_VEC>& trueVals, const Vec<INS_VEC>& falseVals)
{
	check_pair(trueVals, falseVals);
	check_pair_different_type(condition, trueVals);

	Vec<INS_VEC> result(trueVals.size());
	auto pRes = result.start();
	auto pCond = condition.start();
	auto pLhs = trueVals.start();
	auto pRhs = falseVals.start();
	int sz = trueVals.paddedSize();
	Unroll_Select< INS_VEC, UnUsed>::apply_4(sz, pCond, pLhs, pRhs, pRes);
	return result;
}

/*
* selecting between constant values from vectors
Applies the boolean lambda to the data at point [i] in the testData vector, if the result is true it sets the ith element of the result  to the value 
trueVals[i] otherwise falseVals[i]
*/
template< typename INS_VEC, typename BOOL_OPER>
Vec<INS_VEC> ApplySelectionOperationC(BOOL_OPER& COND, const Vec<INS_VEC>& testData, const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	check_pair(lhs, rhs);
	check_pair(lhs, testData);


	if (isScalar(testData) &&isScalar(rhs) && isScalar(lhs)  )
	{
		INS_VEC RHS = testData.getScalarValue();
		INS_VEC TRU = lhs.getScalarValue();
		INS_VEC FLS = rhs.getScalarValue();
		INS_VEC RES = select(COND(RHS), TRU, FLS);
		return Vec<INS_VEC>(RES[0]);
	}


	Vec<INS_VEC> result(rhs.size());
	auto pRes = result.start();
	auto pLhs = lhs.start();
	auto pRhs = rhs.start();
	auto pX = testData.start();
	int sz = lhs.paddedSize();
	Unroll_Select< INS_VEC, BOOL_OPER>::apply_4( sz, pX, pLhs, pRhs, pRes, COND);
	return result;
}

/*
* selecting between scalar constant values
Applies the boolean lambda to the data at point [i] in the testData vector, if the result is true it sets the ith element of the result  to the value
trueVal otherwise falseVal

//rename as select
*/
template< typename INS_VEC, typename BOOL_OPER>
Vec<INS_VEC> ApplySelectionOperationC(BOOL_OPER& COND, const Vec<INS_VEC>& testData, typename InstructionTraits<INS_VEC>::FloatType trueVal, typename InstructionTraits<INS_VEC>::FloatType falseVal)
{
	check_vector(testData);
	Vec<INS_VEC> result(testData.size());
	auto pRes = result.start();
	auto pX = testData.start();
	int sz = testData.paddedSize();
	Unroll_Select< INS_VEC, BOOL_OPER>::apply_4(sz, pX, trueVal, falseVal, pRes, COND);
	return result;
}


/*
* calculating a value using trueOper or falseOper 
Applies the boolean lambda to the data at point [i] in the testData vector, if the result is trueOper(testData[i]) otherwise falseOper(testData[i])

rename as selectTransform
*/
template< typename INS_VEC, typename BOOL_OPER, typename TRUE_OPER, typename FALSE_OPER>
Vec<INS_VEC> ApplySelectionOperationFunc(BOOL_OPER& COND, const Vec<INS_VEC>& testData, TRUE_OPER& trueOper, FALSE_OPER& falseOper)
{
	check_vector(testData);
	Vec<INS_VEC> result(testData.size());
	auto pRes = result.start();
	auto pX = testData.start();
	int sz = testData.paddedSize();
	Unroll_Select< INS_VEC, BOOL_OPER, TRUE_OPER, FALSE_OPER>::apply_4(sz, pX, trueOper, falseOper, pRes, COND);
	return result;
}

