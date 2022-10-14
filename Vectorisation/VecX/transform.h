#pragma once

#include "vec.h"
#include "binary_unitary_operations.h"
#include "conditional_select_eval.h"

//Vec  applyTransform(lambda, const inputVec&, int UR);//unroll
//void  applyTransform(lambda inputVec&, int UR); //unroll

//Unitary lambdas
//unrolled version defaults to unroll x4
template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  ApplyTransform(LAMBDA& lambda, const Vec<INS_VEC>& inputVec)
{
	return ApplyUnitaryOperation(inputVec, lambda); 
}

//not unrolled x1
template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  ApplyTransform1(LAMBDA& lambda, const Vec<INS_VEC>& inputVec)// , int UR)
{
	return ApplyUnitaryOperation1(inputVec, lambda);
}


//inplace transforms

template<typename LAMBDA, typename INS_VEC>
void ApplyTransformM(LAMBDA& lambda,  Vec<INS_VEC>& inputVec)
{
	ApplyUnitaryOperationM(inputVec, lambda);
}

//not unrolled x1
template<typename LAMBDA, typename INS_VEC>
void  ApplyTransform1(LAMBDA& lambda, Vec<INS_VEC>& inputVec)// , int UR)
{
	ApplyUnitaryOperation1(inputVec, lambda);
}


//binary lambdas

//Unitary lambdas
//unrolled version defaults to unroll x4
template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  ApplyTransform(LAMBDA& lambda, const Vec<INS_VEC>& inputVecLHS, const Vec<INS_VEC>& inputVecRHS)
{
	return ApplyBinaryOperation(inputVecLHS, inputVecRHS, lambda);
}

template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  ApplyTransform(LAMBDA& lambda, typename InstructionTraits<INS_VEC>::FloatType LHS, const Vec<INS_VEC>& inputVecRHS)
{
	return ApplyBinaryOperation(LHS, inputVecRHS, lambda);
}

template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  ApplyTransform(LAMBDA& lambda, const Vec<INS_VEC>& inputVecLHS, typename InstructionTraits<INS_VEC>::FloatType RHS)
{
	return ApplyBinaryOperation(inputVecLHS, RHS, lambda);
}

// not unrolled
template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  ApplyTransform1(LAMBDA& lambda, const Vec<INS_VEC>& inputVecLHS, const Vec<INS_VEC>& inputVecRHS)
{
	return ApplyBinaryOperation1(inputVecLHS, inputVecRHS, lambda);
}

//conversion 
template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  ApplyTransform1(LAMBDA& lambda, typename InstructionTraits<INS_VEC>::FloatType LHS, const Vec<INS_VEC>& inputVecRHS)
{
	return ApplyBinaryOperation1(LHS, inputVecRHS, lambda);
}

template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  ApplyTransform1(LAMBDA& lambda, const Vec<INS_VEC>& inputVecLHS, typename InstructionTraits<INS_VEC>::FloatType RHS)
{
	return ApplyBinaryOperation1(inputVecLHS, RHS, lambda);
}

/////////////////////////////// in place binary unroll x 4 ////////////////

template<typename LAMBDA, typename INS_VEC>
void ApplyTransformM(LAMBDA& lambda,  Vec<INS_VEC>& inputVecLHS, const Vec<INS_VEC>& inputVecRHS)
{
	ApplyBinaryOperationMMXY<INS_VEC, LAMBDA>(inputVecLHS, inputVecRHS , lambda);
}

template<typename LAMBDA, typename INS_VEC>
void ApplyTransformM(LAMBDA& lambda, Vec<INS_VEC>& inputVecLHS, typename InstructionTraits<INS_VEC>::FloatType RHS)
{
	ApplyBinaryOperationMMMX<INS_VEC, LAMBDA>(inputVecLHS, RHS, lambda);
}


template<typename LAMBDA, typename INS_VEC>
void ApplyTransformM(LAMBDA& lambda, typename InstructionTraits<INS_VEC>::FloatType LHS,  Vec<INS_VEC>& inputVecRHS)
{
	ApplyBinaryOperationMMM<INS_VEC, LAMBDA>(LHS, inputVecRHS , lambda);
}


////////////////////////////////////////////


//    ApplySparseTransform  takes a boolean condition lambda to determine if the lambda should be used
//    to calculate a given value in the vector
//    if so it applies the transform lambda and then blends it into the target vector 
template< typename LAMBDA, typename INS_VEC, typename CONDITION_LAMBDA>
void ApplySparseTransform(const Vec<INS_VEC>& inputVec, Vec<INS_VEC>& updateResult, LAMBDA& oper, CONDITION_LAMBDA& selectionOp)
{
	ApplySparseUnitaryOperationU(inputVec, updateResult,oper, selectionOp);
}


///////////////////// TO DO  make this the interface /////////////////////

template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType ApplyReduce(const Vec<INS_VEC>& rhs1, OP& oper, typename InstructionTraits<INS_VEC>::FloatType initVal, bool singularInit = true)
{
	return ApplyAccumulate2(rhs1, oper, initVal, singularInit);

}


template< typename INS_VEC, typename OP, typename OPT>
typename InstructionTraits<INS_VEC>::FloatType ApplyTransformReduce(const Vec<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{
	return ApplyTransformAccumulateUR(rhs1, operTransform, operAcc, initVal, singularInit);
}


/////////////////////////////////


template< typename INS_VEC, typename BOOL_OPER>
Vec<INS_VEC> ApplySelection(BOOL_OPER& COND, const Vec<INS_VEC>& testData, const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplySelectionOperationC<INS_VEC, BOOL_OPER>(COND, testData, lhs, rhs);
} 


template< typename INS_VEC, typename BOOL_OPER>
Vec<INS_VEC> ApplySelection(BOOL_OPER& COND, const Vec<INS_VEC>& testData, typename InstructionTraits<INS_VEC>::FloatType trueVal, typename InstructionTraits<INS_VEC>::FloatType falseVal)
{
	return ApplySelectionOperationC<INS_VEC, BOOL_OPER > (COND, testData, trueVal, falseVal);
}

template< typename INS_VEC, typename BOOL_OPER, typename TRUE_OPER, typename FALSE_OPER>
Vec<INS_VEC> ApplySelectionF(BOOL_OPER& COND, const Vec<INS_VEC>& testData, TRUE_OPER& trueOper, FALSE_OPER& falseOper)
{
	return ApplySelectionOperationFunc<INS_VEC, BOOL_OPER, TRUE_OPER, FALSE_OPER >(COND, testData, trueOper, falseOper);
}


/*
 applies the tesFunc to the vector val, if the func returns true  it applies the trueLambda to the value otherwise it applies the falseLambda
 can be slow if true/false Lambda functions are not heavyweight
*/
template< typename INS_VEC, typename  BOOL_TEST_OP, typename  TRUE_LAMBDA, typename  FALSE_LAMBDA>
Vec<INS_VEC> ApplySplitCalculate( BOOL_TEST_OP& testFunc, const Vec<INS_VEC>& val, TRUE_LAMBDA& trueLambda, FALSE_LAMBDA& falseLambda)
{
	return splitConditionalCalculate(val, testFunc, trueLambda, falseLambda);
}


/*
 applies the tesFunc to the view val, if the func returns true  it applies the trueLambda to the value otherwise it applies the falseLambda
 can be slow if true/false Lambda functions are not heavyweight
*/
template< typename INS_VEC, typename  BOOL_TEST_OP, typename  TRUE_LAMBDA, typename  FALSE_LAMBDA>
VecView<INS_VEC> ApplySplitCalculate( BOOL_TEST_OP& testFunc, const VecView<INS_VEC>& val, TRUE_LAMBDA& trueLambda, FALSE_LAMBDA& falseLambda)
{
	return splitConditionalCalculate(val, testFunc, trueLambda, falseLambda);
}


/////////////////////////// filters //////////////////////


template<typename LAMBDA, typename INS_VEC>
void ApplyTransformM(LAMBDA& lambda, VecView<INS_VEC>& inputVec)
{
	return ApplyUnitaryOperation(lambda, inputVec);
}

template<typename LAMBDA, typename INS_VEC>
VecView<INS_VEC> ApplyTransform(LAMBDA& lambda, const VecView<INS_VEC>& inputVec)
{
	return ApplyUnitaryOperation(lambda, inputVec);
}


template<typename LAMBDA, typename INS_VEC>
VecView<INS_VEC> ApplyTransformV(LAMBDA& lambda, const Vec<INS_VEC>& inputVec)
{
	return ApplyUnitaryOperation(lambda, inputVec);
}


/*
applies the OP to the view in and  scatter,  writes the results to the corresponding elements  of the result vector.
*/
template< typename INS_VEC, typename OP>
void ApplyTransformWrite(OP& oper, const VecView<INS_VEC>& view, Vec<INS_VEC>& out)
{
	ApplyUnitaryOperationWrite(oper, view, out);
}
