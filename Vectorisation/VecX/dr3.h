#pragma once
#include "vec.h"
#include "binary_unitary_operations.h"
#include "conditional_select_eval.h"
#include "filter_pipe_and_join.h"
#include "filter_select.h"
#include "alloc_policy_imp.h"

#include "operations.h"
#include "apply_operation.h"

#include "vec.h"
#include "vec_d.h"
#include "vec_bool.h"
#include "vec_view.h"
#include "transform.h"

#include "target_name_space.h"


#ifdef _MSC_VER
	#define _VC_PERF_REG_
#else
	#undef _VC_PERF_REG_
#endif

//////////  transform //////////////

//Unitary lambdas
//unrolled version defaults to unroll x4
template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  transform(LAMBDA& lambda, const Vec<INS_VEC>& inputVec)
{

#ifdef _VC_PERF_REG_
	return  ApplyTransformUR_X(inputVec, lambda);
#else
	return ApplyUnitaryOperation(inputVec, lambda);
#endif

}


//x8 unrolled version 
template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  transformXX(LAMBDA& lambda, const Vec<INS_VEC>& inputVec)
{
	return  ApplyTransformUR_XX(inputVec, lambda);
}


template<typename LAMBDA, typename INS_VEC>
void  transform(LAMBDA& lambda, const Vec<INS_VEC>& inputVec,  Vec<INS_VEC>& outVec)
{
	return  ApplyUnitaryOperation(inputVec, outVec, lambda);
}




//not unrolled x1
template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  transform1(LAMBDA& lambda, const Vec<INS_VEC>& inputVec)
{
	return ApplyUnitaryOperation1(inputVec, lambda);
}


//inplace transforms
template<typename LAMBDA, typename INS_VEC>
void transformM(LAMBDA& lambda, Vec<INS_VEC>& inputVec)
{
	ApplyUnitaryOperationM(inputVec, lambda);
}


//// experimental inplace transforms
// Sampler loads multiple offset INS_VECs so that lagged and advanced values can be used by the lambda oper,
// oper takes the sampler as an input argument, iteration and transformation of input vector   
// is confined to startPos and endPos, subject to not going out of range.
template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP, typename SAMPLER >
void transform(VEC_TYPE<INS_VEC>& inputVec, VEC_TYPE<INS_VEC>& result, OP& oper, SAMPLER& sampler, int startPos = 0, int endPos = -1)
{
	ApplyTransformUR_X_Impl_EX(inputVec, result, oper,  sampler, startPos, endPos);
}



//// experimental inplace transforms binary transform version of above
// Sampler loads multiple offset INS_VECs so that lagged and advanced values can be used by the lambda oper,
// oper takes the sampler as an input argument, iteration and transformation of input vector   
// is confined to startPos and endPos, subject to not going out of range.
template<  template <class> typename VEC_TYPE, typename INS_VEC, typename OP, typename SAMPLER >
void transform(const VEC_TYPE<INS_VEC>& inputVec, const VEC_TYPE<INS_VEC>& aux_inputVec, VEC_TYPE<INS_VEC>& result, OP& oper, SAMPLER& sampler, int startPos = 0, int endPos = -1)
{
	ApplyTransformUR_X_Impl_EX(inputVec, aux_inputVec, result, oper, sampler, startPos, endPos);
}



//not unrolled x1
template<typename LAMBDA, typename INS_VEC>
void  transform1(LAMBDA& lambda, Vec<INS_VEC>& inputVec)
{
	inputVec =ApplyUnitaryOperation1(inputVec, lambda);
}


//binary lambdas

//Unitary lambdas
//unrolled version defaults to unroll x4
template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  transform(LAMBDA& lambda, const Vec<INS_VEC>& inputVecLHS, const Vec<INS_VEC>& inputVecRHS)
{

#if defined (_VC_PERF_REG_)
	return ApplyBinaryTransformUR_X(inputVecLHS, inputVecRHS, lambda);
#else
	return ApplyBinaryOperation(inputVecLHS, inputVecRHS, lambda);
#endif
	
}

template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  transform(LAMBDA& lambda, typename InstructionTraits<INS_VEC>::FloatType LHS, const Vec<INS_VEC>& inputVecRHS)
{
#ifdef _VC_PERF_REG_
	return ApplyBinaryTransformUR_X(LHS, inputVecRHS, lambda);
#else
	return ApplyBinaryOperation(LHS, inputVecRHS, lambda);
#endif

}

template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  transform(LAMBDA& lambda, const Vec<INS_VEC>& inputVecLHS, typename InstructionTraits<INS_VEC>::FloatType RHS)
{
#ifdef _VC_PERF_REG_
	return ApplyBinaryTransformUR_X(inputVecLHS, RHS, lambda);
#else
	return ApplyBinaryOperation(inputVecLHS, RHS, lambda);
#endif

}

// not unrolled
template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  transform1(LAMBDA& lambda, const Vec<INS_VEC>& inputVecLHS, const Vec<INS_VEC>& inputVecRHS)
{
	return ApplyBinaryOperation1(inputVecLHS, inputVecRHS, lambda);
}

//conversion 
template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  transform1(LAMBDA& lambda, typename InstructionTraits<INS_VEC>::FloatType LHS, const Vec<INS_VEC>& inputVecRHS)
{
	return ApplyBinaryOperation1(LHS, inputVecRHS, lambda);
}

template<typename LAMBDA, typename INS_VEC>
Vec<INS_VEC>  transform1(LAMBDA& lambda, const Vec<INS_VEC>& inputVecLHS, typename InstructionTraits<INS_VEC>::FloatType RHS)
{
	return ApplyBinaryOperation1(inputVecLHS, RHS, lambda);
}

/////////////////////////////// in place binary unroll x 4 ////////////////

template<typename LAMBDA, typename INS_VEC>
void transformM(LAMBDA& lambda, Vec<INS_VEC>& inputVecLHS, const Vec<INS_VEC>& inputVecRHS)
{
	ApplyBinaryOperationMMXY<INS_VEC, LAMBDA>(inputVecLHS, inputVecRHS, lambda);
}

template<typename LAMBDA, typename INS_VEC>
void transformM(LAMBDA& lambda, Vec<INS_VEC>& inputVecLHS, typename InstructionTraits<INS_VEC>::FloatType RHS)
{
	ApplyBinaryOperationMMMX<INS_VEC, LAMBDA>(inputVecLHS, RHS, lambda);
}


template<typename LAMBDA, typename INS_VEC>
void transformM(LAMBDA& lambda, typename InstructionTraits<INS_VEC>::FloatType LHS, Vec<INS_VEC>& inputVecRHS)
{
	ApplyBinaryOperationMMM<INS_VEC, LAMBDA>(LHS, inputVecRHS, lambda);
}





///////////////////////////////////////////////////////////////////////
// selection operations generate a vector (or view)  from an input view and a conditional. 
// At their simplest they take a boolean vector and a select between  values in the case of true or false values
// this can become more complex by having a lambda determine the true or false values in situ
// and by calculating the true or false value in situ

template< typename INS_VEC>
Vec<INS_VEC> select(const VecBool<INS_VEC>& condition, const Vec<INS_VEC>& trueVals, const Vec<INS_VEC>& falseVals)
{
	return ApplySelectionOperation(condition, trueVals, falseVals);
}


template< typename INS_VEC>
Vec<INS_VEC> select(const VecBool<INS_VEC>& condition, typename InstructionTraits<INS_VEC>::FloatType trueVal, typename InstructionTraits<INS_VEC>::FloatType falseVal)
{
	return ApplySelectionOperation(condition, trueVal, falseVal);

}

template< typename INS_VEC>
Vec<INS_VEC> select(const VecBool<INS_VEC>& condition, typename InstructionTraits<INS_VEC>::FloatType trueVal, const Vec<INS_VEC>& falseVals)
{
	return ApplySelectionOperation(condition, trueVal, falseVals);

}

template< typename INS_VEC>
Vec<INS_VEC> select(const VecBool<INS_VEC>& condition, const Vec<INS_VEC>& trueVals, typename InstructionTraits<INS_VEC>::FloatType falseVal)
{
	return ApplySelectionOperation(condition, trueVals, falseVal);
}



template< typename INS_VEC, typename BOOL_OPER>
Vec<INS_VEC> select(BOOL_OPER& COND, const Vec<INS_VEC>& testData, const Vec<INS_VEC>& trueVals, const Vec<INS_VEC>& falseVals)
{
	return ApplySelectionOperationC<INS_VEC, BOOL_OPER>(COND, testData, trueVals, falseVals);
}


template< typename INS_VEC, typename BOOL_OPER>
Vec<INS_VEC> select(BOOL_OPER& COND, const Vec<INS_VEC>& testData, typename InstructionTraits<INS_VEC>::FloatType trueVal, typename InstructionTraits<INS_VEC>::FloatType falseVal)
{
#ifdef _VC_PERF_REG_
	return  ApplySelectTransformUR_XC(testData, COND, trueVal, falseVal);
#else
	return ApplySelectionOperationC<INS_VEC, BOOL_OPER >(COND, testData, trueVal, falseVal);
#endif
}


template< typename INS_VEC, typename BOOL_OPER, typename TRUE_OPER, typename FALSE_OPER>
Vec<INS_VEC> selectTransform(BOOL_OPER& COND, const Vec<INS_VEC>& testData, TRUE_OPER& trueOper, FALSE_OPER& falseOper)
{

#ifdef _VC_PERF_REG_
	return ApplySelectionOperationFuncUR_X<INS_VEC, BOOL_OPER, TRUE_OPER, FALSE_OPER >(COND, testData, trueOper, falseOper);
#else
	return ApplySelectionOperationFunc<INS_VEC, BOOL_OPER, TRUE_OPER, FALSE_OPER >(COND, testData, trueOper, falseOper);
#endif

}


/*
 applies the tesFunc to the vector val, if the func returns true  it applies the trueLambda to the value otherwise it applies the falseLambda
 can be slow if true/false Lambda functions are not heavyweight
*/
template< typename INS_VEC, typename  BOOL_TEST_OP, typename  TRUE_LAMBDA, typename  FALSE_LAMBDA>
Vec<INS_VEC> filterTransform(BOOL_TEST_OP& testFunc, const Vec<INS_VEC>& val, TRUE_LAMBDA& trueLambda, FALSE_LAMBDA& falseLambda)
{
#ifdef _VC_PERF_REG_
	return splitConditionalCalculate_X(val, testFunc, trueLambda, falseLambda);
#else
	return splitConditionalCalculate(val, testFunc, trueLambda, falseLambda);
#endif
}


/*
 applies the tesFunc to the view val, if the func returns true  it applies the trueLambda to the value otherwise it applies the falseLambda
 can be slow if true/false Lambda functions are not heavyweight
*/
template< typename INS_VEC, typename  BOOL_TEST_OP, typename  TRUE_LAMBDA, typename  FALSE_LAMBDA>
VecView<INS_VEC> filterTransform(BOOL_TEST_OP& testFunc, const VecView<INS_VEC>& val, TRUE_LAMBDA& trueLambda, FALSE_LAMBDA& falseLambda)
{
	return splitConditionalCalculate(val, testFunc, trueLambda, falseLambda);
}


///////////////////////////  transforming  views values //////////////////////
//why is this needed ??
// so that we get a view not a vector created
template<typename LAMBDA, typename INS_VEC>
VecView<INS_VEC> transformV(LAMBDA& lambda, const Vec<INS_VEC>& inputVec)
{
	return ApplyUnitaryOperationV(inputVec,lambda);
}


/*
Modify  transforms the values held in the VecView object with the lambda. Has unrolled bersion for VC++
performance regression
in place modify avoids the cost of setting up indexes again in a separate view object
*/
template<typename LAMBDA, typename INS_VEC>
void transformM(LAMBDA& lambda, VecView<INS_VEC>& inputVec)
{

#ifdef _VC_PERF_REG_
	ApplyTransformUR_X(inputVec, lambda);
#else
	ApplyUnitaryOperation(inputVec,lambda);
#endif

}

/*
creates a new view indexed onto the same source vector but has values whicvh are transforms the values held in 
input VecView object with the lambda. Has unrolled version for VC++  performance regression
*/
template<typename LAMBDA, typename INS_VEC>
VecView<INS_VEC> transform(LAMBDA& lambda, const VecView<INS_VEC>& inputVec)
{
#ifdef _VC_PERF_REG_
	return ApplyTransformUR_X(inputVec, lambda);
#else
	return ApplyUnitaryOperation( inputVec, lambda);
#endif

}



//input should be const
template<typename LAMBDA, typename INS_VEC>
void transform(LAMBDA& lambda, StridedSpan<INS_VEC>& inputVec, Span<INS_VEC>& outVec)
{
		
	StridedSampler<INS_VEC> strided_sampler(inputVec.stride());

	auto wrappedLambda = [&](StridedSampler<INS_VEC>& sampler)
	{
		auto x = sampler.X_0.value; 
		return lambda(x);
	};

	ApplyTransformUR_X_Impl_EX_STRD(inputVec, outVec, wrappedLambda, strided_sampler, 0, int(inputVec.paddedSize()));

}


//input should be const
template<typename LAMBDA, typename INS_VEC>
void transform(LAMBDA& lambda,   Span<INS_VEC>& inputVec, Span<INS_VEC>& outVec)
{
	UnitarySampler<INS_VEC> identity_sampler;

	auto wrappedLambda = [&](UnitarySampler<INS_VEC>& sampler)
	{
		auto x = sampler.X_0.value;
		return lambda(x);
	};

	ApplyTransformUR_X_Impl_EX(inputVec, outVec, wrappedLambda, identity_sampler, 0, int(inputVec.paddedSize()));

}


////////
template<typename LAMBDA, typename INS_VEC>
void transformM(LAMBDA& lambda, Span<INS_VEC>& inputOutputSpan)
{
	transform(lambda, inputOutputSpan, inputOutputSpan);
}





/*
applies the OP to the view in and  scatter,  writes the results to the corresponding elements  of the result vector.
*/
template< typename INS_VEC, typename OP>
void transformWrite(OP& oper, const VecView<INS_VEC>& view, Vec<INS_VEC>& out)
{
	ApplyUnitaryOperationWrite(oper, view, out);
}


///// filters ///

/*
returns a view containing eleements from the vector which have corresponding true values in the VectorBool condition input
*/
template< typename INS_VEC>
VecView<INS_VEC>  filterB(const VecBool<INS_VEC>& condition, const Vec<INS_VEC>& lhs)
{
	return ApplyFilterBImpl< Vec, INS_VEC>(condition, lhs);
};

/*
returns a view containing eleements from the vector which have corresponding true values in the VectorBool condition input
*/
template< typename INS_VEC>
VecView<INS_VEC>  filterB(const VecBool<INS_VEC>& condition, const VecView<INS_VEC>& lhs)
{
	return ApplyFilterBImpl< VecView, INS_VEC>(condition, lhs);
}


/*
applys filter  boloean lammda  (condition) to a  values held in VectorView  if result is true values
are passed on to the result vvector view which is returned.
*/
template< typename INS_VEC, typename OP>
VecView<INS_VEC>  filter(OP& condition, const VecView<INS_VEC>& lhs)
{
	return ApplyFilterImpl< VecView, INS_VEC, OP>(condition, lhs);
};

/*
* tests if any elements satisfy the condition before trying to do the conditional element copy
*/
template< typename INS_VEC, typename OP>
VecView<INS_VEC>  filter(OP& condition, const Vec<INS_VEC>& lhs)
{
	return ApplyFilterImpl< Vec, INS_VEC, OP>(condition, lhs);
};


template< typename INS_VEC, typename OP>
VecView<INS_VEC> filter(OP& condition, const Span<INS_VEC>& lhs)
{
	return ApplyFilter(condition, lhs);
}


template< typename INS_VEC, typename OP>
VecView<INS_VEC> filter(OP& condition, const StridedSpan<INS_VEC>& lhs)
{

	StridedSampler<INS_VEC> strided_sampler(lhs.stride());
	auto wrappedLambda = [&](StridedSampler<INS_VEC>& sampler)
	{
		auto x = sampler.X_0.value;
		return condition(x);
	};

	assert(lhs.stride() % InstructionTraits<INS_VEC>::width == 0);

	return ApplyFilterImpl_EXt_STRD(wrappedLambda, lhs, strided_sampler);
}


/*
returns a view of the first N elements of the vector that satisfy the condition
*/
template< typename INS_VEC, typename OP>
VecView<INS_VEC>  countedFilter(OP& condition, const Vec<INS_VEC>& lhs, int N)
{
	return CountedFilterImpl(condition, lhs, N);
};

/*
returns a view of the first N elements of the view that satisfy the condition
*/
template< typename INS_VEC, typename OP>
VecView<INS_VEC>  countedFilter(OP& condition, const VecView<INS_VEC>& lhs, int N)
{
	return CountedFilterImpl(condition, lhs, N);
};


/*
returns a tuple of views on the input vector, the first view is for elements satisfying the condition, the second for remaining elements
*/
template< typename INS_VEC, typename OP>
std::tuple<VecView<INS_VEC>, VecView<INS_VEC> >  binaryFilter(OP& condition, const Vec<INS_VEC>& lhs)
{
	return  BinaryFilterImpl(condition, lhs);
};

/*
returns a tuple of views on the input view, the first view is for elements satisfying the condition, the second for remaining elements
*/
template< typename INS_VEC, typename OP>
std::tuple<VecView<INS_VEC>, VecView<INS_VEC> >  binaryFilter(OP& condition, const VecView<INS_VEC>& lhs)
{
	return  BinaryFilterImpl(condition, lhs);
};




//    ApplySparseTransform  takes a boolean condition lambda to determine if the lambda should be used
//    to calculate a given value in the vector
//    if so it applies the transform lambda and then blends it into the target vector 
template< typename LAMBDA, typename INS_VEC, typename CONDITION_LAMBDA>
void sparseTransform(const Vec<INS_VEC>& inputVec, Vec<INS_VEC>& updateResult, LAMBDA& oper, CONDITION_LAMBDA& selectionOp)
{
	ApplySparseUnitaryOperationU(inputVec, updateResult, oper, selectionOp);
}


template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType reduce(OP& oper, const Vec<INS_VEC>& rhs1, typename InstructionTraits<INS_VEC>::FloatType initVal, bool singularInit = true)
{
	return ApplyAccumulate2(rhs1, oper, initVal, singularInit);
}



/*
not unrolled version of reuction applying oper over vector rhs1
*/
template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType reduce1(const Vec<INS_VEC>& rhs1, OP& oper, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue,  bool singularInit = true)
{
	return ApplyAccumulate2(rhs1,oper, initVal,  singularInit);
}


//unroll version
template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType reduce(const Vec<INS_VEC>& rhs1, OP& oper, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue,  bool singularInit = true)
{
	ignore(initVal);
	ignore(singularInit);	
#ifdef _VC_PERF_REG_
	return ApplyAccumulate2UR_X(rhs1, oper);


#else
	return ApplyAccumulate2UR(rhs1, oper, initVal, singularInit);
#endif
	
}


//unroll version
template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType reduce(const Span<INS_VEC>& rhs1, OP& oper, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{
	return ApplyAccumulate2UR_X(rhs1, oper);
	ignore(initVal);
	ignore(singularInit);
}

//unroll version
template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType reduce(const VecView<INS_VEC>& rhs1, OP& oper, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{

	return ApplyAccumulate2UR_X(rhs1, oper);
	ignore(initVal);
	ignore(singularInit);

}


//////experimental unrolled  double accumulation
template< typename INS_VEC, typename OP, typename OP1>
typename std::tuple<typename InstructionTraits<INS_VEC>::FloatType, typename InstructionTraits<INS_VEC>::FloatType>
reduceM(const Vec<INS_VEC>& rhs, OP& oper, OP1& oper1)
{
	return ApplyAccumulate2UR_X2(rhs, oper, oper1);

}




//unitary transform
template< typename INS_VEC, typename OPT, typename OP>
typename InstructionTraits<INS_VEC>::FloatType transformReduce1(const Vec<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{
	return ApplyTransformAccumulate(rhs1, operTransform, operAcc, initVal, singularInit);
}



template< typename INS_VEC, typename OP, typename OPT>
typename InstructionTraits<INS_VEC>::FloatType transformReduce(OPT& operTransform, OP& operAcc, const Vec<INS_VEC>& rhs1, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{
	return ApplyTransformAccumulateUR(rhs1, operTransform, operAcc, initVal, singularInit);
}




//unitary transform unrolled
template< typename INS_VEC, typename OP, typename OPT>
typename InstructionTraits<INS_VEC>::FloatType transformReduce(const Vec<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{
	ignore(initVal);
	ignore(singularInit);
	return ApplyTransformAccumulate2UR_X(rhs1, operTransform, operAcc);
/*
#ifdef _VC_PERF_REG_
	return ApplyTransformAccumulate2UR_X(rhs1, operTransform, operAcc);
	ignore(initVal);
	ignore(singularInit);
#else
	return ApplyTransformAccumulateUR(rhs1, operTransform, operAcc, initVal, singularInit);
#endif
	*/
}


//unitary transform unrolled
template< typename INS_VEC, typename OP, typename OPT>
typename InstructionTraits<INS_VEC>::FloatType transformReduce(const VecView<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{
	return ApplyTransformAccumulate2UR_X(rhs1, operTransform, operAcc);
}


template< typename INS_VEC, typename OP, typename OPT>
typename InstructionTraits<INS_VEC>::FloatType transformReduce(const Span<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{
	return ApplyTransformAccumulate2UR_X_Impl(rhs1, operTransform, operAcc);
}



//////experimental unrolled  double transform accumulation
template< typename INS_VEC, typename TF1, typename RED1,  typename TF2,  typename RED2>
typename std::tuple<typename InstructionTraits<INS_VEC>::FloatType, typename InstructionTraits<INS_VEC>::FloatType>
transformReduceM(const Vec<INS_VEC>& rhs, TF1& transform1, RED1& reduce1, TF2& transform2, RED2& reduce2)
{
	return ApplyTransformAccumulate2UR_X2(rhs, transform1, reduce1, transform2, reduce2);
}




//binary transform  unrolled
template< typename INS_VEC, typename OP, typename OPT>
typename InstructionTraits<INS_VEC>::FloatType transformReduce(const Vec<INS_VEC>& lhs1, const Vec<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{

#if defined( _VC_PERF_REG_)
	return ApplyTransformAccumulate2UR_XBin(lhs1, rhs1, operTransform, operAcc);
	ignore(initVal);
	ignore(singularInit);
#else
	return ApplyTransformAccumulateUR(lhs1, rhs1, operTransform, operAcc, initVal, singularInit);
#endif
	
}


//binary transform  unrolled
template< typename INS_VEC, typename OP, typename OPT>
typename InstructionTraits<INS_VEC>::FloatType transformReduce(const VecView<INS_VEC>& lhs1, const VecView<INS_VEC>& rhs1, OPT& operTransform, OP& operAcc, typename InstructionTraits<INS_VEC>::FloatType initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
{

//#if defined( _VC_PERF_REG_)
	return ApplyTransformAccumulate2UR_XBin(lhs1, rhs1, operTransform, operAcc);
	ignore(initVal);
	ignore(singularInit);
//#else
	//return ApplyTransformAccumulateUR(lhs1, rhs1, operTransform, operAcc, initVal, singularInit);
//#endif

}


///////////////////////////adjacent diff ////////////

template< typename INS_VEC, typename OP>
Vec<INS_VEC> adjacentDiff(const Vec<INS_VEC>& rhs1, OP& oper)
{
	return  ApplyAdjacentDiff(rhs1, oper);
};


//
//inclusive scan
template< typename INS_VEC, typename OP>
Vec<INS_VEC> scan(const Vec<INS_VEC>& rhs1, OP& oper)
{
	return ApplyScan(rhs1, oper);
}
