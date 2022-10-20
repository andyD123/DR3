/****************************  binary_unitary_operations.h  *******************************
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
#include "vec_bool.h"
#include "instruction_traits.h"
#include "unroll_operators.h"
#include "apply_operation.h"
#include "error_utils.h"

template<typename LAMBDA>
struct CustomUnitary
{
	CustomUnitary(LAMBDA& lamb) :myLambda(lamb) {}
	
	template<typename INS_VEC>
	inline INS_VEC operator()(const INS_VEC& X)
	{
		return myLambda(X);
	}

	LAMBDA myLambda;
};


template<typename LAMBDA>
struct CustomUnitaryBool
{
	CustomUnitaryBool(LAMBDA& lamb) :myLambda(lamb) {}
	template<typename INS_VEC>
	inline typename InstructionTraits<INS_VEC>::BoolType operator()(const INS_VEC& X)
	{
		return myLambda(X);
	}
	LAMBDA myLambda;
};


template<typename  LAMBDA>
constexpr CustomUnitary< LAMBDA>  getLambda(LAMBDA& lda)
{
	CustomUnitary< LAMBDA> lam(lda);
	return lam;
}

template<typename  LAMBDA>
constexpr CustomUnitaryBool< LAMBDA>  getLambdaBool(LAMBDA& lda)
{
	CustomUnitaryBool< LAMBDA> lam(lda);
	return lam;
}

template<typename LAMBDA>
struct CustomBinary
{
	CustomBinary(LAMBDA& lamb) :myLambda(lamb) {}
	template<typename INS_VEC>
	inline INS_VEC operator()(const INS_VEC& X, const INS_VEC& Y)
	{
		return myLambda(X, Y);
	}
	LAMBDA myLambda;
};


template<typename LAMBDA>
struct CustomBinaryBool
{
	CustomBinaryBool(LAMBDA& lamb) :myLambda(lamb) {}
	template<typename INS_VEC>
	inline typename InstructionTraits<INS_VEC>::BoolType operator()(const INS_VEC& X, const INS_VEC& Y)
	{
		return myLambda(X, Y);
	}

	template<typename INS_VEC>
	inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& X, const INS_VEC& Y)
	{
		return myLambda(X,Y);
	}


	template<typename INS_VEC>
	inline typename InstructionTraits<INS_VEC>::BoolType apply(const typename InstructionTraits<INS_VEC>::FloatType& X, const INS_VEC& Y)
	{
		return myLambda(X, Y);
	}

	template<typename INS_VEC>
	inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& X, const typename InstructionTraits<INS_VEC>::FloatType& Y)
	{
		return myLambda(X, Y);
	}



	LAMBDA myLambda;
};


template<typename LAMBDA>
struct CustomSparseBinary
{
	CustomSparseBinary(LAMBDA& lamb) :myLambda(lamb) {}
	template<typename INS_VEC>
	inline INS_VEC operator()(const INS_VEC& X, const INS_VEC& Y)
	{
		return myLambda(X, Y);
	}
	LAMBDA myLambda;
};

template< typename INS_VEC, typename LAMBDA, typename LAMBDA_D>
VecD<INS_VEC> ApplyLambdaD(const VecD<INS_VEC>& in, LAMBDA& lamb, LAMBDA_D& lamb_D)
{
	check_vector(in);

	CustomUnitary<LAMBDA> oper(lamb);
	auto val = ApplyUnitaryOperation(in.value(), oper);
	CustomUnitary<LAMBDA_D> operD(lamb_D);
	auto deriv = ApplyUnitaryOperation(in.value(), operD);
	deriv *= in.derivative();
	return  VecD<INS_VEC>(std::move(val), std::move(deriv));
}

template< typename INS_VEC, typename LAMBDA>
Vec<INS_VEC> ApplyLambda(const Vec<INS_VEC>& in, LAMBDA& lamb)
{
	check_vector(in);
	CustomUnitary<LAMBDA> oper(lamb);
	return ApplyUnitaryOperation(in, oper);
}

template< typename INS_VEC, typename LAMBDA>
Vec<INS_VEC> ApplyLambda1(const Vec<INS_VEC>& in, LAMBDA& lamb)
{
	check_vector(in);
	CustomUnitary<LAMBDA> oper(lamb);
	return ApplyUnitaryOperation1(in, oper);
}

template< typename INS_VEC, typename LAMBDA>
Vec<INS_VEC> ApplyLambda2(const Vec<INS_VEC>& in1, const Vec<INS_VEC>& in2, LAMBDA& lamb)
{
	check_pair(in2, in1);
	CustomBinary<LAMBDA> oper(lamb);
	return ApplyBinaryOperation1(in1, in2, oper);
}

template< typename INS_VEC, typename LAMBDA>
VecBool<INS_VEC> ApplyBoolLambda2(const Vec<INS_VEC>& in1, const Vec<INS_VEC>& in2, LAMBDA& lamb)
{
	check_pair(in2, in1);
	CustomBinaryBool<LAMBDA> oper(lamb);
	return ApplyBoolBinaryOperation1(in1, in2, oper);
}


//////////////////////////sparse update lambdas /////////////////////////////////

template< typename INS_VEC, typename LAMBDA>
void SparseUpdateWithLambda1(Vec<INS_VEC>& out, const Vec<INS_VEC>& in, LAMBDA& lamb)
{
	check_pair(out, in);
	CustomSparseBinary<LAMBDA> oper(lamb);
	return ApplySparseUnitaryOperation(out, in, oper);
}

template< typename INS_VEC, typename LAMBDA>
void SparseUpdateWithLambdaT(Vec<INS_VEC>& out, const Vec<INS_VEC>& in, LAMBDA& lamb)
{
	check_pair(out, in);
	return ApplySparseUnitaryOperation(out, in, lamb);
}

//we apply the conditional, if doesnt apply to any we move on 
//otherwise we calculate and over write out
/*
we need to have two lambdas one for the condition and one for the evaluation
the single argument version has this embedded

the test lambda also drives the blend between old values and  to be updated values
see sparse update FMA version of inverse cum normal distribution
*/
template< typename INS_VEC, typename OP>
void ApplySparseUnitaryOperation(Vec<INS_VEC>& rhs1, const Vec<INS_VEC>& rhs2, OP& oper)
{
	check_pair(rhs1, rhs2);
	//assert equality of size ?
	if (rhs1.isScalar() && rhs2.isScalar())
	{
		INS_VEC rhsdat1(rhs1.getScalarValue());
		INS_VEC rhsdat2(rhs2.getScalarValue());
		Vec<INS_VEC> vecRes = oper(rhsdat1, rhsdat2)[0];
		rhs1 = vecRes;
		return;
	}

	Vec<INS_VEC> result(rhs1.size());
	auto pRes = result.start();
	auto pRhs1 = rhs1.start();
	auto pRhs2 = rhs2.start();
	int sz = rhs1.paddedSize();
	Unroll_Binary<INS_VEC, OP>::apply_1(sz, pRhs1, pRhs2, pRes, oper);
	rhs1 = std::move(result); 
}
////////////////////////////////////////////////////////////

template< typename INS_VEC, typename OP, typename BOOL_OP>
void ApplySparseUnitaryOperationU(const Vec<INS_VEC>& lhs1, Vec<INS_VEC>& rhs1, OP& oper, BOOL_OP& selectionOp)
{
	check_pair(lhs1, rhs1);
	//assert equality of size ?
	if (lhs1.isScalar() && rhs1.isScalar())
	{
		INS_VEC rhsdat1(lhs1.getScalarValue());
		INS_VEC rhsdat2(rhs1.getScalarValue());
		INS_VEC result = oper(rhsdat1);
		typename InstructionTraits<INS_VEC>::FloatType scalarResult = result[0]; //first element of register
		Vec<INS_VEC> vecRes = scalarResult; //construct scalar vector from unitary op on lhs1
		rhs1 = vecRes; //wrong i think const
		return;
	}

	auto plhs1 = lhs1.start();
	auto pRhs1 = rhs1.start();
	int sz = rhs1.paddedSize();
	Unroll_Binary<INS_VEC, OP>::apply_1(sz, plhs1, pRhs1,  oper, selectionOp);
}

///////////////////////////////////////  new unitary ops ////////////////////
template< typename INS_VEC, typename OP>
const Vec<INS_VEC>& ApplyUnitaryOperationM(Vec<INS_VEC>& lhs, OP& oper)
{
	check_vector(lhs);
	if (lhs.isScalar())
	{
		lhs =ApplyUnitaryOperationM<INS_VEC, OP>(lhs.getScalarValue(),oper);
		return lhs;
	}

	auto pLhs = lhs.start();
	int sz = lhs.paddedSize();

	Unroll_Unitary< INS_VEC, OP>::apply_4(sz, pLhs, pLhs, oper);
	return lhs;
}


template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyUnitaryOperationM(typename InstructionTraits< INS_VEC>::FloatType lhs, OP& oper)
{
	INS_VEC res = oper(INS_VEC(lhs));
	return Vec<INS_VEC>(res[0]);
}


////////////////////////////////mutating operators ///////////////////////////////////////////////
//operation M == mutators
//used for +=,*=  etc
//doesnt make too much sense for vec oper+= (float&,float) , would need to be allfloat so its a 
//built in type
template< typename INS_VEC, typename OP>
const Vec<INS_VEC>& ApplyBinaryOperationM(Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	check_pair(lhs,rhs);
	if (rhs.isScalar())
	{
		return ApplyBinaryOperationM<INS_VEC, OP>(lhs, rhs.getScalarValue());
	}

	//if (lhs.isScalar()) //makes no sense  doube+=vec ->double
	auto pLhs = lhs.start();
	auto pRhs = rhs.start();
	int sz = lhs.paddedSize();
	OP oper;

	Unroll_Binary<INS_VEC, OP>::apply_4(sz, pLhs, pRhs, pLhs, oper);
	return lhs;
}

template< typename INS_VEC, typename OP>
const Vec<INS_VEC>& ApplyBinaryOperationMMMX(Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs, OP& oper)
{
	check_vector(lhs);
	INS_VEC RHS(rhs);
	auto pLhs = lhs.start();

	int sz = lhs.paddedSize();

	if (lhs.isScalar())
	{
		INS_VEC LHS(lhs.getScalarValue());
		INS_VEC res = oper(LHS, RHS);
		// lhs is a scalar so has no vec data member so no asign
		lhs = res[0];
		return lhs;
	}
	Unroll_Binary<INS_VEC, OP>::apply_4(sz, pLhs, RHS, pLhs, oper);
	return lhs;
}


template< typename INS_VEC, typename OP>
const Vec<INS_VEC>& ApplyBinaryOperationMMXY(Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs ,OP& oper )
{
	check_pair(lhs, rhs);
	if (rhs.isScalar())
	{
		return ApplyBinaryOperationMMMX<INS_VEC, OP>(lhs, rhs.getScalarValue(),oper);
	}

	//if (lhs.isScalar()) //makes no sense  doube+=vec ->double
	auto pLhs = lhs.start();
	auto pRhs = rhs.start();
	int sz = lhs.paddedSize();
	

	Unroll_Binary<INS_VEC, OP>::apply_4(sz, pLhs, pRhs, pLhs, oper);
	return lhs;

}


template< typename INS_VEC, typename OP>
const Vec<INS_VEC>& ApplyBinaryOperationMMM(typename InstructionTraits<INS_VEC>::FloatType lhs, Vec<INS_VEC>& rhs, OP& oper)
{
	check_vector(rhs);
	auto pRhs = rhs.start();

	int sz = rhs.paddedSize();
	INS_VEC LHS(lhs);
	if (rhs.isScalar())
	{
		INS_VEC RHS(rhs.getScalarValue());
		INS_VEC res = oper(LHS, RHS);
		//we dont store since its a scalar
		rhs = res[0];
		return rhs;
	}

	Unroll_Binary<INS_VEC, OP>::apply_4(sz, LHS, pRhs, pRhs, oper);
	return rhs;
}

template< typename INS_VEC, typename OP>
const Vec<INS_VEC>& ApplyBinaryOperationM(Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	check_vector(lhs);
	INS_VEC RHS(rhs);
	auto pLhs = lhs.start();
	OP oper;
	int sz = lhs.paddedSize();

	if (lhs.isScalar())
	{
		INS_VEC LHS(lhs.getScalarValue());
		INS_VEC res = oper(LHS, RHS);
		// lhs is a scalar so has no vec data member so no asign
		lhs = res[0];
		return lhs;
	}
	Unroll_Binary<INS_VEC, OP>::apply_4(sz, pLhs, RHS, pLhs, oper);
	return lhs;
}

template< typename INS_VEC, typename OP>
const Vec<INS_VEC>& ApplyBinaryOperationM(typename InstructionTraits<INS_VEC>::FloatType lhs, Vec<INS_VEC>& rhs)
{
	check_vector(rhs);
	auto pRhs = rhs.start();
	OP oper;
	int sz = rhs.paddedSize();
	INS_VEC LHS(lhs);
	if (rhs.isScalar())
	{
		INS_VEC RHS(rhs.getScalarValue());
		INS_VEC res = oper(LHS, RHS);
		//we dont store since its a scalar
		rhs = res[0];
		return rhs;
	}

	Unroll_Binary<INS_VEC, OP>::apply_4(sz, LHS, pRhs, pRhs, oper);
	return rhs;
}


//////////////////////////  unroll 4  createw new result vector /////////////////// 
template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyBinaryOperation(typename InstructionTraits<INS_VEC>::FloatType  lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	INS_VEC LHS(lhs);
	INS_VEC RHS(rhs);
	OP oper;
	INS_VEC res = oper(LHS, RHS);
	typename InstructionTraits<INS_VEC>::FloatType scalarVal = res[0];
	return Vec<INS_VEC>(scalarVal);
}

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyBinaryOperation(typename InstructionTraits<INS_VEC>::FloatType  lhs, typename InstructionTraits<INS_VEC>::FloatType rhs, OP& oper)
{
	INS_VEC LHS(lhs);
	INS_VEC RHS(rhs);
	INS_VEC res = oper(LHS, RHS);
	typename InstructionTraits<INS_VEC>::FloatType scalarVal = res[0];
	return Vec<INS_VEC>(scalarVal);
}


template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyBinaryOperation(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	check_pair(lhs, rhs);

	if (rhs.isScalar())
	{
		return ApplyBinaryOperation<INS_VEC, OP>(lhs, rhs.getScalarValue());
	}

	if (lhs.isScalar())
	{
		return ApplyBinaryOperation<INS_VEC, OP>(lhs.getScalarValue(), rhs);
	}

	Vec<INS_VEC> result(rhs.size());
	OP oper;

	int sz = lhs.paddedSize();
	auto pRes = result.start();
	auto pLhs = lhs.start();
	auto pRhs = rhs.start();
	Unroll_Binary<INS_VEC, OP>::apply_4(sz, pLhs, pRhs, pRes, oper);
	return result;
}


template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyBinaryOperation(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs,  OP& oper)
{
	check_pair(lhs, rhs);

	if (rhs.isScalar())
	{
		return ApplyBinaryOperation<INS_VEC, OP>(lhs, rhs.getScalarValue(), oper);
	}

	if (lhs.isScalar())
	{
		return ApplyBinaryOperation<INS_VEC, OP>(lhs.getScalarValue(), rhs, oper);
	}

	Vec<INS_VEC> result(rhs.size());
	
	int sz = lhs.paddedSize();
	auto pRes = result.start();
	auto pLhs = lhs.start();
	auto pRhs = rhs.start();
	Unroll_Binary<INS_VEC, OP>::apply_4(sz, pLhs, pRhs, pRes, oper);
	return result;
}


template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyBinaryOperation(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	check_vector(lhs);
	if (lhs.isScalar())
	{
		return ApplyBinaryOperation<INS_VEC, OP>(lhs.getScalarValue(), rhs);
	}

	Vec<INS_VEC> result(lhs.size());
	INS_VEC RHS(rhs);
	auto pRes = result.start();
	auto pLhs = lhs.start();
	OP oper;
	int sz = lhs.paddedSize();
	Unroll_Binary<INS_VEC, OP>::apply_4(sz, pLhs, RHS, pRes, oper);
	return result;
}


template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyBinaryOperation(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs,OP& oper)
{
	check_vector(lhs);
	if (lhs.isScalar())
	{
		return ApplyBinaryOperation<INS_VEC, OP>(lhs.getScalarValue(), rhs, oper);
	}

	Vec<INS_VEC> result(lhs.size());
	INS_VEC RHS(rhs);
	auto pRes = result.start();
	auto pLhs = lhs.start();
	
	int sz = lhs.paddedSize();
	Unroll_Binary<INS_VEC, OP>::apply_4(sz, pLhs, RHS, pRes, oper);
	return result;
}



template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyBinaryOperation(typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	check_vector(rhs);
	if (rhs.isScalar())
	{
		return ApplyBinaryOperation<INS_VEC, OP>(lhs, rhs.getScalarValue());
	}
	Vec<INS_VEC> result(rhs.size());
	auto pRes = result.start();
	auto pRhs = rhs.start();
	INS_VEC LHS(lhs);
	OP oper;
	int sz = rhs.paddedSize();
	Unroll_Binary<INS_VEC, OP>::apply_4(sz, LHS, pRhs, pRes, oper);
	return result;
}




template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyBinaryOperation(typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs, OP& oper)
{
	check_vector(rhs);
	if (rhs.isScalar())
	{
		return ApplyBinaryOperation<INS_VEC, OP>(lhs, rhs.getScalarValue(),oper);
	}
	Vec<INS_VEC> result(rhs.size());
	auto pRes = result.start();
	auto pRhs = rhs.start();
	INS_VEC LHS(lhs);
	int sz = rhs.paddedSize();
	Unroll_Binary<INS_VEC, OP>::apply_4(sz, LHS, pRhs, pRes, oper);
	return result;
}



/////////////////////   no unroll //////////////////

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyBinaryOperation1(const Vec<INS_VEC>& rhs1, const Vec<INS_VEC>& rhs2, OP& oper)
{
	check_pair(rhs1, rhs2);
	//assert equality of size ?
	if (rhs1.isScalar() )
	{
		return ApplyBinaryOperation1<INS_VEC, OP>(rhs1.getScalarValue(), rhs2, oper);
	}
	if (rhs2.isScalar())
	{
		return ApplyBinaryOperation1<INS_VEC, OP>(rhs1, rhs2.getScalarValue(), oper);
	}

	Vec<INS_VEC> result(rhs1.size());
	auto pRes = result.start();
	auto pRhs1 = rhs1.start();
	auto pRhs2 = rhs2.start();
	int sz = rhs1.paddedSize();

	Unroll_Binary<INS_VEC, OP>::apply_1(sz, pRhs1, pRhs2, pRes, oper);
	return result;
}


template< typename INS_VEC, typename OP>
VecBool<INS_VEC> ApplyBoolBinaryOperation1(const Vec<INS_VEC>& rhs1, const Vec<INS_VEC>& rhs2, OP& oper)
{
	check_pair(rhs1, rhs2);
	//assert equality of size ?
	if (rhs1.isScalar())
	{
		return ApplyBoolBinaryOperation1<INS_VEC, OP>(rhs1.getScalarValue(), rhs2, oper);
	}
	if (rhs2.isScalar())
	{
		return ApplyBoolBinaryOperation1<INS_VEC, OP>(rhs1, rhs2.getScalarValue(), oper);
	}

	VecBool<INS_VEC> result(rhs1.size());
	auto pRes = result.start();
	auto pRhs1 = rhs1.start();
	auto pRhs2 = rhs2.start();
	int sz = rhs1.paddedSize();

	BinaryBoolNumericUnroll<INS_VEC, OP>::apply_4(sz, pRhs1, pRhs2, pRes, oper);
	return result;
}


template< typename INS_VEC, typename OP>
VecBool<INS_VEC> ApplyBoolBinaryOperation1(typename InstructionTraits<INS_VEC>::FloatType lhs, typename InstructionTraits<INS_VEC>::FloatType rhs, OP& oper)
{
	INS_VEC LHS(lhs);
	INS_VEC RHS(rhs);

	auto res = oper(LHS, RHS);
	bool val = res[0];
	return VecBool<INS_VEC>(val);
}



template< typename INS_VEC, typename OP>
VecBool<INS_VEC> ApplyBoolBinaryOperation1(typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs, OP& oper)
{
	//check_pair(lhs, rhs);
	//assert equality of size ?

	if (rhs.isScalar())
	{
		return ApplyBoolBinaryOperation1<INS_VEC, OP>(lhs, rhs.getScalarValue(), oper);
	}

	VecBool<INS_VEC> result(rhs.size());
	auto pRes = result.start();
	INS_VEC LHS(lhs);
	auto pRhs = rhs.start();
	int sz = rhs.paddedSize();

	BinaryBoolNumericUnroll<INS_VEC, OP>::apply_4(sz, LHS, pRhs, pRes, oper);
	return result;
}


template< typename INS_VEC, typename OP>
VecBool<INS_VEC> ApplyBoolBinaryOperation1( const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs, OP& oper)
{
	//check_pair(lhs, rhs);
	//assert equality of size ?

	if (lhs.isScalar())
	{
		return ApplyBoolBinaryOperation1<INS_VEC, OP>(lhs.getScalarValue(), rhs, oper);
	}

	VecBool<INS_VEC> result(lhs.size());
	auto pRes = result.start();
	INS_VEC RHS(rhs);
	auto pLhs = lhs.start();
	int sz = lhs.paddedSize();

	BinaryBoolNumericUnroll<INS_VEC, OP>::apply_4(sz, pLhs, RHS, pRes, oper);
	return result;
}




template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyBinaryOperation1(typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs, OP& oper)
{
	check_vector(rhs);
	if (rhs.isScalar())
	{
		return ApplyBinaryOperation<INS_VEC, OP>(lhs, rhs.getScalarValue(), oper);
	}
	Vec<INS_VEC> result(rhs.size());
	auto pRes = result.start();
	auto pRhs = rhs.start();
	INS_VEC LHS(lhs);
	int sz = rhs.paddedSize();
	Unroll_Binary<INS_VEC, OP>::apply_1(sz, LHS, pRhs, pRes, oper);
	return result;
}

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyBinaryOperation1(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs, OP& oper)
{
	check_vector(lhs);
	if (lhs.isScalar())
	{
		return ApplyBinaryOperation<INS_VEC, OP>(lhs.getScalarValue(), rhs, oper);
	}

	Vec<INS_VEC> result(lhs.size());
	INS_VEC RHS(rhs);
	auto pRes = result.start();
	auto pLhs = lhs.start();

	int sz = lhs.paddedSize();
	Unroll_Binary<INS_VEC, OP>::apply_1(sz, pLhs, RHS, pRes, oper);
	return result;
}


template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType
ApplyBinaryOperation1(const typename InstructionTraits<INS_VEC>::FloatType& lhs, const typename InstructionTraits<INS_VEC>::FloatType& rhs, OP& oper)
{
	INS_VEC RHS(rhs);
	INS_VEC LHS(lhs);
	INS_VEC res = oper(lhs, rhs);
	return res[0];
}

///////////////// unitary operations /////////////

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyUnitaryOperation(typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	INS_VEC RHS(rhs);
	OP oper;
	INS_VEC res = oper(RHS);
	return Vec<INS_VEC>(res[0]);
}

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyUnitaryOperation(typename InstructionTraits<INS_VEC>::FloatType rhs, OP& oper)
{
	INS_VEC RHS(rhs);
	INS_VEC res = oper(RHS);
	return Vec<INS_VEC>(res[0]);
}

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyUnitaryOperation(const Vec<INS_VEC>& rhs, OP& oper)
{
	check_vector(rhs);
	if (rhs.isScalar())
	{
		return ApplyUnitaryOperation<INS_VEC, OP>(rhs.getScalarValue(), oper);
	}

	Vec<INS_VEC> result(rhs.size());
	auto pRes = result.start();
	auto pRhs = rhs.start();
	int sz = rhs.paddedSize();

	Unroll_Unitary<INS_VEC, OP>::apply_4(sz, pRhs, pRes, oper);
	return result;
}


template< typename INS_VEC, typename OP>
void ApplyUnitaryOperation(const Vec<INS_VEC>& rhs, Vec<INS_VEC>& res ,OP& oper)
{
	assert(res.size() == rhs.size());
	if (res.size() != rhs.size() || res.isScalar() != rhs.isScalar())
	{
		throw std::exception("ApplyUnitaryOperation ::size mismatch");
	}

	check_vector(res);
	check_vector(rhs);
	if (rhs.isScalar())
	{
		res = ApplyUnitaryOperation<INS_VEC, OP>(rhs.getScalarValue(), oper);
		return;
	}

	auto pRes = res.start();
	auto pRhs = rhs.start();
	int sz = rhs.paddedSize();

	Unroll_Unitary<INS_VEC, OP>::apply_4(sz, pRhs, pRes, oper);

}






template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyUnitaryOperation1(const Vec<INS_VEC>& rhs, OP& oper)
{
	check_vector(rhs);
	if (rhs.isScalar())
	{
		return ApplyUnitaryOperation<INS_VEC, OP>(rhs.getScalarValue(), oper);
	}
	Vec<INS_VEC> result(rhs.size());
	auto pRes = result.start();
	auto pRhs = rhs.start();
	int sz = rhs.paddedSize();
	Unroll_Unitary<INS_VEC, OP>::apply_1(sz, pRhs, pRes, oper);
	return result;
}


template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyUnitaryOperation(const Vec<INS_VEC>& rhs)
{
	check_vector(rhs);
	if (rhs.isScalar())
	{
		return ApplyUnitaryOperation<INS_VEC, OP>(rhs.getScalarValue());
	}
	Vec<INS_VEC> result(rhs.size());
	auto pRes = result.start();
	auto pRhs = rhs.start();
	OP oper;
	int sz = rhs.paddedSize();
	Unroll_Unitary<INS_VEC,OP>::apply_4(sz, pRhs, pRes, oper);
	return result;
}



///////////////  ternary operations  eg FMA ////////////
template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyTernaryOperation(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& mid, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	check_pair(lhs, mid);
	if (isScalar(lhs)  )
	{
		return ApplyTernaryOperation<INS_VEC, OP>(lhs.getScalarValue(), mid, rhs);
	}

	Vec<INS_VEC> result(lhs.size());
	auto pRes = result.start();
	auto pLhs = lhs.start();
	auto pMid = mid.start();
	OP oper;
	int sz = lhs.paddedSize();
	INS_VEC RHS(rhs);
	Unroll_Ternary<INS_VEC, OP>::apply_4(sz, pLhs, pMid, RHS, pRes,oper);
	return result;
}

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyTernaryOperation(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType mid, const Vec<INS_VEC>& rhs)
{
	check_pair(rhs, lhs);
	if ( isScalar(lhs))
	{
		return ApplyTernaryOperation<INS_VEC, OP>(lhs.getScalarValue(), mid, rhs);
	}
	Vec<INS_VEC> result(lhs.size());
	auto pRes = result.start();
	auto pLhs = lhs.start();
	INS_VEC MID(mid);
	auto pRhs = rhs.start();
	OP oper;
	int sz = lhs.paddedSize();
	Unroll_Ternary<INS_VEC, OP>::apply_4(sz, pLhs, MID, pRhs, pRes, oper);
	return result;
}

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyTernaryOperation(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& mid, const Vec<INS_VEC>& rhs)
{
	check_pair(lhs, mid);
	check_pair(rhs, mid);


	if (isScalar(lhs))
	{
		return ApplyTernaryOperation<INS_VEC, OP>(lhs.getScalarValue(), mid, rhs);
    }
	if (isScalar(mid))
	{
		return ApplyTernaryOperation<INS_VEC, OP>(lhs, mid.getScalarValue(), rhs);
	}
	if (isScalar(rhs))
	{
		return ApplyTernaryOperation<INS_VEC, OP>(lhs, mid, rhs.getScalarValue());
	}

	Vec<INS_VEC> result(lhs.size());
	auto pRes = result.start();
	auto pLhs = lhs.start();
	auto pMid = mid.start();
	auto pRhs = rhs.start();
	OP oper;
	int sz = lhs.paddedSize();
	Unroll_Ternary<INS_VEC, OP>::apply_4(sz, pLhs, pMid, pRhs, pRes, oper);
	return result;
}

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyTernaryOperation(typename InstructionTraits<INS_VEC>::FloatType lhs, typename InstructionTraits<INS_VEC>::FloatType mid, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	INS_VEC LHS(lhs);
	INS_VEC MID(mid);
	INS_VEC RHS(rhs);
	OP oper;
	INS_VEC RES = oper(LHS, MID, RHS);
	//scalar constructor of vec
	return Vec<INS_VEC>(RES[0]);
}

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyTernaryOperation(typename InstructionTraits<INS_VEC>::FloatType lhs,  typename InstructionTraits<INS_VEC>::FloatType mid , const Vec<INS_VEC>& rhs )
{
	check_vector(rhs);

	if (isScalar(rhs) )
	{
		return ApplyTernaryOperation<INS_VEC, OP>(lhs, mid, rhs.getScalarValue());
	}
	Vec<INS_VEC> result(rhs.size());
	auto pRes = result.start();
	auto pRhs = rhs.start();
	OP oper;
	int sz = rhs.paddedSize();
	INS_VEC LHS(lhs);
	INS_VEC MID(mid);
	Unroll_Ternary<INS_VEC, OP>::apply_4(sz, LHS, MID, pRhs, pRes, oper);
	return result;
}

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyTernaryOperation(typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& mid, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	check_vector(mid);

	if ( isScalar(mid) )
	{
		return ApplyTernaryOperation<INS_VEC, OP>(lhs, mid.getScalarValue(), rhs);
	}

	Vec<INS_VEC> result(mid.size());
	auto pRes = result.start();
	auto pMid = mid.start();
	OP oper;
	int sz = mid.paddedSize();
	INS_VEC LHS(lhs);
	INS_VEC RHS(rhs);
	Unroll_Ternary<INS_VEC, OP>::apply_4(sz, LHS, pMid, RHS, pRes, oper);
	return result;
}

template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyTernaryOperation(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType mid, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	check_vector(lhs);

	if (lhs.isScalar() )
	{
		return ApplyTernaryOperation<INS_VEC, OP>(lhs.getScalarValue(), mid, rhs);
	}

	Vec<INS_VEC> result(lhs.size());
	auto pRes = result.start();
	auto pLhs = lhs.start();
	OP oper;
	int sz = lhs.paddedSize();
	INS_VEC MID(mid);
	INS_VEC RHS(rhs);
	Unroll_Ternary<INS_VEC, OP>::apply_4(sz, pLhs, MID, RHS, pRes, oper);
	return result;
}


template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyTernaryOperation(typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& mid, const Vec<INS_VEC>& rhs)
{

	check_pair(rhs, mid);

	if (isScalar(mid))
	{
		return ApplyTernaryOperation<INS_VEC, OP>(lhs, mid.getScalarValue(), rhs);
	}
	if (isScalar(rhs))
	{
		return ApplyTernaryOperation<INS_VEC, OP>(lhs, mid, rhs.getScalarValue());
	}

	Vec<INS_VEC> result(rhs.size());
	auto pRes = result.start();
	auto pMid = mid.start();
	auto pRhs = rhs.start();
	OP oper;
	int sz = rhs.paddedSize();
	Unroll_Ternary<INS_VEC, OP>::apply_4(sz, lhs, pMid, pRhs, pRes, oper);
	return result;
}



