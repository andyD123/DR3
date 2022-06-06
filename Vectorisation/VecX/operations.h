/****************************  operations.h   *******************************
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
#include "apply_operation.h"
#include "math_ops.h"


//forward declarations

template<typename INS_VEC> 
class  Vec;

template<typename INS_VEC> 
class  VecBool;


template<typename INS_VEC>
class  VecD;


/******************************************  
                 operations 
*******************************************/

template<typename INS_VEC>
Vec<INS_VEC> operator -( const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC,  DR_CUBED::_unitaryMinus<INS_VEC> >(rhs);
}


template<typename INS_VEC>
Vec<INS_VEC> cdfnorm( const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_cdfnorm<INS_VEC> >(rhs);
}


template<typename INS_VEC>
Vec<INS_VEC> cdfnormD(const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_cdfnormD<INS_VEC> >(rhs);
}


template<typename INS_VEC>
Vec<INS_VEC>fabs(const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_fabs<INS_VEC> >(rhs);
}


template<typename INS_VEC>
Vec<INS_VEC>sqrt( const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_sqrt<INS_VEC> >(rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>exp( const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_exp<INS_VEC> >(rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>log( const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_log<INS_VEC> >(rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>pow(const Vec<INS_VEC>& rhs, typename InstructionTraits<INS_VEC>::FloatType rhs2)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_pow<INS_VEC> >(rhs,rhs2);
}

template<typename INS_VEC>
Vec<INS_VEC>pow(const Vec<INS_VEC>& rhs, const Vec<INS_VEC>& rhs2)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_pow<INS_VEC> >(rhs, rhs2);
}


template<typename INS_VEC>
Vec<INS_VEC>abs( const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_abs<INS_VEC> >(rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>floor( const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_floor<INS_VEC> >(rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>ceil( const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_ceil<INS_VEC> >(rhs);
}


template<typename INS_VEC>
Vec<INS_VEC>sin(const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_sin<INS_VEC> >(rhs);
}


template<typename INS_VEC>
Vec<INS_VEC>cos(const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_cos<INS_VEC> >(rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>tan(const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_tan<INS_VEC> >(rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>asin(const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_asin<INS_VEC> >(rhs);
}


template<typename INS_VEC>
Vec<INS_VEC>acos(const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_acos<INS_VEC> >(rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>atan(const Vec<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation<INS_VEC, DR_CUBED::_atan<INS_VEC> >(rhs);
}

//  Binary Operations ///

template<typename INS_VEC>
Vec<INS_VEC>max(const  Vec<INS_VEC>& lhs,const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_max<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>max(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_max<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>max(typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_max<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>min(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_min<INS_VEC> >(lhs,rhs);
}
template<typename INS_VEC>
Vec<INS_VEC>min(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_min<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>min(typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_min<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC> operator + (const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_plus<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>operator + (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_plus<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>operator + (typename InstructionTraits<INS_VEC>::FloatType lhs,const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_plus<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
const Vec<INS_VEC>& operator += (  Vec<INS_VEC>& lhs,const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperationM<INS_VEC, DR_CUBED::_plus<INS_VEC> >(lhs, rhs);
}

template<typename INS_VEC>
const Vec<INS_VEC>& operator += ( Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBinaryOperationM<INS_VEC, DR_CUBED::_plus<INS_VEC> >(lhs, rhs);
}

template<typename INS_VEC>
Vec<INS_VEC> operator - (const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_minus<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC> operator - (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_minus<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC> operator - (typename InstructionTraits<INS_VEC>::FloatType lhs,const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_minus<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
const Vec<INS_VEC> & operator -= (  Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperationM<INS_VEC, DR_CUBED::_minus<INS_VEC> >(lhs, rhs);
}

template<typename INS_VEC>
const Vec<INS_VEC>&  operator -= (  Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBinaryOperationM<INS_VEC, DR_CUBED::_minus<INS_VEC> >(lhs, rhs);
}


template<typename INS_VEC>
Vec<INS_VEC>operator * (const Vec<INS_VEC>& lhs,const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_times<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>operator * (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_times<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>operator * (typename InstructionTraits<INS_VEC>::FloatType lhs,const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_times<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
const Vec<INS_VEC> & operator *= ( Vec<INS_VEC>& lhs,const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperationM<INS_VEC, DR_CUBED::_times<INS_VEC> >(lhs, rhs);
}

template<typename INS_VEC>
const Vec<INS_VEC>&   operator *= ( Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBinaryOperationM<INS_VEC, DR_CUBED::_times<INS_VEC> >(lhs, rhs);
}


template<typename INS_VEC>
Vec<INS_VEC>operator / (const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_divide<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>operator / (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_divide<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC>operator / (typename InstructionTraits<INS_VEC>::FloatType lhs,const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_divide<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
const Vec<INS_VEC>& operator /= (Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBinaryOperationM<INS_VEC, DR_CUBED::_divide<INS_VEC> >(lhs, rhs);
}

template<typename INS_VEC>
const Vec<INS_VEC>& operator /= (Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	  return ApplyBinaryOperationM<INS_VEC, DR_CUBED::_divide<INS_VEC> >(lhs, rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator > (const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_greaterThan<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator > (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_greaterThan<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator > (typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_greaterThan<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator >= (const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_greaterThanEqual<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator >= (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_greaterThanEqual<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator >= (typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_greaterThanEqual<INS_VEC> >(lhs,rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator < (const Vec<INS_VEC>& lhs,const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_lessThan<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator < (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_lessThan<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator < (typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_lessThan<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator <= (const Vec<INS_VEC>& lhs,const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_lessThanEqual<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator <= (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_lessThanEqual<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator <= (typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_lessThanEqual<INS_VEC> >(lhs,rhs);

}


template<typename INS_VEC>
VecBool<INS_VEC> operator == (const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_Equal<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator == (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_Equal<INS_VEC> >(lhs,rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator == (typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_Equal<INS_VEC> >(lhs,rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator != (const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_NotEqual<INS_VEC> >(lhs,rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator & (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_BitwiseAnd<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator & (typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_BitwiseAnd<INS_VEC> >(lhs, rhs);
}
template<typename INS_VEC>
VecBool<INS_VEC> operator & (const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_BitwiseAnd<INS_VEC> >(lhs, rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator | (const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_BitwiseOr<INS_VEC> >(lhs, rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator | (typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_BitwiseOr<INS_VEC> >(lhs, rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator | (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_BitwiseOr<INS_VEC> >(lhs, rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator ^ (const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_BitwiseXOr<INS_VEC> >(lhs, rhs);
}





template<typename INS_VEC>
VecBool<INS_VEC> operator != (const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_NotEqual<INS_VEC> >(lhs,rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator != (typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_NotEqual<INS_VEC> >(lhs,rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator ! (const VecBool<INS_VEC>& rhs)
{
	return ApplyBooleanUnitaryOperation<INS_VEC, DR_CUBED::_NOT<INS_VEC> >(rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator || (const VecBool<INS_VEC> & lhs, const VecBool<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_OR<INS_VEC> >(lhs,rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator && (const VecBool<INS_VEC> & lhs, const VecBool<INS_VEC>& rhs)
{
	return ApplyBooleanBinaryOperation<INS_VEC, DR_CUBED::_AND<INS_VEC> >(lhs,rhs);
}


template<typename INS_VEC>
Vec<INS_VEC> iff( const VecBool<INS_VEC>& cond ,const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	return ApplySelectionOperation<INS_VEC>(cond,lhs,rhs);
}

template<typename INS_VEC>
Vec<INS_VEC> iff(const VecBool<INS_VEC>& cond, const Vec<INS_VEC>& data, typename InstructionTraits<INS_VEC>::FloatType lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplySelectionOperation<INS_VEC>(cond, data, lhs, rhs);
}


template<typename INS_VEC,typename COND >
Vec<INS_VEC> iff(COND& cond, typename InstructionTraits<INS_VEC>::FloatType lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplySelectionOperation<INS_VEC>(cond, lhs, rhs);
}


template<typename INS_VEC>
Vec<INS_VEC> FMA(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& mid,const Vec<INS_VEC>& rhs)
{
	return ApplyTernaryOperation<INS_VEC, DR_CUBED::_FMA<INS_VEC> >(lhs,mid, rhs);
}



template<typename INS_VEC>
Vec<INS_VEC> FMA(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& mid, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyTernaryOperation<INS_VEC, DR_CUBED::_FMA<INS_VEC> >(lhs, mid, rhs);
}

template<typename INS_VEC>
Vec<INS_VEC> FMA(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType mid, const Vec<INS_VEC>& rhs)
{
	return ApplyTernaryOperation<INS_VEC, DR_CUBED::_FMA<INS_VEC> >(lhs, mid, rhs);
}


template<typename INS_VEC>
Vec<INS_VEC> FMA(typename InstructionTraits<INS_VEC>::FloatType lhs, const Vec<INS_VEC>& mid, typename InstructionTraits<INS_VEC>::FloatType  rhs)
{
	return ApplyTernaryOperation<INS_VEC, DR_CUBED::_FMA<INS_VEC> >(lhs, mid, rhs);
}

template<typename INS_VEC>
Vec<INS_VEC> FMA(const Vec<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType mid, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return ApplyTernaryOperation<INS_VEC, DR_CUBED::_FMA<INS_VEC> >(lhs, mid, rhs);
}


/////////////////////////////////////////////////////////// AAD operations ///////////////////////////////////////////////
template<typename INS_VEC>
VecBool<INS_VEC> operator > (const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs.value() > rhs.value());
}

template<typename INS_VEC>
VecBool<INS_VEC> operator > (const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return (lhs.value() > rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator > (typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs > rhs.value());
}


template<typename INS_VEC>
VecBool<INS_VEC> operator >= (const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs.value() >= rhs.value());
}

template<typename INS_VEC>
VecBool<INS_VEC> operator >= (const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return (lhs.value() >= rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator >= (typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs >= rhs.value());
}


template<typename INS_VEC>
VecBool<INS_VEC> operator < (const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs.value() < rhs.value() );
}

template<typename INS_VEC>
VecBool<INS_VEC> operator < (const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return  (lhs.value() < rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator < (typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs < rhs.value() );
}


template<typename INS_VEC>
VecBool<INS_VEC> operator <= (const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs.value() <= rhs.value() );
}

template<typename INS_VEC>
VecBool<INS_VEC> operator <= (const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return (lhs.value() <= rhs);
}

template<typename INS_VEC>
VecBool<INS_VEC> operator <= (typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs <= rhs.value());
}


template<typename INS_VEC>
VecBool<INS_VEC> operator == (const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs.value() == rhs.value() );
}

template<typename INS_VEC>
VecBool<INS_VEC> operator == (const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return (lhs.value() == rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator == (typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs== rhs.value() );
}


template<typename INS_VEC>
VecBool<INS_VEC> operator != (const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs.value() != rhs.value() );
}

template<typename INS_VEC>
VecBool<INS_VEC> operator != (const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return (lhs.value() != rhs);
}


template<typename INS_VEC>
VecBool<INS_VEC> operator != (typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return (lhs != rhs.value());
}

template<typename INS_VEC>
VecD<INS_VEC> operator + (const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(lhs.value()+ rhs.value(),lhs.derivative()+rhs.derivative() );
}


template<typename INS_VEC>
VecD<INS_VEC>operator + (const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return VecD<INS_VEC>(lhs.value()+ rhs, lhs.derivative());
}

template<typename INS_VEC>
VecD<INS_VEC>operator + (typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return rhs + lhs;
}

template<typename INS_VEC>
const VecD<INS_VEC>& operator += ( VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return lhs = lhs + rhs;
	
}

template<typename INS_VEC>
const VecD<INS_VEC>& operator += (VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return lhs = lhs + rhs;
}

template<typename INS_VEC>
VecD<INS_VEC> operator - (const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(lhs.value() - rhs.value(),lhs.derivative()- rhs.derivative());

}

template<typename INS_VEC>
VecD<INS_VEC> operator - (const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return VecD<INS_VEC>(lhs.value()- rhs, lhs.derivative());
}

template<typename INS_VEC>
VecD<INS_VEC> operator - (typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(lhs -rhs.value() , rhs.derivative());
}

template<typename INS_VEC>
const VecD<INS_VEC> & operator -= (VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return lhs = lhs - rhs;
}

template<typename INS_VEC>
const VecD<INS_VEC>&  operator -= (VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return lhs = lhs - rhs;
}


template<typename INS_VEC>
VecD<INS_VEC>operator * (const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(lhs.value()* rhs.value(), 
		lhs.value()* rhs.derivative() + rhs.value()* lhs.derivative() 	);
}



template<typename INS_VEC>
VecD<INS_VEC>operator * (const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return VecD<INS_VEC>(lhs.value()* rhs, rhs* lhs.derivative());
}



template<typename INS_VEC>
VecD<INS_VEC>operator * (typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return rhs*lhs;
}



template<typename INS_VEC>
const VecD<INS_VEC> & operator *= (VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return lhs = lhs * rhs;
}

template<typename INS_VEC>
const VecD<INS_VEC>&   operator *= (VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return lhs = lhs * rhs;
}


template<typename INS_VEC>
VecD<INS_VEC>operator / (const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(lhs.value()/ rhs.value(),
		((lhs.derivative() - lhs.value()* rhs.derivative() / rhs.value()) / rhs.value() )  );
}

template<typename INS_VEC>
VecD<INS_VEC>operator / (const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return VecD<INS_VEC>(lhs.value()/ rhs,  lhs.derivative()/rhs);
}

template<typename INS_VEC>
VecD<INS_VEC>operator / (typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(lhs/rhs.value(), -lhs*rhs.derivative()/(rhs.value()*rhs.value()) );
}

template<typename INS_VEC>
const VecD<INS_VEC>& operator /= (VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return lhs = lhs / rhs;
}

template<typename INS_VEC>
const VecD<INS_VEC>& operator /= (VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return lhs = lhs / rhs;
}

//////////////////// aad //unitary ops

template<typename INS_VEC>
VecD<INS_VEC> operator -(const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(-rhs.value(), -rhs.derivative());
}

template<typename INS_VEC>
VecD<INS_VEC>sqrt(const VecD<INS_VEC>& rhs)
{
	auto SQ_val = sqrt(rhs.value());
	return VecD<INS_VEC>(SQ_val, 0.5*rhs.derivative()/ SQ_val);
}

template<typename INS_VEC>
VecD<INS_VEC>exp(const VecD<INS_VEC>& rhs)
{
	auto expVal = exp(rhs.value());
	return VecD<INS_VEC>(expVal, expVal* rhs.derivative() );
}

template<typename INS_VEC>
VecD<INS_VEC>log(const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(log(rhs.value()), rhs.derivative()/rhs.value() );
}

template<typename INS_VEC>
Vec<INS_VEC>abs(const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(abs(rhs.value()), abs(rhs.derivative() ) );
}

template<typename INS_VEC>
VecD<INS_VEC>floor(const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(floor(rhs.value()), floor(rhs.derivative()));
}

template<typename INS_VEC>
VecD<INS_VEC>ceil(const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(floor(rhs.value()), floor(rhs.derivative()));
}


template<typename INS_VEC>
VecD<INS_VEC>max(const  VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	const VecBool<INS_VEC> cond = (lhs.value() > rhs.value());
	const Vec<INS_VEC> lhsD = lhs.derivative();
	const Vec<INS_VEC> rhsD  = rhs.derivative();
	Vec<INS_VEC> deriv =  select(cond, lhsD,rhsD);
	Vec<INS_VEC> val = max(lhs.value(), rhs.value());
	VecD<INS_VEC> retVal(val, deriv);
	return retVal;
}

template<typename INS_VEC>
VecD<INS_VEC>max(const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return VecD<INS_VEC>(max(lhs.value(), rhs), select((lhs > rhs), (lhs.derivative(), 0.0)));
}

template<typename INS_VEC>
VecD<INS_VEC>max(typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(max(lhs, rhs.value), select( (lhs > rhs), 0.0,rhs.derivative()));
}


template<typename INS_VEC>
VecD<INS_VEC>min(const VecD<INS_VEC>& lhs, const VecD<INS_VEC>& rhs)
{
	return ApplyBinaryOperation<INS_VEC, DR_CUBED::_min<INS_VEC> >(lhs, rhs);

	const VecBool<INS_VEC> cond = (lhs.value() < rhs.value());
	const Vec<INS_VEC> lhsD = lhs.derivative();
	const Vec<INS_VEC> rhsD = rhs.derivative();
	Vec<INS_VEC> deriv = select(cond, lhsD, rhsD);
	Vec<INS_VEC> val = min(lhs.value(), rhs.value());
	VecD<INS_VEC> retVal(val, deriv);
	return retVal;
}

template<typename INS_VEC>
VecD<INS_VEC>min(const VecD<INS_VEC>& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
{
	return VecD<INS_VEC>(min(lhs.value(), rhs), select((lhs > rhs), (lhs.derivative(), 0.0)));
}

template<typename INS_VEC>
VecD<INS_VEC>min(typename InstructionTraits<INS_VEC>::FloatType lhs, const VecD<INS_VEC>& rhs)
{
	return VecD<INS_VEC>(min(lhs.value(), rhs), select((lhs < rhs), (lhs.derivative(), 0.0)));
}



