/****************************  math_ops.h   *******************************
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
#include "instruction_traits.h"
#include "apply_operation.h" // probably needs to be called extra operations

template<typename INS_VEC>
static  INS_VEC select(typename InstructionTraits<INS_VEC>::BoolType const& s, INS_VEC const& a, INS_VEC const& b)
{
	return select(s, a, b);
}

template<typename INS_VEC>
static  INS_VEC iff(typename InstructionTraits<INS_VEC>::BoolType const& s, INS_VEC const& a, INS_VEC const& b)
{
	return select(s, a, b);
}


template<typename INS_VEC, typename COND>
static  INS_VEC iff(COND& cond, const INS_VEC& condInput, INS_VEC const& a, INS_VEC const& b)
{
	return select(cond(condInput), a, b);
}


template<typename INS_VEC, typename COND>
static  INS_VEC iff(typename InstructionTraits<INS_VEC>::BoolType const& s, typename InstructionTraits<INS_VEC>::FloatType a, typename InstructionTraits<INS_VEC>::FloatType  b)
{
	return select(s, INS_VEC(a), INS_VEC(b));
}






namespace DR_CUBED
{

	template < typename INS_VEC>
	struct _FMA
	{
		inline INS_VEC operator()(const INS_VEC& lhs, const INS_VEC& mid, const INS_VEC& rhs)
		{
			return mul_add(lhs, mid, rhs);
		}

		inline INS_VEC operator()(const INS_VEC& lhs, const INS_VEC& mid, typename InstructionTraits<INS_VEC>::FloatType  val)
		{
			INS_VEC rhs(val);
			return mul_add(lhs, mid, rhs);
		}

		inline INS_VEC operator()(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType val, const INS_VEC& rhs)
		{
			INS_VEC mid(val);
			return mul_add(lhs, mid, rhs);
		}

		inline INS_VEC operator()(typename InstructionTraits<INS_VEC>::FloatType  lhsF, const INS_VEC& mid, const INS_VEC& rhs)
		{
			INS_VEC lhs(lhsF);
			return mul_add(lhs, mid, rhs);
		}

		inline INS_VEC operator()(typename InstructionTraits<INS_VEC>::FloatType lhsF, const INS_VEC& mid, typename InstructionTraits<INS_VEC>::FloatType  val)
		{
			INS_VEC lhs(lhsF);
			INS_VEC rhs(val);
			return mul_add(lhs, mid, rhs);
		}

		inline INS_VEC operator()(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType mid_val, typename InstructionTraits<INS_VEC>::FloatType  val)
		{
			INS_VEC rhs(val);
			INS_VEC mid(mid_val);
			return mul_add(lhs, mid, rhs);
		}

		inline INS_VEC operator()(typename InstructionTraits<INS_VEC>::FloatType lhsF, typename InstructionTraits<INS_VEC>::FloatType midF, typename InstructionTraits<INS_VEC>::FloatType  val)
		{
			INS_VEC lhs(lhsF);
			INS_VEC mid(midF);
			INS_VEC rhs(val);
			return mul_add(lhs, mid, rhs);
		}
	};

	template < typename INS_VEC>
	struct _plus
	{
		inline INS_VEC operator()(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs + rhs;
		}
	};

	template < typename INS_VEC>
	struct _minus
	{
		inline INS_VEC operator()(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs - rhs;
		}
	};

	template < typename INS_VEC>
	struct _times
	{
		inline INS_VEC operator()(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs * rhs;
		}
	};

	template < typename INS_VEC>
	struct _divide
	{
		inline INS_VEC operator()(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs / rhs;
		}
	};

	template < typename INS_VEC>
	struct _log
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return log(rhs);
		}
	};

	template < typename INS_VEC>
	struct _exp
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return exp(rhs);
		}
	};

	template < typename INS_VEC>
	struct _abs
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return abs(rhs);
		}
	};

	template < typename INS_VEC>
	struct _floor
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return floor(rhs);
		}
	};

	template < typename INS_VEC>
	struct _cdfnorminv
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return cdfnorminv(rhs);
		}
	};

	template < typename INS_VEC>
	struct _cdfnorm
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return cdfnorm(rhs);
		}
	};

	template < typename INS_VEC>
	struct  _cdfnormD
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return cdfnormD(rhs);
		}
	};

	template < typename INS_VEC>
	struct _unitaryMinus
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return (-rhs);
		}
	};

	template < typename INS_VEC>
	struct _ceil
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return ceil(rhs);
		}
	};

	template < typename INS_VEC>
	struct _fabs
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return fabs(rhs);
		}
	};

	template < typename INS_VEC>
	struct _sqrt
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return sqrt(rhs);
		}
	};

	template < typename INS_VEC>
	struct _sin
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return sin(rhs);
		}
	};

	template < typename INS_VEC>
	struct _cos
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return cos(rhs);
		}
	};

	template < typename INS_VEC>
	struct _tan
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return tan(rhs);
		}
	};

	template < typename INS_VEC>
	struct _asin
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return asin(rhs);
		}
	};

	template < typename INS_VEC>
	struct _acos
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return acos(rhs);
		}
	};

	template < typename INS_VEC>
	struct _atan
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return atan(rhs);
		}
	};

	template < typename INS_VEC>
	struct _sinh
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return sinh(rhs);
		}
	};

	template < typename INS_VEC>
	struct _cosh
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return cosh(rhs);
		}
	};

	template < typename INS_VEC>
	struct _tanh
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return tanh(rhs);
		}
	};

	template < typename INS_VEC>
	struct _asinh
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return asinh(rhs);
		}
	};

	template < typename INS_VEC>
	struct _acosh
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return acosh(rhs);
		}
	};

	template < typename INS_VEC>
	struct _atanh
	{
		inline INS_VEC operator()(const INS_VEC& rhs)
		{
			return atanh(rhs);
		}
	};

	template < typename INS_VEC>
	struct _max
	{
		inline INS_VEC operator()(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return max(lhs, rhs);
		}
	};

	template < typename INS_VEC>
	struct _min
	{
		inline INS_VEC operator()(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return min(lhs, rhs);
		}
	};

	template < typename INS_VEC>
	struct _pow
	{
		inline INS_VEC operator()(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return pow(lhs, rhs);
		}

		inline INS_VEC operator()(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			return pow(lhs, rhs);
		}

		inline INS_VEC operator()(typename InstructionTraits<INS_VEC>::FloatType lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			return pow(lhs, rhs);
		}
	};

	/////////////////////Bool /////////////

	template < typename INS_VEC>
	struct _lessThan
	{
		inline typename  InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs < rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			INS_VEC vRhs(rhs);
			return lhs < rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(typename InstructionTraits<INS_VEC>::FloatType lhs, const INS_VEC& rhs)
		{
			INS_VEC vLhs(lhs);
			return vLhs < rhs;
		}


		inline  bool apply(typename InstructionTraits<INS_VEC>::FloatType lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			return lhs < rhs;
		}

	};

	template < typename INS_VEC>
	struct _lessThanEqual
	{
		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs <= rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			INS_VEC vRhs(rhs);
			return lhs <= rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(typename InstructionTraits<INS_VEC>::FloatType lhs, const INS_VEC& rhs)
		{
			INS_VEC vLhs(lhs);
			return vLhs <= rhs;
		}

		inline  bool apply(typename InstructionTraits<INS_VEC>::FloatType lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			return lhs <= rhs;
		}

	};

	template < typename INS_VEC>
	struct _greaterThan
	{
		inline typename  InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs > rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			INS_VEC vRhs(rhs);
			return lhs > rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(typename InstructionTraits<INS_VEC>::FloatType lhs, const INS_VEC& rhs)
		{
			INS_VEC vLhs(lhs);
			return vLhs > rhs;
		}

		inline  bool apply(typename InstructionTraits<INS_VEC>::FloatType lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			return lhs > rhs;
		}

	};

	template < typename INS_VEC>
	struct _greaterThanEqual
	{
		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs >= rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			INS_VEC vRhs(rhs);
			return lhs >= rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(typename InstructionTraits<INS_VEC>::FloatType lhs, const INS_VEC& rhs)
		{
			INS_VEC vLhs(lhs);
			return vLhs >= rhs;
		}

		inline  bool apply(typename InstructionTraits<INS_VEC>::FloatType lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			return lhs >= rhs;
		}


	};

	template < typename INS_VEC>
	struct _Equal
	{
		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs == rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			INS_VEC vRhs(rhs);
			return lhs == rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(typename InstructionTraits<INS_VEC>::FloatType lhs, const INS_VEC& rhs)
		{
			INS_VEC vLhs(lhs);
			return vLhs == rhs;
		}

		inline  bool apply(typename InstructionTraits<INS_VEC>::FloatType lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			return lhs == rhs;
		}

	};

	template < typename INS_VEC>
	struct _NotEqual
	{
		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs != rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			INS_VEC vRhs(rhs);
			return lhs != rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(typename InstructionTraits<INS_VEC>::FloatType lhs, const INS_VEC& rhs)
		{
			INS_VEC vLhs(lhs);
			return vLhs != rhs;
		}

		inline  bool apply(typename InstructionTraits<INS_VEC>::FloatType lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			return lhs != rhs;
		}

	};

	template < typename INS_VEC>
	struct _AND
	{
		inline typename InstructionTraits<INS_VEC>::BoolType apply(const typename InstructionTraits<INS_VEC>::BoolType& lhs, const  typename InstructionTraits<INS_VEC>::BoolType& rhs)
		{
			return lhs && rhs;
		}

	};

	template < typename INS_VEC>
	struct _OR
	{
		inline typename InstructionTraits<INS_VEC>::BoolType apply(const typename InstructionTraits<INS_VEC>::BoolType& lhs, const typename InstructionTraits<INS_VEC>::BoolType& rhs)
		{
			return lhs || rhs;
		}

	};

	template < typename INS_VEC>
	struct _NOT
	{
		inline typename InstructionTraits<INS_VEC>::BoolType apply(const typename InstructionTraits<INS_VEC>::BoolType& rhs)
		{
			return !rhs;
		}

	};

	template < typename INS_VEC>
	struct _BitwiseOr
	{
		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs | rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			INS_VEC vRhs(rhs);
			return lhs | rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(typename InstructionTraits<INS_VEC>::FloatType lhs, const INS_VEC& rhs)
		{
			INS_VEC vLhs(lhs);
			return vLhs | rhs;
		}

	};

	template < typename INS_VEC>
	struct _BitwiseAnd
	{
		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs & rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			INS_VEC vRhs(rhs);
			return lhs & rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(typename InstructionTraits<INS_VEC>::FloatType lhs, const INS_VEC& rhs)
		{
			INS_VEC vLhs(lhs);
			return vLhs & rhs;
		}

	};

	template < typename INS_VEC>
	struct _BitwiseNOT
	{
		inline typename InstructionTraits<INS_VEC>::BoolType apply(const typename InstructionTraits<INS_VEC>::BoolType& rhs)
		{
			return ~rhs;
		}

	};



	template < typename INS_VEC>
	struct _BitwiseXOr
	{
		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, const INS_VEC& rhs)
		{
			return lhs^ rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(const INS_VEC& lhs, typename InstructionTraits<INS_VEC>::FloatType rhs)
		{
			INS_VEC vRhs(rhs);
			return lhs ^ rhs;
		}

		inline typename InstructionTraits<INS_VEC>::BoolType apply(typename InstructionTraits<INS_VEC>::FloatType lhs, const INS_VEC& rhs)
		{
			INS_VEC vLhs(lhs);
			return vLhs ^ rhs;
		}

	};



}