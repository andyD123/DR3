/****************************  filter_pipe_and_join.h   *******************************
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
#include "filter_select.h"


/*
 Use "|"  for joining filters and ">" for joining operations
 use braces around sets of filters to control evaluation order

*/

namespace PIPE
{

template< typename INS_VEC>
VecView<INS_VEC> operator |(const Vec<INS_VEC>& rhs, const VecBool<INS_VEC>& condition)
{
	return ApplyFilter(condition, rhs);
}


template< typename INS_VEC, typename OP>
VecView<INS_VEC> operator |(const VecView<INS_VEC>& rhs, OP& condition)
{
	return ApplyFilter(condition, rhs);
}

template< typename INS_VEC, typename OP>
VecView<INS_VEC> operator |(const Vec<INS_VEC>& rhs, OP& condition)
{
	return ApplyFilter(condition, rhs);
}

template< typename INS_VEC, typename OP>
VecView<INS_VEC> operator |(Vec<INS_VEC>& rhs, OP& condition)
{
	return ApplyFilter(condition, rhs);
}



////////////////////
// vector checks are applied inside ApplyUnitaryOperation
template< typename INS_VEC, typename OP>
VecView<INS_VEC> operator > (  VecView<INS_VEC> rhs, OP& oper)
{
	ApplyUnitaryOperation( rhs, oper);
	return rhs;
}


template< typename INS_VEC, typename OP>
VecView<INS_VEC>& operator > ( OP& oper, VecView<INS_VEC>& rhs)
{
	ApplyUnitaryOperation(oper, rhs);
	return rhs;
}

template< typename INS_VEC, typename OP>
VecView<INS_VEC> operator > (OP& oper, const VecView<INS_VEC>& rhs)
{
	return ApplyUnitaryOperation(oper, rhs);
}


template< typename INS_VEC, typename OP>
VecView<INS_VEC>& operator > (Vec<INS_VEC>& rhs, OP& oper)
{
	ApplyUnitaryOperation(oper, rhs);
	return rhs;
}


template< typename INS_VEC, typename OP>
VecView<INS_VEC> operator > (const Vec<INS_VEC>& rhs, OP& oper)
{
	return ApplyUnitaryOperation(oper, rhs);
}


template< typename INS_VEC, typename OP>
VecView<INS_VEC> operator > ( OP& oper , const Vec<INS_VEC>& rhs )
{
	return ApplyUnitaryOperation(oper, rhs);
}


template< typename INS_VEC>
Vec<INS_VEC> operator |(const VecView<INS_VEC>& rhs, Vec<INS_VEC>& out)
{
	auto outRes(out);
	rhs.writeView(outRes);
	return outRes;
}


template< typename INS_VEC>
Vec<INS_VEC> operator |(VecView<INS_VEC>& rhs, const Vec<INS_VEC>& out)
{
	auto outRes(out);
	rhs.writeView(outRes);
	return outRes;
}


struct WriteOut
{};

template< typename INS_VEC>
void operator |(const VecView<INS_VEC>& rhs, WriteOut& out)
{
	//writes back to source to do
	rhs.writeView(out);
}

}// namespace PIPE

/*
 These expression templates are for use at register level combinations of operations
*/
namespace  JOIN
{


	template< typename LHS, typename RHS>
	struct CatOperation
	{
		CatOperation(const LHS& lhs, const RHS& rhs) :m_lhs(lhs), m_rhs(rhs) {}

		template <typename X>
		inline  auto operator()(const X& val) noexcept
		{
			return m_lhs(m_rhs(val));
		}

		//for use with accumulate
		template <typename X>
		inline  auto operator()(const X& lhs_arg, const X& rhs_arg) noexcept
		{
			return m_rhs(rhs_arg, m_lhs(lhs_arg));
		}
		LHS m_lhs;
		RHS m_rhs;
	};

	
	template<typename LHS, typename RHS>
	CatOperation< RHS, LHS> operator | (const LHS& lhs, const RHS& rhs)
	{
		return CatOperation<RHS, LHS>(rhs,lhs);
	}

	/*

		Boolean expression template conjuction for boolean lambdas

		auto isLessThanMinus10 = [](auto x) { return x < -10 };
		auto isGreaterThan10 = [](auto x) { return x > 10 };
		auto isLessThan20 = [](auto x) { return x < 20 };

		we can create simple logical conjunctions of boolean lambdas

		auto betweenTenAndTwenty = isGreaterThan10 && isLessThan20;

		auto isOutsideTenTwenty = !betweenTenAndTwenty;

		auto hasAbsGreaterThanTen = isLessThanMinus10 || isGreaterThan10;

	*/

	template<  typename RHS>
	struct NegateOperation
	{
		NegateOperation(const RHS& rhs) : m_rhs(rhs) {}

		template <typename INS_VEC>
		inline  auto operator()(const INS_VEC& val) noexcept
		{
			return !m_rhs(val);
		}
		RHS m_rhs;
	};


	template<  typename RHS>
	NegateOperation<RHS> operator ! (const RHS& rhs)
	{
		return NegateOperation< RHS>(rhs);
	}





	template< typename LHS, typename RHS>
	struct OROperation
	{
		OROperation(const LHS& lhs, const RHS& rhs) :m_lhs(lhs), m_rhs(rhs) {}

		template <typename X>
		inline  auto operator()(const X& val) noexcept
		{
			return m_lhs(val) || m_rhs(val);
		}
		LHS m_lhs;
		RHS m_rhs;
	};


	template< typename LHS, typename RHS>
	OROperation<  LHS, RHS> operator || (const LHS& lhs, const RHS& rhs)
	{
		return OROperation<LHS, RHS>(lhs, rhs);
	}



	template< typename LHS, typename RHS>
	struct AndOperation
	{
		AndOperation(const LHS& lhs, const RHS& rhs) :m_lhs(lhs), m_rhs(rhs) {}

		template <typename X>
		inline  auto operator()(const X& val) noexcept
		{
			return m_lhs(val) && m_rhs(val);
		}

		LHS m_lhs;
		RHS m_rhs;

	};


	template< typename LHS, typename RHS>
	AndOperation<  LHS, RHS> operator && (const LHS& lhs, const RHS& rhs)
	{
		return AndOperation<LHS, RHS>(lhs, rhs);
	}


}//namespace JOIN

