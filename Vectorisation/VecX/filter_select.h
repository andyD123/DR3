/****************************  filter_select.h   *******************************
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
#include "vec_double.h"
#include "instruction_traits.h"
#include "boolean_operations.h"
#include "accumulate_transform.h"
#include "binary_unitary_operations.h"
#include "unroll_operators.h"
#include "math_ops.h"
#include "filter_pipe_and_join.h" // remove
#include <tuple>


//helper overloads for functions for index into view or vector
template <typename INS_VEC>
inline unsigned int getIndex(const VecView<INS_VEC>& lhs, int i, int j)
{
	return static_cast<unsigned int>(lhs.idxStart()[i + j]);//index from view
}

template <typename INS_VEC>
inline unsigned int getIndex(const Vec<INS_VEC>& lhs, int i, int j)
{
	lhs.start();// to get rid of unused
	return static_cast<unsigned int>(i + j);
}


/*
common implementation for filtering a vec or a view using a vector of pre calulated boolean condition values
*/
template< template <class> typename VEC_TYPE, typename INS_VEC  >
VecView<INS_VEC>  ApplyFilterBImpl(const VecBool<INS_VEC>& condition, const VEC_TYPE<INS_VEC>& lhs)
{
	check_pair_different_type(condition,lhs);

	VecView<INS_VEC> vw(static_cast<size_t>(lhs.size()));
	auto pRes = vw.start();
	auto pIdx = vw.idxStart();
	auto pLhs = lhs.start();
	auto pCond = condition.start();
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 1 * width;
	int sz = static_cast<int>(lhs.size());
	int psn = 0;
	int WDTH = std::min(sz, width);
	for (int i = 0; i < sz; i += step)
	{
		for (int j = 0; (j < WDTH) && ((i + j) < sz); j++)
		{
			if (pCond[i + j])
			{
				pIdx[psn] = getIndex<INS_VEC>(lhs, i, j);
				pRes[psn] = pLhs[i + j];
				++psn;
			}
		}
	}
	vw.setSizeAndPad(psn);
	return vw;
};

/*
returns a view containing eleements from the vector which have corresponding true values in the VectorBool condition input
*/
template< typename INS_VEC>
VecView<INS_VEC>  ApplyFilterB(const VecBool<INS_VEC>& condition, const Vec<INS_VEC>& lhs)
{
	return ApplyFilterBImpl< Vec, INS_VEC>(condition, lhs);
};

/*
returns a view containing eleements from the vector which have corresponding true values in the VectorBool condition input
*/
template< typename INS_VEC>
VecView<INS_VEC>  ApplyFilterB(const VecBool<INS_VEC>& condition, const VecView<INS_VEC>& lhs)
{
	return ApplyFilterBImpl< VecView, INS_VEC>(condition, lhs);
}


/*
common implementation for filtering a vec or a view using a boolean lambda function
*/
template< template <class> typename VEC_TYPE, typename INS_VEC, typename OP  >
VecView<INS_VEC>  ApplyFilterImpl(OP& condition, const VEC_TYPE<INS_VEC>& lhs)
{
	check_vector(lhs);

	VecView<INS_VEC> vw(static_cast<size_t>(lhs.size())); 
	auto pRes = vw.start();
	auto pIdx = vw.idxStart();
	auto pLhs = lhs.start();

	int sz = lhs.size();// should be size() which goes to m_last  element
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 1 * width;
	INS_VEC LHS;

	int psn = 0;
	int WDTH = std::min(sz, width);
	for (int i = 0; i < sz; i += step)
	{
		LHS.load(pLhs + i);
		using boolVType = typename InstructionTraits<INS_VEC>::BoolType;
		boolVType		 COND;

		COND = condition(LHS) ;
		if (horizontal_or(COND))
		{
			for (int j = 0; j < WDTH && ((i + j) < sz); j++)
			{
				if (COND[j])
				{
					pIdx[psn] = getIndex<INS_VEC>(lhs, i, j);
					pRes[psn] = pLhs[i + j];
					++psn;
				}
			}
		}
	}

	vw.setSizeAndPad(psn);
	return vw;
};

/*
applys filter  boloean lammda  (condition) to a  values held in VectorView  if result is true values
are passed on to the result vvector view which is returned.
*/
template< typename INS_VEC, typename OP>
VecView<INS_VEC>  ApplyFilter(OP& condition, const VecView<INS_VEC>& lhs)
{
  return ApplyFilterImpl< VecView, INS_VEC, OP>(condition, lhs);

};

/*
* tests if any elements satisfy the condition before trying to do the conditional element copy
*/
template< typename INS_VEC, typename OP>
VecView<INS_VEC>  ApplyFilter(OP& condition, const Vec<INS_VEC>& lhs)
{
	return ApplyFilterImpl< Vec, INS_VEC, OP>(condition, lhs);	
};

/*
common implementation for filtering the first N elements a vec or a view that satisfy the boolean lambda function "condition"
*/
template< template <class> typename VEC_TYPE, typename INS_VEC, typename OP  >
VecView<INS_VEC>  CountedFilterImpl(OP& condition, const VEC_TYPE<INS_VEC>& lhs, int N)
{

	check_vector(lhs); //should just test that not scalar mode vecx

	VecView<INS_VEC> vw(static_cast<size_t>(lhs.size())); 
	auto pRes = vw.start();
	auto pIdx = vw.idxStart();
	auto pLhs = lhs.start();
	int sz = static_cast<int>(lhs.size());//reduction so we use size 
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 1 * width;
	INS_VEC LHS;
	int psn = 0;
	int WDTH = std::min(sz, width);
	for (int i = 0; i < sz; i += step)
	{
		LHS.load(pLhs + i);
		auto COND = (condition(LHS));
		if (horizontal_or(COND))
		{
			for (int j = 0; j < WDTH && ((i + j) < sz); j++)
			{
				if (COND[j])
				{
					pIdx[psn] = getIndex(lhs, i, j); 
					pRes[psn] = pLhs[i + j];
					++psn;
				}
			}
			// we have found first N of
			if (psn >= N)
			{
				vw.setSizeAndPad(N);
				return vw;
			}
		}
	}
	vw.setSizeAndPad(psn);
	return vw;
};

/*
returns a view of the first N elements of the vector that satisfy the condition 
*/
template< typename INS_VEC, typename OP>
VecView<INS_VEC>  ApplyCountedFilter(OP& condition, const Vec<INS_VEC>& lhs,int N)
{
	return CountedFilterImpl(condition, lhs, N);
};

/*
returns a view of the first N elements of the view that satisfy the condition
*/
template< typename INS_VEC, typename OP>
VecView<INS_VEC>  ApplyCountedFilter(OP& condition, const VecView<INS_VEC>& lhs, int N)
{
	return CountedFilterImpl(condition, lhs, N);
};

/*
common implementation of a function which splits values in a vec or view  into a tuple of views, the first view is for elements satisfying the condition,
the second for remaining elements
*/
template< template <class> typename VEC_TYPE  ,typename INS_VEC , typename OP  >
std::tuple<VecView<INS_VEC>, VecView<INS_VEC> >  BinaryFilterImpl(OP& condition, const  VEC_TYPE< INS_VEC>& lhs)  
{
	VecView<INS_VEC> vwTrue(static_cast<size_t>(lhs.size()));
	VecView<INS_VEC> vwFalse(static_cast<size_t>(lhs.size() ));

	auto pRes = vwTrue.start();
	unsigned int* pIdx = vwTrue.idxStart();
	auto pResFalse = vwFalse.start();
	auto pIdxFalse = vwFalse.idxStart();
	auto pLhs = lhs.start();
	const int width = InstructionTraits<INS_VEC>::width;
	int step = 1 * width;
	int sz = static_cast<int>(lhs.size());

	INS_VEC LHS;
	int psn = 0;
	int psnFalse = 0;
	int WDTH = std::min(sz, width);
	for (int i = 0; i < sz; i += step)
	{
		LHS.load(pLhs + i);
		auto COND(condition(LHS));

		for (int j = 0; j < WDTH && ((i + j) < sz); j++)
		{
			int posnOfIdx = getIndex(lhs, i, j);
			auto val = pLhs[i + j];

			if (COND[j])
			{
				pIdx[psn] = posnOfIdx;
				pRes[psn] = val;
				++psn;
			}
			else
			{
				pIdxFalse[psnFalse] = posnOfIdx;
				pResFalse[psnFalse] = val;
				++psnFalse;
			}
		}
	}

	vwTrue.setSizeAndPad(psn);
	vwFalse.setSizeAndPad(psnFalse);
	return std::make_tuple(std::move(vwTrue), std::move(vwFalse));
};

/*
returns a tuple of views on the input vector, the first view is for elements satisfying the condition, the second for remaining eleemnts
*/
template< typename INS_VEC, typename OP>
std::tuple<VecView<INS_VEC>, VecView<INS_VEC> >  ApplyBinaryFilter(OP& condition, const Vec<INS_VEC>& lhs)
{
	return  BinaryFilterImpl(condition, lhs);
};

/*
returns a tuple of views on the input view, the first view is for elements satisfying the condition, the second for remaining elements
*/
template< typename INS_VEC, typename OP>
std::tuple<VecView<INS_VEC>, VecView<INS_VEC> >  ApplyBinaryFilter(OP& condition, const VecView<INS_VEC>& lhs)
{
	return  BinaryFilterImpl(condition, lhs);
};




////////////////  applying operations to views ///////////////////////

/*
applies the OP to the view updating input in situ.
*/
template< typename INS_VEC, typename OP>
void ApplyUnitaryOperation(VecView<INS_VEC>& rhs, OP& oper)
{
	auto pRhs = rhs.start();
	int sz = rhs.fillSize();
	Unroll_Unitary<INS_VEC, OP>::apply_4(sz, pRhs, pRhs, oper);
}

/*
applies the OP and creates a new view.
*/
template< typename INS_VEC, typename OP>
VecView<INS_VEC> ApplyUnitaryOperation(const VecView<INS_VEC>& rhs, OP& oper)
{
	check_vector_for_filter(rhs);
	VecView<INS_VEC> result(rhs);
	auto pRes = result.start();
	auto pRhs = rhs.start();
	int sz = rhs.fillSize();
	Unroll_Unitary<INS_VEC, OP>::apply_4(sz, pRhs, pRes, oper);
	return result;
}



/*
applies the operation to the input view to generate a new view which it returns
*/

template< typename INS_VEC, typename OP>
VecView<INS_VEC> ApplyUnitaryOperationV(const  Vec<INS_VEC>& rhs, OP& oper)
{

	if (!rhs.isScalar())
	{
		check_vector_for_filter(rhs);
		VecView<INS_VEC> result(rhs);
		auto pRes = result.start();
		auto pRhs = rhs.start();
		int sz = rhs.paddedSize();
		Unroll_Unitary<INS_VEC, OP>::apply_4(sz, pRhs, pRes, oper);
		result.set_last(static_cast<int>(rhs.size()));
		return result;
	
	}
	else
	{
		auto val = rhs.getScalarValue();
		auto scalarRes = oper(INS_VEC(val))[0];
		VecView<INS_VEC> result;
		result =scalarRes;
		return result;
	}
	
}



/*
applies the binary OP to the two view inputs and returns a view object
*/
//optimistic assumes both lhs are same size 
template< typename INS_VEC, typename OP>
VecView<INS_VEC> ApplyBinaryOperation(OP& oper, const  VecView<INS_VEC>& lhs, const  VecView<INS_VEC>& rhs)
{
	check_view_pair(lhs, rhs);
	
	// to do checdks assert  or throw and overloads with scalar and Vec
	VecView<INS_VEC> result(rhs);// 
	auto pRes = result.start();
	auto pRhs = rhs.start();
	auto pLhs = lhs.start();
	int sz = rhs.fillSize();// should be fillsize,  ie padded last

	Unroll_Binary<INS_VEC, OP>::apply_4( sz,  pLhs,  pRhs, pRes, oper);
	//assert ( last i ssame for rhs and lhs ?
	result.set_last(static_cast<int>(rhs.last()) );
	return result;
}



///////////////////  combined operations and filters /////////////////////


// compounding  operations and filtering or writing back to results.
// combined  operation and filter to view 
// combined  operation on a view and write back to source

/*
applies the OP to the input vector in and also applies the (unitary boolean lambda) filter  to populate a
view with values satisfying the filter. The result vector and view are returned in a pair.
*/
template< typename INS_VEC, typename OP, typename FILT>
std::pair< Vec<INS_VEC>, VecView<INS_VEC> > ApplyOperationAndFilter(OP& oper, FILT& filt, const Vec<INS_VEC>& in)
{
	check_vector_for_filter(in);
	VecView<INS_VEC> vw(static_cast<size_t>(in.size()));
	auto pVw = vw.start();
	auto pIdx = vw.idxStart();
	Vec<INS_VEC> result(in.size());
	auto pRes = result.start();
	auto pRhs = in.start();
	//int sz = static_cast<int>(in.size());
	int sz = static_cast<int>(in.paddedSize());

	auto psn = Unroll_Unitary<INS_VEC, OP>::apply_4_filter(sz, pRhs, pRes, oper, pVw, filt, pIdx);
	vw.setSizeAndPad(psn);
	return std::make_pair(std::move(result), std::move(vw));
};


/*
applies the OP to the input vector in and also applies the (unitary boolean lambda) filter  to populates a pair of
view with values satisfying the filter conditions filt and filt2. The result vector and view are returned in a tuple.
*/
template< typename INS_VEC, typename OP, typename FILT1, typename FILT2>
std::tuple< Vec<INS_VEC>, VecView<INS_VEC>, VecView<INS_VEC> > ApplyOperationAndFilter(OP& oper, FILT1& filt1, FILT2& filt2, const Vec<INS_VEC>& in)
{
	check_vector_for_filter(in);
	VecView<INS_VEC> vw(static_cast<size_t>(in.size())); 
	auto pVw = vw.start();
	auto pIdx = vw.idxStart();

	VecView<INS_VEC> vw2(static_cast<size_t>(in.size()));
	auto pVw2 = vw2.start();
	auto pIdx2 = vw2.idxStart();


	Vec<INS_VEC> result(static_cast<int>(in.size()));
	auto pRes = result.start();
	auto pRhs = in.start();
	//int sz = static_cast<int>(in.size());
	int sz = static_cast<int>(in.paddedSize());

	int  psn1 = 0;
	int  psn2 = 0;
	Unroll_Unitary<INS_VEC, OP>::apply_4_filter(sz, pRhs, pRes, oper, pVw, filt1, pIdx, psn1, pVw2, filt2, pIdx2, psn2);
	vw.setSizeAndPad(psn1);
	vw2.setSizeAndPad(psn2);
	return std::make_tuple(std::move(result), std::move(vw), std::move(vw2));

};


/*
applies the OP to the view in and  scatter,  writes the results to the corresponding elements  of the result vector.
*/
template< typename INS_VEC, typename OP>
void ApplyUnitaryOperationWrite(OP& oper, const VecView<INS_VEC>& rhs, Vec<INS_VEC>& res)
{
	if (!rhs.isScalar() && !res.isScalar())
	{
		check_vector(rhs); //checks the view
		auto pRhs = rhs.start();
		int sz = rhs.last(); //because we are writing not necessarily full register
		auto pIdx = rhs.idxStart();
		auto pWrite = res.start();
		Unroll_Unitary<INS_VEC, OP>::apply_4_andWrite(sz, pRhs, pRhs, oper, pWrite, &pIdx[0]);
	}
	else
	{
		throw std::exception("ApplyUnitaryOperationWrite not suppported for scalars ");
	}

}

