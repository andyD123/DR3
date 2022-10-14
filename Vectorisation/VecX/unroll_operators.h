/****************************  unroll_operators.h   *******************************
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
#include "boolean_operations.h"


template< typename INS_VEC, typename OP>
typename InstructionTraits<INS_VEC>::FloatType
ApplyBinaryOperationVec(const typename InstructionTraits<INS_VEC>::FloatType& rhs1, const typename InstructionTraits<INS_VEC>::FloatType& rhs2, OP& oper)
{
	INS_VEC lhs(rhs1);
	INS_VEC rhs(rhs2);
	auto res = oper(lhs, rhs)[0];
	return res;
}


// OpElement's several different use cases for loading and applying operations
// these are register that will be loaded, operated on and writen from
// they can be used in a single case or in a block to give aloop unroll effect
template<typename INS_VEC, int OFFSET>
struct RegisterElement 
{
	RegisterElement():value(InstructionTraits<INS_VEC>::nullValue)
	{
	}

	INS_VEC value;
	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	static const int width = InstructionTraits< INS_VEC>::width;
	static const int relativeOffset = InstructionTraits< INS_VEC>::width * OFFSET;

	 inline void load(Float* pData)
	{
		 value.load_a(pData + relativeOffset);
	}


	 inline void store(Float* pData)
	{
		 value.store_a(pData + relativeOffset);
	}

};


template<typename INS_VEC, typename OP, int OFFSET>
struct  UnitaryOpElement
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	RegisterElement<INS_VEC, OFFSET> elem_lhs;

	using  IdxType = typename InstructionTraits< INS_VEC>::IdxType;
	IdxType idx;

	uint32_t limit = 100000; //maxs size for scatter TO DO 


	UnitaryOpElement()
	{}


	inline void apply(Float* plhs, int i, Float* pOut, OP oper)
	{
		elem_lhs.load(plhs + i); 
		INS_VEC res = oper(elem_lhs.value);
		res.store_a(pOut + i + elem_lhs.relativeOffset);
	}

	//does not store back result but writes to source vector
	inline void calcAndWrite(Float* plhs, int i, Float* pOut, OP oper, Float* pRes, unsigned int* pIdx)
	{
		pOut = 0;
		idx.load(pIdx + i + elem_lhs.relativeOffset);
		elem_lhs.load(plhs + i);
		INS_VEC res = oper(elem_lhs.value);
		if constexpr (InstructionTraits< INS_VEC>::useScatter)
		{
			scatter(idx, limit, res, pRes);
		}
		else
		{
			for (int j =0; j < InstructionTraits< INS_VEC>::width; ++j)
			{
				pRes[j + idx[j]] = res[j];
			}
		}
	}




	template<typename FILT>
	inline void filter(Float* plhs, Float* pFilter, int i, FILT condition, int& psn, unsigned  int* pIdx)
	{
		// i is loop offset == number of elements in loop  * number of times round loop
		using boolVType = typename InstructionTraits<INS_VEC>::BoolType;
		const int width = InstructionTraits<INS_VEC>::width;
		boolVType COND;
		//assumed that element has been already loaded in apply 
		//so no offset 
		COND = condition(elem_lhs.value); //= from pLhs  load  pLhs + i;

		if (horizontal_or(COND))
		{
			for (int j = 0; j < width; j++) //across register
			{
				if (COND[j] != 0)
				{
					pIdx[psn] = (i + elem_lhs.relativeOffset + j); 
					pFilter[psn] = plhs[i + j + elem_lhs.relativeOffset];
					++psn;
				}
			}
		}
	}

};

template<typename INS_VEC, typename OP, int OFFSET>
struct  BinaryOpElement
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	RegisterElement<INS_VEC, OFFSET> elem_lhs;
	RegisterElement<INS_VEC, OFFSET> elem_rhs;


	//takes a unitary op for transform and a boolean unitary for condition 
	//plhs is input
	//prhs is existing results vector
	template<typename BOOL_OP> 
	inline void apply(Float* plhs, Float* prhs, int i,  OP oper, BOOL_OP condition)
	{
		using boolVType = typename InstructionTraits<INS_VEC>::BoolType;
		boolVType COND;

		elem_lhs.load(plhs + i);
		COND = condition(elem_lhs.value); //= from pLhs  load  pLhs + i;

		if (!horizontal_or(COND)) //if any of the test show we need a calculation
		{
			//common case do nothing
		}
		else
		{
			elem_rhs.load(prhs + i);
			INS_VEC res = oper(elem_lhs.value); //apply unitary to input value
			INS_VEC blended = select(COND, res, elem_rhs.value); //blend existing results 
			blended.store_a(prhs + i + elem_lhs.relativeOffset); //write back
		}

	}


	inline void apply(Float* plhs, Float* prhs, int i, Float* pOut, OP oper)
	{
		elem_lhs.load(plhs+i);
		elem_rhs.load(prhs+i);
		INS_VEC res = oper(elem_lhs.value, elem_rhs.value);
		res.store_a(pOut + i+ elem_lhs.relativeOffset);
	}



	 inline void apply(INS_VEC LHS, Float* prhs, int i, Float* pOut, OP oper)
	 {
		 elem_rhs.load(prhs+i);
		 INS_VEC res = oper(LHS, elem_rhs.value);
		 res.store_a(pOut + i+elem_rhs.relativeOffset);
	 }

	 inline void apply(Float* plhs, INS_VEC RHS, int i, Float* pOut, OP oper)
	 {
		 elem_lhs.load(plhs+i);
		 INS_VEC res = oper(elem_lhs.value, RHS);
		 res.store_a(pOut + i+ elem_lhs.relativeOffset);
	 }


};


template<typename INS_VEC, typename OP, int OFFSET>
struct  TernaryOpElement
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	RegisterElement<INS_VEC, OFFSET> elem_lhs;
	RegisterElement<INS_VEC, OFFSET> elem_mid;
	RegisterElement<INS_VEC, OFFSET> elem_rhs;



	inline void apply(Float* plhs, Float* pmid, Float* prhs, int i, Float* pOut, OP oper)
	{
		elem_lhs.load(plhs + i);
		elem_mid.load(pmid + i);
		elem_rhs.load(prhs + i);
		INS_VEC res = oper(elem_lhs.value, elem_mid.value, elem_rhs.value);
		res.store_a(pOut + i + elem_lhs.relativeOffset);
	}



	inline void apply(INS_VEC LHS, Float* pmid, Float* prhs, int i, Float* pOut, OP oper)
	{
		elem_rhs.load(prhs + i);
		elem_mid.load(pmid + i);
		INS_VEC res = oper(LHS, elem_mid.value, elem_rhs.value);
		res.store_a(pOut + i + elem_rhs.relativeOffset);
	}


	inline void apply(Float* plhs, INS_VEC MID, Float* prhs, int i, Float* pOut, OP oper)
	{
		elem_lhs.load(plhs + i);
		elem_rhs.load(prhs + i);
		INS_VEC res = oper(elem_lhs.value, MID, elem_rhs.value);
		res.store_a(pOut + i + elem_lhs.relativeOffset);
	}

	inline void apply(Float* plhs, Float* pmid, INS_VEC RHS, int i, Float* pOut, OP oper)
	{
		elem_lhs.load(plhs + i);
		elem_mid.load(pmid + i);
		INS_VEC res = oper(elem_lhs.value, elem_mid.value, RHS);
		res.store_a(pOut + i + elem_lhs.relativeOffset);
	}


	inline void apply(Float* plhs, INS_VEC MID, INS_VEC RHS, int i, Float* pOut, OP oper)
	{
		elem_lhs.load(plhs + i);
		INS_VEC res = oper(elem_lhs.value, MID, RHS);
		res.store_a(pOut + i + elem_lhs.relativeOffset);
	}


	inline void apply(INS_VEC LHS, Float* pmid, INS_VEC RHS, int i, Float* pOut, OP oper)
	{
		elem_mid.load(pmid + i);
		INS_VEC res = oper(LHS, elem_mid.value, RHS);
		res.store_a(pOut + i + elem_mid.relativeOffset);
	}

	inline void apply(INS_VEC LHS, INS_VEC MID, Float* prhs, int i, Float* pOut, OP oper)
	{
		elem_rhs.load(prhs + i);
		INS_VEC res = oper(LHS, MID, elem_rhs.value);
		res.store_a(pOut + i + elem_rhs.relativeOffset);
	}


};


struct UnUsed {};
template<typename INS_VEC, typename OP, int OFFSET>
struct  SelectOpElement
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	RegisterElement<INS_VEC, OFFSET> elem_test;
	RegisterElement<INS_VEC, OFFSET> elem_true;
	RegisterElement<INS_VEC, OFFSET> elem_false;
		
	inline void applyBool(Float* ptest, Float* plhs, Float* prhs, int i, Float* pOut)
	{
		elem_test.load(ptest + i);
		elem_true.load(plhs + i);
		elem_false.load(prhs + i);
	
		typename InstructionTraits<INS_VEC>::BoolType selectCond = boolConvert<INS_VEC>(elem_test.value);
		INS_VEC res = select(selectCond, elem_true.value, elem_false.value);
		res.store_a(pOut + i + elem_test.relativeOffset);
	}
	

	inline void applyBool(Float* ptest, INS_VEC LHS, Float* prhs, int i, Float* pOut)
	{
		elem_test.load(ptest + i);
		elem_false.load(prhs + i);
		typename InstructionTraits<INS_VEC>::BoolType selectCond = boolConvert<INS_VEC>(elem_test.value);
		INS_VEC res = select(selectCond, LHS, elem_false.value);
		res.store_a(pOut + i + elem_test.relativeOffset);
	}


	inline void applyBool(Float* ptest, Float* plhs, INS_VEC RHS, int i, Float* pOut)
	{
		elem_test.load(ptest + i);
		elem_true.load(plhs + i);
		typename InstructionTraits<INS_VEC>::BoolType selectCond = boolConvert<INS_VEC>(elem_test.value);
		INS_VEC res = select(selectCond, elem_true.value, RHS);
		res.store_a(pOut + i + elem_test.relativeOffset);
	}

	inline void applyBool(Float* ptest, INS_VEC LHS, INS_VEC RHS, int i, Float* pOut)
	{
		elem_test.load(ptest + i);
		typename InstructionTraits<INS_VEC>::BoolType selectCond = boolConvert<INS_VEC>(elem_test.value);
		INS_VEC res = select(selectCond, LHS, RHS);
		res.store_a(pOut + i + elem_test.relativeOffset);
	}



	inline void apply(Float* ptest, Float* plhs, Float* prhs, int i, Float* pOut, OP cond)
	{
		elem_test.load(ptest + i);
		elem_true.load(plhs + i);
		elem_false.load(prhs + i);
		INS_VEC res = select(cond(elem_test.value), elem_true.value, elem_false.value);
		res.store_a(pOut + i + elem_test.relativeOffset);
	}

	inline void apply(Float* ptest, INS_VEC LHS, INS_VEC RHS, int i, Float* pOut, OP cond)
	{
		elem_test.load(ptest + i);
		INS_VEC res = select(cond(elem_test.value), LHS, RHS);
		res.store_a(pOut + i + elem_test.relativeOffset);
	}

	inline void apply(Float* ptest, Float* plhs, INS_VEC RHS, int i, Float* pOut, OP cond)
	{

		elem_test.load(ptest + i);
		elem_true.load(plhs + i);
		INS_VEC res = select(cond(elem_test.value), elem_true.value, RHS);
		res.store_a(pOut + i + elem_test.relativeOffset);
	}

	inline void apply(Float* ptest, INS_VEC LHS, Float* prhs, int i, Float* pOut, OP cond)
	{
		elem_test.load(ptest + i);
		elem_false.load(prhs + i);
		INS_VEC res = select(cond(elem_test.value), LHS, elem_false.value);
		res.store_a(pOut + i + elem_test.relativeOffset);
	}

	template< typename TRUE_OP, typename FALSE_OP >
	inline void apply_f(Float* ptest, TRUE_OP trueOper, FALSE_OP falseOper, int i, Float* pOut, OP cond)
	{
		elem_test.load(ptest + i);
		auto tru = trueOper(elem_test.value);
		auto fls = falseOper(elem_test.value);
		auto res = select(cond(elem_test.value),tru , fls );
		res.store_a(pOut + i + elem_test.relativeOffset);
	}

};


template<typename INS_VEC, typename OP, int OFFSET>
struct  AccumulateOpElement
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	RegisterElement<INS_VEC, OFFSET> elem_rhs;
	RegisterElement<INS_VEC, OFFSET> elem_res;

	inline AccumulateOpElement(Float initVal = InstructionTraits<INS_VEC>::nullValue)
	{
		elem_res.value = initVal;
	}

	inline void setInitialAccumValue(Float initVal)
	{
		elem_res.value = initVal;
	}

	inline void apply(Float* prhs, int i, OP& oper)
	{
		elem_rhs.load(prhs + i);
		elem_res.value = oper(elem_res.value, elem_rhs.value);
	}

};

template<typename INS_VEC, typename OP, int OFFSET>
struct  AccumulateTransformOpElement
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	RegisterElement<INS_VEC, OFFSET> elem_rhs;
	RegisterElement<INS_VEC, OFFSET> elem_res;


	AccumulateTransformOpElement(Float initVal = InstructionTraits<INS_VEC>::nullValue)
	{
		elem_res.value = initVal;
	}

	void setInitialAccumValue(Float initVal)
	{
		elem_res.value = initVal;
	}

	template< typename OP_TRANS >
	inline void apply(Float* prhs, int i, OP& operAcc, OP_TRANS& operTransform)
	{
		elem_rhs.load(prhs + i);
		elem_res.value = operAcc(elem_res.value, operTransform(elem_rhs.value));
	}

};


template<typename INS_VEC, typename OP, int OFFSET> 
struct  AccumulateTransformBinaryOpElement
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;
	RegisterElement<INS_VEC, OFFSET> elem_lhs;
	RegisterElement<INS_VEC, OFFSET> elem_rhs;
	RegisterElement<INS_VEC, OFFSET> elem_res;

	AccumulateTransformBinaryOpElement(Float initVal = InstructionTraits<INS_VEC>::nullValue)
	{
		elem_res.value = initVal;
	}

	inline void setInitialAccumValue(Float initVal)
	{
		elem_res.value = initVal;
	}

	template< typename OP_TRANS>
	inline void apply(Float* plhs, Float* prhs, int i, OP& operAcc, OP_TRANS& operTransform)
	{
		elem_lhs.load(plhs + i);
		elem_rhs.load(prhs + i);
		elem_res.value = operAcc(elem_res.value, operTransform(elem_lhs.value, elem_rhs.value));
	}

};


////////////////////////////////////////////////////////////////////////////////

template< template<typename , typename, int > typename OP_ELEMENT, typename INS_VEC, typename OP>
struct UnrollBinaryBody_4
{
	const int width = InstructionTraits<INS_VEC>::width;
	const int UNROLL = 4;
	const int step = UNROLL * width;

	OP_ELEMENT<INS_VEC, OP, 0>  elem_0;
	OP_ELEMENT<INS_VEC, OP, 1>  elem_1;
	OP_ELEMENT<INS_VEC, OP, 2>  elem_2;
	OP_ELEMENT<INS_VEC, OP, 3>  elem_3;

	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	/////////////// ternary version//////////////////

	inline void apply(int sz, Float* pLhs, Float* pMid,Float* pRhs, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, pMid, pRhs, i, pOut, oper);
			elem_1.apply(pLhs, pMid, pRhs, i, pOut, oper);
			elem_2.apply(pLhs, pMid, pRhs, i, pOut, oper);
			elem_3.apply(pLhs, pMid, pRhs, i, pOut, oper);
		}


		for (; i <= sz - width; i += width)
		{
			elem_0.apply(pLhs, pMid, pRhs, i, pOut, oper);
		}
	}

	inline void apply(int sz, INS_VEC LHS, Float* pMid, Float* pRhs, Float* pOut, OP oper) //not called
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(LHS, pMid, pRhs, i, pOut, oper);
			elem_1.apply(LHS, pMid, pRhs, i, pOut, oper);
			elem_2.apply(LHS, pMid, pRhs, i, pOut, oper);
			elem_3.apply(LHS, pMid, pRhs, i, pOut, oper);

		}

		for (; i <= sz - width; i += width)
		{
			elem_0.apply(LHS, pMid, pRhs, i, pOut, oper);
		}

	}

	inline void apply_4_bool(int sz, Float* pTest, Float* pTrueVals, Float* pFalseVals, Float* pOut)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.applyBool(pTest, pTrueVals, pFalseVals, i, pOut);
			elem_1.applyBool(pTest, pTrueVals, pFalseVals, i, pOut);
			elem_2.applyBool(pTest, pTrueVals, pFalseVals, i, pOut);
			elem_3.applyBool(pTest, pTrueVals, pFalseVals, i, pOut);
		}


		for (; i <= sz - width; i += width)
		{
			elem_0.applyBool(pTest, pTrueVals, pFalseVals, i, pOut);
		}

	}



	inline void apply_4_bool(int sz, Float* pTest, Float* pTrueVals, INS_VEC RHS, Float* pOut)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.applyBool(pTest, pTrueVals, RHS, i, pOut);
			elem_1.applyBool(pTest, pTrueVals, RHS, i, pOut);
			elem_2.applyBool(pTest, pTrueVals, RHS, i, pOut);
			elem_3.applyBool(pTest, pTrueVals, RHS, i, pOut);
		}


		for (; i <= sz - width; i += width)
		{
			elem_0.applyBool(pTest, pTrueVals, RHS, i, pOut);
		}

	}



	inline void apply_4_bool(int sz, Float* pTest, INS_VEC LHS , Float* pFalseVals, Float* pOut)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.applyBool(pTest, LHS, pFalseVals, i, pOut);
			elem_1.applyBool(pTest, LHS, pFalseVals, i, pOut);
			elem_2.applyBool(pTest, LHS, pFalseVals, i, pOut);
			elem_3.applyBool(pTest, LHS, pFalseVals, i, pOut);
		}


		for (; i <= sz - width; i += width)
		{
			elem_0.applyBool(pTest, LHS, pFalseVals, i, pOut);
		}
	}


	inline void apply_4_bool(int sz, Float* pTest, INS_VEC LHS, INS_VEC RHS, Float* pOut)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.applyBool(pTest, LHS, RHS, i, pOut);
			elem_1.applyBool(pTest, LHS, RHS, i, pOut);
			elem_2.applyBool(pTest, LHS, RHS, i, pOut);
			elem_3.applyBool(pTest, LHS, RHS, i, pOut);
		}


		for (; i <= sz - width; i += width)
		{
			elem_0.applyBool(pTest, LHS, RHS, i, pOut);
		}
	}




	inline void apply(int sz, Float* pLhs, INS_VEC MID, Float* pRhs, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, MID, pRhs, i, pOut, oper);
			elem_1.apply(pLhs, MID, pRhs, i, pOut, oper);
			elem_2.apply(pLhs, MID, pRhs, i, pOut, oper);
			elem_3.apply(pLhs, MID, pRhs, i, pOut, oper);
		}


		for (; i <= sz - width; i += width)
		{
			elem_0.apply(pLhs, MID, pRhs, i, pOut, oper);
		}
	}

	inline void apply(int sz, Float* pLhs, Float* pMid, INS_VEC RHS, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, pMid, RHS, i, pOut, oper);
			elem_1.apply(pLhs, pMid, RHS, i, pOut, oper);
			elem_2.apply(pLhs, pMid, RHS, i, pOut, oper);
			elem_3.apply(pLhs, pMid, RHS, i, pOut, oper);
		}

		for (; i <= sz - width; i += width)
		{
			elem_0.apply(pLhs, pMid, RHS, i, pOut, oper);
		}
	}

	template<typename TRUE_OP, typename FALSE_OP>
	inline void apply_f(int sz, Float* pLhs, TRUE_OP& trueOp, FALSE_OP& falseOp, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply_f(pLhs, trueOp, falseOp, i, pOut, oper);
			elem_1.apply_f(pLhs, trueOp, falseOp, i, pOut, oper);
			elem_2.apply_f(pLhs, trueOp, falseOp, i, pOut, oper);
			elem_3.apply_f(pLhs, trueOp, falseOp, i, pOut, oper);
		}

		for (; i <= sz - width; i += width)
		{
			elem_0.apply_f(pLhs, trueOp, falseOp, i, pOut, oper);
		}
	}

	inline void apply(int sz, INS_VEC LHS, INS_VEC MID, Float* pRhs, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(LHS, MID, pRhs, i, pOut, oper);
			elem_1.apply(LHS, MID, pRhs, i, pOut, oper);
			elem_2.apply(LHS, MID, pRhs, i, pOut, oper);
			elem_3.apply(LHS, MID, pRhs, i, pOut, oper);
		}


		for (; i <= sz - width; i += width)
		{
			elem_0.apply(LHS, MID, pRhs, i, pOut, oper);
		}
	}

	inline void apply(int sz, INS_VEC LHS,  Float* pMid, INS_VEC RHS, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(LHS, pMid, RHS, i, pOut, oper);
			elem_1.apply(LHS, pMid, RHS, i, pOut, oper);
			elem_2.apply(LHS, pMid, RHS, i, pOut, oper);
			elem_3.apply(LHS, pMid, RHS, i, pOut, oper);
		}


		for (; i <= sz - width; i += width)
		{
			elem_0.apply(LHS, pMid, RHS, i, pOut, oper);
		}
	}

	inline void apply(int sz, Float* pLhs, INS_VEC MID, INS_VEC RHS, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, MID, RHS, i, pOut, oper);
			elem_1.apply(pLhs, MID, RHS, i, pOut, oper);
			elem_2.apply(pLhs, MID, RHS, i, pOut, oper);
			elem_3.apply(pLhs, MID, RHS, i, pOut, oper);
		}


		for (; i <= sz - width; i += width)
		{
			elem_0.apply(pLhs, MID, RHS, i, pOut, oper);
		}
	}

	/////////////// binary versions //////////////////

	inline void apply(int sz, Float* pLhs, Float* pRhs, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, pRhs, i, pOut, oper);
			elem_1.apply(pLhs, pRhs, i, pOut, oper);
			elem_2.apply(pLhs, pRhs, i, pOut, oper);
			elem_3.apply(pLhs, pRhs, i, pOut, oper);
		}


		for (; i <= sz - width; i += width)
		{
			elem_0.apply(pLhs, pRhs, i, pOut, oper);
		}
	}

	inline void apply(int sz, INS_VEC LHS, Float* pRhs, Float* pOut, OP oper)
	{

		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(LHS, pRhs, i, pOut, oper);
			elem_1.apply(LHS, pRhs, i, pOut, oper);
			elem_2.apply(LHS, pRhs, i, pOut, oper);
			elem_3.apply(LHS, pRhs, i, pOut, oper);

		}

		for (; i <= sz - width; i += width)
		{
			elem_0.apply(LHS, pRhs, i, pOut, oper);
		}

	}

	inline void apply(int sz, Float* pLhs, INS_VEC RHS, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, RHS, i, pOut, oper);
			elem_1.apply(pLhs, RHS, i, pOut, oper);
			elem_2.apply(pLhs, RHS, i, pOut, oper);
			elem_3.apply(pLhs, RHS, i, pOut, oper);

		}

		for (; i <= sz - width; i += width)
		{
			elem_0.apply(pLhs, RHS, i, pOut, oper);

		}
	}

	UnrollBinaryBody_4() {}

	UnrollBinaryBody_4(bool singularInit, Float initVal) 
		: elem_0(singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue),
		elem_1(singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue),
		elem_2(singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue),
		elem_3(singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue)
	{};



	template<typename  OPER_TRANSFORM>
	inline void apply(int sz, Float* pLhs, Float* pRhs, OP oper,OPER_TRANSFORM operTransform,int& i)
	{
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, pRhs, i, oper, operTransform);
			elem_1.apply(pLhs, pRhs, i, oper, operTransform);
			elem_2.apply(pLhs, pRhs, i, oper, operTransform);
			elem_3.apply(pLhs, pRhs, i, oper, operTransform);
		}

		if (i > 0)
		{
			// consolidate across registers
			elem_0.elem_res.value = oper(elem_0.elem_res.value, elem_1.elem_res.value);
			elem_2.elem_res.value = oper(elem_2.elem_res.value, elem_3.elem_res.value);
			elem_0.elem_res.value = oper(elem_0.elem_res.value, elem_2.elem_res.value);
		}
	}

};


template< template<typename, typename, int > typename OP_ELEMENT, typename INS_VEC, typename OP >
struct UnrollBinaryBody_1
{
	const int width = InstructionTraits<INS_VEC>::width;
	const int UNROLL = 1;
	const int step = UNROLL * width;

	OP_ELEMENT<INS_VEC, OP, 0>  elem_0;

	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	/////////////////////ternary ////////////////////

	inline void apply(int sz, Float* pLhs, Float* pMid, Float* pRhs, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, pMid, pRhs, i, pOut, oper);
		}
	}

	inline void apply(int sz, INS_VEC LHS, Float* pMid, Float* pRhs, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(LHS, pMid, pRhs, i, pOut, oper);
		}
	}

	inline void apply(int sz, Float* pLhs, INS_VEC MID, Float* pRhs, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, MID, pRhs, i, pOut, oper);
		}
	}

	inline void apply(int sz, Float* pLhs, Float* pMid, INS_VEC RHS, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, pMid, RHS, i, pOut, oper);
		}
	}



	///////////////////sparse update //////////////////////////

	template <typename BOOL_OP>
	inline void apply(int sz, Float* pLhs, Float* pRhs, /*Float* pOut,*/ OP oper, BOOL_OP selectOp)
	{
		int i = 0;
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, pRhs, i,/* pOut,*/ oper, selectOp);
		}
	}
	///////////////////////////

	inline void apply(int sz, Float* pLhs, Float* pRhs, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, pRhs, i, pOut, oper);
		}
	}

	inline void apply(int sz, INS_VEC LHS, Float* pRhs, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(LHS, pRhs, i, pOut, oper);
		}
	}


	inline void apply(int sz, Float* pLhs, INS_VEC RHS, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, RHS, i, pOut, oper);
		}

	}


};


template< template<typename, typename, int > typename OP_ELEMENT, typename INS_VEC, typename OP >
struct UnrollUnitaryBody_1
{
	const int width = InstructionTraits<INS_VEC>::width;
	const int UNROLL = 1;
	const int step = UNROLL * width;

	OP_ELEMENT<INS_VEC, OP, 0>  elem_0;

	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	inline void apply(int sz,  Float* pRhs, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply( pRhs, i, pOut, oper);
		}
	}

	inline void apply(int sz,  INS_VEC RHS, Float* pOut, OP oper)
	{
		int i = 0;
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply( RHS, i, pOut, oper);
		}
	}


	UnrollUnitaryBody_1() {}

	UnrollUnitaryBody_1(bool singularInit, Float initVal): elem_0(singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue) {};

	inline void apply(int sz, Float* pRhs, OP oper, int& i)
	{
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pRhs, i, oper);
		}
	}

	template<typename OPER_TRANSFORM>
	inline void apply(int sz, Float* pRhs, OP oper, OPER_TRANSFORM operTransform, int& i)
	{
		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pRhs, i, oper, operTransform);
		}
	}


};


template< template<typename, typename, int > typename OP_ELEMENT, typename INS_VEC, typename OP >
struct UnrollUnitaryBody_4
{
	const int width = InstructionTraits<INS_VEC>::width;
	const int UNROLL = 4;
	const int step = UNROLL * width;

	OP_ELEMENT<INS_VEC, OP, 0>  elem_0;
	OP_ELEMENT<INS_VEC, OP, 1>  elem_1;
	OP_ELEMENT<INS_VEC, OP, 2>  elem_2;
	OP_ELEMENT<INS_VEC, OP, 3>  elem_3;

	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	inline void apply(int sz, INS_VEC LHS, Float* pOut, OP oper)
	{

		int i = 0;

		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(LHS, i, pOut, oper);
			elem_1.apply(LHS, i, pOut, oper);
			elem_2.apply(LHS, i, pOut, oper);
			elem_3.apply(LHS, i, pOut, oper);
		}

		for (; i <= sz - width; i += width)
		{
			elem_0.apply(LHS, i, pOut, oper);
		}
	}

	inline void apply(int sz, Float* pLhs, Float* pOut, OP oper)
	{

		int i = 0;

		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, i, pOut, oper);
			elem_1.apply(pLhs, i, pOut, oper);
			elem_2.apply(pLhs, i, pOut, oper);
			elem_3.apply(pLhs, i, pOut, oper);
		}

		for (; i <= sz - width; i += width)
		{
			elem_0.apply(pLhs, i, pOut, oper);
		}
	}


	UnrollUnitaryBody_4() {}

	UnrollUnitaryBody_4(bool singularInit, Float initVal) : elem_0(singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue), 
		elem_1(singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue),
		elem_2(singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue),
		elem_3(singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue)
	{};

	inline void apply(int sz, Float* pRhs, OP oper, int& i)
	{
		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pRhs, i, oper);
			elem_1.apply(pRhs, i, oper);
			elem_2.apply(pRhs, i, oper);
			elem_3.apply(pRhs, i, oper);
		}

		// consolidate across registers
		elem_0.elem_res.value = oper(elem_0.elem_res.value, elem_1.elem_res.value);
		elem_2.elem_res.value = oper(elem_2.elem_res.value, elem_3.elem_res.value);
		elem_0.elem_res.value = oper(elem_0.elem_res.value, elem_2.elem_res.value);
	}

	template<typename OPER_TRANSFORM>
	inline void apply(int sz, Float* pRhs, OP oper, OPER_TRANSFORM operTransform, int& i)
	{

		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pRhs, i, oper, operTransform);
			elem_1.apply(pRhs, i, oper, operTransform);
			elem_2.apply(pRhs, i, oper, operTransform);
			elem_3.apply(pRhs, i, oper, operTransform);
		}

		if (i > 0)
		{
			// consolidate across registers
			elem_0.elem_res.value = oper(elem_0.elem_res.value, elem_1.elem_res.value);
			elem_2.elem_res.value = oper(elem_2.elem_res.value, elem_3.elem_res.value);
			elem_0.elem_res.value = oper(elem_0.elem_res.value, elem_2.elem_res.value);
		}
	}



	//filter to view after transform 
	template<typename FILT>
	inline int apply_filter(int sz, Float* pLhs, Float* pOut, OP oper, Float* pFilt, FILT filt, unsigned int* pIdx)
	{
		int psn = 0;
		int i = 0;

		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, i, pOut, oper);
			elem_0.filter(pLhs, pFilt, i, filt, psn, pIdx);
			elem_1.apply(pLhs, i, pOut, oper);
			elem_1.filter(pLhs, pFilt, i, filt, psn, pIdx);
			elem_2.apply(pLhs, i, pOut, oper);
			elem_2.filter(pLhs, pFilt, i, filt, psn, pIdx);
			elem_3.apply(pLhs, i, pOut, oper);
			elem_3.filter(pLhs, pFilt, i, filt, psn, pIdx);

		}

		for (; i <= sz - width; i += width)
		{
			elem_0.apply(pLhs, i, pOut, oper);
			elem_0.filter(pLhs, pFilt, i, filt, psn, pIdx);
		}

		return psn;
	}


	template<typename FILT1, typename FILT2>
	inline void apply_4_filter(int sz, Float* pLhs, Float* pOut, OP oper,
		Float* pVw1, FILT1 filt1, unsigned int* pIdx1, int& psn1,
		Float* pVw2, FILT2 filt2, unsigned int* pIdx2, int& psn2)
	{
		psn1 = 0;
		psn2 = 0;
		int i = 0;

		for (; i < sz - (UNROLL * width); i += step)
		{
			elem_0.apply(pLhs, i, pOut, oper);
			elem_0.filter(pLhs, pVw1, i, filt1, psn1, pIdx1);
			elem_0.filter(pLhs, pVw2, i, filt2, psn2, pIdx2);

			elem_1.apply(pLhs, i, pOut, oper);
			elem_1.filter(pLhs, pVw1, i, filt1, psn1, pIdx1);
			elem_1.filter(pLhs, pVw2, i, filt2, psn2, pIdx2);

			elem_2.apply(pLhs, i, pOut, oper);
			elem_2.filter(pLhs, pVw1, i, filt1, psn1, pIdx1);
			elem_2.filter(pLhs, pVw2, i, filt2, psn2, pIdx2);

			elem_3.apply(pLhs, i, pOut, oper);
			elem_3.filter(pLhs, pVw1, i, filt1, psn1, pIdx1);
			elem_3.filter(pLhs, pVw2, i, filt2, psn2, pIdx2);

		}

		for (; i <= sz - width; i += width)
		{
			elem_0.apply(pLhs, i, pOut, oper);
			elem_0.filter(pLhs, pVw1, i, filt1, psn1, pIdx1);
			elem_0.filter(pLhs, pVw2, i, filt2, psn2, pIdx2);
		}
	}


	inline void apply_andWrite(int sz, Float* pLhs, Float* pOut, OP oper, Float* pWrite, unsigned int* pIdx)
	{

		int i = 0;

		for (; i <= sz - (UNROLL * width); i += step)
		{
			elem_0.calcAndWrite(pLhs, i, pOut, oper, pWrite, pIdx);
			elem_1.calcAndWrite(pLhs, i, pOut, oper, pWrite, pIdx);
			elem_2.calcAndWrite(pLhs, i, pOut, oper, pWrite, pIdx);
			elem_3.calcAndWrite(pLhs, i, pOut, oper, pWrite, pIdx);

		}

		for (; i <= sz - width; i += width)
		{
			elem_0.calcAndWrite(pLhs, i, pOut, oper, pWrite, pIdx);
		}

		//last element calc all in register and then write out
		if (i < sz)
		{
			elem_0.apply(pLhs, i, pOut, oper);

			for (; i < sz; i++)
			{
				pWrite[pIdx[i]] = pOut[i];
			}
		}
	}

};


/*
sz is the size of the padded data 
last elenent is register width away from the padded end

sizes arguments should be mod register width if all elements are to get calc results
*/
template<typename INS_VEC, typename OP>
struct Unroll_Binary
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	static inline void apply_4(int sz, Float* pLhs, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< BinaryOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, pLhs, pRhs, pOut, oper);
	}

	static inline void apply_4(int sz, INS_VEC LHS, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< BinaryOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, LHS, pRhs, pOut, oper);
	}

	static inline void apply_4(int sz, Float* pLhs, INS_VEC RHS, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< BinaryOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, pLhs, RHS, pOut, oper);
	}

	template <typename BOOL_OP>
	static inline void apply_1(int sz, Float* pLhs, Float* pRhs, /*Float* pOut,*/ OP oper, BOOL_OP selectionOp)
	{
		UnrollBinaryBody_1< BinaryOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, pLhs, pRhs, /*pOut,*/ oper, selectionOp);
	}

	static inline void apply_1(int sz, Float* pLhs, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_1< BinaryOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, pLhs, pRhs, pOut, oper);		
	}

	static inline void apply_1(int sz, INS_VEC LHS, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_1< BinaryOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, LHS, pRhs, pOut, oper);
	}

	static inline void apply_1(int sz, Float* pLhs, INS_VEC RHS, Float* pOut, OP oper)
	{
		UnrollBinaryBody_1< BinaryOpElement, INS_VEC, OP>  binaryUnrolled;
		binaryUnrolled.apply(sz, pLhs, RHS, pOut, oper);
	}

};



//  apply_4_filter 
//  transforms then filters
//  we can load  data into register perform a calculation
//  follow this by writing out results and also performing a filter 
//  to create a view
//  used to calculate common case  and split out rare , saving a traversal
template<typename INS_VEC, typename OP>
struct Unroll_Unitary
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	template<typename FILT1, typename FILT2>
	static void apply_4_filter(int sz, Float* pLhs, Float* pOut, OP oper,
					Float*  pVw1, FILT1 filt1, unsigned int*  pIdx1,int& psn1,
					Float*  pVw2, FILT2 filt2, unsigned int* pIdx2, int& psn2 )
	{
		UnrollUnitaryBody_4< UnitaryOpElement, INS_VEC, OP > unitaryUnrolled;
		unitaryUnrolled.apply_4_filter(sz, pLhs, pOut, oper,
										pVw1, filt1, pIdx1, psn1,
										pVw2, filt2, pIdx2, psn2);

	}

		
	template<typename FILT>
	static inline int apply_4_filter(int sz, Float* pLhs,  Float* pOut, OP oper,  Float* pFilt, FILT filt,unsigned int* pIdx)
	{
		UnrollUnitaryBody_4< UnitaryOpElement, INS_VEC, OP > unitaryUnrolled;
		return unitaryUnrolled.apply_filter(sz, pLhs, pOut, oper, pFilt, filt, pIdx);
	}


	static inline void apply_1(int sz, Float* pLhs, Float* pOut, OP oper)
	{
		UnrollUnitaryBody_1< UnitaryOpElement, INS_VEC, OP > unitaryUnrolled;
		unitaryUnrolled.apply( sz, pLhs,  pOut, oper);
	}

	// sz  is expected to be padded length
	// to end of padded array.
	static inline void apply_4(int sz, Float* pLhs, Float* pOut, OP oper)
	{
		UnrollUnitaryBody_4< UnitaryOpElement, INS_VEC, OP > unitaryUnrolled;
		unitaryUnrolled.apply( sz, pLhs,  pOut, oper);
	}

	// sz  is expected to be padded length
	// to end of padded array.
	static inline void apply_4_andWrite(int sz, Float* pLhs, Float* pOut, OP oper, Float* pWrite,  unsigned int* pIdx)
	{
		// NEED tests
		UnrollUnitaryBody_4< UnitaryOpElement, INS_VEC, OP > unitaryUnrolled;
		unitaryUnrolled.apply_andWrite( sz,  pLhs,  pOut, oper, pWrite, pIdx);
	}

};



/*
sz is the size of the padded data
last element is register width away from the padded end

sizes arguments should be mod register width if all elements are to get calc results
*/
template<typename INS_VEC, typename OP>
struct Unroll_Ternary
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	static inline void apply_4(int sz, Float* pLhs, Float* pMid, Float* pRhs, Float* pOut, OP oper)
	{	
		UnrollBinaryBody_4< TernaryOpElement, INS_VEC, OP>  ternaryUnrolled;
		ternaryUnrolled.apply(sz, pLhs, pMid, pRhs, pOut, oper);
	}


	static inline void apply_4(int sz, INS_VEC LHS,  Float* pMid, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< TernaryOpElement, INS_VEC, OP>  ternaryUnrolled;
		ternaryUnrolled.apply(sz, LHS, pMid, pRhs, pOut, oper);
	}


	static inline void apply_4(int sz, Float* pLhs, INS_VEC MID, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< TernaryOpElement, INS_VEC, OP>  ternaryUnrolled;
		ternaryUnrolled.apply(sz, pLhs, MID, pRhs, pOut, oper);
	}


	static inline void apply_4(int sz, Float* pLhs,  Float* pMid,  INS_VEC RHS, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< TernaryOpElement, INS_VEC, OP>  ternaryUnrolled;
		ternaryUnrolled.apply(sz, pLhs, pMid, RHS, pOut, oper);
	}

	static inline void apply_4(int sz, Float* pLhs, INS_VEC MID, INS_VEC RHS, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< TernaryOpElement, INS_VEC, OP>  ternaryUnrolled;
		ternaryUnrolled.apply(sz, pLhs, MID, RHS, pOut, oper);
	}

	static inline void apply_4(int sz, INS_VEC LHS ,Float* pMid,  INS_VEC RHS, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< TernaryOpElement, INS_VEC, OP>  ternaryUnrolled;
		ternaryUnrolled.apply(sz, LHS, pMid, RHS, pOut, oper);
	}

	static inline void apply_4(int sz, INS_VEC LHS, INS_VEC MID, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< TernaryOpElement, INS_VEC, OP>  ternaryUnrolled;
		ternaryUnrolled.apply(sz, LHS, MID, pRhs, pOut, oper);
	}

	static inline void apply_1(int sz, Float* pLhs, Float* pMid, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_1< TernaryOpElement, INS_VEC, OP>  ternaryUnrolled;
		ternaryUnrolled.apply(sz, pLhs, pMid, pRhs, pOut, oper);
	}

	static inline void apply_1(int sz, INS_VEC LHS, Float* pMid, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_1< TernaryOpElement, INS_VEC, OP>  ternaryUnrolled;
		ternaryUnrolled.apply(sz, LHS, pMid, pRhs, pOut, oper);
	}

	static inline void apply_1(int sz, Float* pLhs, INS_VEC MID, Float* pRhs, Float* pOut, OP oper)
	{
		UnrollBinaryBody_1< TernaryOpElement, INS_VEC, OP>  ternaryUnrolled;
		ternaryUnrolled.apply(sz, pLhs, MID, pRhs, pOut, oper);
	}

	static inline void apply_1(int sz, Float* pLhs, Float* pMid, INS_VEC RHS, Float* pOut, OP oper)
	{
		UnrollBinaryBody_1< TernaryOpElement, INS_VEC, OP>  ternaryUnrolled;
		ternaryUnrolled.apply(sz, pLhs, pMid, RHS, pOut, oper);
	}

};

//////////////select unroll /////////////////
template<typename INS_VEC, typename OP, typename TRUE_OP = UnUsed, typename FALSE_OP = UnUsed >
struct Unroll_Select
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;

	//takes boolean values in pTest and switches PTrueVals and pFalseVals appropriately
	static inline void apply_4(int sz, Float* pTest, Float* pTrueVals, Float* pFalseVals, Float* pOut)
	{
		UnrollBinaryBody_4< SelectOpElement, INS_VEC, OP>  selectUnrolled;
		selectUnrolled.apply_4_bool(sz, pTest, pTrueVals, pFalseVals, pOut);
	}

	static inline void apply_4(int sz, Float* pTest, Float* pTrueVals, INS_VEC FalseVal, Float* pOut)
	{
		UnrollBinaryBody_4< SelectOpElement, INS_VEC, OP>  selectUnrolled;
		selectUnrolled.apply_4_bool(sz, pTest, pTrueVals, FalseVal, pOut);
	}

	static inline void apply_4(int sz, Float* pTest, INS_VEC TrueVal, Float* pFalseVals, Float* pOut)
	{
		UnrollBinaryBody_4< SelectOpElement, INS_VEC, OP>  selectUnrolled;
		selectUnrolled.apply_4_bool(sz, pTest, TrueVal, pFalseVals, pOut);
	}
	

	static inline void apply_4(int sz, Float* pTest, INS_VEC TrueVal, INS_VEC FalseVal, Float* pOut)
	{
		UnrollBinaryBody_4< SelectOpElement, INS_VEC, OP>  selectUnrolled;
		selectUnrolled.apply_4_bool(sz, pTest, TrueVal, FalseVal, pOut);
	}
	

	static inline void apply_4(int sz, Float* pTest, Float* pTrueVals, Float* pFalseVals, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< SelectOpElement, INS_VEC, OP>  selectUnrolled;
		selectUnrolled.apply(sz, pTest, pTrueVals, pFalseVals, pOut, oper); 
	}


	static inline void apply_4(int sz, Float* pTest, Float true_scalar, Float false_scalar, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< SelectOpElement, INS_VEC, OP>  selectUnrolled;
		INS_VEC trueVals(true_scalar);
		INS_VEC falseVals(false_scalar);
		selectUnrolled.apply(sz, pTest, trueVals, falseVals, pOut, oper);
	}


	static inline void apply_4(int sz, Float* pTest, TRUE_OP& true_oper, FALSE_OP& false_oper, Float* pOut, OP oper)
	{
		UnrollBinaryBody_4< SelectOpElement, INS_VEC, OP>  selectUnrolled;
		selectUnrolled.apply_f(sz, pTest, true_oper, false_oper, pOut, oper);
	}

};


//NB apply functions return a float
template<typename INS_VEC, typename OP>
struct Unroll_Accumulate
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;


	static inline Float apply_1(int sz, Float* pRhs,  OP oper, Float initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
	{
		UnrollUnitaryBody_1< AccumulateOpElement, INS_VEC, OP>  accumulateUnrolled(singularInit, initVal);
		int i = 0;
		accumulateUnrolled.apply(sz, pRhs, oper,i);
		return accumulateAcrossLanesAndEnd( sz,i, singularInit,  pRhs, oper, initVal, accumulateUnrolled.elem_0);
	}


	static inline Float apply_4(int sz, Float* pRhs, OP oper, Float initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
	{
		UnrollUnitaryBody_4< AccumulateOpElement, INS_VEC, OP>  accumulateUnrolled(singularInit, initVal);
		int i = 0;
		accumulateUnrolled.apply(sz, pRhs, oper,i);
		return accumulateAcrossLanesAndEnd(sz, i, singularInit, pRhs, oper, initVal, accumulateUnrolled.elem_0);
	}

	template <typename T>
	static inline Float accumulateAcrossLanesAndEnd(int sz, int i, bool singularInit, Float* pRhs, OP& oper, Float initVal, const T&  elem )
	{
		const int width = InstructionTraits<INS_VEC>::width;
		typename InstructionTraits<INS_VEC>::FloatType result = elem.elem_res.value[0];
		int min_wdth = std::min(sz, width);
		//across vectors lanes
		for (int j = 1; j < min_wdth; ++j)
		{
			result = ApplyBinaryOperationVec<INS_VEC, OP>(result, elem.elem_res.value[j], oper);
		}

		//odd bit at end
		for (; i < sz; ++i)
		{
			result = ApplyBinaryOperationVec<INS_VEC, OP>(result, pRhs[i], oper);
		}


		if (!singularInit)
		{
			result = ApplyBinaryOperationVec<INS_VEC, OP>(result, initVal, oper);
		}

		return result;
	}

};


//NB apply functions return a float
template<typename INS_VEC, typename OP, typename OP_TFORM>
struct Unroll_TransformAccumulate
{
	using Float = typename InstructionTraits< INS_VEC>::FloatType;


	static inline Float apply_1(int sz, Float* pRhs, OP oper, OP_TFORM operTransform,Float initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
	{
		UnrollUnitaryBody_1< AccumulateTransformOpElement, INS_VEC, OP>  accumulateTransformUnrolled(singularInit, initVal);
		int i = 0;
		accumulateTransformUnrolled.apply(sz, pRhs, oper, operTransform,i);
		return accumulateTransformAcrossLanesAndEnd(sz, i, singularInit, pRhs, oper, operTransform, initVal, accumulateTransformUnrolled.elem_0);

	}

	static inline Float apply_4(int sz, Float* pRhs, OP oper, OP_TFORM operTransform, Float initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
	{
		UnrollUnitaryBody_4< AccumulateTransformOpElement, INS_VEC, OP>  accumulateTransformUnrolled(singularInit, initVal);
		int i = 0;
		accumulateTransformUnrolled.apply(sz, pRhs, oper, operTransform, i);
		return accumulateTransformAcrossLanesAndEnd(sz, i, singularInit, pRhs, oper, operTransform, initVal, accumulateTransformUnrolled.elem_0);
	}



	// binary transform
	static inline Float apply_4(int sz, Float* pLhs, Float* pRhs, OP oper, OP_TFORM operTransform, Float initVal = InstructionTraits<INS_VEC>::nullValue, bool singularInit = true)
	{
		UnrollBinaryBody_4< AccumulateTransformBinaryOpElement, INS_VEC, OP > accumulateTransformUnrolled(singularInit, initVal);
		int i = 0;
		accumulateTransformUnrolled.apply(sz, pLhs, pRhs, oper, operTransform, i);
		return accumulateBinTransformAcrossLanesAndEnd(sz, i, singularInit, pLhs,pRhs, oper, operTransform, initVal, accumulateTransformUnrolled.elem_0);
	}


	template<typename T>
	static inline Float accumulateBinTransformAcrossLanesAndEnd(int sz, int i, bool singularInit, Float* pLhs, Float* pRhs, OP& oper, OP_TFORM& operTransform, Float initVal, const T& elem)
	{
		const int width = InstructionTraits<INS_VEC>::width;
		int min_wdth = std::min(sz, width);

		auto result = ( i>0 )?elem.elem_res.value[0]: (singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue);

		for (int j = 1; j < min_wdth; ++j)
		{
			result = ApplyBinaryOperationVec<INS_VEC, OP>(result, elem.elem_res.value[j], oper); //oper is accumulate oper
		}

		//odd bit at end
		for (; i < sz; ++i)
		{
			INS_VEC transformed = operTransform(pLhs[i],pRhs[i]);
			result = ApplyBinaryOperationVec<INS_VEC, OP>(transformed[0], result, oper); //oper is accumulate oper
		}


		if (!singularInit)
		{
			//INS_VEC transformedInitialVal = operTransform(initVal);
			//with binary operations we cant have un transformed single value ,
			// so we take the initial value as if it were transformed so that it can then be used in
			// accumulation
			INS_VEC transformedInitialVal = initVal;
			result = ApplyBinaryOperationVec<INS_VEC, OP>(result, transformedInitialVal[0], oper);
		}

		return result;

	}

	template <typename T>
	static inline Float accumulateTransformAcrossLanesAndEnd(int sz, int i, bool singularInit, Float* pRhs, OP& oper, OP_TFORM& operTransform, Float initVal, const T& elem)
	{
		// when sz < width ?? need to make sure init properly
		typename InstructionTraits<INS_VEC>::FloatType result = elem.elem_res.value[0];
		const int width = InstructionTraits<INS_VEC>::width;
		int min_wdth = std::min(sz, width);
		//accumulate across vectors lanes
		if (i > 0)
		{
			for (int j = 1; j < min_wdth; ++j)
			{
				result = ApplyBinaryOperationVec<INS_VEC, OP>(result, elem.elem_res.value[j], oper);
			}
		}
		else
		{
			result = singularInit ? initVal : InstructionTraits<INS_VEC>::nullValue;
		}

		//odd bit at end
		for (; i < sz; ++i)
		{
			INS_VEC transformed = operTransform(pRhs[i]);//
			result = ApplyBinaryOperationVec<INS_VEC, OP>(transformed[0], result, oper);
		}


		if (!singularInit)
		{
			INS_VEC transformedInitialVal = operTransform(initVal);
			result = ApplyBinaryOperationVec<INS_VEC, OP>(result, transformedInitialVal[0], oper);
		}

		return result;
	}




};
