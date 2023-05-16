#pragma once
#include "vec.h"
#include "vec_bool.h"
#include "instruction_traits.h"
#include "binary_unitary_operations.h"
#include "unroll_operators.h"
#include "error_utils.h"

#include <stdexcept>
#include <tuple>


template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyAdjacentDiff(const Vec<INS_VEC>& rhs1, OP& oper)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		throw std::range_error("ApplyAdjacentDiff called on scalar");
	}


	Vec<INS_VEC> result(rhs1.size());

	auto pRes = result.start();
	auto pRhs1 = rhs1.start();

	const int width = InstructionTraits<INS_VEC>::width;
	int step = 1 * width;


	INS_VEC RHS1;
	INS_VEC RHS2;
	INS_VEC RES;

	int sz = rhs1.paddedSize();

	int i = 0;
	for (; i < (sz - step); i += step)
	{
		RHS1.load_a(pRhs1 + i);
		RHS2.load(pRhs1 + i + 1); // load un-alligned
		RES = oper(RHS1, RHS2);
		RES.store_a(pRes + i);
	}
	for (int j = i; j < (sz - 1); j++)
	{
		result[j] = ApplyBinaryOperation1<INS_VEC, OP>(rhs1[j], rhs1[j + 1], oper);
	}

	return result;

};






template<typename INS_VEC, typename OP>
auto scan16(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC runValue, OP& oper)
{
	INS_VEC RHS;
	RHS.load_a(pRhs);

	INS_VEC RHS1 = blend16<16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14>(RHS, runValue);

	auto sum_1 = oper(RHS, RHS1);

	INS_VEC RHS2 = blend16<16, 16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13>(sum_1, runValue);

	auto sum_2 = oper(RHS2, sum_1);

	INS_VEC RHS3 = blend16<16, 16, 16, 16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11>(sum_2, runValue);

	auto sum_3 = oper(RHS3, sum_2);

	INS_VEC RHS4 = blend16<16, 16, 16, 16, 16, 16, 16, 16, 0, 1, 2, 3, 4, 5, 6, 7>(sum_3, runValue);

	auto sum_4 = oper(RHS4, sum_3);

	return sum_4;

}


template<typename INS_VEC, typename OP>

auto scan8(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC runValue, OP& oper)
{
	INS_VEC RHS;
	RHS.load_a(pRhs);
	INS_VEC RHS1 = blend8<8, 0, 1, 2, 3, 4, 5, 6>(RHS, runValue);

	auto sum_1 = oper(RHS, RHS1);

	INS_VEC RHS2 = blend8<8, 8, 0, 1, 2, 3, 4, 5>(sum_1, runValue);

	auto sum_2 = oper(RHS2, sum_1);

	INS_VEC RHS3 = blend8<8, 8, 8, 8, 0, 1, 2, 3>(sum_2, runValue);

	auto sum_3 = oper(RHS3, sum_2);

	return sum_3;

}


template<typename INS_VEC, typename OP>
auto scan4(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC runValue, OP& oper)
{
	INS_VEC RHS;
	RHS.load_a(pRhs);

	INS_VEC RHS1 = blend4<4, 0, 1, 2>(RHS, runValue);

	auto sum_1 = oper(RHS, RHS1);

	INS_VEC RHS2 = blend4<4, 4, 0, 1>(sum_1, runValue);

	return  oper(RHS2, sum_1);
}




template<typename INS_VEC, typename OP>
auto scan2(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC runValue, OP& oper)
{
	INS_VEC RHS;
	RHS.load_a(pRhs);
	INS_VEC RHS1 = blend2<2, 0>(RHS, runValue);
	return  oper(RHS, RHS1);
	
}




template<typename INS_VEC, typename OP>
auto scanN(typename InstructionTraits<INS_VEC>::FloatType*  pRhs, INS_VEC runValue, OP& oper)
{
	if constexpr (InstructionTraits<INS_VEC>::width == 8)
	{
		return scan8(pRhs, runValue, oper);
	}
	else if constexpr (InstructionTraits<INS_VEC>::width == 4)
	{
		return scan4(pRhs, runValue, oper);
	}
	else if constexpr (InstructionTraits<INS_VEC>::width == 16)
	{
		return scan16(pRhs, runValue, oper);
	}
	else  if constexpr (InstructionTraits<INS_VEC>::width == 2 )
	{
		//static_assert(false);
		return scan2(pRhs, runValue, oper);
	}
}




template< typename INS_VEC, typename OP>
Vec<INS_VEC> ApplyScan(const Vec<INS_VEC>& rhs1, OP& oper)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		throw std::range_error("ApplyScan called on scalar");
	}

	Vec<INS_VEC> result(rhs1.size());

	auto pRes = result.start();
	auto pRhs1 = rhs1.start();

	constexpr int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;


	INS_VEC runValue = 0.0;

	INS_VEC contValue = 0.0;

	int sz = rhs1.size();


	//https://gfxcourses.stanford.edu/cs149/fall20content/media/dataparallel/08_dataparallel.pdf

	constexpr int LAST_ELEM = width - 1;

	


	int i = 0;
	for (; i < (sz - step); i += step)
	{

		INS_VEC Res1 = scanN(pRhs1 + i, runValue, oper);
		INS_VEC Res2 = scanN(pRhs1 + i + width, runValue, oper);
		INS_VEC Res3 = scanN(pRhs1 + i + 2 * width, runValue, oper);
		INS_VEC Res4 = scanN(pRhs1 + i+ 3 * width, runValue, oper);

		INS_VEC runValue1 = Res1[LAST_ELEM];
		INS_VEC runValue2 = Res2[LAST_ELEM];
		INS_VEC runValue3 = Res3[LAST_ELEM];
		INS_VEC runValue4 = Res4[LAST_ELEM];

		auto sv12 = oper(runValue1, runValue2);
		auto sv13 = oper(sv12, runValue3);
		auto s14 = oper(sv13,runValue4);

		//Res1 += contValue;
		Res1 =oper(Res1, contValue);
		Res2 = oper(Res2, oper(runValue1, contValue));
		Res3 = oper(Res3, oper(sv12, contValue));
		Res4 = oper(Res4, oper(sv13, contValue));

		//runValue += s14;
		//contValue += s14;

		contValue = oper(contValue, s14);

		Res1.store_a(pRes + i);
		Res2.store_a(pRes + i + width);
		Res3.store_a(pRes + i + 2 * width);
		Res4.store_a(pRes + i + 3 * width);

	}

	if (i == 0) // small and need to init first element
	{ 
		i++; 
		result[0] = rhs1[0];
	};
	for (int j = i; j < sz ; j++)
	{
		result[j] = ApplyBinaryOperation1<INS_VEC, OP>(rhs1[j], result[j - 1], oper);
	}

	return result;

};

/////////////////////TRANSFORM SCAN  BROKEN/////////////////

/*

template<typename INS_VEC, typename OP, typename TRANSFORM>
auto scan16(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC& runValue, OP& oper, typename TRANSFORM& trfrm)
{
	INS_VEC RHS;
	RHS.load_a(pRhs);
	RHS = trfrm(RHS);

	INS_VEC RHS1 = blend16<16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14>(RHS, runValue);
	auto sum_1 = oper(RHS, RHS1);

	INS_VEC RHS2 = blend16<16, 16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13>(sum_1, runValue);
	auto sum_2 = oper(RHS2, sum_1);

	INS_VEC RHS3 = blend16<16, 16, 16, 16, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11>(sum_2, runValue);
	auto sum_3 = oper(RHS3, sum_2);

	INS_VEC RHS4 = blend16<16, 16, 16, 16, 16, 16, 16, 16, 0, 1, 2, 3, 4, 5, 6, 7>(sum_3, runValue);
	auto sum_4 = oper(RHS4, sum_3);

	return sum_4;

}


template<typename INS_VEC, typename OP, typename TRANSFORM >
auto scan8(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC& runValue, OP& oper ,TRANSFORM& trfrm)
{
	INS_VEC RHS;
	RHS.load_a(pRhs);
	RHS = trfrm(RHS);

	INS_VEC RHS1 = blend8<8, 0, 1, 2, 3, 4, 5, 6>(RHS, runValue);
	auto sum_1 = oper(RHS, RHS1);

	INS_VEC RHS2 = blend8<8, 8, 0, 1, 2, 3, 4, 5>(sum_1, runValue);
	auto sum_2 = oper(RHS2, sum_1);

	INS_VEC RHS3 = blend8<8, 8, 8, 8, 0, 1, 2, 3>(sum_2, runValue);
	auto sum_3 = oper(RHS3, sum_2);

	return sum_3;

}


template<typename INS_VEC, typename OP, typename TRANSFORM >
auto scan4(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC& runValue, OP& oper,  TRANSFORM& trfrm)
{
	INS_VEC RHS;
	RHS = trfrm(RHS);

	RHS.load_a(pRhs);
	INS_VEC RHS1 = blend4<4, 0, 1, 2>(RHS, runValue);

	auto sum_1 = oper(RHS, RHS1);

	INS_VEC RHS2 = blend4<4, 4, 0, 1>(sum_1, runValue);

	return  oper(RHS2, sum_1);
}



template<typename INS_VEC, typename OP, typename TRANSFORM >
auto scanN(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC& runValue, OP& oper, typename TRANSFORM& trfrm)
{
	if constexpr (InstructionTraits<INS_VEC>::width == 8)
	{
		return scan8(pRhs, runValue, oper, trfrm);
	}
	else if constexpr (InstructionTraits<INS_VEC>::width == 4)
	{
		return scan4(pRhs, runValue, oper, trfrm);
	}
	else if constexpr (InstructionTraits<INS_VEC>::width == 16)
	{
		return scan16(pRhs, runValue, oper, trfrm);
	}
	else
	{
		//static_assert(false);
	}
}


// works experimental
// generic scan
template< typename INS_VEC, typename OP, typename TRANSFORM>
Vec<INS_VEC> ApplyTransformScan(const Vec<INS_VEC>& rhs1, OP& oper, TRANSFORM& trfrm)
{
	check_vector(rhs1);
	if (isScalar(rhs1))
	{
		throw std::range_error("ApplyScan called on scalar");
	}

	Vec<INS_VEC> result(rhs1.size());

	auto pRes = result.start();
	auto pRhs1 = rhs1.start();

	constexpr int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;


	INS_VEC runValue = 0.0;

	int sz = rhs1.paddedSize();


	//https://gfxcourses.stanford.edu/cs149/fall20content/media/dataparallel/08_dataparallel.pdf

	constexpr int LAST_ELEM = width - 1;


	int i = 0;
	for (; i < (sz - step); i += step)
	{

		INS_VEC Res1 = scanN(pRhs1, runValue, oper, trfrm);
		INS_VEC Res2 = scanN(pRhs1 + width, runValue, oper, trfrm);
		INS_VEC Res3 = scanN(pRhs1 + 2 * width, runValue, oper, trfrm);
		INS_VEC Res4 = scanN(pRhs1 + 3 * width, runValue, oper, trfrm);

		INS_VEC runValue1 = Res1[LAST_ELEM];
		INS_VEC runValue2 = Res2[LAST_ELEM];
		INS_VEC runValue3 = Res3[LAST_ELEM];
		INS_VEC runValue4 = Res4[LAST_ELEM];

		auto sv12 = runValue1 + runValue2;
		auto sv13 = runValue1 + runValue2 + runValue3;
		auto sv34 = runValue3 + runValue4;

		Res2 += runValue1;
		Res3 += sv12;
		Res4 += sv13;

		runValue += sv34 + sv12;

		Res1.store_a(pRes + i);
		Res2.store_a(pRes + i + width);
		Res3.store_a(pRes + i + 2 * width);
		Res4.store_a(pRes + i + 3 * width);


	}
	for (int j = i; j < (sz - 1); j++)
	{
		result[j] = ApplyBinaryOperation1<INS_VEC, OP>(rhs1[j], rhs1[j + 1], oper);
	}

	return result;

};

*/
