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
auto scan16( INS_VEC RHS, INS_VEC runValue, OP& oper)
{

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
auto scan16(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC runValue, OP& oper)
{
	INS_VEC RHS;
	RHS.load_a(pRhs);
	return scan16(RHS, runValue,  oper);

}




template<typename INS_VEC, typename OP>
auto scan8( INS_VEC RHS, INS_VEC runValue, OP& oper)
{

	INS_VEC RHS1 = blend8<8, 0, 1, 2, 3, 4, 5, 6>(RHS, runValue);

	auto sum_1 = oper(RHS, RHS1);

	INS_VEC RHS2 = blend8<8, 8, 0, 1, 2, 3, 4, 5>(sum_1, runValue);

	auto sum_2 = oper(RHS2, sum_1);

	INS_VEC RHS3 = blend8<8, 8, 8, 8, 0, 1, 2, 3>(sum_2, runValue);

	auto sum_3 = oper(RHS3, sum_2);

	return sum_3;

}


template<typename INS_VEC, typename OP>
auto scan8(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC runValue, OP& oper)
{
	INS_VEC RHS;
	RHS.load_a(pRhs);
	return scan8(RHS, runValue, oper);
}



template<typename INS_VEC, typename OP>
auto scan4( INS_VEC RHS, INS_VEC runValue, OP& oper)
{
	INS_VEC RHS1 = blend4<4, 0, 1, 2>(RHS, runValue);
	auto sum_1 = oper(RHS, RHS1);
	INS_VEC RHS2 = blend4<4, 4, 0, 1>(sum_1, runValue);
	return  oper(RHS2, sum_1);
}



template<typename INS_VEC, typename OP>
auto scan4(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC runValue, OP& oper)
{
	INS_VEC RHS;
	RHS.load_a(pRhs);
	return scan4(RHS, runValue, oper);
}




template<typename INS_VEC, typename OP>
auto scan2( INS_VEC RHS, INS_VEC runValue, OP& oper)
{
	INS_VEC RHS1 = blend2<2, 0>(RHS, runValue);
	return  oper(RHS, RHS1);
}



template<typename INS_VEC, typename OP>
auto scan2(typename InstructionTraits<INS_VEC>::FloatType* pRhs, INS_VEC runValue, OP& oper)
{
	INS_VEC RHS;
	RHS.load_a(pRhs);
	return scan2(RHS, runValue, oper);
	
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
		return scan2(pRhs, runValue, oper);
	}
}




template<typename INS_VEC, typename OP>
auto scanN(  INS_VEC RHS, INS_VEC runValue, OP& oper)
{
	if constexpr (InstructionTraits<INS_VEC>::width == 8)
	{
		return scan8(RHS, runValue, oper);
	}
	else if constexpr (InstructionTraits<INS_VEC>::width == 4)
	{
		return scan4(RHS, runValue, oper);
	}
	else if constexpr (InstructionTraits<INS_VEC>::width == 16)
	{
		return scan16(RHS, runValue, oper);
	}
	else  if constexpr (InstructionTraits<INS_VEC>::width == 2)
	{
		return scan2(RHS, runValue, oper);
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

		Res1 =oper(Res1, contValue);
		Res2 = oper(Res2, oper(runValue1, contValue));
		Res3 = oper(Res3, oper(sv12, contValue));
		Res4 = oper(Res4, oper(sv13, contValue));


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

/////////////////////TRANSFORM SCAN /////////////////
// generic transform scan
template< typename INS_VEC, typename OP, typename TRANSFORM>
Vec<INS_VEC> ApplyTransformScan(const Vec<INS_VEC>& rhs1, OP& oper, TRANSFORM& transform)
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

	INS_VEC Res1;
	INS_VEC Res2;
	INS_VEC Res3;
	INS_VEC Res4;


	int i = 0;
	for (; i < (sz - step); i += step)
	{

		Res1.load_a(pRhs1 + i);
		Res2.load_a(pRhs1 + i + width);
		Res3.load_a(pRhs1 + i + 2 * width);
		Res4.load_a(pRhs1 + i + 3 * width);


		Res1 = scanN(transform(Res1), runValue, oper);
		Res2 = scanN(transform(Res2), runValue, oper);
		Res3 = scanN(transform(Res3), runValue, oper);
		Res4 = scanN(transform(Res4), runValue, oper);

		INS_VEC runValue1 = Res1[LAST_ELEM];
		INS_VEC runValue2 = Res2[LAST_ELEM];
		INS_VEC runValue3 = Res3[LAST_ELEM];
		INS_VEC runValue4 = Res4[LAST_ELEM];

		auto sv12 = oper(runValue1, runValue2);
		auto sv13 = oper(sv12, runValue3);
		auto s14 = oper(sv13, runValue4);

		Res1 = oper(Res1, contValue);
		Res2 = oper(Res2, oper(runValue1, contValue));
		Res3 = oper(Res3, oper(sv12, contValue));
		Res4 = oper(Res4, oper(sv13, contValue));


		contValue = oper(contValue, s14);

		Res1.store_a(pRes + i);
		Res2.store_a(pRes + i + width);
		Res3.store_a(pRes + i + 2 * width);
		Res4.store_a(pRes + i + 3 * width);

	}

	if (i == 0) // small and need to init first element
	{
		i++;
		result[0] = transform(INS_VEC(rhs1[0]))[0];
	};
	for (int j = i; j < sz; j++)
	{
		auto trans_of_rhs1_j = transform(INS_VEC(rhs1[j]))[0];

		result[j] = ApplyBinaryOperation1<INS_VEC, OP>(trans_of_rhs1_j, result[j - 1], oper);
	}

	return result;

};


template< typename INS_VEC, typename OP, typename TRANSFORM>
Vec<INS_VEC> ApplyTransformScan(const Vec<INS_VEC>& lhs1, const Vec<INS_VEC>& rhs1, OP& oper, TRANSFORM& transform)
{


	check_vector(rhs1);
	check_vector(lhs1);
	if (isScalar(rhs1))
	{
		throw std::range_error("ApplyScan called on scalar");
	}

	Vec<INS_VEC> result(rhs1.size());

	auto pRes = result.start();
	auto pLhs1 = lhs1.start();
	auto pRhs1 = rhs1.start();

	constexpr int width = InstructionTraits<INS_VEC>::width;
	int step = 4 * width;


	INS_VEC runValue = 0.0;

	INS_VEC contValue = 0.0;

	int sz = rhs1.size();


	//https://gfxcourses.stanford.edu/cs149/fall20content/media/dataparallel/08_dataparallel.pdf

	constexpr int LAST_ELEM = width - 1;

	INS_VEC LHS1;
	INS_VEC LHS2;
	INS_VEC LHS3;
	INS_VEC LHS4;

	INS_VEC RHS1;
	INS_VEC RHS2; 
	INS_VEC RHS3;
	INS_VEC RHS4;


	int i = 0;
	for (; i < (sz - step); i += step)
	{



		RHS1.load_a(pRhs1 + i);
		RHS2.load_a(pRhs1 + i + width);
		RHS3.load_a(pRhs1 + i + 2 * width);
		RHS4.load_a(pRhs1 + i + 3 * width);

		LHS1.load_a(pLhs1 + i);
		LHS2.load_a(pLhs1 + i + width);
		LHS3.load_a(pLhs1 + i + 2 * width);
		LHS4.load_a(pLhs1 + i + 3 * width);

		auto Res1 = scanN(transform(LHS1,RHS1), runValue, oper);
		auto Res2 = scanN(transform(LHS2,RHS2), runValue, oper);
		auto Res3 = scanN(transform(LHS3, RHS3), runValue, oper);
		auto Res4 = scanN(transform(LHS4, RHS4), runValue, oper);

		INS_VEC runValue1 = Res1[LAST_ELEM];
		INS_VEC runValue2 = Res2[LAST_ELEM];
		INS_VEC runValue3 = Res3[LAST_ELEM];
		INS_VEC runValue4 = Res4[LAST_ELEM];

		auto sv12 = oper(runValue1, runValue2);
		auto sv13 = oper(sv12, runValue3);
		auto s14 = oper(sv13, runValue4);

		Res1 = oper(Res1, contValue);
		Res2 = oper(Res2, oper(runValue1, contValue));
		Res3 = oper(Res3, oper(sv12, contValue));
		Res4 = oper(Res4, oper(sv13, contValue));


		contValue = oper(contValue, s14);

		Res1.store_a(pRes + i);
		Res2.store_a(pRes + i + width);
		Res3.store_a(pRes + i + 2 * width);
		Res4.store_a(pRes + i + 3 * width);

	}

	if (i == 0) // small and need to init first element
	{
		i++;
		result[0] = transform(INS_VEC(lhs1[0]),INS_VEC(rhs1[0]))[0];
	};
	for (int j = i; j < sz; j++)
	{
		auto trans_of_rhs1_j = transform(INS_VEC(lhs1[j]),INS_VEC(rhs1[j]))[0];

		result[j] = ApplyBinaryOperation1<INS_VEC, OP>(trans_of_rhs1_j, result[j - 1], oper);
	}

	return result;

};

