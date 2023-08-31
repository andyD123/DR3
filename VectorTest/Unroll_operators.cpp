#include "pch.h"

#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/vec_bool_d.h"
#include "../Vectorisation/VecX/vec_double.h"
#include  "../Vectorisation/VecX/alloc_policy.h"

#include "../Vectorisation/VecX/vec_d.h"
#include "../Vectorisation/VecX/vec_bool.h"
#include "../Vectorisation/VecX/vec_view.h"

#include "../Vectorisation/VecX/target_name_space.h"
#include "../Vectorisation/VecX/unroll_operators.h"

#include "../Vectorisation/VecX/dr3.h"
#include "testNamespace.h"
#include "dr3TestUtil.h"

#include <numeric>



TEST(TestUnroll, RegisterElement_load) 
{
	const int SZ = 100;
	Numeric alignas(64)  input[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input[i] = static_cast<Numeric>(i);
	}


	RegisterElement<typename VecXX::INS, 0> test;
	test.load(input);

	auto res = test.value;

	constexpr int registerWidth = InstructionTraits<VecXX::INS>::width;	
	for (int j = 0; j < registerWidth; j++)
	{
		EXPECT_NUMERIC_EQ(res[j], input[j]);
	}


	RegisterElement<typename VecXX::INS, 1> test_Offset;
	test_Offset.load(input);

	res = test_Offset.value;

	for (int j = 0; j < registerWidth; j++)
	{
		EXPECT_NUMERIC_EQ(res[j], input[j+ registerWidth]);
	}


}

//assuming things are aligned 
TEST(TestUnroll, RegisterElement_Save)
{
	const int SZ = 100;
	Numeric  alignas(64) input[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input[i] = static_cast<Numeric>(i);
	}


	RegisterElement<typename VecXX::INS, 0> test_zero;
	 
	Numeric  alignas(64) saveValues[] = { asNumber(101.1),
		asNumber(102.1),asNumber(103.1),asNumber(104.1),asNumber(105.1),asNumber(106.1),asNumber(107.1),
		asNumber(108.1),asNumber(109.1),asNumber(110.1),asNumber(111.1),asNumber(112.1),asNumber(113.1),
		asNumber(114.1),asNumber(115.1),asNumber(116.1) };

	test_zero.value.load(saveValues);
	test_zero.store(input);

	constexpr int registerWidth = InstructionTraits<VecXX::INS>::width;
	for (int j = 0; j < registerWidth; j++)
	{
		EXPECT_NUMERIC_EQ(saveValues[j], input[j]);
	}


	RegisterElement<typename VecXX::INS, 3> test_three;
	Numeric  alignas(64) saveValuesAgain[] =	{ asNumber(1101.1),	asNumber(1102.1),asNumber(1103.1),asNumber(1104.1),asNumber(1105.1),
		asNumber(1106.1),asNumber(1107.1),	asNumber(1108.1),asNumber(1109.1),asNumber(1110.1),asNumber(1111.1),asNumber(1112.1),asNumber(1113.1),
	    asNumber(1114.1),asNumber(1115.1),asNumber(1116.1) };


	test_three.value.load(saveValuesAgain);
	test_three.store(input);

	for (int j = 0; j < registerWidth; j++)
	{
		EXPECT_NUMERIC_EQ(saveValuesAgain[j], input[j+ 3*registerWidth]);
	}


}



TEST(TestUnroll, BinaryOpElement_apply)
{
	const int SZ = 100;
	Numeric  alignas(64) input_A[SZ];
	Numeric  alignas(64) input_B[SZ];
	Numeric  alignas(64) OUT[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<Numeric>(i);
		input_B[i] = static_cast<Numeric>(100 -i);
		OUT[i] = 0.0;
	}

	auto addLambda = [](auto x, auto y) { return x + y; };
	int sz = SZ;
	int wdth =InstructionTraits<VecXX::INS>::width;

	for (int sz = 0; sz < SZ; sz+= wdth)
	{
		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_1(sz, input_A, input_B, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_NUMERIC_EQ(OUT[i], input_A[i] + input_B[i]);
		}
	}


	for (int sz = 0; sz < SZ; sz += wdth)
	{
		VecXX::INS fixedValue = 23;
		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_1(sz, fixedValue, input_B, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_NUMERIC_EQ(OUT[i], fixedValue[0] + input_B[i]);
		}
	}


	for (int sz = 0; sz < SZ; sz += wdth)
	{
		VecXX::INS fixedValue = 23;
		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_1(sz, input_A, fixedValue, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_NUMERIC_EQ(OUT[i], input_A[i] + fixedValue[0]);
		}
	}
}




TEST(TestUnroll, BinaryOpElement_XXX_apply4)
{
	const int SZ = 100;
	Numeric  alignas(64) input_A[SZ];
	Numeric  alignas(64) input_B[SZ];
	Numeric  alignas(64) OUT[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<Numeric>(i);
		input_B[i] = static_cast<Numeric>(100 - i);
		OUT[i] = 0.0;
	}

	auto addLambda = [](auto x, auto y) { return x + y; };

	int sz = SZ;
	int wdth = InstructionTraits<VecXX::INS>::width;
	for (int sz = 0; sz < SZ; sz += wdth)
	{
		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_4(sz, input_A, input_B, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_NUMERIC_EQ(OUT[i], input_A[i] + input_B[i]);
		}
	}


	for (int sz = 0; sz < SZ; sz += wdth)
	{
		VecXX::INS fixedValue = 23;
		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_4(sz, fixedValue, input_B, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_NUMERIC_EQ(OUT[i], fixedValue[0] + input_B[i]);
		}
	}


	for (int sz = 0; sz < SZ; sz += wdth)
	{
		VecXX::INS fixedValue = 23;
		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_4(sz, input_A, fixedValue, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_NUMERIC_EQ(OUT[i], input_A[i] + fixedValue[0]);
		}
	}

}



TEST(TestUnroll, UnitaryOpElement_apply_1)
{
	const int SZ = 100;
	Numeric  alignas(64) input_A[SZ];
	Numeric  alignas(64) OUT[SZ];

	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<Numeric>(i);
		OUT[i] = 0.0;
	}

	auto squareLambda = [](auto x) { return x * x; };
	int sz = SZ;
	int wdth = InstructionTraits<VecXX::INS>::width;

	for (int sz = 0; sz < SZ; sz += wdth)
	{
		Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_1(sz, input_A,  OUT, squareLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_NUMERIC_EQ(OUT[i], input_A[i] * input_A[i]);
		}
	}


}


TEST(TestUnroll, UnitaryOpElement_apply_4)
{
	const int SZ = 100;
	Numeric  alignas(64) input_A[SZ];
	Numeric  alignas(64) OUT[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<Numeric>(i);
		OUT[i] = 0.0;
	}

	auto squareLambda = [](auto x) { return x * x; };
	int SZ_MAX = SZ;
	int wdth = InstructionTraits<VecXX::INS>::width;

	for (int sz = 0; sz <= SZ_MAX; sz += wdth)
	{
		Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_NUMERIC_EQ(OUT[i], input_A[i] * input_A[i]);
		}
	}
	
}


TEST(TestUnroll, UnitaryOpElement_apply_4_simple)
{
	const int SZ = 100;
	Numeric  alignas(64) input_A[SZ];
	Numeric  alignas(64) OUT[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<Numeric>(i);
		OUT[i] = 0.0;
	}

	auto squareLambda = [](auto x) { return x * x; };
	int SZ_MAX = SZ;
	int wdth = InstructionTraits<VecXX::INS>::width;
	int sz = 4;
	
	for (auto& x : OUT) { x = 0.0; };
	sz = wdth;
	Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
	for (int i = 0; i < sz; i++)
	{
		EXPECT_NUMERIC_EQ(OUT[i], input_A[i] * input_A[i]);
	}

	for (auto& x : OUT) { x = 0.0; };
	sz =2* wdth;
	Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
	for (int i = 0; i < sz; i++)
	{
		EXPECT_NUMERIC_EQ(OUT[i], input_A[i] * input_A[i]);
	}

	for (auto& x : OUT) { x = 0.0; };
	sz = 3 * wdth;
	Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
	for (int i = 0; i < sz; i++)
	{
		EXPECT_NUMERIC_EQ(OUT[i], input_A[i] * input_A[i]);
	}

	for (auto& x : OUT) { x = 0.0; };
	sz = 4* wdth;
	Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
	for (int i = 0; i < sz; i++)
	{
		EXPECT_NUMERIC_EQ(OUT[i], input_A[i] * input_A[i]);
	}

	for (auto& x : OUT) { x = 0.0; };
	sz = 5 * wdth;
	Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
	for (int i = 0; i < sz; i++)
	{
		EXPECT_NUMERIC_EQ(OUT[i], input_A[i] * input_A[i]);
	}

}



TEST(TestUnroll, UnitaryOpElement_filter)
{
	const int SZ = 100;
	Numeric  alignas(64) input_A[SZ];
	Numeric  alignas(64) OUT[SZ];
	unsigned int   alignas(16) idx[SZ];

	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<Numeric>(i);
		OUT[i] = 0.0;
		idx[i] = 000;
	}


	auto filter_Positive = [](auto x) { return x > asNumber(0.0); };
	int sz = SZ;

	int i = 4;
	int  psn = 0;
	unsigned int* pIdx =idx;

	UnitaryOpElement<typename VecXX::INS, typename decltype(filter_Positive), 0 > elem;
	elem.elem_lhs.load(input_A);
	elem.filter(input_A, OUT, i, filter_Positive, psn, pIdx);
	constexpr int registerWidth = InstructionTraits<VecXX::INS>::width;

	for (int j = 0; j < registerWidth-1; j++) // first element is false so we do wdth -1 tests
	{
		EXPECT_NUMERIC_EQ(OUT[j], input_A[i+1+j]);
		EXPECT_EQ(pIdx[j], i + 1+j);
	}


	{
		int i = 0;
		int  psn = 0;
		UnitaryOpElement<typename VecXX::INS, typename decltype(filter_Positive), 0 > elem;
		elem.elem_lhs.load(input_A);
		elem.filter(input_A, OUT, i, filter_Positive, psn, pIdx);
		for (int j = 0; j < registerWidth-1 ; j++) // first element is false so we do wdth -1 tests
		{
			EXPECT_NUMERIC_EQ(OUT[j], input_A[i + j+1]);
			EXPECT_EQ(pIdx[j], 1 + j);
		}
	}


	{
		auto filter_GT_One = [](auto x) { return x > asNumber(1.0); };
		int i = 0;
		int  psn = 0;
		UnitaryOpElement<typename VecXX::INS, typename decltype(filter_GT_One), 0 > elem;
		elem.elem_lhs.load(input_A);
		elem.filter(input_A, OUT, i, filter_GT_One, psn, pIdx);
		elem.filter(input_A, OUT, i, filter_Positive, psn, pIdx);

		for (int j = 0; j < registerWidth-2; j++) // first and 2nd element is false so we do wdth -2 tests
		{
			EXPECT_NUMERIC_EQ(OUT[j], input_A[i + j+2]);
			EXPECT_EQ(pIdx[j], i + j+2);
		}
	}
}

