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


#include <numeric>

//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;



TEST(TestUnroll, RegisterElement_load) 
{
	const int SZ = 100;
	double alignas(16)  input[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input[i] = static_cast<double>(i);
	}


	RegisterElement<typename VecXX::INS, 0> test;
	test.load(input);

	auto res = test.value;

	EXPECT_EQ(res[0], input[0]);
	EXPECT_EQ(res[1], input[1]);
	EXPECT_EQ(res[2], input[2]);
	EXPECT_EQ(res[3], input[3]);


	RegisterElement<typename VecXX::INS, 1> test_Offset;
	test_Offset.load(input);

	res = test_Offset.value;

	EXPECT_EQ(res[0], input[4]);
	EXPECT_EQ(res[1], input[5]);
	EXPECT_EQ(res[2], input[6]);
	EXPECT_EQ(res[3], input[7]);


}

//assuming things are aligned 
TEST(TestUnroll, RegisterElement_Save)
{
	const int SZ = 100;
	double  alignas(16) input[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input[i] = static_cast<double>(i);
	}


	RegisterElement<typename VecXX::INS, 0> test_zero;
	 
	double alignas(16) saveValues[] = { 101.1,102.1,103.1,104.1 };

	test_zero.value.load(saveValues);
	test_zero.store(input);



	EXPECT_EQ(saveValues[0], input[0]);
	EXPECT_EQ(saveValues[1], input[1]);
	EXPECT_EQ(saveValues[2], input[2]);
	EXPECT_EQ(saveValues[3], input[3]);


	RegisterElement<typename VecXX::INS, 3> test_three;
	double alignas(16) saveValuesAgain[] = { 1101.1,1102.1,1103.1,1104.1 };
	test_three.value.load(saveValuesAgain);
	test_three.store(input);


	EXPECT_EQ(saveValuesAgain[0], input[12]);
	EXPECT_EQ(saveValuesAgain[1], input[13]);
	EXPECT_EQ(saveValuesAgain[2], input[14]);
	EXPECT_EQ(saveValuesAgain[3], input[15]);


}



TEST(TestUnroll, BinaryOpElement_apply)
{
	const int SZ = 100;
	double  alignas(16) input_A[SZ];
	double  alignas(16) input_B[SZ];
	double  alignas(16) OUT[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<double>(i);

		input_B[i] = static_cast<double>(100 -i);

		OUT[i] = 0.0;
	}

	auto addLambda = [](auto x, auto y) { return x + y; };


	int sz = SZ;

	for (int sz = 4; sz < SZ; sz+=4)
	{
		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_1(sz, input_A, input_B, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_EQ(OUT[i], input_A[i] + input_B[i]);
		}
	}


	for (int sz = 4; sz < SZ; sz += 4)
	{
		VecXX::INS fixedValue = 23;

		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_1(sz, fixedValue, input_B, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_EQ(OUT[i], fixedValue[0] + input_B[i]);
		}
	}


	for (int sz = 4; sz < SZ; sz += 4)
	{
		VecXX::INS fixedValue = 23;

		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_1(sz, input_A, fixedValue, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_EQ(OUT[i], input_A[i] + fixedValue[0]);
		}
	}



}




TEST(TestUnroll, BinaryOpElement_XXX_apply4)
{
	const int SZ = 100;
	double  alignas(16) input_A[SZ];
	double  alignas(16) input_B[SZ];
	double  alignas(16) OUT[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<double>(i);

		input_B[i] = static_cast<double>(100 - i);

		OUT[i] = 0.0;
	}

	auto addLambda = [](auto x, auto y) { return x + y; };


	int sz = SZ;

	for (int sz = 4; sz < SZ; sz += 4)
	{
		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_4(sz, input_A, input_B, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_EQ(OUT[i], input_A[i] + input_B[i]);
		}
	}


	for (int sz = 4; sz < SZ; sz += 4)
	{
		VecXX::INS fixedValue = 23;

		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_4(sz, fixedValue, input_B, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_EQ(OUT[i], fixedValue[0] + input_B[i]);
		}
	}


	for (int sz = 4; sz < SZ; sz += 4)
	{
		VecXX::INS fixedValue = 23;

		Unroll_Binary<typename VecXX::INS, typename decltype(addLambda)>::apply_4(sz, input_A, fixedValue, OUT, addLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_EQ(OUT[i], input_A[i] + fixedValue[0]);
		}
	}



}



TEST(TestUnroll, UnitaryOpElement_apply_1)
{
	const int SZ = 100;
	double  alignas(16) input_A[SZ];
	//double  alignas(16) input_B[SZ];
	double  alignas(16) OUT[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<double>(i);

		OUT[i] = 0.0;
	}

	auto squareLambda = [](auto x) { return x * x; };


	int sz = SZ;

	for (int sz = 4; sz < SZ; sz += 4)
	{
		Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_1(sz, input_A,  OUT, squareLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_EQ(OUT[i], input_A[i] * input_A[i]);
		}
	}


}


TEST(TestUnroll, UnitaryOpElement_apply_4)
{
	const int SZ = 100;
	double  alignas(16) input_A[SZ];
	double  alignas(16) OUT[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<double>(i);

		OUT[i] = 0.0;
	}

	auto squareLambda = [](auto x) { return x * x; };


	int SZ_MAX = SZ;

	for (int sz = 4; sz <= SZ_MAX; sz += 4)
	{
		Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_EQ(OUT[i], input_A[i] * input_A[i]);
		}
	}
	
}


TEST(TestUnroll, UnitaryOpElement_apply_4_simple)
{
	const int SZ = 100;
	double  alignas(16) input_A[SZ];
	double  alignas(16) OUT[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<double>(i);

		OUT[i] = 0.0;
	}

	auto squareLambda = [](auto x) { return x * x; };


	int SZ_MAX = SZ;

	//for (int sz = 4; sz <= SZ_MAX; sz += 4)
	//{
	int sz = 4;
	
	for (auto& x : OUT) { x = 0.0; };
	sz = 4;
	Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
	for (int i = 0; i < sz; i++)
	{
		EXPECT_EQ(OUT[i], input_A[i] * input_A[i]);
	}

	for (auto& x : OUT) { x = 0.0; };
	sz = 8;
	Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
	for (int i = 0; i < sz; i++)
	{
		EXPECT_EQ(OUT[i], input_A[i] * input_A[i]);
	}

	for (auto& x : OUT) { x = 0.0; };
	sz = 12;
	Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
	for (int i = 0; i < sz; i++)
	{
		EXPECT_EQ(OUT[i], input_A[i] * input_A[i]);
	}

	for (auto& x : OUT) { x = 0.0; };
	sz = 16;
	Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
	for (int i = 0; i < sz; i++)
	{
		EXPECT_EQ(OUT[i], input_A[i] * input_A[i]);
	}

	for (auto& x : OUT) { x = 0.0; };
	sz = 20;
	Unroll_Unitary<typename VecXX::INS, typename decltype(squareLambda)>::apply_4(sz, input_A, OUT, squareLambda);
	for (int i = 0; i < sz; i++)
	{
		EXPECT_EQ(OUT[i], input_A[i] * input_A[i]);
	}
	

}



TEST(TestUnroll, UnitaryOpElement_filter)
{
	const int SZ = 100;
	double  alignas(16) input_A[SZ];
	double  alignas(16) OUT[SZ];


	unsigned int   alignas(16) idx[SZ];


	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<double>(i);

		OUT[i] = 0.0;
		idx[i] = 000;
	}


	auto filter_Positive = [](auto x) { return x > 0; }; 


	int sz = SZ;

	int i = 4;
	int  psn = 0;
	unsigned int* pIdx =idx;

	UnitaryOpElement<typename VecXX::INS, typename decltype(filter_Positive), 0 > elem;
	elem.elem_lhs.load(input_A);

	
	elem.filter(input_A, OUT, i, filter_Positive, psn, pIdx);
		

	EXPECT_EQ(OUT[0], input_A[4+1]);
	EXPECT_EQ(OUT[1], input_A[4+2]);
	EXPECT_EQ(OUT[2], input_A[4+3]);
	EXPECT_EQ(pIdx[0], 4 + 1);
	EXPECT_EQ(pIdx[1], 4 + 2);
	EXPECT_EQ(pIdx[2], 4 + 3);

	{
		int i = 0;
		int  psn = 0;

		UnitaryOpElement<typename VecXX::INS, typename decltype(filter_Positive), 0 > elem;
		elem.elem_lhs.load(input_A);


		elem.filter(input_A, OUT, i, filter_Positive, psn, pIdx);

		EXPECT_EQ(OUT[0], input_A[1]);
		EXPECT_EQ(OUT[1], input_A[2]);
		EXPECT_EQ(OUT[2], input_A[3]);
		EXPECT_EQ(pIdx[0], 1);
		EXPECT_EQ(pIdx[1], 2);
		EXPECT_EQ(pIdx[2], 3);
	}


	{
		auto filter_GT_One = [](auto x) { return x > 1; };

		int i = 0;
		int  psn = 0;

		UnitaryOpElement<typename VecXX::INS, typename decltype(filter_GT_One), 0 > elem;
		elem.elem_lhs.load(input_A);


		elem.filter(input_A, OUT, i, filter_GT_One, psn, pIdx);
		EXPECT_EQ(OUT[0], input_A[2]);
		EXPECT_EQ(OUT[1], input_A[3]);
		EXPECT_EQ(pIdx[0], 2);
		EXPECT_EQ(pIdx[1], 3);

	}


}





/*
TEST(TestUnroll, Unitary_apply_4_filter)
{
	const int SZ = 100;
	double  alignas(16) input_A[SZ];
	double  alignas(16) OUT_FILT[SZ];
	//double  alignas(16) OUT[SZ];
	for (int i = 0; i < SZ; ++i)
	{
		input_A[i] = static_cast<double>(i - SZ / 2);

	//	OUT[i] = 0.0;
		OUT_FILT[i] = 0;
	}

	//auto squareLambda = [](auto x) { return x * x; };

	auto filter_Positive = [](auto x) { return x > 0; }; //bool return 


	int sz = SZ;

	for (int sz = 4; sz < SZ; sz += 4)
	{
		Unroll_Unitary_Filter<typename VecXX::INS, typename decltype(filter_Positive)>::apply_4(sz, input_A,  filter_Positive, OUT_FILT);
		for (int i = 0; i < sz; i++)
		{
			EXPECT_TRUE(OUT_FILT[i] > 0);
		}
	}


}

*/





/*


TEST(TestViews, ApplyUnitaryOperationX) {


	for (int SZ = 3; SZ < 123; SZ++)
	{

		std::vector<double>  v(SZ, 6.66);
		int i = 0;
		for (auto& x : v) { x -= 500.0 + i; ++i; }
		VecXX test(v);


		auto mySquareItLambda = [](const auto& x) {return x * x;  };

		auto testSquare = ApplyUnitaryOperationX(mySquareItLambda, test);

		for (int i = 0; i < test.size(); ++i)
		{
			EXPECT_EQ(testSquare[i], test[i] * test[i]);
		}
	}



}




TEST(TestViews, ApplyBinaryOperation) {


	for (int SZ = 3; SZ < 123; SZ++)
	{


		std::vector<double>  v(SZ, 6.66);
		int i = 0;
		for (auto& x : v) { x -= 500.0 + i; ++i; }
		VecXX test(v);
		auto test2 = test * 3.0;


		auto myDifItLambda = [](const auto& x, const auto& y) {return x - y;  };

		auto testDiff = ApplyBinaryOperation1(test, test2, myDifItLambda);


		EXPECT_EQ(testDiff[0], test[0] - test2[0]);
		EXPECT_EQ(testDiff[SZ - 1], test[SZ - 1] - test2[SZ - 1]);
		EXPECT_EQ(testDiff[SZ / 2], test[SZ / 2] - test2[SZ / 2]);


		for (int i = 1; i < test.size() - 1; ++i)
		{
			EXPECT_EQ(testDiff[i], test[i] - test2[i]);
		}
	}


}


TEST(TestViews, applyFilter1000)
{


	int SZ = 1000;

	std::vector<double>  v(SZ, 6.66);
	int i = 0;
	for (auto& x : v) { x -= 500.0 + i; ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	initVals += 3000.0;
	auto condTestVals = initVals <= 200.0;

	VecVW vw = ApplyFilter(condTestVals, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals = !condTestVals;

	VecVW vwNegate = ApplyFilter(NotCondTestVals, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE(vwNegate.size() + vw.size(), initVals.size());



}



TEST(TestViews, applyFilter3)
{


	int SZ = 3;

	std::vector<double>  v(SZ, 6.66);
	int i = 0;
	for (auto& x : v) { x -= 500.0 + i; ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	initVals += 3000.0;
	auto condTestVals = initVals <= 200.0;

	VecVW vw = ApplyFilter(condTestVals, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals = !condTestVals;

	VecVW vwNegate = ApplyFilter(NotCondTestVals, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE(vwNegate.size() + vw.size(), initVals.size());

}


TEST(TestViews, applyFilterRange4)
{

	int SZ = 4;

	std::vector<double>  v(SZ, 6.66);
	int i = 0;
	for (auto& x : v) { x -= 500.0 + i; ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	initVals += 3000.0;
	auto condTestVals = initVals <= 200.0;

	VecVW vw = ApplyFilter(condTestVals, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals = !condTestVals;

	VecVW vwNegate = ApplyFilter(NotCondTestVals, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE(vwNegate.size() + vw.size(), initVals.size());

}


void testFilterVec(int SZ)
{

	std::vector<double>  v(SZ, 6.66);
	int i = 0;
	for (auto& x : v) { x -= 500.0 + i; ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	initVals += 3000.0;
	auto condTestVals = initVals <= 200.0;

	VecVW vw = ApplyFilter(condTestVals, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals = !condTestVals;

	VecVW vwNegate = ApplyFilter(NotCondTestVals, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE(vwNegate.size() + vw.size(), initVals.size());
}



TEST(TestViews, applyFilterRangeAll)
{

	for (int SZ = 3; SZ < 153; ++SZ) //int bug
	{
		testFilterVec(SZ);
	}

}


TEST(TestViews, applyFilterSmall)
{


	std::vector<double>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };

	VecXX test(v);

	auto initVals = test;

	//boolean vec
	auto condTestVals = test >= 2.0;

	VecVW vw = ApplyFilter(condTestVals, test);


	EXPECT_EQ(vw[0], 2.0);
	EXPECT_EQ(vw[1], 3.0);
	EXPECT_EQ(vw[2], 5.0);


	EXPECT_EQ(vw.size(), 3);


}



TEST(TestViews, applyFilterSmallLambda)
{


	std::vector<double>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };

	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x >= 2.0; };

	VecVW vw = ApplyFilterW(condTest, test);


	EXPECT_EQ(vw[0], 2.0);
	EXPECT_EQ(vw[1], 3.0);
	EXPECT_EQ(vw[2], 5.0);


	EXPECT_EQ(vw.size(), 3);


}



TEST(TestViews, writeView)
{
	// write values back to source position

	std::vector<double>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };

	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x >= 2.0; };

	VecVW vw = ApplyFilterW(condTest, test);


	EXPECT_EQ(vw[0], 2.0);
	EXPECT_EQ(vw[1], 3.0);
	EXPECT_EQ(vw[2], 5.0);


	EXPECT_EQ(vw.size(), 3);

	VecXX myVec(test); //keep target sdame size as input view
	myVec *= 0.0;

	vw.write(myVec);


	EXPECT_EQ(myVec[1], 2.0);
	EXPECT_EQ(myVec[3], 3.0);
	EXPECT_EQ(myVec[4], 5.0);

}


TEST(TestViews, writeThroView)
{
	// write values back to source position

	std::vector<double>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };

	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x >= 2.0; };

	auto condTest2 = [](auto x) {return x >= 5.0; };

	VecVW vw = ApplyFilterW(condTest, test);


	EXPECT_EQ(vw[0], 2.0);
	EXPECT_EQ(vw[1], 3.0);
	EXPECT_EQ(vw[2], 5.0);


	EXPECT_EQ(vw.size(), 3);


	VecVW vw2 = ApplyFilterZ(condTest2, vw);


	VecXX myVec(test); //keep target sdame size as input view
	myVec *= 0.0;

	vw2.write(myVec);


	//EXPECT_EQ(myVec[1], 2.0);
	//EXPECT_EQ(myVec[3], 3.0);
	EXPECT_EQ(myVec[4], 5.0);

}




TEST(TestViews, writeView2)
{
	// write values back to source position

	std::vector<double>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };

	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x >= 2.0; };

	VecVW vw = ApplyFilterX(condTest, test);

	auto squareIt = [](auto x) {return x * x; };


	ApplyUnitaryOperationX(squareIt, vw);
	for (auto x : vw)
	{
		std::cout << x << std::endl;
	}

	VecXX copyOfTest(test);
	for (auto x : copyOfTest)
	{
		std::cout << x << std::endl;
	}

	vw.write(copyOfTest);
	for (auto x : copyOfTest)
	{
		std::cout << x << std::endl;
	}



	//vw = ApplyLambda(vw, squareIt);


}

*/


