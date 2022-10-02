#include "pch.h"

//#include "../Vectorisation/VecX/instruction_traits.h"


#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/vec_bool_d.h"
#include "../Vectorisation/VecX/vec_double.h"
#include  "../Vectorisation/VecX/alloc_policy.h"
#include "../Vectorisation/VecX/vec_d.h"
#include "../Vectorisation/VecX/vec_bool.h"
#include "../Vectorisation/VecX/vec_view.h"

#include "../Vectorisation/VecX/target_name_space.h"




#include <numeric>

//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;




void testSelectFilterB(int SZ, int Pos)
{
	//int SZ = 21;
	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX data(v);

	auto boolVecAboveTenAndHalf = data < VecXX::scalar(Pos);

	auto filterView = ApplyFilterB(boolVecAboveTenAndHalf, data);

	EXPECT_EQ(filterView.size(), Pos);

	for (int i = 0; i < filterView.size(); ++i)
	{
		EXPECT_DOUBLE_EQ(filterView[i], i );
	}

	auto trueCond = [](auto x) { return 1.0 == 1.0; };

	//repeat with data in a view
	{	auto viewOfData = ApplyFilter(trueCond, data);
		auto filterView = ApplyFilterB(boolVecAboveTenAndHalf, viewOfData);

		EXPECT_EQ(filterView.size(), Pos);

		for (int i = 0; i < filterView.size(); ++i)
		{
			EXPECT_DOUBLE_EQ(filterView[i], i);
		}
	}
}

TEST(TestFilterSelect, ApplyFilterB)
{

	for (int SZ = 3; SZ < 131; SZ++)
	{
		for (int Pos = 0; Pos <= SZ; Pos++)
		{
			testSelectFilterB(SZ, Pos);
		}
	}

}



void testSelectFilter(int SZ, int Pos)
{
	//int SZ = 21;
	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX data(v);

	auto abovePos = [=](auto x) { return  x < VecXX::scalar(Pos); };

	auto filterView = ApplyFilter(abovePos, data);

	EXPECT_EQ(filterView.size(), Pos);

	for (int i = 0; i < filterView.size(); ++i)
	{
		EXPECT_DOUBLE_EQ(filterView[i], i);
	}

	auto trueCond = [](auto x) { return 1.0 == 1.0; };

	//repeat with data in a view
	{	
		auto viewOfData = ApplyFilter(trueCond, data);
		auto filterView = ApplyFilter(abovePos, viewOfData);

		EXPECT_EQ(filterView.size(), Pos);

		for (int i = 0; i < filterView.size(); ++i)
		{
			EXPECT_DOUBLE_EQ(filterView[i], i);
		}
	}
}


TEST(TestFilterSelect, ApplyFilter)
{

	for (int SZ = 3; SZ < 131; SZ++)
	{
		for (int Pos = 0; Pos <= SZ; Pos++)
		{
			testSelectFilter(SZ, Pos);
		}
	}

}



void testCountedFilter(int SZ, int Pos, int count)
{
	//int SZ = 21;
	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX data(v);

	auto abovePos = [=](auto x) { return  x < VecXX::scalar(Pos); };

	auto filterView = ApplyCountedFilter(abovePos, data, count);

	EXPECT_EQ(filterView.size(), std::min(count,Pos) );

	for (int i = 0; i < filterView.size(); ++i)
	{
		EXPECT_DOUBLE_EQ(filterView[i], i);
	}

	auto trueCond = [](auto x) { return 1.0 == 1.0; };

	//repeat with data in a view
	{
		auto viewOfData = ApplyFilter(trueCond, data);
		auto filterView = ApplyCountedFilter(abovePos, viewOfData,count);

		EXPECT_EQ(filterView.size(), std::min(count, Pos));

		for (int i = 0; i < filterView.size(); ++i)
		{
			EXPECT_DOUBLE_EQ(filterView[i], i);
		}
	}
}


TEST(TestFilterSelect, ApplyCountedFilter)
{

	for (int SZ = 3; SZ < 71; SZ++)
	{
		for (int Pos = 0; Pos <= SZ; Pos++)
		{
			for (int count  =0; count < (SZ+5); count++)
			testCountedFilter(SZ, Pos, count);
		}
	}

}



void testBinaryFilter(int SZ, int Pos)
{
	//int SZ = 21;
	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX data(v);

	auto abovePos = [=](auto x) { return  x < VecXX::scalar(Pos); };
	auto filterView = ApplyBinaryFilter(abovePos, data);
	auto trueView = std::get<0>(filterView);
	EXPECT_EQ(trueView.size(),  Pos);

	for (int i = 0; i < trueView.size(); ++i)
	{
		EXPECT_DOUBLE_EQ(trueView[i], i);
	}


	auto falseView = std::get<1>(filterView);
	EXPECT_EQ(falseView.size(), (SZ-Pos));
	for (int i = 0; i < falseView.size(); ++i)
	{
		EXPECT_DOUBLE_EQ(falseView[i], i+Pos);
	}


	auto trueCond = [](auto x) { return 1.0 == 1.0; };

	//repeat with data in a view
	{
		auto viewOfData = ApplyFilter(trueCond, data);
		auto filterView = ApplyBinaryFilter(abovePos, viewOfData);
		auto trueView = std::get<0>(filterView);
		EXPECT_EQ(trueView.size(), Pos);

		for (int i = 0; i < trueView.size(); ++i)
		{
			EXPECT_DOUBLE_EQ(trueView[i], i);
		}

		auto falseView = std::get<1>(filterView);
		EXPECT_EQ(falseView.size(), (SZ - Pos));
		for (int i = 0; i < falseView.size(); ++i)
		{
			EXPECT_DOUBLE_EQ(falseView[i], i + Pos);
		}
	}

}


TEST(TestFilterSelect, ApplyBinaryFilter)
{

	for (int SZ = 3; SZ < 131; SZ++)
	{
		for (int Pos = 0; Pos <= SZ; Pos++)
		{
			testBinaryFilter(SZ, Pos);
		}
	}
}






void testUnitaryOpVecViewRef(int SZ)
{

	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX data(v);

	auto SQR = [](auto x) { return  x * x; };

	auto trueCond = [](auto x) { return 1.0 == 1.0; };

	auto viewOfData = ApplyFilter(trueCond, data);

    ApplyUnitaryOperation(SQR, viewOfData);

	auto resultView = viewOfData; 

	EXPECT_EQ(resultView.size(), SZ);

	for (int i = 0; i < resultView.size(); ++i)
	{
		double ii = i;
		EXPECT_DOUBLE_EQ(resultView[i], i*i);
	}


}


TEST(TestFilterSelect, ApplyUnitaryOperationToView)
{
	for (int SZ = 1; SZ < 131; SZ++)
	{
		testUnitaryOpVecViewRef(SZ);
	}
}



void testUnitaryOpVec(int SZ)
{

	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX data(v);

	auto SQR = [](auto x) { return  x * x; };


	auto resultView = ApplyUnitaryOperationV(SQR, data);

	EXPECT_EQ(resultView.size(), SZ);

	for (int i = 0; i < resultView.size(); ++i)
	{
		double ii = i;
		EXPECT_DOUBLE_EQ(resultView[i], i * i);
	}


}


TEST(TestFilterSelect, ApplyUnitaryOperationToVec)
{
	for (int SZ = 1; SZ < 131; SZ++)
	{
		testUnitaryOpVec(SZ);
	}
}



void testBinaryOpVecViewRef(int SZ)
{

	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX data_Lhs(v);

	std::vector<double>  r(SZ, 0.0);
	i = 0;
	for (auto& x : r) { x = SZ-i; ++i; }
	VecXX data_Rhs(r);


	auto ADD = [](auto x,auto y) { return  x + y; };

	auto trueCond = [](auto x) { return 1.0 == 1.0; };

	auto viewOf_Lhs = ApplyFilter(trueCond, data_Lhs);

	auto viewOf_Rhs = ApplyFilter(trueCond, data_Rhs);

	auto resultView = 	ApplyBinaryOperation(ADD, viewOf_Lhs, viewOf_Rhs);



	EXPECT_EQ(resultView.size(), SZ);

	for (int i = 0; i < resultView.size(); ++i)
	{
		double ii = i;
		EXPECT_DOUBLE_EQ(resultView[i], SZ );
	}


}



TEST(TestFilterSelect, ApplyBinaryOperationToView)
{
	for (int SZ = 1; SZ < 131; SZ++)
	{
		testBinaryOpVecViewRef(SZ);
	}
}



// TO DO TEST COMPLEC COMPOUND OPERATIONS


// TO DO TEST JOINS Logical and  functional

TEST(TestFilterSelect, ApplyJoinedFilterToVec)
{
	const int SZ = 40;

	const double NUMBER = SZ * 0.5;
	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX data(v);


	using namespace JOIN;
	auto lessThanThirty = [](const auto& X) { return X < static_cast<decltype(X)>(30); };
	auto isEven = [](auto X) { return (floor(X) / 2.0 - floor(floor(X) / 2.0)) < 0.0001; };

	auto isEvenLessThanThirty = lessThanThirty && isEven;

	auto res = ApplyFilter(isEvenLessThanThirty, data);

	for (auto x : res)
	{
		EXPECT_TRUE(lessThanThirty(x));
		EXPECT_TRUE(isEven(x));
	}

	EXPECT_TRUE(res.size() == 15);

	auto isThirtyOrMore = !lessThanThirty;

	auto resNegate = ApplyFilter(isThirtyOrMore, data);

	for (auto x : resNegate)
	{
		EXPECT_TRUE( x >= 30.);
	}

	EXPECT_TRUE(resNegate.size() == 10);




	auto isEvenOrLessThanThirty = lessThanThirty || isEven;

	auto resOr = ApplyFilter(isEvenOrLessThanThirty, data);

	for (auto x : resOr)
	{
		EXPECT_TRUE(lessThanThirty(x) || isEven(x));
	}

	EXPECT_TRUE(resOr.size() == 35);


}




TEST(TestFilterSelect, ApplyJoinedFilterToView)
{
	const int SZ = 40;

	const double NUMBER = SZ * 0.5;
	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX dataVec(v);

	auto trueCond = [](auto x) { return 1.0 == 1.0; };

	auto data = ApplyFilter(trueCond, dataVec);


	using namespace JOIN;
	auto lessThanThirty = [](const auto& X) { return X < static_cast<decltype(X)>(30); };
	auto isEven = [](auto X) { return (floor(X) / 2.0 - floor(floor(X) / 2.0)) < 0.0001; };

	auto isEvenLessThanThirty = lessThanThirty && isEven;

	auto res = ApplyFilter(isEvenLessThanThirty, data);

	for (auto x : res)
	{
		EXPECT_TRUE(lessThanThirty(x));
		EXPECT_TRUE(isEven(x));
	}

	EXPECT_TRUE(res.size() == 15);

	auto isThirtyOrMore = !lessThanThirty;

	auto resNegate = ApplyFilter(isThirtyOrMore, data);

	for (auto x : resNegate)
	{
		EXPECT_TRUE(x >= 30.);
	}

	EXPECT_TRUE(resNegate.size() == 10);




	auto isEvenOrLessThanThirty = lessThanThirty || isEven;

	auto resOr = ApplyFilter(isEvenOrLessThanThirty, data);

	for (auto x : resOr)
	{
		EXPECT_TRUE(lessThanThirty(x) || isEven(x));
	}

	EXPECT_TRUE(resOr.size() == 35);


}





TEST(TestFilterSelect, ApplyComplexJoinedFilterToVec)
{
	const int SZ = 40;

	const double NUMBER = SZ * 0.5;
	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX data(v);


	using namespace JOIN;
	auto lessThanThirty = [](const auto& X) { return X < static_cast<decltype(X)>(30); };
	auto isEven = [](auto X) { return (floor(X) / 2.0 - floor(floor(X) / 2.0)) < 0.0001; };

	auto isEvenLessThanThirty = lessThanThirty && isEven;


	auto isGreaterThanTen = [](auto X) { return X > static_cast<decltype(X)>(10); };

	auto isGreaterThanTwenty = [](auto X) { return X > static_cast<decltype(X)>(20); };

	auto isGreaterThanThirty = [](auto X) { return X > static_cast<decltype(X)>(30); };

	auto isLessThanthree = [](auto X) { return X < static_cast<decltype(X)>(3); };


   auto complex = ((!isGreaterThanThirty) && isGreaterThanTwenty && isEven) || isLessThanthree; 
   //even numbers  22 -30    and 0,1,2,3

   auto res = ApplyFilter(complex, data);

   for (auto x : res)
   {
	   if (x < 3)
	   {
		   EXPECT_TRUE(true);
	   }
	   else
	   {
		   	EXPECT_TRUE(x <=30.0);
			EXPECT_TRUE(isEven(x));
			EXPECT_TRUE(x > 20.0);
	   }
   }

   EXPECT_TRUE(res.size() == 8);



   auto trueCond = [](auto x) { return 1.0 == 1.0; };

   auto dataView = ApplyFilter(trueCond, data);


   auto resView = ApplyFilter(complex, dataView);

   for (auto x : resView)
   {
	   if (x < 3)
	   {
		   EXPECT_TRUE(true);
	   }
	   else
	   {
		   EXPECT_TRUE(x <= 30.0);
		   EXPECT_TRUE(isEven(x));
		   EXPECT_TRUE(x > 20.0);
	   }
   }

   EXPECT_TRUE(resView.size() == 8);

}




TEST(TestFilterSelect, ApplyComplexJoinedLambdasToVec)
{
	const int SZ = 40;

	const double NUMBER = SZ * 0.5;
	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX data(v);

	auto SQR = [](auto x) { return (x * x); };

	auto AddTen = [](auto x) { return (x + VecXX::INS(10.0)); };

	auto negateIfOverTen = [](auto x) { return iff(x > VecXX::INS(10.), -x, x); };


	using namespace JOIN;
	
	//CAT joins the lambdas and evaluates  left to right in the sequence feeding 
	//the result of the operation into the input of the next lambda
	auto complex2 = AddTen | negateIfOverTen | SQR;

	auto res = ApplyUnitaryOperationV(complex2, data);



	i = 0;
	for (auto resItem : res)
	{
		auto x = data[i];
		++i;

		x += 10;
		x =(x > 10) ? -x : x;
		auto calcResult = x* x;

		EXPECT_DOUBLE_EQ(resItem, calcResult);
	}

	//DO THIS WITH A VIEW

}




TEST(TestFilterSelect, PipeComplexJoinedLambdasToVec)
{
	//using namespace JOIN;

	using namespace PIPE;

	const int SZ = 40;

	const double NUMBER = SZ * 0.5;
	std::vector<double>  v(SZ, 0.0);
	int i = 0;
	for (auto& x : v) { x = i; ++i; }
	VecXX data(v);

	auto SQR = [](auto x) { return (x * x); };

	auto AddTen = [](auto x) { return (x + VecXX::scalar(10.0)); };

	auto negateIfOverTen = [](auto x) { return iff(x > VecXX::scalar(10), -x, x); };

	auto isEven = [](auto X) { return (floor(X) / 2.0 - floor(floor(X) / 2.0)) < 0.0001; };

	using namespace JOIN;

	//CAT joins the lambdas and evaluates  right to left inb the sequence feeding 
	//the result of the operation intyo the input of the next lambda
	//auto complex2 = SQR | negateIfOverTen | AddTen;
	//auto res = ApplyUnitaryOperation(complex2, data);



	auto complex3 =  (data| isEven)> SQR > negateIfOverTen > AddTen ;

	i = 0;
	for (auto resItem : complex3)
	{
		auto x = data[i];
		++i;

		x *= 2.;  //even
		x *= x;   //SQR
		x= (x > 10) ? -x : x;   //negate
		auto calcResult = x +10.0;

		EXPECT_DOUBLE_EQ(resItem, calcResult);


	}

	// Pipes are go left to right 
	// but > are go lrft to right

	auto complexSingleLambda = SQR | negateIfOverTen | AddTen;

	//apply to filtered
	auto complex4 = (data | isEven) > complexSingleLambda;

	i = 0;
	for (auto resItem : complex4)
	{
		auto x = data[i];
		++i;

		x *= 2.;  //even
		x *= x;   //SQR
		x = (x > 10) ? -x : x;   //negate
		auto calcResult = x + 10.0;

		EXPECT_DOUBLE_EQ(resItem, calcResult);
	}

}

