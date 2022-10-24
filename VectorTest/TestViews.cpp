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


#include "../Vectorisation/VecX/dr3.h"
#include "testNamespace.h"
#include "dr3TestUtil.h"


#include <numeric>

TEST(TestViews, ApplyUnitaryOperation)
{
	for (int SZ = 3; SZ < 123; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;
		for (auto& x : v) { x -= asNumber(500.0 + i); ++i; }
		VecXX test(v);

		auto mySquareItLambda = [](const auto& x) {return x * x;  };
		auto testSquare = ApplyUnitaryOperation(test, mySquareItLambda);

		for (int i = 0; i < test.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(testSquare[i], test[i] * test[i]);
		}
	}


}




TEST(TestViews, ApplyUnitaryOperationXXX)
{

	for (int SZ = 3; SZ < 123; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;
		for (auto& x : v) { x -= asNumber(500.0 + i);; ++i; }
		VecXX test(v);

		auto mySquareItLambda = [](const auto& x) {return x * x;  };
		auto testSquare = ApplyUnitaryOperationV(test, mySquareItLambda);

		for (int i = 0; i < test.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(testSquare[i], test[i] * test[i]);
		}
	}

}




TEST(TestViews, ApplyBinaryOperation)
{

	for (int SZ = 3; SZ < 123; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;
		for (auto& x : v) { x -= asNumber(500.0 + i);; ++i; }
		VecXX test(v);
		auto test2 = test * 3.0;
		auto myDifItLambda = [](const auto& x, const auto& y) {return x-y;  };
		auto testDiff= ApplyBinaryOperation1(test,test2, myDifItLambda);

		EXPECT_NUMERIC_EQ(testDiff[0], test[0] - test2[0]);
		EXPECT_NUMERIC_EQ(testDiff[SZ-1], test[SZ - 1] - test2[SZ - 1]);
		EXPECT_NUMERIC_EQ(testDiff[SZ /2], test[SZ/2] - test2[SZ/2]);


		for (int i = 1; i < test.size()-1; ++i)
		{
			EXPECT_NUMERIC_EQ(testDiff[i], test[i] - test2[i]);
		}
	}


}


TEST(TestViews, applyFilter1000)
{

	int SZ = 1000;

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i);; ++i; }
	VecXX test(v);
	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	initVals += 3000.0;
	auto condTestVals = initVals <= 200.0;

	VecVW vw = ApplyFilterB(condTestVals, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals = !condTestVals;

	VecVW vwNegate = ApplyFilterB(NotCondTestVals, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE( (vwNegate.size() + vw.size() )== initVals.size());



}



TEST(TestViews, applyFilter3) 
{

	int SZ = 3;
	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i);; ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	initVals += 3000.0;
	auto condTestVals = initVals <= 200.0;
	auto initSz = initVals.size();

	VecVW vw = ApplyFilterB(condTestVals, initVals);

	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals = initVals > 200.0;//= !condTestVals;
	auto vwNegate = ApplyFilterB(NotCondTestVals, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE((vwNegate.last() + vw.last())== initSz);
}


TEST(TestViews, applyFilterRange4)
{

	int SZ = 4;

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i);; ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	initVals += 3000.0;
	auto condTestVals = initVals <= 200.0;

	VecVW vw = ApplyFilterB(condTestVals, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals = !condTestVals;

	VecVW vwNegate = ApplyFilterB(NotCondTestVals, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE((vwNegate.size() + vw.size()) == initVals.size());

}


TEST(TestViews, applyFilterRange4LambdaX)
{

	int SZ = 4;

	std::vector<Numeric>  v(SZ, asNumber(6.66)) ;
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i);; ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	initVals += 3000.0;
	//auto condTestVals = initVals <= 200.0;

	auto condTestValsLambda = [](auto x) { return x <= VecXX::INS(200.0); };

	VecVW vw = ApplyFilter(condTestValsLambda, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestValsLambda =  [](auto x) { return x > VecXX::INS(200.0); };

	VecVW vwNegate = ApplyFilter(NotCondTestValsLambda, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE((vwNegate.size() + vw.size()) == initVals.size());

}


TEST(TestViews, applyFilterRange4LambdaZ)
{

	int SZ = 4;

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i);; ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	initVals += 3000.0;

	auto condTestValsLambda = [](auto x) { return x <= VecXX::INS(200.0); };

	VecVW vw = ApplyFilter(condTestValsLambda, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestValsLambda = [](auto x) { return x > VecXX::INS(200.0); };

	VecVW vwNegate = ApplyFilter(NotCondTestValsLambda, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE((vwNegate.size() + vw.size()) == initVals.size());

}


TEST(TestViews, applyFilterRange4LambdaW)
{

	int SZ = 4;

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i);; ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	initVals += 3000.0;

	auto condTestValsLambda = [](auto x) { return x <= VecXX::INS(200.0); };

	VecVW vw = ApplyFilter(condTestValsLambda, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestValsLambda = [](auto x) { return x > VecXX::INS(200.0); };

	VecVW vwNegate = ApplyFilter(NotCondTestValsLambda, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE((vwNegate.size() + vw.size()) == initVals.size());

}






void testFilterVec(int SZ)
{

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i); ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	initVals += 3000.0;
	auto condTestVals = initVals <= 200.0;

	VecVW vw = ApplyFilterB(condTestVals, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals = !condTestVals;

	VecVW vwNegate = ApplyFilterB(NotCondTestVals, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE((vwNegate.last() + vw.last())== initVals.size());
}

TEST(TestViews, applyFilterRangeAll)
{

	for (int SZ = 3; SZ < 153; ++SZ)
	{
		testFilterVec(SZ);
	}

}



void testFilterVecX(int SZ)
{

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i); ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	VecVW vw = ApplyFilter(condTest, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals =  [](auto x) {return x > 200.0; };

	VecVW vwNegate = ApplyFilter(NotCondTestVals, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE((vwNegate.last() + vw.last()) == initVals.size());
}



TEST(TestViews, applyFilterRangeAllX)
{

	for (int SZ = 3; SZ < 153; ++SZ)
	{
		testFilterVecX(SZ);
	}

}


void testFilterVecZ(int SZ)
{

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i); ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };
	VecVW vw = ApplyFilter(condTest, initVals);


	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals = [](auto x) {return x > 200.0; };

	VecVW vwNegate = ApplyFilter(NotCondTestVals, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE((vwNegate.last() + vw.last()) == initVals.size());
}



TEST(TestViews, applyFilterRangeAllZ)
{

	for (int SZ = 3; SZ < 153; ++SZ)
	{
		testFilterVecZ(SZ);
	}

}



void testFilterVecW(int SZ)
{

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i); ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };
	VecVW vw = ApplyFilter(condTest, initVals);

	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals = [](auto x) {return x > 200.0; };

	VecVW vwNegate = ApplyFilter(NotCondTestVals, initVals);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE((vwNegate.last() + vw.last()) == initVals.size());
}



TEST(TestViews, applyFilterRangeAllW)
{

	for (int SZ = 3; SZ < 153; ++SZ)
	{
		testFilterVecW(SZ);
	}

}




void testFilterViewZ(int SZ)
{

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i); ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };

	auto evenTest = [](auto x) {return (floor((x - asNumber(0.00000001))) / asNumber(2.0)) == ( floor( (x- asNumber(0.00000001)) / asNumber(2.0))); };

	VecVW evenView = ApplyFilter(evenTest, initVals);
	VecVW vw = ApplyFilter(condTest, evenView);

	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}

	auto NotCondTestVals = [](auto x) {return x > 200.0; };

	VecVW vwNegate = ApplyFilter(NotCondTestVals, evenView);

	for (auto x : vwNegate)
	{
		EXPECT_TRUE(x > 200.0);
	}

	EXPECT_TRUE((vwNegate.last() + vw.last()) == evenView.size());
}



TEST(TestViews, applyFilterRangeAllViewZ)
{

	for (int SZ = 3; SZ < 153; ++SZ)
	{
		testFilterViewZ(SZ);
	}

}


void testCountedFilterZ(int SZ, int count)
{

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i); ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };
	VecVW vw = ApplyCountedFilter(condTest, test, count);

	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}
	EXPECT_TRUE(vw.size() <= count);
}


void testCountedFilterZAllTrue(int SZ, int count)
{

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i); ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x > -10000.0; };

	VecVW vw = ApplyCountedFilter(condTest, test, count);

	for (auto x : vw)
	{
		EXPECT_TRUE(x > -10000.0);
	}
	EXPECT_TRUE(vw.size() == count);
}



TEST(TestViews, applyCountedFilterZRangeAll)
{
	for (int SZ = 3; SZ < 153; ++SZ)
	{
		for (int count = 0; count <SZ;count++ )
			testCountedFilterZ(SZ, count);
	}
}


TEST(TestViews, applyCountedFilterZRangeAllTrue)
{
	for (int SZ = 3; SZ < 153; ++SZ)
	{
		for (int count = 0; count < SZ; count++)
			testCountedFilterZ(SZ, count);
	}
}


void testCountedFilterViewZ(int SZ, int count)
{

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i); ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x <= 200.0; };
	auto evenTest = [](auto x) {return (floor((x - asNumber( 0.00000001))) / asNumber(2.0)) == (floor((x - asNumber(0.00000001)) / asNumber(2.0))); };

	VecVW evenView = ApplyFilter(evenTest, initVals);
	VecVW vw = ApplyCountedFilter(condTest, evenView, count);

	for (auto x : vw)
	{
		EXPECT_TRUE(x <= 200.0);
	}
	EXPECT_TRUE(vw.size() <= count);
}




TEST(TestViews, applyCountedFilterViewZRangeAll)
{
	for (int SZ = 3; SZ < 153; ++SZ)
	{
		for (int count = 0; count < SZ; count++)
			testCountedFilterViewZ(SZ, count);
	}
}




void testCountedFilterViewZAllTrue(int SZ, int count)
{

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i); ++i; }
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x > VecXX::INS(-1200.0); };
	auto trueFilter = [](auto x) {return x==x; };

	VecVW trueView = ApplyFilter(trueFilter, initVals);
	VecVW vw = ApplyCountedFilter(condTest, trueView, count);

	for (auto x : vw)
	{
		EXPECT_TRUE(x > -1200.0);
	}
	EXPECT_TRUE(vw.size() == count);
}



TEST(TestViews, applyCountedFilterViewZRangeAllTrue)
{
	for (int SZ = 3; SZ < 153; ++SZ)
	{
		for (int count = 0; count < SZ; count++)
			testCountedFilterViewZAllTrue(SZ, count);
	}
}



TEST(TestViews, ApplyOperationAndFilter)
{

	int SZ = 21;
	std::vector<Numeric>  v(SZ, asNumber(0.0));
	int i = 0;
	for (auto& x : v) { x = asNumber(i); ++i; }
	VecXX data(v);

	auto SQR = [](auto x) { return x * x; };
	auto aboveFive = [](auto x) { return x > VecXX::scalar(5.0); };
	auto lessThanFifteen = [](auto x) { return x < VecXX::scalar(15.0); };

	using namespace JOIN;
	auto betweenFiveAndFiffteen = aboveFive && lessThanFifteen;
	auto res = ApplyOperationAndFilter(SQR, betweenFiveAndFiffteen, data);

	auto& vec = res.first;
	
	for (int ii = 0; ii < SZ-1; ++ii)
	{
		double i = ii;
		size_t pos = size_t(i);
		auto val = vec[pos];
		auto expected = i * i;
		EXPECT_DOUBLE_EQ(val, expected);
	}

	auto& view = res.second;
	for (int jj = 0; jj < 9; ++jj)
	{
		Numeric j = asNumber(jj);
		EXPECT_DOUBLE_EQ(view[jj] , 6.0+j);
	}

}


TEST(TestViews, ApplyOperationAndDoubleFilter)
{

	int SZ = 21;
	std::vector<Numeric>  v(SZ, asNumber(0.0));
	int i = 0;
	for (auto& x : v) { x = asNumber(i); ++i; }
	VecXX data(v);

	auto SQR = [](auto x) { return x * x; };

	auto aboveFive = [](auto x) { return x > VecXX::scalar(5.0); };

	auto lessThanFifteen = [](auto x) { return x < VecXX::scalar(15.0); };

	using namespace JOIN;
	auto betweenFiveAndFiffteen = aboveFive && lessThanFifteen;

	auto res = ApplyOperationAndFilter(SQR, aboveFive, lessThanFifteen, data);

	auto& vec = std::get<0>(res);

	for (int ii = 0; ii < SZ - 1; ++ii)
	{
		double i = ii;
		size_t pos = size_t(i);
		auto val = vec[pos];
		auto expected = i * i;
		EXPECT_DOUBLE_EQ(val, expected);
	}


	auto& view1 = std::get<1>(res);  //above 5
	for (int jj = 0; jj < 14; ++jj)
	{
		double j = jj;
		EXPECT_DOUBLE_EQ(view1[jj], 6.0 + j);
	}

	auto& view2 = std::get<2>(res); 
	for (int jj = 0; jj < 15; ++jj) //less than 15
	{
		double j = jj;
		EXPECT_DOUBLE_EQ(view2[jj],  j);
	}

}



TEST(TestViews, ApplyOperationAndWrite)
{

	int SZ = 21;
	std::vector<Numeric>  v(SZ, asNumber(0.0));
	int i = 0;
	for (auto& x : v) { x = asNumber(i); ++i; }
	VecXX data(v);

	auto SQR = [](auto x) { return x * x; };
	auto aboveFive = [](auto x) { return x > VecXX::scalar(5.0); };
	auto lessThanFifteen = [](auto x) { return x < VecXX::scalar(15.0); };

	using namespace JOIN;
	
	auto betweenFiveAndFiffteen = aboveFive && lessThanFifteen;
	auto view_WithElementsAboveFive = ApplyFilter(aboveFive, data);

	std::vector<Numeric> valss = view_WithElementsAboveFive;
	std::vector<int> index = view_WithElementsAboveFive.getIndex();
	auto copyOfData = data;

    ApplyUnitaryOperationWrite(SQR, view_WithElementsAboveFive, copyOfData);
	//squares allvalues above five and writes back to copy of source data;

	int j = 0;
	for (; j < 6; ++j) // 5 or less
	{
		Numeric jj = asNumber(j);
		EXPECT_DOUBLE_EQ(copyOfData[j], jj);
	}

	for (; j < 21; ++j) // above 5 is squared values
	{
		Numeric jj = asNumber(j);
		EXPECT_DOUBLE_EQ(copyOfData[j], jj*jj);
	}

}


TEST(TestViews, JoiningLambdas)
{

	int SZ = 21;
	std::vector<Numeric>  v(SZ, asNumber(0.0));
	int i = 0;
	for (auto& x : v) { x = asNumber(i); ++i; }
	VecXX data(v);

	auto SQR = [](auto x) { return x * x; };
	auto SQR_ROOT = [](auto x) { return sqrt(x); };

	auto addFive = [](auto x) {return x + VecXX::scalar(5.0); };
	auto timesTwo = [](auto x) {return x * VecXX::scalar(2.0); };

	auto copyOfData = data;
	auto res_sqr = ApplyUnitaryOperationV(data,SQR);
	using namespace JOIN;
	auto quartic = SQR | SQR;

	auto res_quartic = ApplyUnitaryOperationV(data,quartic);
	auto roundTrip = SQR_ROOT | SQR;
	auto res_quartic_two = ApplyUnitaryOperationV(data,roundTrip);
	auto doIt = addFive | timesTwo;  //2x+5  
	auto res_qtest = ApplyUnitaryOperationV(data,doIt);

	for (int i = 0; i < data.size(); ++i)
	{
		Numeric j = asNumber(i);
		EXPECT_DOUBLE_EQ(res_qtest[i], (5. + j) * 2.0);
	}

}











TEST(TestViews, binaryFilterZVec)
{
	std::vector<Numeric>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };
	VecXX test(v);

	//auto initVals = test;
	auto condTest = [](auto x) {return x >= 2.0; };

	auto tupple_vw = ApplyBinaryFilter(condTest, test);

	auto vw = std::get<0>(tupple_vw);


	EXPECT_NUMERIC_EQ(vw[0], asNumber(2.0));
	EXPECT_NUMERIC_EQ(vw[1], asNumber(3.0));
	EXPECT_NUMERIC_EQ(vw[2], asNumber(5.0) );
	EXPECT_NUMERIC_EQ(vw.size(), 3);


	auto vw_other = std::get<1>(tupple_vw);
	EXPECT_NUMERIC_EQ(vw_other[0], asNumber(1.0));
	EXPECT_NUMERIC_EQ(vw_other[1], asNumber(1.0));
	EXPECT_NUMERIC_EQ(vw_other[2], asNumber(1.0));
	EXPECT_NUMERIC_EQ(vw_other.size(), 3);

	auto addFive = [](auto rhs) { return rhs + VecXX::scalar(5.0); };
	ApplyUnitaryOperation(vw_other,addFive);

	//vw_other.writeView(test);
	vw_other.write(test);
	EXPECT_NUMERIC_EQ(test[0], asNumber(6.0));
	EXPECT_NUMERIC_EQ(test[2], asNumber(6.0));
	EXPECT_NUMERIC_EQ(test[5], asNumber(6.0));

	EXPECT_NUMERIC_EQ(test[1], asNumber(2.0));
	EXPECT_NUMERIC_EQ(test[3], asNumber(3.0));
	EXPECT_NUMERIC_EQ(test[4], asNumber(5.0));

}


TEST(TestViews, binaryFilterZView)
{

	std::vector<Numeric>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };

	VecXX test(v);

	auto trueTest = [](auto x) {return x==x; };

	auto initialVW = ApplyFilter(trueTest, test);

	auto condTest = [](auto x) {return x >= VecXX::INS(2.0); };

	auto tupple_vw = ApplyBinaryFilter(condTest, initialVW);

	auto vw = std::get<0>(tupple_vw);


	EXPECT_NUMERIC_EQ(vw[0], asNumber(2.0));
	EXPECT_NUMERIC_EQ(vw[1], asNumber(3.0));
	EXPECT_NUMERIC_EQ(vw[2], asNumber(5.0));
	EXPECT_NUMERIC_EQ(vw.size(), 3);


	auto vw_other = std::get<1>(tupple_vw);
	EXPECT_NUMERIC_EQ(vw_other[0], asNumber(1.0));
	EXPECT_NUMERIC_EQ(vw_other[1], asNumber(1.0));
	EXPECT_NUMERIC_EQ(vw_other[2], asNumber(1.0));
	EXPECT_NUMERIC_EQ(vw_other.size(), 3);

	auto addFive = [](auto rhs) { return rhs + VecXX::scalar(5.0); };
	ApplyUnitaryOperation(vw_other,addFive);

	vw_other.writeView(test);
	EXPECT_NUMERIC_EQ(test[0], asNumber(6.0));
	EXPECT_NUMERIC_EQ(test[2], asNumber(6.0));
	EXPECT_NUMERIC_EQ(test[5], asNumber(6.0));

	EXPECT_NUMERIC_EQ(test[1], asNumber(2.0));
	EXPECT_NUMERIC_EQ(test[3], asNumber(3.0));
	EXPECT_NUMERIC_EQ(test[4], asNumber(5.0));

}





TEST(TestViews, applyFilterSmall)
{


	std::vector<Numeric>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };

	VecXX test(v);

	auto initVals = test;

	//boolean vec
	auto condTestVals = test >= 2.0;

	VecVW vw = ApplyFilterB(condTestVals, test);


	EXPECT_NUMERIC_EQ(vw[0], asNumber(2.0));
	EXPECT_NUMERIC_EQ(vw[1], asNumber(3.0));
	EXPECT_NUMERIC_EQ(vw[2], asNumber(5.0));
	EXPECT_NUMERIC_EQ(vw.size(), 3);

}



TEST(TestViews, applyFilterSmallLambda)
{


	std::vector<Numeric>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };

	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x >= 2.0; };

	VecVW vw = ApplyFilter(condTest, test );

	EXPECT_NUMERIC_EQ(vw[0], asNumber(2.0));
	EXPECT_NUMERIC_EQ(vw[1], asNumber(3.0));
	EXPECT_NUMERIC_EQ(vw[2], asNumber(5.0));
	EXPECT_NUMERIC_EQ(vw.size(), 3);
}



TEST(TestViews, writeView)
{
	// write values back to source position
	std::vector<Numeric>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };
	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x >= 2.0; };

	VecVW vw = ApplyFilter(condTest, test);

	EXPECT_NUMERIC_EQ(vw[0], asNumber(2.0));
	EXPECT_NUMERIC_EQ(vw[1], asNumber(3.0));
	EXPECT_NUMERIC_EQ(vw[2], asNumber(5.0));
	EXPECT_NUMERIC_EQ(vw.size(), 3);

	VecXX myVec(test); //keep target sdame size as input view
	myVec *= 0.0;
	vw.write(myVec);

	EXPECT_NUMERIC_EQ(vw[0], asNumber(2.0));
	EXPECT_NUMERIC_EQ(vw[1], asNumber(3.0));
	EXPECT_NUMERIC_EQ(vw[2], asNumber(5.0));

}


TEST(TestViews, writeThroView)
{
	// write values back to source position
	std::vector<Numeric>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };

	VecXX test(v);

	auto initVals = test;
	auto condTest = [](auto x) {return x >= 2.0; };
	auto condTest2 = [](auto x) {return x >= 5.0; };

	VecVW vw = ApplyFilter(condTest, test);

	EXPECT_NUMERIC_EQ(vw[0], asNumber(2.0));
	EXPECT_NUMERIC_EQ(vw[1], asNumber(3.0));
	EXPECT_NUMERIC_EQ(vw[2], asNumber(5.0));
	EXPECT_NUMERIC_EQ(vw.size(), 3);


	VecVW vw2 = ApplyFilter(condTest2, vw);


	VecXX myVec(test); //keep target sdame size as input view
	myVec *= 0.0;

	vw2.write(myVec);
	EXPECT_NUMERIC_EQ(myVec[4], asNumber(5.0));

}



TEST(TestViews, writeView2)
{
	// write values back to source position
	std::vector<Numeric>  v = { 1.0,2.0,1.0,3.0,5.0,1.0 };

	VecXX test(v);
	auto initVals = test;
	auto condTest = [](auto x) {return x >= 2.0; };

	VecVW vw = ApplyFilter(condTest, test);
	auto squareIt = [](auto x) {return x * x; };

	 ApplyUnitaryOperation(vw,squareIt );
	 for (const auto& x : vw)
	 {
		 std::cout << x << std::endl;
	 }

	 VecXX copyOfTest(test);
	 for (const auto& x : copyOfTest)
	 {
		 std::cout << x << std::endl;
	 }

	 vw.write(copyOfTest);
	 for (const auto& x : copyOfTest)
	 {
		 std::cout << x << std::endl;
	 }

}


