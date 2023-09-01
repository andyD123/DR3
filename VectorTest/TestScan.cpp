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
#include "../Vectorisation/VecX/scan.h"
#include "../Vectorisation/VecX/instruction_traits.h"


#include "testNamespace.h"
#include "dr3TestUtil.h"

#include <algorithm>

#include <numeric>





void testScan(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);
	auto add = [](auto x, auto y) {return x + y; };


	for (int j = 0; j < SZ; ++j)
	{

		auto res = scan( testVec, add);

		std::vector<Numeric> dbg = res;

		auto expected = testVec[0];

		EXPECT_NUMERIC_EQ(expected, res[0]);

		for (int k = 1; k < SZ; k++)
		{
			expected += testVec[k] ;
			EXPECT_NUMERIC_EQ(expected, res[k]);	
		}
	}




}



long double getErr(long double)
{
	return std::pow(10, 4 - 16);
}

double getErr(double)
{
	return std::pow(10, 4 - 16);
}

double getErr(float)
{
	return std::pow(10, 4 - 8);
}


void testScan1(int SZ ,double start)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(start));

	
	Numeric err = getErr(Numeric(0.));

	VecXX testVec(input);
	auto add = [](auto x, auto y) {return x + y; };


	for (int j = 0; j < SZ; ++j)
	{

		auto res = scan(testVec, add);

		std::vector<Numeric> dbg = res;
	
		std::vector<Numeric> expected;
		std::inclusive_scan(cbegin(input), cend(input), std::back_inserter( expected));

		EXPECT_NEAR(expected[0], res[0], err);

		for (int k = 1; k < SZ; k++)
		{
			auto relErr = err * std::max(Numeric(1.), std::abs(Numeric(expected[k])));
			EXPECT_NEAR(expected[k], res[k], relErr);
			
		}
	}

}






TEST(TestScan, scanShortVector)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testScan(SZ);
	}

	for (int SZ = 3; SZ < 133; SZ++)
	{
		testScan1(SZ,3.14);
	}

}





void testTransformScan1(int SZ, double start)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(start));


	Numeric err = getErr(Numeric(0.));

	VecXX testVec(input);
	auto SQR = [](auto x) { return x * x; };

	auto sqrVec = transform( [](auto x) {return x * x; }, testVec);
	std::vector< Numeric> sq = sqrVec;
	auto add = [](auto x, auto y) {return x + y; };


	for (int j = 0; j < SZ; ++j)
	{

		auto res = ApplyTransformScan(testVec, add, SQR);

		std::vector<Numeric> dbg = res;

		std::vector<Numeric> expected;
		std::inclusive_scan(cbegin(sq), cend(sq), std::back_inserter(expected));

		EXPECT_NEAR(expected[0], res[0], err);

		for (int k = 1; k < SZ; k++)
		{
			auto relErr = err * std::max(Numeric(1.), std::abs(Numeric(expected[k])));
			EXPECT_NEAR(expected[k], res[k], relErr);

		}
	}

}






TEST(TestTransformScanTransform, transformScanShortVector)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransformScan1(SZ,0);
	}

	for (int SZ = 3; SZ < 133; SZ++)
	{
		testTransformScan1(SZ, 3.14);
	}

}







void testTransformScan2(int SZ, double start)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(start));


	Numeric err = getErr(Numeric(0.));

	VecXX testVec(input);

	VecXX testVec1 = testVec + 1.0;

	auto MULT = [](auto x,auto y) { return x * y; };

	auto multVec = testVec * testVec1;


	std::vector< Numeric> prod = multVec;
	auto add = [](auto x, auto y) {return x + y; };


	for (int j = 0; j < SZ; ++j)
	{

		auto res = ApplyTransformScan(testVec, testVec1, add, MULT);

		std::vector<Numeric> dbg = res;

		std::vector<Numeric> expected;
		std::inclusive_scan(cbegin(prod), cend(prod), std::back_inserter(expected));

		EXPECT_NEAR(expected[0], res[0], err);

		for (int k = 1; k < SZ; k++)
		{
			auto relErr = err * std::max(Numeric(1.), std::abs(Numeric(expected[k])));
			EXPECT_NEAR(expected[k], res[k], relErr);

		}
	}

}






TEST(TestTransformScanTransform, transformScanShortVectorBinary)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransformScan2(SZ, 0);
	}

	for (int SZ = 3; SZ < 133; SZ++)
	{
		testTransformScan2(SZ, 3.14);
	}

}


