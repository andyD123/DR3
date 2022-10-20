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
#include "dr3TestUtil.h"

#include <numeric>
#include "testNamespace.h"


void testFilterTransform(int SZ )
{

	auto onlyJlambda = [=](auto x) { return (asNumber(j) > (x - asNumber(0.0001)) && (asNumber(j) < x + asNumber(0.00001))); };
	std::vector<Numeric> input(SZ,asNumber( 0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);
	auto trueLambdaS = [&](auto x) { return x; };
	auto falseLambdaS = [&](auto x) { return -x; };


	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return (asNumber(j) > (x - asNumber(0.0001)) && (asNumber(j) < x + asNumber(0.00001))); };
		VecXX res =  filterTransform(onlyJlambda, testVec, trueLambdaS, falseLambdaS);

		for (int k = 0; k < SZ; k++)
		{
			if( k==j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber( k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k],  asNumber(-k));
			}
		}		
	}

}




TEST(TestFilterTransform, testTransformEachPoint)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testFilterTransform(SZ);
	}


	testFilterTransform(34);
	testFilterTransform(65);
	testFilterTransform(63);
	testFilterTransform(64);

}
