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


#include "../Vectorisation/VecX/dr3.h"




#include <numeric>

//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;




void testFilterTransform(int SZ )
{

	std::vector<double> input(SZ, 0.0);
	std::iota(begin(input), end(input), 0.0);

	VecXX testVec(input);
	auto trueLambdaS = [&](auto x) { return x; };
	auto falseLambdaS = [&](auto x) { return -x; };


	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return (j > (x - 0.0001) && (j < x + 0.00001)); };
		VecXX res =  filterTransform(onlyJlambda, testVec, trueLambdaS, falseLambdaS);

		for (int k = 0; k < SZ; k++)
		{
			if( k==j)
			{
				EXPECT_DOUBLE_EQ(res[k],  k);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -1.0 * k);
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
