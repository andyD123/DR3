#include "pch.h"


#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/vec_bool_d.h"
#include "../Vectorisation/VecX/vec_double.h"
#include  "../Vectorisation/VecX/alloc_policy.h"
#include  "../Vectorisation/VecX/accumulate_transform.h"
#include "../Vectorisation/VecX/target_name_space.h"

#include "../Vectorisation/VecX/dr3.h"
#include "testNamespace.h"
#include "dr3TestUtil.h"

#include <numeric>


auto  getVecBig(int SZ, std::vector<Numeric>& stl)
{
	std::vector<Numeric>  v(SZ, asNumber(1.0/3.0));
	int i = 0;
	
	VecXX test(v);
	stl = v;
	return  test;

}


void evalPrecAccumulate(int startLen, int endLen)
{

	Numeric testEpsilon = 1e-10;
	
	for (int SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<Numeric> v;
		VecXX test = getVecBig(SZ, v);
		using BINNED_ACCUMULATOR = BinsT<VecXX::INS>;
		auto binned_Sum = reduce< BINNED_ACCUMULATOR >(test, BinnedAdd);

		EXPECT_NEAR(double(SZ / 3.0), binned_Sum, testEpsilon);
	}

}

TEST(TestPreciseAccumulator, simpleSummation) 
{
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	//eval over multiple lengths
	evalPrecAccumulate(957, 1043);

	//eval over very small lengths
	evalPrecAccumulate(3, 23);

}


