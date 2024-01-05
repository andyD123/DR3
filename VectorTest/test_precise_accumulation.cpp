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


TEST(TestBin, simpleSummation2)
{
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);


	BinsT<VecXX::INS> bin; 



    EXPECT_EQ(bin.veryBigSummV.extract(0), 0.0);
	EXPECT_EQ(bin.bigSummV.extract(0), 0.0);
	EXPECT_EQ(bin.smallSumV.extract(0), 0.0);
	EXPECT_EQ(bin.tinyV.extract(0), 0.0);


	VecXX::INS testVal =1.0e-16;
	bin += testVal;

	EXPECT_EQ(bin.veryBigSummV.extract(0), 0.0);
	EXPECT_EQ(bin.bigSummV.extract(0), 0.0);
	EXPECT_EQ(bin.smallSumV.extract(0), 0.0);
	EXPECT_EQ(bin.tinyV.extract(0), 1.0e-16);

	bin += testVal;

	EXPECT_EQ(bin.veryBigSummV.extract(0), 0.0);
	EXPECT_EQ(bin.bigSummV.extract(0), 0.0);
	EXPECT_EQ(bin.smallSumV.extract(0), 0.0);
	EXPECT_EQ(bin.tinyV.extract(0), 2.0e-16);


	bin += testVal / 2;

	EXPECT_EQ(bin.veryBigSummV.extract(0), 0.0);
	EXPECT_EQ(bin.bigSummV.extract(0), 0.0);
	EXPECT_EQ(bin.smallSumV.extract(0), 0.0);
	EXPECT_EQ(bin.tinyV.extract(0), 2.5e-16);

	//further tests for the other bins

	testVal = 1.0;

	bin += testVal;

	EXPECT_EQ(bin.veryBigSummV.extract(0), 0.0);
	EXPECT_EQ(bin.bigSummV.extract(0), 1.0);
	EXPECT_EQ(bin.smallSumV.extract(0), 0.0);
	EXPECT_EQ(bin.tinyV.extract(0), 2.5e-16);

	bin += testVal;
	EXPECT_EQ(bin.veryBigSummV.extract(0), 0.0);
	EXPECT_EQ(bin.bigSummV.extract(0), 2.0);
	EXPECT_EQ(bin.smallSumV.extract(0), 0.0);
	EXPECT_EQ(bin.tinyV.extract(0), 2.5e-16);

	bin += testVal / 2;
	
	EXPECT_EQ(bin.veryBigSummV.extract(0), 0.0);
	EXPECT_EQ(bin.bigSummV.extract(0), 2.0);
	EXPECT_EQ(bin.smallSumV.extract(0), 0.5);
	EXPECT_EQ(bin.tinyV.extract(0), 2.5e-16);

	



	BinsT<VecXX::INS> bin2;

	auto oneThird = 1.0 / 3.0;


	bin2 += 1.0e-3 * oneThird;


	EXPECT_EQ(bin2.veryBigSummV.extract(0), 0.0);
	EXPECT_EQ(bin2.bigSummV.extract(0), 0.0);
//	EXPECT_EQ(bin2.smallSumV.extract(0), 1.0/3.0 *1.0e-3);
//	EXPECT_EQ(bin2.tinyV.extract(0), 0.0);

	auto sum = bin2.hsum();

	bin2 = bin2 *100000.0;

	sum = bin2.hsum();

	/*

	//eval over multiple lengths
	evalPrecAccumulate(957, 1043);

	//eval over very small lengths
	evalPrecAccumulate(3, 23);
	*/
}

