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




void testTransform_1(int SZ)
{

	auto doubleIt = [](auto x) { return 2.0 * x; };
	VecXX scalar = 3.33;
	VecXX res = transform(doubleIt, scalar);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_DOUBLE_EQ(6.66, val);


	std::vector<double> input(SZ, 0.0);
	std::iota(begin(input), end(input), 0.0);

	VecXX testVec(input);
	

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select ((j > (x - 0.0001) && (j < x + 0.00001)) ,x,-x); };
		VecXX res = transform(onlyJlambda, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res[k], k);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -1.0 * k);
			}
		}
	}


}




TEST(TestDR3, testTransform_1)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransform_1(SZ);
	}


	testTransform_1(34);
	testTransform_1(65);
	testTransform_1(63);
	testTransform_1(64);

}



void testTransform_XX(int SZ)
{

	auto doubleIt = [](auto x) { return 2.0 * x; };
	VecXX scalar = 3.33;
	VecXX res = transformXX(doubleIt, scalar);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_DOUBLE_EQ(6.66, val);



	std::vector<double> input(SZ, 0.0);
	std::iota(begin(input), end(input), 0.0);

	VecXX testVec(input);

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((j > (x - 0.0001) && (j < x + 0.00001)), x, -x); };
		VecXX res = transformXX(onlyJlambda, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res[k], k);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -1.0 * k);
			}
		}
	}

}




TEST(TestDR3, testTransformXX)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransform_XX(SZ);
	}


	testTransform_XX(34);
	testTransform_XX(65);
	testTransform_XX(63);
	testTransform_XX(64);


}




void testTransform_2(int SZ)
{


	auto doubleIt = [](auto x) { return 2.0 * x; };
	VecXX scalar = 3.33;
	VecXX resScalar;
	transform(doubleIt, scalar, resScalar);
	auto val = resScalar.getScalarValue();
	EXPECT_TRUE(resScalar.isScalar());
	EXPECT_DOUBLE_EQ(6.66, val);


	std::vector<double> input(SZ, 0.0);
	std::iota(begin(input), end(input), 0.0);

	VecXX testVec(input);

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((j > (x - 0.0001) && (j < x + 0.00001)), x, -x); };
		VecXX res = testVec;
		transform(onlyJlambda, testVec,res);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res[k], k);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -1.0 * k);
			}
		}
	}

}




TEST(TestDR3, testTransform_2)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransform_2(SZ);
	}


	testTransform_2(34);
	testTransform_2(65);
	testTransform_2(63);
	testTransform_2(64);


}


void testTransform_3(int SZ)
{

	auto doubleIt = [](auto x) { return 2.0 * x; };
	const VecXX scalar = 3.33;
	VecXX res = transform1(doubleIt, scalar);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_DOUBLE_EQ(6.66, val);



	std::vector<double> input(SZ, 0.0);
	std::iota(begin(input), end(input), 0.0);

	const VecXX testVec(input);

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((j > (x - 0.0001) && (j < x + 0.00001)), x, -x); };
		VecXX res = transform1(onlyJlambda, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res[k], k);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -1.0 * k);
			}
		}
	}

}




TEST(TestDR3, test_transform1)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransform_3(SZ);
	}

	testTransform_3(34);
	testTransform_3(65);
	testTransform_3(63);
	testTransform_3(64);
}



void testTransform_M(int SZ)
{

	auto doubleIt = [](auto x) { return 2.0 * x; };
	VecXX scalar = 3.33;
	VecXX res = scalar;
	transformM(doubleIt, res);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_DOUBLE_EQ(6.66, val);



	std::vector<double> input(SZ, 0.0);
	std::iota(begin(input), end(input), 0.0);

	VecXX testVec(input);

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((j > (x - 0.0001) && (j < x + 0.00001)), x, -x); };
		VecXX res = testVec;
		transformM(onlyJlambda, res);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res[k], k);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -1.0 * k);
			}
		}
	}

}




TEST(TestDR3, test_transformM)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransform_M(SZ);
	}

	testTransform_M(34);
	testTransform_M(65);
	testTransform_M(63);
	testTransform_M(64);
}


void testTransform_M1(int SZ)
{

	auto doubleIt = [](auto x) { return 2.0 * x; };
	VecXX scalar = 3.33;
	VecXX res = scalar;
	transform1(doubleIt, res);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_DOUBLE_EQ(6.66, val);



	std::vector<double> input(SZ, 0.0);
	std::iota(begin(input), end(input), 0.0);

	VecXX testVec(input);

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((j > (x - 0.0001) && (j < x + 0.00001)), x, -x); };
		VecXX res = testVec;
		transform1(onlyJlambda, res);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res[k], k);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -1.0 * k);
			}
		}
	}

}




TEST(TestDR3, test_transform1M)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransform_M1(SZ);
	}

	testTransform_M1(34);
	testTransform_M1(65);
	testTransform_M1(63);
	testTransform_M1(64);
}
