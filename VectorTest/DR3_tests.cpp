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



void testBinaryTransform(int SZ)
{


	auto sumIt = [](auto x,auto y) { return y+ x; };
	VecXX scalar = 3.33;
	// = scalar;
	const VecXX scalarPlusTwo = 2.0 + scalar;
	const VecXX resScalar  =transform(sumIt, scalar, scalarPlusTwo);
	auto val = resScalar.getScalarValue();
	EXPECT_TRUE(resScalar.isScalar());
	EXPECT_DOUBLE_EQ(8.66, val);



	std::vector<double> input(SZ, 0.0);
	std::iota(begin(input), end(input), 0.0);

	const VecXX testVec(input);
	const VecXX testVecPlusTwo = testVec + 2.0;

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJAddlambda = [=](auto x,auto y) { return select((j > (x - 0.0001) && (j < x + 0.00001)), x+y, -(x+y)); };
		//auto onlyJAddlambda_dash = [=](auto x, auto y) { return select((j > (y - 0.0001) && (j < y + 0.00001)), x + y, -(x + y)); };

		VecXX res  = transform(onlyJAddlambda, testVec, testVecPlusTwo);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res[k], 2.0*( k +1.0) );
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -2.0 * (k + 1.0));
			}
		}
	
	}


	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJAddlambda = [=](auto x, auto y) { return select((j > (x - 0.0001) && (j < x + 0.00001)), x + y, -(x + y)); };
		auto onlyJAddlambda_dash = [=](auto x, auto y) { return select((j > (y- 0.0001) && (j < y + 0.00001)), x + y, -(x + y)); };

		const VecXX scalarPlusTwo = 2.0;

		VecXX res = transform(onlyJAddlambda, testVec, scalarPlusTwo);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res[k], k + 2.0);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -1.0 * (k + 2.0));
			}
		}


		VecXX res2 =  transform(onlyJAddlambda_dash, scalarPlusTwo, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res2[k], k + 2.0);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res2[k], -1.0 * (k + 2.0));
			}
		}

		
		VecXX res3 = transform(onlyJAddlambda_dash, 2.0, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res3[k], k + 2.0);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res3[k], -1.0 * (k + 2.0));
			}
		}
		

		VecXX res4 = transform(onlyJAddlambda, testVec,2.0);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res4[k], k + 2.0);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res4[k], -1.0 * (k + 2.0));
			}
		}


	}


}




TEST(TestDR3, test_transform_binary)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testBinaryTransform(SZ);
	}
}




void testBinaryTransform1(int SZ)
{


	auto sumIt = [](auto x, auto y) { return y + x; };
	VecXX scalar = 3.33;
	// = scalar;
	const VecXX scalarPlusTwo = 2.0 + scalar;
	const VecXX resScalar = transform1(sumIt, scalar, scalarPlusTwo);
	auto val = resScalar.getScalarValue();
	EXPECT_TRUE(resScalar.isScalar());
	EXPECT_DOUBLE_EQ(8.66, val);



	std::vector<double> input(SZ, 0.0);
	std::iota(begin(input), end(input), 0.0);

	const VecXX testVec(input);
	const VecXX testVecPlusTwo = testVec + 2.0;

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJAddlambda = [=](auto x, auto y) { return select((j > (x - 0.0001) && (j < x + 0.00001)), x + y, -(x + y)); };
		//auto onlyJAddlambda_dash = [=](auto x, auto y) { return select((j > (y - 0.0001) && (j < y + 0.00001)), x + y, -(x + y)); };

		VecXX res = transform1(onlyJAddlambda, testVec, testVecPlusTwo);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res[k], 2.0 * (k + 1.0));
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -2.0 * (k + 1.0));
			}
		}

	}


	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJAddlambda = [=](auto x, auto y) { return select((j > (x - 0.0001) && (j < x + 0.00001)), x + y, -(x + y)); };
		auto onlyJAddlambda_dash = [=](auto x, auto y) { return select((j > (y - 0.0001) && (j < y + 0.00001)), x + y, -(x + y)); };

		const VecXX scalarPlusTwo = 2.0;

		VecXX res = transform1(onlyJAddlambda, testVec, scalarPlusTwo);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res[k], k + 2.0);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -1.0 * (k + 2.0));
			}
		}


		VecXX res2 = transform1(onlyJAddlambda_dash, scalarPlusTwo, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res2[k], k + 2.0);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res2[k], -1.0 * (k + 2.0));
			}
		}


		VecXX res3 = transform1(onlyJAddlambda_dash, 2.0, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res3[k], k + 2.0);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res3[k], -1.0 * (k + 2.0));
			}
		}


		VecXX res4 = transform1(onlyJAddlambda, testVec, 2.0);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res4[k], k + 2.0);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res4[k], -1.0 * (k + 2.0));
			}
		}


	}


}




TEST(TestDR3, test_transform_binary1)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testBinaryTransform1(SZ);
	}
}





void testBinaryTransformM(int SZ)
{


	auto sumIt = [](auto x, auto y) { return y + x; };
	VecXX scalar = 3.33;
	const VecXX two = 2.0;
	auto resScalar = scalar;
	//const VecXX scalarPlusTwo = 2.0 + scalar;
	transformM(sumIt, resScalar, 2.0);
	auto val = resScalar.getScalarValue();
	EXPECT_TRUE(resScalar.isScalar());
	EXPECT_DOUBLE_EQ(5.33, val);



	std::vector<double> input(SZ, 0.0);
	std::iota(begin(input), end(input), 0.0);


	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJAddlambda = [=](auto x, auto y) { return select((j > (x - 0.0001) && (j < x + 0.00001)), x + y, -(x + y)); };
		auto onlyJAddlambda_dash = [=](auto x, auto y) { return select((j > (y - 0.0001) && (j < y + 0.00001)), x + y, -(x + y)); };

		VecXX res;
		VecXX testVec(input);
		const VecXX testVecPlusTwo = testVec + 2.0;

		transformM(onlyJAddlambda, testVec, testVecPlusTwo);
		res = testVec;

		//std::vector<double> stl = res;

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res[k], 2.0 + 2. *k );
			}
			else
			{
				EXPECT_DOUBLE_EQ(res[k], -(2.0 +2. *k));
			}
		}


		VecXX res1;
		VecXX testVec1(input);
		const VecXX testTwo =  2.0;

		transformM(onlyJAddlambda, testVec1, testTwo);
		res1 = testVec1;

		//std::vector<double> stl = res;

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res1[k], 2.0 +  k);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res1[k], -(2.0 + k));
			}
		}


		VecXX res2;
		VecXX testVec2(input);

		transformM(onlyJAddlambda, testVec2, 2.0);
		res2 = testVec2;

		//std::vector<double> stl = res;

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res2[k], 2.0 + k);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res2[k], -(2.0 + k));
			}
		}


	
		VecXX res3;
		VecXX testVec3(input);

		transformM(onlyJAddlambda_dash,  2.0, testVec3);
		res3 = testVec3;

		std::vector<double> stl = res;

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_DOUBLE_EQ(res3[k], 2.0 + k);
			}
			else
			{
				EXPECT_DOUBLE_EQ(res3[k], -(2.0 + k));
			}
		}
	}
}




TEST(TestDR3, test_transform_binaryM)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testBinaryTransformM(SZ);
	}
}
