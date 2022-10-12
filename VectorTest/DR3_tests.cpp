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

//#include "dr3TestUtil.h"

#include "testNamespace.h"




#include <numeric>
#include <algorithm>
#include <random>



void testTransform_1(int SZ)
{

	auto doubleIt = [](auto x) { return 2.0 * x; };
	VecXX scalar = asNumber(3.33);
	VecXX res = transform(doubleIt, scalar);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(6.66), val);


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);

	
	

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [&](auto x) {
			return select(
			( 
				((asNumber(j) - asNumber(0.0001) ) < x) 
			&&  ((asNumber(j) + asNumber(0.0001) ) > x )
				),x,-x); };
		VecXX res = transform(onlyJlambda, testVec);

		std::vector<Numeric> inspect = res;

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber( -1.0 * k) );
			}
		}
	}


}




TEST(TestDR3, testTransform_1)
{

	testTransform_1(5);


	//for (int SZ = 3; SZ < 5; SZ++)
	//{
	//	testTransform_1(SZ);
	//}


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransform_1(SZ);
	}


	testTransform_1(34);
	testTransform_1(65);
	testTransform_1(63);
	testTransform_1(64);

	/*	*/

}



void testTransform_XX(int SZ)
{

	auto doubleIt = [](auto x) { return 2.0 * x; };
	VecXX scalar = asNumber(3.33);
	VecXX res = transformXX(doubleIt, scalar);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(6.66), val);



	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((asNumber(j) > (x - asNumber(0.001)) && (asNumber(j) < x + asNumber(0.001))), x, -x); };
		VecXX res = transformXX(onlyJlambda, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(-1.0 * k));
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
	VecXX scalar = asNumber(3.33);
	VecXX resScalar;
	transform(doubleIt, scalar, resScalar);
	auto val = resScalar.getScalarValue();
	EXPECT_TRUE(resScalar.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(6.66), val);


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((asNumber(j) > (x - asNumber(0.001)) && (asNumber(j) < x + asNumber(0.001))), x, -x); };
		VecXX res = testVec;
		transform(onlyJlambda, testVec,res);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber( -1.0 * k) );
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
	const VecXX scalar =asNumber(3.33);
	VecXX res = transform1(doubleIt, scalar);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(6.66), val);



	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	const VecXX testVec(input);

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((asNumber(j) > (x - asNumber(0.001)) && (asNumber(j) < x + asNumber(0.001))), x, -x); };
		VecXX res = transform1(onlyJlambda, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(-1.0 * k)) ;
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
	VecXX scalar = asNumber(3.33);
	VecXX res = scalar;
	transformM(doubleIt, res);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(6.66), val);



	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((asNumber(j) > (x - asNumber(0.001)) && (asNumber(j) < x + asNumber(0.001))), x, -x); };
		VecXX res = testVec;
		transformM(onlyJlambda, res);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(-1.0 * k) );
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
	VecXX scalar = asNumber(3.33);
	VecXX res = scalar;
	transform1(doubleIt, res);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(6.66), val);



	std::vector<Numeric> input(SZ, asNumber( 0.0) );
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((asNumber(j) > (x - asNumber(0.001)) && (asNumber(j) < x + asNumber(0.001))), x, -x); };
		VecXX res = testVec;
		transform1(onlyJlambda, res);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(-1.0 * k));
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
	VecXX scalar = asNumber(3.33);
	// = scalar;
	const VecXX scalarPlusTwo = 2.0 + scalar;
	const VecXX resScalar  =transform(sumIt, scalar, scalarPlusTwo);
	auto val = resScalar.getScalarValue();
	EXPECT_TRUE(resScalar.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(8.66), val);



	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	const VecXX testVec(input);
	const VecXX testVecPlusTwo = testVec + 2.0;

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJAddlambda = [=](auto x,auto y) { return select((asNumber(j) > (x - asNumber(0.001)) && (asNumber(j) < x + asNumber(0.001))), x+y, -(x+y)); };
		//auto onlyJAddlambda_dash = [=](auto x, auto y) { return select((j > (y - 0.0001) && (j < y + 0.00001)), x + y, -(x + y)); };

		VecXX res  = transform(onlyJAddlambda, testVec, testVecPlusTwo);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(2.0*( k +1.0)) );
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(-2.0 * (k + 1.0)));
			}
		}
	
	}


	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJAddlambda_dash = [=](auto x, auto y) { return select((asNumber(j) > (y - asNumber(0.001)) && (asNumber(j) < y + asNumber(0.001))), x + y, -(x + y)); };
		auto onlyJAddlambda = [=](auto x, auto y)	   { return select((asNumber(j) > (x - asNumber(0.001)) && (asNumber(j) < x + asNumber(0.001))), x + y, -(x + y)); };


		const VecXX scalarPlusTwo = 2.0;

		VecXX res = transform(onlyJAddlambda, testVec, scalarPlusTwo);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k + 2.0));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(-1.0 * (k + 2.0)));
			}
		}


		VecXX res2 =  transform(onlyJAddlambda_dash, scalarPlusTwo, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res2[k], asNumber( k + 2.0));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res2[k], asNumber( -1.0 * (k + 2.0)));
			}
		}

		
		VecXX res3 = transform(onlyJAddlambda_dash, 2.0, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res3[k], asNumber( k + 2.0));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res3[k], asNumber (-1.0 * (k + 2.0)));
			}
		}
		

		VecXX res4 = transform(onlyJAddlambda, testVec,2.0);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res4[k], asNumber( k + 2.0));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res4[k], asNumber (-1.0 * (k + 2.0)));
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
	VecXX scalar = asNumber(3.33);
	// = scalar;
	const VecXX scalarPlusTwo = 2.0 + scalar;
	const VecXX resScalar = transform1(sumIt, scalar, scalarPlusTwo);
	auto val = resScalar.getScalarValue();
	EXPECT_TRUE(resScalar.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(8.66), val);



	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	const VecXX testVec(input);
	const VecXX testVecPlusTwo = testVec + 2.0;

	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJAddlambda = [=](auto x, auto y) { return select((asNumber(j) > (x - asNumber(0.001)) && (asNumber(j) < x + asNumber(0.001))), x + y, -(x + y)); };
		//auto onlyJAddlambda_dash = [=](auto x, auto y) { return select((j > (y - 0.0001) && (j < y + 0.00001)), x + y, -(x + y)); };

		VecXX res = transform1(onlyJAddlambda, testVec, testVecPlusTwo);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(2.0 * (k + 1.0)));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber (-2.0 * (k + 1.0)));
			}
		}

	}


	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJAddlambda = [=](auto x, auto y) { return select((asNumber(j) > (x - asNumber(0.001)) && (asNumber(j) < x + asNumber(0.001))), x + y, -(x + y)); };
		auto onlyJAddlambda_dash = [=](auto x, auto y) { return select((asNumber(j) > (y - asNumber(0.001)) && (asNumber(j) < y + asNumber(0.001))), x + y, -(x + y)); };

		const VecXX scalarPlusTwo = 2.0;

		VecXX res = transform1(onlyJAddlambda, testVec, scalarPlusTwo);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k + 2.0));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber (-1.0 * (k + 2.0)));
			}
		}


		VecXX res2 = transform1(onlyJAddlambda_dash, scalarPlusTwo, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res2[k], asNumber(k + 2.0));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res2[k], asNumber (-1.0 * (k + 2.0)));
			}
		}


		VecXX res3 = transform1(onlyJAddlambda_dash, 2.0, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res3[k], asNumber( k + 2.0));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res3[k], asNumber (-1.0 * (k + 2.0)));
			}
		}


		VecXX res4 = transform1(onlyJAddlambda, testVec, 2.0);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res4[k], asNumber(k + 2.0));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res4[k], asNumber (-1.0 * (k + 2.0)));
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
	VecXX scalar = asNumber(3.33);
	const VecXX two = asNumber(2.0);
	auto resScalar = scalar;
	
	transformM(sumIt, resScalar, asNumber(2.0) );
	auto val = resScalar.getScalarValue();
	EXPECT_TRUE(resScalar.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(5.33), val);



	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));


	for (int j = 0; j < SZ; ++j)
	{
		auto onlyJAddlambda = [=](auto x, auto y) { return select((asNumber(j) > (x - asNumber(0.001)) && (asNumber(j) < x + asNumber(0.001))), x + y, -(x + y)); };
		auto onlyJAddlambda_dash = [=](auto x, auto y) { return select((asNumber(j) > (y - asNumber(0.001)) && (asNumber(j) < y + asNumber(0.001))), x + y, -(x + y)); };

		VecXX res;
		VecXX testVec(input);
		const VecXX testVecPlusTwo = testVec + 2.0;

		transformM(onlyJAddlambda, testVec, testVecPlusTwo);
		res = testVec;

		
		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber( 2.0 + 2. *k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber (-(2.0 +2. *k)));
			}
		}


		VecXX res1;
		VecXX testVec1(input);
		const VecXX testTwo =  2.0;

		transformM(onlyJAddlambda, testVec1, testTwo);
		res1 = testVec1;

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res1[k], asNumber(2.0 +  k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res1[k], asNumber (-(2.0 + k)));
			}
		}


		VecXX res2;
		VecXX testVec2(input);

		transformM(onlyJAddlambda, testVec2, 2.0);
		res2 = testVec2;


		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res2[k], asNumber(2.0 + k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res2[k], asNumber (-(2.0 + k)));
			}
		}


	
		VecXX res3;
		VecXX testVec3(input);

		transformM(onlyJAddlambda_dash,  2.0, testVec3);
		res3 = testVec3;

		std::vector<Numeric> stl = res;

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res3[k], asNumber( 2.0 + k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res3[k], asNumber (-(2.0 + k)));
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







void testSelect(int SZ)
{

	Numeric j = 666.;

	auto equalsJay = [&](auto x) { return (asNumber(j) > (x - asNumber(0.001))) && (asNumber(j) < (x + asNumber(0.001))); };
	const VecXX trueValue = asNumber(1.);
	const VecXX falseValue = asNumber(-1.);
	const VecXX testValue = asNumber(666.);
	VecXX resScalar = select(equalsJay, testValue, trueValue, falseValue);

	auto val = resScalar.getScalarValue();
	EXPECT_TRUE(resScalar.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(1.0), val);

	const VecXX testValueFalse = asNumber(999.);
	VecXX resScalarFlse = select(equalsJay, testValueFalse, trueValue, falseValue);

	val = resScalarFlse.getScalarValue();
	EXPECT_TRUE(resScalarFlse.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(-1.0), val);

	auto SQR = [](auto x) { return x * x; };
	auto CUBE = [](auto x) { return x * x * x; };

	VecXX resSelScalarLambda = selectTransform(equalsJay, testValue, SQR, CUBE);
	val = resSelScalarLambda.getScalarValue();
	EXPECT_TRUE(resSelScalarLambda.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(666.*666.), val);


	VecXX resSelScalarLambdaFlse = selectTransform(equalsJay, falseValue, SQR, CUBE);
	val = resSelScalarLambdaFlse.getScalarValue();
	EXPECT_TRUE(resSelScalarLambdaFlse.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(-1.), val);



	for (int j = 0; j < SZ; ++j)
	{

		std::vector<Numeric> input(SZ, asNumber(0.0));
		std::iota(begin(input), end(input), asNumber(0.0));

		auto equalsJay = [&](auto x) { return (asNumber(j) > (x - asNumber(0.001))) && (asNumber(j) < (x + asNumber(0.001))); };

		const VecXX trueValues(input);
		const VecXX falseValues = -trueValues;
		const VecXX testValues(input);
		
		VecXX res = select(equalsJay, testValues, trueValues, falseValues);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], -asNumber(k));
			}
		}
	

		VecXX resInv = select(equalsJay, testValues,  falseValues, trueValues);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(resInv[k], asNumber(-k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(resInv[k], asNumber(k));
			}
		}


		VecXX resScalarChoice = select(equalsJay, testValues, asNumber(-666.66), asNumber(999.99) );
		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(resScalarChoice[k], asNumber (-666.66) );
			}
			else
			{
				EXPECT_NUMERIC_EQ(resScalarChoice[k], asNumber( 999.99));
			}
		}


		auto SQR = [](auto x) { return x * x; };
		auto CUBE = [](auto x) { return x * x * x; };

		VecXX resSelLambda = selectTransform(equalsJay, testValues, SQR, CUBE );
		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(resSelLambda[k], SQR(asNumber(k)));
			}
			else
			{
				EXPECT_NUMERIC_EQ(resSelLambda[k], CUBE(asNumber(k)));
			}
		}

		
	}
	
}




TEST(TestDR3, test_select)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testSelect(SZ);
	}
}








void test_FilterTransform(int SZ)
{

	Numeric j = asNumber(666.);

	auto equalsJay = [&](auto x) { return (asNumber(j) > (x - asNumber(0.001))) && (asNumber(j) < (x + asNumber(0.001))); };
	const VecXX trueValue = 1.;
	const VecXX falseValue = -1.;
	const VecXX testValue = 666.;

	auto SQR = [](auto x) { return x * x; };
	auto CUBE = [](auto x) { return x * x * x; };

	VecXX resSelScalarLambda = filterTransform(equalsJay, testValue, SQR, CUBE);
	auto val = resSelScalarLambda.getScalarValue();
	EXPECT_TRUE(resSelScalarLambda.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(666. * 666.), val);


	VecXX resSelScalarLambdaFlse = filterTransform(equalsJay, falseValue, SQR, CUBE);
	val = resSelScalarLambdaFlse.getScalarValue();
	EXPECT_TRUE(resSelScalarLambdaFlse.isScalar());
	EXPECT_NUMERIC_EQ(asNumber (-1.), val);



	for (int j = 0; j < SZ; ++j)
	{

		std::vector<Numeric> input(SZ, asNumber(0.0));
		std::iota(begin(input), end(input), asNumber(0.0));

		auto equalsJay = [&](auto x) { return (asNumber(j) > (x - asNumber(0.001) )) && (asNumber(j) < (x + asNumber(0.001))); };

		const VecXX testValues(input);

		auto SQR = [](auto x) { return x * x; };
		auto CUBE = [](auto x) { return x * x * x; };

		VecXX resSelLambda = filterTransform(equalsJay, testValues, SQR, CUBE);
		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(resSelLambda[k], SQR(asNumber(k)));
			}
			else
			{
				EXPECT_NUMERIC_EQ(resSelLambda[k], CUBE(asNumber(k)));
			}
		}
	}

}




TEST(TestDR3, test_filterTransform)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		test_FilterTransform(SZ);
	}

}





void testTransform_V(int SZ)
{
	
	auto doubleIt = [](auto x) { return 2.0 * x; };
	VecXX scalar = asNumber(3.33);
	VecVW res = transformV(doubleIt, scalar);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(6.66), val);


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);


	for (Numeric j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((asNumber(j) > (x - VecXX::INS(asNumber(0.001)) ) && (asNumber(j) < x + VecXX::INS(asNumber(0.001)) )), x, -x); };
		VecVW res = transformV(onlyJlambda, testVec);

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber (-1.0 * asNumber(k)));
			}
		}
	}


}




TEST(TestDR3, testTransform_view)
{

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransform_V(SZ);
	}


	testTransform_V(34);
	testTransform_V(63);
	testTransform_V(64);


}




void test_FilterTransform_View(int SZ)
{
	auto copyVecToView = [](auto x) {return x; };
	Numeric j = asNumber(666.);

	auto equalsJay = [&](auto x) { return (asNumber(j) > (x - asNumber(0.001))) && (asNumber(j) < (x + asNumber(0.001))); };
	const VecXX trueValue = 1.;
	const VecXX falseValue = -1.;
	const VecXX testValue = 666.;

	auto SQR = [](auto x) { return x * x; };
	auto CUBE = [](auto x) { return x * x * x; };
	VecVW testVwScal = transformV(copyVecToView, testValue);

	try
	{
		VecVW resSelScalarLambda = filterTransform(equalsJay, testVwScal, SQR, CUBE);
		
	}
	catch (const std::exception& ex)
	{
		bool throws = true;
		EXPECT_TRUE(throws);
		ignore(ex);
	}
	


	for (int j = 0; j < SZ; ++j)
	{

		std::vector<Numeric> input(SZ, asNumber(0.0));
		std::iota(begin(input), end(input), asNumber(0.0));
		auto equalsJay = [&](auto x) { return (asNumber(j) > (x - asNumber(0.001))) && (asNumber(j) < (x + asNumber(0.001))); };

		const VecXX testValues(input);
		VecVW testVw = transformV(copyVecToView, testValues);

		std::vector<Numeric> stdVec = testVw;
		auto idz = testVw.getIndex();

		auto SQR = [](auto x) { return x * x; };
		auto CUBE = [](auto x) { return x * x * x; };

		VecVW resSelLambda = filterTransform(equalsJay, testVw, SQR, CUBE);
		stdVec = resSelLambda;

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(resSelLambda[k], SQR(asNumber(k)));
			}
			else
			{
				EXPECT_NUMERIC_EQ(resSelLambda[k], CUBE(asNumber(k)));
			}
		}
	}

}




TEST(TestDR3, test_filterTransformView)
{

	//test_FilterTransform_View(3); //ok
	//test_FilterTransform_View(4); //ok
	//test_FilterTransform_View(5); //ok
	//test_FilterTransform_View(6); //ok
	//test_FilterTransform_View(7); //ok
	//test_FilterTransform_View(8); //ok
	//test_FilterTransform_View(9); // ok
   // test_FilterTransform_View(10); // XXXX
	//test_FilterTransform_View(11); // XXXX
	//test_FilterTransform_View(12); // XXXX

	/*
	for (int SZ = 9; SZ < 11; SZ++)
	{
		test_FilterTransform_View(SZ);
	}
	*/

	test_FilterTransform_View(16);
	test_FilterTransform_View(34);
	test_FilterTransform_View(63);
	test_FilterTransform_View(64);

	for (int SZ = 3; SZ < 33; SZ++)
	{
		test_FilterTransform_View(SZ);
	}
	
}








void test_Transform_V(int SZ)
{


	auto doubleIt = [](auto x) { return 2.0 * x; };
	VecXX scalar = asNumber(3.33);
	VecVW res = transformV(doubleIt, scalar);
	auto val = res.getScalarValue();
	EXPECT_TRUE(res.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(6.66), val);


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);
	auto copyLambda = [&](auto x) {return x; };
	VecVW inputView = transformV(copyLambda, testVec);	




	for (Numeric j = 0; j < SZ; ++j)
	{
		auto onlyJlambda = [=](auto x) { return select((asNumber(j) > (x - VecXX::INS(asNumber(0.01))) && (asNumber(j) < x + VecXX::INS(asNumber(0.01)))), x, -x); };


		auto res  =transform(onlyJlambda, inputView);
		

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], -asNumber(k));
			}
		}
	}


}




TEST(TestDR3, testTransformA_view)
{

	test_Transform_V(3);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		test_Transform_V(SZ);
	}

	test_Transform_V(34);
	test_Transform_V(63);
	test_Transform_V(64);
}




//transform and modify a view 
void testTransformM_V(int SZ)
{

	auto doubleIt = [](auto x) { return 2.0 * x; };
	VecXX scalar = asNumber(3.33);
	VecVW doubled = transformV(doubleIt, scalar);
	transformM(doubleIt,doubled );
	auto val = doubled.getScalarValue();
	EXPECT_TRUE(doubled.isScalar());
	EXPECT_NUMERIC_EQ(asNumber(6.66*2.0), val);
	

	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);

	auto copyLambda = [&](auto x) {return x; };
	VecVW inputView = transformV(copyLambda, testVec);


	for (Numeric j = 0; j < SZ; ++j)
	{

		auto onlyJlambda = [=](auto x) { return select((j == x), x, x* x); };

	
		auto copyMod = inputView;

		transformM(onlyJlambda, copyMod);
		auto res = copyMod;

		auto index2 = res.getIndex();
		std::vector<Numeric>  st2 = res;


		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k));
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k* k));
			}
		}
	}


}


TEST(TestDR3, testTransformM_view)
{

	testTransformM_V(3);
	testTransformM_V(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransformM_V(SZ);
	}

	testTransformM_V(34);
	testTransformM_V(63);
	testTransformM_V(64);
	

}







//transform and modify a view 
void testTransformWrite_Vw(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber( 0.0));

	VecXX testVec(input);

	auto copyLambda = [&](auto x) {return x; };
	VecVW inputView = transformV(copyLambda, testVec);


	for (Numeric  j = 0; j < SZ; ++j)
	{

		auto onlyJlambda = [=](auto x) { return j == x; }; 
		auto SQR = [](auto x) {return x * x; };
		auto fltrVwOnlyJ = filter(onlyJlambda, inputView);
		auto copyInputVec = testVec;
		transformWrite(SQR, fltrVwOnlyJ, copyInputVec);
		auto res = copyInputVec;

		for (int k = 0; k < SZ; k++)
		{
			if (k == j)
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber( k*k) );
			}
			else
			{
				EXPECT_NUMERIC_EQ(res[k], asNumber(k));
			}
		}
	}
}




TEST(TestDR3, testTransformWrite_view)
{

	testTransformWrite_Vw(3);
	testTransformWrite_Vw(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransformWrite_Vw(SZ);
	}

	testTransformWrite_Vw(34);
	testTransformWrite_Vw(63);
	testTransformWrite_Vw(64);
	testTransformWrite_Vw(65);


}



// TO DO filters

// excluding filter by vec of booles 
//



void testFilterVw(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);

	auto copyLambda = [&](auto x) {return x; };
	VecVW inputView = transformV(copyLambda, testVec);

	auto allowAll = [](auto x) {return x==x; };


	auto resAll = filter(allowAll, inputView);
	
	
	for (int i =0; i < resAll.size();++i)
	{
		EXPECT_NUMERIC_EQ(resAll[i], inputView[i]);
	}
	
	auto allowNone = [](auto x) {return x != x; };
	
	auto resNone = filter(allowNone, inputView);

	EXPECT_EQ(0, resNone.size());



	for (Numeric j = 0; j < SZ; ++j)
	{

		auto onlyJlambda = [=](auto x) { return j == x; };
		auto res = filter(onlyJlambda, inputView);

		EXPECT_NUMERIC_EQ(res[0], j);

	}
}




TEST(TestDR3, testFilter_view)
{

	testFilterVw(3);
	testFilterVw(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testFilterVw(SZ);
	}

	testFilterVw(34);
	testFilterVw(63);
	testFilterVw(64);
	testFilterVw(65);

}



void testFilterVecXX(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);

	auto copyLambda = [&](auto x) {return x; };
	//VecVW inputView = transformV(copyLambda, testVec);

	auto allowAll = [](auto x) {return x == x; };


	auto resAll = filter(allowAll, testVec);


	for (int i = 0; i < resAll.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(resAll[i], testVec[i]);
	}

	auto allowNone = [](auto x) {return x != x; };

	auto resNone = filter(allowNone, testVec);

	EXPECT_EQ(0, resNone.size());



	for (Numeric j = 0; j < SZ; ++j)
	{

		auto onlyJlambda = [=](auto x) { return j == x; };
		auto res = filter(onlyJlambda, testVec);

		EXPECT_NUMERIC_EQ(res[0], j);

	}
}




TEST(TestDR3, testFilter_vec)
{

	testFilterVecXX(3);
	testFilterVecXX(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testFilterVecXX(SZ);
	}

	testFilterVecXX(34);
	testFilterVecXX(63);
	testFilterVecXX(64);
	testFilterVecXX(65);

}












void testCountedFilterVw(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber( 0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);

	auto copyLambda = [&](auto x) {return x; };
	VecVW inputView = transformV(copyLambda, testVec);

	auto allowAll = [](auto x) {return x == x; };


	auto resAll = countedFilter(allowAll, inputView,SZ+100);
	auto resAll1 = countedFilter(allowAll, inputView, SZ );
	auto resAll2 = countedFilter(allowAll, inputView, SZ-1);


	for (int i = 0; i < resAll.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(resAll[i], inputView[i]);
	}

	for (int i = 0; i < resAll1.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(resAll1[i], inputView[i]);
	}


	for (int i = 0; i < resAll2.size()-1; ++i)
	{
		EXPECT_NUMERIC_EQ(resAll2[i], inputView[i]);
	}



	auto allowNone = [](auto x) {return x != x; };

	auto resNone = countedFilter(allowNone, inputView,SZ);

	EXPECT_EQ(0, resNone.size());



	for (int j = 0; j < SZ; ++j)
	{

		auto res = countedFilter(allowAll, inputView,j);
		EXPECT_EQ(res.size(), j);

		for (int i = 0; i < res.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(res[i], asNumber(i));
		}

	}
}

TEST(TestDR3, testCountedFilter_view)
{

	testCountedFilterVw(3);
	testCountedFilterVw(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testCountedFilterVw(SZ);
	}

	testCountedFilterVw(34);
	testCountedFilterVw(63);
	testCountedFilterVw(64);
	testCountedFilterVw(65);

}





void testCountedFilterVecXX(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);

	//auto copyLambda = [&](auto x) {return x; };
	//VecVW inputView = transformV(copyLambda, testVec);

	auto allowAll = [](auto x) {return x == x; };


	auto resAll = countedFilter(allowAll, testVec, SZ + 100);
	auto resAll1 = countedFilter(allowAll, testVec, SZ);
	auto resAll2 = countedFilter(allowAll, testVec, SZ - 1);


	for (int i = 0; i < resAll.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(resAll[i], testVec[i]);
	}

	for (int i = 0; i < resAll1.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(resAll1[i], testVec[i]);
	}


	for (int i = 0; i < resAll2.size() - 1; ++i)
	{
		EXPECT_NUMERIC_EQ(resAll2[i], testVec[i]);
	}



	auto allowNone = [](auto x) {return x != x; };

	auto resNone = countedFilter(allowNone, testVec, SZ);

	EXPECT_EQ(0, resNone.size());



	for (int j = 0; j < SZ; ++j)
	{

		auto res = countedFilter(allowAll, testVec, j);
		EXPECT_EQ(res.size(), j);

		for (int i = 0; i < res.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(res[i], asNumber( i));
		}

	}
}

TEST(TestDR3, testCountedFilter_vec)
{

	testCountedFilterVecXX(3);
	testCountedFilterVecXX(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testCountedFilterVecXX(SZ);
	}

	testCountedFilterVecXX(34);
	testCountedFilterVecXX(63);
	testCountedFilterVecXX(64);
	testCountedFilterVecXX(65);

}



//binary filter to do ??


void testBimaryFilterVecXX(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0) );

	VecXX testVec(input);

	
	auto allowAll = [](auto x) {return x == x; };

	auto  tpl = binaryFilter(allowAll, testVec);
	
	auto trueVw = std::get<0>(tpl);
	auto falseVw = std::get<1>(tpl);

	for (int i = 0; i < trueVw.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(trueVw[i], testVec[i]);
	}

	EXPECT_EQ(falseVw.size(), 0);
	EXPECT_EQ(trueVw.size(), SZ);

	auto allowNone = [](auto x) {return x != x; };
	auto  tplNone = binaryFilter(allowNone, testVec);

	trueVw = std::get<0>(tplNone);
	falseVw = std::get<1>(tplNone);


	for (int i = 0; i < falseVw.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(falseVw[i], testVec[i]);
	}

	EXPECT_EQ(trueVw.size(), 0);
	EXPECT_EQ(falseVw.size(), SZ);


	for (Numeric j = 0; j < SZ; ++j)
	{
		auto allowj = [&](auto x) {return x == j; };
		auto res = binaryFilter(allowj, testVec);
	
		auto trueVw = std::get<0>(res);
		auto falseVw = std::get<1>(res);
		EXPECT_EQ(trueVw.size(), 1);

		for (int i = 0; i < testVec.size()-1; ++i)
		{
			if (i < j)
			{
				EXPECT_NUMERIC_EQ(falseVw[i], asNumber(i));
			}
			else if (i > j)
			{
				EXPECT_NUMERIC_EQ(falseVw[i], asNumber(i + 1.));
			}
			else
			{
				EXPECT_NUMERIC_EQ(trueVw[0], asNumber(i) );
			}
		}
	}
}

TEST(TestDR3, testBinaryFilter_vec)
{

	testBimaryFilterVecXX(3);
	testBimaryFilterVecXX(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testBimaryFilterVecXX(SZ);
	}

	testBimaryFilterVecXX(34);
	testBimaryFilterVecXX(63);
	testBimaryFilterVecXX(64);
	testBimaryFilterVecXX(65);

}






void testBinaryFilterVecView(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);

	auto copyLambda = [&](auto x) {return x; };
	VecVW inputView = transformV(copyLambda, testVec);

	auto allowAll = [](auto x) {return x == x; };



	auto  tpl = binaryFilter(allowAll, inputView);

	auto trueVw = std::get<0>(tpl);
	auto falseVw = std::get<1>(tpl);

	for (int i = 0; i < trueVw.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(trueVw[i], testVec[i]);
	}

	EXPECT_EQ(falseVw.size(), 0);
	EXPECT_EQ(trueVw.size(), SZ);

	auto allowNone = [](auto x) {return x != x; };
	auto  tplNone = binaryFilter(allowNone, inputView);

	trueVw = std::get<0>(tplNone);
	falseVw = std::get<1>(tplNone);


	for (int i = 0; i < falseVw.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(falseVw[i], testVec[i]);
	}

	EXPECT_EQ(trueVw.size(), 0);
	EXPECT_EQ(falseVw.size(), SZ);


	for (Numeric j = 0; j < SZ; ++j)
	{
		auto allowj = [&](auto x) {return x == j; };
		auto res = binaryFilter(allowj, inputView);

		auto trueVw = std::get<0>(res);
		auto falseVw = std::get<1>(res);
		EXPECT_EQ(trueVw.size(), 1);

		for (int i = 0; i < testVec.size() - 1; ++i)
		{
			if (i < j)
			{
				EXPECT_NUMERIC_EQ(falseVw[i], asNumber( i) );
			}
			else if (i > j)
			{
				EXPECT_NUMERIC_EQ(falseVw[i], asNumber(i + 1.));
			}
			else
			{
				EXPECT_NUMERIC_EQ(trueVw[0], asNumber( i));
			}
		}
	}
}

TEST(TestDR3, testBinaryFilter_Vw)
{

	testBinaryFilterVecView(3);
	testBinaryFilterVecView(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testBinaryFilterVecView(SZ);
	}

	testBinaryFilterVecView(34);
	testBinaryFilterVecView(63);
	testBinaryFilterVecView(64);
	testBinaryFilterVecView(65);
}





void testSparseTransform_Vec(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	VecXX testVec(input);
	VecXX resultVec(input);
	auto resultVec1 = resultVec;


	auto allowAll = [](auto x) {return x == x; };
	auto allowANone = [](auto x) {return x != x; };

	auto SQR = [](auto x) { return x * x; };


	sparseTransform(testVec, resultVec, SQR, allowAll);

	for (int i = 0; i < resultVec.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(resultVec[i], testVec[i]* testVec[i]);
	}



	sparseTransform(testVec, resultVec1, SQR, allowANone);

	for (int i = 0; i < resultVec.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(resultVec1[i], testVec[i]);
	}


	//transform only one point

	for (Numeric j = 0; j < SZ; ++j)
	{
		auto updateAt_j = [&](auto x) {return x == j; };

		auto updateVec = testVec;

		sparseTransform(testVec, updateVec, SQR, updateAt_j);


		for (int i = 0; i < testVec.size() ; ++i)
		{
			if (i < j)
			{
				EXPECT_NUMERIC_EQ(updateVec[i], testVec[i]);
			}
			else if (i > j)
			{
				EXPECT_NUMERIC_EQ(updateVec[i], testVec[i]);
			}
			else
			{
				EXPECT_NUMERIC_EQ(updateVec[i], testVec[i] * testVec[i]);
			}
		}
	}

}

TEST(TestDR3, testSparseTransform_VecXX)
{

	testSparseTransform_Vec(3);
	testSparseTransform_Vec(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testSparseTransform_Vec(SZ);
	}

	testSparseTransform_Vec(34);
	testSparseTransform_Vec(63);
	testSparseTransform_Vec(64);
	testSparseTransform_Vec(65);
}


/////////////////////  test reductions //////////////////

void testReduce_Vec(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(begin(input), end(input), g);
	VecXX testVec(input);
	VecXX negTestVec = -testVec + asNumber( SZ / 2);


	auto SUM = [](auto x, auto y) { return x + y; };
	auto MAX = [](auto x, auto y) { return iff((x > y), x, y); };
	auto MIN = [](auto x, auto y) { return iff((x < y), x, y); };


	auto resSUM =reduce(testVec, SUM);
	auto expectedSum = asNumber((SZ - 1) * SZ / 2.0);
	EXPECT_NUMERIC_EQ(resSUM, asNumber(expectedSum));

	auto resMAX = reduce(testVec, MAX);
	EXPECT_NUMERIC_EQ(resMAX, asNumber(SZ-1));

	auto resMIN = reduce(negTestVec, MIN);
	EXPECT_NUMERIC_EQ(resMIN, asNumber (-1.0 * (SZ-1) + SZ / 2));






}

TEST(TestDR3, testReduce_VecXX)
{

	testReduce_Vec(3);
	testReduce_Vec(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testReduce_Vec(SZ);
	}

	testReduce_Vec(34);
	testReduce_Vec(63);
	testReduce_Vec(64);
	testReduce_Vec(65);
}




void testReduce1_Vec(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(begin(input), end(input), g);
	VecXX testVec(input);
	VecXX negTestVec = asNumber(SZ / 2) -testVec;


	auto SUM = [](auto x, auto y) { return x + y; };
	auto MAX = [](auto x, auto y) { return iff((x > y), x, y); };
	auto MIN = [](auto x, auto y) { return iff((x < y), x, y); };


	auto resSUM = reduce1(testVec, SUM);
	auto expectedSum = asNumber((SZ - 1) * SZ / 2.0);
	EXPECT_NUMERIC_EQ(resSUM, expectedSum);

	auto resMAX = reduce1(testVec, MAX);
	EXPECT_NUMERIC_EQ(resMAX, asNumber(SZ - 1));

	auto resMIN = reduce1(negTestVec, MIN);
	EXPECT_NUMERIC_EQ(resMIN, asNumber (-1.0 * (SZ - 1) + SZ / 2));






}

TEST(TestDR3, testReduce1_VecXX)
{

	testReduce1_Vec(3);
	testReduce1_Vec(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testReduce1_Vec(SZ);
	}

	testReduce1_Vec(34);
	testReduce1_Vec(63);
	testReduce1_Vec(64);
	testReduce1_Vec(65);
}




void testReduce_Vw(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(begin(input), end(input), g);
	VecXX testVec(input);
	VecXX negTestVec = -testVec + asNumber(SZ / 2);


	auto SUM = [](auto x, auto y) { return x + y; };
	auto MAX = [](auto x, auto y) { return iff((x > y), x, y); };
	auto MIN = [](auto x, auto y) { return iff((x < y), x, y); };


	auto copyLambda = [&](auto x) {return x; };
	VecVW testVw = transformV(copyLambda, testVec);
	VecVW negTestVw = transformV(copyLambda, negTestVec);


	auto resSUM = reduce(testVw, SUM);
	auto expectedSum = asNumber((SZ - 1) * SZ / 2.0);
	EXPECT_NUMERIC_EQ(resSUM, expectedSum);

	auto resMAX = reduce(testVw, MAX);
	EXPECT_NUMERIC_EQ(resMAX, asNumber(SZ - 1));

	auto resMIN = reduce(negTestVw, MIN);
	EXPECT_NUMERIC_EQ(resMIN, asNumber (-1.0 * (SZ - 1) + SZ / 2));






}

TEST(TestDR3, testReduce_Views)
{

	testReduce_Vw(3);
	testReduce_Vw(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testReduce_Vw(SZ);
	}

	testReduce_Vw(34);
	testReduce_Vw(63);
	testReduce_Vw(64);
	testReduce_Vw(65);
}




//////////////////////////////////////////////

/////////////////////  test reductions //////////////////

void testTransformReduceUnitary_Vec(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(begin(input), end(input), g);
	VecXX testVec(input);
	VecXX negTestVec = -testVec + asNumber(SZ / 2);


	auto SUM = [](auto x, auto y) { return x + y; };
	auto MAX = [](auto x, auto y) { return iff((x > y), x, y); };
	auto MIN = [](auto x, auto y) { return iff((x < y), x, y); };

	auto CPY = [](auto x) { return x; };

	


	auto resSUM = transformReduce(testVec, CPY,SUM);
	auto expectedSum = asNumber((SZ - 1) * SZ / 2.0);
	EXPECT_NUMERIC_EQ(resSUM, expectedSum);

	auto resMAX = transformReduce(testVec, CPY, MAX);
	EXPECT_NUMERIC_EQ(resMAX, asNumber(SZ - 1));

	auto resMIN = transformReduce(negTestVec, CPY,MIN);
	EXPECT_NUMERIC_EQ(resMIN, asNumber (-1.0 * (SZ - 1) + SZ / 2));


	auto NEGATE = [](auto x) { return -x; };




	auto resSUMNeg = transformReduce(testVec, NEGATE, SUM);
	expectedSum = asNumber(-(SZ - 1) * SZ / 2.0);
	EXPECT_DOUBLE_EQ(resSUMNeg, expectedSum);

	auto resMAXNeg = transformReduce(testVec, NEGATE, MAX);
	EXPECT_DOUBLE_EQ(resMAXNeg, 0);

	auto resMINNeg = transformReduce(negTestVec, NEGATE, MIN);
	//EXPECT_DOUBLE_EQ(resMIN, -1.0 * (SZ - 1) + SZ / 2);

	/*
	auto SQR = [](auto x) { return x*x; };

	auto resSUMNeg = transformReduce(testVec, SQR, SUM);
	expectedSum = -(SZ - 1) * SZ / 2.0;
	EXPECT_DOUBLE_EQ(resSUMNeg, expectedSum);

	auto resMAXNeg = transformReduce(testVec, SQR, MAX);
	EXPECT_DOUBLE_EQ(resMAXNeg, 0);

	auto resMINNeg = transformReduce(negTestVec, SQR, MIN);
	//EXPECT_DOUBLE_EQ(resMIN, -1.0 * (SZ - 1) + SZ / 2);

	*/

}

TEST(TestDR3, testTransformReduceUnitary_VecXX)
{

	testTransformReduceUnitary_Vec(3);
	testTransformReduceUnitary_Vec(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransformReduceUnitary_Vec(SZ);
	}

	testTransformReduceUnitary_Vec(34);
	testTransformReduceUnitary_Vec(63);
	testTransformReduceUnitary_Vec(64);
	testTransformReduceUnitary_Vec(65);
}




void testTransformReduceUnitary_Vw(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(begin(input), end(input), g);
	VecXX testVec(input);
	VecXX negTestVec = -testVec + asNumber(SZ /2);


	auto SUM = [](auto x, auto y) { return x + y; };
	auto MAX = [](auto x, auto y) { return iff((x > y), x, y); };
	auto MIN = [](auto x, auto y) { return iff((x < y), x, y); };

	auto CPY = [](auto x) { return x; };


	auto testVw = transformV(CPY, testVec);
	auto negTestVw = transform(CPY, negTestVec);

	auto resSUM = transformReduce(testVw, CPY, SUM);
	auto expectedSum = asNumber((SZ - 1) * SZ / 2.0);
	EXPECT_NUMERIC_EQ(resSUM, expectedSum);

	auto resMAX = transformReduce(testVw, CPY, MAX);
	EXPECT_NUMERIC_EQ(resMAX, asNumber(SZ - 1));

	auto resMIN = transformReduce(negTestVw, CPY, MIN);
	EXPECT_NUMERIC_EQ(resMIN, asNumber (-1.0 * (SZ - 1) + SZ / 2));


	auto NEGATE = [](auto x) { return -x; };




	auto resSUMNeg = transformReduce(testVw, NEGATE, SUM);
	expectedSum = asNumber(-(SZ - 1) * SZ / 2.0);
	EXPECT_NUMERIC_EQ(resSUMNeg, expectedSum);

	auto resMAXNeg = transformReduce(testVw, NEGATE, MAX);
	EXPECT_NUMERIC_EQ(resMAXNeg, asNumber(0));

	auto resMINNeg = transformReduce(negTestVw, NEGATE, MIN);
	//EXPECT_DOUBLE_EQ(resMIN, -1.0 * (SZ - 1) + SZ / 2);

	/*
	auto SQR = [](auto x) { return x*x; };

	auto resSUMNeg = transformReduce(testVec, SQR, SUM);
	expectedSum = -(SZ - 1) * SZ / 2.0;
	EXPECT_DOUBLE_EQ(resSUMNeg, expectedSum);

	auto resMAXNeg = transformReduce(testVec, SQR, MAX);
	EXPECT_DOUBLE_EQ(resMAXNeg, 0);

	auto resMINNeg = transformReduce(negTestVec, SQR, MIN);
	//EXPECT_DOUBLE_EQ(resMIN, -1.0 * (SZ - 1) + SZ / 2);

	*/

}

TEST(TestDR3, testTransformReduceUnitary_View)
{

	testTransformReduceUnitary_Vw(3);
	testTransformReduceUnitary_Vw(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransformReduceUnitary_Vw(SZ);
	}

	testTransformReduceUnitary_Vw(34);
	testTransformReduceUnitary_Vw(63);
	testTransformReduceUnitary_Vw(64);
	testTransformReduceUnitary_Vw(65);
}

/////////////////////////////////


void testTransformReduceBinary_Vw(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(begin(input), end(input), g);
	VecXX testVec(input);

	auto testVecPlusTwo = testVec + asNumber(2.0);

	auto SUM = [](auto x, auto y) { return x + y; };
	auto MULT = [](auto x, auto y) { return x * y; };
	
	auto CPY = [](auto x) { return x; };


	auto testVw = transformV(CPY, testVec);
	auto testPlusTwoVw = transformV(CPY, testVecPlusTwo);

	auto inputPlusTwo = input;
	for(auto & x: inputPlusTwo)
	{
		x += asNumber(2.0);
	}

	auto resSUM = transformReduce(testVw, testPlusTwoVw, MULT, SUM);
	auto ipRes = std::inner_product(begin(input), end(input), begin(inputPlusTwo), asNumber(.0));

	EXPECT_NUMERIC_EQ(resSUM, ipRes);

}

TEST(TestDR3, testTransformReduceBinary_View)
{

	testTransformReduceBinary_Vw(3);
	testTransformReduceBinary_Vw(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransformReduceBinary_Vw(SZ);
	}

	testTransformReduceBinary_Vw(34);
	testTransformReduceBinary_Vw(63);
	testTransformReduceBinary_Vw(64);
	testTransformReduceBinary_Vw(65);
}




void testTransformReduceBinary_Vec(int SZ)
{


	std::vector<Numeric> input(SZ, asNumber( 0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(begin(input), end(input), g);
	VecXX testVec(input);

	auto testVecPlusTwo = testVec + asNumber(2.0);

	auto SUM = [](auto x, auto y) { return x + y; };
	auto MULT = [](auto x, auto y) { return x * y; };

	auto CPY = [](auto x) { return x; };


	//auto testVw = transformV(CPY, testVec);
	//auto testPlusTwoVw = transformV(CPY, testVecPlusTwo);

	auto inputPlusTwo = input;
	for (auto& x : inputPlusTwo)
	{
		x += asNumber( 2.0);
	}

	auto resSUM = transformReduce(testVec, testVecPlusTwo, MULT, SUM);
	auto ipRes = std::inner_product(begin(input), end(input), begin(inputPlusTwo), asNumber(.0) );

	EXPECT_NUMERIC_EQ(resSUM, ipRes);

}

TEST(TestDR3, testTransformReduceBinary_Vec)
{

	testTransformReduceBinary_Vec(3);
	testTransformReduceBinary_Vec(4);


	for (int SZ = 3; SZ < 33; SZ++)
	{
		testTransformReduceBinary_Vec(SZ);
	}

	testTransformReduceBinary_Vec(34);
	testTransformReduceBinary_Vec(63);
	testTransformReduceBinary_Vec(64);
	testTransformReduceBinary_Vec(65);
}



//broken on AVX512 

///*

void testFilterBoolVecXX(int SZ)
{
	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	const VecXX testVec(input);

	auto copyLambda = [&](auto x) {return x; };
	

	auto allowAll = [](auto x) {return x == x; };

	VecBL vBoolAllTrue = testVec == testVec;
	VecBL vBoolAllFalse = testVec != testVec;


	 auto resAll = filterB(vBoolAllTrue, testVec);


	for (int i = 0; i < resAll.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(resAll[i], testVec[i]);
	}


	auto allowNone = [](auto x) {return x != x; };

	auto resNone = filterB(vBoolAllFalse, testVec);

	EXPECT_EQ(0, resNone.size());



	for (Numeric j = 0; j < SZ; ++j)
	{
		auto testVecB = vBoolAllFalse;
		testVecB.setAt(j, true);
		auto res = filterB(testVecB, testVec);
		EXPECT_NUMERIC_EQ(res[0], j);

	}
		 /**/

}




TEST(TestDR3, testFilterBool_vec)
{

	testFilterBoolVecXX(3);
	testFilterBoolVecXX(4);

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testFilterBoolVecXX(SZ);
	}

	testFilterBoolVecXX(34);
	testFilterBoolVecXX(63);
	testFilterBoolVecXX(64);
	testFilterBoolVecXX(65);

}



///*

void testFilterBool_vw(int SZ)
{
	std::vector<Numeric> input(SZ, asNumber(0.0));
	std::iota(begin(input), end(input), asNumber(0.0));

	const VecXX testVec(input);

	auto copyLambda = [&](auto x) {return x; };
	VecVW inputView = transformV(copyLambda, testVec);

	auto allowAll = [](auto x) {return x == x; };

	VecBL vBoolAllTrue = testVec == testVec;
	VecBL vBoolAllFalse = testVec != testVec;


	auto resAll = filterB(vBoolAllTrue, inputView);


	for (int i = 0; i < resAll.size(); ++i)
	{
		EXPECT_NUMERIC_EQ(resAll[i], testVec[i]);
	}

	auto allowNone = [](auto x) {return x != x; };

	auto resNone = filterB(vBoolAllFalse, inputView);

	EXPECT_EQ(0, resNone.size());



	for (Numeric j = 0; j < SZ; ++j)
	{
		auto testVecB = vBoolAllFalse;
		testVecB.setAt(j, true);
		auto res = filterB(testVecB, inputView);
		EXPECT_NUMERIC_EQ(res[0], j);

	}

}




TEST(TestDR3, testFilterBool_vw)
{

	testFilterBool_vw(3);
	testFilterBool_vw(4);

	for (int SZ = 3; SZ < 33; SZ++)
	{
		testFilterBool_vw(SZ);
	}

	testFilterBool_vw(34);
	testFilterBool_vw(63);
	testFilterBool_vw(64);
	testFilterBool_vw(65);

}


//*/


