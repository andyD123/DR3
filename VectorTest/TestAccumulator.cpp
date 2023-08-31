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
#include <algorithm>


void setEpsilonTransformAccumulate(long double& val)
{
	const long double testEpsilon = 0.000001;
	val = testEpsilon;
}


void setEpsilonTransformAccumulate(double& val)
{
	const double testEpsilon = 0.000001;
	val= testEpsilon;
}

void setEpsilonTransformAccumulate(float& val)
{
	const float testEpsilon = 100.f;
	val = testEpsilon;
}

void setEpsilonAccumulate(long double& val)
{
	const long double testEpsilon = 0.000000004;
	val = testEpsilon;
}


void setEpsilonAccumulate(double& val)
{
	const double testEpsilon = 0.000000004;
	val = testEpsilon;
}

void setEpsilonAccumulate(float& val)
{
	const float testEpsilon = 0.4f;
	val = testEpsilon;
}


//const Numeric testEpsilon = 0.000001;
const Numeric testEpsilon = 100.;



auto  getVec(int SZ, std::vector<Numeric>& stl)
{

	std::vector<Numeric>  v(SZ, asNumber(6.66));
	int i = 0;
	for (auto& x : v) { x -= asNumber(500.0 + i); ++i; }
	VecXX test(v);
	stl = v;
	return  test;

}


void evalAccumulate(int startLen, int endLen)
{
	
	Numeric testEpsilon;
	setEpsilonTransformAccumulate(testEpsilon);


	for (int  SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<Numeric> v;
		VecXX test = getVec(SZ,v);
		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };

		Numeric resSTL = std::accumulate(v.begin(), v.end(), asNumber(0.0));
		Numeric resAcc = ApplyAccumulate2(test, Sum, Sum, asNumber(0.0));
		Numeric resAcc2 = ApplyAccumulate(test, Sum, asNumber(0.0));
		Numeric resAcc3 = ApplyAccumulate2(test, Sum, asNumber(0.0));

		EXPECT_NEAR(resSTL, resAcc,  testEpsilon);
		EXPECT_NEAR(resSTL, resAcc2, testEpsilon);
		EXPECT_NEAR(resSTL, resAcc3, testEpsilon);
	}

}

TEST(TestAccumulator, simpleSummation) {
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	//eval over multiple3 lengths
	evalAccumulate(957, 1043);

	//eval over very small lengths
	evalAccumulate(3, 23);

}





void evalAccumulateUR(int startLen, int endLen)
{

	Numeric testEpsilon;
	setEpsilonAccumulate(testEpsilon);

	for (int SZ = startLen; SZ <= endLen; SZ++)
	{
		
		std::vector<Numeric>  v(SZ, asNumber(0.0));
		for (int i = 0; i < SZ; i++) { v[i] -= asNumber(500.0 - i); }
		VecXX test(v);
	
		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
		Numeric resSTL = std::accumulate(v.begin(), v.end(), asNumber(0.0));
		Numeric resAcc = ApplyAccumulate(test, Sum, asNumber(0.0));
		Numeric resAcc1 =  ApplyAccumulate2UR(test, Sum, asNumber(0.0));

		EXPECT_NEAR(resSTL, resAcc, testEpsilon);
		EXPECT_NEAR(resSTL, resAcc1, testEpsilon);
	}

}


void evalAccumulateUR_X(int startLen, int endLen)
{
	Numeric testEpsilon;
	setEpsilonAccumulate(testEpsilon);


	for (int SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(0.0));
		for (int i = 0; i < SZ; i++) { v[i] -= asNumber(500.0- i); }
		VecXX test(v);

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
		Numeric resSTL = std::accumulate(v.begin(), v.end(), asNumber(0.0));
		Numeric resAcc = ApplyAccumulate(test, Sum, asNumber(0.0));
		Numeric resAcc1 = ApplyAccumulate2UR_X(test, Sum);

		EXPECT_NEAR(resSTL, resAcc, testEpsilon);
		EXPECT_NEAR(resSTL, resAcc1, testEpsilon);
	}

}

TEST(TestAccumulator, simpleSummationUnRolled)
{
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	//eval over multiple3 lengths
	evalAccumulateUR(957, 1043);

	//eval over very small lengths  OK
	evalAccumulateUR(3, 23);

}

//extended function hand unrolled
TEST(TestAccumulator, simpleSummationUnRolled_X)
{
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	//eval over multiple3 lengths
	evalAccumulateUR_X(957, 1043);

	//eval over very small lengths  OK
	evalAccumulateUR_X(3, 23);

}




void evalTransformAccumulate(int startLen, int endLen)
{
	Numeric testEpsilon;
	setEpsilonTransformAccumulate(testEpsilon);
	//testEpsilon *= 5.;

	for (int SZ = startLen; SZ <= endLen; SZ++)
	{
		
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		for (int i = 0; i < SZ; i++) { v[i] -= asNumber(500.0) - i; }
		VecXX test(v);

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
		auto SQR = [](auto lhs ) { return lhs *lhs; };

		Numeric S = 0.;
		for (Numeric x : v)
		{
			S += x * x;
		}
		
		Numeric resAcc = ApplyTransformAccumulate(test,  SQR, SQR, Sum ,Sum, asNumber( 0.0));
		Numeric resAcc2 = ApplyTransformAccumulate(test, SQR, Sum, asNumber(0.0));
		Numeric resAcc3 = ApplyTransformAccumulateUR(test, SQR, Sum, asNumber(0.0));

		EXPECT_NEAR(S, resAcc, testEpsilon);
		EXPECT_NEAR(S, resAcc2, testEpsilon);
		EXPECT_NEAR(S, resAcc3, testEpsilon);
	}

}


TEST(TestAccumulator, simpleTransformAccumulate)
{
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	//eval over very small lengths  OK
	evalTransformAccumulate(3, 23);

	//eval over multiple lengths
	evalTransformAccumulate(957, 1043);
}




void evalTransformAccumulate2UR_X(int startLen, int endLen)
{
	Numeric testEpsilon;
	setEpsilonTransformAccumulate(testEpsilon);

	for (int SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		for (int i = 0; i < SZ; i++) { v[i] -= asNumber(500.0) - i; }
		VecXX test(v);

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
		auto SQR = [](auto lhs) { return lhs * lhs; };

		Numeric S = 0.;
		for (Numeric x : v)
		{
			S += x * x;
		}

		Numeric resAcc = ApplyTransformAccumulate(test, SQR, SQR, Sum, Sum, asNumber(0.0));
		Numeric resAcc2 = ApplyTransformAccumulate(test, SQR, Sum, asNumber(0.0));
		Numeric resAcc3 = ApplyTransformAccumulate2UR_X(test, SQR, Sum);

		EXPECT_NEAR(S, resAcc, testEpsilon);
		EXPECT_NEAR(S, resAcc2, testEpsilon);
		EXPECT_NEAR(S, resAcc3, testEpsilon);
	}

}




TEST(TestAccumulator, simpleTransformAccumulate_X) {
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	//eval over very small lengths  OK
    evalTransformAccumulate2UR_X(3, 23); //OK

	//eval over multiple lengths
	evalTransformAccumulate2UR_X(957, 1043);

}







void evalTransformAccumulateMax(int startLen, int endLen)
{
	Numeric testEpsilon;
	setEpsilonTransformAccumulate(testEpsilon);

	for (int SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<Numeric>  v;
		VecXX test = getVec(SZ, v);

		auto MAX = [](auto lhs, auto rhs) {  return max(VecXX::INS( lhs), VecXX::INS( rhs)); };
		auto MAX_dbl = [](Numeric lhs, Numeric rhs) {  return std::max(lhs,rhs); };
		auto SQR = [](auto lhs) { return lhs * lhs; };

		Numeric S = v[0] * v[0];
		for (Numeric x : v)
		{
			S = std::max(S, x * x);
		}

		Numeric resAcc = ApplyTransformAccumulate(test, SQR, SQR, MAX, MAX_dbl, test[0]);
		Numeric resAcc2 = ApplyTransformAccumulate(test, SQR, MAX, test[0]);
		Numeric resAcc3 = ApplyTransformAccumulateUR(test, SQR, MAX, test[0]);

		EXPECT_NEAR(S, resAcc, testEpsilon);
		EXPECT_NEAR(S, resAcc2, testEpsilon);
		EXPECT_NEAR(S, resAcc3, testEpsilon);
	}

}



void evalTransformAccumulateMax_X(int startLen, int endLen)
{
	Numeric testEpsilon;
	setEpsilonTransformAccumulate(testEpsilon);

	for (int SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<Numeric>  v;
		VecXX test = getVec(SZ, v);

		auto MAX = [](auto lhs, auto rhs) {  return max(VecXX::INS(lhs), VecXX::INS(rhs)); };
		auto MAX_dbl = [](Numeric lhs, Numeric rhs) {  return std::max(lhs, rhs); };
		auto SQR = [](auto lhs) { return lhs * lhs; };

		Numeric S = v[0] * v[0];
		for (Numeric x : v)
		{
			S = std::max(S, x * x);
		}

		Numeric resAcc = ApplyTransformAccumulate(test, SQR, SQR, MAX, MAX_dbl, test[0]);
		Numeric resAcc2 = ApplyTransformAccumulate(test, SQR, MAX, test[0]);
		Numeric resAcc3 = ApplyTransformAccumulate2UR_X(test, SQR, MAX);

		EXPECT_NEAR(S, resAcc, testEpsilon);
		EXPECT_NEAR(S, resAcc2, testEpsilon);
		EXPECT_NEAR(S, resAcc3, testEpsilon);
	}

}



TEST(TestAccumulator, simpleTransformAccumulate_Max) {
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	//eval over very small lengths  OK
	evalTransformAccumulate(3, 23);

	//eval over multiple lengths
	 evalTransformAccumulateMax(957, 1043);

}


TEST(TestAccumulator, simpleTransformAccumulate_Max_X) {
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	//eval over very small lengths  OK
	evalTransformAccumulate2UR_X(3, 23);

	//eval over multiple lengths
	evalTransformAccumulateMax_X(957, 1043);

}



void evalTransformAccumulateInnerProduct_X(int startLen, int endLen)
{

	Numeric testEpsilon;
	setEpsilonTransformAccumulate(testEpsilon);

	for (int SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		auto w = v;
		for (int i = 0; i < SZ; i++) { v[i] -= asNumber(500.0) - i; w[i] = v[i] + asNumber(1.0); }
		VecXX test(v);
		VecXX test2(w);
		auto MULT = [](auto lhs, auto rhs) {  return lhs * rhs; };
		auto SUM = [](auto lhs, auto rhs) {  return lhs + rhs; };
		auto S = std::inner_product(v.begin(), v.end(), w.begin(), asNumber(0.0));
		Numeric resAcc3 = ApplyTransformAccumulate2UR_XBin(test, test2, MULT, SUM);

		EXPECT_NEAR(S, resAcc3, testEpsilon);
	}

}



TEST(TestAccumulator, simpleTransformAccumulate_BinaryFunc_X)
{
	//eval over very small lengths  OK
	evalTransformAccumulateInnerProduct_X(3, 23);

	//eval over multiple lengths
	evalTransformAccumulateInnerProduct_X(957, 1043);

}


////////////////////////////////////////


void evalAccumulateUR_X_pairwise(int startLen, int endLen)
{
	Numeric testEpsilon;
	setEpsilonAccumulate(testEpsilon);


	for (int SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(0.0));
		for (int i = 0; i < SZ; i++)
		{ 
			v[i] -= asNumber(500.0 - i);
		}

		VecXX test(v);
		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
		Numeric resSTL = std::accumulate(v.begin(), v.end(), asNumber(0.0));
		Numeric resAcc1 = ApplyAccumulate2UR_X_pairwise(test, Sum);
		EXPECT_NEAR(resSTL, resAcc1, testEpsilon);
	}

}

TEST(TestAccumulator, simple_pairwise)
{
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	evalAccumulateUR_X_pairwise(32, 32);
	evalAccumulateUR_X_pairwise(64, 64);
	evalAccumulateUR_X_pairwise(128, 128);
	evalAccumulateUR_X_pairwise(256, 256);
	evalAccumulateUR_X_pairwise(512, 512);
	evalAccumulateUR_X_pairwise(1024, 1024); 

	evalAccumulateUR_X_pairwise(1, 80);

	evalAccumulateUR_X_pairwise(957, 1043);


	evalAccumulateUR_X_pairwise(159, 159);

}



////////////////////////////////////////////////////

void evalTransformAccumulateUR_X_pairwise(int startLen, int endLen)
{
	Numeric testEpsilon;
	setEpsilonAccumulate(testEpsilon);


	for (int SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(0.0));
		for (int i = 0; i < SZ; i++)
		{
			v[i] -= asNumber(500.0 - i);
		}

		VecXX test(v);
		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
		auto SQR = [](auto lhs) { return lhs * lhs; };
		Numeric resSTL = std::inner_product(v.begin(), v.end(), v.begin(),  0.0);
		Numeric resAcc1 = ApplyTransformAccumulate2UR_X_pairwise(test, SQR,Sum);
		EXPECT_NEAR(resSTL, resAcc1, testEpsilon);
	}

}

TEST(TestAccumulator, simple_pairwise_transformReduce)
{
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	evalTransformAccumulateUR_X_pairwise(32, 32);
	evalTransformAccumulateUR_X_pairwise(64, 64);
	evalTransformAccumulateUR_X_pairwise(128, 128);
	evalTransformAccumulateUR_X_pairwise(256, 256);
	evalTransformAccumulateUR_X_pairwise(512, 512);
	evalTransformAccumulateUR_X_pairwise(1024, 1024);

	evalTransformAccumulateUR_X_pairwise(1, 31);

	evalTransformAccumulateUR_X_pairwise(32, 80);

	evalTransformAccumulateUR_X_pairwise(957, 1043);


	evalTransformAccumulateUR_X_pairwise(159, 159);

}



void evalTransformAccumulateUR_X_bivariate_pairwise(int startLen, int endLen)
{
	Numeric testEpsilon;
	setEpsilonAccumulate(testEpsilon);

	 //rounding errors in stl 


	for (int SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(0.0));

		auto v2 = v;
		for (int i = 0; i < SZ; i++)
		{
			v[i] -= asNumber(500.0 - i) / 3.0;
			v2[i] = 2 * v[i];
		}

		VecXX test(v);
		VecXX test2 = asNumber(2.0) * test;

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
		auto mult = [](auto lhs,auto rhs) { return lhs * rhs; };
		Numeric resSTL = std::inner_product(v.begin(), v.end(), v2.begin(), 0.0);
		Numeric resAcc1 = ApplyTransformAccumulate2UR_X_pairwise(test, test2, mult, Sum);

		testEpsilon *= (abs(resSTL) + abs(resAcc1));
		EXPECT_NEAR(resSTL, resAcc1, testEpsilon);
	}

}

TEST(TestAccumulator, bivariate_pairwise_transformReduce)
{
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	evalTransformAccumulateUR_X_bivariate_pairwise(32, 32);
	evalTransformAccumulateUR_X_bivariate_pairwise(64, 64);
	evalTransformAccumulateUR_X_bivariate_pairwise(128, 128);
	evalTransformAccumulateUR_X_bivariate_pairwise(256, 256);
	evalTransformAccumulateUR_X_bivariate_pairwise(512, 512);
	evalTransformAccumulateUR_X_bivariate_pairwise(1024, 1024);

	evalTransformAccumulateUR_X_bivariate_pairwise(1, 31);
	evalTransformAccumulateUR_X_bivariate_pairwise(32, 80);

	evalTransformAccumulateUR_X_bivariate_pairwise(957, 1043);
	evalTransformAccumulateUR_X_bivariate_pairwise(159, 159);

}


