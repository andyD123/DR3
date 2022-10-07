#include "pch.h"


#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/vec_bool_d.h"
#include "../Vectorisation/VecX/vec_double.h"
#include  "../Vectorisation/VecX/alloc_policy.h"
#include  "../Vectorisation/VecX/accumulate_transform.h"


#include "../Vectorisation/VecX/target_name_space.h"



#include <numeric>

using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
//using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;


void evalAccumulate(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 6.66);
		for (int i = 0; i < SZ; i++) { v[i] -= 500.0 - i; }
		VecXX test(v);

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };


		double resSTL = std::accumulate(v.begin(), v.end(), 0.0);
		double resAcc = ApplyAccumulate2(test, Sum, Sum, 0.0);

		double resAcc2 = ApplyAccumulate(test, Sum,  0.0);

		double resAcc3 = ApplyAccumulate2(test, Sum, 0.0);

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





void evalAccumulateUR(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000000004;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 6.66);
		for (int i = 0; i < SZ; i++) { v[i] -= 500.0 - i; }
		VecXX test(v);

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };


		double resSTL = std::accumulate(v.begin(), v.end(), 0.0);
		//UR unrolled
		//double resAcc = ApplyAccumulate2UR(test, Sum, Sum, 0.0);
		double resAcc = ApplyAccumulate(test, Sum,  0.0);

		double resAcc1 =  ApplyAccumulate2UR(test, Sum,  0.0);

		EXPECT_NEAR(resSTL, resAcc, testEpsilon);
		EXPECT_NEAR(resSTL, resAcc1, testEpsilon);
	}

}


void evalAccumulateUR_X(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000000004;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 6.66);
		for (int i = 0; i < SZ; i++) { v[i] -= 500.0 - i; }
		VecXX test(v);

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };


		double resSTL = std::accumulate(v.begin(), v.end(), 0.0);
		//UR unrolled
		//double resAcc = ApplyAccumulate2UR(test, Sum, Sum, 0.0);
		double resAcc = ApplyAccumulate(test, Sum, 0.0);

		double resAcc1 = ApplyAccumulate2UR_X(test, Sum);// , 0.0);

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




void evalTransformAccumulate(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 6.66);
		for (int i = 0; i < SZ; i++) { v[i] -= 500.0 - i; }
		VecXX test(v);

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };

		auto SQR = [](auto lhs ) { return lhs *lhs; };


		//double resSTL = std::accumulate(v.begin(), v.end(), 0.0);
		double S = 0.;
		for (double x : v)
		{
			S += x * x;
		}
		
		double resAcc = ApplyTransformAccumulate(test,  SQR, SQR, Sum ,Sum, 0.0);

		double resAcc2 = ApplyTransformAccumulate(test, SQR, Sum,  0.0);

		double resAcc3 = ApplyTransformAccumulateUR(test, SQR, Sum, 0.0);


		EXPECT_NEAR(S, resAcc, testEpsilon);
		EXPECT_NEAR(S, resAcc2, testEpsilon);
		EXPECT_NEAR(S, resAcc3, testEpsilon);
	}

}


TEST(TestAccumulator, simpleTransformAccumulate) {
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	//eval over very small lengths  OK
	evalTransformAccumulate(3, 23);

	//eval over multiple lengths
	evalTransformAccumulate(957, 1043);



}




void evalTransformAccumulate2UR_X(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 6.66);
		for (int i = 0; i < SZ; i++) { v[i] -= 500.0 - i; }
		VecXX test(v);

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };

		auto SQR = [](auto lhs) { return lhs * lhs; };


		//double resSTL = std::accumulate(v.begin(), v.end(), 0.0);
		double S = 0.;
		for (double x : v)
		{
			S += x * x;
		}

		double resAcc = ApplyTransformAccumulate(test, SQR, SQR, Sum, Sum, 0.0);

		double resAcc2 = ApplyTransformAccumulate(test, SQR, Sum, 0.0);

		double resAcc3 = ApplyTransformAccumulate2UR_X(test, SQR, Sum);// , 0.0);


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







void evalTransformAccumulateMax(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 6.66);
		for (int i = 0; i < SZ; i++) { v[i] -= 500.0 - i; }
		VecXX test(v);

		auto MAX = [](auto lhs, auto rhs) {  return max(VecXX::INS( lhs), VecXX::INS( rhs)); };
		auto MAX_dbl = [](double lhs, double rhs) {  return std::max(lhs,rhs); };

		auto SQR = [](auto lhs) { return lhs * lhs; };


		//double resSTL = std::accumulate(v.begin(), v.end(), 0.0);
		double S = v[0] * v[0];
		for (double x : v)
		{
			S = std::max(S, x * x);
		}

		double resAcc = ApplyTransformAccumulate(test, SQR, SQR, MAX, MAX_dbl, test[0]);

		double resAcc2 = ApplyTransformAccumulate(test, SQR, MAX, test[0]);

		double resAcc3 = ApplyTransformAccumulateUR(test, SQR, MAX, test[0]);


		EXPECT_NEAR(S, resAcc, testEpsilon);
		EXPECT_NEAR(S, resAcc2, testEpsilon);
		EXPECT_NEAR(S, resAcc3, testEpsilon);
	}

}



void evalTransformAccumulateMax_X(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 6.66);
		for (int i = 0; i < SZ; i++) { v[i] -= 500.0 - i; }
		VecXX test(v);

		auto MAX = [](auto lhs, auto rhs) {  return max(VecXX::INS(lhs), VecXX::INS(rhs)); };
		auto MAX_dbl = [](double lhs, double rhs) {  return std::max(lhs, rhs); };

		auto SQR = [](auto lhs) { return lhs * lhs; };


		//double resSTL = std::accumulate(v.begin(), v.end(), 0.0);
		double S = v[0] * v[0];
		for (double x : v)
		{
			S = std::max(S, x * x);
		}

		double resAcc = ApplyTransformAccumulate(test, SQR, SQR, MAX, MAX_dbl, test[0]);

		double resAcc2 = ApplyTransformAccumulate(test, SQR, MAX, test[0]);

		double resAcc3 = ApplyTransformAccumulate2UR_X(test, SQR, MAX);// , test[0]);


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



void evalTransformAccumulateInnerProduct_X(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 6.66);
		auto w = v;
		for (int i = 0; i < SZ; i++) { v[i] -= 500.0 - i; w[i] = v[i] + 1.0; }
		VecXX test(v);
		VecXX test2(w);


		auto MULT = [](auto lhs, auto rhs) {  return lhs * rhs; };
		auto SUM = [](auto lhs, auto rhs) {  return lhs + rhs; };

		auto S = std::inner_product(v.begin(), v.end(), w.begin(), 0.0);
		double resAcc3 = ApplyTransformAccumulate2UR_XBin(test, test2, MULT, SUM);

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