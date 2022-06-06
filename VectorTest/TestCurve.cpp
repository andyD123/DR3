#include "pch.h"


#include "../../Vectorisation/ExampleVectors/curve.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/vec_bool_d.h"
#include "../Vectorisation/VecX/vec_double.h"
#include  "../Vectorisation/VecX/alloc_policy.h"

typedef VecD<VecDouble>  VecxD;
typedef VecD<VecDouble>  Vecx;
typedef Vec<VecDouble>  VecXX;

TEST(TestCaseCurve, Test1) {
	EXPECT_EQ(1, 1);
	EXPECT_TRUE(true);

	std::vector<double>  values{ 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 };
	std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10 };
	std::vector <double>   datesD = { 0,1,2,3,4,5,6,7,8,9,10 };




	Curve< double, double> testCurve;
	testCurve.setValues(begin(values), end(values), begin(datesD), end(datesD)); //wrong way round

	auto val = testCurve.valueAt(0.0);

	EXPECT_EQ(val, 0.0);
	val = testCurve.valueAt(0.5);
	EXPECT_EQ(val, 0.5);


	///////////////////////////
	std::vector< VecXX>  vecVals;
	for (int i = 0; i < 11; i++)
	{
		VecXX vals(i * 0.001 + 0.06, 100);
		vecVals.push_back(vals);

	}


	{
		using ZeroCrv = Curve< double, VecXX, ZeroInterp<double, VecXX> >;

		ZeroCrv testCurve2;
		testCurve2.setValues(begin(datesD), end(datesD), begin(vecVals), end(vecVals));

		auto valV = testCurve2.valueAt(0.0);

		auto valV2 = testCurve2.valueAt(0.5);
	}


	{
		std::vector<double>  values{ 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 };
		std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10 };
		std::vector <double>   datesD = { 0,1,2,3,4,5,6,7,8,9,10 };

		std::vector< VecXX>  vecVals;
		for (int i = 0; i < 11; i++)
		{
			VecXX vals(i * 0.001 + 0.06, 100);
			vecVals.push_back(vals);

		}


		Curve2<double, VecXX, ZeroInterpCached<double, VecXX>>  testCurve2(10);
		testCurve2.setValues(begin(datesD), end(datesD), begin(vecVals), end(vecVals));

		auto valV = testCurve2.valueAt(0.0);

		auto valV2 = testCurve2.valueAt(0.5);


		for (long l = 0; l < 10000; l++)
		{
			auto valV3 = testCurve2.valueAt(0.5);
		}
	}

	//EXPECT_EQ(val, 0.0);
	//val = testCurve.valueAt(0.5);
	//EXPECT_EQ(val, 0.5);
}