#include "pch.h"


#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/vec_bool_d.h"
#include "../Vectorisation/VecX/vec_double.h"
#include  "../Vectorisation/VecX/alloc_policy.h"
#include "../Vectorisation/VecX/vec_d.h"
#include "../Vectorisation/VecX/vec_bool.h"
#include "../Vectorisation/VecX/apply_operation.h"



#include "../Vectorisation/VecX/target_name_space.h"
#include "../Vectorisation/VecX/dr3.h"
#include "testNamespace.h"
#include "dr3TestUtil.h"


#include <numeric>



TEST(TestCaseName, TestMakeVec)
{

	std::vector<Numeric> three(3, asNumber(42.0));
	//std::vector<Numeric> mix{ 1.0,2.0,3.0 };
	std::vector<Numeric> mix{ asNumber(1.0), asNumber(2.0),asNumber(3.0) };
	VecXX Vec2(three);
	VecXX Vec1(mix);

	EXPECT_EQ(Vec2.size(), 3);
	EXPECT_EQ(Vec1[0], 1.0);
	EXPECT_EQ(Vec1[1], 2.0);
	EXPECT_EQ(Vec1[2], 3.0);

	EXPECT_NUMERIC_EQ(Vec1[0], asNumber(1.0));
	EXPECT_NUMERIC_EQ(Vec1[1], asNumber(2.0));
	EXPECT_NUMERIC_EQ(Vec1[2], asNumber(3.0));

}

TEST(TestCaseName, TestAddVec)
{

	std::vector<Numeric> three(3, asNumber(42.0));
	//std::vector<Numeric> mix{ 1.0,2.0,3.0 };
	std::vector<Numeric> mix{ asNumber(1.0), asNumber(2.0),asNumber(3.0) };
	VecXX Vec2(three);
	VecXX Vec1(mix);

	auto added = Vec2 + Vec1;
	EXPECT_EQ(added.size(), 3);
	//EXPECT_EQ(added[0], 43.0);
	//EXPECT_EQ(added[1], 44.0);
	//EXPECT_EQ(added[2], 45.0);

	EXPECT_NUMERIC_EQ(added[0], asNumber(43.0));
	EXPECT_NUMERIC_EQ(added[1], asNumber(44.0));
	EXPECT_NUMERIC_EQ(added[2], asNumber(45.0));

	auto added1 = Vec1 + Vec2;
	EXPECT_EQ(added1.size(), 3);
	//EXPECT_EQ(added1[0], 43.0);
	//EXPECT_EQ(added1[1], 44.0);
	//EXPECT_EQ(added1[2], 45.0);

	EXPECT_NUMERIC_EQ(added1[0], asNumber(43.0));
	EXPECT_NUMERIC_EQ(added1[1], asNumber(44.0));
	EXPECT_NUMERIC_EQ(added1[2], asNumber(45.0));

	auto added2 = Vec1;
	added2 += Vec2;
	EXPECT_EQ(added2.size(), 3);
	//EXPECT_EQ(added2[0], 43.0);
	//EXPECT_EQ(added2[1], 44.0);
	//EXPECT_EQ(added2[2], 45.0);

	EXPECT_NUMERIC_EQ(added2[0], asNumber(43.0));
	EXPECT_NUMERIC_EQ(added2[1], asNumber(44.0));
	EXPECT_NUMERIC_EQ(added2[2], asNumber(45.0));

	Numeric val = 42.0;
	auto added3 = Vec1 + val;
	EXPECT_EQ(added3.size(), 3);
	//EXPECT_EQ(added3[0], 43.0);
	//EXPECT_EQ(added3[1], 44.0);
	//EXPECT_EQ(added3[2], 45.0);

	EXPECT_NUMERIC_EQ(added3[0], asNumber(43.0));
	EXPECT_NUMERIC_EQ(added3[1], asNumber(44.0));
	EXPECT_NUMERIC_EQ(added3[2], asNumber(45.0));


	auto added4 = val + Vec1;
	EXPECT_EQ(added4.size(), 3);
	//EXPECT_EQ(added4[0], 43.0);
	//EXPECT_EQ(added4[1], 44.0);
	//EXPECT_EQ(added4[2], 45.0);

	EXPECT_NUMERIC_EQ(added4[0], asNumber(43.0));
	EXPECT_NUMERIC_EQ(added4[1], asNumber(44.0));
	EXPECT_NUMERIC_EQ(added4[2], asNumber(45.0));



	//VecXX added5( val);
	//added5 +=Vec1;
	Vec1 += val;
	EXPECT_EQ(Vec1.size(), 3);
	//EXPECT_EQ(Vec1[0], 43.0);
	//EXPECT_EQ(Vec1[1], 44.0);
	//EXPECT_EQ(Vec1[2], 45.0);

	EXPECT_NUMERIC_EQ(Vec1[0], asNumber(43.0));
	EXPECT_NUMERIC_EQ(Vec1[1], asNumber(44.0));
	EXPECT_NUMERIC_EQ(Vec1[2], asNumber(45.0));

}


TEST(TestCaseName, TestSubVec)
{

	std::vector<Numeric> three(3, asNumber(42.0));
	///std::vector<Numeric> mix{ 1.0,2.0,3.0 };

	std::vector<Numeric> mix{ asNumber(1.0), asNumber(2.0),asNumber(3.0) };
	VecXX Vec2(three);
	VecXX Vec1(mix);

	auto result = Vec2 - Vec1;
	EXPECT_EQ(result.size(), 3);

	//EXPECT_EQ(result[0], 41.0);
	//EXPECT_EQ(result[1], 40.0);
	//EXPECT_EQ(result[2], 39.0);

	EXPECT_NUMERIC_EQ(result[0], asNumber(41.0));
	EXPECT_NUMERIC_EQ(result[1], asNumber(40.0));
	EXPECT_NUMERIC_EQ(result[2], asNumber(39.0));

	auto result1 = Vec1 - Vec2;
	EXPECT_EQ(result1.size(), 3);
	//EXPECT_EQ(result1[0], -41.0);
	//EXPECT_EQ(result1[1], -40.0);
	//EXPECT_EQ(result1[2], -39.0);

	EXPECT_NUMERIC_EQ(result1[0], asNumber(-41.0));
	EXPECT_NUMERIC_EQ(result1[1], asNumber(-40.0));
	EXPECT_NUMERIC_EQ(result1[2], asNumber(-39.0));


	auto result2 = Vec1;
	result2 -= Vec2;
	EXPECT_EQ(result2.size(), 3);
	//EXPECT_EQ(result2[0], -41.0);
	//EXPECT_EQ(result2[1], -40.0);
	//EXPECT_EQ(result2[2], -39.0);

	EXPECT_NUMERIC_EQ(result2[0], asNumber(-41.0));
	EXPECT_NUMERIC_EQ(result2[1], asNumber(-40.0));
	EXPECT_NUMERIC_EQ(result2[2], asNumber(-39.0));


	Numeric val = asNumber(42.0);
	auto result3 = Vec1 - val;
	EXPECT_EQ(result3.size(), 3);
	/*
	EXPECT_EQ(result3[0], -41.0);
	EXPECT_EQ(result3[1], -40.0);
	EXPECT_EQ(result3[2], -39.0);
	*/

	EXPECT_NUMERIC_EQ(result3[0], asNumber(-41.0));
	EXPECT_NUMERIC_EQ(result3[1], asNumber(-40.0));
	EXPECT_NUMERIC_EQ(result3[2], asNumber(-39.0));


	auto result4 = val - Vec1;
	EXPECT_EQ(result4.size(), 3);
/*
	EXPECT_EQ(result4[0], 41.0);
	EXPECT_EQ(result4[1], 40.0);
	EXPECT_EQ(result4[2], 39.0);
*/
	EXPECT_NUMERIC_EQ(result4[0], asNumber(41.0));
	EXPECT_NUMERIC_EQ(result4[1], asNumber(40.0));
	EXPECT_NUMERIC_EQ(result4[2], asNumber(39.0));

	//VecXX result5(val);
	Vec1 -= val;
	EXPECT_EQ(Vec1.size(), 3);
	/*
	EXPECT_EQ(Vec1[0], -41.0);
	EXPECT_EQ(Vec1[1], -40.0);
	EXPECT_EQ(Vec1[2], -39.0);
	*/
	EXPECT_NUMERIC_EQ(Vec1[0], asNumber(-41.0));
	EXPECT_NUMERIC_EQ(Vec1[1], asNumber(-40.0));
	EXPECT_NUMERIC_EQ(Vec1[2], asNumber(-39.0));

}


TEST(TestCaseName, TestMultVec)
{

	std::vector<Numeric> three(3, asNumber(42.0));
	//std::vector<Numeric> mix{ 1.0,2.0,3.0 };
	std::vector<Numeric> mix{ asNumber(1.0), asNumber(2.0),asNumber(3.0) };
	VecXX Vec2(three);
	VecXX Vec1(mix);

	auto result = Vec2 * Vec1;
	EXPECT_EQ(result.size(), 3);
	/*
	EXPECT_EQ(added[0], 42.0);
	EXPECT_EQ(added[1], 84.0);
	EXPECT_EQ(added[2], 126.0);
	*/

	EXPECT_NUMERIC_EQ(result[0], asNumber(42.0));
	EXPECT_NUMERIC_EQ(result[1], asNumber(84.0));
	EXPECT_NUMERIC_EQ(result[2], asNumber(126.0));



	auto result1 = Vec1 * Vec2;
	EXPECT_EQ(result1.size(), 3);
	EXPECT_NUMERIC_EQ(result1[0], asNumber(42.0));
	EXPECT_NUMERIC_EQ(result1[1], asNumber(84.0));
	EXPECT_NUMERIC_EQ(result1[2], asNumber(126.0));

	auto result2 = Vec1;
	result2 *= Vec2;

/*
	EXPECT_EQ(result2.size(), 3);
	EXPECT_EQ(added2[0], 42.0);
	EXPECT_EQ(added2[1], 84.0);
	EXPECT_EQ(added2[2], 126.0);
*/
	EXPECT_EQ(result2.size(), 3);
	EXPECT_NUMERIC_EQ(result2[0], asNumber(42.0));
	EXPECT_NUMERIC_EQ(result2[1], asNumber(84.0));
	EXPECT_NUMERIC_EQ(result2[2], asNumber(126.0));

	Numeric val = asNumber(42.0);
	auto result3 = Vec1 * val;
	EXPECT_EQ(result3.size(), 3);
	/*
	EXPECT_EQ(added3[0], 42.0);
	EXPECT_EQ(added3[1], 84.0);
	EXPECT_EQ(added3[2], 126.0);
	*/
	EXPECT_NUMERIC_EQ(result3[0], asNumber(42.0));
	EXPECT_NUMERIC_EQ(result3[1], asNumber(84.0));
	EXPECT_NUMERIC_EQ(result3[2], asNumber(126.0));


	auto result4 = val * Vec1;
	EXPECT_EQ(result4.size(), 3);
	EXPECT_NUMERIC_EQ(result4[0], asNumber(42.0));
	EXPECT_NUMERIC_EQ(result4[1], asNumber(84.0));
	EXPECT_NUMERIC_EQ(result4[2], asNumber(126.0));


	/*
	EXPECT_EQ(added4[0], 42.0);
	EXPECT_EQ(added4[1], 84.0);
	EXPECT_EQ(added4[2], 126.0);
	*/

	//VecXX added5(val);
	Vec1 *= val;
	EXPECT_EQ(Vec1.size(), 3);

	/*
	EXPECT_EQ(Vec1[0], 42.0);
	EXPECT_EQ(Vec1[1], 84.0);
	EXPECT_EQ(Vec1[2], 126.0);
	*/
	EXPECT_NUMERIC_EQ(Vec1[0], asNumber(42.0));
	EXPECT_NUMERIC_EQ(Vec1[1], asNumber(84.0));
	EXPECT_NUMERIC_EQ(Vec1[2], asNumber(126.0));

}


TEST(TestCaseName, TestDivVec)
{

	std::vector<Numeric> three(3, asNumber(4.0));
	std::vector<Numeric> mix{ asNumber(10.0),asNumber(20.),asNumber(30.0) };
	VecXX Vec2(mix);
	VecXX Vec1(three);

	auto added = Vec2 / Vec1;
	EXPECT_EQ(added.size(), 3);
	EXPECT_NUMERIC_EQ(added[0], asNumber(2.5));
	EXPECT_NUMERIC_EQ(added[1], asNumber(5.0));
	EXPECT_NUMERIC_EQ(added[2], asNumber(7.5));

	auto added1 = Vec1 / Vec2;
	EXPECT_EQ(added1.size(), 3);
	EXPECT_NUMERIC_EQ(added1[0], asNumber(0.40));
	EXPECT_NUMERIC_EQ(added1[1], asNumber(0.20));
	EXPECT_NUMERIC_EQ(added1[2], asNumber(4.0 / 30.0));

	auto added2 = Vec1;
	added2 /= Vec2;
	EXPECT_EQ(added2.size(), 3);
	EXPECT_NUMERIC_EQ(added2[0], asNumber(0.40));
	EXPECT_NUMERIC_EQ(added2[1], asNumber(0.20));
	EXPECT_NUMERIC_EQ(added2[2], asNumber(4.0 / 30.0));


	Numeric val = 4.0;
	auto added3 = Vec2 / val;
	EXPECT_EQ(added3.size(), 3);
	EXPECT_NUMERIC_EQ(added3[0], asNumber(2.5));
	EXPECT_NUMERIC_EQ(added3[1], asNumber(5.0));
	EXPECT_NUMERIC_EQ(added3[2], asNumber(7.5));


	auto added4 = val / Vec1;
	EXPECT_EQ(added4.size(), 3);
	EXPECT_NUMERIC_EQ(added4[0], asNumber(1.));
	EXPECT_NUMERIC_EQ(added4[1], asNumber(1.));
	EXPECT_NUMERIC_EQ(added4[2], asNumber(1.));


	Vec2 /= val;
	EXPECT_EQ(Vec2.size(), 3);
	EXPECT_NUMERIC_EQ(Vec2[0], asNumber(2.5));
	EXPECT_NUMERIC_EQ(Vec2[1], asNumber(5.0));
	EXPECT_NUMERIC_EQ(Vec2[2], asNumber(7.5));

}


TEST(TestFuncions, TestExp)
{

	std::vector<Numeric> three(3, asNumber(4.0));
	std::vector<Numeric> mix{ asNumber(1.0),asNumber(2.),asNumber(3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(2.2));

//	std::vector<double> three(3, 4.0);
//	std::vector<double> mix{ 1.0,2.0,3.0 };
//	//VecXX Vec2(mix);
//	VecXX Vec1(mix);
//	VecXX vecScalar(2.2);

	auto res = exp(Vec1);
	EXPECT_EQ(res.size(), 3);

	EXPECT_NUMERIC_EQ(asNumber( res[0] ), asNumber( exp(mix[0]) )  );
	EXPECT_NUMERIC_EQ(asNumber( res[1] ), asNumber( exp(mix[1]) )  );
	EXPECT_NUMERIC_EQ(asNumber (res[2] ), asNumber( exp(mix[2]) )  );

	auto res2 = exp(vecScalar);
	EXPECT_NUMERIC_EQ(res2.getScalarValue(), exp(vecScalar.getScalarValue()));
}


TEST(TestFuncions, TestLog)
{
/*
	std::vector<double> three(3, 4.0);
	std::vector<double> mix{ 1.0,2.0,3.0 };
	//VecXX Vec2(mix);
	VecXX Vec1(mix);
	VecXX vecScalar(2.2);
*/
	std::vector<Numeric> three(3, asNumber(4.0));
	//std::vector<Numeric> mix{ 1.0,2.0,3.0 };
	std::vector<Numeric> mix{ asNumber(1.0),asNumber(2.),asNumber(3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(2.2));

	auto res = log(Vec1);
	EXPECT_EQ(res.size(), 3);
	/*
	EXPECT_DOUBLE_EQ(res[0], log(mix[0]));
	EXPECT_DOUBLE_EQ(res[1], log(mix[1]));
	EXPECT_DOUBLE_EQ(res[2], log(mix[2]));
	*/
	EXPECT_NUMERIC_EQ(asNumber(res[0]), asNumber(log(mix[0])));
	EXPECT_NUMERIC_EQ(asNumber(res[1]), asNumber(log(mix[1])));
	EXPECT_NUMERIC_EQ(asNumber(res[2]), asNumber(log(mix[2])));

	auto res2 = log(vecScalar);
	//EXPECT_DOUBLE_EQ(res2.getScalarValue(), log(vecScalar.getScalarValue()));
	EXPECT_NUMERIC_EQ(res2.getScalarValue(), log(vecScalar.getScalarValue()));
}

TEST(TestFuncions, Testabs)
{
/*
	std::vector<double> three(3, 4.0);
	std::vector<double> mix{ 1.0,2.0,3.0 };
	//VecXX Vec2(mix);
	VecXX Vec1(mix);
	VecXX vecScalar(2.2);
*/
	//std::vector<Numeric> three(3, asNumber(4.0));
	//std::vector<Numeric> mix{ 1.0,2.0,3.0 };
	std::vector<Numeric> mix{ asNumber(1.0),asNumber(2.),asNumber(3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(2.2));


	auto res = abs(Vec1);
	EXPECT_EQ(res.size(), 3);
	/*
	EXPECT_DOUBLE_EQ(res[0], abs(mix[0]));
	EXPECT_DOUBLE_EQ(res[1], abs(mix[1]));
	EXPECT_DOUBLE_EQ(res[2], abs(mix[2]));
	*/

	EXPECT_NUMERIC_EQ(asNumber(res[0]), asNumber(abs(mix[0])));
	EXPECT_NUMERIC_EQ(asNumber(res[1]), asNumber(abs(mix[1])));
	EXPECT_NUMERIC_EQ(asNumber(res[2]), asNumber(abs(mix[2])));

	auto res2 = abs(vecScalar);
	//EXPECT_DOUBLE_EQ(res2.getScalarValue(), abs(vecScalar.getScalarValue()));
	EXPECT_NUMERIC_EQ(res2.getScalarValue(), abs(vecScalar.getScalarValue()));
}


TEST(TestFuncions, Testfloor)
{
	//std::vector<Numeric> three(3, asNumber(4.0));
	//std::vector<Numeric> mix{ 1.0,2.0,3.0 };
	std::vector<Numeric> mix{ asNumber(1.0),asNumber(2.),asNumber(3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(2.2));
/*

	std::vector<double> three(3, 4.0);
	std::vector<double> mix{ 1.0,2.0,3.0 };
	//VecXX Vec2(mix);
	VecXX Vec1(mix);
	VecXX vecScalar(2.2);
*/
	auto res = floor(Vec1);
	EXPECT_EQ(res.size(), 3);
/*
	EXPECT_DOUBLE_EQ(res[0], floor(mix[0]));
	EXPECT_DOUBLE_EQ(res[1], floor(mix[1]));
	EXPECT_DOUBLE_EQ(res[2], floor(mix[2]));
	*/

	EXPECT_NUMERIC_EQ(asNumber(res[0]), asNumber(floor(mix[0])));
	EXPECT_NUMERIC_EQ(asNumber(res[1]), asNumber(floor(mix[1])));
	EXPECT_NUMERIC_EQ(asNumber(res[2]), asNumber(floor(mix[2])));


	auto res2 = floor(vecScalar);
	EXPECT_DOUBLE_EQ(res2.getScalarValue(), floor(vecScalar.getScalarValue()));
}

TEST(TestFuncions, cdfnorminv)
{
	
	std::vector<Numeric> mix{ asNumber(1.0 / 4.0), asNumber(1. / 2.0), asNumber(1./ 3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(1.0 / 2.2));


	auto res = cdfnorminv<VecXX::INS>(Vec1);
	EXPECT_EQ(res.size(), 3);
	EXPECT_NUMERIC_EQ(asNumber(res[0]), asNumber(cdfnorminv(mix[0])));
	EXPECT_NUMERIC_EQ(asNumber(res[1]), asNumber(cdfnorminv(mix[1])));
	EXPECT_NUMERIC_EQ(asNumber(res[2]), asNumber(cdfnorminv(mix[2])));

	auto res2 = cdfnorminv(vecScalar);
	EXPECT_NUMERIC_EQ(res2.getScalarValue(), asNumber( cdfnorminv(vecScalar.getScalarValue()) ));


}

TEST(TestFuncions, cdfnorm)
{
	std::vector<Numeric> mix{ asNumber(1.0 / 4.0), asNumber(1./ 2.0), asNumber(1./ 3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(1.0 / 2.2));

	auto res = cdfnorm(Vec1);
	EXPECT_EQ(res.size(), 3);
/*
	EXPECT_DOUBLE_EQ(res[0], cdfnorm(mix[0]));
	EXPECT_DOUBLE_EQ(res[1], cdfnorm(mix[1]));
	EXPECT_DOUBLE_EQ(res[2], cdfnorm(mix[2]));

	auto res2 = cdfnorm(vecScalar);
	EXPECT_DOUBLE_EQ(res2.getScalarValue(), cdfnorm(vecScalar.getScalarValue()));

	*/

	EXPECT_NUMERIC_EQ(asNumber(res[0]), asNumber(cdfnorm(mix[0])));
	EXPECT_NUMERIC_EQ(asNumber(res[1]), asNumber(cdfnorm(mix[1])));
	EXPECT_NUMERIC_EQ(asNumber(res[2]), asNumber(cdfnorm(mix[2])));

	auto res2 = cdfnorm(vecScalar);
	EXPECT_NUMERIC_EQ(res2.getScalarValue(), asNumber(cdfnorm(vecScalar.getScalarValue())));

}


TEST(TestFuncions, cdfnormD)
{
	/*
	std::vector<double> mix{ 1.0 / 4.0,1 / 2.0,1 / 3.0 };
	//VecXX Vec2(mix);
	VecXX Vec1(mix);
	VecXX vecScalar(1.0 / 2.2);
	*/

	std::vector<Numeric> mix{ asNumber(1.0 / 4.0), asNumber(1. / 2.0), asNumber(1. / 3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(1.0 / 2.2));

	auto res = cdfnormD(Vec1);
	EXPECT_EQ(res.size(), 3);

	EXPECT_DOUBLE_EQ(res[0], cdfnormD(mix[0]));
	EXPECT_DOUBLE_EQ(res[1], cdfnormD(mix[1]));
	EXPECT_DOUBLE_EQ(res[2], cdfnormD(mix[2]));

	auto res2 = cdfnormD(vecScalar);
	EXPECT_DOUBLE_EQ(res2.getScalarValue(), cdfnormD(vecScalar.getScalarValue()));
}


TEST(TestFuncions, unitaryMinus)
{

/*
	std::vector<double> mix{ 1.0 / 4.0,1 / 2.0,1 / 3.0 };
	//VecXX Vec2(mix);
	VecXX Vec1(mix);
	VecXX vecScalar(1.0 / 2.2);
*/
	std::vector<Numeric> mix{ asNumber(1.0 / 4.0), asNumber(1. / 2.0), asNumber(1. / 3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(1.0 / 2.2));

	auto res = -Vec1;
	EXPECT_EQ(res.size(), 3);
/*
	EXPECT_EQ(res[0], -(mix[0]));
	EXPECT_EQ(res[1], -(mix[1]));
	EXPECT_EQ(res[2], -(mix[2]));
	*/

	EXPECT_NUMERIC_EQ(asNumber(res[0]), asNumber(-(mix[0])));
	EXPECT_NUMERIC_EQ(asNumber(res[1]), asNumber(-(mix[1])));
	EXPECT_NUMERIC_EQ(asNumber(res[2]), asNumber(-(mix[2])));

	auto res2 = -(vecScalar);
	EXPECT_EQ(res2.getScalarValue(), -(vecScalar.getScalarValue()));
}

TEST(TestFuncions, ceil)
{
/*
	std::vector<double> mix{ 1.0 / 4.0,1 / 2.0,1 / 3.0 };
	//VecXX Vec2(mix);
	VecXX Vec1(mix);
	VecXX vecScalar(1.0 / 2.2);
*/
	std::vector<Numeric> mix{ asNumber(1.0/4.0), asNumber(1.0/2.0), asNumber(1.0/3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(1.0 / 2.2));


	auto res = ceil(Vec1);
	EXPECT_EQ(res.size(), 3);
	/*
	EXPECT_EQ(res[0], ceil(mix[0]));
	EXPECT_EQ(res[1], ceil(mix[1]));
	EXPECT_EQ(res[2], ceil(mix[2]));
	*/
	EXPECT_NUMERIC_EQ(asNumber(res[0]), asNumber(ceil(mix[0])));
	EXPECT_NUMERIC_EQ(asNumber(res[1]), asNumber(ceil(mix[1])));
	EXPECT_NUMERIC_EQ(asNumber(res[2]), asNumber(ceil(mix[2])));


	auto res2 = ceil(vecScalar);
	EXPECT_NUMERIC_EQ(res2.getScalarValue(), ceil(vecScalar.getScalarValue()));
}


TEST(TestFuncions, sqrt)
{
/*
	std::vector<double> mix{ 1.0 / 4.0,1 / 2.0,1 / 3.0 };
	//VecXX Vec2(mix);
	VecXX Vec1(mix);
	VecXX vecScalar(1.0 / 2.2);
*/
	std::vector<Numeric> mix{ asNumber(1.0 / 4.0), asNumber(1.0 / 2.0), asNumber(1.0 / 3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(1.0 / 2.2));


	auto res = sqrt(Vec1);
	EXPECT_EQ(res.size(), 3);
	/*
	EXPECT_EQ(res[0], sqrt(mix[0]));
	EXPECT_EQ(res[1], sqrt(mix[1]));
	EXPECT_EQ(res[2], sqrt(mix[2]));
	*/
	EXPECT_NUMERIC_EQ(asNumber(res[0]), asNumber(sqrt(mix[0])));
	EXPECT_NUMERIC_EQ(asNumber(res[1]), asNumber(sqrt(mix[1])));
	EXPECT_NUMERIC_EQ(asNumber(res[2]), asNumber(sqrt(mix[2])));


	auto res2 = sqrt(vecScalar);
	EXPECT_EQ(res2.getScalarValue(), sqrt(vecScalar.getScalarValue()));
}


TEST(TestFuncions, pow)
{
/*
	std::vector<double> mix{ 1.0 / 4.0,1 / 2.0,1 / 3.0 };
	//VecXX Vec2(mix);
	VecXX Vec1(mix);
	VecXX vecScalar(1.0 / 2.0);
*/


	std::vector<Numeric> mix{ asNumber(1.0/4.0), asNumber(1.0/2.0), asNumber(1.0/3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(1.0/2.2));

	auto res = pow(Vec1, asNumber(0.5) );
	EXPECT_EQ(res.size(), 3);
/*
	EXPECT_DOUBLE_EQ(res[0], sqrt(mix[0]));
	EXPECT_DOUBLE_EQ(res[1], sqrt(mix[1]));
	EXPECT_DOUBLE_EQ(res[2], sqrt(mix[2]));
*/

	EXPECT_NUMERIC_EQ(res[0], asNumber(sqrt(mix[0])));
	EXPECT_NUMERIC_EQ(res[1], asNumber(sqrt(mix[1])));
	EXPECT_NUMERIC_EQ(res[2], asNumber(sqrt(mix[2])));

	//std::vector<Numeric> mix1{ 1.0, 2.0, 3.0 };
	std::vector<Numeric> mix1{ asNumber(1.0 ), asNumber( 2.0), asNumber(3.0) };
	VecXX powers(mix1);


	auto res1 = pow(Vec1, powers);

	EXPECT_EQ(res1.size(), 3);

	EXPECT_NUMERIC_EQ(res1[0], asNumber( Vec1[0]));
	EXPECT_NUMERIC_EQ(res1[1], asNumber(Vec1[1] * Vec1[1]));
	EXPECT_NUMERIC_EQ(res1[2], asNumber (Vec1[2] * Vec1[2]* Vec1[2]));


}



TEST(TestFuncions, max)
{
/*
	std::vector<double> mix{ 1.0 / 4.0,1 / 2.0,1 / 3.0 };
	//VecXX Vec2(mix);
	VecXX Vec1(mix);
	VecXX vecScalar(1.0 / 2.2);
*/

	std::vector<Numeric> mix{ asNumber(1.0 / 4.0), asNumber(1.0/ 2.0), asNumber(1.0/ 3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(1.0 / 2.2));

	VecXX Vec2 = Vec1 + 1.0;
	VecXX Vec3 = Vec1 - 1.0;

	auto res = max(Vec1, Vec2);
	EXPECT_EQ(res.size(), 3);

/*
	EXPECT_EQ(res[0], Vec2[0]);
	EXPECT_EQ(res[1], Vec2[1]);
	EXPECT_EQ(res[2], Vec2[2]);
*/
	EXPECT_NUMERIC_EQ(asNumber(res[0]), asNumber(Vec2[0]));
	EXPECT_NUMERIC_EQ(asNumber(res[1]), asNumber(Vec2[1]));
	EXPECT_NUMERIC_EQ(asNumber(res[2]), asNumber(Vec2[2]));



	auto res3 = max(vecScalar, asNumber( -100.0) );
	EXPECT_NUMERIC_EQ(res3.getScalarValue(),  asNumber( 1.0 / 2.2) );

	auto res2 = max(Vec1, Vec3);
	EXPECT_EQ(res.size(), 3);

	//EXPECT_EQ(res2[0], Vec1[0]);
	//EXPECT_EQ(res2[1], Vec1[1]);
	//EXPECT_EQ(res2[2], Vec1[2]);

	EXPECT_NUMERIC_EQ(asNumber(res2[0]), asNumber(Vec1[0]));
	EXPECT_NUMERIC_EQ(asNumber(res2[1]), asNumber(Vec1[1]));
	EXPECT_NUMERIC_EQ(asNumber(res2[2]), asNumber(Vec1[2]));



	//VecxD testVals(Vec2, Vec3);
	//VecxD testVals2 = testVals * 10.0;

	//auto res4 = max(testVals, testVals2);


}



TEST(TestFuncions, min)
{
/*
	std::vector<double> mix{ 1.0 / 4.0,1 / 2.0,1 / 3.0 };
	//VecXX Vec2(mix);
	VecXX Vec1(mix);
	VecXX vecScalar(1.0 / 2.2);
*/
	std::vector<Numeric> mix{ asNumber(1.0 / 4.0), asNumber(1.0 / 2.0), asNumber(1.0 / 3.0) };
	VecXX Vec1(mix);
	VecXX vecScalar(asNumber(1.0 / 2.2));



	VecXX Vec2 = Vec1 + asNumber(1.0);
	VecXX Vec3 = Vec1 + asNumber(2.0);

	auto res = min(Vec1, Vec2);
	EXPECT_EQ(res.size(), 3);

/*
	EXPECT_EQ(res[0], Vec1[0]);
	EXPECT_EQ(res[1], Vec1[1]);
	EXPECT_EQ(res[2], Vec1[2]);
*/
	EXPECT_NUMERIC_EQ(asNumber(res[0]), asNumber(Vec1[0]));
	EXPECT_NUMERIC_EQ(asNumber(res[1]), asNumber(Vec1[1]));
	EXPECT_NUMERIC_EQ(asNumber(res[2]), asNumber(Vec1[2]));

	res = min(Vec2, Vec1);
	EXPECT_NUMERIC_EQ(asNumber(res[0]), asNumber(Vec1[0]));
	EXPECT_NUMERIC_EQ(asNumber(res[1]), asNumber(Vec1[1]));
	EXPECT_NUMERIC_EQ(asNumber(res[2]), asNumber(Vec1[2]));


	auto res3 = min(vecScalar, asNumber (-100.0));
	EXPECT_EQ(res3.getScalarValue(), asNumber(-100.0));

	auto res2 = min(Vec3, Vec1);
	EXPECT_EQ(res.size(), 3);
	/*
	EXPECT_EQ(res2[0], Vec2[0]);
	EXPECT_EQ(res2[1], Vec2[1]);
	EXPECT_EQ(res2[2], Vec2[2]);
	*/
	EXPECT_NUMERIC_EQ(asNumber(res2[0]), asNumber(Vec1[0]));
	EXPECT_NUMERIC_EQ(asNumber(res2[1]), asNumber(Vec1[1]));
	EXPECT_NUMERIC_EQ(asNumber(res2[2]), asNumber(Vec1[2]));


}


TEST(TestFuncions, multAdd)
{
	//std::vector<Numeric> mix{ 1.0,  2.0, 3.0 };
	std::vector<Numeric> mix{ asNumber(1.0 ), asNumber(2.0), asNumber(3.0) };
	VecXX Vec2(mix);

	std::vector<Numeric> one{ asNumber(1.0), asNumber(1.0), asNumber(1.0) };

	VecXX Vec1(one);

	VecXX Vec3 = Vec1 + asNumber(1.0);

	auto res = FMA(Vec2, Vec1, Vec3);
	EXPECT_EQ(res.size(), 3);

	EXPECT_NUMERIC_EQ(res[0], Vec1[0] * Vec2[0] + Vec3[0]);
	EXPECT_NUMERIC_EQ(res[1], Vec1[1] * Vec2[1] + Vec3[1]);
	EXPECT_NUMERIC_EQ(res[2], Vec1[2] * Vec2[2] + Vec3[2]);

	auto res1 = FMA(Vec1, Vec2, asNumber(1.0));
	EXPECT_NUMERIC_EQ(res1[0], Vec1[0] * Vec2[0] + asNumber(1.0));
	EXPECT_NUMERIC_EQ(res1[1], Vec1[1] * Vec2[1] + asNumber(1.0));
	EXPECT_NUMERIC_EQ(res1[2], Vec1[2] * Vec2[2] + asNumber(1.0));

	auto res2 = FMA(Vec2, asNumber(1.0), 3.0);
	EXPECT_EQ(res2.size(), 3);

	EXPECT_NUMERIC_EQ(res2[0], Vec2[0] * asNumber(1.0) + asNumber(3.0));
	EXPECT_NUMERIC_EQ(res2[1], Vec2[1] * asNumber(1.0) + asNumber(3.0));
	EXPECT_NUMERIC_EQ(res2[2], Vec2[2] * asNumber(1.0) + asNumber(3.0));


}


TEST(TestFuncions, testLamdaOneD)
{
	//std::vector<Numeric> mix{ 1.0,  2.0, 3.0 };
	std::vector<Numeric> mix{ asNumber(1.0), asNumber(2.0), asNumber(3.0) };
	VecXX Vec2(mix);

	//std::vector<Numeric> one{ 1.0,  1.0, 1.0 };
	std::vector<Numeric> one{ asNumber(1.0), asNumber(1.0), asNumber(1.0) };

	VecXX Vec1(one);

	auto myFunc = [](auto X) {return (X * X); };

	auto res = ApplyLambda(Vec2, myFunc);
	EXPECT_EQ(res.size(), 3);


	EXPECT_EQ(res[0], Vec2[0] * Vec2[0]);
	EXPECT_EQ(res[1], Vec2[1] * Vec2[1]);
	EXPECT_EQ(res[2], Vec2[2] * Vec2[2]);


}


TEST(TestFuncions, testLamdaTwoD)
{
	//std::vector<Numeric> mix{ 1.0,  2.0, 3.0 };
	std::vector<Numeric> mix{ asNumber(1.0), asNumber(2.0), asNumber(3.0) };
	VecXX Vec2(mix);

	std::vector<Numeric> one{ 1.0,  1.0, 1.0 };

	VecXX Vec1(one);

	auto myFunc = [](auto X, auto Y) {return (X * X + Y * Y + X * Y); };

	auto res = ApplyLambda2(Vec1, Vec2, myFunc);
	EXPECT_EQ(res.size(), 3);


	EXPECT_EQ(res[0], 3.0);
	EXPECT_EQ(res[1], 7.0);
	EXPECT_EQ(res[2], 13.0);


}

/*
TEST(TestFuncions, testBoolLamdaTwoD)
{
	std::vector<Numeric> mix{ 1.0,  2.0, 0.1 };
	VecXX Vec2(mix);

	std::vector<Numeric> one{ 1.0,  1.0, 1.0 };

	VecXX Vec1(one);
	auto myFunc = [](const auto& X, const auto& Y) {return ((X * X - Y * Y + X * Y) <  X * X); };

//	auto wrap = CustomBinaryBool(myFunc);
	auto res = ApplyBoolLambda2(Vec1, Vec2, myFunc);

	
	EXPECT_EQ(res.size(), 3);

	auto res0 = res[0];
	auto res1 = res[1];
	auto res2 = res[2];

	EXPECT_FALSE(res0);
	EXPECT_TRUE(res1);
	EXPECT_FALSE(res2);

	//reverse arguments
	res = ApplyBoolLambda2(Vec2, Vec1, myFunc);


	EXPECT_EQ(res.size(), 3);

	res0 = res[0];
	res1 = res[1];
	res2 = res[2];

	EXPECT_FALSE(res0);
	EXPECT_FALSE(res1);
	EXPECT_TRUE(res2);



	///  //////////////////////////////

	auto vcl_lte = (Vec1 <= Vec2);

	EXPECT_EQ(vcl_lte.size(), 3);
	res0 = vcl_lte[0];
	res1 = vcl_lte[1];
	res2 = vcl_lte[2];

	EXPECT_TRUE(res0);
	EXPECT_TRUE(res1);
	EXPECT_FALSE(res2);

	vcl_lte = (Vec2 <= Vec1);

	EXPECT_EQ(vcl_lte.size(), 3);
	res0 = vcl_lte[0];
	res1 = vcl_lte[1];
	res2 = vcl_lte[2];

	EXPECT_TRUE(res0);
	EXPECT_FALSE(res1);
	EXPECT_TRUE(res2);


	auto vcl_eq = (Vec1 == Vec2);

	EXPECT_EQ(vcl_lte.size(), 3);
	res0 = vcl_eq[0];
	res1 = vcl_eq[1];
	res2 = vcl_eq[2];

	EXPECT_TRUE(res0);
	EXPECT_FALSE(res1);
	EXPECT_FALSE(res2);


	vcl_eq = (Vec2 == Vec1);

	EXPECT_EQ(vcl_lte.size(), 3);
	res0 = vcl_eq[0];
	res1 = vcl_eq[1];
	res2 = vcl_eq[2];

	EXPECT_TRUE(res0);
	EXPECT_FALSE(res1);
	EXPECT_FALSE(res2);


	VecXX scalar1(asNumber(100.));
	VecXX scalar2(asNumber(100.));

	

	auto vc2a = (scalar2 == scalar1);



	EXPECT_TRUE(vc2a.isScalar()); // wrong always false, do we allow scala bool vectors

	//vc2a.

	auto vc3 = (Vec2 <= asNumber(0.2));

	EXPECT_EQ(vc3.size(), 3);

	res0 = vc3[0];
	res1 = vc3[1];
	res2 = vc3[2];

	EXPECT_FALSE(res0);
	EXPECT_FALSE(res1);
	EXPECT_TRUE(res2);


	auto vc4 = (Vec2 >= asNumber(0.2));

	EXPECT_EQ(vc3.size(), 3);

	res0 = vc4[0];
	res1 = vc4[1];
	res2 = vc4[2];

	EXPECT_TRUE(res0);
	EXPECT_TRUE(res1);
	EXPECT_FALSE(res2);

	/////////////////////////////////

	res = ApplyBoolLambda2(scalar1, Vec1, myFunc);
	EXPECT_EQ(res.size(), 3);

	res0 = res[0];
	res1 = res[1];
	res2 = res[2];


	EXPECT_FALSE(res0);
	EXPECT_FALSE(res1);
	EXPECT_FALSE(res2);


	res = ApplyBoolLambda2( Vec1, scalar1, myFunc);
	EXPECT_EQ(res.size(), 3);

	res0 = res[0];
	res1 = res[1];
	res2 = res[2];

	EXPECT_TRUE(res0);
	EXPECT_TRUE(res1);
	EXPECT_TRUE(res2);
		


}

*/

TEST(TestFuncions, testIff)
{
	//std::vector<Numeric>  mix{ 1.0,  2.0, 3.0 };
	std::vector<Numeric> mix{ asNumber(1.0), asNumber(2.0), asNumber(3.0) };
	VecXX Vec2(mix);

	std::vector<Numeric>  one{ asNumber(1.0),asNumber(1.0),asNumber(1.0) };// 1.0, 1.0, 1.0

	std::vector<Numeric>  other{ asNumber(8.0),asNumber(7.0),asNumber(6.0) }; // other{ 8.0,  7.0, 6.0 };

	VecXX scalar = asNumber(999.);
	VecXX scalar2 = asNumber(222.);

	VecXX Vec1(one);
	VecXX otherFalse(other);

	auto boolCond = Vec2 > (Vec1 + asNumber(1.0) );

	auto res = select(boolCond, Vec1, -Vec1);

	EXPECT_EQ(res.size(), 3);
	//EXPECT_EQ(res[0], -1.0);
	//EXPECT_EQ(res[1], -1.0);
	//EXPECT_EQ(res[2], 1.0);

	EXPECT_NUMERIC_EQ(res[0], asNumber(-1.0));
	EXPECT_NUMERIC_EQ(res[1], asNumber(-1.0));
	EXPECT_NUMERIC_EQ(res[2], asNumber(1.0) );



	res = select(boolCond, Vec1, otherFalse);

	EXPECT_EQ(res.size(), 3);
	//EXPECT_EQ(res[0], 8.0);
	//EXPECT_EQ(res[1], 7.0);
	//EXPECT_EQ(res[2], 1.0);

	EXPECT_NUMERIC_EQ(res[0], asNumber(8.0));
	EXPECT_NUMERIC_EQ(res[1], asNumber(7.0));
	EXPECT_NUMERIC_EQ(res[2], asNumber(1.0));


	res = select(boolCond, scalar, scalar2);

	EXPECT_EQ(res.size(), 3);
//	EXPECT_EQ(res[0], 222.0);
//	EXPECT_EQ(res[1], 222.0);
//	EXPECT_EQ(res[2], 999.0);

	EXPECT_NUMERIC_EQ(res[0], asNumber(222.0));
	EXPECT_NUMERIC_EQ(res[1], asNumber(222.0));
	EXPECT_NUMERIC_EQ(res[2], asNumber(999.0));




	res = select(boolCond, asNumber( 999.0), asNumber(222.0));

	EXPECT_EQ(res.size(), 3);
	//EXPECT_EQ(res[0], 222.0);
	//EXPECT_EQ(res[1], 222.0);
	//EXPECT_EQ(res[2], 999.0);

	EXPECT_NUMERIC_EQ(res[0], asNumber(222.0));
	EXPECT_NUMERIC_EQ(res[1], asNumber(222.0));
	EXPECT_NUMERIC_EQ(res[2], asNumber(999.0));

	res = select(boolCond, Vec2, asNumber(333.1));

	EXPECT_EQ(res.size(), 3);
	//EXPECT_EQ(res[0], 333.1);
	//EXPECT_EQ(res[1], 333.1);
	//EXPECT_EQ(res[2], 3.0);

	EXPECT_NUMERIC_EQ(res[0], asNumber(333.1));
	EXPECT_NUMERIC_EQ(res[1], asNumber(333.1));
	EXPECT_NUMERIC_EQ(res[2], asNumber(3.0));


	res = select(boolCond, asNumber(333.1), Vec2);

	EXPECT_EQ(res.size(), 3);
	//EXPECT_EQ(res[0], 1.0);
	//EXPECT_EQ(res[1], 2.0);
	//EXPECT_EQ(res[2], asNumber(333.1));

	EXPECT_NUMERIC_EQ(res[0], asNumber(1.));
	EXPECT_NUMERIC_EQ(res[1], asNumber(2.));
	EXPECT_NUMERIC_EQ(res[2], asNumber(333.1));





	//auto res = iff(Vec2 > Vec1 + 1.0, Vec1, -Vec1);
	res = iff(boolCond, Vec1, -Vec1);
	//res = iff(Vec2 > Vec1 + 1.0, Vec1, -Vec1);
	EXPECT_EQ(res.size(), 3);
	//EXPECT_EQ(res[0], -1.0);
	//EXPECT_EQ(res[1], -1.0);
	//EXPECT_EQ(res[2], 1.0);

	EXPECT_NUMERIC_EQ(res[0], asNumber(-1.));
	EXPECT_NUMERIC_EQ(res[1], asNumber(-1.));
	EXPECT_NUMERIC_EQ(res[2], asNumber(1.0));
		


	//res = iff(boolCond, Vec1, -Vec1);
	res = iff(Vec2 > (Vec1 + asNumber(1.0)), Vec1, -Vec1);
	EXPECT_EQ(res.size(), 3);
	//EXPECT_EQ(res[0], -1.0);
	//EXPECT_EQ(res[1], -1.0);
	//EXPECT_EQ(res[2], 1.0);

	EXPECT_NUMERIC_EQ(res[0], asNumber(-1.));
	EXPECT_NUMERIC_EQ(res[1], asNumber(-1.));
	EXPECT_NUMERIC_EQ(res[2], asNumber(1.0));


}
	

Numeric getNull(Numeric)
{
	return 0.0;
}


VecXX::INS getNull(VecXX::INS)
{
	return VecXX::INS(0.0);
}

TEST(TestFuncions, accum1)
{
	std::vector<Numeric> oneOone(101,asNumber(1.0) );
	VecXX Vec(oneOone);

	auto sum = [](auto lhs, auto rhs) {return rhs + lhs; };
	Numeric sumRes = ApplyAccumulate(Vec, sum,0.0);



	auto KhanAddV = [c = getNull(VecXX::INS(0.0)), sum = getNull(VecXX::INS(0.0))](auto lhs, auto rhs) mutable
	{
		auto y = rhs - c;
		auto t = sum + y;
		c = (t - sum);
		c = c - y;
		sum = t;
		return t;
	};

	auto KhanAddS = [c = getNull(0.0), sum = getNull(0.0)](auto lhs, auto rhs) mutable
	{
		auto y = rhs - c;
		auto t = sum + y;
		c = (t - sum);
		c = c - y;
		sum = t;
		return t;
	};


	auto Add = [](auto lhs, auto rhs) { return lhs + rhs; };

	auto sum1 = ApplyAccumulate2(Vec, KhanAddV, KhanAddS, asNumber(0.0));
	auto sum2 = ApplyAccumulate2(Vec, Add, Add, asNumber(0.0));



	/*
	std::vector<double> one{ 1.0,  1.0, 1.0 };

	VecXX Vec1(one);

	VecXX Vec3 = Vec1 + 1.0;

	auto res = FMA(Vec2, Vec1, Vec3);
	EXPECT_EQ(res.size(), 3);

	EXPECT_EQ(res[0], Vec1[0] * Vec2[0] + Vec3[0]);
	EXPECT_EQ(res[1], Vec1[1] * Vec2[1] + Vec3[1]);
	EXPECT_EQ(res[2], Vec1[2] * Vec2[2] + Vec3[2]);

	auto res1 = FMA(Vec1, Vec2, 1.0);
	EXPECT_EQ(res1[0], Vec1[0] * Vec2[0] + 1.0);
	EXPECT_EQ(res1[1], Vec1[1] * Vec2[1] + 1.0);
	EXPECT_EQ(res1[2], Vec1[2] * Vec2[2] + 1.0);

	auto res2 = FMA(Vec2, 1.0, 3.0);
	EXPECT_EQ(res2.size(), 3);

	EXPECT_EQ(res2[0], Vec2[0] * 1.0 + 3.0);
	EXPECT_EQ(res2[1], Vec2[1] * 1.0 + 3.0);
	EXPECT_EQ(res2[2], Vec2[2] * 1.0 + 3.0);

	*/


}
