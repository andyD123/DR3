#include "pch.h"


#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/vec_bool_d.h"
#include "../Vectorisation/VecX/vec_double.h"
#include  "../Vectorisation/VecX/alloc_policy.h"
#include "../Vectorisation/VecX/binary_unitary_operations.h"
#include "../Vectorisation/VecX/target_name_space.h"

#include <numeric>



using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
//using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;







template <typename INS_VEC>
void evalApplyTernaryOperation_VFF(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] =  i; }
		VecXX test(v);

		auto resAcc = FMA(test, VecXX::scalar(3.0), VecXX::scalar(1.0));

		double ii = 0.;
		for (double x : resAcc)
		{
			EXPECT_NEAR(x,ii * 3.0 + 1.0, testEpsilon);
			ii++;
		}	
	}

}


template <typename INS_VEC>
void evalApplyTernaryOperation_VVF(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }
		VecXX test(v);

		for (int i = 0; i < SZ; i++) { v[i] = i+1.1; }
		VecXX mult(v);

		auto resAcc = FMA(test, mult, VecXX::scalar(1.0));

		double ii = 0.;
		int i = 0;
		for (double x : resAcc)
		{
			EXPECT_NEAR(x, ii * mult[i] + 1.0, testEpsilon);
			ii++;
			i++;
		}	

	}

}



template <typename INS_VEC>
void evalApplyTernaryOperation_VVV(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 0.0);
		std::vector<double>  z(SZ, 0.0);
		std::vector<double>  y(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }
		VecXX test(v);

		for (int i = 0; i < SZ; i++) { z[i] = i + 1/6.1; }
		VecXX mult(z);

		for (int i = 0; i < SZ; i++) { y[i] = i + 1.1; }
		VecXX addd(v);

		auto resAcc = FMA(test, mult, addd);

		double ii = 0.;
		int i = 0;
		for (double x : resAcc)
		{
			EXPECT_NEAR(x, ii * mult[i] + addd[i], testEpsilon);
			ii++;
			i++;
		}	

	}

}


template <typename INS_VEC>
void evalApplyTernaryOperation_VFV(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }
		VecXX test(v);

		for (int i = 0; i < SZ; i++) { v[i] = i + 1.1; }
		VecXX addd(v);

		auto resAcc = FMA(test, VecXX::scalar(1.1), addd);

		double ii = 0.;
		int i = 0;
		for (double x : resAcc)
		{
			EXPECT_NEAR(x, ii * 1.1 + addd[i], testEpsilon);
			ii++;
			i++;
		}

	}

}

//////////////////// start with scalar ////////


template <typename INS_VEC>
void evalApplyTernaryOperation_FFF(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }
		VecXX test(v);

		test = 1.123; //scalar

		auto resAcc = FMA(test, VecXX::scalar(3.0), VecXX::scalar(1.0));

		double ii = 0.;
		for (double x : resAcc)
		{
			EXPECT_NEAR(x, 1.123 * 3.0 + 1.0, testEpsilon);
			ii++;
		}
	}

}


template <typename INS_VEC>
void evalApplyTernaryOperation_FVF(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }
		VecXX test(v);
		test = 1.123; //scalar

		for (int i = 0; i < SZ; i++) { v[i] = i + 1.1; }
		VecXX mult(v);

		auto resAcc = FMA(test, mult, VecXX::scalar(1.0));

		double ii = 0.;
		int i = 0;
		for (double x : resAcc)
		{
			EXPECT_NEAR(x, 1.123 * mult[i] + 1.0, testEpsilon);
			ii++;
			i++;
		}

	}

}



template <typename INS_VEC>
void evalApplyTernaryOperation_FVV(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 0.0);
		std::vector<double>  z(SZ, 0.0);
		std::vector<double>  y(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }
		VecXX test(v);
		test = 1.123; //scalar


		for (int i = 0; i < SZ; i++) { z[i] = i + 1 / 6.1; }
		VecXX mult(z);
		

		for (int i = 0; i < SZ; i++) { y[i] = i + 1.1; }
		VecXX addd(v);

		auto resAcc = FMA(test, mult, addd);

		double ii = 0.;
		int i = 0;
		for (double x : resAcc)
		{
			EXPECT_NEAR(x, 1.123 * mult[i] + addd[i], testEpsilon);
			ii++;
			i++;
		}

	}

}


template <typename INS_VEC>
void evalApplyTernaryOperation_FFV(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;

	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{
		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }
		VecXX test(v);
		test = 1.123; //scalar

		for (int i = 0; i < SZ; i++) { v[i] = i + 1.1; }
		VecXX addd(v);

		auto resAcc = FMA(test, VecXX::scalar(1.1), addd);

		double ii = 0.;
		int i = 0;
		for (double x : resAcc)
		{
			EXPECT_NEAR(x, 1.123 * 1.1 + addd[i], testEpsilon);
			ii++;
			i++;
		}
	}

}







TEST(TestBinaryAndUnitaryFunctions, ApplyTernaryOperation)
{

	for (int j = 1; j < 102; ++j)
	{
		evalApplyTernaryOperation_VFF<typename VecXX::INS>(1, 1+j);
	}

	for (int j = 1; j < 102; ++j)
	{
		evalApplyTernaryOperation_VVF<typename VecXX::INS>(1, 1 + j);
	}


	for (int j = 1; j < 102; ++j)
	{
		evalApplyTernaryOperation_VFV<typename VecXX::INS>(1, 1 + j);
	}

	for (int j = 1; j < 102; ++j)
	{
		evalApplyTernaryOperation_VVV<typename VecXX::INS>(1, 1 + j);
	}
	

	//
	
	for (int j = 1; j < 102; ++j)
	{
		evalApplyTernaryOperation_FFF<typename VecXX::INS>(1, 1 + j);
	}

	for (int j = 1; j < 102; ++j)
	{
		evalApplyTernaryOperation_FVF<typename VecXX::INS>(1, 1 + j);
	}


	for (int j = 1; j < 102; ++j)
	{
		evalApplyTernaryOperation_FFV<typename VecXX::INS>(1, 1 + j);
	}

	for (int j = 1; j < 102; ++j)
	{
		evalApplyTernaryOperation_FVV<typename VecXX::INS>(1, 1 + j);
	}


	
}


template <typename INS_VEC>
void evalApplyUnitaryOperation(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;
	//int SZ = 1000;
	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{	

		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = (0==i%2)? -i:i;  }
		VecXX test(v);

		//calling unitary minus
		auto negativeRes = -test;

		int i = 0;
		for (double x : negativeRes)
		{
			EXPECT_NEAR(x, -test[i], testEpsilon);
			i++;
		}
	}
}


template <typename INS_VEC>
void evalApplyUnitaryOperation1(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;
	//int SZ = 1000;
	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{

		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = (0 == i % 2) ? -i : i; }
		VecXX test(v);

		//calling unitary minus
		auto negativeRes = ApplyUnitaryOperation1<  VecXX::INS, DR_CUBED::_unitaryMinus<VecXX::INS > >(test, DR_CUBED::_unitaryMinus<VecXX::INS >());

		int i = 0;
		for (double x : negativeRes)
		{
			EXPECT_NEAR(x, -test[i], testEpsilon);
			i++;
		}
	}
}


template <typename INS_VEC>
void evalApplyUnitaryOperation_4(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;
	//int SZ = 1000;
	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{

		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = (0 == i % 2) ? -i : i; }
		VecXX test(v);

		//calling unitary minus
		auto negativeRes = ApplyUnitaryOperation<  VecXX::INS, DR_CUBED::_unitaryMinus<VecXX::INS > >(test);

		int i = 0;
		for (double x : negativeRes)
		{
			EXPECT_NEAR(x, -test[i], testEpsilon);
			i++;
		}
	}
}




TEST(TestBinaryAndUnitaryFunctions, ApplyUnitaryOperation)
{
	const double testEpsilon = 0.000001;
	{
		auto  testValue = VecXX::scalar(1.1);
		VecXX res = ApplyUnitaryOperation<  VecXX::INS, DR_CUBED::_unitaryMinus<VecXX::INS > >(testValue);
		EXPECT_TRUE(res.isScalar());
		EXPECT_NEAR(res.getScalarValue(), -1.1, testEpsilon);

		for (int j = 1; j < 102; ++j)
		{
			evalApplyUnitaryOperation< VecXX::INS>(1, 102 + j);
		}
	}

	{
		VecXX  testValue = VecXX::scalar(1.1);
		VecXX res = ApplyUnitaryOperation<  VecXX::INS, DR_CUBED::_unitaryMinus<VecXX::INS > >(testValue, DR_CUBED::_unitaryMinus<VecXX::INS >());
		EXPECT_TRUE(res.isScalar());
		EXPECT_NEAR(res.getScalarValue(), -1.1, testEpsilon);
		/*
		for (int j = 1; j < 102; ++j)
		{
			evalApplyUnitaryOperation< VecXX::INS>(1, 102 + j);
		}
		*/
	}


	{

		VecXX  testValue = VecXX::scalar(1.1);
		VecXX res = ApplyUnitaryOperation1<  VecXX::INS, DR_CUBED::_unitaryMinus<VecXX::INS > >(testValue, DR_CUBED::_unitaryMinus<VecXX::INS >());
		EXPECT_TRUE(res.isScalar());
		EXPECT_NEAR(res.getScalarValue(), -1.1, testEpsilon);

		for (int j = 1; j < 102; ++j)
		{
			evalApplyUnitaryOperation1< VecXX::INS>(1, 102 + j);
		}
	}


	{

		VecXX  testValue = VecXX::scalar(1.1);
		VecXX res = ApplyUnitaryOperation<  VecXX::INS, DR_CUBED::_unitaryMinus<VecXX::INS > >(testValue);
		EXPECT_TRUE(res.isScalar());
		EXPECT_NEAR(res.getScalarValue(), -1.1, testEpsilon);

		for (int j = 1; j < 102; ++j)
		{
			evalApplyUnitaryOperation_4< VecXX::INS>(1, 102 + j);
		}
	}


}




template <typename INS_VEC>
void evalApplyBinaryOperationVV(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;
	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{

		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }// (0 == i % 2) ? -i : i;

		VecXX Lhs(v);

		VecXX Rhs = Lhs * 3.2 - 1.0;

		auto res =ApplyBinaryOperation<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(Lhs,Rhs);

		int i = 0;
		for (double x : res)
		{
			EXPECT_NEAR(x, v[i]*(1+3.2) -1.0   , testEpsilon);
			i++;
		}
	}
}


template <typename INS_VEC>
void evalApplyBinaryOperationVF(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;
	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{

		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }

		VecXX Lhs(v);

		VecXX Rhs = 10.1;

		auto res = ApplyBinaryOperation<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(Lhs, Rhs);

		int i = 0;
		for (double x : res)
		{
			EXPECT_NEAR(x, v[i] +10.1, testEpsilon);
			i++;
		}
	}
}


template <typename INS_VEC>
void evalApplyBinaryOperationFV(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;
	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{

		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }

		VecXX Lhs(v);

		VecXX Rhs = 10.1;

		auto res = ApplyBinaryOperation<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(Rhs,Lhs);

		int i = 0;
		for (double x : res)
		{
			EXPECT_NEAR(x, v[i] + 10.1, testEpsilon);
			i++;
		}
	}
}




TEST(TestBinaryAndUnitaryFunctions, ApplyBinaryOperationVV)
{
	const double testEpsilon = 0.000001;
	{
		VecXX  LHS = VecXX::scalar(1.1);
		VecXX  RHS = VecXX::scalar(3.1);



		auto res = ApplyBinaryOperation<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(LHS, RHS);
		EXPECT_TRUE(res.isScalar());
		EXPECT_NEAR(res.getScalarValue(), 4.2, testEpsilon);
		////////////////// both inputs scalar /////////////

		int SZ = 13;
		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }
		VecXX smallVec(v);

		res = ApplyBinaryOperation<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(smallVec, RHS);
		EXPECT_TRUE(!res.isScalar());

		int i = 0;
		for (double x : res)
		{
			EXPECT_NEAR(x, smallVec[i] + 3.1, testEpsilon);
			i++;
		}
		////////////////////////////

		res = ApplyBinaryOperation<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(LHS, smallVec);
		EXPECT_TRUE(!res.isScalar());
		i = 0;
		for (double x : res)
		{
			EXPECT_NEAR(x, smallVec[i] + 1.1, testEpsilon);
			i++;
		}

		//////////////////////////////////


		for (int j = 1; j < 102; ++j)
		{
			evalApplyBinaryOperationVV< VecXX::INS>(1, 102 + j);
		}
	}
}






TEST(TestBinaryAndUnitaryFunctions, ApplyBinaryOperationVF)
{
	const double testEpsilon = 0.000001;

	{

		VecXX  LHS = VecXX::scalar(1.1);
		//auto  RHS = VecXX::scalar(3.1);

		int SZ = 13;
		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }
		VecXX smallVec(v);

		auto floatVal = VecXX::scalar(11.1);

		auto res = ApplyBinaryOperation<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(smallVec, floatVal);
		EXPECT_TRUE(!res.isScalar());
		int i = 0;
		for (double x : res)
		{
			EXPECT_NEAR(x, smallVec[i] + 11.1, testEpsilon);
			i++;
		}
		
		/////////////////////////////
		res = ApplyBinaryOperation<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(LHS, floatVal);
		EXPECT_TRUE(res.isScalar());
		EXPECT_NEAR(res.getScalarValue(), 12.2, testEpsilon);
		
		for (int j = 1; j < 102; ++j)
		{
			evalApplyBinaryOperationVF< VecXX::INS>(1, 102 + j);
		}
	}


	{

		VecXX  RHS = VecXX::scalar(1.1);

		int SZ = 13;
		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }
		VecXX smallVec(v);

		auto floatVal = VecXX::scalar(11.1);

		auto res = ApplyBinaryOperation<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >( floatVal, smallVec);
		EXPECT_TRUE(!res.isScalar());
		int i = 0;
		for (double x : res)
		{
			EXPECT_NEAR(x, smallVec[i] + 11.1, testEpsilon);
			i++;
		}

		/////////////////////////////
		res = ApplyBinaryOperation<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >( floatVal, RHS);
		EXPECT_TRUE(res.isScalar());
		EXPECT_NEAR(res.getScalarValue(), 12.2, testEpsilon);

		for (int j = 1; j < 102; ++j)
		{
			evalApplyBinaryOperationFV< VecXX::INS>(1, 102 + j);
		}
	}

	{
		// float float 
		auto floatVal1 = VecXX::scalar(11.1);
		auto floatVal2 = VecXX::scalar(13.1);

		auto res = ApplyBinaryOperation<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(floatVal1, floatVal2);
		EXPECT_TRUE(res.isScalar());
		EXPECT_NEAR(res.getScalarValue(), 24.2, testEpsilon);

	}


}



template <typename INS_VEC>
void evalApplyBinaryOperation1VV(size_t startLen, size_t endLen)
{
	const double testEpsilon = 0.000001;
	for (size_t SZ = startLen; SZ <= endLen; SZ++)
	{

		std::vector<double>  v(SZ, 0.0);
		for (int i = 0; i < SZ; i++) { v[i] = i; }// (0 == i % 2) ? -i : i;

		VecXX Lhs(v);

		VecXX Rhs = Lhs * 3.2 - 1.0;

		auto res = ApplyBinaryOperation1<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(Lhs, Rhs, DR_CUBED::_plus<VecXX::INS >());

		int i = 0;
		for (double x : res)
		{
			EXPECT_NEAR(x, v[i] * (1 + 3.2) - 1.0, testEpsilon);
			i++;
		}
	}
}




TEST(TestBinaryAndUnitaryFunctions, ApplyBinaryOperation1VV)
{
	const double testEpsilon = 0.000001;
	{
		VecXX  LHS = VecXX::scalar(1.1);
		VecXX  RHS = VecXX::scalar(3.1);

		auto res = ApplyBinaryOperation1(LHS, RHS, DR_CUBED::_plus < VecXX::INS>());
		EXPECT_TRUE(res.isScalar());
		EXPECT_NEAR(res.getScalarValue(), 4.2, testEpsilon);
		////////////////// both inputs scalar /////////////

		auto  scalarLHS = VecXX::scalar(1.1);
		auto  scalarRHS = VecXX::scalar(3.1);

		res = ApplyBinaryOperation1<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(scalarLHS, scalarRHS, DR_CUBED::_plus < VecXX::INS>());
		EXPECT_TRUE(res.isScalar());
		EXPECT_NEAR(res.getScalarValue(), 4.2, testEpsilon);
		////////////////// both inputs scalar ////////////


		for (int j = 1; j < 102; ++j)
		{
			evalApplyBinaryOperation1VV< VecXX::INS>(1, 102 + j);
		}
	}
}




//TO DO make these run over variable length arrays on test
TEST(TestBinaryAndUnitaryFunctions, ApplyBinaryOperationsMutating)
{
	const double testEpsilon = 0.000001;

	int SZ = 13;
	std::vector<double>  v(SZ, 0.0);
	for (int i = 0; i < SZ; i++) { v[i] = i; }
	VecXX smallVec(v);
	int i = 0;

	/******* this is OK */
	//float+=vec ??
	{
		// float float 
		auto floatVal1 = VecXX::scalar(11.1);
		//auto floatVal2 = VecXX::scalar(13.1);
		//vec += flt , when vec
		auto res = ApplyBinaryOperationM<VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(floatVal1, smallVec);
		EXPECT_TRUE(!res.isScalar());
		i = 0;
		for (double x : res)
		{
			EXPECT_NEAR(x, v[i] + 11.1, testEpsilon);
			i++;
		}

	
		VecXX scalarVec = VecXX::scalar(13.1);

		res = ApplyBinaryOperationM<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(floatVal1, scalarVec);
		EXPECT_TRUE(res.isScalar());  //not a scalar a vec with one element
		EXPECT_TRUE(res.getScalarValue() ==13.1+11.1);

	}

	
	/////////////////////////////////
	//vec += float
	{
		// float float 
		auto floatVal1 = VecXX::scalar(11.1);
		//vec += flt , when vec
		VecXX smallVec1(v);
		auto res = ApplyBinaryOperationM<VecXX::INS, DR_CUBED::_plus<VecXX::INS > >( smallVec1, floatVal1);
		EXPECT_TRUE(!res.isScalar());
		i = 0;
		for (double x : res)
		{
			EXPECT_NEAR(x, v[i] + 11.1, testEpsilon);
			i++;
		}

		VecXX scalarVec = VecXX::scalar(13.1);
		res = ApplyBinaryOperationM<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(scalarVec,floatVal1);
		EXPECT_TRUE(res.isScalar());  //not a scalar a vec with one element
		EXPECT_TRUE(res.getScalarValue() == 13.1 + 11.1);

	}

	
	//vec += vec
	{
		// float float 
		auto floatVal1 = VecXX::scalar(11.1);
		VecXX smallVec1(v);
		VecXX scalarVec = floatVal1; // this is a VecXX initialised as scalar  mode

		//vec += flt , when flt is a single element vec
		auto res = ApplyBinaryOperationM<VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(smallVec1, scalarVec);
		EXPECT_TRUE(!res.isScalar());
		i = 0;
		for (double x : res)
		{
			EXPECT_NEAR(x, v[i] + 11.1, testEpsilon);
			i++;
		}
		
		VecXX smallVec2(v);
		VecXX anotherVec = smallVec2 * 1.1;
		res = ApplyBinaryOperationM<  VecXX::INS, DR_CUBED::_plus<VecXX::INS > >(smallVec2, anotherVec);
		EXPECT_TRUE(!res.isScalar());  //not a scalar a vec with one element
		int i = 0;
	    for (double x : res)
		{
			double expected = v[i];
			expected *= 2.1;
			EXPECT_NEAR(x, expected, testEpsilon);
			i++;
		}

	}



}


void testApplySparse(int start, int SZ, int mod)
{

	const double testEpsilon = 0.000001;

	SZ = 130;
	std::vector<double>  v(SZ, 0.0);
	for (int i = 0; i < SZ; i++) { v[i] = i; }
	VecXX smallVec(v);


	auto sinOper = [](auto x) {return sin( VecXX::INS(0.25) * x); };

	auto twenty = VecXX::scalar(20.);
	auto sparseFunc = [=](auto x) { return ((x / twenty - floor(x / twenty)) < 0.000001); };

	auto bigOPer = [=](auto x) {return twenty; };

	auto res = ApplyUnitaryOperation1(smallVec, sinOper);

	int i = 0;
	for (double x : res)
	{
		double expected = v[i];
		expected = sin(0.25 * i);
		EXPECT_NEAR(x, expected, testEpsilon);
		i++;
	}

	
	ApplySparseUnitaryOperationU(smallVec, res, bigOPer, sparseFunc);


	i = 0;
	for (double x : res)
	{
		if (i % 20 == 0)
		{
			EXPECT_NEAR(x, 20.0, testEpsilon);
		}
		else
		{
			double expected = v[i];
			expected = sin(0.25 * i);
			EXPECT_NEAR(x, expected, testEpsilon);
		}
		i++;
	}
}


TEST(TestBinaryAndUnitaryFunctions, ApplySparseUnitaryOperation)
{
	//TO DO LOOKS LIKE ITS NOT inmplemented as remembered, is the sparseness down to the user
	// to use an iff or select

	for (int mod = 1; mod < 35; ++mod)
	{
		for (int j = 1; j < 102; ++j)
		{
			testApplySparse(1, 102 + j, mod);
		}
	}
}

////////////////////////////////////////////////