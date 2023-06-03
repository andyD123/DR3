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
#include "testNamespace.h"
#include "dr3TestUtil.h"


#include <numeric>

TEST(TestSpan, ApplyUnitaryOperation)
{
	
	for (int SZ = 3; SZ < 125; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;
		
		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		SpanXX test(testV.begin(), SZ);

		auto mySquareItLambda = [](auto x) {return x * x;  };

		auto testSquare = test;

		std::vector<double> std_vec_in = test;

		ApplyUnitaryOperation(testSquare, mySquareItLambda);


		for (int i = 0; i < test.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(testSquare[i], std_vec_in[i] * std_vec_in[i]);
		}
		
	}

}

TEST(TestSpan, ApplyUnitaryOperation_shifting_start)
{
	
	for (int start = 0; start < 125; start++)
	{

		for (int SZ = 1; SZ < 125; SZ++)
		{
			std::vector<Numeric>  v(SZ+200, asNumber(6.66));
			int i = 0;

			for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
			const VecXX testV(v);

			SpanXX test(testV.begin() + start, SZ);

			auto mySquareItLambda = [](auto x) {return x * x;  };

			auto testSquare = test;

			std::vector<double> std_vec_in = test;

			ApplyUnitaryOperation(testSquare, mySquareItLambda);

			for (int i = 0; i < test.size(); ++i)
			{
				EXPECT_NUMERIC_EQ(testSquare[i], std_vec_in[i] * std_vec_in[i]);
			}

		}
	}

}

TEST(TestSpan, transform)
{

	for (int SZ = 3; SZ < 125; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		SpanXX test(testV.begin(), SZ);

		auto mySquareItLambda = [](auto x) {return x * x;  };

		std::vector<double> std_vec_in = test;

		VecXX outVec(v);
		SpanXX result(outVec.begin(), SZ);

		transform( mySquareItLambda, test, result);


		for (int i = 0; i < test.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(result[i], std_vec_in[i] * std_vec_in[i]);
		}

	}

}

TEST(TestSpan, transforM)
{

	for (int SZ = 3; SZ < 125; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		SpanXX test(testV.begin(), SZ);

		auto mySquareItLambda = [](auto x) {return x * x;  };

		std::vector<double> std_vec_in = test;


		transformM(mySquareItLambda, test);


		for (int i = 0; i < test.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(test[i], std_vec_in[i] * std_vec_in[i]);
		}

	}


}

TEST(TestSpan, transforM_shifting_start)
{

	for (int start = 0; start < 125; start++)
	{

		for (int SZ = 1; SZ < 125; SZ++)
		{
			std::vector<Numeric>  v(SZ + 200, asNumber(6.66));
			int i = 0;

			for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
			const VecXX testV(v);

			SpanXX test(testV.begin() + start, SZ);

			auto mySquareItLambda = [](auto x) {return x * x;  };

			auto testSquare = test;

			std::vector<double> std_vec_in = test;

			transformM(mySquareItLambda, testSquare);

			for (int i = 0; i < test.size(); ++i)
			{
				EXPECT_NUMERIC_EQ(testSquare[i], std_vec_in[i] * std_vec_in[i]);
			}
		}
	}

}

TEST(TestSpan, filter)
{

	for (int SZ = 3; SZ < 125; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		SpanXX testSpan(testV.begin(), SZ);

		std::vector<double> expected;

		for (const auto& X : testSpan)
		{
			if ((2.0 * std::floor(X * 0.5) - X) > -0.0001)
			{
				expected.emplace_back(X);
			}

		}

		auto isEven = [](auto x) { return !((x - VecXX::INS((2.0) * floor(x * VecXX::INS(0.5)))) > VecXX::INS(0.)); };
		

		auto evenView = filter(isEven, testSpan);

		std::vector<double> dbgVw = evenView;


		for (int i = 0; i < evenView.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(expected[i], evenView[i]);
		}

	}


}

TEST(TestSpan, filter_shifting_start)
{
	for (int start = 0; start < 37; start++)
	{
		for (int SZ = 1; SZ < 125; SZ++)
		{
			std::vector<Numeric>  v(600, asNumber(6.66));
			int i = 0;

			for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
			const VecXX testV(v);

			SpanXX testSpan(testV.begin() + start, SZ);

			std::vector<double> expected;

			for (const auto& X : testSpan)
			{
				if ((2.0 * std::floor(X * 0.5) - X) > -0.0001)
				{
					expected.emplace_back(X);
				}

			}

			auto isEven = [](auto x) { return !((x - VecXX::INS((2.0) * floor(x * VecXX::INS(0.5)))) > VecXX::INS(0.)); };


			auto evenView = filter(isEven, testSpan);

			std::vector<double> dbgVw = evenView;


			for (int i = 0; i < evenView.size(); ++i)
			{
				EXPECT_NUMERIC_EQ(expected[i], evenView[i]);
			}
		}
	}

}

TEST(TestSpan, reduce)
{

	for (int SZ = 3; SZ < 125; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		SpanXX testSpan(testV.begin(), SZ);

		auto SUM = [](auto x, auto y) { return x + y; };

		auto result = reduce(testSpan, SUM);
		double expected = 0.;

		for (const auto& X : testSpan)
		{
			expected += X;
		}
		
		EXPECT_NUMERIC_EQ(expected, result);	

	}

}

TEST(TestSpan, reduce_shifting_start)
{
	for (int start = 0; start < 37; start++)
	{

		for (int SZ = 3; SZ < 125; SZ++)
		{
			std::vector<Numeric>  v(SZ+300, asNumber(6.66));
			int i = 0;

			for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
			const VecXX testV(v);

			SpanXX testSpan(testV.begin() + start, SZ);

			auto SUM = [](auto x, auto y) { return x + y; };

			auto result = reduce(testSpan, SUM);
			double expected = 0.;

			for (const auto& X : testSpan)
			{
				expected += X;
			}

			EXPECT_NUMERIC_EQ(expected, result);

		}
	}

}




TEST(TestSpan, transform_reduce)
{

	for (int SZ = 3; SZ < 125; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		SpanXX testSpan(testV.begin(), SZ);

		auto SUM = [](auto x, auto y) { return x + y; };
		auto SQR = [](auto x) { return x * x; };

		auto result = transformReduce(testSpan, SQR, SUM);
		double expected = 0.;

		for (const auto& X : testSpan)
		{
			expected += X*X;
		}

		EXPECT_NUMERIC_EQ(expected, result);

	}

}



TEST(TestSpan, transform_reduce_shifting_start)
{

	for (int start = 0; start < 37; start++)
	{

		for (int SZ = 3; SZ < 125; SZ++)
		{
			std::vector<Numeric>  v(SZ + 300, asNumber(6.66));
			int i = 0;

			for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
			const VecXX testV(v);

			SpanXX testSpan(testV.begin(), SZ);

			auto SUM = [](auto x, auto y) { return x + y; };
			auto SQR = [](auto x) { return x * x; };

			auto result = transformReduce(testSpan, SQR, SUM);
			double expected = 0.;

			for (const auto& X : testSpan)
			{
				expected += X * X;
			}

			EXPECT_NUMERIC_EQ(expected, result);
		}

	}

}