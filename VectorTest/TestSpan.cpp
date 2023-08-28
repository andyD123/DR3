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

		std::vector<Numeric> std_vec_in = test;

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

			std::vector<Numeric> std_vec_in = test;

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

		std::vector<Numeric> std_vec_in = test;

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

		std::vector<Numeric> std_vec_in = test;


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

			std::vector<Numeric> std_vec_in = test;

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

		std::vector<Numeric> expected;

		for (const auto& X : testSpan)
		{
			if ((2.0 * std::floor(X * 0.5) - X) > -0.0001)
			{
				expected.emplace_back(X);
			}

		}

		auto isEven = [](auto x) { return !((x - VecXX::INS((2.0) * floor(x * VecXX::INS(0.5)))) > VecXX::INS(0.)); };
		

		auto evenView = filter(isEven, testSpan);

		std::vector<Numeric> dbgVw = evenView;


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

			std::vector<Numeric> expected;

			std::vector<Numeric> inspan = testSpan;

			expected = testSpan;
			expected.clear();

			auto isEvenScalar = [](auto x) { return !(((x - 2.0 * floor(x * 0.5)) > 0.)); };

			for (const auto& X : testSpan)
			{
				if ( isEvenScalar(X) )   
				{
					expected.emplace_back(X);
				}

			}

			auto isEven = [](auto x) { return !((x - VecXX::INS((2.0) * floor(x * VecXX::INS(0.5)))) > VecXX::INS(0.)); };


			auto evenView = filter(isEven, testSpan);

			std::vector<Numeric> dbgVw = evenView;


			for (int i = 0; i < evenView.size(); ++i)
			{
				EXPECT_NUMERIC_EQ(expected[i], evenView[i]);
			}
		}
	}

}


TEST(TestSpan, reduce)
{
	//this appears to work on skylake, but loads are not unaligned

	for (int SZ = 3; SZ < 125; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		SpanXX testSpan(testV.begin(), SZ);

		auto SUM = [](auto x, auto y) { return x + y; };

		auto result = reduce(testSpan, SUM);
		Numeric expected = 0.;

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
			Numeric expected = 0.;

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
		Numeric expected = 0.;

		for (const auto& X : testSpan)
		{
			expected += X*X;
		}

		EXPECT_NUMERIC_EQ(expected, result);

	}

}



TEST(TestSpan, transform_reduce_binary)
{
	
	for (int SZ = 3; SZ < 125; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;
		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		SpanXX testSpan(testV.begin(), SZ);

		SpanXX testSpan2(testV.begin(), SZ);

		auto SUM = [](auto x, auto y) { return x + y; };
		auto MULT = [](auto x, auto y) { return x * y; };

		auto result = transformReduce(testSpan, testSpan2, MULT, SUM);
		Numeric expected = 0.;

		for (int i = 0; i < SZ; ++i)
		{
			expected += testV[i] * testV[i];
		}

		EXPECT_NUMERIC_EQ(expected, result);

	}

	// with offsets
	for (int offset = 0; offset < 20; ++offset)
	{
		for (int SZ = 3; SZ < 125; SZ++)
		{
			std::vector<Numeric>  v(SZ+30, asNumber(6.66));
			int i = 0;
			for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
			const VecXX testV(v);

			SpanXX testSpan(testV.begin()+ offset, SZ);

			SpanXX testSpan2(testV.begin(), SZ);

			auto SUM = [](auto x, auto y) { return x + y; };
			auto MULT = [](auto x, auto y) { return x * y; };

			auto result = transformReduce(testSpan, testSpan2, MULT, SUM);
			Numeric expected = 0.;

			for (int i = 0; i < SZ; ++i)
			{
				expected += testV[i+ offset] * testV[i];
			}

			EXPECT_NUMERIC_EQ(expected, result);

		}
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
			Numeric expected = 0.;

			for (const auto& X : testSpan)
			{
				expected += X * X;
			}

			EXPECT_NUMERIC_EQ(expected, result);
		}

	}

}



TEST(TestSpan, strided_transform)
{

	{
		int stride = 8;
		int SZ = 320;
		std::vector<Numeric>  v(SZ + 300, asNumber(0.0));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		StrdSpanXX test(testV.begin(), SZ, stride);

		auto mySquareItLambda = [](auto x) {return x *x;  };

		auto sameVal = [](auto x) {return x;  };


		auto testSquare = test;

		std::vector<Numeric> std_vec_in = test;

		
		VecXX outVec(v);
		outVec *= 0.;
		SpanXX result(outVec.begin(), SZ);

		transform(sameVal, test, result);

		std::vector<Numeric> dbug = result;

		std::vector<Numeric> expected;

		for (const auto& x : test)
		{
			int testVal = static_cast<int>(x);
			if (testVal % stride == 0)
				expected.emplace_back(x);
		}

	

		for (int i = 0; i < expected.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(result[i], expected[i] );
		}
	}

}


TEST(TestSpan, strided_transform_withOffset)
{

	{
		int offset = 1;

		int stride = 8;
		int SZ = 320;
		std::vector<Numeric>  v(SZ + 300, asNumber(0.0));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		StrdSpanXX test(testV.begin()+ offset, SZ, stride);

		auto mySquareItLambda = [](auto x) {return x * x;  };

		auto sameVal = [](auto x) {return x;  };


		auto testSquare = test;

		std::vector<Numeric> std_vec_in = test;


		VecXX outVec(v);
		outVec *= 0.;
		SpanXX result(outVec.begin(), SZ);

		transform(sameVal, test, result);

		std::vector<Numeric> dbug = result;

		std::vector<Numeric> expected;

		for (const auto& x : test)
		{
			int testVal = static_cast<int>(x);
			if (testVal % stride == offset % stride)
				expected.emplace_back(x);
		}



		for (int i = 0; i < expected.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(result[i], expected[i]);
		}


	}

}




TEST(TestSpan, strided_transform_2)
{

	for (int SZ = 3; SZ < 125; SZ++)
	{
		int stride = 8;

		std::vector<Numeric>  v(SZ + 300, asNumber(0.0));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		StrdSpanXX test(testV.begin(), SZ, stride);

		auto mySquareItLambda = [](auto x) {return x * x;  };

		auto sameVal = [](auto x) {return x;  };


		auto testSquare = test;

		std::vector<Numeric> std_vec_in = test;


		VecXX outVec(v);
		outVec *= 0.;
		SpanXX result(outVec.begin(), SZ);

		transform(sameVal, test, result);

		std::vector<Numeric> dbug = result;

		std::vector<Numeric> expected;

		for (const auto& x : test)
		{
			int testVal = static_cast<int>(x);
			if (testVal % stride == 0)
				expected.emplace_back(x);
		}


		for (int i = 0; i < expected.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(result[i], expected[i]);
		}
	}

}



TEST(TestSpan, strided_transform_withOffset_2)
{

	int stride = 8;

	for (int offset = 0; offset < 66; offset++)
	{

		for (int SZ = 3; SZ < 125; SZ++)
		{

			std::vector<Numeric>  v(SZ + 300, asNumber(0.0));
			int i = 0;

			for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
			const VecXX testV(v);

			StrdSpanXX test(testV.begin() + offset, SZ, stride);

			auto mySquareItLambda = [](auto x) {return x * x;  };

			auto sameVal = [](auto x) {return x;  };


			auto testSquare = test;

			std::vector<Numeric> std_vec_in = test;


			VecXX outVec(v);
			outVec *= 0.;
			SpanXX result(outVec.begin(), SZ);

			transform(sameVal, test, result);

			std::vector<Numeric> dbug = result;

			std::vector<Numeric> expected;

			for (const auto& x : test)
			{
				int testVal = static_cast<int>(x);
				if (testVal % stride == offset % stride)
					expected.emplace_back(x);
			}

			for (int i = 0; i < expected.size(); ++i)
			{
				EXPECT_NUMERIC_EQ(result[i], expected[i]);
			}


		}
	}

}



TEST(TestSpan, strided_filter)
{

	int stride = 8;


	int offset = 0;


	//span size
	for (int SZ = 3; SZ < 125; SZ++)
	{


			std::vector<Numeric>  v(160 * stride, asNumber(0.));
			int i = 0;

			for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
			const VecXX testV(v);

			StrdSpanXX testSpan(testV.begin() + offset, SZ, stride);

			std::vector<Numeric> expected;

			auto  SZ_MAX = std::min(v.size(), static_cast<size_t>(SZ * stride));

			for (size_t jj = offset; jj < (SZ+ offset); jj += stride)
			{
				const auto& X = testV[jj];

				if ((2.0 * std::floor(X * 0.5) - X) > -0.0001)
				{
					expected.emplace_back(X);
				}
			}

			auto isEven = [](auto x) { return !((x - VecXX::INS((asNumber(2.0)) * floor(x * VecXX::INS(asNumber(0.5))))) > VecXX::INS(asNumber(0.))); };
			auto evenView = filter(isEven, testSpan);

			std::vector<Numeric> dbgSpan = testSpan;
			std::vector<Numeric> dbgVw = evenView;


			for (int i = 0; i < evenView.size(); ++i)
			{
				EXPECT_NUMERIC_EQ(expected[i], evenView[i]);
			}

	}


}


TEST(TestSpan, strided_filter_vary_offsets)
{

	int stride = 8; // stride  is multiple of width of instruction  set  i/e we assume padded and stride =  INT * width

	for (int offset = 0; offset < 45; ++offset) //offset
	{
		for (int SZ =1; SZ < 125; SZ++) //span size 
		{

			std::vector<Numeric>  v(160 * stride, asNumber(0.));
			int i = 0;

			for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
			const VecXX testV(v);

			StrdSpanXX testSpan(testV.begin() + offset, SZ, stride);

			std::vector<Numeric> expected;

			auto  SZ_MAX = std::min(v.size(), static_cast<size_t>(SZ * stride));

			for (size_t jj = offset; jj < (SZ + offset); jj += stride)
			{
				const auto& X = testV[jj];

				if ((2.0 * std::floor(X * 0.5) - X) > -0.0001)
				{
					expected.emplace_back(X);
				}
			}

			auto isEven = [](auto x) { return !((x - VecXX::INS((asNumber(2.0)) * floor(x * VecXX::INS(asNumber(0.5))))) > VecXX::INS(asNumber(0.))); };
			auto evenView = filter(isEven, testSpan);

			std::vector<Numeric> dbgSpan = testSpan;
			std::vector<Numeric> dbgVw = evenView;


			EXPECT_NUMERIC_EQ(static_cast<int>(expected.size()), evenView.size());

			for (int i = 0; i < evenView.size(); ++i)
			{
				EXPECT_NUMERIC_EQ(expected[i], evenView[i]);
			}
			
		}

	}

}





TEST(TestSpan, strided_reduce)
{
	int offset = 0;

	
	for (int SZ = 3; SZ < 125; SZ++)
	{
		int stride = 8;
		
		std::vector<Numeric>  v(SZ + 300, asNumber(0.0));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		StrdSpanXX test(testV.begin(), SZ, stride);
		
		auto testSquare = test;

		std::vector<Numeric> std_vec_in = test;


		auto SUM = [](auto x, auto y) { return x + y; };
		

		auto result = reduce(test,  SUM);
		Numeric expected = 0.;

		for (int i = 0; i < SZ;i+= stride)
		{
			expected += v[i + offset];
		}

		EXPECT_NUMERIC_EQ(expected, result);

	}

}




TEST(TestSpan, strided_reduce_large)
{
	int offset = 0;

	int SZ = 4003;
	{
		int stride = 8;
		
		std::vector<Numeric>  v(SZ + 300, asNumber(0.0));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		StrdSpanXX test(testV.begin(), SZ, stride);

		auto SUM = [](auto x, auto y) { return x + y; };
		auto result = reduce(test, SUM);
		Numeric expected = 0.;

		for (int i = 0; i < SZ; i += stride)
		{
			expected += v[i + offset];
		}

		EXPECT_NUMERIC_EQ(expected, result);
	}
}



void testStridedSpanReduce(int offset, int startSZ, int endSZ)
{
	for (int SZ = startSZ; SZ <= endSZ; SZ++)
	{
		int stride = 8;

		std::vector<Numeric>  v(SZ + 300, asNumber(0.0));
		int i = 0;

		for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
		const VecXX testV(v);

		StrdSpanXX test(testV.begin() + offset, SZ, stride);

		auto SUM = [](auto x, auto y) { return x + y; };
		auto result = reduce(test, SUM);
		Numeric expected = 0.;

		for (int i = 0; i < SZ; i += stride)
		{
			expected += v[i + offset];
		}

		EXPECT_NUMERIC_EQ(expected, result);
	}

}


TEST(TestSpan, strided_reduce_large_2)
{
	
	int offset = 1;

	for (int offset = 0; offset < 31; offset++)
	{

		testStridedSpanReduce( offset, 3, 33);

		testStridedSpanReduce( offset, 250, 303 );

		testStridedSpanReduce( offset,600 , 649);
	}

}



void testTransformReduceStrided(int startSz, int maxSz, int offsetStart, int offsetEnd)
{

	for (int offset = offsetEnd; offset < offsetEnd; offset++)
	{


		for (int SZ = startSz; SZ < maxSz; SZ++)
		{
			int stride = 8;

			std::vector<Numeric>  v(SZ + 300, asNumber(0.0));
			int i = 0;

			for (auto& x : v) { x = asNumber(0.0 + i); ++i; }
			const VecXX testV(v);

			StrdSpanXX test(testV.begin(), SZ, stride);

			auto testSquare = test;

			std::vector<Numeric> std_vec_in = test;

			auto SQR = [](auto x) { return x * x; };
			auto SUM = [](auto x, auto y) { return x + y; };


			auto result = transformReduce(test, SUM, SQR);
			Numeric expected = 0.;

			for (int i = 0; i < SZ; i += stride)
			{
				expected += v[i + offset] * v[i + offset];
			}

			EXPECT_NUMERIC_EQ(expected, result);

		}
	}

}

TEST(TestSpan, strided_transform_reduce)
{
	//int offset = 0;

	testTransformReduceStrided(3, 125, 0, 1);

	//small but with offset range greater than stride
	testTransformReduceStrided(3, 125, 0, 31);


	testTransformReduceStrided(250, 300, 0, 31);


	testTransformReduceStrided(600, 650, 0, 31);


}



TEST(TestSpan, binary_transform_simple)
{
	
	for (int SZ = 3; SZ < 125; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;

		for (auto& x : v) { x = asNumber(1.0 + i); ++i; }
		const VecXX testV(v);

		VecXX testY =  testV * 2.0;
		testY += 1.0;

		SpanXX spn1(testV.begin(), SZ);

		SpanXX spn2(testY.begin(), SZ);

		VecXX testV2(0.0, SZ);

		SpanXX outSpan(testV2.begin(), SZ);

		auto DIV = [](auto x, auto y) {return x / y;  };

		std::vector<Numeric> std_vec_in = spn1;
		std::vector<Numeric> std_vec_in2 = spn2;

		transform(DIV, spn1, spn2, outSpan);

		std::vector<Numeric> outDbg = outSpan;

		for (int i = 1; i < spn2.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(outSpan[i], spn1[i] / spn2[i]);
		}

	}

}



TEST(TestSpan, binary_transform_simple_with_offset)
{
	

	for (int SZ = 3; SZ < 125; SZ++)
	{
		for (int offset = 0; offset < SZ; offset++)
		{

			std::vector<Numeric>  v(SZ * 2, asNumber(6.66));
			int i = 0;

			for (auto& x : v) { x = asNumber(1.0 + i); ++i; }
			const VecXX testV(v);

			VecXX testY = testV * 2.0;
			testY += 1.0;

			SpanXX spn1(testV.begin() + offset, SZ);

			SpanXX spn2(testY.begin() + offset, SZ);

			VecXX testV2(0.0, SZ * 2);

			SpanXX outSpan(testV2.begin(), SZ);

			auto DIV = [](auto x, auto y) {return x / y;  };

			//std::vector<double> std_vec_in = spn1;
			//std::vector<double> std_vec_in2 = spn2;

			transform(DIV, spn1, spn2, outSpan);

			//std::vector<double> outDbg = outSpan;

			for (int i = 1; i < spn2.size(); ++i)
			{
				EXPECT_NUMERIC_EQ(outSpan[i ], spn1[i ] / spn2[i ]);
			}
		}

	}


	for (int SZ = 3; SZ < 125; SZ++)
	{
		//int offset = 0;
		for (int offset = 0; offset < SZ - 5; offset++)
		{
			//int offset2 = 0;
			for (int offset2 = 0; offset2 < 5; offset2++)
			{

				std::vector<Numeric>  v(SZ * 2, asNumber(6.66));
				int i = 0;

				for (auto& x : v) { x = asNumber(1.0 + i); ++i; }
				const VecXX testV(v);

				VecXX testY = testV * 2.0;
				testY += 1.0;

				SpanXX spn1(testV.begin() + offset + offset2, SZ);

				SpanXX spn2(testY.begin() + offset, SZ);

				VecXX testV2(0.0, SZ * 2);

				SpanXX outSpan(testV2.begin(), SZ);

				auto DIV = [](auto x, auto y) {return x / y;  };

				std::vector<Numeric> std_vec_in = spn1;
				std::vector<Numeric> std_vec_in2 = spn2;

				transform(DIV, spn1, spn2, outSpan);

				std::vector<Numeric> outDbg = outSpan;

				for (int i = 1; i < spn2.size(); ++i)
				{
					EXPECT_NUMERIC_EQ(outSpan[i], spn1[i ] / spn2[i]);
				}
			}
		}
		
	}

}





TEST(TestSpan, binary_transform_simple_with_vector)
{

	for (int SZ = 3; SZ < 125; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;

		for (auto& x : v) { x = asNumber(1.0 + i); ++i; }
		const VecXX testV(v);

		VecXX testY = testV * 2.0;
		testY += 1.0;

		SpanXX spn1(testV.begin(), SZ);

		SpanXX spn2(testY.begin(), SZ);

		VecXX testV2(0.0, SZ);

		SpanXX outSpan(testV2.begin(), SZ);

		auto DIV = [](auto x, auto y) {return x / y;  };

		std::vector<Numeric> std_vec_in = spn1;
		std::vector<Numeric> std_vec_in2 = spn2;

		//transform(DIV, spn1, spn2, outSpan);
		transform(DIV, spn1, testY, outSpan);

		std::vector<Numeric> outDbg = outSpan;

		for (int i = 1; i < spn2.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(outSpan[i], spn1[i] / testY[i]);
		}

	}


	for (int SZ = 3; SZ < 125; SZ++)
	{
		std::vector<Numeric>  v(SZ, asNumber(6.66));
		int i = 0;

		for (auto& x : v) { x = asNumber(1.0 + i); ++i; }
		const VecXX testV(v);

		VecXX testY = testV * 2.0;
		testY += 1.0;

		SpanXX spn2(testY.begin(), SZ);

		VecXX testV2(0.0, SZ);

		SpanXX outSpan(testV2.begin(), SZ);

		auto DIV = [](auto x, auto y) {return x / y;  };

		std::vector<Numeric> std_vec_in = testY;
		std::vector<Numeric> std_vec_in2 = spn2;

		transform(DIV, testY, spn2, outSpan);
		
		std::vector<Numeric> outDbg = outSpan;

		for (int i = 1; i < spn2.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(outSpan[i], testY[i] / spn2[i]);
		}

		// check that its also  written to testV2
		for (int i = 1; i < spn2.size(); ++i)
		{
			EXPECT_NUMERIC_EQ(testV2[i], testY[i] / spn2[i]);
		}

	}


}






TEST(TestSpan, basic_2D_SPAN)
{

	//create a vector of 1D data 160 in length
	VecXX owningVec(0.1, 16 * 10);

	auto val = 0.0;
	for (auto& x : owningVec)
	{
		x = val;
		val++;
	}



	using MAT1 = Layout2D<Numeric, 8, 0>;

	//auto pDat =
	MDSpan<Numeric, MAT1> mat(owningVec.data(), 10, 10);

	/*
	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			std::cout << mat(i, j) << ",";
		}
		std::cout << "\n";
	}

	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			mat(i, j) = i * 100 + j;
		}

	}

	std::cout << "\n";
	std::cout << "\n";

	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			std::cout << mat(i, j) << ",";
		}
		std::cout << "\n";
	}
	*/
}


