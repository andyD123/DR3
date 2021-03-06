// accumulateExample.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <algorithm>
#include <random>
#include <numeric>
#include <iterator>
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>  



#include "../Vectorisation/VecX/dr3.h"

#include "norm.h"


/*
Uncomment one of the Using namespace lines below to select the instruction set that you wish to run  
Those ending in F have float type as underlying, those ending with D have a double.

The project is set to compile using the AVX512  enhanced instruction set. The namespace selection 
choses the type of the intrinsics that are used to instantiate lambdas.

If your hardware does not support AVX512 chose the next level down AVX2 and avoid using namespaces 
DRC::VecD8D or DRC::VecF8F which will cause generation of code with instructions that your computer doesn't support. 

check device manager/processor to determine what processor you have and check against web site 
https://ark.intel.com/content/www/us/en/ark/products/123550/intel-xeon-silver-4114-processor-13-75m-cache-2-20-ghz.html
or
https://www.intel.com/content/www/us/en/products/details/processors/xeon/scalable.html

//Next uncomment example that you wish to run.
You can change projects setting for the compiler debug/release  are for VS2019, clang  and ICC2022

precision of results compare for float types is incorrect.

*/

//pick an instruction set for intrinsics by selecting a name space

//Using namespace DRC::VecDb;
//using namespace DRC::VecD2D;  //sse2   double
using namespace DRC::VecD4D;	//avx2   double
//using namespace DRC::VecF8F;	// avx2  float
//using namespace DRC::VecD8D;  //avx512 double
//using namespace DRC::VecF16F; //avx512   float




const double billion = 1000000000.0;


template<typename T>
bool vectorsEqual(const std::vector<T>& C1, const std::vector<T>& C2, const std::vector<T>& C3)
{
	bool  testOK = true;
	const double ERR = 1e-13; //for examples
	if (C1.size() != C2.size())
	{
		return false;
	}

	if (C3.size() != C2.size())
	{
		return false;
	}

	for (size_t i = 0; i < C3.size(); i++)
	{
		auto err1 = fabs((C1[i] - C2[i]) / (C2[i] + C1[i]));
		auto err2 = fabs((C1[i] - C3[i]) / (C1[i] + C3[i]));

		if ((err1 > ERR) || (err2 > ERR))
		{
			testOK = false;
			std::cout << "\n err diff@ " << i << " err1 =" << err1 << ", err2 = " << err2 << "\n";
			break;
		}
	}

	return testOK;

}


template<typename T>
bool vectorsEqualD(const std::vector<T>& C1, const std::vector<T>& C2, const std::vector<T>& C3, const std::vector<T>& input)
{
	bool  testOK = true;
	const double ERR = 1e-14;
	if (C1.size() != C2.size())
	{
		return false;
	}

	if (C3.size() != C2.size())
	{
		return false;
	}

	for (int i = 0; i < C3.size(); i++)
	{
		auto err1 = fabs((C1[i] - C2[i]) / (C2[i] + C1[i]));
		auto err2 = fabs((C1[i] - C3[i]) / (C1[i] + C3[i]));

		if ((err1 > ERR) || (err2 > ERR))
		{
			testOK = false;
			std::cout << "\n err diff@ " << i << " err1 =" << err1 << ", err2 = " << err2 << "\n";
			std::cout << "\n val @ " << i << " C1[i] =" << C1[i] << ", C3[i] = " << C3[i] << "input val=" << input[i] <<"\n";
			break;
		}
	}

	return testOK;

}


auto getRandomShuffledVector(int SZ)
{
	using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
	std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
	for (int i = 0; i < SZ; i++) { v[i] += FloatType(SZ / 2) + i; }//{ v[i] += FloatType(SZ/2) - i; }
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(begin(v), end(v), g);
	return v;
}


auto numOps = [](int TEST_LOOP_SZ, int SZ) { return  static_cast<int>(double(TEST_LOOP_SZ) * double(SZ)); };


double getnull(double)
{
	return 0.0;
}


//example functions fwd decl

void	testMemCpy2();
void    doMax();
void 	doSum(); 
void    doInnerProd();
void	doTransform();
void 	doSumSqrs();
void    khanAccumulation();
void	testBinarySelection();
void	testBinarySelection1();
void	testBinarySelection2();
void    testBinarySelection3();
void	doCountIf();







int main()
{

 //Uncomment  a function to play with

	//	testMemCpy2(); 
	 doMax();
	//	doSum(); // stl slower with intel  stl slower
	//    doInnerProd();
	//	doTransform();
	// 	doSumSqrs();
	//  khanAccumulation();
	//	testBinarySelection(); //select between constants
	//	testBinarySelection1();//select between light functions
	//	testBinarySelection2(); //medium
	//    testBinarySelection3(); //heavy	
	//	doCountIf();



	return 0;
}




void testMemCpy()
{

	using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
	std::vector<FloatType>  v(400, VecXX::SCALA_TYPE(0.5));
	std::vector<FloatType>  c(400, VecXX::SCALA_TYPE(0.5));
	VecXX test(v);
	VecXX test2(v);
	test2 *= 0.0;

	{
		auto startTme = std::chrono::high_resolution_clock::now();


		for (int i = 0; i < 100000000; ++i)
		{
			std::memcpy(&c[0], &v[0], 400 * sizeof(double));
		}

		auto endTime = std::chrono::high_resolution_clock::now();
		auto runtime = endTime - startTme;
		std::cout << "MemCpy run time mem 10 bill doubles copy  320GB=" << runtime.count() / 1000000000.0 << " seconds  = " << 320. / (runtime.count() / 1000000000.0) << " GB per sec" << std::endl;

	}

	{
		auto startTme = std::chrono::high_resolution_clock::now();
		auto cpyLambda = [](const auto& rhs) { return rhs; };

		for (int i = 0; i < 100000000; ++i)
		{
			auto res = ApplyUnitaryOperation(test, cpyLambda);
		}

		auto endTime = std::chrono::high_resolution_clock::now();
		auto runtime = endTime - startTme;
		std::cout << "Unitary lambda copy run time mem 40 bill doubles copy  320GB=" << runtime.count() / 1000000000.0 << " seconds  = " << 320. / (runtime.count() / 1000000000.0) << " GB per sec" << std::endl;
	}
}



void testMemCpy2()
{
	for (int SZ = 100; SZ < 10000; SZ += 100)
	{
		using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
		std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(0.5));
		std::vector<FloatType>  c(SZ, VecXX::SCALA_TYPE(0.5));
		VecXX test(v);
		VecXX test2(v);
		test2 *= 0.0;


		std::cout << " \n Buffers number of doubles = , " << SZ;

		long loop = 100000;

		{
			auto startTme = std::chrono::high_resolution_clock::now();

			for (long i = 0; i < loop; ++i) //40 bil copy
			{
				std::memcpy(&c[0], &v[0], SZ * sizeof(FloatType));
			}

			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << " , MemCpy run , " << SZ * sizeof(FloatType) * loop / (1000000000.0) * 1.0 / (runtime.count() / 1000000000.0) << " ,  GB per sec";
		}

		{

			auto startTme = std::chrono::high_resolution_clock::now();

			auto cpyLambda = [](const auto& rhs) { return rhs; };

			for (long i = 0; i < loop; ++i) //40 bil copy
			{
				//auto res = ApplyUnitaryOperation(test, cpyLambda);
				auto res = transform( cpyLambda, test);
			}

			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << "			Unitary  run  , " << SZ * sizeof(FloatType) * loop / (1000000000.0) * 1.0 / (runtime.count() / 1000000000.0) << " ,  GB per sec";
		}
	}
}

class TimerGuard
{
	double& m_runTime;
	std::chrono::steady_clock::time_point  m_startTme;

public:
	TimerGuard(double& runTime) : m_runTime(runTime), m_startTme(std::chrono::high_resolution_clock::now()) { runTime = 0.; }

	~TimerGuard()
	{
		auto endTime = std::chrono::high_resolution_clock::now();
		auto runtime = endTime - m_startTme;
		m_runTime = runtime.count() / billion;
	}
};


//finds the max value in a vector
void doMax()
{
	const int TEST_LOOP_SZ = 1000;

	for (long int SZ = 200; SZ < 10000; SZ += 200)
	{
		auto v = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX test(v);

		// generic lambda for max either calling a max instruction or doing a selection with iff
		// auto mxDbl = [](auto lhs, auto rhs) { return max(lhs, rhs); };
		   auto mxDbl = [](auto lhs, auto rhs) { return iff(lhs > rhs, lhs, rhs); }; //using iff fastest 

		double time = 0.;
		auto runName = "";
		auto writeResults = [&](auto res) {std::cout << "size" << SZ << "," << runName << " result =, " << res << ", " << " Number of operations, " << numOps(TEST_LOOP_SZ, SZ) << ", run time  =, " << time << ", rate  =, " << numOps(TEST_LOOP_SZ, SZ) / time << ", , "; };

		{	runName = "STL  max elemen";
			volatile  double res = 0.;
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res = *std::max_element(cbegin(v), cend(v));
					}
				}
			}
			writeResults(res);
		}


		{
			runName = "DR3 reduce   max element";
			volatile  double res = 0.;
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res = reduce(test,mxDbl );
					}
				}
			}
			writeResults(res);
			std::cout << "\n";
		}
	}
}

//sums all elements in a vector
void doSum()
{
	const auto zero = InstructionTraits<VecXX::INS>::nullValue;
	const int TEST_LOOP_SZ = 1000;




	for (int SZ = 200; SZ < 60000; SZ += 200)
	{
		auto v = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX test(v);

		double time = 0.;
		auto runName = "";
		volatile  double res = 0.;

		auto writeResults = [&](auto res) {std::cout << "size" << SZ << "," << runName << " result =, " << res << ", " << " Number of operations, " << numOps(TEST_LOOP_SZ, SZ) << ", run time  =, " << time << ", rate  =, " << numOps(TEST_LOOP_SZ, SZ) / time << ", , "; };

		{	runName = "\n For loop accumulate";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res = 0.;
						for (auto x : v)
						{
							res += x;
						}
					}
				}
			}
			writeResults(res);
		}

		{	runName = "\n std::accumulate";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res = std::accumulate(v.begin(), v.end(), zero);
					}
				}
			}
			writeResults(res);
		}

		{	runName = "\n DR3 reduce";
			{  	auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res = reduce(test, Sum, zero, true);
					}
				}
			}
			writeResults(res);
			std::cout << "\n";
		}

	}
}


void doTransform()
{

	const int TEST_LOOP_SZ = 1000;

	for (long SZ = 200; SZ < 60000; SZ += 200)
	{
		auto v = getRandomShuffledVector(SZ); // std stl vector double or float 
		auto targetVec = v;
		auto transformVec = v;
		VecXX test(v);
		auto SQR = [](auto rhs) { return rhs * rhs; };

		double time = 0.;
		auto runName = "";

		auto writeResults = [&]() {std::cout << "size" << SZ << "," << runName << " Number of operations, " << numOps(TEST_LOOP_SZ, SZ) << ", run time  =, " << time << ", rate  =, " << numOps(TEST_LOOP_SZ, SZ) / time << ", , "; };

		{	runName = "sqr std::transform";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						std::transform(v.begin(), v.end(), targetVec.begin(), SQR);
					}
				}
			}
			writeResults();

		}

		{	runName = "DR3 transform";
			VecXX res;

			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res = transform(SQR, test);
					}
				}
			}
			writeResults();
			transformVec = res;
			
			bool testOK = vectorsEqual(transformVec, transformVec, targetVec);
			if (testOK)
			{
				std::cout << "Matching results";
			}
			else
			{
				std::cout << " FAIL results dont match";
			}
			std::cout << "\n";
		}
	}
}

void doInnerProd()
{

	const int TEST_LOOP_SZ = 1000;
	auto zero = InstructionTraits<VecXX::INS>::nullValue;

	for (long SZ = 200; SZ < 60000; SZ += 200)
	{

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		auto v2 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX t1(v1);
		VecXX t2(v2);



		double time = 0.;
		auto runName = "";
		volatile  double res = 0.;

		auto writeResults = [&](auto res) {std::cout << "size" << SZ << "," << runName << " result =, " << res << ", " << " Number of operations, " << numOps(TEST_LOOP_SZ, SZ) << ", run time  =, " << time << ", rate  =, " << numOps(TEST_LOOP_SZ, SZ) / time << ", , "; };


		{	runName = "std::inner_product";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res = inner_product(v1.cbegin(), v1.cend(), v2.cbegin(), zero);
					}
				}
			}
			writeResults(res);
		}


		{	runName = "DR3 transformReduce";
			{	auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
			    auto Mult = [](auto X, auto Y) { return X * Y; };
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res =transformReduce(t1, t2, Mult, Sum);
					}
				}
			}
			writeResults(res);
			std::cout << "\n";

		}

	}
}

void doSumSqrs()
{
	const int TEST_LOOP_SZ = 1000;
	const auto zero = InstructionTraits<VecXX::INS>::nullValue;

	for (long SZ = 200; SZ < 60000; SZ += 200)  //vector size
	{

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX t1(v1);



		double time = 0.;
		auto runName = "";
		volatile  double res = 0.;

		auto writeResults = [&](auto res) {std::cout << "size" << SZ << "," << runName << " result =, " << res << ", " << " Number of operations, " << numOps(TEST_LOOP_SZ, SZ) << ", run time  =, " << time << ", rate  =, " << numOps(TEST_LOOP_SZ, SZ) / time << ", , "; };

		{	runName = "sqrs std::inner_product";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res = std::inner_product(v1.cbegin(), v1.cend(), v1.cbegin(), zero);
					}
				}
			}
			writeResults(res);
		}


		{	runName = "DR3 sum sqrs transformReduce";
			{
				auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
				auto SQR = [](auto X) { return X * X; };
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res = transformReduce(t1, SQR, Sum);
					}
				}
			}
			writeResults(res);
			std::cout << "\n";
		}
	}
}

void khanAccumulation()
{
	long loop = 10000;

	auto NULL_Vec = VecXX::INS(0.0);
	const auto zero = InstructionTraits<VecXX::INS>::nullValue;
	//lambdas
	auto KhanAddV = [c = NULL_Vec](auto sum, auto rhs) mutable
	{
		auto y = rhs - c;
		auto t = sum + y;
		c = (t - sum);
		c = c - y;
		sum = t;
		return t;
	};

	for (long SZ = 200; SZ < 100000; SZ += 500)
	{
		VecXX test( VecXX::scalar(1.0 / 6.0), SZ);
		
		std::vector<InstructionTraits<VecXX::INS>::FloatType> v = test;

		double time = 0.;
		auto runName = "";
		double err = 0.;
		volatile  double res = 0.;
		auto writeResults = [&](auto res) {std::cout << "size" << SZ << "," << runName << " result =, " << res << ",  error =," << err << " Number of operations, " << numOps(loop, SZ) << ", run time  =, " << time << ", rate  =, " << numOps(loop, SZ) / time << ", , "; };

		{	runName = "std accumulate Sum";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < loop; l++)
					{
						res = std::accumulate(begin(v), end(v), zero);
					}
				}
				err = (res - SZ / 6.0) / res;
			}
			writeResults(res);
		}


		{	runName = "DR3 reduce Khan Sum";
			{	
				TimerGuard timer(time);
				{
					for (long l = 0; l < loop; l++)
					{
						res = reduce(test, KhanAddV, zero);
					}
				}
				err = (res - SZ / 6.0) / res;
			}
			writeResults(res);
			std::cout << "\n";
		}
	}
}


void testBinarySelection()
{
	
	const int TEST_LOOP_SZ = 1000;

	for (long SZ = 200; SZ < 60000; SZ += 200)
	{
		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX testVec(v1);

		auto C = v1; //copy of STL vector to write vakues to
		double time = 0.;
		auto runName = "";

		//use auto on scalars so we can switch between float and double instruction sets too
		auto one = VecXX::scalar(1.0);
		auto two = VecXX::scalar(2.0);
		auto half = VecXX::scalar(0.5);

		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };


		auto truVal = one;
		auto falseVal = two;

		auto writeResults = [&]() {std::cout << "size" << SZ << "," << runName <<  ", " << " Number of operations, " << numOps(TEST_LOOP_SZ, SZ) << ", run time  =, " << time << ", rate  =, " << numOps(TEST_LOOP_SZ, SZ) / time << ", , "; };

		{	runName = "for loop";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						for (int k = 0; k < SZ; k++)
						{
							auto x = v1[k];
							if ((x - two * floor(x * half)) >= one)
							{
								C[k] = one;
							}
							else
							{
								C[k] = two;
							}
						}
					}
				}
			}
			writeResults();
		}

		// selects true or false value into result res on the basisi of the lambda MyOddLmbda
		{	runName = "DR3 select constants boolean lambda";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						volatile auto res = select(MyOddLmbda, testVec, truVal, falseVal); //was ApplySelectionOperationC
					}
				}
			}
			writeResults();
			std::cout << "\n";
		}

	}
}


//selecting between very light weight functions
void testBinarySelection1()
{
	const int TEST_LOOP_SZ = 1000;

	for (long SZ = 200; SZ < 60000; SZ += 200)
	{

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX testVec(v1);

		auto C = v1; //copy of STL vector to write vakues to
		auto C1 = v1;
		auto C2 = v1;

		double time = 0.;
		auto runName = "";
		volatile  double res = 0.;
		res *= res;

		//use auto on scalars so we can switch between float and double instruction sets too
		const auto one = VecXX::scalar(1.0);
		const auto two = VecXX::scalar(2.0);
		const auto half = VecXX::scalar(0.5);

		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

		// making sure scalar constants are coerced to type of  float associated with chosen instruction set for ease opf switching sets
		auto trueLambda = [](auto x) { return  (x * VecXX::scalar(3.7) + VecXX::scalar(16.2));  };
		auto falseLambda = [](auto x) { return  x * ((x * VecXX::scalar(3.5) ) + VecXX::scalar(12.2) );  };


		auto writeResults = [&]() {std::cout << "size" << SZ << "," << runName <<  ", " << " Number of operations, " << numOps(TEST_LOOP_SZ, SZ) << ", run time  =, " << time << ", rate  =, " << numOps(TEST_LOOP_SZ, SZ) / time << ", , "; };

		{
			runName = "for loop";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{

						for (int k = 0; k < SZ; k++)
						{
							auto x = v1[k];
							if ((x - two * floor(x * half)) >= one )
							{
								C[k] = (x * VecXX::scalar(3.7) + VecXX::scalar(16.2));
							}
							else
							{
								C[k] = x * ((x * VecXX::scalar(3.5)) + VecXX::scalar(12.2));
							}
						}
					}
				}
			}
			writeResults();
		}

		{	runName = "std::transform conditional";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						std::transform(v1.begin(), v1.end(), C1.begin(),
										[&](auto x)
										{
											if (MyOddLmbda(x))
											{
												return (x * VecXX::scalar(3.7) + VecXX::scalar(16.2));
											}
											return x * ((x * VecXX::scalar(3.5)) + VecXX::scalar(12.2));
										}
									);
					}
				}
			}
			writeResults();

		}


		{	runName = "DR3 selectTransform ";
		    auto resDr = testVec;
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						resDr = selectTransform(MyOddLmbda, testVec, trueLambda, falseLambda);
					}
				}
			}
			writeResults();
			
			C2 = resDr; //fill up stl vector for compare
		}

		//////////compare ///////////
		bool testOK = vectorsEqual(C, C1, C2);
		if (testOK)
		{
			std::cout << "Matching results";
		}
		else
		{
			std::cout << " FAIL results dont match";
		}
		std::cout << "\n";

	}


}


//selecting between middle weight functions
void testBinarySelection2()
{
	using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
	const int TEST_LOOP_SZ = 1000;

	for (long SZ = 200; SZ < 60000; SZ += 200)
	{

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float
		VecXX testVec(v1);

		auto C = v1; //copy of STL vector to write values to
		auto C1 = v1;
		auto C2 = v1;

		double time = 0.;
		auto runName = "";
		const auto one = VecXX::scalar(.0);
		const auto two = VecXX::scalar(2.0);
		const auto half = VecXX::scalar(0.5);

		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

		/// from acklams inverse cdf normal
		static FloatType a[] = { 0.0,  -3.969683028665376e+01, 2.209460984245205e+02,-2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00 };
		static FloatType b[] = { 0.0, -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,  6.680131188771972e+01, -1.328068155288572e+01 };
		static FloatType c[] = { 0.0,-7.784894002430293e-03,-3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
		static FloatType d[] = { 0.0,  7.784695709041462e-03, 3.224671290700398e-01,  2.445134137142996e+00, 3.754408661907416e+00 };


		auto trueLambda = [&](auto q)
		{
			auto X = (((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
				((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + VecXX::scalar(1.0));
			return X;
		};


		auto falseLambda = [&](auto q)
		{
			auto X = (((((a[1] * q + a[2]) * q + a[3]) * q + a[4]) * q + a[5]) * q + a[6]) /
				((((b[1] * q + b[2]) * q + b[3]) * q + b[4]) * q + VecXX::scalar(1.0));
			return X;
		};

		auto writeResults = [&]() { std::cout << "size" << SZ << "," << runName <<  ", " << " # calcs, " << numOps(TEST_LOOP_SZ , SZ) << ", run time  =, " << time << ", rate  =, " << numOps(TEST_LOOP_SZ, SZ) / time << ", , "; };

		
		{	runName = "for loop";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{

						for (int k = 0; k < SZ; k++)
						{
							auto x = v1[k];
							C[k] = MyOddLmbda(x) ? trueLambda(x) : falseLambda(x);
						}
					}
				}
			}
			writeResults();
		}

		//calculate both true and false lambdas for all values and select/blend values together
		{	runName = "DR3 selectTransform ";
			VecXX res;
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res = selectTransform(MyOddLmbda, testVec, trueLambda, falseLambda);
					}
				}
			}
			writeResults();
			C1 = res;
		}
		

	
		//binary filter into concrete views  of true and false elements according to MyOddLambda. Then apply  true and false lambdas and finally  merge values together
		{	runName = "DR3 filterTransform";
			VecXX res;
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						res = filterTransform(MyOddLmbda, testVec, trueLambda, falseLambda);
					}
				}
			}
			writeResults();
			C2 = res;
		}

		//compare resulyts of calcs
		bool testOK = vectorsEqual(C1, C2, C);

		if (testOK)
		{
			std::cout << "Matching results";
		}
		else
		{
			std::cout << " FAIL results dont match";
		}
		std::cout << "\n";
	}

}



//heavy weight one odd lambda  functions
void testBinarySelection3()
{
	using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
	const int TEST_LOOP_SZ = 1000;

	for (long SZ = 200; SZ < 60000; SZ += 200)
	{

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float
		VecXX testVec(v1);

		auto C = v1; //copy of STL vector to write vakues to
		auto C1 = v1;
		auto C2 = v1;


		double time = 0.;
		auto runName = "";
		volatile  double res = 0.;
		res *= res;

		/// acklams inverse cdf normal style coefficients
		static FloatType c[] = { 0.0,-7.784894002430293e-03,-3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
		static FloatType d[] = { 0.0,  7.784695709041462e-03, 3.224671290700398e-01,  2.445134137142996e+00, 3.754408661907416e+00 };

		const auto one = VecXX::scalar(.0);
		const auto two = VecXX::scalar(2.0);
		const auto half = VecXX::scalar(0.5);

		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

		auto trueLambda = [&](auto q)
		{
			auto X = (((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
				((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + VecXX::scalar(1.0));
			return X;
		};


		//auto falseLambda = [](auto x) { return exp(-sin(x / VecXX::INS(20.))) + VecXX::INS(1.0) / (VecXX::INS(1.0) + exp(x * (sin(x * VecXX::INS(3.7) + VecXX::INS(12.2) / x))));  };
		auto falseLambdaS = [](auto x) { return exp(-sin(x / VecXX::scalar( 20.))) + VecXX::scalar(1.0) / (VecXX::scalar(1.0) + exp(x * (sin(x * VecXX::scalar(3.7) + VecXX::scalar(12.2) / x))));  };
		auto writeResults = [&]() { std::cout << "size" << SZ << "," << runName << ", " << " # ops, " << numOps(TEST_LOOP_SZ, SZ) << ", run time  =, " << time << ", rate  =, " << numOps(TEST_LOOP_SZ, SZ) / time << ", , "; };
		

		{	runName = "for loop";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						for (int k = 0; k < SZ; k++)
						{
							auto x = v1[k];
							C[k] = ( MyOddLmbda(x)) ? trueLambda(x) : falseLambdaS(x);

						}
					}
				}
			}
			writeResults();
		}



		{	//performs both true and false lambda and then selects/blends the results together.
			runName = "DR3 selectTransform";
		    VecXX resDr;
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						resDr = selectTransform(MyOddLmbda, testVec, trueLambda, falseLambdaS);
					}
				}
			}
			writeResults();
			C1 = resDr;// copy to STL vector
		}


		{	runName = "DR3 filterTransform";
			VecXX resDr;
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						 resDr = filterTransform(MyOddLmbda, testVec, trueLambda, falseLambdaS);
					}
				}
			}
			writeResults();
			C2 = resDr;
		}


		//////////compare calc results///////////
		bool testOK = vectorsEqual(C, C1, C2);


		if (testOK)
		{
			std::cout << "Matching results";
		}
		else
		{
			std::cout << " FAIL results dont match";
		}
		std::cout << "\n";

	}

}


void doCountIf()
{

	const int TEST_LOOP_SZ = 1000;
	//auto zero = InstructionTraits<VecXX::INS>::nullValue;

	for (long SZ = 200; SZ < 60000; SZ += 200)
	{

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX test(v1);


		double time = 0.;
		auto runName = "";
		volatile  double res = 0.;

		auto writeResults = [&](auto res) {std::cout << "size" << SZ << "," << runName << " result =, " << res << ", " << " # ops, " << numOps(TEST_LOOP_SZ, SZ) << ", run time  =, " << time << ", rate  =, " << numOps(TEST_LOOP_SZ, SZ) / time << ", , "; };

		{	runName = "STL  count_if";
			long long  resStl = 0;
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{

						auto isOdd = [](auto x) {
							auto y = (x - (2. * floor(x * 0.5)));
							if (y >= 1.0)
								return false;
							else return true;
						};
						resStl = std::count_if(begin(v1), end(v1), isOdd);
					}
				}	
			}
			writeResults(resStl);
		}
		
		{	runName = "DR3 transformReduce";
			{
				TimerGuard timer(time);
				{
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						auto oneIfOddLmbda = [](auto x) { return  iff((x - (VecXX::INS(2.0) * floor(x * VecXX::INS(0.5)))) >= VecXX::INS(1.0), VecXX::INS(1.0), VecXX::INS(0.0)); };
						auto Add = [](auto x, auto y) { return  x + y;  };
						res = transformReduce(test, oneIfOddLmbda, Add);
					}
				}
			}
			writeResults(res);
			std::cout << "\n";
		}
	}
}










