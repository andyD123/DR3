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
#include <thread>
#include <map>



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

//using namespace DRC::VecDb;
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

bool valuesAreEqual(double x, double y,double tol =  1e-14)
{
	auto err1 = fabs((x - y) / (x + y));

	return (err1 > tol) ? false : true;
}


auto getRandomShuffledVector(int SZ, int instance_number=0)
{
	using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;


	static std::map<int, std::vector<FloatType> > vectors;


	int key = 10 * SZ + instance_number;
	//store vectors with key 10 times size  and add on 0-9 integer for instance of different random vector

	if (SZ < 0)
	{
		vectors.clear();
		SZ = 0;
	}


	if (vectors.find(key) != vectors.end())
	{
		return vectors[key];
	}
	else
	{
		std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
		for (int i = 0; i < SZ; i++) { v[i] += FloatType(SZ / 2) + i; }
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(begin(v), end(v), g);
		vectors[key] = v;
		return v;
	}
}




auto numOps = [](int TEST_LOOP_SZ, int SZ) { return  static_cast<int>(double(TEST_LOOP_SZ) * double(SZ)); };


double getnull(double)
{
	return 0.0;
}

using Calc_Values = std::map<int, double>;
using  Mapped_Performance_Results = std::map<int, std::vector<double> >; // array size  v vector<throughput for runs>
using Mapped_Stats = std::map<int, std::pair<double, double> >; // size -.pair ( throughput ,  std dev of through put)

struct RunResults
{
	Mapped_Performance_Results m_raw_results;
	Calc_Values  m_calc_results;
};


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


auto runFunctionOverDifferentSize = [](int testRepeats, int vec_start_size, int vec_stepSZ, int vec_maxSZ, const auto& func, long testLoopSZ)
{

	RunResults results;

	for (int j = 0; j < testRepeats; ++j)
	{
		int VEC_SZ = vec_start_size;
		for (; VEC_SZ < vec_maxSZ; VEC_SZ += vec_stepSZ)
		{
			auto res = func(VEC_SZ, testLoopSZ);
			auto calculation_rate = res.second;
			auto calc_value = res.first;
			results.m_raw_results[VEC_SZ].push_back(calculation_rate);
			if (j == 0)
			{
				results.m_calc_results[VEC_SZ] = calc_value;
			}
		}
	}
	return results;
};


auto runFunctionOverDifferentSizeVec = [](int testRepeats, int vec_start_size, int vec_stepSZ, int vec_maxSZ, const auto& func, long testLoopSZ)
{

	RunResults results;

	for (int j = 0; j < testRepeats; ++j)
	{
		int VEC_SZ = vec_start_size;
		for (; VEC_SZ < vec_maxSZ; VEC_SZ += vec_stepSZ)
		{
			auto res = func(VEC_SZ, testLoopSZ);
			auto calculation_rate = res.second;
			auto calc_value = res.first;
			results.m_raw_results[VEC_SZ].push_back(calculation_rate);
			//if (j == 0)
			//{
			//	results.m_calc_results[VEC_SZ] = calc_value;
			//}
		}
	}
	return results;
};



auto performanceStats = [](const Mapped_Performance_Results& raw_results)
{

	Mapped_Stats stats;

	for (const auto& item : raw_results)
	{
		double sum = 0;
		double sum_sqrd = 0;
		double N = 0.0;
		for (const auto run_rate : item.second)
		{
			sum += run_rate;
			sum_sqrd += (run_rate * run_rate);
			N++;
		}

		double avg = sum / N;
		double varSqrd = sum_sqrd + (avg * avg * N) - (2.0 * avg * sum);
		double var = std::sqrt(varSqrd / (N - 1.));

		stats[item.first] = { avg ,var };

	}
	return stats;
};




//example functions fwd decl

void	testMemCpy2();
void    doMax();
void 	doSum(); 
void    doInnerProd();
void	doTransform();
void 	doSumSqrs();
void    khanAccumulation();

void	binarySelectionBetweenConst();
void	binarySelectionBetweenLinearFunction(); // y= mx + c    a couple of of operations
void	testBinarySelection2();
void    binarySelectionBetweenMiddleWeightFunction();
void    testBinarySelection3();
void	binarySelectionBetweenHeavyWeightFunction();
void	doCountIf();


void doAVXMax512Dance();







int main()
{

 //Uncomment  a function to play with
 // 
 //  transform 
	//	doTransform();

	std::cout << "\n \n \n \n testMemCpy2() \n" << std::endl;
		testMemCpy2(); 
//accumulate 
// 
		std::cout << "\n \n \n \n doSum() \n" << std::endl;
    	doSum(); // stl slower with intel  stl slower

		std::cout << "\n \n \n \n doMax() \n"  << std::endl;
	    doMax();

//	    doAVXMax512Dance();

//transform accum
		std::cout << "\n \n \n \n doInnerProd() \n" << std::endl;
	  doInnerProd();


	  std::cout << "\n \n \n \n doSumSqrs() \n" << std::endl;
		doSumSqrs();
// lambda capture

	 std::cout << "\n \n \n \n khanAccumulation() \n" << std::endl;
	 khanAccumulation();

// branching
	 std::cout << "\n \n \n \n binarySelectionBetweenConst() \n" << std::endl;

	binarySelectionBetweenConst(); //select between constants


	std::cout << "\n \n \n \n binarySelectionBetweenLinearFunction() \n" << std::endl;
	//	testBinarySelection1();//select between light functions
	binarySelectionBetweenLinearFunction();

	std::cout << "\n \n \n \n binarySelectionBetweenMiddleWeightFunction() \n" << std::endl;
	//	testBinarySelection2(); //medium
	binarySelectionBetweenMiddleWeightFunction();

	std::cout << "\n \n \n \n binarySelectionBetweenHeavyWeightFunction() \n" << std::endl;
	//  testBinarySelection3(); //heavy	
	binarySelectionBetweenHeavyWeightFunction();
	
	std::cout << "\n \n \n \n doCountIf() \n" << std::endl;
	doCountIf();



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



//finds the max value in a vector
/*
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

*/

/*
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


*/



/*
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
//
*/


void doSum()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	auto zero = InstructionTraits<VecXX::INS>::nullValue;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto accumulate_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ, 0); 

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = std::accumulate(v1.begin(), v1.end(), zero);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = std::accumulate(v1.begin(), v1.end(), zero);
				}
			}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};


	auto accumulate_for_each = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ, 0);

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = 0.;
				for (auto x : v1)
				{
					res += x;
				}

			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = 0.;
					for (auto x : v1)
					{
						res += x;
					}
				}
			}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_accumulate = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(SZ, 0); // std stl vector double or float 
		VecXX t1(v1);
		{
			auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = reduce(t1, Sum, zero, true);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = reduce(t1, Sum, zero, true);
				}
			}
		}
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};

	auto run_res_for_each = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, accumulate_for_each, TEST_LOOP_SZ);
	auto stats_foreach = performanceStats(run_res_for_each.m_raw_results);


	auto run_res_stl = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, accumulate_run, TEST_LOOP_SZ);
	auto stats_stl = performanceStats(run_res_stl.m_raw_results);


	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_accumulate, TEST_LOOP_SZ);
	auto stats_DR3_inner_prod = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_stl)
	{
		auto  forLoop = run_res_for_each.m_calc_results[elem.first]; 
		auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = run_res_stl.m_calc_results[elem.first];
		auto strMatch = valuesAreEqual(valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "forloop , size " << elem.first << " , " << stats_foreach[elem.first].first << " +- " << stats_foreach[elem.first].second << "\t \t  std::accumulate , size " << elem.first << " , " << elem.second.first << " +- " << elem.second.second << "\t \t DR3 accumulate , size " << elem.first << " , " << stats_DR3_inner_prod[elem.first].first << " +- " << stats_DR3_inner_prod[elem.first].second << ", numerical check: " << strMatch << "\n";
	}

}



void doMax()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	/*
	const int maxVectorSize = 4400;
	const int minVectorSize = 3800;
	const long TEST_LOOP_SZ = 100000;
	const int vectorStepSize = 8;
	const int repeatRuns = 10;
	*/
	auto zero = InstructionTraits<VecXX::INS>::nullValue;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto accumulate_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		//auto v = getRandomShuffledVector(SZ); // std stl vector double or float 
		auto v1 = getRandomShuffledVector(VEC_SZ, 0);

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = *std::max_element(v1.begin(), v1.end());
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = *std::max_element(v1.begin(), v1.end());
				}
			}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};



	auto DR3_accumulate = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;


		// generic lambda for max either calling a max instruction or doing a selection with iff
	//	 auto mxDbl = [](auto lhs, auto rhs) { return max(lhs, rhs); };
		auto mxDbl = [](auto lhs, auto rhs) { return iff(lhs > rhs, lhs, rhs); }; //using iff fastest 

		auto v1 = getRandomShuffledVector(SZ, 0); // std stl vector double or float 
		VecXX vec(v1);
		{
			
			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = reduce(vec, mxDbl);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = reduce(vec, mxDbl);
				}
			}
		}
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};



	auto run_res_stl = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, accumulate_run, TEST_LOOP_SZ);
	auto stats_stl = performanceStats(run_res_stl.m_raw_results);


	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_accumulate, TEST_LOOP_SZ);
	auto stats_DR3_perf = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& perf_stl : stats_stl)
	{
		
		auto  valDr3 = dr3_raw_results.m_calc_results[perf_stl.first];
		auto  valStl = run_res_stl.m_calc_results[perf_stl.first];
		auto strMatch = valuesAreEqual(valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout <<"  std::max_element, size " << perf_stl.first << ", " << perf_stl.second.first << " + -" << perf_stl.second.second << "\t \t DR3 reduce, size " << perf_stl.first << ", " << stats_DR3_perf[perf_stl.first].first << " + -" << stats_DR3_perf[perf_stl.first].second << ", numerical check : " << strMatch << "\n";
	}


}









/*
void doTransform()
{

	const int TEST_LOOP_SZ = 1000;

	for (long SZ = 200; SZ < 60000; SZ += 200)
	{
		auto v = getRandomShuffledVector(SZ); // std stl vector double or float 
		auto targetVec = v;
		auto transformVec = v;
		const VecXX test(v);
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
						//res = transformXX(SQR, test);
						//res = transform1(SQR, test);
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
*/

//

void doTransform() //best with avx2 
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 400;
	const int maxVectorSize = 60000;
	const int minVectorSize = 400;

	auto zero = InstructionTraits<VecXX::INS>::nullValue;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto stl_transform = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		auto v = getRandomShuffledVector(VEC_SZ, 0); // std stl vector double or float 
		auto targetVec = v;
		volatile auto transformVec = v;
		VecXX test(v);
		auto SQR = [](auto rhs) { return rhs * rhs; };

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				std::transform(v.begin(), v.end(), targetVec.begin(), SQR);
				//res = inner_product(v1.cbegin(), v1.cend(), v2.cbegin(), zero);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					std::transform(v.begin(), v.end(), targetVec.begin(), SQR);

				}
			}
		}
		return  std::make_pair(targetVec, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_transform = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		auto v = getRandomShuffledVector(VEC_SZ, 0); // std stl vector double or float 
		auto targetVec = v;
		auto res = v;
		const VecXX test(v);
		auto SQR = [](auto rhs) { return rhs * rhs; };


		{
			auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
			auto Mult = [](auto X, auto Y) { return X * Y; };

			//warm up
			for (long l = 0; l < 100; l++)
			{
				//res = transformXX(SQR, test);
				res = transform(SQR, test);

			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = transform(SQR, test);
					//res = transformXX(SQR, test);
				}
			}
		}
		return std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);

	};

	auto dr3_raw_results = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_transform, TEST_LOOP_SZ);
	auto stats_DR3 = performanceStats(dr3_raw_results.m_raw_results);

	auto run_res_innerProd = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, stl_transform, TEST_LOOP_SZ);
	auto stats_inner_prod = performanceStats(run_res_innerProd.m_raw_results);





	//print out results

	for (const auto& elem : stats_inner_prod)
	{
		auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = run_res_innerProd.m_calc_results[elem.first];
		auto strMatch = "TO DO";// vectorsEqual(valDr3, valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "STL transform , size " << elem.first << " , " << stats_inner_prod[elem.first].first << " +- " << stats_inner_prod[elem.first].second << "\t \t DR3 transform , size " << elem.first << " , " << stats_DR3[elem.first].first << " +- " << stats_DR3[elem.first].second << ", numerical check: " << strMatch << "\n";
	}
	
	



}
//





void doInnerProd()
{

	const long TEST_LOOP_SZ =  1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	auto zero = InstructionTraits<VecXX::INS>::nullValue;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto inner_prod_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ,0); // std stl vector double or float 
		auto v2 = getRandomShuffledVector(VEC_SZ,1); // std stl vector double or float 

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = inner_product(v1.cbegin(), v1.cend(), v2.cbegin(), zero);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = inner_product(v1.cbegin(), v1.cend(), v2.cbegin(), zero);
				}
			}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_inner_prod = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(SZ,0); // std stl vector double or float 
		auto v2 = getRandomShuffledVector(SZ,1); // std stl vector double or float 
		VecXX t1(v1);
		VecXX t2(v2);

		{	
			auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
			auto Mult = [](auto X, auto Y) { return X * Y; };

			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = transformReduce(t1, t2, Mult, Sum);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = transformReduce(t1, t2, Mult, Sum);
				}
			}
		}
		return std::make_pair(res , numOps(TEST_LOOP_SZ, SZ) / time );

	};


	
	auto run_res_innerProd = runFunctionOverDifferentSize(repeatRuns, minVectorSize,vectorStepSize, maxVectorSize, inner_prod_run, TEST_LOOP_SZ);
	auto stats_inner_prod = performanceStats(run_res_innerProd.m_raw_results);


	auto dr3_raw_results  =  runFunctionOverDifferentSize(repeatRuns, minVectorSize,vectorStepSize, maxVectorSize, DR3_inner_prod, TEST_LOOP_SZ);
	auto stats_DR3_inner_prod = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_inner_prod)
	{
		auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = run_res_innerProd.m_calc_results[elem.first];
		auto strMatch = valuesAreEqual(valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "STL inner product , size " <<  elem.first << " , " << elem.second.first << " +- " << elem.second.second << "\t \t DR3 inner product , size " << elem.first << " , " << stats_DR3_inner_prod[elem.first].first << " +- " << stats_DR3_inner_prod[elem.first].second << ", numerical check: " << strMatch << "\n";
	}

	

}

/*
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

*/


void doSumSqrs()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	auto zero = InstructionTraits<VecXX::INS>::nullValue;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto inner_prod_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		//auto v2 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = inner_product(v1.cbegin(), v1.cend(), v1.cbegin(), zero);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = inner_product(v1.cbegin(), v1.cend(), v1.cbegin(), zero);
				}
			}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_inner_prod = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX t1(v1);


		{
			auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
			auto Mult = [](auto X, auto Y) { return X * Y; };

			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = transformReduce(t1, t1, Mult, Sum);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = transformReduce(t1, t1, Mult, Sum);
				}
			}
		}
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};



	auto run_res_innerProd = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, inner_prod_run, TEST_LOOP_SZ);
	auto stats_inner_prod = performanceStats(run_res_innerProd.m_raw_results);


	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_inner_prod, TEST_LOOP_SZ);
	auto stats_DR3_inner_prod = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_inner_prod)
	{
		auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = run_res_innerProd.m_calc_results[elem.first];
		auto strMatch = valuesAreEqual(valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "STL inner product sum sqrs , size " << elem.first << " , " << elem.second.first << " +- " << elem.second.second << "\t \t DR3 inner product sum sqrs , size " << elem.first << " , " << stats_DR3_inner_prod[elem.first].first << " +- " << stats_DR3_inner_prod[elem.first].second << ", numerical check: " << strMatch << "\n";
	}



}


/*
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

*/




void khanAccumulation()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	auto zero = InstructionTraits<VecXX::INS>::nullValue;

	auto stl_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		double relative_error = 0.;
		VecXX test(VecXX::scalar(1.0 / 6.0), VEC_SZ);
		std::vector<InstructionTraits<VecXX::INS>::FloatType> v1 = test;

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				
				res = std::accumulate(begin(v1), end(v1), zero);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = std::accumulate(begin(v1), end(v1), zero);
				}
			}

			relative_error = (VEC_SZ / 6. - res) / res;
		}
		return  std::make_pair(relative_error, numOps(TEST_LOOP_SZ, VEC_SZ) / time);

	};

	auto DR3_reduce = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double relative_error = 0.;
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


		double time = 0.;
		volatile  double res = 0.;
		VecXX t1(VecXX::scalar(1.0 / 6.0), VEC_SZ);

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = reduce(t1, KhanAddV, zero);

			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = reduce(t1, KhanAddV, zero);
				}

			}
		}

		relative_error = (VEC_SZ / 6. - res) / res;
		return std::make_pair(relative_error, numOps(TEST_LOOP_SZ, VEC_SZ) / time);

	};



	auto run_res_innerProd = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, stl_run, TEST_LOOP_SZ);
	auto stats_inner_prod = performanceStats(run_res_innerProd.m_raw_results);


	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_reduce, TEST_LOOP_SZ);
	auto stats_DR3_inner_prod = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_inner_prod)
	{
		auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = run_res_innerProd.m_calc_results[elem.first];
		//auto strMatch = valuesAreEqual(valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "STL accumulate , size " << elem.first << " , " << elem.second.first << " +- " << elem.second.second <<"relative err"<< valStl << "\t \t DR3 khan accumulate , size " << elem.first << " , " << stats_DR3_inner_prod[elem.first].first << " +- " << stats_DR3_inner_prod[elem.first].second << ", relative err " << valDr3 << "\n";
	}



}




/*

void binarySelectionBetweenConst()
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
*/


void binarySelectionBetweenConst()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	auto zero = InstructionTraits<VecXX::INS>::nullValue;
	auto one = VecXX::scalar(1.0);
	auto two = VecXX::scalar(2.0);
	auto half = VecXX::scalar(0.5);
	auto truVal = one;
	auto falseVal = two;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto for_loop_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		auto C = v1;
		//use auto on scalars so we can switch between float and double instruction sets too


		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				for (int k = 0; k < VEC_SZ; k++)
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

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					for (int k = 0; k < VEC_SZ; k++)
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
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_select = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		auto one = VecXX::scalar(1.0);
		auto two = VecXX::scalar(2.0);

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX testVec(v1);


		{
			auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

			//warm up
			for (long l = 0; l < 100; l++)
			{
				auto res = select(MyOddLmbda, testVec, truVal, falseVal); //was ApplySelectionOperationC
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					auto res = select(MyOddLmbda, testVec, truVal, falseVal);
				}
			}
		}
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};



	auto run_res_innerProd = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, for_loop_run, TEST_LOOP_SZ);
	auto stats_inner_prod = performanceStats(run_res_innerProd.m_raw_results);


	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_select, TEST_LOOP_SZ);
	auto stats_DR3_inner_prod = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_inner_prod)
	{
		auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = run_res_innerProd.m_calc_results[elem.first];
		auto strMatch = "TO DO";//"valuesAreEqual(valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "for loop binarySelectionBetweenConst , size " << elem.first << " , " << elem.second.first << " +- " << elem.second.second << "\t \t DR3 binarySelectionBetweenConst , size " << elem.first << " , " << stats_DR3_inner_prod[elem.first].first << " +- " << stats_DR3_inner_prod[elem.first].second << ", numerical check: " << strMatch << "\n";
	}



}
////////




/*
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
*/

//binarySelectionBetweenLinearFunction

void binarySelectionBetweenLinearFunction()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	auto zero = InstructionTraits<VecXX::INS>::nullValue;
	auto one = VecXX::scalar(1.0);
	auto two = VecXX::scalar(2.0);
	auto half = VecXX::scalar(0.5);
	auto truVal = one;
	auto falseVal = two;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto for_loop_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		auto C = v1;
		//use auto on scalars so we can switch between float and double instruction sets too


		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				for (int k = 0; k < VEC_SZ; k++)
				{
					auto x = v1[k];
					if ((x - two * floor(x * half)) >= one)
					{
						C[k] = (x * VecXX::scalar(3.7) + VecXX::scalar(16.2));
					}
					else
					{
						C[k] = x * ((x * VecXX::scalar(3.5)) + VecXX::scalar(12.2));
					}
				}
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					for (int k = 0; k < VEC_SZ; k++)
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
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};



	auto stl_transform_with_branch_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		auto C = v1;
		//use auto on scalars so we can switch between float and double instruction sets too


		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				std::transform(v1.begin(), v1.end(), C.begin(),
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

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					std::transform(v1.begin(), v1.end(), C.begin(),
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
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};




	auto DR3_select = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		auto one = VecXX::scalar(1.0);
		auto two = VecXX::scalar(2.0);
		// making sure scalar constants are coerced to type of  float associated with chosen instruction set for ease opf switching sets
		auto trueLambda = [](auto x) { return  (x * VecXX::scalar(3.7) + VecXX::scalar(16.2));  };
		auto falseLambda = [](auto x) { return  x * ((x * VecXX::scalar(3.5)) + VecXX::scalar(12.2));  };

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX testVec(v1);


		{
			auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

			//warm up
			for (long l = 0; l < 100; l++)
			{
				auto res = selectTransform(MyOddLmbda, testVec, trueLambda, falseLambda);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					auto res = selectTransform(MyOddLmbda, testVec, trueLambda, falseLambda);
				}
			}
		}
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};



	auto run_res_innerProd = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, for_loop_run, TEST_LOOP_SZ);
	auto stats_inner_prod = performanceStats(run_res_innerProd.m_raw_results);

	auto run_transformBranchy = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, stl_transform_with_branch_run, TEST_LOOP_SZ);
	auto stats_transform_branchy = performanceStats(run_transformBranchy.m_raw_results);


	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_select, TEST_LOOP_SZ);
	auto stats_DR3 = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_inner_prod)
	{
		//auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		//auto  valStl = run_res_innerProd.m_calc_results[elem.first];
		auto strMatch = "TO DO";//"valuesAreEqual(valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "for loop binarySelectionBetweenSimpleFunctions , size " << elem.first << " , " << elem.second.first << " +- " << elem.second.second << "\t \t stl transform branchy lambda  , size " << elem.first << " , " << stats_transform_branchy[elem.first].first << " +- " << stats_transform_branchy[elem.first].second << "\t \t DR3 binarySelectionBetweenVSimpleLambda , size " << elem.first << " , " << stats_DR3[elem.first].first << " +- " << stats_DR3[elem.first].second << ", numerical check: " << strMatch << "\n";
	}



}

//selecting between middle weight functions
void testBinarySelection2()
{
	using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
	const int TEST_LOOP_SZ = 1000;

	for (long SZ = 200; SZ < 20000; SZ += 200)
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

void binarySelectionBetweenMiddleWeightFunction()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	auto zero = InstructionTraits<VecXX::INS>::nullValue;
	auto one = VecXX::scalar(1.0);
	auto two = VecXX::scalar(2.0);
	auto half = VecXX::scalar(0.5);


	auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };
	using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
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



	getRandomShuffledVector(-1); // reset  random input vectors

	auto for_loop_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		auto C = v1;
		//use auto on scalars so we can switch between float and double instruction sets too


		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				for (int k = 0; k < VEC_SZ; k++)
				{
					auto x = v1[k];
					if ((x - two * floor(x * half)) >= one)
					{
						C[k] = trueLambda(x);
					}
					else
					{
						C[k] = falseLambda(x);
					}
				}
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					for (int k = 0; k < VEC_SZ; k++)
					{
						auto x = v1[k];
						if ((x - two * floor(x * half)) >= one)
						{
							C[k] = trueLambda(x);
						}
						else
						{
							C[k] = falseLambda(x);
						}
					}
				}
			}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};



	auto stl_transform_with_branch_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		auto C = v1;
		//use auto on scalars so we can switch between float and double instruction sets too


		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				std::transform(v1.begin(), v1.end(), C.begin(),
					[&](auto x)
					{
						return MyOddLmbda(x) ? trueLambda(x) : falseLambda(x);
					}
				);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					std::transform(v1.begin(), v1.end(), C.begin(),
						[&](auto x)
						{
							return MyOddLmbda(x) ? trueLambda(x) : falseLambda(x);
						}
					);
				}
			}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};




	auto DR3_select = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		auto one = VecXX::scalar(1.0);
		auto two = VecXX::scalar(2.0);
		// making sure scalar constants are coerced to type of  float associated with chosen instruction set for ease opf switching sets
		//auto trueLambda = [](auto x) { return  (x * VecXX::scalar(3.7) + VecXX::scalar(16.2));  };
		//auto falseLambda = [](auto x) { return  x * ((x * VecXX::scalar(3.5)) + VecXX::scalar(12.2));  };

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX testVec(v1);


		{
			auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

			//warm up
			for (long l = 0; l < 100; l++)
			{
				auto res = selectTransform(MyOddLmbda, testVec, trueLambda, falseLambda);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					auto res = selectTransform(MyOddLmbda, testVec, trueLambda, falseLambda);
				}
			}
		}
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};



	auto DR3_filterTransdform = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		auto one = VecXX::scalar(1.0);
		auto two = VecXX::scalar(2.0);
		// making sure scalar constants are coerced to type of  float associated with chosen instruction set for ease opf switching sets
		//auto trueLambda = [](auto x) { return  (x * VecXX::scalar(3.7) + VecXX::scalar(16.2));  };
		//auto falseLambda = [](auto x) { return  x * ((x * VecXX::scalar(3.5)) + VecXX::scalar(12.2));  };

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX testVec(v1);


		{
			auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

			//warm up
			for (long l = 0; l < 100; l++)
			{
				auto res = filterTransform(MyOddLmbda, testVec, trueLambda, falseLambda);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					auto res = filterTransform(MyOddLmbda, testVec, trueLambda, falseLambda);
				}
			}
		}
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};


	auto run_res_innerProd = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, for_loop_run, TEST_LOOP_SZ);
	auto stats_inner_prod = performanceStats(run_res_innerProd.m_raw_results);

	//auto run_transformBranchy = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, stl_transform_with_branch_run, TEST_LOOP_SZ);
	//auto stats_transform_branchy = performanceStats(run_transformBranchy.m_raw_results);

	auto dr3_raw_resultsFilter = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_filterTransdform, TEST_LOOP_SZ);
	auto stats_DR3_filter = performanceStats(dr3_raw_resultsFilter.m_raw_results);

	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_select, TEST_LOOP_SZ);
	auto stats_DR3 = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_inner_prod)
	{
		//auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		//auto  valStl = run_res_innerProd.m_calc_results[elem.first];
		auto strMatch = "TO DO";//"valuesAreEqual(valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "for loop binarySelectionBetweenSimpleFunctions , size " << elem.first << " , " << elem.second.first << " +- " << elem.second.second << "\t \t DR3 filter_transform medium weight  , size " << elem.first << " , " << stats_DR3_filter[elem.first].first << " +- " << stats_DR3_filter[elem.first].second << "\t \t DR3 binarySelection medium Weight , size " << elem.first << " , " << stats_DR3[elem.first].first << " +- " << stats_DR3[elem.first].second << ", numerical check: " << strMatch << "\n";
	}



}


//heavy weight one odd lambda  functions
void testBinarySelection3()
{
	using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
	const int TEST_LOOP_SZ = 1000;

	for (long SZ = 200; SZ < 20000; SZ += 200)
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


void binarySelectionBetweenHeavyWeightFunction()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	auto zero = InstructionTraits<VecXX::INS>::nullValue;
	auto one = VecXX::scalar(1.0);
	auto two = VecXX::scalar(2.0);
	auto half = VecXX::scalar(0.5);


	auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };
	using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
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

	/*
	auto falseLambda = [&](auto q)
	{
		auto X = (((((a[1] * q + a[2]) * q + a[3]) * q + a[4]) * q + a[5]) * q + a[6]) /
			((((b[1] * q + b[2]) * q + b[3]) * q + b[4]) * q + VecXX::scalar(1.0));
		return X;
	};
	*/

	auto falseLambdaS = [](auto x) { return exp(-sin(x / VecXX::scalar(20.))) + VecXX::scalar(1.0) / (VecXX::scalar(1.0) + exp(x * (sin(x * VecXX::scalar(3.7) + VecXX::scalar(12.2) / x))));  };


	getRandomShuffledVector(-1); // reset  random input vectors

	auto for_loop_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		auto C = v1;
		//use auto on scalars so we can switch between float and double instruction sets too


		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				for (int k = 0; k < VEC_SZ; k++)
				{
					auto x = v1[k];
					if ((x - two * floor(x * half)) >= one)
					{
						C[k] = trueLambda(x);
					}
					else
					{
						C[k] = falseLambdaS(x);
					}
				}
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					for (int k = 0; k < VEC_SZ; k++)
					{
						auto x = v1[k];
						if ((x - two * floor(x * half)) >= one)
						{
							C[k] = trueLambda(x);
						}
						else
						{
							C[k] = falseLambdaS(x);
						}
					}
				}
			}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};



	auto stl_transform_with_branch_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		auto C = v1;
		//use auto on scalars so we can switch between float and double instruction sets too


		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				std::transform(v1.begin(), v1.end(), C.begin(),
					[&](auto x)
					{
						return MyOddLmbda(x) ? trueLambda(x) : falseLambdaS(x);
					}
				);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					std::transform(v1.begin(), v1.end(), C.begin(),
						[&](auto x)
						{
							return MyOddLmbda(x) ? trueLambda(x) : falseLambdaS(x);
						}
					);
				}
			}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};




	auto DR3_select = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		auto one = VecXX::scalar(1.0);
		auto two = VecXX::scalar(2.0);
		// making sure scalar constants are coerced to type of  float associated with chosen instruction set for ease opf switching sets
		//auto trueLambda = [](auto x) { return  (x * VecXX::scalar(3.7) + VecXX::scalar(16.2));  };
		//auto falseLambda = [](auto x) { return  x * ((x * VecXX::scalar(3.5)) + VecXX::scalar(12.2));  };

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX testVec(v1);


		{
			auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

			//warm up
			for (long l = 0; l < 100; l++)
			{
				auto res = selectTransform(MyOddLmbda, testVec, trueLambda, falseLambdaS);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					auto res = selectTransform(MyOddLmbda, testVec, trueLambda, falseLambdaS);
				}
			}
		}
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};



	auto DR3_filterTransdform = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		auto one = VecXX::scalar(1.0);
		auto two = VecXX::scalar(2.0);
		// making sure scalar constants are coerced to type of  float associated with chosen instruction set for ease opf switching sets
		//auto trueLambda = [](auto x) { return  (x * VecXX::scalar(3.7) + VecXX::scalar(16.2));  };
		//auto falseLambda = [](auto x) { return  x * ((x * VecXX::scalar(3.5)) + VecXX::scalar(12.2));  };

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX testVec(v1);


		{
			auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

			//warm up
			for (long l = 0; l < 100; l++)
			{
				auto res = filterTransform(MyOddLmbda, testVec, trueLambda, falseLambdaS);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					auto res = filterTransform(MyOddLmbda, testVec, trueLambda, falseLambdaS);
				}
			}
		}
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};


	auto run_res_innerProd = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, for_loop_run, TEST_LOOP_SZ);
	auto stats_inner_prod = performanceStats(run_res_innerProd.m_raw_results);

	//auto run_transformBranchy = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, stl_transform_with_branch_run, TEST_LOOP_SZ);
	//auto stats_transform_branchy = performanceStats(run_transformBranchy.m_raw_results);

	auto dr3_raw_resultsFilter = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_filterTransdform, TEST_LOOP_SZ);
	auto stats_DR3_filter = performanceStats(dr3_raw_resultsFilter.m_raw_results);

	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_select, TEST_LOOP_SZ);
	auto stats_DR3 = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_inner_prod)
	{
		//auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		//auto  valStl = run_res_innerProd.m_calc_results[elem.first];
		auto strMatch = "TO DO";//"valuesAreEqual(valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "for loop binarySelectionBetweenSimpleAndHeavyFunctions , size " << elem.first << " , " << elem.second.first << " +- " << elem.second.second << "\t \t DR3 filter_transform heavy weight  , size " << elem.first << " , " << stats_DR3_filter[elem.first].first << " +- " << stats_DR3_filter[elem.first].second << "\t \t DR3 binarySelection heavy Weight , size " << elem.first << " , " << stats_DR3[elem.first].first << " +- " << stats_DR3[elem.first].second << ", numerical check: " << strMatch << "\n";
	}



}




/*
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

*/


void doCountIf()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	auto zero = InstructionTraits<VecXX::INS>::nullValue;
	auto one = VecXX::scalar(1.0);
	auto two = VecXX::scalar(2.0);
	auto half = VecXX::scalar(0.5);


	auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };
	using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
	/// from acklams inverse cdf normal



	getRandomShuffledVector(-1); // reset  random input vectors

	auto std_count_if = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		auto C = v1;
		//use auto on scalars so we can switch between float and double instruction sets too


		auto MyOddLmbda = [&](auto x) { return  (x - two * floor(x * half)) >= one;  };

		volatile auto resStl = 0.0;

		{

			auto isOdd = [](auto x) {
				auto y = (x - (2. * floor(x * 0.5)));
				if (y >= 1.0)
					return false;
				else return true;
			};

			//warm up
			for (long l = 0; l < 100; l++)
			{
				resStl = std::count_if(begin(v1), end(v1), isOdd);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					resStl = std::count_if(begin(v1), end(v1), isOdd);
				}
			}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};





	auto DR3_count_if = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		auto one = VecXX::scalar(1.0);
		auto two = VecXX::scalar(2.0);
		// making sure scalar constants are coerced to type of  float associated with chosen instruction set for ease opf switching sets
		//auto trueLambda = [](auto x) { return  (x * VecXX::scalar(3.7) + VecXX::scalar(16.2));  };
		//auto falseLambda = [](auto x) { return  x * ((x * VecXX::scalar(3.5)) + VecXX::scalar(12.2));  };

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX test(v1);

		auto Add = [](auto x, auto y) { return  x + y;  };
		

		{

			auto oneIfOddLmbda = [](auto x) { return  iff((x - (VecXX::INS(2.0) * floor(x * VecXX::INS(0.5)))) >= VecXX::INS(1.0), VecXX::INS(1.0), VecXX::INS(0.0)); };
			auto Add = [](auto x, auto y) { return  x + y;  };
			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = transformReduce(test, oneIfOddLmbda, Add);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = transformReduce(test, oneIfOddLmbda, Add);
				}
			}
		}
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};


	auto stl_run_res = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, std_count_if, TEST_LOOP_SZ);
	auto stats_stl_count = performanceStats(stl_run_res.m_raw_results);

	//auto run_transformBranchy = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, stl_transform_with_branch_run, TEST_LOOP_SZ);
	//auto stats_transform_branchy = performanceStats(run_transformBranchy.m_raw_results);

	//auto dr3_raw_resultsFilter = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_filterTransdform, TEST_LOOP_SZ);
	//auto stats_DR3_filter = performanceStats(dr3_raw_resultsFilter.m_raw_results);

	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_count_if, TEST_LOOP_SZ);
	auto stats_DR3 = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_stl_count)
	{
		//auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		//auto  valStl = run_res_innerProd.m_calc_results[elem.first];
		auto strMatch = "TO DO";//"valuesAreEqual(valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "std count if , size " << elem.first << " , " << elem.second.first << " +- " << elem.second.second <<  "\t \t DR3 count if using transform reduce , size " << elem.first << " , " << stats_DR3[elem.first].first << " +- " << stats_DR3[elem.first].second << ", numerical check: " << strMatch << "\n";
	}



}







