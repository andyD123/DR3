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
#include <cstring>


#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/accumulate_transform.h"
#include "../Vectorisation/VecX/error_utils.h"


#include "../Vectorisation/VecX/zip_utils.h"

#include "norm.h"


/*
Uncomment one of the Using namespace lines below to select the instruction set that you wish to run  
Those ending in F have float type as underlying, those ending with D have a double.

The project is set to compile using the AVX512  enhanced instruction set. The namespace selection 
choses the type of the intrinsics that are used to instantiate lambdas.

If your hardware does not support AVX512 chose the next level down AVX2 and avoid using namespaces 
DRC::VecD8D or DRC::VecF16F which will cause generation of code with instructions that your computer doesn't support. 

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
//using namespace DRC::VecD4D;	//avx2   double
//using namespace DRC::VecF8F;	// avx2  float
using namespace DRC::VecD8D;  //avx512 double
//using namespace DRC::VecF16F; //avx512   float


using FLOAT = typename InstructionTraits<VecXX::INS>::FloatType;


//AllAllocatorsGuard<FLOAT> allocGuard;

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


double getErr(const  std::vector<double>& t) 
{
	ignore(t);
	return 2e-11;
}


float getErr(const  std::vector<float>& t)
{
	ignore(t);
	return 1.0;
}


template<typename T>
bool vectorsEqualD(const std::vector<T>& C1, const std::vector<T>& C2, const std::vector<T>& C3, const std::vector<T>& input,  T ERR = 1e-13 )
{

	
	ERR = getErr(C1);

	bool  testOK = true;
	
	if (C1.size() != C2.size())
	{
		std::cout << "wrong size C1,C2"<< C1.size() << ", " << C2.size() <<std::endl;
		return false;
	}

	if (C3.size() != C2.size())
	{
		std::cout << "wrong size C2,C3" << C2.size() << ", " << C3.size() << std::endl;
		return false;
	}

	for (size_t i = 0; i < C3.size(); i++)
	{
		auto err1 = fabs((C1[i] - C2[i]) / (fabs(C2[i]) + fabs(C1[i])));
		auto err2 = fabs((C1[i] - C3[i]) / (fabs(C1[i])+ fabs(C3[i])));

		if ((err1 > ERR) || (err2 > ERR))
		{
			testOK = false;
			std::cout << "\n err diff@ " << i << " err1 =" << err1 << ", err2 = " << err2 << "\n";
			std::cout << "\n val @ " << i << " C1[i] =" << C1[i] << ", C2[i] = " << C2[i] << ", C3[i] = " << C3[i] << "input val=" << input[i] <<"\n";
			std::cout << std::endl;
			break;
		}
	}

	return testOK;

}

bool valuesAreEqual(double x, double y,double tol =  1e-14)
{

	auto err1 = fabs((x - y) / (fabs(x) + fabs(y)));

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
		for (int i = 0; i < SZ; i++) { v[i] += /*FloatType(SZ / 2)*/ + i; }
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

using Calc_Values = std::map<int, FLOAT>;
using Calc_Values_V = std::map<int,std::vector<FLOAT> >;
using  Mapped_Performance_Results = std::map<int, std::vector<double> >; // array size  v vector<throughput for runs>
using Mapped_Stats = std::map<int, std::pair<double, double> >; // size -.pair ( throughput ,  std dev of through put)

struct RunResults
{
	Mapped_Performance_Results m_raw_results;
	Calc_Values  m_calc_results;
};

struct RunResultsVec
{
	Mapped_Performance_Results m_raw_results;
	Calc_Values_V  m_calc_results;
};

class TimerGuard
{
	double& m_runTime;
	std::chrono::steady_clock::time_point  m_startTme;

public:
	TimerGuard(double& runTime) : m_runTime(runTime), m_startTme(std::chrono::steady_clock::now()) { runTime = 0.; }

	~TimerGuard()
	{
		auto endTime = std::chrono::steady_clock::now();
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
				results.m_calc_results[VEC_SZ] = static_cast<FLOAT>(calc_value);
			}
		}
	}
	return results;
};


auto runFunctionOverDifferentSizeVec = [](int testRepeats, int vec_start_size, int vec_stepSZ, int vec_maxSZ, const auto& func, long testLoopSZ)
{

	RunResultsVec results;

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
				std::vector<FLOAT> tmp = res.first;
				results.m_calc_results[VEC_SZ] = tmp;
			}
			
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
void    doInnerProd();
void	doTransform();
void 	doSumSqrs();
void    khanAccumulation();
void	binarySelectionBetweenConst();
void	binarySelectionBetweenLinearFunction(); // y= mx + c    a couple of of operations
void    binarySelectionBetweenMiddleWeightFunction();
void	binarySelectionBetweenHeavyWeightFunction();
void	doCountIf(); 
void    doMinMax();
void    doSumSqrsMulti();



void doAVXMax512Dance();



int main()
{


	std::cout << "\n \n \n \n testMemCpy2() \n" << std::endl;
	testMemCpy2(); 

	return 0;

	//accumulate 
	std::cout << "\n \n \n \n doMax() \n"  << std::endl;
	doMax();

	//multi reduce
	std::cout << "\n \n \n \n doMinMax() \n" << std::endl;
	doMinMax();

//transform accum
	std::cout << "\n \n \n \n doInnerProd() \n" << std::endl;
	doInnerProd();


	std::cout << "\n \n \n \n doSumSqrs() \n" << std::endl;
	doSumSqrs();

	// lambda capture   sum with  error correction
	std::cout << "\n \n \n \n khanAccumulation() \n" << std::endl;
	khanAccumulation();
	
// branching
	std::cout << "\n \n \n \n binarySelectionBetweenConst() \n" << std::endl;
	binarySelectionBetweenConst(); //select between constants


	std::cout << "\n \n \n \n binarySelectionBetweenLinearFunction() \n" << std::endl;
	binarySelectionBetweenLinearFunction();

	std::cout << "\n \n \n \n binarySelectionBetweenMiddleWeightFunction() \n" << std::endl;
	binarySelectionBetweenMiddleWeightFunction();

	std::cout << "\n \n \n \n binarySelectionBetweenHeavyWeightFunction() \n" << std::endl;
	binarySelectionBetweenHeavyWeightFunction();
	
	std::cout << "\n \n \n \n doCountIf() \n" << std::endl;
	doCountIf();

	// experimental
	doMinMax();
	doSumSqrsMulti();
	//return 0;

	
	// use namespace DRC::VecD8D  run this and watch power consumption
	// switches between AVX2 and AVX512  implementations
	// AVX512 uses less energy in this case
	// doAVXMax512Dance();

	return 0;
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
		FloatType* cp = nullptr;
		FloatType* vp = nullptr;
		auto bufferSz = SZ * sizeof(FloatType);


		loop = 500000;
		cp = new FloatType[SZ];
		vp = new FloatType[SZ];
		bufferSz = SZ * sizeof(FloatType);
		auto startTme = std::chrono::high_resolution_clock::now();

		for (long i = 0; i < loop; ++i)
		{
			volatile auto xx =std::memcpy(cp, vp, bufferSz);

			ignore(xx);
		}
			
		auto endTime = std::chrono::high_resolution_clock::now();
		auto runtime = endTime - startTme;
		std::cout << " , MemCpy run , " << SZ * sizeof(FloatType) * loop / (1000000000.0) * 1.0 / (runtime.count() / 1000000000.0) << " ,  GB per sec";	

		delete[] cp;
		delete[] vp;



		startTme = std::chrono::high_resolution_clock::now();
		auto cpyLambda = [](const auto& rhs) { return rhs; };

		for (long i = 0; i < loop; ++i) 
		{
			auto res = transform( cpyLambda, test);
		}

		endTime = std::chrono::high_resolution_clock::now();
		runtime = endTime - startTme;
		std::cout << "			Unitary  run  , " << SZ * sizeof(FloatType) * loop / (1000000000.0) * 1.0 / (runtime.count() / 1000000000.0) << " ,  GB per sec";
	
	}
}


void doMax()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto accumulate_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		auto v1 = getRandomShuffledVector(VEC_SZ, 0);

		//warm up
		for (long l = 0; l < 100; l++)
		{
			res = *std::max_element(v1.begin(), v1.end());
		}


		{   TimerGuard timer(time);
		for (long l = 0; l < TEST_LOOP_SZ; l++)
		{
			res = *std::max_element(v1.begin(), v1.end());
		}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};



	auto DR3_accumulate = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		// generic lambda for max either calling a max instruction or doing a selection with iff
	   // auto mxDbl = [](auto lhs, auto rhs) { return max(lhs, rhs); };
		auto mxDbl = [](auto lhs, auto rhs) { return iff(lhs > rhs, lhs, rhs); }; //using iff fastest 

		auto v1 = getRandomShuffledVector(SZ, 0); // std stl vector double or float 
		VecXX vec(v1);

		//warm up
		for (long l = 0; l < 100; l++)
		{
			res = reduce(vec, mxDbl);
		}


		{   TimerGuard timer(time);
		for (long l = 0; l < TEST_LOOP_SZ; l++)
		{
			res = reduce(vec, mxDbl);
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
		std::cout << "  std::max_element, size " << perf_stl.first << ", " << perf_stl.second.first << ", + - ," << perf_stl.second.second << "\t \t DR3 reduce, size " << perf_stl.first << ", " << stats_DR3_perf[perf_stl.first].first << ",  + - ," << stats_DR3_perf[perf_stl.first].second << ", numerical check : " << strMatch << "\n";
	}
}





//example doing  multi reduction operation after load
//should be faster 
void doMinMax()
{

	const long TEST_LOOP_SZ = 100;// 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto accumulate_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		volatile  double res_min = 0.;
		auto v1 = getRandomShuffledVector(VEC_SZ, 0);

		//warm up
		for (long l = 0; l < 100; l++)
		{
			res = *std::max_element(v1.begin(), v1.end());

			res_min = *std::min_element(v1.begin(), v1.end());
		}


		{   TimerGuard timer(time);
		for (long l = 0; l < TEST_LOOP_SZ; l++)
		{
			res = *std::max_element(v1.begin(), v1.end());
			res_min = *std::min_element(v1.begin(), v1.end());
		}
		}
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);

		ignore(res_min);
	};



	auto DR3_accumulate = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		// generic lambda for max either calling a max instruction or doing a selection with iff
	   // auto mxDbl = [](auto lhs, auto rhs) { return max(lhs, rhs); };
		auto mxDbl = [](auto lhs, auto rhs) { return iff(lhs > rhs, lhs, rhs); }; //using iff fastest 
		auto minDbl = [](auto lhs, auto rhs) { return iff(lhs < rhs, lhs, rhs); }; //using iff fastest 

		auto v1 = getRandomShuffledVector(SZ, 0); // std stl vector double or float 
		VecXX vec(v1);

		auto ress = reduceM(vec, mxDbl, minDbl);
		ignore(ress);
		//warm up
		for (long l = 0; l < 100; l++)
		{
			 res = reduce(vec, mxDbl);
			 double mnn = reduce(vec, minDbl);
			 ignore(mnn);
		}


		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
		  		 res = reduce(vec, mxDbl);
	   			 double mnn = reduce(vec, minDbl);
				 ignore(mnn);
			}

		}

		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};



	auto DR3_accumulate_multi = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		// generic lambda for max either calling a max instruction or doing a selection with iff
	   // auto mxDbl = [](auto lhs, auto rhs) { return max(lhs, rhs); };
		auto mxDbl = [](auto lhs, auto rhs) { return iff(lhs > rhs, lhs, rhs); }; //using iff fastest 
		auto minDbl = [](auto lhs, auto rhs) { return iff(lhs < rhs, lhs, rhs); }; //using iff fastest 

		auto v1 = getRandomShuffledVector(SZ, 0); // std stl vector double or float 
		VecXX vec(v1);

		auto ress = reduceM(vec, mxDbl, minDbl);
		//warm up
		for (long l = 0; l < 100; l++)
		{
			ress = reduceM(vec, mxDbl, minDbl);
		}


		{   TimerGuard timer(time);
		for (long l = 0; l < TEST_LOOP_SZ; l++)
		{
			ress = reduceM(vec, mxDbl, minDbl);
			double mmmm = std::get<0>(ress);
			double mnn = std::get<1>(ress);
			ignore(mnn);
			ignore(mmmm);
		}
		}

		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};


	auto run_res_stl = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, accumulate_run, TEST_LOOP_SZ);
	auto stats_stl = performanceStats(run_res_stl.m_raw_results);


	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_accumulate, TEST_LOOP_SZ);
	auto stats_DR3_perf = performanceStats(dr3_raw_results.m_raw_results);

	auto dr3_raw_results_mult = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_accumulate_multi, TEST_LOOP_SZ);
	auto stats_DR3_perf_mult = performanceStats(dr3_raw_results_mult.m_raw_results);


	//print out results
	for (const auto& perf_stl : stats_stl)
	{
		auto  valDr3 = dr3_raw_results.m_calc_results[perf_stl.first];
		auto  valStl = run_res_stl.m_calc_results[perf_stl.first];
		auto  valDr3Mult = dr3_raw_results_mult.m_calc_results[perf_stl.first];
		auto strMatch = valuesAreEqual(valDr3, valStl, valDr3Mult) ? "calcs match" : "cal difference";
		std::cout << "  std::max_element, size " << perf_stl.first << ", " << perf_stl.second.first << ", + - ," << perf_stl.second.second 
			<< "\t \t DR3 reduce, size " << perf_stl.first << ", " << stats_DR3_perf[perf_stl.first].first << ",  + - ," << stats_DR3_perf[perf_stl.first].second 
			<< "\t \t DR3 reduce_mult, size " << perf_stl.first << ", " << stats_DR3_perf_mult[perf_stl.first].first << ",  + - ," << stats_DR3_perf_mult[perf_stl.first].second
			<< ", numerical check : " << strMatch << "\n";
	}


	
}






void doTransform() 
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 400;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto stl_transform = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		auto v = getRandomShuffledVector(VEC_SZ, 0);
		auto targetVec = v;
		volatile auto transformVec = v;
		VecXX test(v);
		auto SQR = [](auto rhs) { return rhs * rhs; };
		
		//warm up
		for (long l = 0; l < 100; l++)
		{
			std::transform(v.begin(), v.end(), targetVec.begin(), SQR);
		}

		TimerGuard timer(time);
		{
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				std::transform(v.begin(), v.end(), targetVec.begin(), SQR);

			}
		}
		
		return  std::make_pair(targetVec, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};



	auto DR3_transform = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		auto v = getRandomShuffledVector(VEC_SZ, 0);
		auto targetVec = v;
		VecXX res = v;
		const VecXX test(v);
		auto SQR = [](auto rhs) { return rhs * rhs; };
		
		//warm up
		for (long l = 0; l < 100; l++)
		{
			res = transform(SQR, test);
		}

		TimerGuard timer(time);
		{
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				res = transform(SQR, test);
			}
		}
		
		std::vector<FLOAT> vvv = res;
		return std::make_pair(vvv, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto dr3_raw_results = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_transform, TEST_LOOP_SZ);
	auto stats_DR3 = performanceStats(dr3_raw_results.m_raw_results);

	auto run_res_stl = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, stl_transform, TEST_LOOP_SZ);
	auto stats_stl = performanceStats(run_res_stl.m_raw_results);


	//print out results

	for (const auto& elem : stats_stl)
	{
		auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = run_res_stl.m_calc_results[elem.first];
		auto strMatch =  vectorsEqual(valDr3, valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "STL transform , size " << elem.first << " , " << stats_stl[elem.first].first << ", +- ," << stats_stl[elem.first].second << "\t \t DR3 transform , size " << elem.first << " , " << stats_DR3[elem.first].first << ", +- ," << stats_DR3[elem.first].second << ", numerical check: " << strMatch << "\n";
	}
	
}






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

		auto v1 = getRandomShuffledVector(VEC_SZ,0); 
		auto v2 = getRandomShuffledVector(VEC_SZ,1); 
		
		//warm up
		for (long l = 0; l < 100; l++)
		{
			res = inner_product(v1.cbegin(), v1.cend(), v2.cbegin(), zero);
		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				res = inner_product(v1.cbegin(), v1.cend(), v2.cbegin(), zero);
			}
		}
		
		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_inner_prod = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(SZ,0); 
		auto v2 = getRandomShuffledVector(SZ,1); 
		VecXX t1(v1);
		VecXX t2(v2);

		
		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
		auto Mult = [](auto X, auto Y) { return X * Y; };

		//warm up
		for (long l = 0; l < 100; l++)
		{
			res = transformReduce(t1, t2, Mult, Sum);
		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				res = transformReduce(t1, t2, Mult, Sum);
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
		std::cout << "STL inner product , size " <<  elem.first << " , " << elem.second.first << ", +- ," << elem.second.second << "\t \t DR3 inner product , size " << elem.first << " , " << stats_DR3_inner_prod[elem.first].first << ", +- ," << stats_DR3_inner_prod[elem.first].second << ", numerical check: " << strMatch << "\n";
	}

	

}


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

		auto v1 = getRandomShuffledVector(VEC_SZ); 

		//warm up
		for (long l = 0; l < 100; l++)
		{
			res = inner_product(v1.cbegin(), v1.cend(), v1.cbegin(), zero);
		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				res = inner_product(v1.cbegin(), v1.cend(), v1.cbegin(), zero);
			}
		}

		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_inner_prod = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(SZ); 
		VecXX t1(v1);

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
		auto Mult = [](auto X, auto Y) { return X * Y; };

		//warm up
		for (long l = 0; l < 100; l++)
		{
			res = transformReduce(t1, t1, Mult, Sum);
		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				res = transformReduce(t1, t1, Mult, Sum);
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
		std::cout << "STL inner product sum sqrs , size " << elem.first << " , " << elem.second.first << ", +- ," << elem.second.second << "\t \t DR3 inner product sum sqrs , size " << elem.first << " , " << stats_DR3_inner_prod[elem.first].first << ", +- ," << stats_DR3_inner_prod[elem.first].second << ", numerical check: " << strMatch << "\n";
	}



}



/*
applies multiple transforms and reductions after the load operation 
reducing memorry traversal

here we do squaring an just return a copy as transform
and sum as the reduction, usefuk for stats to get sum of all values
and sum of all squares of values
*/
void doSumSqrsMulti()
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
		volatile  double res2 = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ);

		//warm up
		for (long l = 0; l < 100; l++)
		{
			res = inner_product(v1.cbegin(), v1.cend(), v1.cbegin(), zero);
		}


		{   TimerGuard timer(time);
		for (long l = 0; l < TEST_LOOP_SZ; l++)
		{
			res = inner_product(v1.cbegin(), v1.cend(), v1.cbegin(), zero);

			res2 =std::accumulate(v1.cbegin(), v1.cend(),  zero);
		}
		}

		return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);

		ignore(res2);
	};

	auto DR3_inner_prod = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		volatile  double res2 = 0.;

		auto v1 = getRandomShuffledVector(SZ);
		VecXX t1(v1);

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
		//auto Mult = [](auto X, auto Y) { return X * Y; };
		auto Unit = [](auto X) { return X; };
		auto SQR = [](auto X) { return X*X; };

		//warm up
		for (long l = 0; l < 100; l++)
		{

			res = transformReduce(t1, SQR, Sum);
			res2 = transformReduce(t1, Unit, Sum);

		}


		{   TimerGuard timer(time);
		for (long l = 0; l < TEST_LOOP_SZ; l++)
		{
		 	res = transformReduce( t1, SQR, Sum);
			res2 = transformReduce(t1, Unit, Sum);	
		}
		}

		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

		ignore(res2);

	};

	auto DR3_inner_prod_multi_transformReduce = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		volatile  double res2 = 0.;

		auto v1 = getRandomShuffledVector(SZ);
		VecXX t1(v1);

		auto Sum = [](auto lhs, auto rhs) { return lhs + rhs; };
		//auto Mult = [](auto X, auto Y) { return X * Y; };
		auto Unit = [](auto X) { return X; };
		auto SQR = [](auto X) { return X * X; };

		//warm up
		for (long l = 0; l < 100; l++)
		{
	
			auto tpl = transformReduceM(t1, Unit, Sum, SQR, Sum);
			//auto total 
			res2 = std::get<0>(tpl);
			//auto total_sqr
			res = std::get<1>(tpl);
		}


		{   TimerGuard timer(time);
		for (long l = 0; l < TEST_LOOP_SZ; l++)
		{

			auto tpl = transformReduceM(t1, Unit, Sum, SQR, Sum);
			res2 = std::get<0>(tpl);
			auto res0 = std::get<1>(tpl);
			ignore(res0);
	
		}
		}

		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

		ignore(res2);

	};





	auto run_res_innerProd = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, inner_prod_run, TEST_LOOP_SZ);
	auto stats_inner_prod = performanceStats(run_res_innerProd.m_raw_results);


	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_inner_prod, TEST_LOOP_SZ);
	auto stats_DR3_inner_prod = performanceStats(dr3_raw_results.m_raw_results);


	auto dr3_raw_results_multiReduce = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_inner_prod_multi_transformReduce, TEST_LOOP_SZ);
	auto stats_DR3_inner_prod_multiReduce = performanceStats(dr3_raw_results_multiReduce.m_raw_results);


	//print out results
	for (const auto& elem : stats_inner_prod)
	{
		auto  valDr3 =  dr3_raw_results.m_calc_results[elem.first];
		auto  valDr3_multi =  dr3_raw_results_multiReduce.m_calc_results[elem.first];
		auto  valStl = run_res_innerProd.m_calc_results[elem.first];
		auto strMatch = valuesAreEqual(valDr3, valStl, valDr3_multi) ? "calcs match" : "cal difference";
		std::cout << "STL inner product sum sqrs , size " << elem.first << " , " << elem.second.first << ", +- ," << elem.second.second 
			<< "\t \t DR3 inner product sum sqrs , size " << elem.first << " , " << stats_DR3_inner_prod[elem.first].first << ", +- ," << stats_DR3_inner_prod[elem.first].second 
			<< "\t \t DR3 Multi reduce , size " << elem.first << " , " << stats_DR3_inner_prod_multiReduce[elem.first].first << ", +- ," << stats_DR3_inner_prod_multiReduce[elem.first].second
			<< ", numerical check: " << strMatch << "\n";
	}



}



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
		
		//warm up
		for (long l = 0; l < 100; l++)
		{
				
			res = std::accumulate(begin(v1), end(v1), zero);
		}

	
		{	TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				res = std::accumulate(begin(v1), end(v1), zero);
			}
		}

		relative_error = (VEC_SZ / 6. - res) / res;
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

		//warm up
		for (long l = 0; l < 100; l++)
		{
			res = reduce(t1, KhanAddV, zero);

		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				res = reduce(t1, KhanAddV, zero);
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
		std::cout << "STL accumulate , size " << elem.first << " , " << elem.second.first << ", +- ," << elem.second.second <<"relative err"<< valStl << "\t \t DR3 khan accumulate , size " << elem.first << " , " << stats_DR3_inner_prod[elem.first].first << ", +- ," << stats_DR3_inner_prod[elem.first].second << ", relative err " << valDr3 << "\n";
	}



}




void binarySelectionBetweenConst()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	auto one = VecXX::scalar(1.0);
	auto two = VecXX::scalar(2.0);
	auto truVal = one;
	auto falseVal = two;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto for_loop_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		//volatile  double res = 0.;

		auto v1 = getRandomShuffledVector(VEC_SZ); 
		auto C = v1;
		double halfSize = VEC_SZ * 0.5;	
		//warm up
		for (long l = 0; l < 100; l++)
		{
			for (int k = 0; k < VEC_SZ; k++)
			{
				auto x = v1[k];
				if (x > halfSize) 
				{
					C[k] = one;
				}
				else
				{
					C[k] = two;
				}
			}
		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				for (int k = 0; k < VEC_SZ; k++)
				{
					auto x = v1[k];
					if (x > halfSize) 
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
		
		return  std::make_pair(C, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_select = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX testVec(v1);
		VecXX res;

		FLOAT halfSize = SZ * static_cast<FLOAT>(0.5); 
		auto upperHalfLmbda = [&](auto x) { return x > halfSize; };
		//warm up
		for (long l = 0; l < 100; l++)
		{
			auto res1 = select(upperHalfLmbda, testVec, truVal, falseVal);
		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				res = select(upperHalfLmbda, testVec, truVal, falseVal);
			}
		}
		
		std::vector<FLOAT> vec = res;
		return std::make_pair(vec, numOps(TEST_LOOP_SZ, SZ) / time);

	};



	auto run_res_for_loop = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, for_loop_run, TEST_LOOP_SZ);
	auto stats_loop = performanceStats(run_res_for_loop.m_raw_results);


	auto dr3_raw_results = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_select, TEST_LOOP_SZ);
	auto stats_DR3_inner_prod = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_loop)
	{
		auto  valDr3_select = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = run_res_for_loop.m_calc_results[elem.first];
		bool VecsOK = vectorsEqualD(valDr3_select, valStl, valStl, valStl);
		auto strMatch = VecsOK ? "calcs match" : "cal difference";
		std::cout << "for loop binarySelectionBetweenConst , size " << elem.first << " , " << elem.second.first << ", +- ," << elem.second.second << "\t \t DR3 binarySelectionBetweenConst , size " << elem.first << " , " << stats_DR3_inner_prod[elem.first].first << ", +- , " << stats_DR3_inner_prod[elem.first].second << ", numerical check: " << strMatch << "\n";
	}



}


void binarySelectionBetweenLinearFunction()
{

	const long TEST_LOOP_SZ = 100;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto for_loop_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		auto C = v1;
		double halfSize = VEC_SZ * 0.5;
		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				for (int k = 0; k < VEC_SZ; k++)
				{
					auto x = v1[k];
					if (x > halfSize)
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
						if (x > halfSize) 
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
		return  std::make_pair(C, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto stl_transform_with_branch_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		auto C = v1;
		double halfSize = VEC_SZ * 0.5;
		auto MyOddLmbda = [&](auto x) { return x > halfSize; };

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
		return  std::make_pair(C, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_select = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		// making sure scalar constants are coerced to type of  float associated with chosen instruction set for ease opf switching sets
		auto trueLambda = [](auto x) { return  (x * VecXX::scalar(3.7) + VecXX::scalar(16.2));  };
		auto falseLambda = [](auto x) { return  x * ((x * VecXX::scalar(3.5)) + VecXX::scalar(12.2));  };

		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX testVec(v1);
		VecXX res;
		FLOAT halfSize = SZ * static_cast<FLOAT>(0.5);
		auto halfOfLmbda = [&](auto x) { return x > halfSize; };

		//warm up
		for (long l = 0; l < 100; l++)
		{
				res = selectTransform(halfOfLmbda, testVec, trueLambda, falseLambda);
		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
					res = selectTransform(halfOfLmbda, testVec, trueLambda, falseLambda);
			}
		}
		
		std::vector<FLOAT> v = res;
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};



	auto run_res_for_loop = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, for_loop_run, TEST_LOOP_SZ);
	auto stats_for_loop = performanceStats(run_res_for_loop.m_raw_results);

	auto run_transformBranchy = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, stl_transform_with_branch_run, TEST_LOOP_SZ);
	auto stats_transform_branchy = performanceStats(run_transformBranchy.m_raw_results);

	auto dr3_raw_results = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_select, TEST_LOOP_SZ);
	auto stats_DR3 = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_for_loop)
	{
		auto  valDr3_select = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = run_res_for_loop.m_calc_results[elem.first];
		auto  val_branchy = run_transformBranchy.m_calc_results[elem.first];
		bool VecsOK = vectorsEqualD(valDr3_select, valStl, val_branchy, valStl);
		auto strMatch = VecsOK ? "calcs match" : "cal difference";
		std::cout << "for loop binarySelectionBetweenSimpleFunctions , size " << elem.first << " , " << elem.second.first << ", +- ," << elem.second.second << "\t \t stl transform branchy lambda  , size " << elem.first << " , " << stats_transform_branchy[elem.first].first << ", +- ," << stats_transform_branchy[elem.first].second << "\t \t DR3 binarySelectionBetweenVSimpleLambda , size " << elem.first << " , " << stats_DR3[elem.first].first << ", +- ," << stats_DR3[elem.first].second << ", numerical check: " << strMatch << "\n";
	}



}


void binarySelectionBetweenMiddleWeightFunction()
{

	const long TEST_LOOP_SZ = 100;// 0;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

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

		auto v1 = getRandomShuffledVector(VEC_SZ); 
		auto C = v1;
		double halfSize = VEC_SZ * 0.5;
		{
			//warm up
			for (long l = 0; l < 100; l++)
			{
				for (int k = 0; k < VEC_SZ; k++)
				{
					auto x = v1[k];
					if (x > halfSize ) 
					{
						C[k] = trueLambda(x);
					}
					else
					{
						C[k] = falseLambda(x);
					}
				}
			}

			
			{   TimerGuard timer(time);
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					for (int k = 0; k < VEC_SZ; k++)
					{
						auto x = v1[k];
						if (x > halfSize ) 
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
		return  std::make_pair(C, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_select = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		auto v1 = getRandomShuffledVector(SZ);
		VecXX testVec(v1);
		VecXX res;	
		FLOAT halfSize = SZ * static_cast<FLOAT>(0.5);
		auto halfOfLmbda = [&](auto x) { return x > halfSize; };
			
		for (long l = 0; l < 100; l++)//warm up
		{
				res = selectTransform(halfOfLmbda, testVec, trueLambda, falseLambda);
		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
					res = selectTransform(halfOfLmbda, testVec, trueLambda, falseLambda);
			}
		}
		
		std::vector<FLOAT> v = res;
		return std::make_pair(v, numOps(TEST_LOOP_SZ, SZ) / time);

	};

	auto DR3_filterTransdform = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		auto v1 = getRandomShuffledVector(SZ); 
		VecXX testVec(v1);
		VecXX res;
		{
			FLOAT halfSize = SZ * static_cast<FLOAT>(0.5);
			auto MyOddLmbda = [&](auto x) { return x > halfSize; };

			//warm up
			for (long l = 0; l < 100; l++)
			{
				res = filterTransform(MyOddLmbda, testVec, trueLambda, falseLambda);
			}

			TimerGuard timer(time);
			{
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					res = filterTransform(MyOddLmbda, testVec, trueLambda, falseLambda);

				}
			}
		}
		std::vector<FLOAT> v = res;
		return std::make_pair(v, numOps(TEST_LOOP_SZ, SZ) / time);

	};

	auto dr3_raw_resultsFilter = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_filterTransdform, TEST_LOOP_SZ);
	auto stats_DR3_filter = performanceStats(dr3_raw_resultsFilter.m_raw_results);


	auto run_res_for_loop = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, for_loop_run, TEST_LOOP_SZ);
	auto stats_for_loop = performanceStats(run_res_for_loop.m_raw_results);


	auto dr3_raw_results = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_select, TEST_LOOP_SZ);
	auto stats_DR3 = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_for_loop)
	{
		auto  valDr3_select = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = run_res_for_loop.m_calc_results[elem.first];
		auto  valDr3_filterTransform = dr3_raw_resultsFilter.m_calc_results[elem.first];

		bool VecsOK =vectorsEqualD(valDr3_select, valDr3_filterTransform, valStl, valStl);
		auto strMatch = VecsOK ? "calcs match" : "cal difference";
		std::cout << "for loop binarySelectionBetweenMiddleWeightFunctions , size " << elem.first << " , " << elem.second.first << ", +-  " << elem.second.second << "\t \t DR3 filter_transform medium weight  , size " << elem.first << " , " << stats_DR3_filter[elem.first].first << ", +- ," << stats_DR3_filter[elem.first].second << "\t \t DR3 binarySelection medium Weight , size " << elem.first << " , " << stats_DR3[elem.first].first << ", +- ," << stats_DR3[elem.first].second << ", numerical check: " << strMatch << "\n";
	}



}


//heavy weight one odd lambda  functions

void binarySelectionBetweenHeavyWeightFunction()
{
	const int WARM_UP_LOOP = 10;
	const long TEST_LOOP_SZ = 10;// 00;
	const int repeatRuns =  20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;


	auto falseLambdaS = [](auto x) { return exp(-sin(x / VecXX::scalar(20.1))) + VecXX::scalar(1.1) / (VecXX::scalar(1.03) + exp(x * (sin(x * VecXX::scalar(3.9) + VecXX::scalar(12.4) / x))));  };
	auto trueLambdaS = [](auto x) { return exp(-sin(x / VecXX::scalar(20.))) + VecXX::scalar(1.0) / (VecXX::scalar(1.0) + exp(x * (sin(x * VecXX::scalar(3.7) + VecXX::scalar(12.2) / x))));  };


	getRandomShuffledVector(-1); // reset  random input vectors

	auto for_loop_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		auto v1 = getRandomShuffledVector(VEC_SZ); // std stl vector double or float 
		auto C = v1;
		double halfSize = VEC_SZ * 0.5;

		
		//warm up
		for (long l = 0; l < WARM_UP_LOOP; l++)
		{
			for (int k = 0; k < VEC_SZ; k++)
			{
				auto x = v1[k];
				if ( x > halfSize) 
				{
					C[k] = trueLambdaS(x);
				}
				else
				{
					C[k] = falseLambdaS(x);
				}
			}
		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				for (int k = 0; k < VEC_SZ; k++)
				{
					auto x = v1[k];
					if ( x > halfSize) 
					{
						C[k] = trueLambdaS(x);
					}
					else
					{
						C[k] = falseLambdaS(x);
					}
				}
			}
		}
		
		return  std::make_pair(C, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_select = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		VecXX  res;
		auto v1 = getRandomShuffledVector(SZ); 
		VecXX testVec(v1);
		FLOAT halfSize = SZ * static_cast<FLOAT>(0.5);
		auto halfOfLmbda = [&](auto x) { return x > halfSize; };

		//warm up
		for (long l = 0; l < WARM_UP_LOOP; l++)
		{
			res = selectTransform(halfOfLmbda, testVec, trueLambdaS, falseLambdaS);
		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				res = selectTransform(halfOfLmbda, testVec, trueLambdaS, falseLambdaS);
			}
		}
		
		std::vector<FLOAT> v = res;
		return std::make_pair(v, numOps(TEST_LOOP_SZ, SZ) / time);

	};

	auto DR3_filterTransdform = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		VecXX  res;
		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX testVec(v1);
		FLOAT halfSize = SZ * static_cast<FLOAT>(0.5);
		auto halfOfLmbda = [&](auto x) { return x > halfSize; };

		//warm up
		for (long l = 0; l < WARM_UP_LOOP; l++)
		{
			res = filterTransform(halfOfLmbda, testVec, trueLambdaS, falseLambdaS);
		}

		
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				res = filterTransform(halfOfLmbda, testVec, trueLambdaS, falseLambdaS);
			}
		}
		
		std::vector<FLOAT> v = res;
		return std::make_pair(v, numOps(TEST_LOOP_SZ, SZ) / time);

	};


	auto run_res_for_loop = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, for_loop_run, TEST_LOOP_SZ);
	auto stats_for_loop = performanceStats(run_res_for_loop.m_raw_results);

	auto dr3_raw_resultsFilter = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_filterTransdform, TEST_LOOP_SZ);
	auto stats_DR3_filter = performanceStats(dr3_raw_resultsFilter.m_raw_results);

	auto dr3_raw_results = runFunctionOverDifferentSizeVec(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_select, TEST_LOOP_SZ);
	auto stats_DR3 = performanceStats(dr3_raw_results.m_raw_results);


	//print out results
	for (const auto& elem : stats_for_loop)
	{
		auto  valDr3_select = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = run_res_for_loop.m_calc_results[elem.first];
		auto  valDr3_filterTransform = dr3_raw_resultsFilter.m_calc_results[elem.first];

		bool VecsOK = vectorsEqualD(valDr3_select, valStl, valDr3_filterTransform, valDr3_filterTransform,  getErr(valDr3_select));
		auto strMatch = VecsOK ? "calcs match" : "cal difference";

		std::cout << "for loop binarySelectionBetweenHeavyFunctions , size " << elem.first << " , " << elem.second.first << ", +- ," << elem.second.second << "\t \t DR3 filter_transform heavy weight  , size " << elem.first << " , " << stats_DR3_filter[elem.first].first << ", +- ," << stats_DR3_filter[elem.first].second << "\t \t DR3 binarySelection heavy Weight , size " << elem.first << " , " << stats_DR3[elem.first].first << " , +- , " << stats_DR3[elem.first].second << ", numerical check: " << strMatch << "\n";
	}
}



void doCountIf()
{

	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	getRandomShuffledVector(-1); // reset  random input vectors

	auto std_count_if = [&](int VEC_SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
	
		auto v1 = getRandomShuffledVector(VEC_SZ); 
		auto C = v1;
	
		volatile long resStl = 0;

		double halfSize = VEC_SZ * 0.5;
		auto isOverHalf= [&](auto x) { return x > halfSize; };

		for (long l = 0; l < 100; l++)
		{
			resStl = static_cast< long>( std::count_if(begin(v1), end(v1), isOverHalf) );
		}

			
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				resStl = static_cast<long>(std::count_if(begin(v1), end(v1), isOverHalf));
			}
		}
		
		return  std::make_pair((double) resStl, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
	};

	auto DR3_count_if = [&](int SZ, long TEST_LOOP_SZ)
	{
		double time = 0.;
		volatile  double res = 0.;
		auto v1 = getRandomShuffledVector(SZ); // std stl vector double or float 
		VecXX test(v1);

		//double halfSize = SZ * 0.5;
		FLOAT halfSize = SZ * static_cast<FLOAT>(0.5);

		auto oneIfOverHalfLambda = [&](auto x) 
		{
			return  iff((x > VecXX::INS(halfSize)),VecXX::INS(1.0),	VecXX::INS(0.0));
		};

		auto Add = [](auto x, auto y) { return  x + y;  };

		for (long l = 0; l < 10; l++)
		{
			res = transformReduce(test, oneIfOverHalfLambda, Add);
		}

			
		{   TimerGuard timer(time);
			for (long l = 0; l < TEST_LOOP_SZ; l++)
			{
				res = transformReduce(test, oneIfOverHalfLambda, Add);
			}
		}
		
		return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};



	auto stl_run_res = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, std_count_if, TEST_LOOP_SZ);
	auto stats_stl_count = performanceStats(stl_run_res.m_raw_results);

	auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_count_if, TEST_LOOP_SZ);
	auto stats_DR3 = performanceStats(dr3_raw_results.m_raw_results);



	//print out results
	for (const auto& elem : stats_stl_count)
	{
		auto  valDr3 = dr3_raw_results.m_calc_results[elem.first];
		auto  valStl = stl_run_res.m_calc_results[elem.first];
		auto strMatch = valuesAreEqual(valDr3, valStl) ? "calcs match" : "cal difference";
		std::cout << "std count if , size " << elem.first << " , " << elem.second.first << ", +-, " << elem.second.second <<  "\t \t DR3 count if using transform reduce , size " << elem.first << " , " << stats_DR3[elem.first].first << ", +- ," << stats_DR3[elem.first].second << ", numerical check: " << strMatch << "\n";
	}

}







