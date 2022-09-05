

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

//pick an instruction set for intrinsics by selecting a name space

//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;  //sse2   double
//using namespace DRC::VecD4D;	//avx2   double
//using namespace DRC::VecF8F;	// avx2  float
//using namespace DRC::VecD8D;  //avx512 double
//using namespace DRC::VecF16F; //avx512   float
const double billion = 1000000000.0;


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


auto getRandomShuffledVectorxxx(int SZ, int instance_number = 0)
{
	using FloatType = double; // typename InstructionTraits<VecXX::INS>::FloatType;


	static std::map<int, std::vector<double> > vectors;


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
		std::vector<double>  v(SZ, double(6.66));
		for (int i = 0; i < SZ; i++) { v[i] += double(SZ / 2) + i; }
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(begin(v), end(v), g);
		vectors[key] = v;
		return v;
	}
}


auto runFunctionOverDifferentSize = [](int testRepeats, int vec_start_size, int vec_stepSZ, int vec_maxSZ, const auto& func, long testLoopSZ)
{

	RunResults results;

	for (int j = 0; j < testRepeats; ++j)
	{
		int VEC_SZ = vec_start_size;
		for (; VEC_SZ < vec_maxSZ; VEC_SZ += vec_stepSZ)
		{
			func(VEC_SZ, testLoopSZ);
		}
	}
	return results;
};




void doAVXMax512Dance()
{
	/*
	const long TEST_LOOP_SZ = 1000;
	const int repeatRuns = 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;
	*/

	const int sleepTime = 10000.;/// 10 seconds


	const int maxVectorSize = 4400;
	const int minVectorSize = 3800;
	const long TEST_LOOP_SZ = 100000;
	const int vectorStepSize = 8;
	const int repeatRuns = 10;

	auto zero = 0.0;// InstructionTraits<VecXX::INS>::nullValue;

	getRandomShuffledVectorxxx(-1); // reset  random input vectors

								 /*
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

*/

	auto DR3_accumulate = [&](int SZ, long TEST_LOOP_SZ)
	{
		using namespace DRC::VecD8D;

		double time = 0.;
		volatile  double res = 0.;

		// generic lambda for max either calling a max instruction or doing a selection with iff
	//	 auto mxDbl = [](auto lhs, auto rhs) { return max(lhs, rhs); };
		auto mxDbl = [](auto lhs, auto rhs) { return iff(lhs > rhs, lhs, rhs); }; //using iff fastest 

		auto v1 = getRandomShuffledVectorxxx(SZ, 0); // std stl vector double or float 
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
		//return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};


	auto DR3_accumulate2 = [&](int SZ, long TEST_LOOP_SZ)
	{
		using namespace DRC::VecD4D;

		double time = 0.;
		volatile  double res = 0.;

		// generic lambda for max either calling a max instruction or doing a selection with iff
	//	 auto mxDbl = [](auto lhs, auto rhs) { return max(lhs, rhs); };
		auto mxDbl = [](auto lhs, auto rhs) { return iff(lhs > rhs, lhs, rhs); }; //using iff fastest 

		auto v1 = getRandomShuffledVectorxxx(SZ, 0); // std stl vector double or float 
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
		//return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

	};



	for (;;)
	{

		double time = 0.0;
		{
			
			for (int K = 0; K < 4; K++)
			{
				time = 0.;
				std::cout << "AVX 512 " << K + 1 << "of 4" << std::endl;
				auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_accumulate, TEST_LOOP_SZ);
				std::cout << "AVX 512 " << K + 1 << "of 4" << time << "seconds   now sleep" << std::endl;

				using namespace std::chrono_literals;
				std::this_thread::sleep_for(15000ms);

			}
		}
		using namespace std::chrono_literals;
		std::this_thread::sleep_for(15000ms);

		{	
			
			for (int K = 0; K < 4; K++)
			{
				time = 0.;
				std::cout << "AVX 2 " << K + 1 << "of 4" << std::endl;
				auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_accumulate2, TEST_LOOP_SZ);
				std::cout << "AVX 2 " << K + 1 << "of 4" << time << "seconds   now sleep" << std::endl;
				using namespace std::chrono_literals;
				std::this_thread::sleep_for(15000ms);

			}
		}

	}

}


