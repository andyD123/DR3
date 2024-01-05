

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


const double billion = 1000000000.0;


using Calc_Values = std::map<int, double>;
using  Mapped_Performance_Results = std::map<int, std::vector<double> >; // array size  v vector<throughput for runs>
using Mapped_Stats = std::map<int, std::pair<double, double> >; // size -.pair ( throughput ,  std dev of through put)

struct RunResults
{
	Mapped_Performance_Results m_raw_results;
	Calc_Values  m_calc_results;
	double time;
};

class TimerGuard
{
	double& m_runTime;
	std::chrono::high_resolution_clock::time_point  m_startTme;

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

	const int maxVectorSize = 4400;
	const int minVectorSize = 800;
	const long TEST_LOOP_SZ = 10000;
	const int vectorStepSize = 8;
	const int repeatRuns = 13;

	getRandomShuffledVectorxxx(-1); // reset  random input vectors


	//avx512 lambda
	auto DR3_avx512 = [&](int SZ, long TEST_LOOP_SZ)
	{
		using namespace DRC::VecD8D;

		double time = 0.;
		volatile  double res = 0.;

		auto mxDbl = [](auto lhs, auto rhs) { return iff(lhs > rhs, lhs, rhs); };

		auto v1 = getRandomShuffledVectorxxx(SZ, 0);
		VecXX vec(v1);


		for (long l = 0; l < TEST_LOOP_SZ; l++)
		{
			res = reduce(vec, mxDbl);
		}

		return std::make_pair(res, time);
	};


	auto DR3_avx2 = [&](int SZ, long TEST_LOOP_SZ)
	{
		using namespace DRC::VecD4D;

		double time = 0.;
		volatile  double res = 0.;

		auto mxDbl = [](auto lhs, auto rhs) { return iff(lhs > rhs, lhs, rhs); };

		auto v1 = getRandomShuffledVectorxxx(SZ, 0);
		VecXX vec(v1);


		for (long l = 0; l < TEST_LOOP_SZ; l++)
		{
			res = reduce(vec, mxDbl);
		}


		return std::make_pair(res, time);

	};




	auto DR3_sse2 = [&](int SZ, long TEST_LOOP_SZ)
	{
		using namespace DRC::VecD2D;

		double time = 0.;
		volatile  double res = 0.;

		auto mxDbl = [](auto lhs, auto rhs) { return iff(lhs > rhs, lhs, rhs); };

		auto v1 = getRandomShuffledVectorxxx(SZ, 0);
		VecXX vec(v1);


		for (long l = 0; l < TEST_LOOP_SZ; l++)
		{
			res = reduce(vec, mxDbl);
		}

		return std::make_pair(res, time);

	};



	auto DR3_stl = [&](int SZ, long TEST_LOOP_SZ)
	{
		using namespace DRC::VecD2D;

		double time = 0.;
		volatile  double res = 0.;

		auto v1 = getRandomShuffledVectorxxx(SZ, 0);


		for (long l = 0; l < TEST_LOOP_SZ; l++)
		{
			res = *std::max_element(begin(v1), end(v1));
		}

		return std::make_pair(res, time);
	};




	using namespace std::chrono_literals;

	for (;;)
	{

		double time = 0.0;
		constexpr int NUM_BURSTS = 3;
		constexpr auto SLEEP_TIME = 20000ms;





		//STL
		for (int K = 0; K < NUM_BURSTS; K++)
		{
			time = 0.;
			std::cout << "1/3rd the work using STL max " << K + 1 << "of " << NUM_BURSTS << std::endl;
			{	TimerGuard timer(time);
			auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_stl, TEST_LOOP_SZ / 3);
			}
			std::cout << "STL " << K + 1 << " of  " << NUM_BURSTS << "    " << time << " seconds   now sleep" << std::endl;
			std::this_thread::sleep_for(SLEEP_TIME);
		}

		std::this_thread::sleep_for(SLEEP_TIME);



		//SSE2
		for (int K = 0; K < NUM_BURSTS; K++)
		{
			time = 0.;
			std::cout << "SSE2 " << K + 1 << " of  " << NUM_BURSTS << std::endl;
			{	TimerGuard timer(time);
			auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_sse2, TEST_LOOP_SZ);
			}
			std::cout << "SSE2 " << K + 1 << " of  " << NUM_BURSTS << "    " << time << " seconds   now sleep" << std::endl;
			std::this_thread::sleep_for(SLEEP_TIME);
		}

		std::this_thread::sleep_for(SLEEP_TIME);


		//AVX2 
		for (int K = 0; K < NUM_BURSTS; K++)
		{
			time = 0.;
			std::cout << "AVX2 " << K + 1 << " of  " << NUM_BURSTS << std::endl;
			{	TimerGuard timer(time);
			auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_avx2, TEST_LOOP_SZ);
			}
			std::cout << "AVX2 " << K + 1 << " of  " << NUM_BURSTS << "    " << time << " seconds   now sleep" << std::endl;
			std::this_thread::sleep_for(SLEEP_TIME);

		}

		std::this_thread::sleep_for(SLEEP_TIME);

		//AVX512
		for (int K = 0; K < NUM_BURSTS; K++)
		{
			time = 0.;
			std::cout << "AVX512 " << K + 1 << " of  " << NUM_BURSTS << std::endl;
			{
				TimerGuard timer(time);
				auto dr3_raw_results = runFunctionOverDifferentSize(repeatRuns, minVectorSize, vectorStepSize, maxVectorSize, DR3_avx512, TEST_LOOP_SZ);
			}
			std::cout << "AVX512 " << K + 1 << " of  " << NUM_BURSTS << "    " << time << " seconds   now sleep" << std::endl;
			std::this_thread::sleep_for(SLEEP_TIME);
		}

		std::this_thread::sleep_for(SLEEP_TIME);





	}

}




