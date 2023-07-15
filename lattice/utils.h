#pragma once


#include "../Vectorisation/VecX/error_utils.h"
#include <vector>
#include <chrono>
#include <map>
#include <iostream>

#include <algorithm>
#include <random>
#include <numeric>
#include <iterator>
#include <iomanip>  
#include <cstring>


//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;  //sse2   double
//using namespace DRC::VecD4D;	//avx2   double
//using namespace DRC::VecF8F;	// avx2  float
using namespace DRC::VecD8D;  //avx512 double
//using namespace DRC::VecF16F; //avx512   float

using FLOAT = typename InstructionTraits<VecXX::INS>::FloatType;


//using FLOAT = typename InstructionTraits<VecXX::INS>::FloatType;

constexpr double billion = 1000000000.0;


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


static double getErr(const  std::vector<double>& t)
{
	ignore(t);
	return 2e-11;
}


static float getErr(const  std::vector<float>& t)
{
	ignore(t);
	return 1.0;
}


template<typename T>
bool vectorsEqualD(const std::vector<T>& C1, const std::vector<T>& C2, const std::vector<T>& C3, const std::vector<T>& input, T ERR = 1e-13)
{


	ERR = getErr(C1);

	bool  testOK = true;

	if (C1.size() != C2.size())
	{
		std::cout << "wrong size C1,C2" << C1.size() << ", " << C2.size() << std::endl;
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
		auto err2 = fabs((C1[i] - C3[i]) / (fabs(C1[i]) + fabs(C3[i])));

		if ((err1 > ERR) || (err2 > ERR))
		{
			testOK = false;
			std::cout << "\n err diff@ " << i << " err1 =" << err1 << ", err2 = " << err2 << "\n";
			std::cout << "\n val @ " << i << " C1[i] =" << C1[i] << ", C2[i] = " << C2[i] << ", C3[i] = " << C3[i] << "input val=" << input[i] << "\n";
			std::cout << std::endl;
			break;
		}
	}

	return testOK;

}

static bool valuesAreEqual(double x, double y, double tol = 1e-14)
{

	auto err1 = fabs((x - y) / (fabs(x) + fabs(y)));

	return (err1 > tol) ? false : true;
}


static auto getRandomShuffledVector(int SZ, int instance_number = 0)
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
		for (int i = 0; i < SZ; i++) { v[i] += /*FloatType(SZ / 2)*/ +i; }
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(begin(v), end(v), g);
		vectors[key] = v;
		return v;
	}
}







static auto numOps = [](int TEST_LOOP_SZ, int SZ) { return  static_cast<int>(double(TEST_LOOP_SZ) * double(SZ)); };


static double getnull(double)
{
	return 0.0;
}

using Calc_Values = std::map<int, FLOAT>;
using Calc_Values_V = std::map<int, std::vector<FLOAT> >;
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


static auto runFunctionOverDifferentSize = [](int testRepeats, int vec_start_size, int vec_stepSZ, int vec_maxSZ, const auto& func, long testLoopSZ)
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


static auto runFunctionOverDifferentSizeVec = [](int testRepeats, int vec_start_size, int vec_stepSZ, int vec_maxSZ, const auto& func, long testLoopSZ)
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



static auto performanceStats = [](const Mapped_Performance_Results& raw_results)
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


