// VariadicReduction.cpp : This file contains the 'main' function. Program execution begins and ends there.
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
#include <execution>


#include <iomanip>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_map>



#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/accumulate_transform.h"
#include "../Vectorisation/VecX/error_utils.h"


#include "../Vectorisation/VecX/zip_utils.h"

//#include "norm.h"


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
//using namespace DRC::VecLDb;  //long double 
//using namespace DRC::VecDb; //double
//using namespace DRC::VecD2D;  //sse2   double
//using namespace DRC::VecD4D;	//avx2   double
//using namespace DRC::VecF8F;	// avx2  float

using namespace DRC::VecD8D;  //avx512 double
//using namespace DRC::VecF16F; //avx512   float


using FLOAT = typename InstructionTraits<VecXX::INS>::FloatType;



const double billion = 1000000000.0;





// Helper function to call a single lambda with the current state and value
template <std::size_t Index, typename Tuple, typename State, typename Value>
auto applyLambdaX(const Tuple& lambdas, State state, Value value) {
	return std::get<Index>(lambdas)(state, value);
}

// Recursive case: apply the lambdas to each value in the sequence
template <typename Tuple, typename InIt, typename... States,
	std::size_t... Indices>
auto applyLambdasImplX(const Tuple& lambdas, InIt first, InIt last,
	std::tuple<States...> states,
	std::index_sequence<Indices...>) {
	for (; first != last; ++first) {
		states = std::make_tuple(applyLambdaX<Indices>(
			lambdas, std::get<Indices>(states), *first)...);
	}
	return states;
}

// Main function to apply multiple lambdas to each value in the sequence
template <typename InIt, typename... Lambdas, typename... States>
auto reduceVariadic(InIt first, InIt last, const std::tuple<Lambdas...>& lambdas,
	std::tuple<States...> initStates) {
	if (first == last) {
		throw std::invalid_argument("Input range cannot be empty");
	}
	return applyLambdasImplX(lambdas, first, last, initStates,
		std::index_sequence_for<Lambdas...>{});
}

// timing size constants
constexpr unsigned TIMING_SZ = 10;// 1000;  // number of loops used for timing
constexpr const long SZ = 1000000l;   // size of data

// container utils
template <typename T>
void update(T& map) {
	for (long l = 0; l < SZ; ++l) {
		map.emplace(l, static_cast<double>(l));
	}
}

void updateVec(std::vector<double>& vec) {
	for (long l = 0; l < SZ; ++l) {
		vec.emplace_back(static_cast<double>(l));
	}
}

// Example 1: variadic reduction over an unordered map

void runMap() {
	std::unordered_map<long, double> testMap;
	// set up map
	update(testMap);

	// define reduction lambdas
	auto sum = [](auto lhs, auto rhs) {
		lhs.second += rhs.second;
		return lhs;
		};

	auto sum1 = [](auto lhs, auto rhs) {
		lhs.second += rhs.second * rhs.second;
		return lhs;
		};

	auto sum2 = [](auto lhs, auto rhs) {
		lhs.second += rhs.second * rhs.second * rhs.second;
		return lhs;
		};

	auto sum3 = [](auto lhs, auto rhs) {
		lhs.second += rhs.second * rhs.second * rhs.second * rhs.second;
		return lhs;
		};

	auto sum4 = [](auto lhs, auto rhs) {
		lhs.second +=
			rhs.second * rhs.second * rhs.second * rhs.second * rhs.second;
		return lhs;
		};

	std::pair<long, double> nil = { 0L, 0.0 };

	auto result = std::reduce(begin(testMap), end(testMap), nil, sum);
	auto result1 = std::reduce(begin(testMap), end(testMap), nil, sum1);
	auto result2 = std::reduce(begin(testMap), end(testMap), nil, sum2);
	auto result3 = std::reduce(begin(testMap), end(testMap), nil, sum3);
	auto result4 = std::reduce(begin(testMap), end(testMap), nil, sum4);

	volatile auto db = result1.second;

	auto start = std::chrono::system_clock::now();  // start timer

	// reduce over the data vector separately for each lambda
	for (int u = 0; u < TIMING_SZ; u++) {
		result = std::reduce(begin(testMap), end(testMap), nil, sum);
		result1 = std::reduce(begin(testMap), end(testMap), nil, sum1);
		result2 = std::reduce(begin(testMap), end(testMap), nil, sum2);
		result3 = std::reduce(begin(testMap), end(testMap), nil, sum3);
		result4 = std::reduce(begin(testMap), end(testMap), nil, sum4);
	}

	auto last = std::chrono::system_clock::now();  // stop timer
	auto elapsed_seconds =
		std::chrono::duration_cast<std::chrono::microseconds>(last - start)
		.count();

	std::cout << std::setprecision(12);
	std::cout << "\n running reducing over map with sum and sum1, sum2, sum3 "
		"sum 4 time = "
		<< elapsed_seconds / 1000000.0 << "\n";
	std::cout << "\n";

	std::cout << result.second << "   ,  " << result1.second << "   ,  "
		<< result2.second << "   ,  " << result3.second << "   ,  "
		<< result4.second << "\n ";

	// VARIADIC reduction
	// create lambda tuples
	auto lambdas = std::make_tuple(sum, sum1, sum2, sum3, sum4);
	// create initial value tuples
	auto init_values = std::make_tuple(nil, nil, nil, nil, nil);

	// call variadic reduce
	auto tuple_result =
		reduceVariadic(begin(testMap), end(testMap), lambdas, init_values);

	auto start2 = std::chrono::system_clock::now();

	// run the timing loop
	for (int u = 0; u < TIMING_SZ; u++) {
		tuple_result =
			reduceVariadic(begin(testMap), end(testMap), lambdas, init_values);
	}
	auto last2 = std::chrono::system_clock::now();

	auto elapsed_seconds2 =
		std::chrono::duration_cast<std::chrono::microseconds>(last2 - start2)
		.count();

	std::cout << std::setprecision(12);
	std::cout << "\n running variadic reduce over map with sum and sum1, sum2, "
		"sum3 sum4 time = "
		<< elapsed_seconds2 / 1000000.0 << "\n";
	std::cout << "\n";

	std::cout << std::get<0>(tuple_result).second << ",   ";
	std::cout << std::get<1>(tuple_result).second << ",   ";
	std::cout << std::get<2>(tuple_result).second << ",   ";
	std::cout << std::get<3>(tuple_result).second << ",   ";
	std::cout << std::get<4>(tuple_result).second << ",   ";

	std::cout << "\n variadic reduce is faster by "
		<< (double)elapsed_seconds / elapsed_seconds2 << "\n";
}

// Example 2: variadic reduction over a vector
void runVector() {
	std::vector<double> vec;
	vec.reserve(SZ);

	for (long ll = 0l; ll < SZ; ++ll) {
		double val = static_cast<double>(ll);
		vec.push_back(val);
	}

	auto sum = [](auto lhs, auto rhs) { return lhs + rhs; };

	auto sum1 = [](double lhs, double rhs) { return lhs + (rhs * rhs); };

	auto sum2 = [](auto lhs, auto rhs) -> decltype(lhs) {
		return lhs + rhs * rhs * rhs;
		};

	auto sum3 = [](auto lhs, auto rhs) -> decltype(lhs) {
		return lhs + rhs * rhs * rhs * rhs;
		};

	auto sum4 = [](auto lhs, auto rhs) -> decltype(lhs) {
		return lhs + rhs * rhs * rhs * rhs * rhs;
		};

	const double nil = { 0.0 };

	auto start = std::chrono::system_clock::now();

	double result = .0;
	double result1 = 0.;
	double result2 = 0.0;
	double result3 = 0.0;
	double result4 = 0.0;

	for (int u = 0; u < TIMING_SZ; u++) {
		result = std::accumulate(begin(vec), end(vec), nil, sum);
		result1 = std::accumulate(begin(vec), end(vec), nil, sum1);
		result2 = std::accumulate(begin(vec), end(vec), nil, sum2);
		result3 = std::accumulate(begin(vec), end(vec), nil, sum3);
		result4 = std::accumulate(begin(vec), end(vec), nil, sum4);
	}

	auto last = std::chrono::system_clock::now();
	auto elapsed_seconds =
		std::chrono::duration_cast<std::chrono::microseconds>(last - start)
		.count();

	std::cout << std::setprecision(12);
	std::cout << "\n running reducing over vector with sum and sum1, sum2, "
		"sum3 sum 4 time = "
		<< elapsed_seconds / 1000000.0 << "\n";
	std::cout << "\n";

	std::cout << result << "   ,  " << result1 << "   ,  " << result2
		<< "   ,  " << result3 << "   ,  " << result4 << "\n ";

	auto lambdas = std::make_tuple(sum, sum1, sum2, sum3, sum4);
	const auto init_values = std::make_tuple(nil, nil, nil, nil, nil);

	auto tuple_result = reduceVariadic(begin(vec), end(vec), lambdas, init_values);

	auto start2 = std::chrono::system_clock::now();

	for (int u = 0; u < TIMING_SZ; u++) {
		tuple_result = reduceVariadic(begin(vec), end(vec), lambdas, init_values);
	}

	auto last2 = std::chrono::system_clock::now();
	auto elapsed_seconds2 =
		std::chrono::duration_cast<std::chrono::microseconds>(last2 - start2)
		.count();

	std::cout << std::setprecision(12);

	std::cout << "\n running variadic reduce over vector with sum and sum1, "
		"sum2, sum3 sum4 time = "
		<< elapsed_seconds2 / 1000000.0 << "\n";
	std::cout << "\n";

	std::cout << std::get<0>(tuple_result) << ",   ";
	std::cout << std::get<1>(tuple_result) << ",   ";
	std::cout << std::get<2>(tuple_result) << ",   ";
	std::cout << std::get<3>(tuple_result) << ",   ";
	std::cout << std::get<4>(tuple_result) << ",   ";

	// std::cout << db << "\n";

	std::cout << "variadic reduce is faster by "
		<< (double)elapsed_seconds / elapsed_seconds2 << "\n";
}





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


long double getErr(const  std::vector<long double>& t)
{
	ignore(t);
	return 2e-11;
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

bool valuesAreEqual(double x, double y, double tol = 1e-14)
{

	auto err1 = fabs((x - y) / (fabs(x) + fabs(y)));

	return (err1 > tol) ? false : true;
}


auto getRandomShuffledVector(int SZ, int instance_number = 0)
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







auto numOps = [](int TEST_LOOP_SZ, int SZ) { return  static_cast<int>(double(TEST_LOOP_SZ) * double(SZ)); };


double getnull(double)
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



//example doing  multi reduction operation after load
//should be faster 
void doMultiStats(int runType =7)
{

	const long TEST_LOOP_SZ = 10;// 100;// 100;//400;// 1000;
	const int repeatRuns = 20;// 20;
	const int vectorStepSize = 200;
	const int maxVectorSize = 20000;
	const int minVectorSize = 400;

	getRandomShuffledVector(-1); // reset  random input vectors


		auto accumulate_run = [&](int VEC_SZ, long TEST_LOOP_SZ)
		{
			double time = 0.;
			double res = 0.;
			double res_min = 0.;
			volatile  double total = 0.;  //sum to prevent optimizing out
		
			auto v1 = getRandomShuffledVector(VEC_SZ, 0);

			auto sum_vals = [](double accumulate, double val) {return accumulate + val; };
			auto sumSqr = [](double accumulate, double val) {return accumulate + val * val; };
			auto sumCube = [](double accumulate, double val) {return accumulate + val * val * val; };
			auto sumQuart = [](double accumulate, double val) {return accumulate + val * val * val * val; };

			//warm up
			for (long l = 0; l < 100; l++)
			{
				auto sum_val = std::reduce( v1.begin(), v1.end(), 0.0, sum_vals);
				auto sum_val1 = std::reduce( v1.begin(), v1.end(), 0.0, sumSqr);
				auto sum_val2 = std::reduce( v1.begin(), v1.end(), 0.0, sumCube);
				auto sum_val3 = std::reduce( v1.begin(), v1.end(), 0.0, sumQuart);
				
				total = sum_val + sum_val1 + sum_val2 + sum_val3;
			}


			{
				TimerGuard timer(time);
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{

					auto sum_val = std::reduce( v1.begin(), v1.end(), 0.0, sum_vals);
					auto sum_val1 = std::reduce( v1.begin(), v1.end(), 0.0, sumSqr);
					auto sum_val2 = std::reduce( v1.begin(), v1.end(), 0.0, sumCube);
					auto sum_val3 = std::reduce( v1.begin(), v1.end(), 0.0, sumQuart);


					total = sum_val + sum_val1 + sum_val2 + sum_val3; 
					//compute combined value over all results  to check tests and stop short circuit
				}
			}

			ignore(total);
			ignore(res_min);

			return  std::make_pair(res, numOps(TEST_LOOP_SZ, VEC_SZ) / time);
		};



		auto DR3_accumulate = [&](int SZ, long TEST_LOOP_SZ)
		{
			double time = 0.;
			volatile  double res = 0.;
		

			auto sum = [](auto accumulate, auto val) {return accumulate + val; };
			auto sumSqr = [](auto accumulate, auto val) {return accumulate + val * val; };
			auto sumCube = [](auto accumulate, auto val) {return accumulate + val * val * val; };
			auto sumQuart = [](auto accumulate, auto val) {return accumulate + val * val * val * val; };

			auto v1 = getRandomShuffledVector(SZ, 0); // std stl vector double or float 
			VecXX vec(v1);


			for (long l = 0; l < 100; l++)
			{
				double sum_mnn = reduce(vec, sum);
				double sumSqr_mnn = reduce(vec, sumSqr);
				double sum_Cube = reduce(vec, sumCube);
				double sum_Qrt = reduce(vec, sumQuart);

				res = 0.5 * (sum_mnn + sumSqr_mnn + sum_Cube + sum_Qrt);
			}


			{
				TimerGuard timer(time);
				for (long l = 0; l < TEST_LOOP_SZ; l++)
				{
					double sum_mnn = reduce(vec, sum);
					double sumSqr_mnn = reduce(vec, sumSqr);
					double sum_Cube = reduce(vec, sumCube);
					double sum_Qrt = reduce(vec, sumQuart);

					res = 0.5 * (sum_mnn + sumSqr_mnn + sum_Cube + sum_Qrt);
				}
			}

			return std::make_pair(res, numOps(TEST_LOOP_SZ, SZ) / time);

		};

		auto DR3_accumulate_multi = [&](int SZ, long TEST_LOOP_SZ)
			{
				double time = 0.;
				volatile  double res = 0.;

				auto sum = [](auto accumulate, auto val) {return accumulate + val; };
				auto sumSqr = [](auto accumulate, auto val) {return accumulate + val * val; };
				auto sumCube = [](auto accumulate, auto val) {return accumulate + val * val * val; };
				auto sumQuart = [](auto accumulate, auto val) {return accumulate + val * val * val * val; };

				auto v1 = getRandomShuffledVector(SZ, 0); // std stl vector double or float 
				VecXX vec(v1);

				auto ress = reduceM(vec, sum, sumSqr, sumCube, sumQuart);
				//warm up
				for (long l = 0; l < 100; l++)
				{
					ress = reduceM(vec, sum, sumSqr, sumCube, sumQuart);
				}


				{
					TimerGuard timer(time);
					for (long l = 0; l < TEST_LOOP_SZ; l++)
					{
						ress = reduceM(vec, sum, sumSqr, sumCube, sumQuart);

						double sum_mnn = std::get<0>(ress);
						double sumSqr_mnn = std::get<1>(ress);
						double sum_Cube = std::get<2>(ress);
						double sum_Qrt = std::get<3>(ress);

						res = 0.5 * (sum_mnn + sumSqr_mnn + sum_Cube + sum_Qrt);
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

			std::cout << std::setprecision(8);
			//print out results
			for (const auto& perf_stl : stats_stl)
			{
				auto  valDr3 = dr3_raw_results.m_calc_results[perf_stl.first];
				auto  valStl = run_res_stl.m_calc_results[perf_stl.first];
				auto  valDr3Mult = dr3_raw_results_mult.m_calc_results[perf_stl.first];
				auto strMatch = valuesAreEqual(valDr3, valStl, valDr3Mult) ? "calcs match" : "cal difference";
				std::cout << "  std::reduce, size " << perf_stl.first << ", " << perf_stl.second.first << ", + - ," << perf_stl.second.second
					<< "\t \t DR3 reduce, size " << perf_stl.first << ", " << stats_DR3_perf[perf_stl.first].first << ",  + - ," << stats_DR3_perf[perf_stl.first].second
					<< "\t \t DR3 reduce_mult, size " << perf_stl.first << ", " << stats_DR3_perf_mult[perf_stl.first].first << ",  + - ," << stats_DR3_perf_mult[perf_stl.first].second
					<< ", numerical check : " << strMatch << "\n";
			}
	
}




////variadic reduce driver//////////////////

int main()
{
	// variadic reduction over map
	runMap();

	//variadic reduction over std vector
	runVector();

	//multiple reductions over SIMD vector VecXX
	doMultiStats(7);
	return 0;
}


