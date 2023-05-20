// lattice.cpp : This file contains the 'main' function. Program execution begins and ends there.
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

//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;  //sse2   double
//using namespace DRC::VecD4D;	//avx2   double
//using namespace DRC::VecF8F;	// avx2  float
using namespace DRC::VecD8D;  //avx512 double
//using namespace DRC::VecF16F; //avx512   float


using FLOAT = typename InstructionTraits<VecXX::INS>::FloatType;


void testSampler()
{

	FLOAT sampleData[] = { 0., 1.0,2.,3.,4.,5.,6.,7.,8.,9.,10. };

	//TriSampler<DRC::VecD4D::VecXX::INS>  zzz;
	TrinomialSampler<VecXX::INS>  zzz;

	zzz.load(&sampleData[1]);

	auto x = zzz.get<-1>();
	auto y = zzz.get<0>();
	auto z = zzz.get<1>();

	ignore(x);
	ignore(y);
	ignore(z);
}




using FLOAT = InstructionTraits<VecXX::INS>::FloatType;

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



/*
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
*/

void doScan()
{

	auto v1 = std::vector<VecXX::SCALA_TYPE>(4000, 0.0);
	FLOAT last = 0.0;// -1.0;
	for (auto& x : v1)
	{
		x = last + 1.0;
		last = x;
	}
	// getRandomShuffledVector(200);
	VecXX a(v1);
	std::vector<FLOAT> dbg = v1;

	//auto sum = [](auto X, auto Y) { return X + Y; };

	auto sum = [](auto X, auto Y) { return  X + Y; };


	double time = 0;
	{
		TimerGuard timer(time);

		//a += 1.0;

		VecXX result;
		for (long l = 0; l < 1000000; ++l)
		{
			result = ApplyScan(a, sum);
			//dbg = a;

		}

		dbg = result;
	}
	std::cout << "run time = " << time << "\n";

	///*
	auto dbg_cpy = dbg;
	time = 0;
	{
		TimerGuard timer(time);

		//a += 1.0;

		for (long l = 0; l < 1000000; ++l)
		{
			//VecXX result = ApplyScan8(a, sum);
			std::inclusive_scan(begin(v1), end(v1), dbg.begin());
			//	dbg = a;
			//	dbg = result;


		}
	}
	std::cout << "run time std::scan = " << time << "\n";

	//	/*
	for (size_t i = 0; i < dbg_cpy.size(); ++i)
	{
		std::cout << i << "  ,  " << dbg[i] << "   , " << dbg_cpy[i] << std::endl;
	}

	//	*/



}



void doZipping()
{

	auto v1 = getRandomShuffledVector(200);
	VecXX a(v1);
	VecXX b = 2. * a;
	VecXX c = 3. * a;
	VecXX d = 4. * a;
	VecXX e = 5. * a;

	auto zipped = make_Zipped<VecXX::INS>(a, b, c);
	auto it = make_Zipped_itr<VecXX::INS>(a, b, c);


	//auto zipped =
	const VecXX& aa = a;
	const VecXX& bb = b;
	const VecXX& cc = c;
	const VecXX& dd = d;
	const VecXX& ee = e;

	auto zp5 = make_Zipped<VecXX::INS>(aa, bb, cc, dd, ee);
	auto zp5_itr = make_Zipped_itr<VecXX::INS>(aa, bb, cc, dd, ee);

	ignore(zp5);
	ignore(zp5_itr);


	//it.inc(1);

	auto addingLambda = [](auto& zpped)
	{
		const auto& X = std::get<0>(zpped.m_registers);
		const auto& Y = std::get<1>(zpped.m_registers);
		const auto& Z = std::get<2>(zpped.m_registers);

		return X + Y + Z;
	};


	VecXX::INS X = 0.;
	for (int i = 0; i < 200 / 8; ++i)
	{
		zipped.load(it);
		auto sum = addingLambda(zipped);
		X += sum;
		it.inc(1);
	}

	std::cout << X[0];

	zipped.load(it);


	auto addingLambda5 = [](auto& zpped)
	{
		const auto& X = std::get<0>(zpped.m_registers);
		const auto& Y = std::get<1>(zpped.m_registers);
		const auto& Z = std::get<2>(zpped.m_registers);
		const auto& L = std::get<3>(zpped.m_registers);
		const auto& M = std::get<4>(zpped.m_registers);

		return X + Y + Z + L + M;
	};



	//Vec<INS_VEC> res =   transform(OP & oper, const Zipped_ITR< INS_VEC, N>&zip)

	auto res = transform(addingLambda5, zp5_itr);

	std::vector<FLOAT> vdb = res;

	enum class Access { down = -1, mid = 0, up = 1 };
	//	using Mysample = Named_Zip_iter < VecXX::INS, 3, Access > ;

		//NAMED_INDEX

		//Mysample instance;


	//	Named_Zipped_Reg < VecXX::INS, 3, Access > myReg;
	Zipped_Reg < VecXX::INS, 3 > myReg;
	//	using reg = Zipped_Reg < VecXX::INS, 3 >;

		//auto zz = Access::down;
	auto yy = myReg.get<int(Access::up)>();
	ignore(yy);



	//////////////////////////////


	VecXX out1 = a;
	VecXX out2 = b;

	VecXX& out1_r = out1;
	VecXX& out2_r = out2;


	auto zp_out = make_Zipped_itr_ref<VecXX::INS>(out1_r, out2_r);

	auto addingLambdaIO = [&](const auto& zpped, auto& out_zipp)
	{
		const auto& X = std::get<0>(zpped.m_registers);
		const auto& Y = std::get<1>(zpped.m_registers);
		const auto& Z = std::get<2>(zpped.m_registers);

		std::get<0>(out_zipp.m_registers) = X + Y + Z;
		std::get<1>(out_zipp.m_registers) = X + Y;
	};


	transform(addingLambdaIO, zp5_itr, zp_out);


	std::vector<FLOAT> vdb1 = out1;
	std::vector<FLOAT> vdb2 = out2;

}



double europeanBinomialPricer(double S, double K, double sig, double r, double T, int N)
{



	//bool bisOdd = (N & 1);

	VecXX terminalAssetPrices(1.0, N + 1);

	double Dt = T / N;
	double  u = std::exp(sig * std::sqrt(Dt));
	double d = 1. / u;


	VecXX::INS pu = (exp(r * Dt) - d) / (u - d);
	VecXX::INS oneMinusP = (1.0 - pu);
	VecXX::INS disc = exp(-r * Dt);

	BinomialSampler<VecXX::INS> sampler;

	auto binomialRollBack = [=](BinomialSampler<VecXX::INS>& sampler)
	{
		auto X1 = sampler.get<1>();
		auto X0 = sampler.get<0>();
		return disc * (X1 * pu + X0 * oneMinusP);
	};


	auto payOffFunc = [=](auto X) { return select(X > K, X - K, 0.0); };
	//	auto payOffFunc = [=](auto X) { return select(X < K, K - X, 0.0); };

		//set up underlying asset prices at maturity
	double last = S * std::pow(d, N + 2);
	for (auto& el : terminalAssetPrices)
	{
		last *= (u * u);
		el = last;
	}


	auto odd_slice = transform(payOffFunc, terminalAssetPrices);
	auto even_slice = odd_slice;

	int j = N + 1;
	for (int i = 0; i < N / 2; ++i)
	{
		transform(odd_slice, even_slice, binomialRollBack, sampler, 0, j);
		transform(even_slice, odd_slice, binomialRollBack, sampler, 0, j - 1);
		j -= 2;

	}
	return odd_slice[0];;
}



void doBinomialPricer()
{


	double res = 0.0;
	double time = 0;
	{
		TimerGuard timer(time);



		double S = 100.;
		double K = 100.;
		double vol = 0.2;
		double rate = 0.06;
		double T = 1;
		int N = 10000;

		res = europeanBinomialPricer(S, K, vol, rate, T, N);

	}
	std::cout << std::setprecision(12) << "price " << res << "\n";
	std::cout << " takes secs" << time << "\n";
}


double europeanTrinomialPricer(double S, double K, double sig, double r, double T, int N)
{

	double y = 0.0;// 0.03; //div yield


	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	//double Dx = sig * std::sqrt(3.0* Dt);
	double Dx = sig * std::sqrt(2.0 * Dt);

	double v = r - y - 0.5 * sig * sig;


	//double  u = Dx;// std::exp(sig * std::sqrt(Dt));
	//double d = 1. / u;


	VecXX::INS pu = 0.5 * ((Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx) + (v * Dt) / Dx);

	VecXX::INS pd = 0.5 * ((Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx) - (v * Dt) / Dx);

	VecXX::INS pm = 1. - (Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx);


	VecXX::INS disc = exp(-r * Dt);

	TrinomialSampler<VecXX::INS> sampler;

	auto trinomialRollBack = [=](TrinomialSampler<VecXX::INS>& sampler)
	{
		auto X1 = sampler.get<1>();
		auto X0 = sampler.get<0>();
		auto X_1 = sampler.get<-1>();
		return disc * (X1 * pu + X0 * pm + X_1 * pd);
	};


	auto payOffFunc = [=](auto X) { return select(X > K, X - K, 0.0); };

	//set up underlying asset prices at maturity
	double last = S * exp(-(N + 1) * Dx);
	double edx = exp(Dx);
	for (auto& el : terminalAssetPrices)
	{
		last *= edx;
		el = last;
	}

	auto odd_slice = transform(payOffFunc, terminalAssetPrices);
	auto even_slice = odd_slice;

	int j = 2 * N + 1 - 1;
	int i = 0;
	for (; i < N; i += 2)
	{
		transform(odd_slice, even_slice, trinomialRollBack, sampler, i, j);
		transform(even_slice, odd_slice, trinomialRollBack, sampler, i + 1, j - 1);
		j -= 2;
	}

	return odd_slice[N];
}




double americanTrinomialPricer(double S, double K, double sig, double r, double T, int N)
{

	double y = 0.0;// 0.03; //div yield


	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	double Dx = sig * std::sqrt(2.0 * Dt);
	double v = r - y - 0.5 * sig * sig;

	//double  u = Dx;
	//double d = 1. / u;


	VecXX::INS pu = 0.5 * ((Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx) + (v * Dt) / Dx);
	VecXX::INS pd = 0.5 * ((Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx) - (v * Dt) / Dx);
	VecXX::INS pm = 1. - (Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx);


	VecXX::INS disc = exp(-r * Dt);
	TrinomialSampler<VecXX::INS> sampler;

	auto trinomialRollBack = [=](TrinomialSampler<VecXX::INS>& sampler)
	{
		auto X1 = sampler.get<1>();
		auto X0 = sampler.get<0>();
		auto X_1 = sampler.get<-1>();
		return disc * (X1 * pu + X0 * pm + X_1 * pd);
	};

	//call
	auto payOffFunc = [=](auto X) { return select(X > K, X - K, 0.0); };

	//put
	//auto payOffFunc = [=](auto X) { return select(X < K, K -X , 0.0); };

	//set up underlying asset prices at maturity
	double last = S * exp(-(N + 1) * Dx);
	double edx = exp(Dx);
	for (auto& el : terminalAssetPrices)
	{
		last *= edx;
		el = last;
	}

	auto excerciseValue = transform(payOffFunc, terminalAssetPrices);
	auto odd_slice = excerciseValue;

	UnitarySampler<VecXX::INS> identity_sampler; //identity just 

	auto applyEarlyExcercise = [=](UnitarySampler<VecXX::INS>& sampler, auto excercisePrice)
	{
		auto optPrice = sampler.get<0>();
		return max(optPrice, excercisePrice);
	};


	auto even_slice = odd_slice;

	int j = 2 * N + 1 - 1;
	int i = 0;
	for (; i < N; i += 2)
	{
		transform(odd_slice, even_slice, trinomialRollBack, sampler, i, j);
		// transform to get early excercise for american bit , iderntity sampler just passes values straight through
		transform(even_slice, excerciseValue, even_slice, applyEarlyExcercise, identity_sampler, i, j);

		transform(even_slice, odd_slice, trinomialRollBack, sampler, i + 1, j - 1);
		transform(odd_slice, excerciseValue, odd_slice, applyEarlyExcercise, identity_sampler, i + 1, j - 1);

		j -= 2;
	}

	return odd_slice[N];
}




void doTrinomialPricer()
{

	double res = 0.0;
	double time = 0;
	{
		TimerGuard timer(time);

		double S = 100.;
		double K = 100.;
		double vol = 0.2;
		double rate = 0.06;
		double T = 1;
		int N = 1000;// 500;
		res = americanTrinomialPricer(S, K, vol, rate, T, N);


	}
	std::cout << std::setprecision(12) << "price " << res << "\n";
	std::cout << " takes secs" << time << "\n";
}




double americanFiniteDiffPricer(double S, double K, double sig, double r, double T, int N)
{

	//dividend yield
	double y = 0.0;// 0.03;// 0.03; //div yield
	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	double Dx = sig * std::sqrt(2.0 * Dt);
	double v = r - y - 0.5 * sig * sig;

	//double  u = Dx;
	//double d = 1. / u;


	VecXX::INS pu = 0.5 * Dt * ((sig * sig) / (Dx * Dx) + v / Dx);
	VecXX::INS pd = 0.5 * Dt * ((sig * sig) / (Dx * Dx) - v / Dx);
	VecXX::INS pm = 1. - Dt * (sig * sig) / (Dx * Dx) - r * Dt;


	//VecXX::INS disc = exp(-r * Dt);
	TrinomialSampler<VecXX::INS> sampler;

	auto trinomialRollBack = [=](TrinomialSampler<VecXX::INS>& sampler)
	{
		auto X1 = sampler.get<1>();
		auto X0 = sampler.get<0>();
		auto X_1 = sampler.get<-1>();
		return  (X1 * pu + X0 * pm + X_1 * pd);
	};


	//Pay off functions
	//call
	auto payOffFunc = [=](auto X) { return select(X > K, X - K, 0.0); };

	//put
	//auto payOffFunc = [=](auto X) { return select(X < K, K - X, 0.0); };

	//set up underlying asset prices at maturity
	double last = S * exp(-(N + 1) * Dx);
	double edx = exp(Dx);
	for (auto& el : terminalAssetPrices)
	{
		last *= edx;
		el = last;
	}

	auto excerciseValue = transform(payOffFunc, terminalAssetPrices);
	auto odd_slice = excerciseValue;

	UnitarySampler<VecXX::INS> identity_sampler; //identity  

	auto applyEarlyExcercise = [=](UnitarySampler<VecXX::INS>& sampler, auto excercisePrice)
	{
		auto optPrice = sampler.get<0>();
		return max(optPrice, excercisePrice);
	};

	auto even_slice = odd_slice * 0.0;

	int J = 2 * N;

	int k = 0;
	for (; k < N; k += 2)
	{

		transform(odd_slice, even_slice, trinomialRollBack, sampler, 0, J);

		//boundary condition
		even_slice[0] = even_slice[1] + terminalAssetPrices[1] - terminalAssetPrices[0];
		even_slice[J] = even_slice[J - 1];

		// transform to get early excercise for american bit , iderntity sampler just passes values straight through
		transform(even_slice, excerciseValue, even_slice, applyEarlyExcercise, identity_sampler, 0, J);


		transform(even_slice, odd_slice, trinomialRollBack, sampler, 0, J);
		//boundary condition
		odd_slice[0] = odd_slice[1] + terminalAssetPrices[1] - terminalAssetPrices[0];
		odd_slice[J] = odd_slice[J - 1];
		transform(odd_slice, excerciseValue, odd_slice, applyEarlyExcercise, identity_sampler, 0, J);

	}

	return odd_slice[N];
}




void doAmericanFiniteDiff()
{
	double time = 0;
	double res = 0.0;
	{
		TimerGuard timer(time);



		std::cout << std::setprecision(12) << " start \n";
		double S = 100.;
		double K = 100.;
		double vol = 0.2;
		double rate = 0.06;
		double T = 1;
		int N = 10000;
		res = americanFiniteDiffPricer(S, K, vol, rate, T, N);


	}
	std::cout << std::setprecision(12) << "price " << res << "\n";
	std::cout << " takes secs" << time << "\n";
}




double americanImplicitFiniteDiffPricerFast(double S, double K, double sig, double r, double T, int N)
{

	//dividend yield
	double y = 0.0;// 0.03;// 0.03; //div yield
	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	double Dx = sig * std::sqrt(2.0 * Dt);
	double v = r - y - 0.5 * sig * sig;

	//double  u = Dx;
	//double d = 1. / u;


	VecXX::SCALA_TYPE pu = -0.5 * Dt * ((sig * sig) / (Dx * Dx) + v / Dx);
	VecXX::SCALA_TYPE pd = -0.5 * Dt * ((sig * sig) / (Dx * Dx) - v / Dx);
	VecXX::SCALA_TYPE pm = 1. + Dt * (sig * sig) / (Dx * Dx) + r * Dt;


	std::vector<FLOAT> vdbg;
	//Pay off functions

	//call
	auto payOffFunc = [=](auto X) { return select(X > K, X - K, 0.0); };

	//put
	//auto payOffFunc = [=](auto X) { return select(X < K, K - X, 0.0); };

	//set up underlying asset prices at maturity
	double last = S * exp(-(N + 1) * Dx);
	double edx = exp(Dx);
	for (auto& el : terminalAssetPrices)
	{
		last *= edx;
		el = last;
	}

	//option value at maturity

	auto excerciseValue = transform(payOffFunc, terminalAssetPrices);


	auto american = [](auto X, auto Y) { return select(X > Y, X, Y); };


	//derivative boundary condition
	double  lambda_L = -1. * (terminalAssetPrices[1] - terminalAssetPrices[0]);
	double  lambda_U = 0.0;

	auto odd_slice = excerciseValue;
	//	vdbg = odd_slice;


	auto even_slice = odd_slice * 0.0;

	int J = 2 * N;
	int k = 0;

	VecXX pmp(1.0, J + 1);
	VecXX pp(1.0, J + 1);

	////////////
	// SOLVE IMPLICIT TRIDIAGONAL  IN LINE 
	//SUB BOUNDARY CONDITION AT J = -n INTO  J = -n+1
	pmp[1] = pm + pd;
	pp[1] = odd_slice[1] + pd * lambda_L;


	auto pu_pd = pu * pd;

	// eliminate upper diagonal
	for (int j = 2; j < J; ++j)
	{
		pmp[j] = pm - pu_pd / pmp[j - 1];
	}

	auto inv_pmp = 1.0 / pmp;

	auto pd_inv_pmp = pd * inv_pmp;

	/////////////

	for (; k < N; k += 2)
	{

		// SOLVE IMPLICIT TRIDIAGONAL  IN LINE 
		//SUB BOUNDARY CONDITION AT J = -n INTO  J = -n+1
		//pmp[1] = pm + pd;
		pp[1] = odd_slice[1] + pd * lambda_L;

		// eliminate upper diagonal
		for (int j = 2; j < J; ++j)
		{

			pp[j] = odd_slice[j] - pp[j - 1] * pd_inv_pmp[j - 1];
		}

		even_slice[1] = (pp[J - 1] + pmp[J - 1] * lambda_U) / (pu + pmp[J - 1]);
		even_slice[J - 1] = even_slice[J] - lambda_U;


		// back substitution
		for (int j = J - 2; j != 0; j--)
		{
			even_slice[j] = (pp[j] - pu * even_slice[j + 1]) * inv_pmp[j];
		}


		even_slice = transform(american, even_slice, (const VecXX&)excerciseValue);


		// calc odd slice now
//////////////////////////////////////////////////////////////////////


		pp[1] = even_slice[1] + pd * lambda_L;


		// eliminate upper diagonal
		for (int j = 2; j < J; ++j)
		{
			pp[j] = even_slice[j] - pp[j - 1] * pd_inv_pmp[j - 1];
		}

		odd_slice[1] = (pp[J - 1] + pmp[J - 1] * lambda_U) / (pu + pmp[J - 1]);
		odd_slice[J - 1] = odd_slice[J] - lambda_U;


		// back substitution
		for (int j = J - 2; j != 0; j--)
		{
			odd_slice[j] = (pp[j] - pu * odd_slice[j + 1]) * inv_pmp[j];
		}


		odd_slice = transform(american, odd_slice, (const VecXX&)excerciseValue);


	}

	return odd_slice[N];
}





double americanImplicitFiniteDiffPricer(double S, double K, double sig, double r, double T, int N)
{

	//dividend yield
	double y = 0.0;// 0.03;// 0.03; //div yield
	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	double Dx = sig * std::sqrt(2.0 * Dt);
	double v = r - y - 0.5 * sig * sig;

	//double  u = Dx;
	//double d = 1. / u;


	VecXX::SCALA_TYPE pu = -0.5 * Dt * ((sig * sig) / (Dx * Dx) + v / Dx);
	VecXX::SCALA_TYPE pd = -0.5 * Dt * ((sig * sig) / (Dx * Dx) - v / Dx);
	VecXX::SCALA_TYPE pm = 1. + Dt * (sig * sig) / (Dx * Dx) + r * Dt;


	std::vector<FLOAT> vdbg;
	//Pay off functions

	//call
	auto payOffFunc = [=](auto X) { return select(X > K, X - K, 0.0); };

	//put
	//auto payOffFunc = [=](auto X) { return select(X < K, K - X, 0.0); };

	//set up underlying asset prices at maturity
	double last = S * exp(-(N + 1) * Dx);
	double edx = exp(Dx);
	for (auto& el : terminalAssetPrices)
	{
		last *= edx;
		el = last;
	}

	//option value at maturity

	auto excerciseValue = transform(payOffFunc, terminalAssetPrices);

	//derivative boundary condition
	double  lambda_L = -1. * (terminalAssetPrices[1] - terminalAssetPrices[0]);
	double  lambda_U = 0.0;

	auto odd_slice = excerciseValue;
	vdbg = odd_slice;


	auto even_slice = odd_slice * 0.0;

	int J = 2 * N;
	int k = 0;

	VecXX pmp(1.0, J + 1);
	VecXX pp(1.0, J + 1);

	for (; k < N; k += 2)
	{

		// SOLVE IMPLICIT TRIDIAGONAL  IN LINE 
		//SUB BOUNDARY CONDITION AT J = -n INTO  J = -n+1
		pmp[1] = pm + pd;
		pp[1] = odd_slice[1] + pd * lambda_L;



		// eliminate upper diagonal
		for (int j = 2; j < J; ++j)
		{
			pmp[j] = pm - pu * pd / pmp[j - 1];
			pp[j] = odd_slice[j] - pp[j - 1] * pd / pmp[j - 1];
		}

		even_slice[1] = (pp[J - 1] + pmp[J - 1] * lambda_U) / (pu + pmp[J - 1]);
		even_slice[J - 1] = even_slice[J] - lambda_U;


		// back substitution
		for (int j = J - 2; j != 0; j--)
		{
			even_slice[j] = (pp[j] - pu * even_slice[j + 1]) / pmp[j];
		}

		/*

		// american
		for (int j = 0; j < (J+1); j++)
		{
			even_slice[j] = std::max(even_slice[j], excerciseValue[j]);
		}

*/
//even_slice =	transform(even_slice, excerciseValue);

// calc odd slice now
//////////////////////////////////////////////////////////////////////

		pmp[1] = pm + pd;
		pp[1] = even_slice[1] + pd * lambda_L;


		// eliminate upper diagonal
		for (int j = 2; j < J; ++j)
		{
			pmp[j] = pm - pu * pd / pmp[j - 1];
			pp[j] = even_slice[j] - pp[j - 1] * pd / pmp[j - 1];
		}

		odd_slice[1] = (pp[J - 1] + pmp[J - 1] * lambda_U) / (pu + pmp[J - 1]);
		odd_slice[J - 1] = odd_slice[J] - lambda_U;


		// back substitution
		for (int j = J - 2; j != 0; j--)
		{
			odd_slice[j] = (pp[j] - pu * odd_slice[j + 1]) / pmp[j];
		}

		/*
		//american
		for (int j = 0; j < (J + 1); j++)
		{
			odd_slice[j] = std::max(odd_slice[j], excerciseValue[j]);
		}
		*/

	}

	return odd_slice[N];
}





void doAmericanImplicitFiniteDiff()
{
	double time = 0;
	double res = 0.0;
	{
		TimerGuard timer(time);

		std::cout << std::setprecision(12) << " start \n";
		double S = 100.;
		double K = 100.;
		double vol = 0.2;
		double rate = 0.06;
		double T = 1;
		int N = 1000;
		//res = americanImplicitFiniteDiffPricer(S, K, vol, rate, T, N);		
		res = americanImplicitFiniteDiffPricerFast(S, K, vol, rate, T, N);
	}
	std::cout << std::setprecision(12) << "price " << res << "\n";
	std::cout << "americanImplicitFiniteDiffPricerFast takes secs" << time << "\n";
}





/////////////////////

//still broken ???
double americanCrankNicholsonPricer(double S, double K, double sig, double r, double T, int N)
{

	//dividend yield
	double y = 0.;// 0.03;// 0.0;// 0.03;// 0.03; //div yield
	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	double Dx = sig * std::sqrt(1.0 * Dt);// 0.2;//
	double v = r - y - 0.5 * sig * sig;

	//double  u = Dx;
	//double d = 1. / u;


	VecXX::SCALA_TYPE pu = -0.25 * Dt * ((sig * sig) / (Dx * Dx) + v / Dx);
	VecXX::SCALA_TYPE pd = -0.25 * Dt * ((sig * sig) / (Dx * Dx) - v / Dx);
	VecXX::SCALA_TYPE pm = 1. + 0.5 * Dt * (sig * sig) / (Dx * Dx) + 0.5 * r * Dt;

	std::vector<FLOAT> vdbg;
	//Pay off functions

	//call
	auto payOffFunc = [=](auto X) { return select(X > K, X - K, 0.0); };

	//put
	//auto payOffFunc = [=](auto X) { return select(X < K, K - X, 0.0); };

	//set up underlying asset prices at maturity
	double last = S * exp(-(N + 1) * Dx);
	double edx = exp(Dx);
	for (auto& el : terminalAssetPrices)
	{
		last *= edx;
		el = last;
	}

	//option vakue at maturity

	auto excerciseValue = transform(payOffFunc, terminalAssetPrices);

	//vdbg = terminalAssetPrices;
	//vdbg = excerciseValue;

	//derivative boundary condition
	double  lambda_L = -1. * (terminalAssetPrices[1] - terminalAssetPrices[0]);
	double  lambda_U = 0.0;

	auto odd_slice = excerciseValue;
	//	vdbg = odd_slice;



	auto even_slice = odd_slice * 0.0;

	even_slice[0] = odd_slice[0];

	int J = 2 * N;

	int k = 0;

	VecXX pmp(1.0, J + 1);
	VecXX pp(1.0, J + 1);

	for (; k <= N; k += 2)
	{

		// SOLVE IMPLICIT TRIDIAGONAL  IN LINE 
		//SUB BOUNDARY CONDITION AT J = -n INTO  J = -n+1
		pmp[1] = pm + pd;
		pp[1] = -pu * odd_slice[2] - (pm - 2.) * odd_slice[1] - pd * odd_slice[0] + pd * lambda_L;



		// eliminate upper diagonal
		for (int j = 2; j < J; ++j)
		{
			pmp[j] = pm - pu * pd / pmp[j - 1];
			pp[j] = -pu * odd_slice[j + 1] - (pm - 2.0) * odd_slice[j] - pd * odd_slice[j - 1] - pp[j - 1] * pd / pmp[j - 1];
		}



		even_slice[J] = (pp[J - 1] + pmp[J - 1] * lambda_U) / (pu + pmp[J - 1]);
		even_slice[J - 1] = even_slice[J] - lambda_U;


		// back substitution
		for (int j = J - 1; j >= 0; j--)
		{
			even_slice[j] = (pp[j] - pu * even_slice[j + 1]) / pmp[j];
		}


		even_slice[0] = odd_slice[0];
		//vdbg = even_slice;
/*
		//american condition
		for (int j = 0; j < (J+1); j++)
		{
			even_slice[j] = std::max(even_slice[j], excerciseValue[j]);
		}
*/

//	vdbg = even_slice;

//	vdbg = odd_slice;
//	vdbg = pmp;
//	vdbg = pp;


	// calc odd slice now
//////////////////////////////////////////////////////////////////////

		pmp[1] = pm + pd;
		pp[1] = -pu * even_slice[2] - (pm - 2.) * even_slice[1] - pd * even_slice[0] + pd * lambda_L;


		// eliminate upper diagonal
		for (int j = 2; j < J; ++j)
		{
			pmp[j] = pm - pu * pd / pmp[j - 1];
			pp[j] = -pu * even_slice[j + 1] - (pm - 2.0) * even_slice[j] - pd * even_slice[j - 1] - pp[j - 1] * pd / pmp[j - 1];
		}

		odd_slice[J] = (pp[J - 1] + pmp[J - 1] * lambda_U) / (pu + pmp[J - 1]);
		odd_slice[J - 1] = odd_slice[J] - lambda_U;


		// back substitution
		for (int j = J - 1; j >= 0; j--)
		{
			odd_slice[j] = (pp[j] - pu * odd_slice[j + 1]) / pmp[j];
		}

		odd_slice[0] = even_slice[0];
		/*
				//american condition
				for (int j = 0; j < (J + 1); j++)
				{
					odd_slice[j] = std::max(odd_slice[j], excerciseValue[j]);
				}
			*/

			//	vdbg = odd_slice;
			//	vdbg = even_slice;
			//	vdbg = pmp;
			//	vdbg = pp;

	}

	return odd_slice[N];
}




void doAmericanCrankNicholson()
{
	double time = 0;
	double res = 0.0;
	{
		TimerGuard timer(time);

		std::cout << std::setprecision(12) << " start \n";
		double S = 100.;
		double K = 100.;
		double vol = 0.2;
		double rate = 0.06;
		double T = 1;
		int N = 10000;
		res = americanCrankNicholsonPricer(S, K, vol, rate, T, N);
	}
	std::cout << std::setprecision(12) << "price " << res << "\n";
	std::cout << " takes secs" << time << "\n";
}




double americanTrinomialPricerUpAndOut(double S, double K, double sig, double r, double T, double H, double rebate, int N)
{

	double y = 0.0;// 0.03; //div yield


	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	double Dx = sig * std::sqrt(2.0 * Dt);
	double v = r - y - 0.5 * sig * sig;

	//double  u = Dx;
	//double d = 1. / u;


	VecXX::INS pu = 0.5 * ((Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx) + (v * Dt) / Dx);
	VecXX::INS pd = 0.5 * ((Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx) - (v * Dt) / Dx);
	VecXX::INS pm = 1. - (Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx);


	VecXX::INS disc = exp(-r * Dt);
	TrinomialSampler<VecXX::INS> sampler;

	auto trinomialRollBack = [=](TrinomialSampler<VecXX::INS>& sampler)
	{
		auto X1 = sampler.get<1>();
		auto X0 = sampler.get<0>();
		auto X_1 = sampler.get<-1>();
		return disc * (X1 * pu + X0 * pm + X_1 * pd);
	};


	//auto payOffFunc = [=](auto X) { return select(X > K, X - K, 0.0); }; //call
	auto payOffFunc = [=](auto X) { return select(X < K, K - X, 0.0); };  //put

	//set up underlying asset prices at maturity
	double last = S * exp(-(N + 1) * Dx);
	double edx = exp(Dx);
	for (auto& el : terminalAssetPrices)
	{
		last *= edx;
		el = last;
	}



	auto excerciseValue = transform(payOffFunc, terminalAssetPrices);
	auto odd_slice = excerciseValue;

	UnitarySampler<VecXX::INS> identity_sampler; //identity just 

	auto applyEarlyExcercise = [=](UnitarySampler<VecXX::INS>& sampler, auto excercisePrice)
	{
		auto optPrice = sampler.get<0>();
		return max(optPrice, excercisePrice);
	};


	auto applyBarrier = [=](UnitarySampler<VecXX::INS>& sampler, auto stockPrice)
	{
		auto optPrice = sampler.get<0>();
		return select(stockPrice < H, optPrice, rebate);
	};


	auto even_slice = odd_slice;

	int j = 2 * N + 1 - 1;
	int i = 0;
	for (; i < N; i += 2)
	{
		transform(odd_slice, even_slice, trinomialRollBack, sampler, i, j);
		// transform to get early excercise for american bit , iderntity sampler just passes values straight through
	//	transform(even_slice, excerciseValue, even_slice, applyEarlyExcercise, identity_sampler, i, j);
		transform(even_slice, terminalAssetPrices, even_slice, applyBarrier, identity_sampler, i, j);

		transform(even_slice, odd_slice, trinomialRollBack, sampler, i + 1, j - 1);
		//	transform(odd_slice, excerciseValue, odd_slice, applyEarlyExcercise, identity_sampler, i + 1, j - 1);
		transform(odd_slice, terminalAssetPrices, odd_slice, applyBarrier, identity_sampler, i + 1, j - 1);

		j -= 2;
	}

	return odd_slice[N];

	ignore(applyEarlyExcercise);
}




void doAmericanTrinomialPricerUpAndOut()
{

	double res = 0.0;
	double time = 0;
	{
		TimerGuard timer(time);

		double S = 100.;
		double K = 100.;
		double vol = 0.2;
		double rate = 0.06;
		double T = 1;
		int N = 1000;// 500;
		double H = 3.;
		double rebate = 0.0;
		res = americanTrinomialPricerUpAndOut(S, K, vol, rate, T, H, rebate, N);


	}
	std::cout << std::setprecision(12) << "price " << res << "\n";
	std::cout << " takes secs" << time << "\n";
}




VecXX  blackScholes(const VecXX& S, const VecXX& K, const VecXX& t, const VecXX& r, const VecXX& sigma)
{
	auto invK = 1.0 / K;
	auto discountedRate = exp(-r * t);
	auto S_invK = S * invK;
	auto log_sK = log(S_invK);
	auto rootT = sqrt(t);
	auto sigmaRootT = rootT * sigma;
	auto invSigmaRootT = 1.0 / sigmaRootT;
	auto halfSigmaSqrd_t = (0.5 * sigma * sigma + r) * t;
	auto d1 = invSigmaRootT * (log_sK + halfSigmaSqrd_t);
	auto d2 = d1 - sigmaRootT;
	auto normD1 = cdfnorm(d1);
	auto normD2 = cdfnorm(d2);
	auto C = S * normD1 - K * discountedRate * normD2;
	// auto delta = normD1;
	return C;
}




double euroTrinomialPricerWithInit(double S, double K, double sig, double r, double T, int N)
{

	double y = 0.0;// 0.03; //div yield


	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	double Dx = sig * std::sqrt(2.0 * Dt);
	double v = r - y - 0.5 * sig * sig;

	//double  u = Dx;
	//double d = 1. / u;


	VecXX::INS pu = 0.5 * ((Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx) + (v * Dt) / Dx);
	VecXX::INS pd = 0.5 * ((Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx) - (v * Dt) / Dx);
	VecXX::INS pm = 1. - (Dt * sig * sig + v * v * Dt * Dt) / (Dx * Dx);


	VecXX::INS disc = exp(-r * Dt);
	TrinomialSampler<VecXX::INS> sampler;

	auto trinomialRollBack = [=](TrinomialSampler<VecXX::INS>& sampler)
	{
		auto X1 = sampler.get<1>();
		auto X0 = sampler.get<0>();
		auto X_1 = sampler.get<-1>();
		return disc * (X1 * pu + X0 * pm + X_1 * pd);
	};

	//call
	auto payOffFunc = [=](auto X) { return select(X > K, X - K, 0.0); };

	//put
	//auto payOffFunc = [=](auto X) { return select(X < K, K -X , 0.0); };

	//set up underlying asset prices at maturity
	double last = S * exp(-(N + 1) * Dx);
	double edx = exp(Dx);
	for (auto& el : terminalAssetPrices)
	{
		last *= edx;
		el = last;
	}

	auto excerciseValue = transform(payOffFunc, terminalAssetPrices);
	auto odd_slice = excerciseValue;


	UnitarySampler<VecXX::INS> identity_sampler; //identity just 

	auto applyEarlyExcercise = [=](UnitarySampler<VecXX::INS>& sampler, auto excercisePrice)
	{
		auto optPrice = sampler.get<0>();
		return max(optPrice, excercisePrice);
	};


	auto even_slice = odd_slice;


	/// blacks initialisation 
	VecXX::INS invK = 1.0 / K;
	VecXX::INS discountedRate = exp(-r * Dt);

	VecXX::INS rootT = sqrt(Dt);
	VecXX::INS sigmaRootT = rootT * sig;
	VecXX::INS invSigmaRootT = 1.0 / sigmaRootT;
	VecXX::INS halfSigmaSqrd_t = (0.5 * sig * sig + r) * Dt;

	VecXX::INS Strike = K;

	auto blackScholeInit = [&](VecXX::INS S)
	{
		VecXX::INS S_invK = S * invK;
		VecXX::INS log_sK = log(S_invK);

		VecXX::INS d1 = invSigmaRootT * (log_sK + halfSigmaSqrd_t);
		VecXX::INS d2 = d1 - sigmaRootT;
		VecXX::INS normD1 = cdfnorm(d1);
		VecXX::INS normD2 = cdfnorm(d2);
		VecXX::INS C = S * normD1 - Strike * discountedRate * normD2;
		return C;

	};

	int j = 2 * N + 1 - 1;
	int i = 0;

	//use BS transform and normal for first pair of slices
	even_slice = transform(blackScholeInit, terminalAssetPrices);

	//std::vector<double> dbg = even_slice;

	transform(even_slice, odd_slice, trinomialRollBack, sampler, i + 1, j - 1);

	i += 2;
	j -= 2;

	for (; i < N; i += 2)
	{
		transform(odd_slice, even_slice, trinomialRollBack, sampler, i, j);
		// transform to get early excercise for american bit , iderntity sampler just passes values straight through
		//transform(even_slice, excerciseValue, even_slice, applyEarlyExcercise, identity_sampler, i, j);

		transform(even_slice, odd_slice, trinomialRollBack, sampler, i + 1, j - 1);
		//transform(odd_slice, excerciseValue, odd_slice, applyEarlyExcercise, identity_sampler, i + 1, j - 1);

		j -= 2;
	}

	return odd_slice[N];

	ignore(applyEarlyExcercise);
}




void doTrinomialPricerWithInit()
{

	double res = 0.0;
	double time = 0;
	{
		TimerGuard timer(time);

		double S = 100.;
		double K = 100.;
		double vol = 0.2;
		double rate = 0.06;
		double T = 1;
		int N = 10000;// 500;

		for (int k = 0; k < 100; k++)
		{
			res = euroTrinomialPricerWithInit(S, K, vol, rate, T, N);
		}


	}
	std::cout << std::setprecision(12) << "price " << res << "\n";
	std::cout << " takes secs" << time << "\n";
}





int main()
{

	try
	{
		doScan();
	}
	catch (...)
	{
	}
	//	return 0;


	///

	doTrinomialPricerWithInit();
	//return 0;


	doAmericanTrinomialPricerUpAndOut();
	//return 0;

	doAmericanImplicitFiniteDiff();
	//return 0;



	doAmericanCrankNicholson();
	//	return 0;



	doAmericanFiniteDiff();
	//	return 0;

	doTrinomialPricer();
	//	return 0;

	doZipping();
	//	return 0;

	doBinomialPricer();
	//	return 0;
	//

}