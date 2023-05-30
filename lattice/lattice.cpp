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
#include "../Vectorisation/VecX/span.h"



                                                                   

/*
* 
Some initial implementations of  lattice models taken from the psudo code in clewlow and strickland
https://www.amazon.co.uk/Implementing-Derivative-Models-Clewlow-1998-04-29/dp/B01K0S2UMC/ref=sr_1_1?crid=1X6ARB2YR7JAG&keywords=strickland+and+clewlow&qid=1685036394&sprefix=stricklamd+and+clewlow%2Caps%2C57&sr=8-1

to see how easy to code up with SIMD lambdas 

Sampler makes it easy to vectorise offset relationships  eg f ( X[i], x[i+1]) )

There is a need for iterating over multiple vexx  quantities at the same time so we need to introduce
so sort of zip iteration


We need a vectorised scan function  or something similar to achieve some sort of  vectorisation in tri diagonal solvers
// eg thomas. And schemes by Gilkes et al .  
// finally leading on to the need for representing multiple dimensional quantities 
// so we need a span and MD span like functionality  to work with PDE implementations
we also need to have some sort of span and padded/alligned  multi dimensional layout ie an MD span for SIMD
and also the layout 

however this is perhaps not the end. polyhedra, hilbert and koch curves 
and of course tiling .

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
//in utils
#include "utils.h"


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

	auto sum = [](auto X, auto Y) { return  X + Y; };

	double time = 0;
	{
		TimerGuard timer(time);

		VecXX result;
		for (long l = 0; l < 1000000; ++l)
		{
			result = ApplyScan(a, sum);
			//dbg = a;
		}

		dbg = result;
	}
	std::cout << "run time = " << time << "\n";



	///////TRANSFORM SCAN //////
	
	//{
	//	TimerGuard timer(time);
	//	auto doubleIt= [](auto X) { return  X*2.; };
	//	VecXX result;
	//	for (long l = 0; l < 1000000; ++l)
	//	{
	//		result = ApplyTransformScan(a, sum, doubleIt);
	//		//dbg = a;
	//	}
	//	dbg = result;
	//}

	//////////


	
	auto dbg_cpy = dbg;
	time = 0;
	{
		TimerGuard timer(time);
		
		for (long l = 0; l < 1000000; ++l)
		{
			std::inclusive_scan(begin(v1), end(v1), dbg.begin());
		}
	}
	std::cout << "run time std::scan = " << time << "\n";

	for (size_t i = 0; i < dbg_cpy.size(); ++i)
	{
		std::cout << i << "  ,  " << dbg[i]  << "   , " << dbg_cpy[i]  << std::endl;
	}

}

void testSampler()
{

	FLOAT sampleData[] = { 0., 1.0,2.,3.,4.,5.,6.,7.,8.,9.,10. };

	TrinomialSampler<VecXX::INS>  zzz;

	zzz.load(&sampleData[1]);

	auto x = zzz.get<-1>();
	auto y = zzz.get<0>();
	auto z = zzz.get<1>();

	ignore(x);
	ignore(y);
	ignore(z);
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


double europeanBinomialPricer(double S, double K, double sig, double r, double T, int N)
{

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

	double S = 100.;
	double K = 100.;
	double vol = 0.2;
	double rate = 0.06;
	double T = 1;
	{
		TimerGuard timer(time);



		int N = 30000;

		res = europeanBinomialPricer(S, K, vol, rate, T, N);

	}
	auto expected = blackScholes(S, K, T, rate, vol);
	std::cout << std::setprecision(12) << "binomial price " << res << "expected =" <<   expected.getScalarValue()  << "\n";
	std::cout << " takes secs" << time << "\n";
}

double europeanTrinomialPricer(double S, double K, double sig, double r, double T, int N)
{

	double y = 0.0;// 0.03; //div yield


	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	double Dx = sig * std::sqrt(3.0* Dt);
	
	double v = r - y - 0.5 * sig * sig;



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

	UnitarySampler<VecXX::INS> identity_sampler; //identity just  passes through

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



double europeanTrinomialPricer1(double S, double K, double sig, double r, double T, int N)
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

//	auto applyEarlyExcercise = [=](UnitarySampler<VecXX::INS>& sampler, auto excercisePrice)
//	{
//		auto optPrice = sampler.get<0>();
//		return max(optPrice, excercisePrice);
//	};


	auto even_slice = odd_slice;

	int j = 2 * N + 1 - 1;
	int i = 0;
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
		int N = 60000;// 500;
		//res = americanTrinomialPricer(S, K, vol, rate, T, N);
		res = europeanTrinomialPricer1(S, K, vol, rate, T, N);


	}
	std::cout << std::setprecision(12) << "price " << res << "\n";
	std::cout << "trinomial takes secs" << time << "\n";
}

double americanFiniteDiffPricer(double S, double K, double sig, double r, double T, int N)
{

	//dividend yield
	double y = 0.0;
	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	double Dx = sig * std::sqrt(2.0 * Dt);
	double v = r - y - 0.5 * sig * sig;


	VecXX::INS pu = 0.5 * Dt * ((sig * sig) / (Dx * Dx) + v / Dx);
	VecXX::INS pd = 0.5 * Dt * ((sig * sig) / (Dx * Dx) - v / Dx);
	VecXX::INS pm = 1. - Dt * (sig * sig) / (Dx * Dx) - r * Dt;



	TrinomialSampler<VecXX::INS> sampler;
	//introduces offset variables so that we can get vectorised versions
	// of X[i+1], X[i] and  x[i-1]
	// under the hood these do unaligned loads into registers so taht we can still
	// apply the vectorised 
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

	UnitarySampler<VecXX::INS> identity_sampler; //identity  pas through

	//this is the american part of the option excercise
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

		//apply boundary condition
		even_slice[0] = even_slice[1] + terminalAssetPrices[1] - terminalAssetPrices[0];
		even_slice[J] = even_slice[J - 1];
		// transform to get early excercise for american exercise , iderntity sampler just passes values straight through
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
		int N = 60000;
		res = americanFiniteDiffPricer(S, K, vol, rate, T, N);


	}
	std::cout << std::setprecision(12) << "price " << res << "\n";
	std::cout << "finite dif takes secs" << time << "\n";
}

double americanImplicitFiniteDiffPricerFast(double S, double K, double sig, double r, double T, int N)
{

	
	double y = 0.0;//dividend yield

	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	double Dx = sig * std::sqrt(2.0 * Dt);
	double v = r - y - 0.5 * sig * sig;


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
	//LOOP HOIST BITS FROM  IMPLICIT TRIDIAGONAL 

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

		// SOLVE IMPLICIT TRIDIAGONAL  IN LINE 	SUB BOUNDARY CONDITION AT J = -n INTO  J = -n+1
	
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

		//american excercise bit
		even_slice = transform(american, even_slice, (const VecXX&)excerciseValue);


		// now calculate the  odd slice 

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

		//american excercise bit
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
		int N = 30000;
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


	//derivative boundary condition
	double  lambda_L = -1. * (terminalAssetPrices[1] - terminalAssetPrices[0]);
	double  lambda_U = 0.0;

	auto odd_slice = excerciseValue;

	//set up slices
	auto even_slice = odd_slice * 0.0;
	even_slice[0] = odd_slice[0];

	int J = 2 * N;
	int k = 0;

	VecXX pmp(1.0, J + 1);
	VecXX pp(1.0, J + 1);

	for (; k <= N; k += 2)
	{

		// SOLVE IMPLICIT TRIDIAGONAL  IN LINE //SUB BOUNDARY CONDITION AT J = -n INTO  J = -n+1
		
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


double euroTrinomialPricerWithInit(double S, double K, double sig, double r, double T, int N)
{

	double y = 0.0;// 0.03; //div yield
	VecXX terminalAssetPrices(1.0, 2 * N + 1);

	double Dt = T / N;
	double Dx = sig * std::sqrt(2.0 * Dt);
	double v = r - y - 0.5 * sig * sig;

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

	//even_slice = transform(payOffFunc, terminalAssetPrices);

	//std::vector<double> dbg = even_slice;

	transform(even_slice, odd_slice, trinomialRollBack, sampler, i + 1, j - 1);

	i += 2;
	j -= 2;

	for (; i < N; i += 2)
	{
		transform(odd_slice, even_slice, trinomialRollBack, sampler, i, j);
		// transform to get early excercise for american bit , identity sampler just passes values straight through
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

	//warm up
	double res = 0.0;
	{
		
		
		{

			double S = 100.;
			double K = 100.;
			double vol = 0.2;
			double rate = 0.06;
			double T = 1;
			int N = 300;// 500;

			for (int k = 0; k < 100; k++)
			{
				res = euroTrinomialPricerWithInit(S, K, vol, rate, T, N);
			}
		}
	}


	//double res = 0.0;
	double time = 0;


	double S = 100.;
	double K = 100.;
	double vol = 0.2;
	double rate = 0.06;
	double T = 1;

	{
		TimerGuard timer(time);


		int N = 3000;// 500;

		for (int k = 0; k < 100; k++)
		{
			res = euroTrinomialPricerWithInit(S, K, vol, rate, T, N);
		}


	}

	auto expected = blackScholes(S, K, T, rate, vol);
	std::cout << std::setprecision(12) << "price " << res << "   expected =" << expected.getScalarValue() << "\n";

	//std::cout << std::setprecision(12) << "price " << res << "\n";
	std::cout << " takes secs" << time << "\n";
}



void doMatrix()
{
	VecXX owningVec(0.1, 16*10);

	auto val = 0.0;
	for (auto& x : owningVec)
	{
		x = val;
		val++;
	}

	using MAT1 =  DR3::Layout2D<double, 8, 0>;

	//auto pDat =
	DR3::MDSpan<double,MAT1> mat(owningVec.data(), 10, 10);

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

}



int main()
{
/*
	doBinomialPricer();
	//	return 0;
	doTrinomialPricer();
	//return 0;
	doAmericanFiniteDiff();

	doTrinomialPricerWithInit();
	//return 0;  

	doAmericanImplicitFiniteDiff();
	//return 0;

	doAmericanCrankNicholson();
	//return 0;

	return 0;



   //doAmericanTrinomialPricerUpAndOut();
	//return 0;


	doZipping();
	//	return 0;
*/
 
	doMatrix();

	return 0;


	try
	{
		doScan();
	}
	catch (...)
	{
	}
		return 0;

	
	


	//

}