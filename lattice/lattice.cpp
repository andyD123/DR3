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


#include "lattice_tools.h"
#include "pricers.h"

                                                                   

/*
* 
Some initial implementations of  lattice models taken from the pseudo code in clewlow and strickland
https://www.amazon.co.uk/Implementing-Derivative-Models-Clewlow-1998-04-29/dp/B01K0S2UMC/ref=sr_1_1?crid=1X6ARB2YR7JAG&keywords=strickland+and+clewlow&qid=1685036394&sprefix=stricklamd+and+clewlow%2Caps%2C57&sr=8-1

to see how easy to code up with SIMD lambdas 

Sampler makes it easy to vectorise offset relationships  eg f ( X[i], x[i+1]) )

There is a need for iterating over multiple vexx  quantities at the same time so we need to introduce
some sort of zip iteration


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


void doBinomialPricer()
{

	double res = 0.0;
	double time = 0;

	double S = 100.;
	double K = 100.;
	double vol = 0.2;
	double rate = 0.06;
	double T = 1;
    int N = 30000;
	{
		TimerGuard timer(time);
		res = europeanBinomialPricer(S, K, vol, rate, T, N);
	}
	auto expected = blackScholes(S, K, T, rate, vol);
	std::cout << std::setprecision(12) << "binomial price " << res << " expected = " <<   expected.getScalarValue()  ;
	std::cout << "steps = "<< N  <<"  takes secs " << time << "\n";
}


void doTrinomialPricer()
{

	double res = 0.0;
	double time = 0;

    double S = 100.;
    double K = 100.;
    double vol = 0.2;
    double rate = 0.06;
    double T = 1;
    int N = 60000;// 500;
    {
        TimerGuard timer(time);

        //res = americanTrinomialPricer(S, K, vol, rate, T, N);
        res = europeanTrinomialPricer1(S, K, vol, rate, T, N);
    }

    auto expected = blackScholes(S, K, T, rate, vol);
    std::cout << std::setprecision(12) << "trinomial price " << res << " expected = " <<   expected.getScalarValue()  ;
    std::cout << "steps = "<< N  <<"  takes secs " << time << "\n";
}


void doAmericanFiniteDiff()
{
	double time = 0;
	double res = 0.0;

    double S = 100.;
    double K = 100.;
    double vol = 0.2;
    double rate = 0.06;
    double T = 1;
    int N = 60000;
    {
        TimerGuard timer(time);
		res = americanFiniteDiffPricer(S, K, vol, rate, T, N);
	}

    auto expected = blackScholes(S, K, T, rate, vol);
    std::cout << std::setprecision(12) << "americanFiniteDiff price " << res << " expected = " <<   expected.getScalarValue()  ;
    std::cout << "steps = "<< N  <<" takes secs " << time << "\n";
}


void doAmericanImplicitFiniteDiff()
{
	double time = 0;
	double res = 0.0;
    double S = 100.;
    double K = 100.;
    double vol = 0.2;
    double rate = 0.06;
    double T = 1;
    int N = 30000;
    {
        TimerGuard timer(time);
        res = americanImplicitFiniteDiffPricerFast(S, K, vol, rate, T, N);
	}
    auto expected = blackScholes(S, K, T, rate, vol);
    std::cout << std::setprecision(12) << "americanImplicitFiniteDiffPricer price " << res << " expected = " <<   expected.getScalarValue()  ;
    std::cout << "steps = "<< N  <<" takes secs " << time << "\n";
}


void doAmericanCrankNicholson()
{
	double time = 0;
	double res = 0.0;
    double S = 100.;
    double K = 100.;
    double vol = 0.2;
    double rate = 0.06;
    double T = 1;
    int N = 10000;
    {
        TimerGuard timer(time);
		res = americanCrankNicholsonPricer(S, K, vol, rate, T, N);
	}

    auto expected = blackScholes(S, K, T, rate, vol);
    std::cout << std::setprecision(12) << "americanCrankNicholsonPricer price " << res << " expected = " <<   expected.getScalarValue()  ;
    std::cout << "steps = "<< N  <<" takes secs " << time << "\n";


}

void doAmericanTrinomialPricerUpAndOut()
{

	double res = 0.0;
	double time = 0;
    double S = 100.;
    double K = 100.;
    double vol = 0.2;
    double rate = 0.06;
    double T = 1;
    int N = 1000;// 500;
    double H = 3.;
    double rebate = 0.0;
    {
        TimerGuard timer(time);
		res = americanTrinomialPricerUpAndOut(S, K, vol, rate, T, H, rebate, N);
	}
    auto expected = blackScholes(S, K, T, rate, vol);
    std::cout << std::setprecision(12) << "americanTrinomialPricerUpAndOut price " << res << " expected = " <<   expected.getScalarValue()  ;
    std::cout << "steps = "<< N  <<" takes secs " << time << "\n";
}


void doTrinomialPricerWithInit()
{
	double res = 0.0;
	double time = 0;
	double S = 100.;
	double K = 100.;
	double vol = 0.2;
	double rate = 0.06;
	double T = 1;

    int N = 3000;// 500;

    for (int k = 0; k < 100; k++)
    {
        TimerGuard timer(time);
        res = euroTrinomialPricerWithInit(S, K, vol, rate, T, N);
    }

    auto expected = blackScholes(S, K, T, rate, vol);
    std::cout << std::setprecision(12) << "euroTrinomialPricerWithInit price " << res << " expected = " <<   expected.getScalarValue()  ;
    std::cout << "steps = "<< N  <<" takes secs " << time << "\n";
}


int main()
{
	doBinomialPricer();

	doTrinomialPricer();

	doAmericanFiniteDiff();

	doTrinomialPricerWithInit();

	doAmericanImplicitFiniteDiff();

	doAmericanCrankNicholson();

	doAmericanTrinomialPricerUpAndOut();

    /*
	//lattice bits
	doZipping();

	doScan();
		
	doTransformWithASpan();

	doStridedSpan();

	doMatrix();
*/
	return 0;
	
	//

}