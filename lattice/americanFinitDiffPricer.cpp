
#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/accumulate_transform.h"
#include "../Vectorisation/VecX/error_utils.h"


#include "../Vectorisation/VecX/zip_utils.h"
#include "../Vectorisation/VecX/span.h"


#include "utils.h"
#include "pricers.h"


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
		const auto& X1 = sampler.X_1.value;
		const auto& X0 = sampler.X_0.value;
		const auto& X_1 = sampler.X_Minus_1.value;
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
		auto optPrice = sampler.X_0.value; //.get<0>();
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
