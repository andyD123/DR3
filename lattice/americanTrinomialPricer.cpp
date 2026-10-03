
#include "../Vectorisation/VecX/dr3.h"
#include "utils.h"
#include "pricers.h"
#include <stdexcept>


double americanTrinomialPricer(double S, double K, double sig, double r, double T, int N)
{
	if (N <= 0 || (N % 2) != 0)
		throw std::invalid_argument("americanTrinomialPricer requires a positive even N");

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
		const auto& X1 = sampler.X_1.value;
		const auto& X0 = sampler.X_0.value;
		const auto& X_1 = sampler.X_Minus_1.value;
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
		auto optPrice = sampler.X_0.value; //.get<0>();
		return max(optPrice, excercisePrice);
	};


	auto even_slice = odd_slice;

	int j = 2 * N + 1;
	int i = 0;
	for (; i < N; i += 2)
	{
		transform(odd_slice, even_slice, trinomialRollBack, sampler, i, j);
		// The trinomial sampler shrinks the valid range by one node at each edge.
		transform(even_slice, excerciseValue, even_slice, applyEarlyExcercise, identity_sampler, i + 1, j - 1);

		transform(even_slice, odd_slice, trinomialRollBack, sampler, i + 1, j - 1);
		transform(odd_slice, excerciseValue, odd_slice, applyEarlyExcercise, identity_sampler, i + 2, j - 2);

		j -= 2;
	}

	return odd_slice[N];
}
