#include "../Vectorisation/VecX/dr3.h"
#include "utils.h"
#include "pricers.h"


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
		const auto& X1 = sampler.X_1.value;
		const auto& X0 = sampler.X_0.value;
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



