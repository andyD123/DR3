#include "../Vectorisation/VecX/dr3.h"

#include "utils.h"
#include "pricers.h"


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


	UnitarySampler<VecXX::INS> identity_sampler; //identity just 

	auto applyEarlyExcercise = [=](UnitarySampler<VecXX::INS>& sampler, auto excercisePrice)
	{
		auto optPrice = sampler.X_0.value; //.get<0>();
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

	ignore(applyEarlyExcercise);

	return odd_slice[N];

}

