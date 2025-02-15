#include "../Vectorisation/VecX/dr3.h"
#include "utils.h"
#include "pricers.h"



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

		auto X1 = sampler.X_Minus_1.value;
		auto X0 = sampler.X_0.value;
		auto X_1 = sampler.X_1.value;

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
		auto optPrice = sampler.X_0.value; //.get<0>();
		return max(optPrice, excercisePrice);
	};


	auto applyBarrier = [=](UnitarySampler<VecXX::INS>& sampler, auto stockPrice)
	{
		auto optPrice = sampler.X_0.value;//.get<0>();
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

	ignore(applyEarlyExcercise);

	return odd_slice[N];

}

