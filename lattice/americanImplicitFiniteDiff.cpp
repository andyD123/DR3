#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/zip_utils.h"
#include "utils.h"
#include "pricers.h"


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



		// american
		for (int j = 0; j < (J + 1); j++)
		{
			even_slice[j] = std::max(even_slice[j], excerciseValue[j]);
		}


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


		//american
		for (int j = 0; j < (J + 1); j++)
		{
			odd_slice[j] = std::max(odd_slice[j], excerciseValue[j]);
		}


	}

	return odd_slice[N];
}
