#include "../Vectorisation/VecX/dr3.h"
#include "utils.h"
#include "pricers.h"


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

