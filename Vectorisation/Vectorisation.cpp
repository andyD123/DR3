// Vectorisation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
/*
//#include "pch.h"
#include <iostream>




#include "../Vectorisation/VecX/Curve.h"
#include "../Vectorisation/VecX/Vec.h"
#include "../Vectorisation/VecX/operations.h"
//#include "operations.h"

#include <chrono>


//typedef Vec<VecDouble>  Vecx;
typedef VecD<VecDouble>  VecxD;

typedef VecD<VecDouble>  Vecx;
typedef Vec<VecDouble>  VecXX;



Vecx blackScholes2(const Vecx& spot, const Vecx& strike, const Vecx& time, const Vecx& rate, const Vecx& vol)
{


	auto K = strike;
	auto S = spot;
	auto sigma = vol;
	auto r = rate;


	//	std::cout << "spot" <<std::endl;
	//	print(spot);
	//	std::cout << "sigma" <<std::endl;
	//	print(sigma);

	auto invK = 1.0 / K;

	auto discountedRate = exp(-rate * time);

	auto log_sK(log(S*invK));
	//auto logSK ( log_s + log(K*discountedRate) );

	auto rootT = sqrt(time);
	auto sigmaRootT = rootT * sigma;
	auto invSigmaRootT =1.0/ sigmaRootT;

	auto d1 = invSigmaRootT * (log_sK + (0.5 *sigma*sigma + r)*time);// / sigmaRootT;// *invSigmaRootT; //logSK/sigmaRootT +0.5f*sigmaRootT;
	auto d2 = d1 - sigmaRootT;

	auto normD1 = cdfnorm(d1);
	auto normD2 = cdfnorm(d2);

	auto C = S * normD1 - K * discountedRate* normD2;

	//	std::cout << "C" <<std::endl;
	//	print(C);
	auto delta = discountedRate * normD1;

	return C;

}



void call_BS()
{

	size_t size = 20;// 1024;
	Vecx spot(size);
	double strike = 40.0;
	double time = 0.5;
	double rate = 0.1;
	Vecx vol(size);


	for (size_t i = 0; i < size; i++)
	{
		vol[i] = 0.20;//+0.001* static_cast<double>(i-500);

		spot[i] = static_cast<double>(42.0);
	}


	auto fixedVol = 0.25;
	auto fixedSpot = 100.0;



	for (long i = 0; i < 1000000; i++)
	{
		//volatile Vecx prices = blackScholes(  spot,  strike,  time,  rate,  vol);
		volatile Vecx prices = blackScholes2(spot, Vecx(strike), Vecx(time), Vecx(rate), vol);//Vecx(fixedVol) );//
	}
	std::cout << "done !\n";

}


void call_BSDeriv()
{

	size_t size = 20;// 1024;
	VecXX spot(size);
	double strike = 40.0;
	double time = 0.5;
	double rate = 0.1;
	VecXX vol(size);

	std::vector<double> ones(20,1.0);

	VecXX Vec1(ones);


	std::vector<double> two(20, 42.0);
	VecXX Vec2(two);


	VecXX strikeV(Vec2);


	std::vector<double> timeV(20, time);
	VecXX timeVec(timeV);

	std::vector<double> rateV(20, rate);
	VecXX rateVec(rateV);

	std::vector<double> volV(20, 0.2);
	VecXX volVec(volV);

	auto spotD = VecD< VecDouble>::makeDVecOnes(Vec2);
	auto strikeD = VecD< VecDouble>::makeDVecZero(strikeV);
	auto timeD = VecD< VecDouble>::makeDVecZero(timeVec);
	auto rateD = VecD< VecDouble>::makeDVecZero(rateVec);

	auto volD = VecD< VecDouble>::makeDVecZero(volVec);// VecXX(vol));


	for (size_t i = 0; i < size; i++)
	{
		vol[i] = 0.20;//+0.001* static_cast<double>(i-500);

		spot[i] = static_cast<double>(42.0);
	}


	auto fixedVol = 0.2;
	auto fixedSpot = 100.0;



//	for (long i = 0; i < 1000000; i++)
//	{
		//volatile Vecx prices = blackScholes(  spot,  strike,  time,  rate,  vol);
	//	volatile Vecx prices = blackScholes2(spot, Vecx(strike), Vecx(time), Vecx(rate), vol);//Vecx(fixedVol) );//
		volatile Vecx prices = blackScholes2(spotD, strikeD, timeD, rateD, volD);//Vecx(fixedVol) );
//	}
	std::cout << "done !\n";

}



void diffFunction()
{
	VecxD testVec(1.0, 0.0);

	std::vector<double> ones;
	std::vector<double> twos;
	double SZ = 20;
	for (double i = 0; i < SZ; ++i)
	{
		ones.push_back(i);
		twos.push_back(2.0);
	}


	VecxD Vec1(ones);
	VecxD Vec2(twos);

	auto res = Vec1 + Vec2;

	VecxD scalar(3.0,11.1);

	auto res2 = Vec1 + scalar;

	auto res3 = Vec1 + 2.666;

	auto res4 = 2.666 + Vec1;

	auto res5 = Vec1 - Vec2;
	auto res6 = Vec1 - scalar;
	auto res7 = Vec1 - 3.142;
	auto res8 = 1.1234 - Vec1;
	res8 += Vec1;
	res8 -= Vec1;


	auto res9 = Vec1 * Vec2;
	auto res10 = Vec1 * 10.1;
	res10 *= Vec1;
	res10 *=10.1;
	//res10


	auto res11 = Vec1 / Vec2;
	auto res12 = Vec1 / 10.1;
	res11 /= Vec1;
	res11 /= 10.1;

	auto negV1 = -Vec1;

	auto root = sqrt(Vec1);

	auto rEX = exp(Vec1);

	auto Lg = log(Vec1);

	std::vector<double> X = { 1.0,2.0 };
	std::vector<double> Y = { 1.0,1.0 };

	VecxD tv(X, Y);


	VecxD r = sqrt(tv);

	r = tv*tv;

	std::cout << r.value()[0] << "," <<  r.value()[1] << "\n";
	std::cout << r.derivative()[0] << ","  << r.derivative()[1] << "\n";


	r = tv *10;

	std::cout << r.value()[0] << "," << r.value()[1] << "\n";
	std::cout << r.derivative()[0] << "," << r.derivative()[1] << "\n";


	r = 10*tv;
	std::cout << r.value()[0] << "," << r.value()[1] << "\n";
	std::cout << r.derivative()[0] << "," << r.derivative()[1] << "\n";

	r = 1 / tv;
	std::cout << r.value()[0] << "," << r.value()[1] << "\n";
	std::cout << r.derivative()[0] << "," << r.derivative()[1] << "\n";



	r =  tv +10;
	std::cout << r.value()[0] << "," << r.value()[1] << "\n";
	std::cout << r.derivative()[0] << "," << r.derivative()[1] << "\n";



	r = exp(tv);
	std::cout << r.value()[0] << "," << r.value()[1] << "\n";
	std::cout << r.derivative()[0] << "," << r.derivative()[1] << "\n";

}



int main()
{

//	call_BSDeriv();

	call_BS();
	diffFunction();

    std::cout << "done \n"; 
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

*/