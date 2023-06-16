// curveExample.cpp : This file contains the 'main' function. Program execution begins and ends there.
//



#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/accumulate_transform.h"
#include "../Vectorisation/VecX/error_utils.h"


#include "../Vectorisation/VecX/zip_utils.h"
#include "../Vectorisation/VecX/span.h"


//#include "lattice_tools.h"
//#include "pricers.h"

/*

#include <iostream>
#include <algorithm>
#include <random>
#include <vector>

#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/apply_operation.h"
#include "../Vectorisation/VecX/vec_d.h"
#include "../Vectorisation/VecX/vec_bool.h"
#include "../Vectorisation/VecX/vec_view.h"


#include "../Vectorisation/VecX/target_name_space.h"
*/
//
#include "../Vectorisation/VecX/target_name_space.h"




#include <chrono>

#include "curve.h"



//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
using namespace DRC::VecD4D;
//using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;



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
#include <vector>


#include <iostream>




void do_curve()
{
	std::vector<double>  values = { 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 };
	std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10 };
	std::vector <double>   datesD = { 0,1,2,3,4,5,6,7,8,9,10 };

	std::vector< VecXX>  vecVals;
	for (int i = 0; i < 11; i++)
	{
		VecXX vals(i * 0.001 + 0.06, 200);
		vecVals.push_back(vals);

	}


	Curve<double, VecXX> forecastCurve; // just linear interpolation

	Curve2<double, VecXX, ZeroInterpCached<double, VecXX>>  discountCurve(200);
	discountCurve.setValues(begin(datesD), end(datesD), begin(vecVals), end(vecVals));

	auto valV = discountCurve.valueAt(0.0);

	auto valV2 = discountCurve.valueAt(0.5);

	std::vector<double> cashFlows(200, 0.5);







	double offsetDay = 0.0;
	for (long ll = 0; ll < 30; ll++) // so month days 
	{
		VecXX prices(0.0, 200);

		for (long l = 0; l < 30000; l++) //do all on the same month
		{
			int jj = 0;
			for (double d = 0.1; d < 10; d += 0.25, ++jj)
			{

				auto DF = discountCurve.valueAt(d + offsetDay);
				prices += DF * cashFlows[jj];
			}
			//auto valV3 = testCurve2.valueAt(0.5);
		}

		offsetDay += 0.003;
	} //ll


}


void doCurve2()
{


	std::vector<double>  values = { 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 };
	std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10 };
	std::vector <double>   datesD = { 0,1,2,3,4,5,6,7,8,9,10 };

	int KRVSZ = 200;// 00;
	std::vector< VecXX>  vecVals;
	for (int i = 0; i < 11; i++)
	{
		VecXX vals(i * 0.001 + 0.06, KRVSZ);
		//if (i > 0)
		//{

		for (int k = 1; k < KRVSZ; k++)
		{
			vals[k] += 0.0001 * k + 0.01 * i;
		}
		//}

		vecVals.push_back(vals);

	}


	Curve2<double, VecXX, ZeroInterpCached<double, VecXX>>  testCurve2(KRVSZ);
	testCurve2.setValues(begin(datesD), end(datesD), begin(vecVals), end(vecVals));

	auto valV = testCurve2.valueAt(0.0);

	auto valV2 = testCurve2.valueAt(0.5);

	//std::vector<double> cashFlows(200, 0.5);


	using cashFlow = std::pair< double, double>;

	using swap = cashFlow[40];

	// using  instruments = std::vector<swap>;

	using PF = std::vector<swap>;

	auto bigPF = PF(2000000);
	// double stepSz = 0.25;

	double count = 0;
	for (auto& swp : bigPF)
	{
		double date = 0.0;
		double period = 0.25;

		double offset = int(count / 30000.0) * 0.03 * 365.0;
		for (auto& cf : swp)
		{
			cf.first = date + offset;
			date += period;
			cf.second = 0.1;
		}
		count++;
	}

	auto results = new double[2000000 * KRVSZ];

	auto startTme = std::chrono::high_resolution_clock::now();


	//auto results=  new double[100000*200];

	//std::vector<VecXX> allPrices[1000000];

	// VecXX prices(0.0, 200);


	auto dfsFixed = testCurve2.valueAt(60);

	int j = 0;
	for (auto& swp : bigPF)
	{

		VecXX prices(0.0, KRVSZ);
		SpanXX prc(prices.begin(), KRVSZ);
		for (auto& cf : swp)
		{
			// auto DF = testCurve2.valueAt(cf.first);
			// prices += cf.second * DF;// testCurve2.valueAt(cf.first);// DF * cf.second + prices;


/*
			VecXX dfs = testCurve2.valueAt(cf.first);
			const  VecXX::SCALA_TYPE  cfs = cf.second;

			prices = FMA(dfs, cfs, prices);

*/
			const VecXX& dfs = testCurve2.valueAt(cf.first);


			SpanXX spnDf(dfs.begin(), dfs.size());

			VecXX::INS  cashflow = cf.second;

			auto priceAndAccumulateCashFlow = [cashflow](auto price, auto  df)
			{
				return mul_add(cashflow, df, price);
			};

			//prices =
		//	transformM(priceAndAccumulateCashFlow, prices, dfs);
			transform(priceAndAccumulateCashFlow, prc, spnDf, prc);



		}
		auto basePrice = prices[0];
		auto risk = prices - basePrice;
		//allPrices[j] =prices;
   //	 count++;

		/*	 */

		for (int kk = 0; kk < KRVSZ; kk++)
		{
			results[j * KRVSZ + kk] = risk[kk];
		}
		j++;




	}

	auto endTime = std::chrono::high_resolution_clock::now();
	auto runtime = endTime - startTme;
	std::cout << "run time 1M swaps " << KRVSZ << " buckets = " << runtime.count() / 1000000000.0 << std::endl;

}

template <typename numType, typename T1, typename T2>
struct MYFUNC : public T1, public T2
{
	numType F(const numType& val) { return T1(val); }
	numType DF(const numType& val) { return T2(val); }
};



void do_curve3()
{
	std::vector<double>  values = { 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 };
	std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10 };
	std::vector <double>   datesD = { 0,1,2,3,4,5,6,7,8,9,10 };

	std::vector< VecXX>  vecVals;
	for (int i = 0; i < 11; i++)
	{
		VecXX vals(i * 0.001 + 0.06, 200);
		vecVals.push_back(vals);

	}


	//Curve<double, VecXX> forecastCurve; // just linear interpolation

	Curve2<double, VecXX, ZeroInterpCached<double, VecXX>>  discountCurve(200);
	discountCurve.setValues(begin(datesD), end(datesD), begin(vecVals), end(vecVals));

	auto valV = discountCurve.valueAt(0.0);

	auto valV2 = discountCurve.valueAt(0.5);

	std::vector<double> cashFlows(200, 0.5);


	double offsetDay = 0.0;
	for (long ll = 0; ll < 30; ll++) // so month days 
	{
		VecXX prices(0.0, 200);

		for (long l = 0; l < 30000; l++) //do all on the same month
		{
			int jj = 0;
			for (double d = 0.1; d < 10; d += 0.25, ++jj)
			{

				auto DF = discountCurve.valueAt(d + offsetDay);
				prices += DF * cashFlows[jj];
			}
			//auto valV3 = testCurve2.valueAt(0.5);
		}

		offsetDay += 0.003;
	} //ll


}

void sensiCurve()
{
	std::vector<double>  values{ 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 };
	std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10 };
	//std::vector <double>   datesD = { 0,1,2,3,4,5,6,7,8,9,10 };

	std::vector< VecxD>  vecVals;
	for (int i = 0; i < 11; i++)
	{

		//VecxD vals((std::vector<double>(i*0.001 + 0.06, 200)), (std::vector<double>(0.0, 200)));
		//if (i > 0)
		auto vals = (i < 2) ? VecxD::makeDVecOnes(0.06, 200) : VecxD::makeDVecOnes(0.06, 200);
		auto valsZ = (i < 2) ? VecxD::makeDVecZero(0.06, 200) : VecxD::makeDVecZero(0.06, 200);

		for (int k = 1; k < 200; k++)
		{
			valsZ[k] += 0.0001 * k;
			vals[k] += 0.0001 * k;
		}
		//}

		//VecxD vals(vals

		if (i == 3)// || i== 4)
		{
			vecVals.push_back(vals);
		}
		else
		{
			vecVals.push_back(valsZ);
		}

	}

	auto datesD = values;
	Curve2<double, VecxD, ZeroInterpCached<double, VecxD>>  testCurve2(200);
	testCurve2.setValues(begin(datesD), end(datesD), begin(vecVals), end(vecVals));

	auto valV = testCurve2.valueAt(0.0);

	auto valV2 = testCurve2.valueAt(0.5);

	auto valV3 = testCurve2.valueAt(1.0);

	auto valV4 = testCurve2.valueAt(1.5);

	auto valV5 = testCurve2.valueAt(2.0);

	auto valV5a = testCurve2.valueAt(2.25);

	auto valV6 = testCurve2.valueAt(2.5);

	auto valV6a = testCurve2.valueAt(2.75);

	auto valV7 = testCurve2.valueAt(3.0);

	auto valV7a = testCurve2.valueAt(3.00001);

	auto valV8 = testCurve2.valueAt(3.5);

	auto valV9 = testCurve2.valueAt(3.75);

	auto valV10 = testCurve2.valueAt(4.0);

}




int main()
{
 //  std::cout << "Hello World!\n";
	doCurve2();
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
