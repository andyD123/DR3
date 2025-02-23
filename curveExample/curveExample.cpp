// curveExample.cpp : This file contains the 'main' function. Program execution begins and ends there.
//



#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/accumulate_transform.h"
#include "../Vectorisation/VecX/error_utils.h"


#include "../Vectorisation/VecX/zip_utils.h"
#include "../Vectorisation/VecX/span.h"


#include "../Vectorisation/VecX/target_name_space.h"


#include <memory>

#include <chrono>

#include <random>
#include <algorithm>  // For std::shuffle

#include "curve.h"



//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
//using namespace DRC::VecD4D;
using namespace DRC::VecD8D;
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





////////////////  cashflow stuff ////////////////////


using namespace std;
using namespace std::chrono;

// Struct to hold the cash flow information: amount and cumulative day
struct CashFlow {
	double amount;  // Coupon (or final payment including principal)
	long day;       // date == days since 2000 //Cumulative days since the bond's issue date
	long date;
};

// Struct to represent a Bond using Array of Structures (AOS)
struct Bond_AOS {
	long m_ID;                          // Unique bond identifier
	system_clock::time_point m_issueDate; // Actual issue date of the bond
	std::vector<CashFlow> cashflows;      // Cash flows (annual coupon payments; first period skipped)

	// Less-than operator: compares bonds by the day-of-year of their last cash flow.
	bool operator<(const Bond_AOS& other) const {
		int thisLast = get_day_of_year(cashflows.back().day, m_issueDate);
		int otherLast = get_day_of_year(other.cashflows.back().day, other.m_issueDate);
		return thisLast < otherLast;
	}

	// Returns the maturity day-of-year for this bond.
	int get_maturity_day_of_year() const {
		return get_day_of_year(cashflows.back().day, m_issueDate);
	}

private:
	// Helper: converts a cumulative day (since issue) and the bond's issue date into a day-of-year.
	int get_day_of_year(long cumulative_day, system_clock::time_point issueDate) const {
		system_clock::time_point target_date = issueDate + days(cumulative_day);
		time_t tt = system_clock::to_time_t(target_date);
		tm local;
		localtime_s(&local, &tt);
		return local.tm_yday + 1; // tm_yday is 0-based
	}
};

// Utility: creates a tm structure from year, month, and day.
tm makeDate(int year, int month, int day) {
	tm date = {};
	date.tm_year = year - 1900;
	date.tm_mon = month - 1;
	date.tm_mday = day;
	return date;
}



// Returns the number of days between two time_points.
int calculate_days_between(system_clock::time_point start, system_clock::time_point end) {
	return (int)duration_cast<days>(end - start).count();
}

static system_clock::time_point baseDate()
{
	static auto date = makeDate(1900, 1, 1);
	static time_t tt = mktime(&date);
	static auto baseDate = std::chrono::system_clock::from_time_t(tt);
	return baseDate;
}

// Generates a Bond_AOS with annual coupon payments.
// The parameter 'maturity_day' is an offset (in days) to set the first payment date.
Bond_AOS generate_bond_cashflows_AOS(double face_value, double coupon_rate, int maturity_years,
	system_clock::time_point issue_date, long bond_id)//, int maturity_day)
{
	Bond_AOS bond;
	bond.m_ID = bond_id;

	bond.m_issueDate = issue_date; // Save the bond's issue date

	double annual_coupon = coupon_rate * face_value;

	auto next_payment_date = issue_date;

	// Skip the first period (period 1) and add cash flows starting from period 2.
	for (int period = 1; period <= maturity_years; ++period)
	{
		int days_to_cashflow = calculate_days_between(issue_date, next_payment_date);
		auto cfDate = bond.m_issueDate + years(period);//bond.m_issueDate + days(days_to_cashflow);
		long cfDate_l = calculate_days_between(baseDate(), cfDate);
		double coupon_payment = annual_coupon; // Payment for this period

		if (period > 1)
		{ // Skip first period
			CashFlow cf;
			cf.amount = (period == maturity_years) ? coupon_payment + face_value : coupon_payment;
			cf.day = days_to_cashflow;
			cf.date = cfDate_l;
			bond.cashflows.push_back(cf);
		}
		//issue_date = next_payment_date;
		next_payment_date = issue_date + years(period);
	}
	return bond;
}

// Creates a portfolio of bonds grouped by maturity day offset.
// For each predetermined maturity day offset, generates multiple bonds.
map<int, vector<Bond_AOS>> createPortfolio(int bondsPerMaturity, int numberOfYears)
{
	map<int, vector<Bond_AOS>> portfolio_by_maturity_day;

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> couponDis(0.02, 0.10);


	int year = 2000;
	int month = 1;
	int day = 1;
	double face_value = 1000;
	int maturity_years = numberOfYears;

	auto baseDate = makeDate(1900, 1, 1);



	tm issue_tm = makeDate(year, month, day);
	system_clock::time_point issue_date = system_clock::from_time_t(mktime(&issue_tm));

	long bondID = 1;
	for (int md = 0; md < 365; md++)
	{


		auto issue_date_this_bond = issue_date + days(md);
		for (int j = 0; j < bondsPerMaturity; ++j)
		{
			double coupon_rate = couponDis(gen);
			Bond_AOS bond = generate_bond_cashflows_AOS(face_value, coupon_rate, maturity_years, issue_date_this_bond, bondID);//, md);
			portfolio_by_maturity_day[md].push_back(bond);
			++bondID;
		}
	}
	return portfolio_by_maturity_day;
}

// Displays only the first and last cash flow for a bond,
// printing the corresponding calendar date (computed from the bond's issue date).
void display_bond_cashflows_AOS(const Bond_AOS& bond) 
{
	// Use the bond's issue date as base.
	time_t base_time = system_clock::to_time_t(bond.m_issueDate);
	tm base_tm;
	localtime_s(&base_tm, &base_time);
	system_clock::time_point base_date = system_clock::from_time_t(mktime(&base_tm));

	cout << fixed << setprecision(2);
	if (bond.cashflows.empty()) {
		cout << "Bond ID: " << bond.m_ID << " has no cash flows." << endl;
	}
	else if (bond.cashflows.size() == 1) {
		const CashFlow& cf = bond.cashflows.front();
		system_clock::time_point cf_date = bond.m_issueDate + days(cf.day);
		time_t tt = system_clock::to_time_t(cf_date);
		tm local_tm;
		localtime_s(&local_tm, &tt);
		cout << "Bond ID: " << bond.m_ID
			<< " | Only Cashflow -> Amount: £" << cf.amount
			<< " | Cumulative Day: " << cf.day
			<< " | Date: " << put_time(&local_tm, "%Y-%m-%d") << endl;
	}
	else {
		// Print first cashflow
		const CashFlow& cf_first = bond.cashflows.front();
		//     system_clock::time_point cf_date_first = bond.m_issueDate + days(cf_first.day);

		system_clock::time_point cf_date_first = baseDate() + days(cf_first.date);

		time_t tt_first = system_clock::to_time_t(cf_date_first);
		tm local_tm_first;
		localtime_s(&local_tm_first, &tt_first);
		cout << "Bond ID: " << bond.m_ID
			<< " | First Cashflow -> Amount: £" << cf_first.amount
			<< " | Cumulative Day: " << cf_first.day
			<< " | Date: " << put_time(&local_tm_first, "%Y-%m-%d") << endl;
		// Print last cashflow
		const CashFlow& cf_last = bond.cashflows.back();
		system_clock::time_point cf_date_last = baseDate() + days(cf_last.date);
		time_t tt_last = system_clock::to_time_t(cf_date_last);
		tm local_tm_last;
		localtime_s(&local_tm_last, &tt_last);
		cout << "Bond ID: " << bond.m_ID
			<< " | Last Cashflow -> Amount: £" << cf_last.amount
			<< " | Cumulative Day: " << cf_last.day
			<< " | Date: " << put_time(&local_tm_last, "%Y-%m-%d") << endl;
	}
	cout << "Maturity Day of the Year: " << bond.get_maturity_day_of_year() << endl;
}



int example_cashflows() 
{
	int bondsPerMaturity = 3;// 3000;// 10000;// 3;  // Number of bonds per maturity day offset.
	int numberOfYears = 20;    // Each bond matures in 20 years.

	// Create the portfolio of bonds grouped by maturity day offset.
	map<int, vector<Bond_AOS>> portfolio = createPortfolio(bondsPerMaturity, numberOfYears);

	cout << "Portfolio grouped by maturity day offset:\n";

	for (const auto& entry : portfolio)
	{
		cout << "\nMaturity Day Offset: " << entry.first << "\n";
		for (const auto& bond : entry.second)
		 {
			   cout << "\nBond ID: " << bond.m_ID << "\n";
		       display_bond_cashflows_AOS(bond);
		 }
	}
	
	return 0;
}



/////////////////

/*

struct CashFlow {
	double amount =1.0;  // Coupon (or final payment including principal)
//	long day = 0;       // Cumulative days since the bond's issue date
	long date = 0;
};
*/


struct OldCurveImpl
{
	virtual ~OldCurveImpl() {};
	virtual double getDiscountFactor(long days) const = 0;
};

class OldCurve
{
public:
	virtual ~OldCurve() {};
	OldCurve( std::unique_ptr<OldCurveImpl>&& newCurve) :m_pImpl(std::move(newCurve)) {};

	// no copy/assign 

	double getDiscountFactor( long days) const
	{
		return m_pImpl->getDiscountFactor(days);
	}


private:

	std::unique_ptr<OldCurveImpl>  m_pImpl;
};


class SimpleInterp : public  OldCurveImpl
{
public:
	~SimpleInterp() = default;

	SimpleInterp(const std::vector<long>& dts, const std::vector<double>& zeros) :zeroRates(zeros), dates(dts) {} //check arrays same size, dates monatonic and rates sensible

	double getDiscountFactor(long days) const override
	{
		if (true) //inRange(days))
		{

			//small number of curve points 
			size_t largestPos = dates.size();
			size_t i = 1;
			for (; i < largestPos; ++i)
			{
				if (dates[i] >= days) break;
			}
			if ( i >= largestPos)
			{
				i = largestPos - 1;
			}

			double fraction =(days - dates[i - 1]) / (dates[i] - dates[i - 1]);

			double zeroRate =(zeroRates[i] - zeroRates[i - 1])* fraction + zeroRates[i - 1];

			double years = days / 365.;

			return exp(-zeroRate * years);

		}

	}

	

private:
	std::vector<double> zeroRates;
	std::vector<long> dates;


};


OldCurve makeCurve(const std::vector<long>& dts, const std::vector<double>& zeros)
{
	auto innerCurve =	std::make_unique<SimpleInterp>(dts, zeros);
	return OldCurve(std::move(innerCurve) );
}

/*
//SOA 
struct CashFlowsSOA
{
	CashFlowsSOA(long start, double repeat_period_days, int number, double amount) : amounts(number, amount)
	{
	  std::vector<double>   days_tmp(number, repeat_period_days);
	  std::inclusive_scan(begin(days_tmp), end(days_tmp), begin(days_tmp));

	 std::transform(begin(days_tmp), end(days_tmp), begin(days_tmp), std::roundf);

	  days.resize(number);
	  for (const auto& dy : days_tmp)
	  {
		  days.emplace_back(static_cast<long>(dy) );
	  }
		
	}

	std::vector<double> amounts;
	
	std::vector<long>   days;
};


CashFlowsSOA makeCashFlowsSOA(long start, double repeat_period_days, int number, double amount)
{
	return CashFlowsSOA(start, repeat_period_days, number, amount);
}

//using CashFlow = std::pair< long, double>;

std::vector<CashFlow>   makeCashFlowsAOS(long start, double repeat_period_days, int number, double amount)
{
	CashFlowsSOA temp(start, repeat_period_days, number, amount);

	std::vector<CashFlow> flows(number);
	
	for (int i = 0; i < number; ++i)
	{
		flows[i].date = temp.days[i];
		flows[i].amount = temp.amounts[i];
	}

	return flows;
	
}
*/




std::vector<CashFlow>   makeCashFlowsAOS(int Y,int M, int D, int maturity_years, double coupon_rate , double face)
{
	auto date = makeDate(Y, M, D);
	time_t tt = mktime(&date);
	auto issueDate = std::chrono::system_clock::from_time_t(tt);
	Bond_AOS bnd = generate_bond_cashflows_AOS(100.0, coupon_rate,  maturity_years, issueDate,  0);

	return bnd.cashflows;
}



void do_basic_curve()
{

	//auto baseDate = makeDate(1900, 1, 1);

	auto refDate2K = makeDate(2000, 1, 1);
	time_t tt = mktime(&refDate2K);
	auto refDate = std::chrono::system_clock::from_time_t(tt);

	long daysOffset = calculate_days_between(baseDate(), refDate);

	std::vector<double>  zero_rates = { 0.05, 0.05, 0.05,0.05, 0.05, 0.05, 0.05,0.05, 0.05, 0.05, 0.05,0.05 };
	std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10,30 };

	for (auto& dt : dates)
	{
		dt *= 365;
	//	dt+=daysOffset;
	}

	auto my_curve = makeCurve(dates, zero_rates);


	for (int j = 0; j < 20; j++)
	{
		auto df = my_curve.getDiscountFactor(j * 365.);// +daysOffset);
		std::cout << "years =" << j << "days =" << j * 365 << "expected df" << exp(-0.05 * j) << "actual df =" << df << std::endl;
	}

	std::cout << "half years \n";
	for (int j = 0; j < 20; j++)
	{
		auto df = my_curve.getDiscountFactor(j * 365. + 182.);// +daysOffset); // half year

		auto df_full_yr = my_curve.getDiscountFactor(j * 365.);// +daysOffset); // full year

		auto df_halfYear = exp(-0.05 * 182 / 365.);
		std::cout << "years =" << j << "days =" << j * 365 +180 << " expected df" << df_halfYear* df_full_yr << " actual df =" << df << std::endl;
	}

	/////////////////////////////////////////////// new curve /////////////////////////////
	std::cout << "curve sloping up \n";

	std::vector<double>  zero_rates_2 = { 0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.11,0.12, 0.13, 0.14, 0.15,0.16 };
	//std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10,30 };
	auto my_curve2 = makeCurve(dates, zero_rates_2);

	for (int j = 0; j < 10; j++)
	{
		auto df = my_curve2.getDiscountFactor(j * 365.);
		auto rate = (0.05 + 0.01 * j);
		std::cout << "years =" << j << "days =" << j * 365 << "expected df" << exp(-(rate* j))<< "actual df =" << df << std::endl;
	}


	//make some cashflows
	/*
	auto  flows = makeCashFlowsAOS(0, 182.5, 20, 10.0);

	double total_price = 0.;

	for (const auto& flow : flows)
	{
		total_price += flow.amount * my_curve2.getDiscountFactor(flow.date);
	}

	*/
	//CashFlowsSOA tenYearAnnual = makeCashFlowsSOA(0, 182.5, 20, 10);

	//double price = 0.;

	
	double total_price = 0.0;
	auto cashFlows = makeCashFlowsAOS(2000, 1, 1, 10, 0.05, 100.0);
	for (const auto& flow : cashFlows)
	{
		total_price += flow.amount * my_curve2.getDiscountFactor(flow.date);
	}


	std::cout << "price =" << total_price << std::endl;

	


}

/**/

OldCurve makeBigCurve()
{

	std::vector<double>  zero_rates = {  0.05, 0.05,0.05, 0.05, 0.05, 0.05,0.05, 0.05, 0.05, 0.05,0.05 };
	std::vector <long>   dates = { 1,2,3,4,5,6,7,8,9,10,30 };

	for (auto& dt : dates)
	{
		dt *= 365;
	}

	std::vector <long> monthly = { 0, 30, 60,90,120,150,180,210,240,270,300,330 };
	std::vector <double> monthly_zero_rates = { 0.02, 0.03, 0.04, 0.04,0.05,0.05,0.05, 0.055,0.056,0.057,0.058, 0.059 };

	for (const auto& dt : dates)
	{
		monthly.emplace_back(dt);
	}

	for (const auto& rate : zero_rates)
	{
		monthly_zero_rates.emplace_back(rate);
	}

	return makeCurve(monthly, monthly_zero_rates);

}

void do_basic_curve_1()
{
	OldCurve testCrv = makeBigCurve();
	auto  flows = makeCashFlowsAOS(2000, 1, 1, 50, 0.05, 100.0);
		//makeCashFlowsAOS(0, 182.5, 50, 1.0); // 50 semi annual

	using BOND = std::vector<CashFlow>;
	long numBonds = 1000000;

	std::vector<double> prc_res(numBonds);
	
	std::vector<BOND>  PF;
	for (long l = 0; l < numBonds; l++)
	{
		PF.emplace_back(flows);
	}

	volatile double price = 0.0;
	int KRVSZ = 1;

	std::cout << "start" << std::endl;
	auto startTme = std::chrono::high_resolution_clock::now();

	int i = 0;
	for( const auto& bnd: PF)
	{
		price = 0.0;
		for (const auto& cf : bnd)
		{
			price += cf.amount*    testCrv.getDiscountFactor(cf.date);
		}
		prc_res[i] = price; //store result
		++i;
	}

	std::cout << "done" << std::endl;

	auto endTime = std::chrono::high_resolution_clock::now();
	auto runtime = endTime - startTme;
	std::cout << "run time " << numBonds << " bonds of " << 50 << " casflows and " << KRVSZ << " buckets = " << runtime.count() / 1000000000.0 << std::endl;


}





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
	try
	{

	std::vector<double>  values = { 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.};
	std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10,50 };
	std::vector <double>   datesD = { 0,1,2,3,4,5,6,7,8,9,10 };

	int KRVSZ = 200;// 00;
	std::vector< VecXX>  vecVals;
	for (int i = 0; i < 12; i++)
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

	auto makeYears = [](auto x) { return x * 365; };
	std::transform(begin(dates), end(dates), begin(dates), makeYears);

	Curve2<long, VecXX, ZeroInterpCached<long, VecXX>>  testCurve2(KRVSZ);
	testCurve2.setValues(begin(dates), end(dates), begin(vecVals), end(vecVals));

	auto valV = testCurve2.valueAt(0.0);

	auto valV2 = testCurve2.valueAt(5);

	//using cashFlow = std::pair< long, double>;

	//using swap = CashFlow[40];

	using swap = std::vector<CashFlow>;

	// using  instruments = std::vector<swap>;

	using PF = std::vector<swap>;

	//auto bigPF = PF(1000000);


	//auto bigPF = PF(999999);// 00000);


	PF bigPF(999999);
	// double stepSz = 0.25;

	//PF bigPF;


	double count = 0;
	for (auto& swp : bigPF)
	{
		swp.reserve(40);

		swp = swap(40);

		long date = 0;
		double period = 0.25 * 365;

		//double offset = int(count / 30000.0) * 0.03 * 365.0;
		long offset = int(count / (999999 / 30));// *30;
		for (auto& cf : swp)
		{
			cf.date = date + offset;
			date += period;
			cf.amount = 0.1;
		}
		count++;
	}
	


	auto results = new double[1000000 * KRVSZ];

	auto startTme = std::chrono::high_resolution_clock::now();


	auto dfsFixed = testCurve2.valueAt(60);

	int j = 0;
	for (auto& swp : bigPF)
	{

		VecXX prices(0.0, KRVSZ);
		SpanXX prc(prices.begin(), KRVSZ);

		auto cf_t = (bigPF[0])[1];

		const VecXX& dfs = testCurve2.valueAt(cf_t.date);

		for (auto& cf : swp)
		{
			// auto DF = testCurve2.valueAt(cf.first);
			// prices += cf.second * DF;// testCurve2.valueAt(cf.first);// DF * cf.second + prices;



			//VecXX dfs = testCurve2.valueAt(cf.first);
			//const  VecXX::SCALA_TYPE  cfs = cf.second;
			//prices = FMA(dfs, cfs, prices);





		//	std::cout << ", " << cf.first;

//			const VecXX& dfs = testCurve2.valueAt(cf.first);


			SpanXX spnDf(dfs.begin(), dfs.size());

			VecXX::INS  cashflow = cf.amount;

			auto priceAndAccumulateCashFlow = [cashflow](auto price, auto  df)
			{
				return mul_add(cashflow, df, price);
			};

			//prices =
			transformM(priceAndAccumulateCashFlow, prices, dfs);
			//	transform(priceAndAccumulateCashFlow, prc, spnDf, prc);



		}

	//	std::cout << "\n";
	//	std::cout << "\n";


		auto basePrice = prices[0];
		auto risk = prices - basePrice;
		//allPrices[j] =prices;
   //	 count++;

		

		for (int kk = 0; kk < KRVSZ; kk++)
		{
			results[j * KRVSZ + kk] = risk[kk];
		}
		j++;

 

	}

	std::cout << "run \n";

	auto endTime = std::chrono::high_resolution_clock::now();
	auto runtime = endTime - startTme;
	std::cout << "run time 1M swaps " << KRVSZ << " buckets = " << runtime.count() / 1000000000.0 << std::endl;

}
catch (std::exception& expt)
{
	std::cout << "exception thrown" << expt.what() << "\n";
}
}




void doCurve2A()
{
	try
	{

		std::vector<double>  values = { 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11. };
		std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10,50 };
		std::vector <double>   datesD = { 0,1,2,3,4,5,6,7,8,9,10 };

		int KRVSZ = 200;// 00;
		std::vector< VecXX>  vecVals;
		for (int i = 0; i < 12; i++)
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

		auto makeYears = [](auto x) { return x * 365; };
		std::transform(begin(dates), end(dates), begin(dates), makeYears);


		//Curve2<double, VecXX, ZeroInterpCached<double, VecXX>>  testCurve2(KRVSZ);
		Curve2<long, VecXX, ZeroInterpCached<long, VecXX>>  testCurve2(KRVSZ);
		//	testCurve2.setValues(begin(datesD), end(datesD), begin(vecVals), end(vecVals));
		testCurve2.setValues(begin(dates), end(dates), begin(vecVals), end(vecVals));

		auto valV = testCurve2.valueAt(0.0);

		//auto valV2 = testCurve2.valueAt(0.5);

		
		using cashFlow = std::pair< long, double>;

		//using swap = cashFlow[40];
		//using swap = cashFlow[40] ;

		using swap = std::vector<cashFlow>;

		
		using PF = std::vector<swap>;

		PF bigPF;// = PF;// (1000000);
		//auto bigPF = PF(999999);// 00000);
		// double stepSz = 0.25;
		//PF bigPF;

		
		double count = 0;
		for (auto& swp : bigPF)
		{
			long date = 0.0;
			double period = 0.25 * 365;

			//double offset = int(count / 30000.0) * 0.03 * 365.0;
			long offset = int(count / (999999 / 365));
			for (auto& cf : swp)
			{
				cf.first = date + offset;
				date += period;
				cf.second = 0.1;
			}
			count++;
		}



		auto results = new double[1000000 * KRVSZ];

		auto startTme = std::chrono::high_resolution_clock::now();


		//auto results=  new double[100000*200];

		//std::vector<VecXX> allPrices[1000000];

		// VecXX prices(0.0, 200);


		auto dfsFixed = testCurve2.valueAt(60);

		int j = 0;

		std::vector<  const VecXX*> last_dfs( 1000, nullptr);
		std::vector< long>  last_dates( 1000, -1);


		for (auto& swp : bigPF)
		{

			VecXX prices(0.0, KRVSZ);
			SpanXX prc(prices.begin(), KRVSZ);
	

			int pos = -1;

			for (auto& cf : swp)
			{
	
				pos++;
				const VecXX* pdfs = nullptr;
				if (last_dates[pos] == cf.first)
				{
					pdfs = last_dfs[pos];
				}
				else
				{
					//local check
					if (last_dates[pos + 1] == cf.first)
					{
						pdfs = last_dfs[pos+1];
					}
					else
					{
						auto pdfs_lu = testCurve2.valueAt(cf.first);

						last_dates.insert(last_dates.begin() + pos, cf.first);
						last_dfs.insert(last_dfs.begin() + pos, &pdfs_lu);
						pdfs = &pdfs_lu;
					}
				}
			

				const VecXX& dfs = *pdfs;
				SpanXX spnDf(dfs.begin(), dfs.size());

				VecXX::INS  cashflow = cf.second;

				auto priceAndAccumulateCashFlow = [cashflow](auto price, auto  df)
				{
					return mul_add(cashflow, df, price);
				};


				transformM(priceAndAccumulateCashFlow, prices, dfs);

			}

			auto basePrice = prices[0];
			auto risk = prices - basePrice;
			//allPrices[j] =prices;
	   //	 count++;

			

			for (int kk = 0; kk < KRVSZ; kk++)
			{
				results[j * KRVSZ + kk] = risk[kk];
			}
			j++;

	


		}

		std::cout << "run \n";

		auto endTime = std::chrono::high_resolution_clock::now();
		auto runtime = endTime - startTme;
		std::cout << "run time 1M swaps " << KRVSZ << " buckets = " << runtime.count() / 1000000000.0 << std::endl;

	}
	catch (std::exception& expt)
	{
		std::cout << "exception thrown" << expt.what() << "\n";
	}
}



/*

template <typename numType, typename T1, typename T2>
struct MYFUNC : public T1, public T2
{
	numType F(const numType& val) { return T1(val); }
	numType DF(const numType& val) { return T2(val); }
};



void do_curve3()
{
	std::vector<double>  values = { 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 ,30.0};
	std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10,30 };
	std::vector <double>   datesD = { 0,1,2,3,4,5,6,7,8,9,10,30 };

	std::vector< VecXX>  vecVals;
	for (int i = 0; i < 12; i++)
	{
		VecXX vals(i * 0.001 + 0.06, 200);
		vecVals.push_back(vals);

	}


	int KRVSZ = 200;

	//Curve<double, VecXX> forecastCurve; // just linear interpolation

	Curve2<double, VecXX, ZeroInterpCached<double, VecXX>>  discountCurve(KRVSZ);
	discountCurve.setValues(begin(datesD), end(datesD), begin(vecVals), end(vecVals));

	auto valV = discountCurve.valueAt(0.0);

	auto valV2 = discountCurve.valueAt(0.5);

	std::vector<double> cashFlows(200, 0.5);

	auto startTme = std::chrono::high_resolution_clock::now();

	double offsetDay = 0.0;

	for (long ll = 0; ll < 30; ll++) // so month days 
	{
		VecXX prices(0.0, 200);

		for (long l = 0; l < 30000; l++) //do all on the same month
		{
			int jj = 0;
			for (double d = 0.1; d < 10; d += 0.25, ++jj)
			{

				const auto& DF = discountCurve.valueAt(d + offsetDay);
				//prices += DF * cashFlows[jj];

				//FMA(DF, cashFlows, prices);

/////////////////////////////////////////////////////////////////
				//const VecXX& dfs = testCurve2.valueAt(cf.first);


				//SpanXX spnDf(DF.begin(), DF.size());

				VecXX::INS  cashflow = cashFlows[jj];// cf.second;

				auto priceAndAccumulateCashFlow = [cashflow](auto price, auto  df)
				{
					return mul_add(cashflow, df, price);
				};

				//prices =
			//	transformM(priceAndAccumulateCashFlow, prices, dfs);
				transformM(priceAndAccumulateCashFlow, prices, DF);// , prices);
			}
			//auto valV3 = testCurve2.valueAt(0.5);
		}

		offsetDay += 1./365;
	} //ll

	auto endTime = std::chrono::high_resolution_clock::now();
	auto runtime = endTime - startTme;
	std::cout << "run time 900k swaps of 40 casflows and  " << KRVSZ << " buckets = " << runtime.count() / 1000000000.0 << std::endl;


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







void doCurve2A_shuffled()
{
	try
	{

		std::vector<double>  values = { 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11. };
		std::vector <long>   dates = { 0,1,2,3,4,5,6,7,8,9,10,50 };
		std::vector <double>   datesD = { 0,1,2,3,4,5,6,7,8,9,10 };

		int KRVSZ = 200;// 00;
		std::vector< VecXX>  vecVals;
		for (int i = 0; i < 12; i++)
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

		auto makeYears = [](auto x) { return x * 365; };
		std::transform(begin(dates), end(dates), begin(dates), makeYears);


		//Curve2<double, VecXX, ZeroInterpCached<double, VecXX>>  testCurve2(KRVSZ);
		Curve2<long, VecXX, ZeroInterpCached<long, VecXX>>  testCurve2(KRVSZ);
		//	testCurve2.setValues(begin(datesD), end(datesD), begin(vecVals), end(vecVals));
		testCurve2.setValues(begin(dates), end(dates), begin(vecVals), end(vecVals));

		auto valV = testCurve2.valueAt(0.0);

		//auto valV2 = testCurve2.valueAt(0.5);


		using cashFlow = std::pair< long, double>;

		using swap = cashFlow[40];


		using PF = std::vector<swap>;

		//auto bigPF = PF(1000000);
	//	auto bigPF = PF(999999);// 00000);

		long BIG_PF_SIZE = 999999;

		BIG_PF_SIZE = 20999;

		//auto bigPF = PF(20999);// 00000);

		//auto bigPF = PF(BIG_PF_SIZE);// 00000);

		PF bigPF;
		// double stepSz = 0.25;

		
		double count = 0;
		for (auto& swp : bigPF)
		{
			long date = 0.0;
			double period = 0.25 * 365;

			//double offset = int(count / 30000.0) * 0.03 * 365.0;
			long offset = int(count / (BIG_PF_SIZE / 365));
			for (auto& cf : swp)
			{
				cf.first = date + offset;
				date += period;
				cf.second = 0.1;
			}
			count++;
		}


		std::cout << "starting shuffle" << std::endl;

		std::random_device rd;
		std::mt19937 g(rd());
		shuffle(bigPF.begin(), bigPF.end(), g);
		//return shuffled_path;


		std::cout << "end shuffle" << std::endl;

		auto results = new double[1000000 * KRVSZ];


		std::cout << "starting pricing" << std::endl;

		auto startTme = std::chrono::high_resolution_clock::now();


		//auto results=  new double[100000*200];

		//std::vector<VecXX> allPrices[1000000];

		// VecXX prices(0.0, 200);


		auto dfsFixed = testCurve2.valueAt(60);

		int j = 0;

		std::vector<  const VecXX*> last_dfs(1000, nullptr);
		std::vector< long>  last_dates(1000, -1);


		for (auto& swp : bigPF)
		{

			VecXX prices(0.0, KRVSZ);
			SpanXX prc(prices.begin(), KRVSZ);


			int pos = -1;

			for (auto& cf : swp)
			{

				pos++;
				const VecXX* pdfs = nullptr;
				if (last_dates[pos] == cf.first)
				{
					pdfs = last_dfs[pos];
				}
				else
				{
					//local check
					if (last_dates[pos + 1] == cf.first)
					{
						pdfs = last_dfs[pos + 1];
					}
					else
					{
						auto pdfs_lu = testCurve2.valueAt(cf.first);

						last_dates.insert(last_dates.begin() + pos, cf.first);
						last_dfs.insert(last_dfs.begin() + pos, &pdfs_lu);
						pdfs = &pdfs_lu;
					}
				}


				const VecXX& dfs = *pdfs;
				SpanXX spnDf(dfs.begin(), dfs.size());

				VecXX::INS  cashflow = cf.second;

				auto priceAndAccumulateCashFlow = [cashflow](auto price, auto  df)
					{
						return mul_add(cashflow, df, price);
					};


				transformM(priceAndAccumulateCashFlow, prices, dfs);

			}

			auto basePrice = prices[0];
			auto risk = prices - basePrice;
			//allPrices[j] =prices;
	   //	 count++;

			

			for (int kk = 0; kk < KRVSZ; kk++)
			{
				results[j * KRVSZ + kk] = risk[kk];
			}
			j++;




		}

		std::cout << "run \n";

		auto endTime = std::chrono::high_resolution_clock::now();
		auto runtime = endTime - startTme;
		std::cout << "run time 1M swaps " << KRVSZ << " buckets = " << runtime.count() / 1000000000.0 << std::endl;

	}
	catch (std::exception& expt)
	{
		std::cout << "exception thrown" << expt.what() << "\n";
	}
}

*/



int main()
{

//	example_cashflows();

///	do_basic_curve();


//	do_basic_curve_1();


 //  std::cout << "Hello World!\n";

	

//	for (int i = 0; i < 15; i++)
//	{

		doCurve2();

//		doCurve2A();

	//	doCurve2A_shuffled();
//	}

	//do_curve3();
	return 0;
}

