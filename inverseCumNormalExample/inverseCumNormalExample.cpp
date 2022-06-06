// accumulateExample.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <algorithm>
#include <random>
#include <numeric>
#include <iterator>
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>  
#include <immintrin.h>
#include <chrono>
#include <iostream>
#include <functional>



#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/apply_operation.h"
#include "../Vectorisation/VecX/vec_d.h"
#include "../Vectorisation/VecX/vec_bool.h"
#include "../Vectorisation/VecX/vec_view.h"

#include "../Vectorisation/VecX/target_name_space.h"



#include "cdfNormalInverse.h"

//set instruction set via namespace in cdfNormalInverse.h
void compareNumerics();
void doPerformanceComparison();
void profileExample();



int main()
{
	doPerformanceComparison();
	//	compareNumerics();
	//	profileRun();
	return 0;
}


void doPerformanceComparison()
{

	long loop = 10000;
	for (int SZ = 100; SZ < 10000; SZ += 100)
	{

		using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
		std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));


		for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }

		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(&v[0], &v[SZ - 1], g);


		VecXX test(v);
		auto test2 = test;

		{
			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto 	res = calcCDFNormsSparseFMA(test);
			}
			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << SZ << " calcCDFNormsSparseFMA ," << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << " , evluations run , " << SZ * loop * 1000000000.0 / (runtime.count()) << ",  evals per sec |,";// << std::endl;
		}



		{
			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto res = calcCDFNormWithViews(test);
			}
			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << " calcCDFNormWithViews ," << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << ", secs , evluations run , " << SZ * loop * 1000000000.0 / (runtime.count()) << ", evals per sec |,";// << std::endl;
		}



		{
			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto res = calcCDFNormsSparseFMAOnePass(test);
			}
			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << " calcCDFNormsSparseFMAOnePass ," << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << ",secs , evluations run , " << SZ * loop * 1000000000.0 / (runtime.count()) << ", evals per sec |, ";// << std::endl;
		}


	
		{
			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto res =  calcCDFNormWichuraViewsAndFMA2splits(test);
			}
			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << " calcCDFNormWichuraViewsAndFMA time ," << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << ", secs , evluations run , " << SZ * loop * 1000000000.0 / (runtime.count()) << ",  evals per sec |,";// << std::endl;

		}

		{
			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto 	res = calcCDFNormWithViewsAndFMA(test);
			}

			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << SZ << " calcCDFNormWithViewsAndFMA time aclam 1 bill = " << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << ",secs , evluations run , " << SZ * loop * 1000000000.0 / (runtime.count()) << ",  evals per sec | ,";// << std::endl;
		}


		{
			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto 	res = calcCDFNormWithViewsAndFMAWriteFromViewCalc(test);
			}
			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << SZ << " calcCDFNormWithViewsAndFMAandViewCalcWrite time aclam 1 bill = " << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << ",secs , evluations run , " << SZ * loop * 1000000000.0 / (runtime.count()) << ",  evals per sec | ,";// << std::endl;
		}

		{

			std::vector<double> res(test.size(), .0);
			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto cdfNormLambda = [](auto x) { return static_cast<double>(qnorm8(x) ); };
				auto itRes = res.begin();
				std::transform(test.begin(), test.end(), itRes, cdfNormLambda);
			}

			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << " transform  with lambda of a double,  time aclam 1 bill = " << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << ",secs , evluations run , " << SZ * loop * 1000000000.0 / (runtime.count()) << "  evals per sec |\n \n" << std::endl;
		}


		/*
		// intel svnl called as lambda over vector, set vcl to lambda calling svml function  for target instruction set 
		{
			//SVML functions
			//  auto vcl = [&](__m512d x) { return _mm512_cdfnorminv_pd(x);  };
			//  auto vcl = [](Vec8d& x) { return _mm512_cdfnorminv_pd(x);  };
			//  auto vcl = [&](Vec4d& x) { return _mm256_cdfnorminv_pd(x); };  
			//	auto vcl = [](Vec4d& x) { return _mm256_cdfnorminv_pd(x); };     
			//	auto vcl = [](Vec16f& x) { return _mm512_cdfnorminv_ps(x);  };
			//  auto vcl = [](Vec8f& x) { return _mm256_cdfnorminv_ps(x);  };

				auto startTme = std::chrono::high_resolution_clock::now();
				for (long l = 0; l < loop; l++)
				{
					auto res = ApplyUnitaryOperation(test, vcl);
				}
				auto endTime = std::chrono::high_resolution_clock::now();
				auto runtime = endTime - startTme;
				std::cout << " svml raw  time aclam ," << runtime.count() / 1000000000.0  ;// << std::endl;
				std::cout << ", secs , evluations run , " << int(SZ * loop) << ", " << (SZ * loop) * 1000000000.0 / (runtime.count()) << " |,";// evals per sec" << std::endl;
		}

		*/


	}

}


//compare different implementations
void compareNumerics()
{

	for (float x = 0.001f; x < 1.0f; x += 0.003f)
	{
		VecXX X(x, 16);
		//auto res =calcCDFNormsSparseFMAOnePass(X)[0]; //OK
		//auto res = calcCDFNormsSparseFMA(X)[0]; // good
		//auto res = calcCDFNormWichuraViewsAndFMA(X)[0]; // superb accuract
		//auto res = calcCDFNormWithViewsAndFMA(X)[0]; // superb accuract
		//auto res = calcCDFNormWichuraViewsAndFMA(X)[0]; // superb accuract
		auto res =  calcCDFNormWichuraViewsAndFMA2splits(X)[0]; // superb accuract
		std::cout << x << "," << res << "expected" << qnorm8(x) << ",err" << res - qnorm8(x) << std::endl;
	}

}


void profileExample()
{
	//	int SZ = 2000;
	long loop = 100000;// 100;
//	for (int SZ = 100; SZ < 10000; SZ += 100)
	int SZ = 10000;
	{

		using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
		std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));


		for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }

		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(&v[0], &v[SZ - 1], g);

		VecXX test(v);
		auto test2 = test;


		{
			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				//auto 	res = calcCDFNormsSparseFMA(test);
				  auto 	res = calcCDFNormsSparseFMA2(test); //uses sparse update proper
				//auto 	res = calcCDFNormWichuraViewsAndFMA(test);
			    //auto 	res = calcCDFNormWichuraViewsAndFMA2splits(test);
				//auto 	res =  calcCDFNormBranchlessFloat(test);
			}

			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << " , evluations run , " << SZ * loop * 1000000000.0 / (runtime.count()) << ",  evals per sec |,";// << std::endl;

		}
	}
}
