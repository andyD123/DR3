// accumulateExample.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "../Vectorisation/VecX/vec.h"
#include "../Vectorisation/VecX/operations.h"
#include "../Vectorisation/VecX/apply_operation.h"
#include "../Vectorisation/VecX/vec_d.h"
#include "../Vectorisation/VecX/vec_bool.h"
#include "../Vectorisation/VecX/vec_view.h"
#include "../Vectorisation/VecX/target_name_space.h"

#include "cdfNormalInverse.h"


#include <algorithm>
#include <random>
#include <numeric>
#include <iterator>
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>  
#include <immintrin.h>



auto numOps = [](int TEST_LOOP_SZ, int SZ) { return  static_cast<int>(double(TEST_LOOP_SZ) * double(SZ)); };



void compareNumerics();
void doPerformanceComparison();
void profileExample();



template<typename T>
bool vectorsEqual(const std::vector<T>& C1, const std::vector<T>& C2, const std::vector<T>& C3)
{
	bool  testOK = true;
	const double ERR = 1e-15; //for examples
	if (C1.size() != C2.size())
	{
		return false;
	}

	if (C3.size() != C2.size())
	{
		return false;
	}

	for (size_t i = 0; i < C3.size(); i++)
	{
		auto err1 = fabs((C1[i] - C2[i]) / (C2[i] + C1[i]));
		auto err2 = fabs((C1[i] - C3[i]) / (C1[i] + C3[i]));

		if ((err1 > ERR) || (err2 > ERR))
		{
			testOK = false;
			std::cout << "\n err diff@ " << i << " err1 =" << err1 << ", err2 = " << err2 << "\n";
			break;
		}
	}

	return testOK;

}

//set instruction set via namespace in cdfNormalInverse.h
int main()
{
      doPerformanceComparison();
	//	compareNumerics();
	//profileExample();
	return 0;
}


void doPerformanceComparison()
{

	long loop = 10000;
	//for (int SZ = 100; SZ < 10000; SZ += 100)
	//{

		using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;
	//	std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));


		//for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }

		std::random_device rd;
		std::mt19937 g(rd());
	//	std::shuffle(&v[0], &v[SZ - 1], g);


		//VecXX test(v);
		//auto test2 = test;

		//VecXX resultsWichura;
		//std::vector<double> resT(test.size(), .0);

		//
		///*
		for (int SZ = 100; SZ < 10000; SZ += 100)
		{
			std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
			for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }
			std::shuffle(&v[0], &v[SZ - 1], g);


			VecXX test(v);
			auto test2 = test;

			//warm
			for (long l = 0; l < loop; l++)
			{
				auto 	res = calcCDFNormsSparseFMA(test);
			}

			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto 	res = calcCDFNormsSparseFMA(test);
			}
			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << SZ  <<  "," << " calcCDFNormsSparseFMA ," << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << " , evluations run , " << numOps(loop,SZ) * 1000000000.0 / (runtime.count()) << ",  evals per sec |, \n";// << std::endl;
		}

		for (int SZ = 100; SZ < 10000; SZ += 100)
		{
			std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
			for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }
			std::shuffle(&v[0], &v[SZ - 1], g);


			VecXX test(v);
			auto test2 = test;
		
			for (long l = 0; l < loop; l++)
			{
				auto res = calcCDFNormsSparseFMAOnePass(test);
			}
			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto res = calcCDFNormsSparseFMAOnePass(test);
			}
			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << SZ << ", calcCDFNormsSparseFMAOnePass ," << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << ",secs , evluations run , " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << ", evals per sec |, \n";// << std::endl;
		}


		for (int SZ = 100; SZ < 10000; SZ += 100)
		{
			std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
			for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }
			std::shuffle(&v[0], &v[SZ - 1], g);


			VecXX test(v);
			auto test2 = test;
			//warm
			for (long l = 0; l < loop; l++)
			{
				auto res = calcCDFNormWithViews(test);
			}

			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto res = calcCDFNormWithViews(test);
			}
			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << SZ << ", calcCDFNormWithViews, " << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << ", secs , evluations run , " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << ", evals per sec |, \n";// << std::endl;
		}



		for (int SZ = 100; SZ < 10000; SZ += 100)
		{
			std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
			for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }
			std::shuffle(&v[0], &v[SZ - 1], g);


			VecXX test(v);
			auto test2 = test;
			for (long l = 0; l < loop; l++)
			{
				auto 	res = calcCDFNormWithViewsAndFMA(test);
			}

			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto 	res = calcCDFNormWithViewsAndFMA(test);
			}

			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << SZ << ",  calcCDFNormWithViewsAndFMA time aclam 1 bill = " << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << ",secs , evluations run , " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << ",  evals per sec | ,\n";// << std::endl;
		}

		

		for (int SZ = 100; SZ < 10000; SZ += 100)
		{
			std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
			for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }
			std::shuffle(&v[0], &v[SZ - 1], g);


			VecXX test(v);
			auto test2 = test;
			for (long l = 0; l < loop; l++)
			{
				auto 	res = calcCDFNormWithViewsAndFMAWriteFromViewCalc(test);
			}

			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto 	res = calcCDFNormWithViewsAndFMAWriteFromViewCalc(test);
			}
			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << SZ << ",  calcCDFNormWithViewsAndFMAandViewCalcWrite time aclam 1 bill = " << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << ",secs , evluations run , " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << ",  evals per sec | ,\n";// << std::endl;
		}
		
		for (int SZ = 100; SZ < 10000; SZ += 100)
		{
			std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
			for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }
			std::shuffle(&v[0], &v[SZ - 1], g);


			VecXX test(v);
			auto test2 = test;

			VecXX resultsWichura;
			std::vector<double> resT(test.size(), .0);


			auto startTme = std::chrono::high_resolution_clock::now();
			for (long l = 0; l < loop; l++)
			{
				auto cdfNormLambda = [](auto x) { return static_cast<double>(qnorm8(x) ); };
				auto itRes = resT.begin();
				std::transform(test.begin(), test.end(), itRes, cdfNormLambda);
			}

			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << SZ << ",  transform  with lambda of a double,  time aclam 1 bill = " << runtime.count() / 1000000000.0;//<< std::endl;
			std::cout << ",secs , evluations run , " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << "  evals per sec |\n";
		}

		
//*/

		for (int SZ = 100; SZ < 10000; SZ += 100)
		{
			std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
			for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }
			std::shuffle(&v[0], &v[SZ - 1], g);


			VecXX test(v);
			auto test2 = test;
			VecXX resultsWichura;
			std::vector<double> resT(test.size(), .0);

			{

				for (long l = 0; l < loop; l++)
				{
					//auto res
					resultsWichura = calcCDFNormWichuraViewsAndFMA2splits(test);
				}

				auto startTme = std::chrono::high_resolution_clock::now();
				for (long l = 0; l < loop; l++)
				{
					//auto res
					resultsWichura = calcCDFNormWichuraViewsAndFMA2splits(test);
				}
				auto endTime = std::chrono::high_resolution_clock::now();
				auto runtime = endTime - startTme;
				std::cout << SZ << ",  calcCDFNormWichuraViewsAndFMA time ," << runtime.count() / 1000000000.0;//<< std::endl;
				std::cout << ", secs , evluations run , " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << ",  evals per sec |,\n";// << std::endl;

			}
		}
		//
		/*
		

		//*/
//
/*
//

		for (int SZ = 100; SZ < 10000; SZ += 100)
		{
			std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
			for (int i = 0; i < SZ; i++) { v[i] = i / (FloatType)(SZ); }
			std::shuffle(&v[0], &v[SZ - 1], g);
			VecXX test(v);


			// intel svml called as lambda over vector, set vcl to lambda calling svml function  for target instruction set 
			{
				//SVML functions
				//  auto vcl = [&](__m512d x) { return _mm512_cdfnorminv_pd(x);  };
				auto vcl = [](Vec8d& x) { return _mm512_cdfnorminv_pd(x);  };
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
				std::cout << " svml raw  time aclam ," << runtime.count() / 1000000000.0;// << std::endl;
				std::cout << ", secs , evluations run , " << numOps(loop, SZ) << ", " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << " |, evals per sec" << std::endl;
			}

		}
//
*/
		/*	*/
		/*
		std::vector<VecXX::SCALA_TYPE> wch = resultsWichura;
		
		if (vectorsEqual(resT, resT, wch))
		{
			std::cout << "\n transform and WS241 MATCH \n";
		}
		else
		{
			std::cout << "\n FAIL transform and WS241 DO NOT MATCH \n";
		}
		std::cout << " | \n \n" << std::endl;
	//}

		*/

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
				//  auto 	res = calcCDFNormsSparseFMA2(test); //uses sparse update proper
				//auto 	res = calcCDFNormWichuraViewsAndFMA(test);
			    auto 	res = calcCDFNormWichuraViewsAndFMA2splits(test);
				//auto 	res =  calcCDFNormBranchlessFloat(test);
			}

			auto endTime = std::chrono::high_resolution_clock::now();
			auto runtime = endTime - startTme;
			std::cout << " , evluations run , " << numOps(loop, SZ) * 1000000000.0 / (runtime.count()) << ",  evals per sec |,";// << std::endl;

		}
	}
}
