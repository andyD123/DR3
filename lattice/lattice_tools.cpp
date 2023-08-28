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


#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/accumulate_transform.h"
#include "../Vectorisation/VecX/error_utils.h"


#include "../Vectorisation/VecX/zip_utils.h"
#include "../Vectorisation/VecX/span.h"

#include "utils.h"
#include "lattice_tools.h"


/*
Experimental features for using SIMD in rangified manner 
lagged/sampled valued being introduced into iterating lambda

zip iterator

strided iteration

span based access 

and abridged version of layout scheme  MDSpan .. not the standard

*/

static auto getRandomShuffledVector(int SZ, int instance_number = 0)
{
	using FloatType = typename InstructionTraits<VecXX::INS>::FloatType;


	static std::map<int, std::vector<FloatType> > vectors;


	int key = 10 * SZ + instance_number;
	//store vectors with key 10 times size  and add on 0-9 integer for instance of different random vector

	if (SZ < 0)
	{
		vectors.clear();
		SZ = 0;
	}


	if (vectors.find(key) != vectors.end())
	{
		return vectors[key];
	}
	else
	{
		std::vector<FloatType>  v(SZ, VecXX::SCALA_TYPE(6.66));
		for (int i = 0; i < SZ; i++) { v[i] += i; }
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(begin(v), end(v), g);
		vectors[key] = v;
		return v;
	}
}



void doScan()
{

	auto v1 = std::vector<VecXX::SCALA_TYPE>(4000, 0.0);
	FLOAT last = 0.0;
	for (auto& x : v1)
	{
		x = last + 1.0;
		last = x;
	}
	
	VecXX a(v1);
	std::vector<FLOAT> dbg = v1;

	auto sum = [](auto X, auto Y) { return  X + Y; };

	double time = 0;
	{
		TimerGuard timer(time);

		VecXX result;
		for (long l = 0; l < 1000000; ++l)
		{
			result = ApplyScan(a, sum);
			//dbg = a;
		}

		dbg = result;
	}
	std::cout << "run time = " << time << "\n";



	auto dbg_cpy = dbg;
	time = 0;
	{
		TimerGuard timer(time);

		for (long l = 0; l < 1000000; ++l)
		{
			std::inclusive_scan(begin(v1), end(v1), dbg.begin());
		}
	}
	std::cout << "run time std::scan = " << time << "\n";

	for (size_t i = 0; i < dbg_cpy.size(); ++i)
	{
		std::cout << i << "  ,  " << dbg[i] << "   , " << dbg_cpy[i] << std::endl;
	}

}

void testSampler()
{

	FLOAT sampleData[] = { 0., 1.0,2.,3.,4.,5.,6.,7.,8.,9.,10. };

	TrinomialSampler<VecXX::INS>  zzz;

	zzz.load(&sampleData[1]);

	auto x = zzz.X_Minus_1.value;
	auto y = zzz.X_0.value;
	auto z = zzz.X_1.value;

	ignore(x);
	ignore(y);
	ignore(z);
}

void doStridedSpan()
{
	using FLOAT = VecXX::SCALA_TYPE;
	auto v1 = std::vector<VecXX::SCALA_TYPE>(4000, 0.0);
	FLOAT last = 0.0;
	for (auto& x : v1)
	{
		x = last + 1.0;
		last = x;
	}

	VecXX a_vec(v1);
	std::vector<FLOAT> dbg = v1;
	VecXX result = a_vec;
	result *= 0.0;

	StridedSampler<VecXX::INS> strided_sampler(2);
	SpanXX spn(a_vec.start(), 800);

	SpanXX spnPlus(result.start(), 800);

	auto SQR = [](StridedSampler<VecXX::INS>& sampler)
	{
		auto x = sampler.X_0.value;
		return x * x;
	};


	ApplyTransformUR_X_Impl_EX(spn, spnPlus, SQR, strided_sampler, 0, int(spn.paddedSize()));

	std::vector<FLOAT> ddbbgg = spn;
	ddbbgg = spnPlus;

	VecXX rootvec = result;

	auto SQRT2 = [](auto& x)
	{
		return sqrt(x);
	};

	auto rooted = spnPlus;
	transform(SQRT2, spnPlus, rooted);

	ddbbgg = rooted;

}

void doTransformWithASpan()
{

	auto v1 = std::vector<VecXX::SCALA_TYPE>(4000, 0.0);
	FLOAT last = 0.0;
	for (auto& x : v1)
	{
		x = last + 1.0;
		last = x;
	}

	VecXX a_vec(v1);
	std::vector<FLOAT> dbg = v1;

	VecXX result = a_vec;

	SpanXX spn(a_vec.start(), 200);

	SpanXX spnPlus(a_vec.start() + 200, 200);

	UnitarySampler<VecXX::INS> identity_sampler;


	auto SQR = [](UnitarySampler<VecXX::INS>& sampler)
	{
		auto x = sampler.X_0.value; 
		return x * x;
	};


	auto SQR2 = [](auto x)
	{
		return x * x;
	};

	//need a ptr diff  betwee start span and start vecXX

	ApplyTransformUR_X_Impl_EX(a_vec, spnPlus, SQR, identity_sampler, 0, int(spn.paddedSize()));

	ApplyTransformUR_X_Impl_EX(spn, spnPlus, SQR, identity_sampler, 0, int(spn.paddedSize()));

	//trapped into using same container , ie vecxx vecxx  , not mixed so add extra type argument 
	std::vector<FLOAT> ddbbgg = spn;

	transform(SQR2, spn, spnPlus);

	ddbbgg = spnPlus;


	auto isEven = [](auto x) { return !((x - VecXX::INS((2.0) * floor(x * VecXX::INS(0.5)))) > VecXX::INS(0.)); };
	auto evenView = filter(isEven, a_vec);

	ddbbgg = evenView;

	ApplyTransformUR_X_Impl_EX(spn, evenView, SQR, identity_sampler, 100, int(spn.paddedSize()));

	ddbbgg = evenView;


	VecXX anotherVector(v1);
	VecXX resvec(v1);


	SpanXX smallSpan(anotherVector.start() + 20, 100);

	auto vwOfSpan = filter(isEven, smallSpan);

	ddbbgg = vwOfSpan;

	auto CUBE_IT = [](auto x) { return x * x * x; };

	auto neVw = vwOfSpan;
	transformM(CUBE_IT, neVw);

	for (auto& x : neVw)
	{
		std::cout << x << std::endl;
	}

	ddbbgg = neVw;

	ddbbgg = vwOfSpan;

	//transforming view as a result of filtering a span

	vwOfSpan = transform(CUBE_IT, vwOfSpan); // broken somehow ?

	ddbbgg = vwOfSpan;

	for (auto& x : vwOfSpan)
	{
		std::cout << x << std::endl;
	}



}

void doZipping()
{

	auto v1 = getRandomShuffledVector(200);
	VecXX a(v1);
	VecXX b = 2. * a;
	VecXX c = 3. * a;
	VecXX d = 4. * a;
	VecXX e = 5. * a;

	auto zipped = make_Zipped<VecXX::INS>(a, b, c);
	auto it = make_Zipped_itr<VecXX::INS>(a, b, c);


	const VecXX& aa = a;
	const VecXX& bb = b;
	const VecXX& cc = c;
	const VecXX& dd = d;
	const VecXX& ee = e;

	auto zp5 = make_Zipped<VecXX::INS>(aa, bb, cc, dd, ee);
	auto zp5_itr = make_Zipped_itr<VecXX::INS>(aa, bb, cc, dd, ee);

	ignore(zp5);
	ignore(zp5_itr);


	auto addingLambda = [](auto& zpped)
	{
		const auto& X = std::get<0>(zpped.m_registers);
		const auto& Y = std::get<1>(zpped.m_registers);
		const auto& Z = std::get<2>(zpped.m_registers);

		return X + Y + Z;
	};


	VecXX::INS X = 0.;
	for (int i = 0; i < 200 / 8; ++i)
	{
		zipped.load(it);
		auto sum = addingLambda(zipped);
		X += sum;
		it.inc(1);
	}

	std::cout << X[0];

	zipped.load(it);


	auto addingLambda5 = [](auto& zpped)
	{
		const auto& X = std::get<0>(zpped.m_registers);
		const auto& Y = std::get<1>(zpped.m_registers);
		const auto& Z = std::get<2>(zpped.m_registers);
		const auto& L = std::get<3>(zpped.m_registers);
		const auto& M = std::get<4>(zpped.m_registers);

		return X + Y + Z + L + M;
	};


	auto res = transform(addingLambda5, zp5_itr);

	std::vector<FLOAT> vdb = res;

	enum class Access { down = -1, mid = 0, up = 1 };

	Zipped_Reg < VecXX::INS, 3 > myReg;

	auto yy = myReg.get((int)Access::up);
	ignore(yy);



	//////////////////////////////

	VecXX out1 = a;
	VecXX out2 = b;

	VecXX& out1_r = out1;
	VecXX& out2_r = out2;


	auto zp_out = make_Zipped_itr_ref<VecXX::INS>(out1_r, out2_r);

	auto addingLambdaIO = [&](const auto& zpped, auto& out_zipp)
	{
		const auto& X = std::get<0>(zpped.m_registers);
		const auto& Y = std::get<1>(zpped.m_registers);
		const auto& Z = std::get<2>(zpped.m_registers);

		std::get<0>(out_zipp.m_registers) = X + Y + Z;
		std::get<1>(out_zipp.m_registers) = X + Y;
	};


	transform(addingLambdaIO, zp5_itr, zp_out);


	std::vector<FLOAT> vdb1 = out1;
	std::vector<FLOAT> vdb2 = out2;

}


void doMatrix()
{
	VecXX owningVec(0.1, 16 * 10);

	auto val = 0.0;
	for (auto& x : owningVec)
	{
		x = val;
		val++;
	}

	using MAT1 = Layout2D<FLOAT, 8, 0>; //row order  pad/simd size = 8


	MDSpan<FLOAT, MAT1> mat8(owningVec.data(), 10, 10);

	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			std::cout << mat8(i, j) << ",";
		}
		std::cout << "\n";
	}


	std::cout << "  8 \n  \n \n  \n ";
		

	using MAT4 = Layout2D<FLOAT, 4, 0>;


	MDSpan<FLOAT, MAT4> mat4(owningVec.data(), 10, 10);

	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			std::cout << mat4(i, j) << ",";
		}
		std::cout << "\n";
	}


	std::cout << " 4 \n  \n \n  \n ";


	using MAT2 = Layout2D<FLOAT, 2, 0>;


	MDSpan<FLOAT, MAT2> mat2(owningVec.data(), 10, 10);

	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			std::cout << mat2(i, j) << ",";
		}
		std::cout << "\n";
	}



	std::cout << "2 \n  \n \n  \n ";


	using MAT16 = Layout2D<FLOAT, 16, 0>; //row order  pad/simd size = 8

	MDSpan<FLOAT, MAT16> mat16(owningVec.data(), 10, 10);

	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			std::cout << mat16(i, j) << ",";
		}
		std::cout << "\n";
	}

	std::cout << "16 \n  \n \n  \n ";


	////////////////////////


	using MAT11 = Layout2D<FLOAT, 8, 1>; //row order 
	MDSpan<FLOAT, MAT11> mat1_8(owningVec.data(), 10, 10);

	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			std::cout << mat1_8(i, j) << ",";
		}
		std::cout << "\n";
	}


	std::cout << "col  8 \n  \n \n  \n ";


	using MAT14 = Layout2D<FLOAT, 4, 1>; //row order  

	MDSpan<FLOAT, MAT14> mat1_4(owningVec.data(), 10, 10);

	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			std::cout << mat1_4(i, j) << ",";
		}
		std::cout << "\n";
	}


	std::cout << "col 4 \n  \n \n  \n ";


	using MAT12 = Layout2D<FLOAT, 2,1>; //row order  
	MDSpan<FLOAT, MAT12> mat1_2(owningVec.data(), 10, 10);

	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			std::cout << mat1_2(i, j) << ",";
		}
		std::cout << "\n";
	}



	std::cout << "col 2 \n  \n \n  \n ";


	using MAT116 = Layout2D<FLOAT, 16, 1>; //row order  
	MDSpan<FLOAT, MAT116> mat116(owningVec.data(), 10, 10);

	for (int i = 0; i < 10; ++i)
	{
		for (int j = 0; j < 10; ++j)
		{
			std::cout << mat116(i, j) << ",";
		}
		std::cout << "\n";
	}

	std::cout << "col 16 \n  \n \n  \n ";





	
	{

		MDSpan<FLOAT, MAT1> mat(owningVec.data(), 10, 10);


		for (int i = 0; i < 10; ++i)
		{
			for (int j = 0; j < 10; ++j)
			{
				mat(i, j) = i * 100 + j;
			}

		}

		std::cout << "\n";
		std::cout << "\n";

		for (int i = 0; i < 10; ++i)
		{
			for (int j = 0; j < 10; ++j)
			{
				std::cout << mat(i, j) << ",";
			}
			std::cout << "\n";
		}


		Span< VecXX::INS>  spn = getSpan<VecXX::INS>(mat, 0);

		std::vector<FLOAT> vdbg = spn;

		Span< VecXX::INS>  spn2 = getSpan<VecXX::INS>(mat, 1);

		std::vector<FLOAT> vdbg2 = spn2;

		auto add = [](auto x, auto y) {return x + y; };
		FLOAT sum = reduce(spn, add);

		FLOAT sum1 = 0;

		for (int k = 0; k < 10; k++)
		{
			auto  strd_spn = getStridedSpan<VecXX::INS>(mat, 1, k);

			std::vector<FLOAT> vdbg_strd = strd_spn;

			sum1 = reduce(strd_spn, add);
		}


		ignore(sum);
		ignore(sum1);

		// TO DO
		//VecXX dotProductResult(0., spn.size());

		//dotProductResult = transformReduce()

	}


}


