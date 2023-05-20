#pragma once

#include "vec.h"
#include "instruction_traits.h"
#include "error_utils.h"


#include <array>
#include <utility>




//usage  read only 

/*

we have a function of several variables to vectorise

first create the zip Reg

auto make_Zipped(Vexx1 Vex2, Vex3 ....)

define lambda 

auto mybigLambda  =[]( auto& Zipped)
{
	auto P = get<Zipped::PRESSURE>(Zipped);
	auto T = get<Zipped::TEMPERATURE>(Zipped);
	auto V = get<Zipped::VOLUME>(Zipped);

	INS_VEC R =1.3111e-8;
	
	return V*P/(R*T);
}

auto my_zipped = make_Zipped( enum{ PRESSURE,TEMP,VOL}, {Vex1, Vex2, Vex3}, );

auto result = transform( gasPVT,mybigLambda);


*/


template <typename INS_VEC, int N>
struct Zipped_ITR
{
	//align to do
	using FT = typename InstructionTraits<INS_VEC>::FloatType;
	std::array< FT*, N> m_itrs;

	int m_offset = 0;

	int  setOffset(int i)
	{
		int prevOffset = m_offset;
		m_offset = i;

		incOffset();
		return prevOffset;
	}


	void incOffset()
	{
		for (size_t i = 0; i < m_itrs.size(); ++i)
		{
			m_itrs[i] += m_offset;
		}
	}



	//increase position whole number of registers
	void inc(int reg)
	{
		for (size_t i = 0; i < m_itrs.size();++i)
		{
			m_itrs[i]+= (reg * InstructionTraits<INS_VEC>::width);
		}
	}

	int m_size = 0;
	int m_paddedSize = 0;

	int paddedSize() const
	{
		return m_paddedSize;
	}

	int size() const
	{
		return m_size;
	}

};

template <typename INS_VEC, int N, typename NAMED_INDEX >
struct  Named_Zip_iter : public Zipped_ITR< INS_VEC,N>, public  NAMED_INDEX
{

};




template <typename INS_VEC, const int N>
struct Zipped_Reg
{
	//align to do
	
	 std::array<INS_VEC, N> m_registers;


	//std::get  ?

	template<int M>
	INS_VEC& get()
	{
		return std::get<M>(m_registers);
	}

	template<int M>
	constexpr INS_VEC get() const
	{
		return  std::get<M>(m_registers);
	}

	//these might need unroll
	void load(Zipped_ITR<INS_VEC, N>& it)
	{
		for( size_t i= 0; i < N;++i  )
		{
			auto& reg = m_registers[i];
			reg.load(it.m_itrs[i]);
		}
	}


	void store(Zipped_ITR<INS_VEC, N>& it)
	{
		for (size_t i = 0; i < N; ++i)
		{
			auto& reg = m_registers[i];
			reg.store(it.m_itrs[i]);
		}
	}

};






/**/

template <typename  SAMPLER, typename INS_VEC,  int N, int M>
struct Sampled_Zipped_Reg
{
	//align to do

	std::array< SAMPLER, N> m_sampler;


	//st//d::array< SAMPLE<INS_VEC> , N> m_registers;




	//TO DO

	template<int M, int N>
	INS_VEC& get()
	{
		return std::get<N>(std::get<M>(m_sampler));
	}

	template<int M, int N>
	constexpr INS_VEC get() const
	{
		return  std::get<N>(std::get<M>(m_sampler));
	}

	//these might need unroll
	void load(Zipped_ITR<INS_VEC, N>& it)
	{
		for (size_t i = 0; i < N; ++i)
		{
			auto& reg = m_sampler[i];
			reg.load(it.m_itrs[i]);
		}
	}


	//dont store sampler
	void store(Zipped_ITR<INS_VEC, N>& it)
	{
		
		//for (size_t i = 0; i < N; ++i)
		//{
		//	auto& reg = m_sampler[i];
		//	//reg.store(it.m_itrs[i]);
		//}
		
	}

};





//template <typename INS_VEC,  int N, class NAMED_INDEX>
//struct Named_Zipped_Reg : public Zipped_Reg<INS_VEC,N>
//{
//	constexpr static NAMED_INDEX idx;
//};
//noyt good idea



template < typename INS_VEC,   typename  ...T>
constexpr auto make_Zipped( T& ...args)
{
	
//	using FT = typename InstructionTraits<INS_VEC>::FloatType;

	constexpr int N = sizeof...(T);
	Zipped_Reg<INS_VEC, N> ret;
	return ret;

	

}


/*
template < typename INS_VEC, typename  ...T>
constexpr auto make_Zipped_itr(typename T& ...args)
{
	constexpr int N = sizeof...(T);
	Zipped_ITR<INS_VEC, N > zpd = { {(std::begin(args), ...) } };
	return zpd;

}
*/


template < typename INS_VEC, typename  ...T>
constexpr auto make_Zipped_itr(  T& ...args)
{
//	using FT = typename InstructionTraits<INS_VEC>::FloatType;

	constexpr int N = sizeof...(T);
	Zipped_ITR<INS_VEC, N > zpd;

	int i = 0;
	int size = 0;
	int padded_size = 0;
	auto getStart = [&](const auto& x)
	{	 
		 zpd.m_itrs[i]= x.begin();

		 assert((padded_size == 0) ||((padded_size > 0) && (padded_size == x.paddedSize()) ) );
		 assert((size == 0) || ((size > 0) && (size == x.size())));
		 size = x.size();
		 padded_size = x.paddedSize();

		++i;
	};

	(getStart(args), ...);

	zpd.m_size = size;
	zpd.m_paddedSize = padded_size;

	return zpd;

}


template < typename INS_VEC, typename  ...T>
 auto make_Zipped_itr_ref(T& ...args)
{
//	using FT = typename InstructionTraits<INS_VEC>::FloatType;

	constexpr int N = sizeof...(T);
	Zipped_ITR<INS_VEC, N > zpd;

	int i = 0;
	int size = 0;
	int padded_size = 0;
	auto getStart = [&]( auto& x)
	{
		zpd.m_itrs[i] = x.begin();

		assert((padded_size == 0) || ((padded_size > 0) && (padded_size == x.paddedSize())));
		assert((size == 0) || ((size > 0) && (size == x.size())));
		size = x.size();
		padded_size = x.paddedSize();

		++i;
	};

	(getStart(args), ...);

	zpd.m_size = size;
	zpd.m_paddedSize = padded_size;

	return zpd;

}




template < typename INS_VEC>
Zipped_ITR<INS_VEC, 3 >  make_Zipped_itr(const Vec<INS_VEC>& v1, const Vec<INS_VEC>&v2, const Vec<INS_VEC>&v3 )
{
	return  { v1.begin(), v2.begin(),v3.begin() };
}




template<typename OP, typename INS_VEC, int N>
Vec<INS_VEC>  transform(OP& oper,  Zipped_ITR< INS_VEC,N>&  zip)
{

	/*
	check_vector(lhs); //TO DO

	INS_VEC LHS;
	INS_VEC RHS;
	INS_VEC RES;
	RHS = rhs;

	if (isScalar(lhs))
	{
		// both are scalar
		LHS = lhs.getScalarValue();
		RES = oper(LHS, RHS);
		return VEC_TYPE<INS_VEC>(RES[0]);
	}
	*/

	Vec<INS_VEC> ret(zip.m_size);


	//int sz = lhs.size();
	//auto pLhs = ret.start();
	//VEC_TYPE<INS_VEC> ret(initTransformer(lhs));// ret(sz);
	auto pRet = ret.start();

	Zipped_Reg< INS_VEC, N> LHS;
	/*
	Zipped_Reg< INS_VEC, N> LHS2;
	Zipped_Reg< INS_VEC, N> LHS3;
	Zipped_Reg< INS_VEC, N> LHS4;
	*/

	INS_VEC RES;
//	INS_VEC RES2;
//	INS_VEC RES3;
//	INS_VEC RES4;


	int i = 0;

	int step = InstructionTraits< INS_VEC>::width;

	int impSZ = zip.m_paddedSize;

	int rhsSZ = impSZ;// -step;
//	int j = 0;
	for (; i < rhsSZ; i += step)// j++)
	{

		LHS.load(zip); //j);
		RES = oper(LHS);
		RES.store_a(pRet + i);
		zip.inc(1);

		/*
		LHS1.load_a(pLhs + i + width);
		RES1 = oper(LHS1, RHS);
		RES1.store_a(pRet + i + width);

		LHS2.load_a(pLhs + i + 2 * width);
		RES2 = oper(LHS2, RHS);
		RES2.store_a(pRet + i + 2 * width);

		LHS3.load_a(pLhs + i + 3 * width);
		RES3 = oper(LHS3, RHS);
		RES3.store_a(pRet + i + 3 * width);
		*/
	}

	/*	
	for (; i <= impSZ - width; i += width)
	{
		LHS.load_a(pLhs + i);
		RES = oper(LHS, RHS);
		RES.store_a(pRet + i);
	}
*/

	//Since vector is padded no odd end bits, just unused end bits 
	return ret;
}



//have samplers to zipped iterator
/*
make zipped
zipped has vector of samplers


SampledZip ; 
 
auto pu = 	get<price,1>->SZ;
auto pm =  	get<price,0>->SZ;
auto pd = 	get<price,-1>->SZ;


*/


template<typename OP, typename INS_VEC, int N,int M>
void  transform(OP& oper, Zipped_ITR< INS_VEC, N> in, Zipped_ITR< INS_VEC, M>& out)  // by value in parameter
{

	//check in and out sizes are same.

	Zipped_Reg< INS_VEC, N> IN;
	Zipped_Reg< INS_VEC, M> OUT;


	int i = 0;

	int step = InstructionTraits< INS_VEC>::width;

	int impSZ = in.m_paddedSize;

	int rhsSZ = impSZ;
	int j = 0;
	for (; i < rhsSZ; i += step, j++)
	{
		IN.load(in); 
		oper(IN, OUT);
		OUT.store(out);
		in.inc(1);
		out.inc(1);
	}


//Since vector is padded no odd end bits, just unused end bits 
//	return ret;
}





// experimental 
// we combine the sampler with constant inputs from in
template<  typename INS_VEC,  int M, int N, typename OP, typename SAMPLER >
void ApplyTransformUR_X_Impl_EX(SAMPLER& sampler,Zipped_ITR< INS_VEC, N> in, Zipped_ITR<INS_VEC, M>& out, OP& oper, int i = 0, int impSZ = -1)
{

	impSZ = (impSZ < 0) ? in.paddedSize() : impSZ;


//	auto pRhs1 = rhs1.start();
//	auto pRet = ret.start();

	const int width = InstructionTraits<INS_VEC>::width;
	
	
	//int step = 4 * width;
	int step =  width;


	//SAMPLER RHS1(sampler);
	//INS_VEC RES;

	//Zipped_Reg< INS_VEC, N> IN;
	Zipped_Reg< INS_VEC, M> OUT;

	Sampled_Zipped_Reg< Sampled_Zipped_Reg,INS_VEC, N,M> IN;


	//we can only get a starting position bigger than  zero when we access points in the 
	// data preceeding the starting point, so we advance to a popint where we sample valid /existing data
	i = i + std::max(0, -sampler.min());


	IN.setOffset(i);
	OUT.setOffset(i);



	//similarly if we are sampling  points beyond current index, we need to reduce maximum value iterated to so
	// that we stay in a valid range 
	impSZ = impSZ - std::max(0, sampler.max());

	int ld_offset = 0;
	int SV_offset = 0;

	int rhsSZ = impSZ - step;
	for (; i < rhsSZ; i += step)
	{
	//	RHS1.load(pRhs1 + i + ld_offset);
	//	RES = oper(RHS1);
	//	RES.store(pRet + i + SV_offset);


		IN.load(in);
		oper(IN, OUT);
		OUT.store(out);
		in.inc(1);
		out.inc(1);


		/*
		RHS2.load(pRhs1 + i + ld_offset + width);
		RES1 = oper(RHS2);
		RES1.store(pRet + i + SV_offset + width);

		RHS3.load(pRhs1 + i + ld_offset + width * 2);
		RES2 = oper(RHS3);
		RES2.store(pRet + i + SV_offset + width * 2);

		RHS4.load(pRhs1 + i + ld_offset + width * 3);
		RES3 = oper(RHS4);
		RES3.store(pRet + i + SV_offset + width * 3);
		*/
	}
	/*
	for (; i <= impSZ - width; i += width)
	{
		RHS1.load(pRhs1 + i + ld_offset);
		RES = oper(RHS1);
		RES.store(pRet + i + SV_offset);
	}
	*/

	

	if (i < (impSZ - width))
	{
		/*
		RHS1.load(pRhs1 + i )//+ ld_offset);
		RES = oper(RHS1);
		RES.store(pRet + i);// +SV_offset);
		*/

		IN.load(in);
		oper(IN, OUT);
		OUT.store(out);
		in.inc(1);
		out.inc(1);
	}

	//move to one register width from last valid
	//point to calculate
	i = impSZ - width;
/*
	RHS1.load(pRhs1 + i);// +ld_offset);
	RES = oper(RHS1);
	RES.store(pRet + i);// +SV_offset);
*/

	IN.setOffset(i);
	OUT.setOffset(i);

	IN.load(in);
	oper(IN, OUT);
	OUT.store(out);
	in.inc(1);
	out.inc(1);

}

