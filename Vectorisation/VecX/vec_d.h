/****************************  vec_d.h  *******************************
* Author:        Andrew Drakeford
* Date created:  2021-04-10
* Last modified: 2021-04-10
* Version:       1.0
* Project:       DR Cubed
* Description:
*
* (c) Copyright 2019 Andrew Drakeford
* Apache License version 2.0 or later.
*****************************************************************************/
#pragma once
#include <algorithm>
#include <iterator>
#include "vec.h"
#include "vec_bool.h"



template <typename INS_VEC>
class VecD
{
private:
public:

	Vec< INS_VEC> val;
	Vec< INS_VEC> deriv;

public:

	VecD()
	{}

	static VecD< INS_VEC> makeDVecZero(const Vec<INS_VEC>&  value)
	{
		if (value.isScalar())
		{
			return  VecD(value.getScalarValue(), InstructionTraits<INS_VEC>::nullValue);
		}
		std::vector< typename InstructionTraits<INS_VEC>::FloatType> zeros(value.size(), InstructionTraits<INS_VEC>::nullValue);
		return VecD(value, Vec< INS_VEC>(zeros));
	}

	static VecD< INS_VEC> makeDVecOnes(const Vec<INS_VEC>&  value)
	{
		if (value.isScalar())
		{
			return  VecD(value.getScalarValue(), InstructionTraits<INS_VEC>::oneValue);
		}
		std::vector< typename InstructionTraits<INS_VEC>::FloatType> ones(value.size(), InstructionTraits<INS_VEC>::oneValue);
		return VecD(value, ones);
	}


	static VecD< INS_VEC> makeDVecOnes(const typename InstructionTraits<INS_VEC>::FloatType&  value, int sz)
	{
		Vec< INS_VEC> values(value, sz);
		Vec< INS_VEC> ones(InstructionTraits<INS_VEC>::oneValue, sz);
		return VecD(values, ones);
	}

	static VecD< INS_VEC> makeDVecZero(const typename InstructionTraits<INS_VEC>::FloatType&  value, int sz)
	{
		Vec< INS_VEC> values(value, sz);
		Vec< INS_VEC> zeros(InstructionTraits<INS_VEC>::nullValue, sz);
		return VecD(values, zeros);
	}

	static VecD< INS_VEC> makeDVecOnesV(const typename InstructionTraits<INS_VEC>::FloatType&  value, int sz)
	{
		Vec< INS_VEC> values(value, sz);
		std::vector< typename InstructionTraits<INS_VEC>::FloatType> ones(value.size(), InstructionTraits<INS_VEC>::oneValue);
		return VecD(values, ones);
	}

	static VecD< INS_VEC> makeDVecZeroV(const typename InstructionTraits<INS_VEC>::FloatType&  value, int sz)
	{
		Vec< INS_VEC> values(value, sz);
		std::vector< typename InstructionTraits<INS_VEC>::FloatType> nulls(value.size(), InstructionTraits<INS_VEC>::nullValue);
		return VecD(values, nulls);
	}

	explicit VecD(typename InstructionTraits<INS_VEC>::FloatType scalarVal)
		:val(scalarVal), deriv(InstructionTraits<INS_VEC>::nullValue)
	{

	}


	VecD(typename InstructionTraits<INS_VEC>::FloatType scalarVal, typename InstructionTraits<INS_VEC>::FloatType derivVal)
		:val(scalarVal), deriv(derivVal)
	{

	}



	VecD(const std::vector< typename InstructionTraits<INS_VEC>::FloatType> & ctr) :val(ctr), deriv(InstructionTraits<INS_VEC>::nullValue, ctr.size())
	{

	}


	VecD(const Vec<INS_VEC>&  value, const Vec<INS_VEC>&  derivative) : val(value), deriv(derivative)
	{}


	VecD(Vec<INS_VEC>&&  value, Vec<INS_VEC>&&  derivative) :
		val(std::forward< Vec<INS_VEC>>(value)),
		deriv(std::forward<Vec<INS_VEC>>(derivative))
	{}


	VecD(Vec<INS_VEC>&&  value) :
		val(std::forward< Vec<INS_VEC>>(value))
	{
		if (!val.isScalar())
		{
			deriv(InstructionTraits<INS_VEC>::nullVal, value.size());
		}
		else
		{
			deriv(InstructionTraits<INS_VEC>::nullVal);
		}
	}


	VecD(Vec<INS_VEC>&&  value, const Vec<INS_VEC>&  d) :
		val(std::forward< Vec<INS_VEC>>(value)), deriv(d)
	{
	}


	//explicit
	VecD(const Vec<INS_VEC>& value) :
		val(value), deriv(value.isScalar()? InstructionTraits<INS_VEC>::nullValue : Vec< INS_VEC>(InstructionTraits<INS_VEC>::nullValue,value.size()  ) )
	{
	}
	

	explicit VecD(size_t sz) :val(sz), deriv(sz)
	{

	}


	typename InstructionTraits<INS_VEC>::FloatType& operator[](size_t pos)
	{
		return val[pos];
	}

	typename InstructionTraits<INS_VEC>::FloatType operator[](size_t pos) const
	{
		return val[pos];
	}


	inline typename InstructionTraits<INS_VEC>::FloatType* start() const
	{
		return val.start();
	}


	inline size_t size() const
	{
		return val.size();
	}



	inline int  paddedSize() const
	{
		return static_cast<int>(val.paddedSize());
	}

	inline bool isScalar() const
	{
		return val.isScalar();
	}

	inline typename InstructionTraits<INS_VEC>::FloatType getScalarValue() const
	{
		return val.getScalarValue();
	}

	inline void setScalarValue(typename InstructionTraits<INS_VEC>::FloatType newVal)
	{
		val.setScalarValue(newVal);
	}


	inline typename InstructionTraits<INS_VEC>::FloatType getScalarDeriv() const
	{
		return deriv.getScalarValue();
	}

	inline void setScalarDeriv(typename InstructionTraits<INS_VEC>::FloatType newVal)
	{
		deriv.setScalarValue(newVal);
	}



	inline const Vec< INS_VEC>& value() const
	{
		return val;
	}

	inline const Vec< INS_VEC>& derivative() const
	{
		return deriv;
	}

	inline  Vec< INS_VEC>& value()
	{
		return val;
	}

	inline  Vec< INS_VEC>& derivative()
	{
		return deriv;
	}


};


template <typename INS_VEC>
VecD<INS_VEC> D(const Vec<INS_VEC>& rhs)
{
	return VecD<INS_VEC>::makeDVecOnes(rhs);
}


template <typename INS_VEC>
VecD<INS_VEC> C(const Vec<INS_VEC>& rhs)
{
	return VecD<INS_VEC>::makeDVecZero(rhs);
}

