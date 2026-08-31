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
#include <stdexcept>
#include <utility>
#include "vec.h"
#include "vec_bool.h"



template <typename INS_VEC>
class VecD
{
private:
    Vec<INS_VEC> m_value;
    Vec<INS_VEC> m_derivative;

    static void validateShape(const Vec<INS_VEC>& value, const Vec<INS_VEC>& derivative)
    {
        if (value.isScalar() != derivative.isScalar()) {
            throw std::invalid_argument(
                "VecD: primal and derivative must both be scalar or both be vectors");
        }
        if (!value.isScalar() && value.size() != derivative.size()) {
            throw std::invalid_argument(
                "VecD: primal and derivative logical sizes must match");
        }
    }

public:

	VecD()
	{}

	VecD(const VecD&) = default;
	VecD(VecD&&) noexcept = default;

	VecD& operator=(const VecD& rhs)
	{
		if (this != &rhs)
		{
			replaceValueAndDerivative(rhs.m_value, rhs.m_derivative);
		}
		return *this;
	}

	VecD& operator=(VecD&&) noexcept = default;

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
		return makeDVecOnes(value, sz);
	}

	static VecD< INS_VEC> makeDVecZeroV(const typename InstructionTraits<INS_VEC>::FloatType&  value, int sz)
	{
		return makeDVecZero(value, sz);
	}

	explicit VecD(typename InstructionTraits<INS_VEC>::FloatType scalarVal)
		:m_value(scalarVal), m_derivative(InstructionTraits<INS_VEC>::nullValue)
	{

	}


	VecD(typename InstructionTraits<INS_VEC>::FloatType scalarVal, typename InstructionTraits<INS_VEC>::FloatType derivVal)
		:m_value(scalarVal), m_derivative(derivVal)
	{

	}



	VecD(const std::vector< typename InstructionTraits<INS_VEC>::FloatType> & ctr)
		:m_value(ctr),
		 m_derivative(InstructionTraits<INS_VEC>::nullValue, ctr.size())
	{

	}


	VecD(const Vec<INS_VEC>&  value, const Vec<INS_VEC>&  derivative)
		:m_value(value), m_derivative(derivative)
    {
        validateShape(m_value, m_derivative);
    }


	VecD(Vec<INS_VEC>&&  value, Vec<INS_VEC>&&  derivative) :
		m_value(std::forward< Vec<INS_VEC>>(value)),
		m_derivative(std::forward<Vec<INS_VEC>>(derivative))
	{
        validateShape(m_value, m_derivative);
    }


	VecD(Vec<INS_VEC>&&  value) :
		m_value(std::forward< Vec<INS_VEC>>(value))
	{
		if (!m_value.isScalar())
		{
			m_derivative = Vec<INS_VEC>(
				InstructionTraits<INS_VEC>::nullValue, m_value.size());
		}
		else
		{
			m_derivative = InstructionTraits<INS_VEC>::nullValue;
		}
	}


	VecD(Vec<INS_VEC>&&  value, const Vec<INS_VEC>&  d) :
		m_value(std::forward< Vec<INS_VEC>>(value)), m_derivative(d)
	{
		validateShape(m_value, m_derivative);
	}


	//explicit
	VecD(const Vec<INS_VEC>& value) :
		m_value(value),
		m_derivative(value.isScalar()
			? Vec<INS_VEC>(InstructionTraits<INS_VEC>::nullValue)
			: Vec<INS_VEC>(InstructionTraits<INS_VEC>::nullValue, value.size()))
	{
	}
	

	explicit VecD(size_t sz) :m_value(sz), m_derivative(sz)
	{

	}


	typename InstructionTraits<INS_VEC>::FloatType& operator[](size_t pos)
	{
		return m_value[pos];
	}

	typename InstructionTraits<INS_VEC>::FloatType operator[](size_t pos) const
	{
		return m_value[pos];
	}


	inline typename InstructionTraits<INS_VEC>::FloatType* start() const
	{
		return m_value.start();
	}


	inline size_t size() const
	{
		return m_value.size();
	}



	inline int  paddedSize() const
	{
		return static_cast<int>(m_value.paddedSize());
	}

	inline bool isScalar() const
	{
		return m_value.isScalar();
	}

	inline typename InstructionTraits<INS_VEC>::FloatType getScalarValue() const
	{
		return m_value.getScalarValue();
	}

	inline void setScalarValue(typename InstructionTraits<INS_VEC>::FloatType newVal)
	{
		m_value.setScalarValue(newVal);
	}


	inline typename InstructionTraits<INS_VEC>::FloatType getScalarDeriv() const
	{
		return m_derivative.getScalarValue();
	}

	inline void setScalarDeriv(typename InstructionTraits<INS_VEC>::FloatType newVal)
	{
		m_derivative.setScalarValue(newVal);
	}



	// Shape-bearing storage is read-only to callers. Use the validated
	// replacement methods below when either side must change.
	inline const Vec<INS_VEC>& value() const
	{
		return m_value;
	}

	inline const Vec<INS_VEC>& derivative() const
	{
		return m_derivative;
	}

	void replaceValue(Vec<INS_VEC> value)
	{
		validateShape(value, m_derivative);
		m_value = std::move(value);
	}

	void replaceDerivative(Vec<INS_VEC> derivative)
	{
		validateShape(m_value, derivative);
		m_derivative = std::move(derivative);
	}

	void replaceValueAndDerivative(
		Vec<INS_VEC> value, Vec<INS_VEC> derivative)
	{
		// A paired replacement may change the logical shape atomically.
		validateShape(value, derivative);
		m_value = std::move(value);
		m_derivative = std::move(derivative);
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
