/****************************  vec.h  *******************************
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
#include "instruction_traits.h"
#include "alloc_policy.h"
#include "apply_operation.h"
#include "vec_view.h"

#include <algorithm>
#include <iterator>
#include <vector>



template <typename INS_VEC>
class Vec
{
public:

	friend class VecView< INS_VEC>;

	typedef INS_VEC INS;
	using SCALA_TYPE = typename InstructionTraits<INS_VEC>::FloatType;

	//convert to scalar
	template<typename T>
	static SCALA_TYPE scalar(const T& val)
	{
		return static_cast<SCALA_TYPE>(val);
	}

	//convert to register 
	template<typename T>
	static INS reg(T& val)
	{
		INS vec(static_cast<SCALA_TYPE>(val));
		return vec;
	}


private:

	typename InstructionTraits<INS_VEC>::FloatType 	 m_scalarVal;
	bool m_isScalar;

	typename InstructionTraits<INS_VEC>::FloatType*  m_pData;

	int m_size; // number of elements represented
	size_t m_implSize;// actual size of allocated block

public:

	Vec():m_scalarVal(0.),m_isScalar(true)
	{
		m_size = 0;
		m_implSize = 0;
		m_pData = nullptr;
	}


	//not explicit  allow conversions
	Vec( typename InstructionTraits<INS_VEC>::FloatType scalarVal):m_scalarVal(scalarVal),m_isScalar(true)
	{
		m_size = 0;
		m_implSize = 0;
		m_pData = nullptr;
	}


	Vec& operator =(typename InstructionTraits<INS_VEC>::FloatType scalarVal) 
	{
		m_isScalar=true;
		m_scalarVal=scalarVal;
		if (m_pData != nullptr)
		{
			free(m_implSize, m_pData);
		}
		m_size = 0;
		m_implSize = 0;
		m_pData = nullptr;
		return *this;
	}
	

	Vec(const std::vector<  typename InstructionTraits<INS_VEC>::FloatType > & ctr)
	{

		int sz = static_cast<int>(std::distance(ctr.begin(), ctr.end()) );
		m_size =sz;
		m_implSize = sz;
		alloc(m_implSize,m_pData);

		auto repeatedPaddingValue = ctr.at(sz-1);
		for(auto s =sz; s < static_cast<int>(m_implSize);s++)
		{
			m_pData[s] =repeatedPaddingValue;
		}

		std::copy(cbegin(ctr),cend(ctr),m_pData);

		m_isScalar = false;
		m_scalarVal = 0.0;
	}

	explicit Vec(int sz) :m_size(sz), m_implSize(sz)
	{
		alloc(m_implSize, m_pData);
		m_isScalar = false;
		m_scalarVal = 0.0;
	}



	Vec(typename InstructionTraits<INS_VEC>::FloatType val, int sz) :m_size(sz), m_implSize(sz)
	{
		alloc(m_implSize, m_pData);
		m_isScalar = false;
		m_scalarVal = 0.0;

		std::fill_n(start(), sz, val);

	}


	~Vec()
	{
		if(m_pData != nullptr)
		{
			free(m_implSize,m_pData);
		}
	}

	Vec(const Vec& rhs):  m_scalarVal(rhs.m_scalarVal), m_isScalar(rhs.m_isScalar), m_size(rhs.m_size), m_implSize(rhs.m_implSize)
	{
		m_pData = nullptr;

		if( !m_isScalar)
		{
			m_implSize = m_size;
			alloc(m_implSize,m_pData);
			std::copy(rhs.m_pData, rhs.m_pData+ m_implSize , m_pData);
		}
	}

	Vec& operator=(const Vec& rhs)
	{
		if (&rhs != this)
		{
			if (m_pData != nullptr)
			{
				free(m_implSize, m_pData);
				m_pData = nullptr;
				m_size = 0;
				m_implSize = 0;
			}

			m_isScalar = rhs.m_isScalar;
			m_scalarVal = rhs.m_scalarVal;

			if( !m_isScalar)
			{
				m_size= rhs.m_size;
				m_implSize = m_size;
				alloc(m_implSize,m_pData);
				std::copy(rhs.m_pData, rhs.m_pData+ m_implSize , m_pData);
			}
		}

	   return *this;
	}


	Vec(Vec&& rhs) noexcept
	{
		m_implSize = 0;
		m_implSize= rhs.m_implSize;
		m_isScalar = false;
		m_scalarVal = InstructionTraits<INS_VEC>::nullValue;
		m_size =rhs.size();
		m_pData = nullptr;
		*this = std::move(rhs);
	}

	Vec& operator=( Vec&& rhs) noexcept
	{
		if (&rhs != this)
		{
			std::swap(m_isScalar , rhs.m_isScalar);
			std::swap(m_scalarVal , rhs.m_scalarVal);
			std::swap( m_implSize , rhs.m_implSize);
			std::swap(m_size, rhs.m_size);
			std::swap(m_pData, rhs.m_pData);
		}
		return *this;
	}

	
	//explicit
	operator std::vector<typename InstructionTraits<INS_VEC>::FloatType>()
	{
		return std::vector<typename InstructionTraits<INS_VEC>::FloatType>(begin(), end());
	}


	typename InstructionTraits<INS_VEC>::FloatType& operator[](size_t pos) 
	{
		return m_pData[pos];
	}

	typename InstructionTraits<INS_VEC>::FloatType operator[](size_t pos) const
	{
		return m_pData[pos];
	}


	inline typename InstructionTraits<INS_VEC>::FloatType* start() const
	{
		return m_pData;
	}


	inline int size() const
	{
		return m_size;
	}

	
	inline int paddedSize() const
	{
		return  static_cast<int>(m_implSize);
	}

	inline bool isScalar() const
	{
		return m_isScalar;
	}

	inline typename InstructionTraits<INS_VEC>::FloatType getScalarValue() const
	{
		return m_scalarVal;
	}

	inline void setScalarValue( typename InstructionTraits<INS_VEC>::FloatType val)
	{
		m_scalarVal = val;
	}

	inline typename InstructionTraits<INS_VEC>::FloatType* begin() const
	{
		return start();
	}
	
	inline typename InstructionTraits<INS_VEC>::FloatType* end() const
	{
		return start() + m_size;
	}

	inline static  INS_VEC reg(typename InstructionTraits<INS_VEC>::FloatType val)
	{
		return INS_VEC(val);
	}

};




template<typename T>
bool isScalar(const Vec<T> & X)
{
	return  X.isScalar();
}

template<typename T>
bool isScalar( Vec<T>& X)
{
	return  X.isScalar();
}


template<typename T>
bool isScalar(const VecView<T>& X)
{
	return  X.isScalar();
}

template<typename T>
bool isScalar( VecView<T>& X)
{
	return  X.isScalar();
}


template<typename T> bool isScalar(const T&)
{
	return true;
}

template<typename T> bool isScalar( T&)
{
	return true;
}

