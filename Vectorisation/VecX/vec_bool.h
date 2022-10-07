/****************************  vec_bool.h  *******************************
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
#include "alloc_policy.h"
#include "apply_operation.h"

template <typename INS_VEC>
class VecBool
{
private:

	typename InstructionTraits<INS_VEC>::FloatType* m_pData;
	size_t m_size;
	size_t m_implSize;

public:
	VecBool(int sz) :m_size(sz), m_implSize(sz)
	{
		alloc(m_implSize, m_pData);
	}


	~VecBool()
	{
		if (m_pData != nullptr)
		{
			free(m_size, m_pData);
		}
	}

	VecBool(const VecBool& rhs)
	{
		m_size = rhs.m_size;
		m_implSize = m_size;
		alloc(m_implSize, m_pData);
		std::copy(rhs.m_pData, rhs.m_pData + m_implSize, m_pData);
	}


	VecBool& operator=(const VecBool& rhs)
	{
		if (&rhs != this)
		{
			m_size = rhs.m_size;
			m_implSize = rhs.m_implSize;
			std::copy(rhs.m_pData, rhs.m_pData + m_implSize, m_pData);
		}
		return *this;
	}


	VecBool(VecBool&& rhs) noexcept
	{
		m_implSize = 0;
		m_implSize= rhs.m_implSize;
		m_size = rhs.size();
		m_pData = nullptr;
		*this = std::move(rhs);
	}

	VecBool& operator=(VecBool&& rhs) noexcept
	{
		if (&rhs != this)
		{
			std::swap(m_implSize, rhs.m_implSize);
			std::swap(m_size, rhs.m_size);
			std::swap(m_pData, rhs.m_pData);
		}
		return *this;
	}


	inline typename InstructionTraits<INS_VEC>::FloatType* start() const
	{
		return m_pData;
	}

	inline size_t size() const
	{
		return m_size;

	}

	inline size_t paddedSize() const
	{
		return m_implSize;
	}


	inline bool isScalar() const
	{
		return false;
	}

	inline bool operator[](int j)const
	{
		return m_pData[j];
	}


	void setAt(int j, bool val)
	{
		m_pData[j] = val; // need set get function for instruction set
	}

};

