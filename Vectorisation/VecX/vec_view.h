/****************************  vec_view.h  *******************************
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
#include "error_utils.h"
#include <vector>


template <typename INS_VEC>
class VecView
{
private:

	typename InstructionTraits<INS_VEC>::FloatType 	 m_scalarVal;
	bool m_isScalar;

	typename InstructionTraits<INS_VEC>::FloatType* m_pData;
	int m_size; // number of elements represented
	size_t	m_implSize;// actual size of allocated block
	size_t m_implSizeIdx;
	unsigned int* m_pIndex;
	int m_fillSize;
	int m_last;


public:

	VecView() : m_scalarVal(0.),m_isScalar(true)
	{
		
		m_last = 0;
		m_fillSize = 0;
		m_size = 0;
		m_implSize = 0;
		m_implSizeIdx = 0;
		m_pData = nullptr;
		m_pIndex = nullptr;
	}

	explicit VecView(size_t sz): m_scalarVal(InstructionTraits<INS_VEC>::nullValue), m_size(static_cast<int>(sz)), m_implSize(sz), m_implSizeIdx(sz)
	{
		alloc(m_implSize, m_pData);
		alloc(m_implSizeIdx, m_pIndex);
		m_isScalar = false;
		m_scalarVal = 0.0;
		m_fillSize = 0;
		m_last = 0;

	}


	explicit VecView(int sz, const VecView& rhs) : m_scalarVal(InstructionTraits<INS_VEC>::nullValue), m_size(static_cast<int>(sz)), m_implSize(sz), m_implSizeIdx(sz)
	{
		alloc(m_implSize, m_pData);
		alloc(m_implSizeIdx, m_pIndex);
		std::copy(rhs.m_pIndex, rhs.m_pIndex + m_implSizeIdx, m_pIndex);
		m_isScalar = false;
		m_scalarVal = 0.0;
		m_fillSize = 0;
		m_last = 0;

	}


	explicit VecView(typename InstructionTraits<INS_VEC>::FloatType scalarVal) 
		: m_scalarVal(scalarVal),m_isScalar(true),m_implSize(0),m_implSizeIdx(0),
		m_pData(nullptr),m_pIndex(nullptr), m_size(0), m_fillSize(0), m_last(0)
	{
	
	}


	VecView& operator =(typename InstructionTraits<INS_VEC>::FloatType scalarVal)
	{
		m_isScalar = true;
		m_scalarVal = scalarVal;
		if (m_pData != nullptr)
		{
			free(m_implSize, m_pData);
			free(m_implSizeIdx, m_pIndex);
		}
		m_size = 0;
		m_implSize = 0;
		m_implSizeIdx = 0;
		m_pData = nullptr;
		m_pIndex = nullptr;
		m_fillSize = 0;
		m_last = 0;
		return *this;
	}



	~VecView()
	{
		if (m_pData != nullptr)
		{
			free(m_implSize, m_pData);
			free(m_implSizeIdx, m_pIndex);
		}
	}

	VecView(const VecView& rhs)
	{
		m_last = rhs.m_last;
		m_fillSize = rhs.m_fillSize;
		m_isScalar = rhs.m_isScalar;
		m_scalarVal = rhs.m_scalarVal;
		m_size = 0;
		m_implSize = 0;
		m_implSizeIdx = 0;
		m_pData = nullptr;
		m_pIndex = nullptr;

		if (!m_isScalar && (rhs.m_size > 0))
		{
			m_size = rhs.m_size;
			m_implSize = m_size;
			alloc(m_implSize, m_pData);
			std::copy(rhs.m_pData, rhs.m_pData + m_implSize, m_pData);
			m_implSizeIdx = m_size;
			alloc(m_implSizeIdx, m_pIndex);
			std::copy(rhs.m_pIndex, rhs.m_pIndex + m_implSizeIdx, m_pIndex);
		}
	}

	VecView& operator=(const VecView& rhs)
	{
		if (&rhs != this)
		{
			if (m_pData != nullptr)
			{
				free(m_implSize, m_pData);
				free(m_implSizeIdx, m_pIndex);
			}
			m_last = rhs.m_last;
			m_fillSize = rhs.m_fillSize;
			m_isScalar = rhs.m_isScalar;
			m_scalarVal = rhs.m_scalarVal;
			m_size = 0;
			m_implSize = 0;
			m_implSizeIdx = 0;
			m_pData = nullptr;
			m_pIndex = nullptr;

			if (!m_isScalar && (rhs.size() > 0))
			{
				m_size = rhs.m_size;
				m_implSize = m_size;
				alloc(m_implSize, m_pData);
				std::copy(rhs.m_pData, rhs.m_pData + m_implSize, m_pData);
				m_implSizeIdx = m_size;
				alloc(m_implSizeIdx, m_pIndex);
				std::copy(rhs.m_pIndex, rhs.m_pIndex + m_implSizeIdx, m_pIndex);
			}
		}

		return *this;
	}

	VecView(VecView&& rhs) noexcept
	{
		m_last = rhs.m_last;
		m_fillSize = rhs.m_fillSize;
		m_isScalar = false;
		m_scalarVal = InstructionTraits<INS_VEC>::nullValue;
		m_size = static_cast<int>(rhs.size());
		m_pData = nullptr;
		m_pIndex = nullptr;
		m_implSize = 0;
		m_implSizeIdx = 0;
		*this = std::move(rhs);
	}

	VecView& operator=(VecView&& rhs) noexcept
	{
		if (&rhs != this)
		{
			std::swap(m_last, rhs.m_last);
			std::swap(m_fillSize, rhs.m_fillSize);
			std::swap(m_isScalar, rhs.m_isScalar);
			std::swap(m_scalarVal, rhs.m_scalarVal);
			std::swap(m_implSize, rhs.m_implSize);
			std::swap(m_implSizeIdx, rhs.m_implSizeIdx);
			std::swap(m_size, rhs.m_size);
			std::swap(m_pData, rhs.m_pData);
			std::swap(m_pIndex, rhs.m_pIndex);
		}
		return *this;
	}

	//copy view but set values from the given vector
	VecView(const VecView& rhs, const Vec< INS_VEC>& rhsVec)
	{
		m_last = rhs.m_last;
		m_fillSize = rhs.m_fillSize;
		m_isScalar = rhs.m_isScalar;
		m_scalarVal = rhs.m_scalarVal;
		m_size = 0;
		m_implSize = 0;
		m_implSizeIdx = 0;
		m_pData = nullptr;
		m_pIndex = nullptr;

		if (!m_isScalar && (rhs.size() > 0))
		{
			m_size = rhs.m_size;
			m_implSize = m_size;
			alloc(m_implSize, m_pData);
			m_implSizeIdx = m_size;
			alloc(m_implSizeIdx, m_pIndex);
			std::copy(rhs.m_pIndex, rhs.m_pIndex + m_implSizeIdx, m_pIndex);

			for (int i = 0; i < static_cast<int>(m_implSizeIdx); ++i)
			{
				m_pData[i] = rhsVec[rhs.m_pIndex[i]];
			}
		}
	}

	VecView( const Vec< INS_VEC>& rhs)
	{
		m_last = rhs.m_size;
		m_fillSize = rhs.m_size;
		m_isScalar = rhs.m_isScalar;
		m_scalarVal = rhs.m_scalarVal;
		m_size = 0;
		m_implSize = 0;
		m_implSizeIdx = 0;
		m_pData = nullptr;
		m_pIndex = nullptr;

		if (!m_isScalar && (rhs.size() > 0))
		{
			m_size = rhs.m_size;
			m_implSize = m_size;
			alloc(m_implSize, m_pData);
			std::copy(rhs.start(), rhs.start() + m_implSize, m_pData);
			
			m_implSizeIdx = m_size;
			alloc(m_implSizeIdx, m_pIndex);

			for (int i = 0; i < static_cast<int>(m_implSizeIdx); ++i)
			{
				m_pIndex[i] = i;
			}
		}
	}

	//explicit
	operator std::vector<typename InstructionTraits<INS_VEC>::FloatType>()
	{
		return std::vector<typename InstructionTraits<INS_VEC>::FloatType>(begin(), end());
	}

	operator std::vector<typename InstructionTraits<INS_VEC>::FloatType>() const
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


	inline unsigned int* idxStart() const
	{
		return m_pIndex;
	}

	inline typename InstructionTraits<INS_VEC>::FloatType* start() const
	{
		return m_pData;
	}


	//for modern iteration
	inline typename InstructionTraits<INS_VEC>::FloatType* begin() const
	{
		return start();
	}

	//for modern iteration
	inline typename InstructionTraits<INS_VEC>::FloatType* end() const
	{
		return start() + m_last;
	}

	inline int size() const
	{
		return m_last;
	}

	inline int srcSize() const
	{
		return m_size;
	}

	inline void setSizeAndPad(size_t SZ)
	{
		int sz = static_cast<int>(SZ);
		typename InstructionTraits<INS_VEC>::FloatType* ttt = nullptr;
		m_fillSize = getAllignedSize(sz, ttt);
		int cnt = sz;
		m_last = sz;
		while ((cnt < m_fillSize))
		{
			m_pData[cnt] = m_pData[sz - 1];
			m_pIndex[cnt] = m_pIndex[sz - 1];
			cnt++;
		}
	}

	inline int last() const
	{
		return m_last;
	}

	inline void set_last( int newVal) 
	{
		 m_last = newVal;
	}

	inline int fillSize() const
	{
		return m_fillSize;
	}
	
	inline int paddedSize() const
	{
		return static_cast<int>(m_implSize);
	}

	inline bool isScalar() const
	{
		return m_isScalar;
	}

	inline typename InstructionTraits<INS_VEC>::FloatType getScalarValue() const
	{
		return m_scalarVal;
	}

	inline void setScalarValue(typename InstructionTraits<INS_VEC>::FloatType val)
	{
		m_scalarVal = val;
	}


	template<typename T>
	bool isScalar(const Vec<T>&)
	{
		return false;
	}

	template<typename T> bool isScalar(const T&)
	{
		return true;
	}

	template<typename T> bool isScalar(T&)
	{
		return true;
	}

	std::vector<int> getIndex() const
	{
		std::vector<int> ret;
		for (int i = 0; i < m_size; ++i)
		{
			ret.push_back(m_pIndex[i]);	
		}
		return ret;
	}

	void writeView(Vec<INS_VEC>& vec) const
	{

		typename InstructionTraits<INS_VEC>::FloatType* pTrgt = vec.start();
		size_t i = 0;
		for (; i < m_size - (size_t)(4); ++i)
		{
			pTrgt[(m_pIndex[i])] = m_pData[i];
			pTrgt[(m_pIndex[i + 1])] = m_pData[i + 1];
			pTrgt[(m_pIndex[i + 2])] = m_pData[i + 2];
			pTrgt[(m_pIndex[i + 3])] = m_pData[i + 3];
		}

		auto resid = (m_size - 1) % 4;
		switch (resid)
		{
		case 3:
			pTrgt[(m_pIndex[i + 3])] = m_pData[i + 3];
			[[fallthrough]];
		case 2:
			pTrgt[(m_pIndex[i + 2])] = m_pData[i + 2];
			[[fallthrough]];
		case 1:
			pTrgt[(m_pIndex[i + 1])] = m_pData[i + 1];
			[[fallthrough]];
		case 0:
			pTrgt[(m_pIndex[i])] = m_pData[i];
		}
	}
	// do using scatter instructions

	void write(Vec<INS_VEC>& rhs) const
	{

		typename InstructionTraits< INS_VEC>::FloatType* pRes = rhs.start();
		int constexpr width = InstructionTraits< INS_VEC>::width;
		int constexpr unrollFactor = 4;
		using  IdxType = typename InstructionTraits< INS_VEC>::IdxType;




		auto pIdx = &m_pIndex[0];
		int i = 0;
		int SZ = (m_last);

		if constexpr (InstructionTraits< INS_VEC>::useScatter)
		{


			INS_VEC r0;
			INS_VEC r1;
			INS_VEC r2;
			INS_VEC r3;

			IdxType idx0;
			IdxType idx1;
			IdxType idx2;
			IdxType idx3;

			uint32_t limit = 100000;



			for (; i < SZ - (int)(width * unrollFactor); i += width * unrollFactor)
			{
				idx0.load(pIdx + i);
				r0.load_a(m_pData + i);
				scatter(idx0, limit, r0, pRes);

				idx1.load(pIdx + i + width);
				r1.load_a(m_pData + i + width);
				scatter(idx1, limit, r1, pRes);

				idx2.load(pIdx + i + 2 * width);
				r2.load_a(m_pData + i + 2 * width);
				scatter(idx2, limit, r2, pRes);

				idx3.load(pIdx + i + 3 * width);
				r3.load_a(m_pData + i + 3 * width);
				scatter(idx3, limit, r3, pRes);

			}
		}
		
		for (; i < SZ; ++i)
		{
			pRes[(m_pIndex[i])] = m_pData[i];
		}

	}

	void write(VecView<INS_VEC>& rhs) const
	{

		typename InstructionTraits< INS_VEC>::FloatType* pRes = rhs.start();
		int constexpr width = InstructionTraits< INS_VEC>::width;
		int constexpr unrollFactor = 4;
		using  IdxType = typename InstructionTraits< INS_VEC>::IdxType;

		INS_VEC r0;
		INS_VEC r1;
		INS_VEC r2;
		INS_VEC r3;

		IdxType idx0;
		IdxType idx1;
		IdxType idx2;
		IdxType idx3;


		auto pIdx = &m_pIndex[0];

		uint32_t limit = 100000;

		int i = 0;
		int SZ = (m_last);
		for (; i < SZ - (int)(width * unrollFactor); i += width * unrollFactor)
		{
			idx0.load(pIdx +i);
			r0.load_a(m_pData + i);
			scatter(idx0, limit, r0, pRes);

			idx1.load(pIdx + i + width);
			r1.load_a(m_pData + i + width);
			scatter(idx1, limit, r1, pRes);

			idx2.load(pIdx + i + 2 * width);
			r2.load_a(m_pData + i + 2 * width);
			scatter(idx2, limit, r2, pRes);

			idx3.load(pIdx + i + 3 * width);
			r3.load_a(m_pData + i + 3 * width);
			scatter(idx3, limit, r3, pRes);

		}

		for (; i < SZ; ++i)
		{
			pRes[(m_pIndex[i])] = m_pData[i];
		}
	}
}; //class



template< typename INS_VEC>
Vec<INS_VEC>  merge(std::tuple<VecView<INS_VEC>, VecView<INS_VEC> >& src)
{

	if (
		std::get<0>(src).isScalar()
		|| std::get<1>(src).isScalar()
		)
	{
		throw std::exception("cant merge scalar views");
	}

	Vec<INS_VEC> ret(std::get<0>(src).srcSize());
	const VecView<INS_VEC>& trueData = std::get<0>(src);
	const VecView<INS_VEC>& falseData = std::get<1>(src);
	falseData.write(ret);
	trueData.write(ret);
	return ret;
}

template< typename INS_VEC>
VecView<INS_VEC>  mergeToViews(VecView<INS_VEC>& lhs, VecView<INS_VEC>& rhs)
{
	VecView<INS_VEC> ret(static_cast<size_t>(lhs.srxSize()));
	lhs.write(ret);
	rhs.write(ret);
	return ret;
}


