#pragma once
#include "instruction_traits.h"

// Span
//  pretty much like std span with a few minor tweaks
// 
//  The  layout and MDSpan are under construction, 
//  efforts to get a simple and useful array function working with SIMD  
//  and the various transforms
// 
// 
//  padded size == actual size  since its not padded
template<typename  INS_VEC>
struct Span
{

	using T = typename InstructionTraits<INS_VEC>::FloatType;

	Span(const T&) {}

	Span(T* pdata, size_t extent) :m_pstart(pdata), m_extent( extent) {}


	//explicit
	operator std::vector<typename InstructionTraits<INS_VEC>::FloatType>()
	{
		return std::vector<typename InstructionTraits<INS_VEC>::FloatType>(start(), start()+ m_extent);
	}

	const T* begin() const
	{
		return m_pstart;
	}

	T* begin() 
	{
		return m_pstart;
	}


	const T* end() const
	{
		return m_pstart + m_extent;
	}


	T& operator[](size_t pos)
	{
		return *(m_pstart + pos);
	}

	const T& operator[](size_t pos)const
	{
		return *(m_pstart + pos);
	}

	size_t size() const
	{
		return m_extent;
	}

	//front()  back()


	template< size_t N>
    Span<T> first() const 
	{
		return span(m_pstart, N);
	}


	template< size_t N>
	Span<T> last() const
	{
		return span(m_pstart + m_extent - N, N);
	}

	bool empty() const
	{
		return m_extent = 0;
	}

	T* start() const
	{
		return m_pstart;
	}


	// spans are user defined in length so we
	// dont have padded size
	// also they dont own memory so cant guarantee to alloc
	// padding space
	size_t paddedSize() const
	{
		return m_extent;
	}


	constexpr bool isScalar() const
	{
		return false;
	}

	typename InstructionTraits<INS_VEC>::FloatType getScalarValue()const
	{
		return InstructionTraits<INS_VEC>::nullValue;
	}


private:
	T* m_pstart;
	size_t m_extent;
};





struct SpanIncrement
{
	long long m_start;// start
	long long  m_end;// end

};


template <typename T>
Span<T> operator + (const Span<T>& span, const SpanIncrement& inc)
{
	return { span.m_pstart + inc.m_start , span.size() + inc.m_end };
}



template <typename T>
Span<T> operator += (const Span<T>& span, const SpanIncrement& inc)
{
	 span.m_pstart += inc.m_start;
	 span.size() += inc.m_end;
	 return span;
}



/////////////////////////////////////////
/// under development //////////////
template<typename  INS_VEC>
struct StridedSpan 
{

	using T = typename InstructionTraits<INS_VEC>::FloatType;

	StridedSpan(const T&) {}

	StridedSpan(T* pdata, size_t extent, int stride) :m_pstart(pdata), m_extent(extent), m_stride(stride){}


	//explicit
	operator std::vector<typename InstructionTraits<INS_VEC>::FloatType>()
	{
		//return std::vector<typename InstructionTraits<INS_VEC>::FloatType>(start(), start() + m_extent);
		using Float = typename InstructionTraits<INS_VEC>::FloatType;
		std::vector<Float> ret;
		ret.reserve((m_extent + 1) / m_stride);

		for (Float* it = start(); it < start() + m_extent; it += m_stride)
		{
			ret.emplace_back(*it);
		}
		return ret;
	}

	const T* begin() const
	{
		return m_pstart;
	}

	T* begin()
	{
		return m_pstart;
	}


	const T* end() const
	{
		return m_pstart + m_extent;
	}


	T& operator[](size_t pos)
	{
		return *(m_pstart + pos* m_stride);
	}

	const T& operator[](size_t pos)const
	{
		return *(m_pstart + pos * m_stride);
	}

	size_t size() const
	{
		return m_extent;
	}

	//front()  back()


	template< size_t N>
	StridedSpan<T> first() const
	{
		return StridedSpan(m_pstart, N * m_stride, m_stride);
	}


	template< size_t N>
	StridedSpan<T> last() const
	{
		return StridedSpan(m_pstart + m_extent - N * m_stride, N * m_stride);
	}

	bool empty() const
	{
		return m_extent = 0;
	}

	T* start() const
	{
		return m_pstart;
	}


	// spans are user defined in length so we
	// dont have padded size
	// also they dont own memory so cant guarantee to alloc
	// padding space
	// does this apply to strided spans ????
	size_t paddedSize() const
	{
		return m_extent;
	}


	constexpr bool isScalar() const
	{
		return false;
	}

	typename InstructionTraits<INS_VEC>::FloatType getScalarValue()const
	{
		return InstructionTraits<INS_VEC>::nullValue;
	}

	int stride() const 
	{
		return m_stride;
	}

private:
	T* m_pstart;
	size_t m_extent;
	int m_stride;
};




constexpr int ROW_LAYOUT = 0;
constexpr int COL_LAYOUT = 1;
// Layout maps user defined index  structure to an offset

template<typename T, size_t SIMD_SZ, size_t aligned_extent =0>
struct Layout2D
{
	T* m_pAlignedStart;
	size_t m_SimdSize;
	size_t numSIMDS;
	size_t m_rows;
	size_t m_cols;
	size_t m_extent;

	bool isRowOrder;

	T* dataRef() 
	{ 
		return m_pAlignedStart;
	}
	
	Layout2D(T* pdata, size_t  rows, size_t cols):m_pAlignedStart(pdata),m_rows(rows),m_cols(cols)
	{
		if constexpr  (aligned_extent == 0)
		{
			isRowOrder = true;
			numSIMDS = static_cast<int>(m_rows / SIMD_SZ);
			if (m_rows % SIMD_SZ > 0) numSIMDS++;
			m_SimdSize = numSIMDS * SIMD_SZ;
			m_extent = m_SimdSize * m_cols;
		}
		else
		{
			isRowOrder = false;
			numSIMDS =  static_cast<int>(m_cols / SIMD_SZ);
			if (m_cols % SIMD_SZ > 0) numSIMDS++;

			m_SimdSize = numSIMDS * SIMD_SZ;
			m_extent = m_SimdSize * m_rows;

		}

		
	}

	inline size_t stride(size_t extent) const
	{
		
		if constexpr (aligned_extent == 0)
		{
			//return (aligned_extent == extent) ? SIMD_SZ : SIMD_SZ * m_cols;
			return 1;
		}
		else
		{
			//return (aligned_extent == extent) ? SIMD_SZ : SIMD_SZ * m_rows;
			return m_extent;
		}
		
		
	}


	 inline size_t getArrayPos(size_t  row, size_t col)
	{
		 if constexpr  (aligned_extent == 0)
		{
			return col + m_SimdSize * row;
		}
		else
		{
			return  row + m_SimdSize * col;
		}
	}
	
	const T& operator() (size_t row, size_t col) const
	{
		return *(m_pAlignedStart + getArrayPos(row, col));
	}

	T& operator() (size_t row, size_t col) 
	{
		return *(m_pAlignedStart + getArrayPos(row, col));
	}



};



//light version
// for 2D simD

// one axis is aligned and the other strided
template <typename T, typename Layout>
struct MDSpan :public Layout
{

	//MDSpan(T* pStart, Layout& lyOut) :Layout(layOut), m_pAlignedStart(pStart) {}
	MDSpan(T* pStart, size_t row, size_t col) :Layout(pStart, row,col), m_pAlignedStart(pStart) {}


	T* data() const
	{
		return m_pAlignedStart;
	}

    /*

	size_t extent(size_t sz) const
	{
		return ::Layout.extent(sz);
	}


	size_t empty() const
	{
		return ::Layout.empty();
	}
    */

private:
	T* m_pAlignedStart;

};




template<typename INS_VEC>
Span<INS_VEC>  getSpan(Layout2D< typename InstructionTraits<INS_VEC>::FloatType, InstructionTraits<INS_VEC>::width >& layout, size_t pos)
{
	return Span<INS_VEC>(layout.dataRef() + layout.getArrayPos(pos,0), layout.isRowOrder? layout.m_rows : layout.m_cols);
}



template<typename INS_VEC>
StridedSpan<INS_VEC>  getStridedSpan(Layout2D< typename InstructionTraits<INS_VEC>::FloatType, InstructionTraits<INS_VEC>::width >& layout, size_t extnt_id, size_t pos)
{


	return StridedSpan<INS_VEC>(layout.dataRef() + layout.getArrayPos(0,pos), layout.m_extent -pos, layout.isRowOrder ? layout.m_SimdSize : layout.m_cols);
}


