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
    Span<INS_VEC> first() const
	{
		const size_t count = std::min(N, m_extent);
		return Span<INS_VEC>(m_pstart, count);
	}


	template< size_t N>
	Span<INS_VEC> last() const
	{
		const size_t count = std::min(N, m_extent);
		return Span<INS_VEC>(m_pstart + m_extent - count, count);
	}

	bool empty() const
	{
		return m_extent == 0;
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

	StridedSpan(T* pdata, size_t extent, size_t stride) :m_pstart(pdata), m_extent(extent), m_stride(stride){}


	//explicit
	operator std::vector<typename InstructionTraits<INS_VEC>::FloatType>()
	{
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

	// Legacy algorithms interpret size()/paddedSize() as the physical extent
	// covered by the strided view. logicalSize() exposes the element count.
	size_t size() const
	{
		return m_extent;
	}

	size_t logicalSize() const
	{
		return m_extent == 0 ? 0 : 1 + (m_extent - 1) / m_stride;
	}

	size_t physicalExtent() const
	{
		return m_extent;
	}

	//front()  back()


	template< size_t N>
	StridedSpan<INS_VEC> first() const
	{
		const size_t count = std::min(N, logicalSize());
		const size_t extent = count == 0 ? 0 : 1 + (count - 1) * m_stride;
		return StridedSpan<INS_VEC>(m_pstart, extent, m_stride);
	}


	template< size_t N>
	StridedSpan<INS_VEC> last() const
	{
		const size_t count = std::min(N, logicalSize());
		const size_t offset = (logicalSize() - count) * m_stride;
		const size_t extent = count == 0 ? 0 : 1 + (count - 1) * m_stride;
		return StridedSpan<INS_VEC>(m_pstart + offset, extent, m_stride);
	}

	bool empty() const
	{
		return m_extent == 0;
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

	size_t stride() const
	{
		return m_stride;
	}

private:
	T* m_pstart;
	size_t m_extent;
	size_t m_stride;
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

	const T* dataRef() const
	{
		return m_pAlignedStart;
	}

	size_t rows() const { return m_rows; }
	size_t cols() const { return m_cols; }
	size_t storageSize() const { return m_extent; }
	size_t rowStride() const { return isRowOrder ? m_SimdSize : 1; }
	size_t columnStride() const { return isRowOrder ? 1 : m_SimdSize; }
	
	Layout2D(T* pdata, size_t rows, size_t cols):m_pAlignedStart(pdata),m_rows(rows),m_cols(cols)
	{
		if constexpr (aligned_extent == 0)
		{
			// Row-major storage pads the contiguous column extent.
			isRowOrder = true;
			numSIMDS = static_cast<int>(m_cols / SIMD_SZ);
			if (m_cols % SIMD_SZ > 0) numSIMDS++;
			m_SimdSize = numSIMDS * SIMD_SZ;
			m_extent = m_SimdSize * m_rows;
		}
		else
		{
			// Column-major storage pads the contiguous row extent.
			isRowOrder = false;
			numSIMDS = static_cast<int>(m_rows / SIMD_SZ);
			if (m_rows % SIMD_SZ > 0) numSIMDS++;
			m_SimdSize = numSIMDS * SIMD_SZ;
			m_extent = m_SimdSize * m_cols;
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
StridedSpan<INS_VEC> makeStridedSpanFromCount(
	typename InstructionTraits<INS_VEC>::FloatType* start,
	size_t count,
	size_t stride)
{
	const size_t extent = count == 0 ? 0 : 1 + (count - 1) * stride;
	return StridedSpan<INS_VEC>(start, extent, stride);
}


// Explicit matrix views. Row-major storage gives contiguous rows and strided
// columns; column-major storage gives the inverse. The numerical code can name
// the axis it wants without knowing the physical layout.
template<typename INS_VEC, size_t SIMD_SZ = InstructionTraits<INS_VEC>::width, size_t EXTENT = 0>
auto getRowSpan(MDSpan<typename InstructionTraits<INS_VEC>::FloatType,
	Layout2D<typename InstructionTraits<INS_VEC>::FloatType, SIMD_SZ, EXTENT> >& layout,
	size_t row)
{
	if constexpr (EXTENT == 0)
		return Span<INS_VEC>(layout.dataRef() + layout.getArrayPos(row, 0), layout.m_cols);
	else
		return makeStridedSpanFromCount<INS_VEC>(
			layout.dataRef() + layout.getArrayPos(row, 0), layout.m_cols, layout.m_SimdSize);
}


template<typename INS_VEC, size_t SIMD_SZ = InstructionTraits<INS_VEC>::width, size_t EXTENT = 0>
auto getColumnSpan(MDSpan<typename InstructionTraits<INS_VEC>::FloatType,
	Layout2D<typename InstructionTraits<INS_VEC>::FloatType, SIMD_SZ, EXTENT> >& layout,
	size_t col)
{
	if constexpr (EXTENT == 0)
		return makeStridedSpanFromCount<INS_VEC>(
			layout.dataRef() + layout.getArrayPos(0, col), layout.m_rows, layout.m_SimdSize);
	else
		return Span<INS_VEC>(layout.dataRef() + layout.getArrayPos(0, col), layout.m_rows);
}


// Compatibility helpers retain the old names. getSpan selects the contiguous
// axis and getStridedSpan selects the orthogonal axis. The legacy extent
// argument is retained for source compatibility.
template<typename INS_VEC, size_t SIMD_SZ = InstructionTraits<INS_VEC>::width, size_t EXTENT = 0>
StridedSpan<INS_VEC> getStridedSpan(MDSpan<typename InstructionTraits<INS_VEC>::FloatType,
	Layout2D<typename InstructionTraits<INS_VEC>::FloatType, SIMD_SZ, EXTENT> >& layout,
	size_t /*extent*/, size_t pos)
{
	if constexpr (EXTENT == 0)
		return getColumnSpan<INS_VEC>(layout, pos);
	else
		return getRowSpan<INS_VEC>(layout, pos);
}


template<typename INS_VEC, size_t SIMD_SZ = InstructionTraits<INS_VEC>::width, size_t EXTENT = 0>
Span<INS_VEC> getSpan(MDSpan<typename InstructionTraits<INS_VEC>::FloatType,
	Layout2D<typename InstructionTraits<INS_VEC>::FloatType, SIMD_SZ, EXTENT> >& layout,
	size_t pos)
{
	if constexpr (EXTENT == 0)
		return getRowSpan<INS_VEC>(layout, pos);
	else
		return getColumnSpan<INS_VEC>(layout, pos);
}

