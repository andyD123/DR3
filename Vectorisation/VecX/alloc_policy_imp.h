/****************************  alloc_policy_imp.h   *******************************
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
#include <vector>
#include <unordered_map>


//need a function to reduce size pools to  a minimum
// get rid of magic numbers of byte sizes etc

const int BytesOnCacheLine = 64;
const long MemPoolInitialIncrement = 16;
const long  MemPoolScaleFactor = 2;
const int  ByteAllignment = 64;

template <typename T>
class PoolStrat
{
public:

	PoolStrat(const PoolStrat&) = delete;
	PoolStrat& operator=(const PoolStrat&) = delete;
	PoolStrat& operator=( PoolStrat&&) = delete;
	PoolStrat(PoolStrat&&) = delete;


	explicit PoolStrat(int vecSz) :m_vecSize(vecSz)
	{
		m_sz = 0;
		m_incrementSize = MemPoolInitialIncrement;// 16;
		addToPool(m_incrementSize);
		m_pos = 0;
	}

	~PoolStrat()
	{
		for (auto vec : m_allocatedVecs)
		{
			delete vec;
		}
	}



	T* alloc()
	{
		if (m_pos < (m_sz - 1))
		{
			T* ret = m_memPool[m_pos];
			m_pos++;
			return ret;
		}
		else
		{
			m_incrementSize *= MemPoolScaleFactor;
			addToPool(m_incrementSize);
			return alloc();
		}
	}


	void free(T* pToFree)
	{
		//typically this should be next one down from top of stack
		if ((m_pos <= 0) || (nullptr == pToFree))
		{
			return;
		}

		if (m_memPool[m_pos - 1] == pToFree)
		{
			//pToFree[0] = 666;
			m_pos--;
			return;
		}

		//search for values of i > 0
		int i = m_pos;
		if (i >= static_cast<int>(m_memPool.size()))
		{
			i = static_cast<int>(m_memPool.size()) - 1;
		}
		int maxPos = i;

		for (; i > -1; i--)
		{
			if (m_memPool[i] == pToFree)
			{
				//swap to be freed with top element and  decrement//bubble to the top
				for (int k = i; k < maxPos - 1; k++)
				{
					std::swap(m_memPool[k], m_memPool[k + 1]);
				}
				//pToFree[0] = 666;
				m_pos--;
				return;
			}
		}

	}


	void addToPool(int numElements)
	{
		//m_vecSize for double 64 byte align ie cache line
		size_t offsetAlgn = ByteAllignment;// 64;//   16;
		std::vector<T>* pVecsMem = new std::vector<T>((long)(numElements)*m_vecSize + offsetAlgn);
		m_allocatedVecs.push_back(pVecsMem);

		T* pstrtPt = &((*pVecsMem)[0]);
		while ((reinterpret_cast<long long>(pstrtPt)) % offsetAlgn) pstrtPt++;

		for (int i = 0; i < numElements; i++)
		{
			m_memPool.push_back(pstrtPt);
			pstrtPt += m_vecSize;
		}

		m_sz += numElements;

	}

	inline long pos() const
	{
		return m_pos;
	}

	inline long size() const
	{
		return m_sz;
	}

	const std::vector<std::vector<T>* >& getAllocVecs() const
	{
		m_allocatedVecs;
	}

private:
	long m_pos;
	long m_sz;
	std::vector<T*>  m_memPool;
	long m_incrementSize; // next number of vectors for allocation
	long m_vecSize; //size of element vector considering allignment and padding 
	std::vector<std::vector<T>* >  m_allocatedVecs;

};


//////////////////////////////////////////


template <typename T>
class AllocPolicy
{
	int m_vec_size;
	PoolStrat<T>* m_pool;
public:
	int size() const
	{
		return m_vec_size;
	}

	AllocPolicy(int size) :m_vec_size(size)
	{
		m_pool = new PoolStrat<T>(size);
	}
	~AllocPolicy()
	{
		delete m_pool;
	}


	inline T* alloc()
	{
		return m_pool->alloc();
	}

	inline void free(T* pElement)
	{
		m_pool->free(pElement);
	}

};



template <typename T = double>
class AllAllocators
{
	static int lastSize_N;
	static AllocPolicy<T>* pAllocPolicy;
	static std::unordered_map<int, AllocPolicy<T>*>  m_map_sizeToAllocPolicy;


	static 	void setUpPolicy(int size_N)
	{
		auto itr = m_map_sizeToAllocPolicy.find(size_N);
		if (m_map_sizeToAllocPolicy.end() == itr)
		{
			pAllocPolicy = new AllocPolicy<T>(size_N);
			m_map_sizeToAllocPolicy[size_N] = pAllocPolicy;
		}
	}



public:

	static 	void removePolicy(int size_N)
	{
		auto itr = m_map_sizeToAllocPolicy.find(size_N);
		if (m_map_sizeToAllocPolicy.end() != itr)
		{
			auto policyPtr = m_map_sizeToAllocPolicy[size_N];
			delete policyPtr;
			m_map_sizeToAllocPolicy.erase(itr);
		}
		
	}

	static 	void freeAll()
	{
		for (auto& item : m_map_sizeToAllocPolicy)
		{
			delete item.second;
		}
		m_map_sizeToAllocPolicy.clear();
	}


	static T* alloc(int size_N)
	{
		if (lastSize_N == size_N)
		{
			return  pAllocPolicy->alloc();
		}

		setUpPolicy(size_N);

		pAllocPolicy = m_map_sizeToAllocPolicy[size_N];
		lastSize_N = size_N;
		return pAllocPolicy->alloc();
	}



	static void  free(size_t size_N, T* pMem)
	{
		int sz_N = static_cast<int>(size_N);

		if (lastSize_N == sz_N)
		{
			return  pAllocPolicy->free(pMem);
		}

		setUpPolicy(sz_N);
		pAllocPolicy = m_map_sizeToAllocPolicy[sz_N];
		lastSize_N = sz_N;
		return pAllocPolicy->free(pMem);

	}


};

template< typename T>
struct NumOnCacheLine
{
	static inline int size()
	{
		return BytesOnCacheLine / sizeof(T);
	}
};


template<typename T>
int  getAllignedSizeT(size_t N, T*)
{
	const int M = NumOnCacheLine<T>::size();
	size_t res = (N % M == 0) ? N : (N / M + 1) * M;
	return static_cast<int>(res);
}



template< typename T>
void allocT(size_t& N, T*& pMem)
{
	int n = getAllignedSize(N, pMem);
	N = static_cast<size_t>(n);
	pMem = AllAllocators<T>::alloc(n);
}

template< typename T>
void freeT(size_t N, T* pOld)
{
	//find element and mark as unused 
	return AllAllocators<T>::free(N, pOld);

}

void freeAllAllocators(double);
void freeAllAllocators(float);
void freeAllAllocators(unsigned int);


template <typename T = double>
class AllAllocatorsGuard
{
public:
	~AllAllocatorsGuard()
	{
		freeAllAllocators(T());
	}

};



