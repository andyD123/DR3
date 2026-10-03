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
#include <algorithm>
#include <memory>
#include <mutex>
#include <stdexcept>
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
		if (m_pos < m_sz)
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
		if (nullptr == pToFree)
		{
			return;
		}

		//typically this should be next one down from top of stack
		if (m_pos <= 0)
		{
			reportInvalidFree(pToFree);
			return;
		}

		if (m_memPool[m_pos - 1] == pToFree)
		{
			//pToFree[0] = 666;
			m_pos--;
			return;
		}

		//search for values of i > 0
		int i = static_cast<int>(m_pos) - 1;
		const int maxPos = static_cast<int>(m_pos);

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

		reportInvalidFree(pToFree);
	}


	void addToPool(int numElements)
	{
		constexpr size_t offsetAlgn = ByteAllignment;
		// Reserve every potentially throwing container growth before publishing
		// the new block into the pool.
		m_allocatedVecs.reserve(m_allocatedVecs.size() + 1);
		m_memPool.reserve(m_memPool.size() + static_cast<size_t>(numElements));
		auto storage = std::make_unique<std::vector<T>>(
			static_cast<size_t>(numElements) * static_cast<size_t>(m_vecSize)
			+ ByteAllignment);

		T* pstrtPt = storage->data();
		while ((reinterpret_cast<long long>(pstrtPt)) % offsetAlgn) pstrtPt++;
		std::vector<T*> newEntries;
		newEntries.reserve(static_cast<size_t>(numElements));

		for (int i = 0; i < numElements; i++)
		{
			newEntries.push_back(pstrtPt);
			pstrtPt += m_vecSize;
		}

		m_allocatedVecs.push_back(storage.get());
		m_memPool.insert(m_memPool.end(), newEntries.begin(), newEntries.end());
		storage.release();
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
		return m_allocatedVecs;
	}

private:
	void reportInvalidFree(T* pointer) const
	{
		if (!dr3AllocatorDiagnosticsEnabled())
		{
			return;
		}
		const auto entry = std::find(m_memPool.begin(), m_memPool.end(), pointer);
		if (entry != m_memPool.end())
		{
			throw std::logic_error("double free returned to allocator pool");
		}
		throw std::logic_error("pointer does not belong to allocator pool");
	}

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

	inline long liveCount() const
	{
		return m_pool->pos();
	}

};



template <typename T = double>
class AllAllocators
{
	static int lastSize_N;
	static AllocPolicy<T>* pAllocPolicy;

	// Heap-allocated so the map outlives AllAllocatorsGuard destructors in
	// other TUs. Static std::unordered_map members are destroyed in an
	// unspecified order relative to those guards and caused heap-use-after-free
	// / segfault after gtest had already reported all tests passed.
	static std::unordered_map<int, AllocPolicy<T>*>& policies()
	{
		static auto* map = new std::unordered_map<int, AllocPolicy<T>*>();
		return *map;
	}

	// Heap allocation deliberately keeps the mutex alive through process
	// teardown, including AllAllocatorsGuard destructors in other TUs.
	static std::mutex& registryMutex()
	{
		static auto* mutex = new std::mutex();
		return *mutex;
	}


	static AllocPolicy<T>* setUpPolicy(int size_N)
	{
		auto& map = policies();
		auto itr = map.find(size_N);
		if (map.end() != itr)
		{
			return itr->second;
		}

		// Do not publish the cache pointer until map insertion succeeds. If
		// allocation or insertion throws, unique_ptr releases the new policy
		// and the registry/cache remain unchanged.
		auto newPolicy = std::make_unique<AllocPolicy<T>>(size_N);
		auto inserted = map.emplace(size_N, newPolicy.get());
		if (inserted.second)
		{
			newPolicy.release();
		}
		return inserted.first->second;
	}



public:

	static 	void removePolicy(int size_N)
	{
		std::lock_guard<std::mutex> lock(registryMutex());
		auto& map = policies();
		auto itr = map.find(size_N);
		if (map.end() != itr)
		{
			auto policyPtr = itr->second;
			if (policyPtr->liveCount() != 0)
			{
				throw std::logic_error("cannot remove allocator pool while blocks are live");
			}
			if (pAllocPolicy == policyPtr)
			{
				pAllocPolicy = nullptr;
				lastSize_N = -1;
			}
			delete policyPtr;
			map.erase(itr);
		}
		
	}

	static 	void freeAll()
	{
		std::lock_guard<std::mutex> lock(registryMutex());
		auto& map = policies();
		for (const auto& item : map)
		{
			if (item.second->liveCount() != 0)
			{
				throw std::logic_error("cannot clean allocator pools while blocks are live");
			}
		}
		for (auto& item : map)
		{
			delete item.second;
		}
		map.clear();
		pAllocPolicy = nullptr;
		lastSize_N = -1;
	}


	static T* alloc(int size_N)
	{
		std::lock_guard<std::mutex> lock(registryMutex());
		if (lastSize_N == size_N)
		{
			return  pAllocPolicy->alloc();
		}

		pAllocPolicy = setUpPolicy(size_N);
		lastSize_N = size_N;
		return pAllocPolicy->alloc();
	}



	static void  free(size_t size_N, T* pMem)
	{
		std::lock_guard<std::mutex> lock(registryMutex());
		if (nullptr == pMem)
		{
			return;
		}
		int sz_N = static_cast<int>(size_N);
		auto& map = policies();
		auto policy = map.find(sz_N);
		if (map.end() == policy)
		{
			if (dr3AllocatorDiagnosticsEnabled())
			{
				throw std::logic_error("no allocator pool exists for returned block size");
			}
			return;
		}

		pAllocPolicy = policy->second;
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
void freeAllAllocators(long double);
void freeAllAllocators(double);
void freeAllAllocators(float);
void freeAllAllocators(unsigned int);


template <typename T = double>
class AllAllocatorsGuard
{
public:
	~AllAllocatorsGuard() noexcept
	{
		// Explicit cleanup reports live blocks to callers, but static teardown
		// must never allow that diagnostic to escape a destructor. The registry
		// and mutex intentionally have process lifetime, so leaving a live pool
		// for the operating system to reclaim is safe at process exit.
		try
		{
			freeAllAllocators(T());
		}
		catch (const std::logic_error&)
		{
		}
	}

};
