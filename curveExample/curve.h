#pragma once
#include <vector>
#include <algorithm>
#include <exception>
#include <deque>
#include <unordered_map>
#include <list>
#include <cstddef>
#include <stdexcept>

//#include "ApplyOperation.h"

template<typename T, typename V>
struct Extrap_end_flat_last
{
	static V extrap(const T& , const std::vector<V>& vals, const  std::vector<T>&  ords)
	{
		if(! ords.empty()   && !vals.empty() )
		{
			return vals.back();
		}
		else
		{
			throw std::exception( " T or V vectors  empty in extrap");
		}
	}
};

template<typename T, typename V>
struct Extrap_start_flat_first
{
	static V extrap(const T& , const std::vector<V>& vals, const  std::vector<T>&  ords)
	{
		
		if(! ords.empty()   && !vals.empty() )
		{
			return vals.front();
		}
		else
		{
			throw std::exception( " T or V vectors  empty in extrap");
		}
	}
};


template<typename T, typename V>
struct LinearInterp
{
	static V calc(const T& dt, const std::vector<V>& vals, const  std::vector<T>&  ords)
	{
		if( (!ords.empty())   &&   (ords.size() == vals.size() )  )
		{
			auto loBound = std::lower_bound(cbegin(ords),cend(ords),dt);
			if( dt != *loBound)
			{
				--loBound;
			}
			auto pos = std::distance(cbegin(ords), loBound);

			//auto fraction = (dt-ords[pos])/static_cast<V>(ords[pos+1] - ords[pos]);
			auto fraction = (dt - ords[pos])*1.0 /( (ords[pos + 1] - ords[pos]) );
			auto& lowerVal = vals[pos];
			auto& upperVal = vals[pos+1];

			return lowerVal + (upperVal -lowerVal) * fraction;

		}
		else
		{
			throw std::exception( " T or V vectors  empty in extrap");
		}
	}

};



template<typename T, typename V>
struct FlatInterp
{
	static V calc(const T& dt, const std::vector<V>& vals, const  std::vector<T>&  ords)
	{
		if ((!ords.empty()) && (ords.size() == vals.size()))
		{
			auto loBound = std::lower_bound(cbegin(ords), cend(ords), dt);
			if (dt != *loBound)
			{
				--loBound;
			}
			auto pos = std::distance(cbegin(ords), loBound);
			return vals[pos];
		}
		else
		{
			throw std::exception(" T or V vectors  empty in extrap");
		}
	}

};




template <typename T , typename V, 
	typename INTERP_POLICY =  LinearInterp<T,V>,
	typename  EXTRAP_END_POLICY =  Extrap_end_flat_last<T,V>,
	typename EXTRAP_START_POLICY  =   Extrap_start_flat_first<T,V> 
	>
class Curve
{
public:
	//to do non copy etc

	Curve(){}
	~Curve(){}

	V  valueAt( T dt)
	{
		if ( ( dt >= m_min) && ( dt <= m_max) )
		{
			return INTERP_POLICY::calc(dt, m_vals, m_ords);
		}

		if( dt < m_min)
		{
			return EXTRAP_START_POLICY::extrap(dt, m_vals, m_ords);
		}
		else //if( dt > m_max)
		{
			return EXTRAP_END_POLICY::extrap(dt, m_vals, m_ords);
		}

		
	}

	template< typename IT_T,typename IT_V>
	void setValues(  IT_T start_T,  IT_T last_T,  IT_V start_V,  IT_V last_V)
	{
		auto T_size = std::distance(start_T, last_T);
		auto V_size = std::distance(start_V, last_V);

		if ( T_size != V_size)
		{
			throw std::exception( " T & V vectors of different sizes");
		}

		m_vals.clear();
		m_vals.resize(V_size);

		std::copy(start_V,last_V, m_vals.begin() );

		m_ords.clear();
		m_ords.resize(T_size);

		std::copy(start_T,last_T, m_ords.begin() );

		m_min = *start_T;
		m_max = *(--last_T);
		//assumes we put in ordered unique  ordinates

	}


private:
	T m_min = 0;
    T m_max = 0;

	std::vector<V> m_vals;
	std::vector<T> m_ords;

};



template<typename T, typename V>
struct ZeroInterp : protected LinearInterp<T,V>
{
	static V calc(const T& dt, const std::vector<V>& vals, const  std::vector<T>&  ords)
	{

		const auto& Yields = LinearInterp<T, V>::calc(dt, vals, ords);

		return exp(-Yields * dt);

	}
};








/// 
//git hub cache
////

template<typename key_t, typename value_t>//, typename CALC = []() {throw std::range_error("There is no such key in cache")} >
class lru_cache {
public:
	typedef typename std::pair<key_t, value_t> key_value_pair_t;
	typedef typename std::list<key_value_pair_t>::iterator list_iterator_t;

	lru_cache(size_t max_size) :
		_max_size(max_size) {
	}

	void put(const key_t& key, const value_t& value) {
		auto it = _cache_items_map.find(key);
		_cache_items_list.push_front(key_value_pair_t(key, value));
		if (it != _cache_items_map.end()) {
			_cache_items_list.erase(it->second);
			_cache_items_map.erase(it);
		}
		_cache_items_map[key] = _cache_items_list.begin();

		if (_cache_items_map.size() > _max_size) {
			auto last = _cache_items_list.end();
			last--;
			_cache_items_map.erase(last->first);
			_cache_items_list.pop_back();
		}
	}

	const value_t& get(const key_t& key)//, CALC& calc) 
	{
		auto it = _cache_items_map.find(key);
		if (it == _cache_items_map.end()) {
			throw std::range_error("There is no such key in cache");
			//auto calcRes = calc(key);
			//put(key, calcRes);
			//return calcRes;
		}
		else {
			_cache_items_list.splice(_cache_items_list.begin(), _cache_items_list, it->second);
			return it->second->second;
		}
	}

	bool exists(const key_t& key) const {
		return _cache_items_map.find(key) != _cache_items_map.end();
	}

	size_t size() const {
		return _cache_items_map.size();
	}

private:
	std::list<key_value_pair_t> _cache_items_list;
	std::unordered_map<key_t, list_iterator_t> _cache_items_map;
	size_t _max_size;
};



template <typename T, typename V,
	typename INTERP_POLICY = LinearInterp<T, V>,
	typename  EXTRAP_END_POLICY = Extrap_end_flat_last<T, V>,
	typename EXTRAP_START_POLICY = Extrap_start_flat_first<T, V>
>
class Curve2
{
public:
	//to do non copy etc

	Curve2(size_t cacheSz): m_min(0),m_max(0), m_interpCache(cacheSz) {}

	~Curve2() {}

	V  valueAt(T dt)
	{
		if ((dt >= m_min) && (dt <= m_max))
		{

			return  m_interpCache.calc(dt, m_vals, m_ords);
			//return INTERP_POLICY::calc(dt, m_vals, m_ords);
		}

		if (dt < m_min)
		{
			return EXTRAP_START_POLICY::extrap(dt, m_vals, m_ords);
		}
		else //if( dt > m_max)
		{
			return EXTRAP_END_POLICY::extrap(dt, m_vals, m_ords);
		}


	}

	template< typename IT_T, typename IT_V>
	void setValues(IT_T start_T, IT_T last_T, IT_V start_V, IT_V last_V)
	{
		auto T_size = std::distance(start_T, last_T);
		auto V_size = std::distance(start_V, last_V);

		if (T_size != V_size)
		{
			throw std::exception(" T & V vectors of different sizes");
		}

		m_vals.clear();
		m_vals.resize(V_size);

		std::copy(start_V, last_V, m_vals.begin());

		m_ords.clear();
		m_ords.resize(T_size);

		std::copy(start_T, last_T, m_ords.begin());

		m_min = *start_T;
		m_max = *(--last_T);
		//assumes we put in ordered unique  ordinates

	}


private:
	T m_min;
	T m_max;

	std::vector<V> m_vals;
	std::vector<T> m_ords;

	INTERP_POLICY  m_interpCache;

};


template<typename T, typename V>
struct MyCalc
{

	MyCalc() {};

	MyCalc(const std::vector<V>&  vals, const  std::vector<T>&  ords) :
		m_vals(vals), m_ords(ords)
	{
		m_deltaTAtOrds.emplace_back(0);
		V discountfactor = vals[0] / vals[0];// { 1.0, 0 };// InstructionTraits<V>::oneValue;//1.0;  
		V fwdRate;
		m_DiscountFactorAtOrd.emplace_back(discountfactor);
		m_Fwds.emplace_back(vals[0]);
		//m_Fwds.emplace_back(vals[0]);// flat  rate
		//discountfactor = exp(-ords[0] * vals[1]);
		for (size_t i = 1; i < ords.size(); ++i)
		{
			auto deltaT = ords[i] - ords[i - 1];
			m_deltaTAtOrds.emplace_back(deltaT);

			auto newDiscountfactor = exp(-m_ords[i] * vals[i]);

			auto dfOverInterval = newDiscountfactor / m_DiscountFactorAtOrd.back();

			auto Yavg = -log(dfOverInterval) / deltaT;

			m_Fwds.emplace_back(2.0*Yavg - m_Fwds.back());
			m_DiscountFactorAtOrd.emplace_back(newDiscountfactor);
		}

	}

	void init(const std::vector<V>&  vals, const  std::vector<T>&  ords) 
	{
		m_vals=vals;
		m_ords =ords;

		m_deltaTAtOrds.emplace_back(0);
		V discountfactor = vals[0] / vals[0];// { 1.0, 0 };// InstructionTraits<V>::oneValue;//1.0;  
		V fwdRate;
		m_DiscountFactorAtOrd.emplace_back(discountfactor);
		m_Fwds.emplace_back(vals[0]);
		//m_Fwds.emplace_back(vals[0]);
		//m_Fwds.emplace_back(vals[0]);// flat  rate
		//discountfactor = exp(-ords[0] * vals[1]);
		for (size_t i = 1; i < ords.size(); ++i)
		{
			auto deltaT = ords[i] - ords[i - 1];
			m_deltaTAtOrds.emplace_back(deltaT);

			auto newDiscountfactor = exp(-m_ords[i] * vals[i]);

			auto dfOverInterval = newDiscountfactor / m_DiscountFactorAtOrd[i-1];

			auto Yavg = -log(dfOverInterval) / deltaT;

			discountfactor *= dfOverInterval;

			i==1?m_Fwds.emplace_back(vals[1]):m_Fwds.emplace_back(2.0*Yavg - m_Fwds[i-1]);
			m_DiscountFactorAtOrd.emplace_back(discountfactor);// newDiscountfactor);
		}

	}



	V operator()(const T& t)
	{
		//return Calc(t, m_vals, m_ords, m_DiscountFactorAtOrd, m_deltaTAtOrds);
		return Calc(t, m_Fwds, m_ords, m_DiscountFactorAtOrd, m_deltaTAtOrds);
	}

	//V Calc(const T& dt, const std::vector<V>& vals, const  std::vector<T>&  ords)
	//{
	//	const auto& Yields = LinearInterp<T, V>::calc(dt, vals, ords);
	//	return exp(-Yields * dt);
	//}


	V Calc(const T& t, const std::vector<V>& vals, const  std::vector<T>&  ords,
		const  std::vector<V>&  discountFactorAtOrd,
		const  std::vector<T>& )// deltaTAtOrd)
	{
	const auto& Yields = LinearInterp<T, V>::calc(t, vals, ords);
	const auto& DFsAtPoints = FlatInterp<T, V>::calc(t, discountFactorAtOrd, ords);
	//const auto& dt = t-FlatInterp<T, T>::calc(t, deltaTAtOrd, ords);
	const auto& dt = t - FlatInterp<T, T>::calc(t, ords, ords);
	//const auto& interval = FlatInterp<T, T>::calc(deltaTAtOrd, ords, ords);
	//const auto& dt = t - FlatInterp<T, T>::calc(deltaTAtOrd, ords, ords);
	return exp(-Yields * dt)* DFsAtPoints;

	//return exp(-Yields * t); zero rate intep
	}


//	const std::vector<V>&  m_vals;
//	const  std::vector<T>&  m_ords;


	std::vector<V>  m_vals;
	std::vector<T>  m_ords;
	std::vector<V>  m_Fwds;
	std::vector<V>  m_DiscountFactorAtOrd;
	std::vector<T>  m_deltaTAtOrds;
};

template<typename T, typename V>
struct ZeroInterpCached  : protected  LinearInterp<T, V>
{
	using LRU_C = lru_cache<T, V>;

	LRU_C m_lru;
	MyCalc<T, V> clc;// (vals, ords);
	bool isInit = false;


	ZeroInterpCached<T, V>(size_t maxCache) : m_lru(maxCache) {};

	 V calc(const T& dt, const std::vector<V>& vals, const  std::vector<T>&  ords)
	{

		 if (!isInit)
		 {
			 clc.init(vals, ords);
			 isInit = true;
		 }

		 if (m_lru.exists(dt))
		 {
			// V ret( m_lru.get(dt) );
			// return ret;

			 return V(m_lru.get(dt));
		 }
		 else
		 {
			// MyCalc<T, V> clc(vals, ords);
			 //V ret = clc(dt);
			 m_lru.put(dt, clc(dt));
			 return V(m_lru.get(dt));
		 }

		//MyCalc<T, V> clc(vals, ords);
		//exists
		//V ret = m_lru.get(dt, clc);
		// V ret; //not reached
		 //return ret;
	}


};







/*
template <typename IDX, typename VAL>
struct CurveCash
{

	struct CashItem
	{
		IDX  idx;
		VAL* pVal;
	};
	int maxCount;
	int count = 0;
	std::deque<CashItem> cash;
	//std::unordererd_map< IDX,

	void addNewItem(VAL* pVal, IDX idx)
	{
		if (count == max)
		{
			auto item = cash.front();
			cash.pop_front();
			delete item->pVal;
			--count;
		}

		cash.emplace_back(idx, pVal);
		++count;
	}

	void updateItem(VAL* pVal, IDX idx)
	{

		if (count == max)
		{
			auto item = cash.front();
			cash.pop_front();
			delete item->pVal;
			--count;
		}

		cash.emplace_back(idx, pVal);
		++count;
	}

	VAL*  update(IDX idx)
	{
		if (isFound(idx))
		{
			moveToFront()
		}
		else
		{
			addNewItem(pVal, idx);
		}
	}
	

};

*/
