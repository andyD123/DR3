/****************************  error_utils.h   *******************************
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
#include "vec.h"
#include "vec_view.h"
#include <exception>
#include <cassert>

//ignore for unused from sutter
template<class T> void ignore(const T&) { }

template<typename VEC>
bool check_vector( const VEC& rhs)
{
	auto rhsSz = rhs.size();

	if ( ( rhsSz > 0) || rhs.isScalar() )
	{
		return true;		
	}
	else
	{
	    //std::
		assert(false);
		throw std::exception("bad vector size of non scalar");
	}
}

template<typename VEC>
bool check_pair(const VEC& lhs, const VEC& rhs)
{
	check_vector(lhs);
	check_vector(rhs);

	if ( (lhs.size() == rhs.size() ) && (lhs.size() > 0 ) )
	{
		return true;
	}
	

	if (rhs.isScalar() || lhs.isScalar())
	{
		return true;
	}
	else
	{
		
		assert(false);
		throw std::exception("bad vector size");
	}
}

template<typename VEC1, typename VEC2>
bool check_pair_different_type(const VEC1& lhs, const VEC2& rhs)
{
	check_vector(lhs);
	check_vector(rhs);

	if (lhs.size() == rhs.size())
		return true;
	if (rhs.isScalar() || lhs.isScalar())
	{
		return true;
	}
	else
	{
		//std::
		assert(false);
		throw std::exception("bad vector size");
	}
}




//////////////  views ////////////////
template<typename INS_VEC>
bool check_vector(const VecView<INS_VEC>& /*rhs*/)
{
	//TO DO
	/*
	auto rhsSz = rhs.size();

	if ((rhsSz > 0) || rhs.isScalar())
	{
		return true;
	}
	else
	{
		//std::assert(false);
		throw std::exception("bad vector size of non scalar");
	}
	*/
	return true;
}

template<typename INS_VEC>
bool check_vector(const Vec<INS_VEC>& rhs)
{
	auto rhsSz = rhs.size();

	if ((rhsSz > 0) || rhs.isScalar())
	{
		return true;
	}
	else
	{
		//std::
		assert(false);
		throw std::exception("bad vector size of non scalar");
	}
}



template<typename INS_VEC>
bool check_vector(const VecD<INS_VEC>& rhs)
{
	//do nothing
	/*
	auto rhsSz = rhs.size();

	if ((rhsSz > 0) || rhs.isScalar())
	{
		return true;
	}
	else
	{
		//std::
		assert(false);
		throw std::exception("bad vector size of non scalar");
	}
	*/
}

template<typename INS_VEC>
bool check_vector_for_filter(const Vec<INS_VEC>& rhs)
{
	auto rhsSz = rhs.size();

	if ((rhsSz > 0) || !rhs.isScalar()) //no scalar vectors for filtering to views
	{
		return true;
	}
	else
	{
		//std::
		assert(false);
		throw std::exception("bad vector size of non scalar");
	}
}


template<typename INS_VEC>
bool check_vector_for_filter(const VecView<INS_VEC>&/* rhs*/)
{
	return true;// views can be empty
	/*
	auto rhsSz = rhs.size();

	if ((rhsSz > 0) || !rhs.isScalar()) //no scalar vectors for filtering to views
	{
		return true;
	}
	else
	{
		//std::assert(false);
		throw std::exception("bad vector size of non scalar");
	}
	*/
}


template<typename INS_VEC>
bool check_view_pair(const Vec<INS_VEC>& lhs, const Vec<INS_VEC>& rhs)
{
	check_vector_for_filter(lhs);
	check_vector_for_filter(rhs);

	if (lhs.size() == rhs.size())
		return true;

	//std::
	assert(false);
	throw std::exception("bad vector size");

}

template<typename INS_VEC>
bool check_view_pair(const VecView<INS_VEC>& lhs, const VecView<INS_VEC>& rhs)
{
	check_vector_for_filter(lhs);
	check_vector_for_filter(rhs);

	if (lhs.size() == rhs.size())
		return true;

	//std::
	assert(false);
	throw std::exception("bad vector size");

}