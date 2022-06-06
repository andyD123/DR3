/****************************  vec_bool_d.h  *******************************
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
#include "vec_double.h"

static inline double asDouble( bool bVal)
{
	if(!bVal) return 0.0;
	return static_cast<double>(0xFFFFFFFFFFFFFFFF);
}


static inline bool asBool( double val)
{
	if( 0.0 == val) return false;
	return true;
	
}


class VecBoolD
{
private:

	 double m_data[2];
public:
	VecBoolD()	{ m_data[0] = asDouble(false);	m_data[1] = asDouble(false);}

	VecBoolD( bool d0,bool d1){ m_data[0] = asDouble(d0);	m_data[1] = asDouble(d1);}

	VecBoolD& operator = ( const VecBoolD& rhs)
	{
		 m_data[0] =rhs.m_data[0];	m_data[1] =rhs.m_data[1];
		 return *this;
	}


	VecBoolD& load_a ( const double* p)
	{
		 m_data[0] =p[0];
		 m_data[1] =p[1];
		 return *this;
	}


   void store_a( double* p)
	{
		p[0] = m_data[0];
		p[1] =m_data[1];
	}

   bool extract(size_t idx) const
   {
	    return  asBool(m_data[idx]);
   }


   void insert(size_t idx,bool val)
   {
	  m_data[idx]= asDouble(val);
   }


    bool operator [] (size_t index) const
	{
        return extract(index);
	}

   static int size() { return 2;}

   inline bool isScalar() const
   {
	   return false;
   }

	
};



static inline bool horizontal_or(VecBoolD const & a)
{
	return a[0] || a[1];
}

static inline bool horizontal_and(VecBoolD const & a)
{
	return a[0] && a[1];
}


static inline VecBoolD operator  &&(VecBoolD const& a, const VecBoolD& b)
{
	return  VecBoolD( a[0] & b[1],a[1] & b[1]);
}


static inline VecBoolD operator  ||(VecBoolD const& a, const VecBoolD& b)
{
	return  VecBoolD(a[0] | b[1], a[1] | b[1]);
}


static inline VecBoolD operator !(VecBoolD const& a)
{
	return  VecBoolD(!a[0] , !a[1]);
}
