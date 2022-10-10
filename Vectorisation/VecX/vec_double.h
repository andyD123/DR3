/****************************  vec_double.h *******************************
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
#include <math.h>
#include <cmath>

#include "vec_bool_d.h"


// https://stackedboxes.org/2017/05/01/acklams-normal-quantile-function/

static double cdfnorminv( double p )
{
	/// acklams inverse cdf normal
	static double a[] = { 0.0,  -3.969683028665376e+01, 2.209460984245205e+02,-2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00 };
	static double b[] = { 0.0, -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,  6.680131188771972e+01, -1.328068155288572e+01};
	static double c[] = { 0.0,-7.784894002430293e-03,-3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
	static double d[] = { 0.0,  7.784695709041462e-03, 3.224671290700398e-01,  2.445134137142996e+00, 3.754408661907416e+00 };


	const double p_low  = 0.02425;
	const double p_high = 1 - p_low;


	//Rational approximation for central region.
	if ( (p_low <= p) && (p <= p_high) )
	{
	  double q = p - 0.5;
	  double r = q*q;
	  double x = (((((a[1]*r+a[2])*r+a[3])*r+a[4])*r+a[5])*r+a[6])*q /
			(((((b[1]*r+b[2])*r+b[3])*r+b[4])*r+b[5])*r+1);
	  return x;
	}


	//Rational approximation for lower region.
	if ( ( 0.0 <= p) &&( p < p_low) )
	{
	   double q = sqrt(-2.0*log(p));
	   double x = (((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6]) /
			 ((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1.0);

	   return x;
	}

	//Rational approximation for upper region.
	if ( (p_high < p) && (p <= 1 ) )
	{
	  double q = sqrt(-2.0*log(1-p));
	  double x = -(((((c[1]*q+c[2])*q+c[3])*q+c[4])*q+c[5])*q+c[6]) /
			  ((((d[1]*q+d[2])*q+d[3])*q+d[4])*q+1.0);

	 return x;
	}

//	throw std::range_error("input out of range");
	return 0.0;
}

////

// https://en.cppreference.com/w/cpp/numeric/math/erfc
static double cdfnorm(double x)
{
	return std::erfc(-x/std::sqrt(2.0))/2.0;
}



static double cdfnormD(double x)
{

//   https://mathworld.wolfram.com/Erfc.html
	constexpr  double invRootPi = 0.564189583547756;
	constexpr double invRootTwo = 0.707106781186548;
	return invRootTwo *invRootPi*(exp(-0.5*x*x));
}




class VecDouble
{
protected:

	 double m_data[2];
public:
	VecDouble() {
		m_data[0] = 0.0; m_data[1] = 0.0;
	}

	VecDouble( double d){ m_data[0] =d;	m_data[1] =d;}
	VecDouble( double d0,double d1){ m_data[0] =d0;	m_data[1] =d1;}

	VecDouble(double const   d[2] ){ m_data[0] =d[0];	m_data[1] =d[1];}


	VecDouble( const VecBoolD& rhs)
	{
		m_data[0] = rhs[0];
		m_data[1] = rhs[1];
	}

	operator	VecBoolD() const
	{
		return	VecBoolD(m_data[0],m_data[1]);
	}

	VecDouble& operator = ( const VecDouble& rhs)
	{
		 m_data[0] =rhs.m_data[0];	m_data[1] =rhs.m_data[1];
		 return *this;
	}


	VecDouble& load_a ( const double* p)
	{
		 m_data[0] =p[0];
		 m_data[1] =p[1];
		 return *this;
	}


	VecDouble& load(const double* p)
	{
		m_data[0] = p[0];
		m_data[1] = p[1];
		return *this;
	}

   void store_a( double* p)
	{
		p[0] = m_data[0];
		p[1] =m_data[1];
	}

   void store(double* p)
   {
	   p[0] = m_data[0];
	   p[1] = m_data[1];
   }

   double extract(size_t idx) const
   {
	    return m_data[idx];
   }


   void insert(size_t idx,double val)
   {
	  m_data[idx]= val;
   }


    double operator [] (size_t index) const
	{
        return extract(index);
	}

   static int size() { return 2;}

	
};



// vector operator + : add element by element
static inline VecDouble operator + (VecDouble const & a, VecDouble const & b) {
	return VecDouble(a[0]+b[0],a[1]+b[1]);
}

// vector operator + : add vector and scalar
static inline VecDouble operator + (VecDouble const & a, double b) {
    return a + VecDouble(b);
}
static inline VecDouble operator + (double a, VecDouble const & b) {
    return VecDouble(a) + b;
}


static inline VecDouble & operator += (VecDouble & a, VecDouble const & b)
{
    a = a + b;
    return a;
}


static inline VecDouble operator ++ (VecDouble & a, int) {
    VecDouble a0 = a;
    a = a + 1.0;
    return a0;
}


static inline VecDouble & operator ++ (VecDouble & a) {
    a = a + 1.0;
    return a;
}


static inline VecDouble operator - (VecDouble const & a, VecDouble const & b)
{
    return  VecDouble(a[0]-b[0],a[1]-b[1]);
}


static inline VecDouble operator - (VecDouble const & a, double b) {
    return a - VecDouble(b);
}
static inline VecDouble operator - (double a, VecDouble const & b) {
    return VecDouble(a) - b;
}


static inline VecDouble operator - (VecDouble const & a)
{
	return VecDouble(-a[0],-a[1]);
}


static inline VecDouble & operator -= (VecDouble & a, VecDouble const & b) {
    a = a - b;
    return a;
}


static inline VecDouble operator -- (VecDouble & a, int) {
    VecDouble a0 = a;
    a = a - 1.0;
    return a0;
}


static inline VecDouble & operator -- (VecDouble & a)
{
    a = a - 1.0;
    return a;
}


static inline VecDouble operator * (VecDouble const & a, VecDouble const & b)
{
    return   VecDouble(a[0]*b[0],a[1]*b[1]);
}


static inline VecDouble operator * (VecDouble const & a, double b) {
    return a * VecDouble(b);
}
static inline VecDouble operator * (double a, VecDouble const & b) {
    return VecDouble(a) * b;
}


static inline VecDouble & operator *= (VecDouble & a, VecDouble const & b)
{
    a = a * b;
    return a;
}


static inline VecDouble operator / (VecDouble const & a, VecDouble const & b)
{
    return VecDouble(a[0]/b[0],a[1]/b[1]);
}


static inline VecDouble operator / (VecDouble const & a, double b)
{
    return a / VecDouble(b);
}
static inline VecDouble operator / (double a, VecDouble const & b)
{
    return VecDouble(a) / b;
}


static inline VecDouble & operator /= (VecDouble & a, VecDouble const & b)
{
    a = a / b;
    return a;
}



static inline VecBoolD operator == (VecDouble const & a, VecDouble const & b) 
{
    return VecBoolD(a[0] == b[0], a[1] ==b[1]);
}


static inline VecBoolD operator != (VecDouble const & a, VecDouble const & b)
{
    return VecBoolD(a[0] != b[0], a[1] !=b[1]);
}


static inline VecBoolD operator < (VecDouble const & a, VecDouble const & b)
{
    return VecBoolD(a[0] < b[0], a[1] < b[1]);
}


static inline VecBoolD operator <= (VecDouble const & a, VecDouble const & b) {
     return VecBoolD(a[0] <= b[0], a[1] <= b[1]);
}


static inline VecBoolD operator > (VecDouble const & a, VecDouble const & b)
{
    return b < a;
}


static inline VecBoolD operator >= (VecDouble const & a, VecDouble const & b) 
{
    return b <= a;
}


static inline VecDouble operator & (const VecDouble  & a, const VecDouble  & b)
{
    return  VecDouble(long(a[0]) & long( b[0]) , long( a[1])& long(b[1]) );
}


static inline VecDouble & operator &= (VecDouble & a, VecDouble const & b) {
    a = a & b;
    return a;
}


static inline VecDouble operator & (VecDouble const & a, VecBoolD const & b)
{
	 return VecDouble(long(a[0])&long(b[0]), long(a[1])&long(b[1]) );
}

static inline VecDouble operator & (VecBoolD const & a, VecDouble const & b)
{
	return VecDouble(long(a[0])&long(b[0]), long(a[1]) &long(b[1]) );
}


static inline VecDouble operator | (VecDouble const & a, VecDouble const & b)
{
	return VecDouble(long(a[0]) | long( b[0]), long( a[1]) | long(b[1]) );
}


static inline VecDouble & operator |= (VecDouble & a, VecDouble const & b) {
    a = a | b;
    return a;
}


static inline VecDouble operator ^ (VecDouble const & a, VecDouble const & b)
{
	return VecDouble(long(a[0]) ^long( b[0]), long(a[1]) ^long( b[1] ));
}


static inline VecDouble & operator ^= (VecDouble & a, VecDouble const & b)
{
    a = a ^ b;
    return a;
}


static inline VecDouble operator ! (VecDouble const & a)
{
	return VecDouble( asDouble( a[0]!=0.0), asDouble( a[1]!=0.0) );
}


static inline VecDouble select ( VecBoolD const & s,  VecDouble const & a,  VecDouble const & b)
{
	 return VecDouble( s[0]?a[0]:b[0], s[1]?a[1]:b[1] );
}


static inline double horizontal_add (VecDouble const & a)
{
	return a[0]+ a[1];
}


static inline VecDouble max(VecDouble const & a, VecDouble const & b)
{
	return VecDouble(a[0] >b[0]? a[0] :b[0] ,a[1] >b[1]? a[1] :b[1] );
}



static inline VecDouble min(VecDouble const & a, VecDouble const & b)
{
	return VecDouble(a[0] <b[0]? a[0] :b[0] ,a[1]<b[1]? a[1] :b[1] );
}


static inline VecDouble mul_add(VecDouble const & a, VecDouble const & b, VecDouble const & c )
{
	return VecDouble(a[0]* b[0] +c[0], a[1] * b[1] +c[1]);
}

static inline VecDouble mul_add(VecDouble const & a, VecDouble const & b, double const & c)
{
	return VecDouble(a[0] * b[0] + c, a[1] * b[1] + c);
}

static inline VecDouble mul_add(VecDouble const & a, double const & b, double const & c)
{
	return VecDouble(a[0] * b + c, a[1] * b + c);
}


static inline VecDouble abs(VecDouble const & a) 
{
	return VecDouble(::fabs(a[0]) , ::fabs(a[1] ) );
}


static inline VecDouble sqrt(VecDouble const & a) 
{
	return VecDouble(::sqrt((a[0])) , ::sqrt(a[1] ) );
}


static inline VecDouble square(VecDouble const & a) 
{
	return VecDouble(a[0] * a[0] , a[1] * a[1]  );
}


static inline VecDouble pow(VecDouble const & a,double b) 
{
	return VecDouble( ::pow(a[0],b) , ::pow(a[1],b)  );
}


static inline VecDouble pow(VecDouble const & a, VecDouble const & b)
{
	return VecDouble(::pow(a[0], b[0]), ::pow(a[1], b[1]));
}

static inline VecDouble exp(VecDouble const & a) 
{
	return VecDouble( ::exp(a[0]) , ::exp(a[1])  );
}


static inline VecDouble log(VecDouble const & a) 
{
	return VecDouble( ::log(a[0]) , ::log(a[1])  );
}


static inline VecDouble ceil(VecDouble const & a) 
{
	return VecDouble( ::ceil(a[0]) , ::ceil(a[1])  );
}


static inline VecDouble floor(VecDouble const & a) 
{
	return VecDouble( ::floor(a[0]) , ::floor(a[1])  );
}


static inline VecDouble sin(VecDouble const& a)
{
	return VecDouble(::sin(a[0]), ::sin(a[1]));
}

static inline VecDouble cos(VecDouble const& a)
{
	return VecDouble(::cos(a[0]), ::cos(a[1]));
}



static inline VecDouble cdfnorm(VecDouble const & a) 
{
	return VecDouble( ::cdfnorm(a[0]) , ::cdfnorm(a[1])  );
}


static inline VecDouble cdfnormD(VecDouble const & a)
{
	return VecDouble(::cdfnormD(a[0]), ::cdfnormD(a[1]));
}


static inline VecDouble cdfnorminv(VecDouble const & a) 
{
	return VecDouble( ::cdfnorminv(a[0]) , ::cdfnorminv(a[1])  );
}
