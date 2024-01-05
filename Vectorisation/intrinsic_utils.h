#pragma once
/****************************  unroll_operators.h   *******************************
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


#include "../Vectorisation/VecX/instruction_traits.h"
#include <immintrin.h>

/*
//_MM_MANTISSA_NORM_ENUM n;
//_MM_MANTISSA_SIGN_ENUM m;

static inline __m256d  getMantissa(__m256d x)// const int n, const int m) 
{

	return _mm256_getmant_pd(x, _MM_MANT_NORM_1_2, _MM_MANT_SIGN_zero);
}

static inline __m512d  getMantissa(__m512d x, int n, const int m)
{
	return _mm512_getmant_pd(x,n,m);
}



static inline __m256d  getExponent(__m256d x)
{

	return _mm256_getexp_pd(x);
}

static inline __m512d  getExponent(__m512d x)
{
	return _mm512_getexp_pd(x);
}

static inline __m128d  getExponent(__m128d x)
{
	return _mm_getexp_pd(x);
}

static inline __m256  getExponent(__m256 x)
{

	return _mm256_getexp_ps(x);
}

static inline __m512  getExponent(__m512 x)
{
	return _mm512_getexp_ps(x);
}

static inline __m128  getExponent(__m128 x)
{
	return _mm_getexp_ps(x);
}

*/