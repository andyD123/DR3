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