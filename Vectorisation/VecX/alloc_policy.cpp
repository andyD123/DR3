/****************************  alloc_policy.cpp   *******************************
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
#include "alloc_policy.h"
#include "alloc_policy_imp.h"

bool dr3AllocatorDiagnosticsEnabled() noexcept
{
#ifdef NDEBUG
	return false;
#else
	return true;
#endif
}

template<>
int AllAllocators<long double>::lastSize_N = -1;
template<>
int AllAllocators<double>::lastSize_N = -1;
template<>
int AllAllocators<float>::lastSize_N = -1;
template<>
int AllAllocators<unsigned int>::lastSize_N = -1;

template<>
AllocPolicy<long double>* AllAllocators<long double>::pAllocPolicy = nullptr;
template<>
AllocPolicy<double>* AllAllocators<double>::pAllocPolicy = nullptr;
template<>
AllocPolicy<float>* AllAllocators<float>::pAllocPolicy = nullptr;
template<>
AllocPolicy<unsigned int>* AllAllocators<unsigned int>::pAllocPolicy = nullptr;



void freePool(size_t N, long double* pOld)
{
	return freeT(N, pOld);
}


void freePool(size_t N, double* pOld)
{
	return freeT(N, pOld);
}

void freePool(size_t N, float* pOld)
{
	return freeT(N, pOld);
}

void freePool(size_t N, unsigned int* pOld)
{
	return freeT(N, pOld);
}

void allocPool(size_t& N, long double*& pMem)
{
	allocT(N, pMem);
}

void allocPool(size_t& N, double*& pMem)
{
	allocT(N, pMem);
}

void allocPool(size_t& N, float*& pMem)
{
	allocT(N, pMem);
}

void allocPool(size_t& N, unsigned int*& pMem)
{
	allocT(N, pMem);
}

int  getAllignedSize(size_t N, long double* pOld)
{
	return getAllignedSizeT(N, pOld);
}

int  getAllignedSize(size_t N, double* pOld)
{
	return getAllignedSizeT(N, pOld);
}

int  getAllignedSize(size_t N, float* pOld)
{
	return getAllignedSizeT(N, pOld);
}

int  getAllignedSize(size_t N, unsigned int* pOld)
{
	return getAllignedSizeT(N, pOld);
}
void freeAllAllocators(long double)
{
	AllAllocators<long double>::freeAll();
}
void freeAllAllocators(double)
{
	AllAllocators<double>::freeAll();
}
void freeAllAllocators(float)
{
	AllAllocators<float>::freeAll();
}
void freeAllAllocators(unsigned int)
{
	AllAllocators<unsigned int>::freeAll();
}
