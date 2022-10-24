/****************************  alloc_policy.h   *******************************
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

void freePool(size_t N, double* pOld);
void freePool(size_t N, float* pOld);
void freePool(size_t N, unsigned int* pOld);


void allocPool(size_t& N, double*& pMem);
void allocPool(size_t& N, float*& pOld);
void allocPool(size_t& N, unsigned int*& pOld);


int  getAllignedSize(size_t N, double* pOld);
int  getAllignedSize(size_t N, float* pOld);
int  getAllignedSize(size_t N, unsigned int* pOld);









