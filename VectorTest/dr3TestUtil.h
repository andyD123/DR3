#pragma once
#include "pch.h"


#include "../Vectorisation/VecX/vec.h"
#include "testNamespace.h"


//using Numeric = InstructionTraits<VecXX::INS>::FloatType;

Numeric asNumber(double x);

Numeric asNumber(float x);

Numeric asNumber(int x);

void EXPECT_NUMERIC_EQ(double x, double y);

void EXPECT_NUMERIC_EQ(float x, float y);

