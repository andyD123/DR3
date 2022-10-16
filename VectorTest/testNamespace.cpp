#include "testNamespace.h"

Numeric asNumber(double x)
{
	return static_cast<Numeric>(x);
}

Numeric asNumber(float x)
{
	return static_cast<Numeric>(x);
}


Numeric asNumber(int x)
{
	return static_cast<Numeric>(x);
}


void EXPECT_NUMERIC_EQ(double x, double y)
{
	EXPECT_DOUBLE_EQ(x, y);
}


void EXPECT_NUMERIC_EQ(float x, float y)
{
	EXPECT_FLOAT_EQ(x, y);
}


void EXPECT_NUMERIC_EQ(int x, int y)
{
	EXPECT_EQ(x, y);
}

