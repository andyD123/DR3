// dancingAVX512.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <algorithm>
#include <random>
#include <numeric>
#include <iterator>
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>  
#include <thread>
#include <map>
#include <cstring>


#include "../Vectorisation/VecX/dr3.h"
#include "../Vectorisation/VecX/accumulate_transform.h"
#include "../Vectorisation/VecX/error_utils.h"


#include "../Vectorisation/VecX/zip_utils.h"

//#include "norm.h"

#include "AVX512Dance.h"

// use namespace DRC::VecD8D  run this and watch power consumption
// switches between AVX2 and AVX512  implementations
// AVX512 uses less energy in this case

int main()
{
  //  std::cout << "Hello World!\n";
    doAVXMax512Dance();
}

