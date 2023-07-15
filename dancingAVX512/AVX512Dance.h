#pragma once


void doAVXMax512Dance();

 //return  iff((x - (VecXX::INS(2.0) * floor(x * VecXX::INS(0.5)))) >= VecXX::INS(1.0),			//auto oneIfOddLmbda = [&](auto x) { iff(x > VecXX::INS(0.5 * SZ), 