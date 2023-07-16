// GettingStarted.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <algorithm>
#include <random>
#include <numeric>
#include <iterator>
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>  


#include "../Vectorisation/VecX/dr3.h"



//using namespace DRC::VecDb;
//using namespace DRC::VecD2D;
//using namespace DRC::VecD4D;
using namespace DRC::VecD8D;
//using namespace DRC::VecF16F;
//using namespace DRC::VecF8F;

//#include <iostream>

void basicVector();
void applyLambdasToVector();
void accumulateVector();
void applyBinaryLambdasToVector();
void applySparseLambdasToVector();
void transformReduce();
void branching();
void filtersAndViews();
void complexExampleOfJoiningfiltersAndViews();
void experimentalAAD();

void optionPricerDriver();
double blackScholes_1(double S, double K, double t, double r, double sigma);
double  blackScholes_2(double S, double K, double t, double r, double sigma);
VecXX  blackScholes_3(const VecXX& S, const VecXX& K, const VecXX& t, const VecXX& r, const VecXX& sigma);



AllAllocatorsGuard<typename VecXX::SCALA_TYPE> allocGuard;



int main()
{
   basicVector();

   applyLambdasToVector();

   applyBinaryLambdasToVector();
   
   applySparseLambdasToVector();

   accumulateVector();

   transformReduce();

   branching();

   filtersAndViews();

   complexExampleOfJoiningfiltersAndViews();

   experimentalAAD();

   optionPricerDriver();
}


VecXX getVector(int SZ=100)
{
    std::vector<double>  stlVec(SZ);
    std::iota(begin(stlVec), end(stlVec), 0.0);
    //creates a vector from STL
    VecXX vec1 = stlVec;
    return vec1;
}






void basicVector()
{
    std::vector<double>  stlVec(100);
    std::iota(begin(stlVec), end(stlVec), 0.0);


    //create a vector from STL
    VecXX vec1 = stlVec;

    //create vecXX from a sum
    VecXX vec2 = vec1 + 2.0;

    //create an STL vector from the VexXX
    std::vector<double> outVec = vec2; // conversion   

    for (auto& yy : vec2)
    {
        yy = 21.1;
    }

    //do some element wise maths on the vectors
    vec2 = vec1 * 3.0; //mult scalar

    vec2 *= 2.0;
    vec2 += 0.1;
    auto vec3 = sqrt(vec2);  // sqrt function , all elements in vec3[i] = sqrt(vec2[i]) 

    
    //with vectors too
    vec2 *= vec2;
    vec2 += vec3;


    std::vector<double> outVec1 = vec3; // a bit easier  to see the elements (std::vector<double>)

    //create a vecXX from a scalar 
    VecXX four = 4.0; // a scalar VecXX

    vec2 = vec2 + four;
    vec2 = four + vec2;

    bool testScalar = four.isScalar();
    bool testScalar2 = vec2.isScalar();

   // do some vector vector maths
    auto elementWiseSum = vec2 + vec3;

    auto elemWiseMult = vec2 * vec3;

    //complex expression
    auto x = vec2;
    auto y = 10.0 -3.2*x + exp(-sqrt(x * x + sin(x)));


    ignore(testScalar);
    ignore(testScalar2);

}




void applyLambdasToVector()
{

    //square all elements 
    auto square = [](auto x) { return x * x; };
    

//  ApplyTransform 
    

    auto  v1 = getVector(); 
    std::vector<double> dbVec = v1;

    auto res = transform(square, v1);
    dbVec = res;

    //transform a constant vector
    const auto  v_const = getVector();
    res = transform(square, v_const);
    dbVec = res;

    // transform modify inplace 
   VecXX updater1 = getVector();
   transformM(square, updater1);
   dbVec = updater1;

   transformM(square, updater1);
   dbVec = updater1;

   // ApplyTransform to scalar vector with Val 3.0
 
   {
       const VecXX v = 3.0;
       auto res1 = transform(square, v);
       dbVec = res1;

       res1 = transform1(square, v);
       dbVec = res1;


       VecXX updater = 3.0;

       transformM(square, updater);
       dbVec = updater;

       transform1(square, updater);
       dbVec = updater;
   }

 

}


void applyBinaryLambdasToVector()
{

    auto add1 = [](auto x, auto y) { return x + y; };


    //default transform unroll 4
    //basic binary op 
    std::vector<double> dbVec;
    const auto  x1 = getVector();
    const auto  y1 = getVector()+2.0;
    auto res1 = transform(add1, x1,y1);
    dbVec = res1;

    auto res2 = transform(add1, x1, 10.0);
    dbVec = res2;

    auto res3 = transform(add1, 100.0,y1);
    dbVec = res3;

   // double unused = 0.0;

    //apply transforms no unroll
    { //only single unroll
 
        std::vector<double> dbVec1;

        auto add = [](auto x, auto y) { return x + y; };
        //basic binary op 
        const auto  x = getVector();
        const auto  y = getVector();// +2.0;
        auto resa = transform1(add, x, y);
        dbVec1 = resa;

        auto resb = transform1(add, x, VecXX(10.0)); //creating  vector containing a scalar
        dbVec1 = resb;

        auto resc = transform1(add, VecXX(100.0), y);
        dbVec1 = resc;


        auto resd = transform1(add, x, 10.0);
        dbVec1 = resd;

  
        auto rese = transform1(add, 100.0, y);
        dbVec1 = rese;

    }

    // binary transform in place  modify  unroll x4
    // transforM  M is for modify 
    {
        auto add_mod = [](auto& x, auto y) {  x += y;  return x; };
        //basic binary op 
        std::vector<double> dbVecA;
        const auto  x = getVector();
        auto xx = x;
        const auto  y = getVector() + 2.0;
        auto yy = y;
        transformM(add_mod, xx, y);
        dbVecA = xx;

        auto add_modxx = [](auto& x, auto y) { x += y; return x; };
        transformM(add_modxx, xx, 10.0);
        dbVecA = xx;

        auto add_mody = [](auto x, auto& y) { y += x; return y; };
        transformM(add_mody, 100.0, yy);
        dbVecA = yy;
    }


}



void applySparseLambdasToVector()
{

    auto poly1 = [](auto x) { return x*x -20.2*x +3.1; };


    //default transform unroll 4
    //basic binary op 
    std::vector<double> dbVec;
    const auto  x = getVector();

    auto resultVec = transform(poly1, x);
    dbVec = resultVec;

    auto conditionLambda = [](auto x) { return (x > VecXX::scalar(16.9)) && (x < VecXX::scalar(20.1)); };//boolean lambda 

   //auto lambda = [](auto x) {return iff(x > VecXX::scalar(30.0), -x * x, x); };
    auto lambda = [](auto x) {x = 0; return 666; };

    sparseTransform(x, resultVec, lambda, conditionLambda); //clamping

    dbVec = resultVec;

}





void accumulateVector()
{
    auto  vec = getVector(10000000);

   //NB accumulate really means reduce  so lets try max value
   auto lambdaMax = [](auto lhs, auto rhs) { return max(lhs, rhs); };
   auto maxVal = reduce(vec, lambdaMax);

    vec += 1.0 / 6.0; 
    auto sum = [](auto x, auto y) {return x + y; };
    auto acc = reduce(vec, sum);


    auto NULL_Vec = (VecXX::INS(0.0));
    auto KhanAdd = [c = NULL_Vec](auto sum, auto rhs) mutable
    {
        auto y = rhs - c;
        auto t = sum + y;
        c = (t - sum);
        c = c - y;
        sum = t;
        return t;
    };
    double KnAdd = reduce(vec, KhanAdd);


   ignore(acc);
   ignore(KnAdd);
   ignore(maxVal);
   
}


void transformReduce()
{  

    auto  test = getVector(1000);

    // we use the cast VecXX::INS  to cast the hard scalar value to a register full of the values of appropriate float/double, 
    // helps to keep the lambda types consistent and manageable
    //also helps if we switch between name spaces and instruction sets
    auto oneIfOddLmbda = [](auto x) { return iff(((x - VecXX::INS((2.0) * floor(x * VecXX::INS(0.5)))) > VecXX::INS(0.)), VecXX::INS(1.0), VecXX::INS(0.0)); };
    auto oneIfEvenLmbda = [](auto x) { return iff(((x - VecXX::INS((2.0) * floor(x * VecXX::INS(0.5)))) <= VecXX::INS(0.)), VecXX::INS(1.0), VecXX::INS(0.0)); };
    auto Add = [](auto x, auto y) { return x + y; };
    auto theCount = transformReduce(oneIfOddLmbda, Add, test);

    auto SQR = [](auto X) {return X * X; };
    auto sumSquares = transformReduce(SQR, Add, test);

    //return val of X or zero  not logical XOR
    auto XorZero = [=](auto x) {  return x * oneIfOddLmbda(x); }; //capturing a lambda 
    auto XorZeroEven = [=](auto x) {  return x * oneIfEvenLmbda(x); }; //capturing a lambda 

    using namespace JOIN;

    auto oddSqr = XorZero | SQR;
    auto evenSqr = XorZeroEven | SQR;

    auto sumOddSquares = transformReduce(oddSqr, Add, test);
    auto sumEvenSquares = transformReduce(evenSqr, Add, test);

    auto sumboth = sumEvenSquares + sumOddSquares;

    ignore(sumSquares);
    ignore(sumboth);
    ignore(theCount);

}



void branching()
{
    //using namespace VecD4D;

    auto test = getVector(100);
    auto vA = test-50.0;
    auto vB = -vA;

    auto ABiggerThanMinusA = vA > vB;  //boolean

    std::vector<double> stlVec;

   //    auto branch = ApplySelection(ABiggerThanMinusA, vA, vB);
    auto branch = select(ABiggerThanMinusA, vA, vB);
 //   auto branchProb =  ApplySelectionOperation(ABiggerThanMinusA, vA, vB); //problem
    auto branchProb = select(ABiggerThanMinusA, vA, vB); //problem
    stlVec = branchProb; // inspect vector in debugger
   auto branchScalar = select(ABiggerThanMinusA, VecXX::scalar(1.0), VecXX::scalar(-1.0));
   // 
   // auto branch = ApplySelectionOperationC(ABiggerThanMinusA, VecXX::scalar(1.0), VecXX::scalar(-1.0) );

    auto ABiggerThan= [](auto A ) { return A > VecXX::reg(20.1); };

   // auto branch = select(ABiggerThan,test, vA, vB);
   // stlVec = branch;

    auto branchC = select(ABiggerThan, test, VecXX::scalar(1.1), VecXX::scalar(1000.1));
    stlVec = branchC;

    auto truFunc = [](auto x) {return x * x + 100.1; };
    auto falseFunc = [](auto x) {return x+21.1; };

    auto branchF = selectTransform(ABiggerThan, test, truFunc, falseFunc);
    stlVec = branchF;

    auto anotherCalc = filterTransform( ABiggerThan, test, truFunc, falseFunc);
    stlVec = anotherCalc;

    // issue with branching from boolean vector
    //todo applyselction if ( VecBOOL, B);
}



void filtersAndViews()
{
    std::vector<double> stlVec;
    auto  test = getVector(1000);

    stlVec = test;
    //must be boolean lambda
    auto evenLmbda = [](auto x) { return ((x - VecXX::INS((2.0) * floor(x * VecXX::INS(0.5)))) == VecXX::INS(0.0)); };

    auto view = filter(evenLmbda, test);
    stlVec = view;
    auto SQR = [](auto x) {return x * x; };

   // ApplyUnitaryOperation(SQR, view);
    transformM(SQR, view);
    stlVec = view;

    view.writeView(test); //over write back to vector
    stlVec = test;

    ///// using a counted if ///
    auto  test1 = getVector(100);
    auto view1 =  countedFilter(evenLmbda, test1,20);
    stlVec = view1;

    transformM(SQR, view1);
    stlVec = view1;

    view1.writeView(test1); //over write back to vector
    stlVec = test1;

    /////////////////////////////

    auto  test2 = getVector(100);
    auto view2 = countedFilter(evenLmbda, test2, 30);

    transformWrite(SQR, view2, test2);
    stlVec = test2;

    ///////////////////////////////////

    //binary filter

    auto  v = getVector(180);
    //returns a tupple
    auto res = binaryFilter(evenLmbda, v);
    auto& evens = std::get<0>(res);
    auto& odds = std::get<1>(res);

    stlVec = evens;
    stlVec = odds;

    auto overHundred = [](auto x) {  return (x >= VecXX::INS(100.0)); };
    auto lessThanOnetwenty = [](auto x) {  return x <= VecXX::INS(120.0); };

    using namespace JOIN;

    auto betweenHundredAndOnetwenty = overHundred && lessThanOnetwenty;

    auto evenFilterView = filter(betweenHundredAndOnetwenty, evens);
    stlVec = evenFilterView;

    auto copyOfEvenView = evenFilterView;

   transformM(SQR, evenFilterView);
   evenFilterView.writeView(v);
   stlVec = v;


    auto newView =transform(SQR, copyOfEvenView);
    newView.writeView(v);
    stlVec = v;


}




void complexExampleOfJoiningfiltersAndViews()
{
 
    std::vector<double> stlVec;
    auto  vec = getVector(1000);
    //must be boolean lambda
    auto SQR = [](auto x) { return (x * x); };

    auto AddTen = [](auto x) { return (x + VecXX::scalar(10.0)); };

    auto negateIfOverTen = [](auto x) { return iff(x > VecXX::scalar(10), -x, x); };

    auto isEven = [](auto X) { return (floor(X) / 2.0 - floor(floor(X) / 2.0)) < 0.0001; };

    using namespace PIPE;
    // PIPE uses operator | to feed vector to filters and views to views via filter
    // it uses operator > to feed transform vector/view into next vector/view
    // 
    //filter vector for even numbers and then SQR all elements , then negate   then add ten
    auto complexResultView = (vec | isEven) > SQR > negateIfOverTen > AddTen;
    //dump to stl vec for debug
    stlVec = complexResultView;


    // JOIN fuses the operations into compound operations
    // so that they operate on the same data sequentially while in registers
    using namespace JOIN;
    //compose the lambdas into one that can be applied sequentially at register level
    //expression template in JOIN
    auto complexSingleLambda = SQR | negateIfOverTen | AddTen;


    auto complexResultView1 = (vec | isEven) > complexSingleLambda;
    //dump to stl vec for debug
    stlVec = complexResultView1;

    // boolean tweaks to filters

    auto isOdd = !isEven;
    auto complexResultView2 = (vec | isOdd) > complexSingleLambda;
    //dump to stl vec for debug
    stlVec = complexResultView2;


    auto isLessThanHundred = [](auto x) { return x < VecXX::scalar(100.); };

    auto isLessThan20 = [](auto x) { return x < VecXX::scalar(20.); };

    auto twenty_to_hundred = (!isLessThan20) && isLessThanHundred && isOdd;


    auto complexResultView3 = (vec | twenty_to_hundred) > complexSingleLambda;
    //dump to stl vec for debug
    stlVec = complexResultView3;
    complexResultView3.writeView(vec);
    stlVec = vec;

  

}






template<typename x>
auto blackScholes(const x& S, const x& K, const x& t, const x& r, const x& sigma)
{
    auto invK = 1.0 / K;
    auto discountedRate = exp(-r * t);
    auto S_invK = S * invK;
    auto log_sK = log(S_invK);
    auto rootT = sqrt(t);
    auto sigmaRootT = rootT * sigma;
    auto invSigmaRootT = 1.0 / sigmaRootT;
    auto halfSigmaSqrd_t = (0.5 * sigma * sigma + r) * t;
    auto d1 = invSigmaRootT * (log_sK + halfSigmaSqrd_t);
    auto d2 = d1 - sigmaRootT;
    auto normD1 = cdfnorm(d1);
    auto normD2 = cdfnorm(d2);
    auto C = S * normD1 - K * discountedRate * normD2;
    auto delta = normD1;
    return C;
}


template<typename INS>
void call_BSDerivAll()
{

    int size = 20;// 1024;
    VecXX spot(size);
    double strike = 40.0;
    double time = 1.0;
    double rate = 0.1;
    VecXX vol(size);

    std::vector<double> ones(20, 1.0);
    VecXX Vec1(ones);

    auto spotV = Vec1 * 42.0;
    std::vector<double> two(20, 42.0);
    VecXX Vec2(two);
    VecXX strikeV(Vec2);


    std::vector<double> timeV(20, time);
    VecXX timeVec(timeV);

    std::vector<double> rateV(20, rate);
    VecXX rateVec(rateV);

    std::vector<double> volV(20, 0.2);
    VecXX volVec(volV);

    for (int i = 0; i < size; i++)
    {
        vol[i] = 0.20;//+0.001* static_cast<double>(i-500);

        spot[i] = static_cast<double>(42.0);
    }


    std::vector<double> dbgVec;

    auto basePrice = blackScholes(spotV, strikeV, timeVec, rateVec, volVec);
    double blip = 0.00000001;
    auto bumpSpot = spotV + blip;
    auto bumpPriceVector = blackScholes(bumpSpot, strikeV, timeVec, rateVec, volVec);
    auto delta = (bumpPriceVector - basePrice) / blip;
    dbgVec = delta;
    auto dbgVecOther = dbgVec;

    // force function  to take double or differential vectors, standard conversion is to make 
    // null derivatives, using D() function creates differential vectors that are differentiable
    // so use  a D() around the parameter that you want to be differentiated against

    // explicit creation of template with derivative vector type causes all the other arguments to need to be converted
    // so we create non differentiable versions of derivative vector
    // wrapping a D() around the vecX input creates a temporary differentiable vector
    // which causes the derivative to be calculated using that vector
    auto pricesDeltaMagic = blackScholes< VecxD>(D(spotV), strikeV, timeVec, rateVec, volVec);
    auto delta1 = pricesDeltaMagic.derivative();
    dbgVecOther = delta1;

    auto pricesVegaMagic = blackScholes< VecxD>(spotV, strikeV, timeVec, rateVec, D(volVec));
    auto vega = pricesVegaMagic.derivative();

    dbgVecOther = vega;

    auto volVecBump = volVec + blip;
    bumpPriceVector = blackScholes(spotV, strikeV, timeVec, rateVec, volVecBump);
    auto sensi = (bumpPriceVector - basePrice) / blip;
    dbgVec = sensi;


    auto pricesthetaMagic = blackScholes< VecxD>(spotV, strikeV, D(timeVec), rateVec, volVec);
    auto theta = pricesthetaMagic.derivative();
    dbgVecOther = theta;

    auto thetaVecBump = timeVec + blip;
    bumpPriceVector = blackScholes(spotV, strikeV, thetaVecBump, rateVec, volVec);
    sensi = (bumpPriceVector - basePrice) / blip;
    dbgVec = sensi;


    ignore(strike);

}





void experimentalAAD()
{
    call_BSDerivAll<VecXX::INS>();
}





double blackScholes_1(double S, double K, double t, double r, double sigma)
{
    double invK = 1.0 / K;
    double discountedRate = exp(-r * t);
    double S_invK = S * invK;
    double log_sK = log(S_invK);
    double rootT = sqrt(t);
    double sigmaRootT = rootT * sigma;
    double invSigmaRootT = 1.0 / sigmaRootT;
    double halfSigmaSqrd_t = (0.5 * sigma * sigma + r) * t;
    double d1 = invSigmaRootT * (log_sK + halfSigmaSqrd_t);
    double d2 = d1 - sigmaRootT;
    double normD1 = cdfnorm(d1);
    double normD2 = cdfnorm(d2);
    double C = S * normD1 - K * discountedRate * normD2;
    //double delta = normD1;
    return C;
}



double  blackScholes_2(double S, double K, double t, double r, double sigma)
{
    auto invK = 1.0 / K;
    auto discountedRate = exp(-r * t);
    auto S_invK = S * invK;
    auto log_sK = log(S_invK);
    auto rootT = sqrt(t);
    auto sigmaRootT = rootT * sigma;
    auto invSigmaRootT = 1.0 / sigmaRootT;
    auto halfSigmaSqrd_t = (0.5 * sigma * sigma + r) * t;
    auto d1 = invSigmaRootT * (log_sK + halfSigmaSqrd_t);
    auto d2 = d1 - sigmaRootT;
    auto normD1 = cdfnorm(d1);
    auto normD2 = cdfnorm(d2);
    auto C = S * normD1 - K * discountedRate * normD2;
  //  auto delta = normD1;
    return C;
}




// make automatic variables auto 
VecXX  blackScholes_3(const VecXX& S, const VecXX& K, const VecXX& t, const VecXX& r, const VecXX& sigma)
{
    auto invK = 1.0 / K;
    auto discountedRate = exp(-r * t);
    auto S_invK = S * invK;
    auto log_sK = log(S_invK);
    auto rootT = sqrt(t);
    auto sigmaRootT = rootT * sigma;
    auto invSigmaRootT = 1.0 / sigmaRootT;
    auto halfSigmaSqrd_t = (0.5 * sigma * sigma + r) * t;
    auto d1 = invSigmaRootT * (log_sK + halfSigmaSqrd_t);
    auto d2 = d1 - sigmaRootT;
    auto normD1 = cdfnorm(d1);
    auto normD2 = cdfnorm(d2);
    auto C = S * normD1 - K * discountedRate * normD2;
   // auto delta = normD1;
    return C;
}


void optionPricerDriver()
{
    double S = 40.0;  // current spot price
    double K = 40;    // strike
    double t = 1.0;   // 1 year maturity
    double r = 0.1;   // 10% , risk free rate of return
    double sigma = 0.1;   //10%  sigma  volaiity of spot price 

    double optionPrice1 = blackScholes_1(S, K, t, r, sigma);

   

    std::vector stl(100, 0.0);
    std::iota(begin(stl), end(stl), 0.0);

    VecXX sigXX = stl;
    sigXX *= 0.001;
    sigXX += sigma;

    std::vector<double> testVec = sigXX;

    auto optionPrice3 = blackScholes_3(S, K, t, r, sigXX);

    std::vector<double> prices = optionPrice3;

    ignore(optionPrice1);


}


