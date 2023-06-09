#pragma once

double europeanBinomialPricer(double S, double K, double sig, double r, double T, int N);

double europeanTrinomialPricer(double S, double K, double sig, double r, double T, int N);

double europeanTrinomialPricer1(double S, double K, double sig, double r, double T, int N);

double americanTrinomialPricer(double S, double K, double sig, double r, double T, int N);

double americanFiniteDiffPricer(double S, double K, double sig, double r, double T, int N);

double americanImplicitFiniteDiffPricer(double S, double K, double sig, double r, double T, int N);

double americanImplicitFiniteDiffPricerFast(double S, double K, double sig, double r, double T, int N);

double euroTrinomialPricerWithInit(double S, double K, double sig, double r, double T, int N);

double americanCrankNicholsonPricer(double S, double K, double sig, double r, double T, int N);

double americanTrinomialPricerUpAndOut(double S, double K, double sig, double r, double T, double H, double rebate, int N);