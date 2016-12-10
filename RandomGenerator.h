#pragma once

#include <functional>
#include <vector>

typedef unsigned __int64 uBigInt;

// Define prototype of random number generator
typedef std::function<double()> RDN_GEN;
//class RDN_GEN {
//	virtual double operator()() = 0;
//};

// Linear Congruential Uniform Generator
RDN_GEN getLinCongUniGen(int x0=1, int a=39373, int c=0, unsigned int k=(1uL<<31)-1uL);

// Inverse Transform Method to generate normal variable
RDN_GEN getInvTransNormalGen(RDN_GEN uniform = getLinCongUniGen());
std::vector<double> getInverselist(std::vector<double> v);

// Acceptance-Rejection normal variable generator
RDN_GEN getARNormalGen(RDN_GEN uniform = getLinCongUniGen());
std::vector<double> getARlist(std::vector<double> v);

// Box-Muller normal variable generator
RDN_GEN getBoxMullerNormalGen(RDN_GEN uniform = getLinCongUniGen());
std::vector<double> getBMlist(std::vector<double> v);
