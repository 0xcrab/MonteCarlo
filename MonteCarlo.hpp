#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include "RandomGenerator.h"
#include "Option.h"
#include <map>

double MonteCarlo_EUdiv(double S0, double r, double vol,
		std::map<double, Divident> divs, size_t size);


#endif
