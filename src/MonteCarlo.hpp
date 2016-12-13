#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include "RandomGenerator.h"
#include "Divident.hpp"
#include <map>
#include <functional>

double MonteCarlo_EUdiv(double S0, double r, double vol, double T,
		std::function<double(double)> payoff,
		std::map<double, Divident*> divs, size_t size);

double MonteCarlo_EUdiv_CV(double S0, double r, double vol, double T,
		std::function<double(double)> payoff,
		std::map<double, Divident*> divs,
		size_t npath, double BS_Price);

#endif
