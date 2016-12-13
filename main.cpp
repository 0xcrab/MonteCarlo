#include <iostream>
#include "src/MonteCarlo.hpp"
#include <cmath>
#include <exception>
#include "src/NumericTools/Midpoint Rule.h"
using namespace std;

double BSPrice(double St, double _K, double _r, double _q, double _sigma, double dt, bool isPut)
{
	double d1 = (log(St / _K) + (_r - _q + _sigma*_sigma / 2.) * dt) / (_sigma*sqrt(dt));
	double d2 = (log(St / _K) + (_r - _q - _sigma*_sigma / 2.) * dt) / (_sigma*sqrt(dt));
	double price;
	// for put option
	if (!isPut)
		price = cum_dist_normal(d1)*St*exp(-_q*dt)
		- cum_dist_normal(d2)*_K*exp(-_r*dt);
	// for call option
	else
		price = cum_dist_normal(-d2)*_K*exp(-_r*dt)
		- cum_dist_normal(-d1)*St*exp(-_q*dt);
	return price;
}

int main(){

	double S0 = 50;
	double K = 55.55;
	double r = 0.02;
	double vol = 0.3;
	double T = 7.0/12.0;
	auto payoff = [=](double St){return max(K-St, 0.0);};
	Divident_Fixed		  div1(2.0/12.0, .5, r);
	Divident_Proportional div2(4.0/12.0, 0.01);
	Divident_Fixed		  div3(6.0/12.0, .75, r);
	map<double, Divident*> divs{{div1.getTime(), &div1},
								{div2.getTime(), &div2},
								{div3.getTime(), &div3},
	};
	double bs = BSPrice(S0, K, r, 0, vol, T, true);
	cout << "bs = " << bs << endl;
	cout << MonteCarlo_EUdiv_CV(S0, r, vol, T, payoff, divs, 5000, bs) << endl;

	auto tryMethod = [=](int npath, auto method){
		cout << npath << '\t' << method(S0, r, vol, T, payoff, divs, npath) << endl;
	};

	for(int k=0; k<8; k++)
		tryMethod(10000 * (1 << k), MonteCarlo_EUdiv);
	cout << "--------------------------------------------------" << endl;
	auto tryMethod2 = [=](int npath, auto method){
		cout << npath << '\t' << method(S0, r, vol, T, payoff, divs, npath, bs) << endl;
	};
	for(int k=0; k<7; k++)
		tryMethod2(10000 * (1 << k), MonteCarlo_EUdiv_CV);

	return 0;
}
