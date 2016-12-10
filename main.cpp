#include <iostream>
#include "src/MonteCarlo.hpp"
#include <cmath>
#include <exception>
using namespace std


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
	cout << MonteCarlo_EUdiv_CV(S0, r, vol, T, payoff, divs, 5000) << endl;

	auto tryMethod = [=](int npath, auto method){
		cout << npath << '\t' << method(S0, r, vol, T, payoff, divs, npath) << endl;
	};

	for(int k=0; k<8; k++)
		tryMethod(10000 * (1 << k), MonteCarlo_EUdiv);
	cout << "--------------------------------------------------" << endl;
	for(int k=0; k<7; k++)
		tryMethod(10000 * (1 << k), MonteCarlo_EUdiv_CV);

	return 0;
}
