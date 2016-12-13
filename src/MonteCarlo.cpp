#include "MonteCarlo.hpp"
#include <vector>
#include <algorithm>
#include "RandomGenerator.h"
#include <cmath>
#include <functional>
#include <numeric>
#include <iostream>
using namespace std;

vector<double> MonteCarlo_EUdiv_helper(double S0, double r, double vol, double T,
		function<double(double)> payoff,
		std::map<double, Divident*> divs,
		size_t npath){

	size_t size = npath * (divs.size() + 1);
	vector<double> normal(size);
	generate(normal.begin(), normal.end(), getBoxMullerNormalGen());
	auto iter_normal = normal.begin();
	auto mc_price = [=](double z, double dt) {return exp((r - 0.5*pow(vol, 2)) * dt + vol*sqrt(dt)*z); };
	double St;
	double old_t;
	vector<double> price;
	while(npath--){
		St = S0;
		old_t = 0;
		for(auto& i : divs){
			auto&& div = i.second;
			St *= mc_price(*iter_normal++, i.first-old_t);
			//cout << St << "\t->";
			St -= div->getDivident(St, i.first);
			//cout << St << "\t";
			old_t = i.first;
		}
		St *= mc_price(*iter_normal++, T-old_t);
		//cout << "     ==>     " << St << "\t" << payoff(St) << "\t" << exp(-r*T) << "\t" << payoff(St)*exp(-r*T) << endl;
		price.push_back(payoff(St) * exp(-r*T));
	}
	//double sum = 0;
	//for(size_t i = 0; i<price.size(); i++){
		//sum += price[i];
		//cout << i << '\t' << sum  << "\t" << price[i] << endl;
	//}
	//cout << accumulate(price.begin(), price.end(), 0.0) << " " <<  static_cast<double>(price.size()) << endl;
	return price;
}

double MonteCarlo_EUdiv(double S0, double r, double vol, double T,
		function<double(double)> payoff,
		std::map<double, Divident*> divs,
		size_t npath){

	auto price = MonteCarlo_EUdiv_helper(S0, r, vol, T, payoff, divs, npath);
	return accumulate(price.begin(), price.end(), 0.0) / static_cast<double>(price.size());
}

double MonteCarlo_EUdiv_CV(double S0, double r, double vol, double T,
		function<double(double)> payoff,
		std::map<double, Divident*> divs,
		size_t npath, double BS_Price){
	
	auto price_div = MonteCarlo_EUdiv_helper(S0, r, vol, T, payoff, divs, npath);
	map<double, Divident*> zero_divs;
	Divident_Proportional div(0.0, 0.0);
	for(auto& i : divs) zero_divs.insert({i.second->getTime(), &div});
	auto price_nondiv = MonteCarlo_EUdiv_helper(S0, r, vol, T, payoff, zero_divs, npath);

	auto price_div_mean = std::accumulate(price_div.begin(), price_div.end(), 0.0) / price_div.size();
	auto price_nondiv_mean = std::accumulate(price_nondiv.begin(), price_nondiv.end(), 0.0) / price_nondiv.size();
	//cout << "nondiv mean = " << price_nondiv_mean << endl;
	//cout << "div mean = " << price_div_mean << endl;
	double up = 0.0, down = 0.0;
	for(size_t i=0; i<price_div.size(); i++){
		up += (price_div[i] - price_div_mean) * (price_nondiv[i] - price_nondiv_mean);
		down += pow(price_nondiv[i] - price_nondiv_mean, 2.0);
	}
	auto b = up / down;
	//cout << "b=" << b << endl;
	vector<double> cv(price_div.size());
	for(size_t i=0; i<price_div.size(); i++){
		cv[i] = price_div[i] - b * (price_nondiv[i] - BS_Price);
	}
	return accumulate(cv.begin(), cv.end(), 0.0) / static_cast<double>(cv.size());
}

