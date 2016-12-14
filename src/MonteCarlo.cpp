#include "MonteCarlo.hpp"
#include <vector>
#include <algorithm>
#include "RandomGenerator.h"
#include <cmath>
#include <functional>
#include <numeric>
#include <iostream>
using namespace std;

vector<double> MonteCarlo_EUdiv_getSt(double S0, double r, double vol, double T,
		std::map<double, Divident*> divs,
		size_t npath){

	size_t size = npath * (divs.size() + 1);
	vector<double> normal(size);
	generate(normal.begin(), normal.end(), getBoxMullerNormalGen());
	auto iter_normal = normal.begin();
	auto mc_price = [=](double z, double dt) {return exp((r - 0.5*pow(vol, 2)) * dt + vol*sqrt(dt)*z); };
	double St;
	double old_t;
	vector<double> St_list;
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
		St_list.push_back(St);
		//cout << "     ==>     " << St << "\t" << payoff(St) << "\t" << exp(-r*T) << "\t" << payoff(St)*exp(-r*T) << endl;
	}
	//double sum = 0;
	//for(size_t i = 0; i<price.size(); i++){
		//sum += price[i];
		//cout << i << '\t' << sum  << "\t" << price[i] << endl;
	//}
	//cout << accumulate(price.begin(), price.end(), 0.0) << " " <<  static_cast<double>(price.size()) << endl;
	return St_list;
}

vector<double> MonteCarlo_EUdiv_helper(double S0, double r, double vol, double T,
		function<double(double)> payoff,
		std::map<double, Divident*> divs,
		size_t npath){

	vector<double> price;
	auto St_list = MonteCarlo_EUdiv_getSt(S0, r, vol, T, divs, npath);
	double disc = exp(-r*T);
	for_each(St_list.begin(), St_list.end(), 
			[&price, &payoff, disc](auto& St){
			price.push_back(payoff(St) * disc);
			});
	
	return price;
}

double MonteCarlo_EUdiv(double S0, double r, double vol, double T,
		function<double(double)> payoff,
		std::map<double, Divident*> divs,
		size_t npath){

	auto price = MonteCarlo_EUdiv_helper(S0, r, vol, T, payoff, divs, npath);
	return accumulate(price.begin(), price.end(), 0.0) / static_cast<double>(price.size());
}

vector<double> Control_Variance(const vector<double>& price_div, const vector<double>& price_nondiv, double BS_Price){

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
	return cv;
}

double MonteCarlo_EUdiv_CV(double S0, double r, double vol, double T,
		function<double(double)> payoff,
		std::map<double, Divident*> divs,
		size_t npath, double BS_Price, 
		double BS_delta, double *output_delta, double *output_delta_cv){
	
	auto price_div = MonteCarlo_EUdiv_helper(S0, r, vol, T, payoff, divs, npath);
	map<double, Divident*> zero_divs;
	Divident_Proportional div(0.0, 0.0);
	for(auto& i : divs) zero_divs.insert({i.second->getTime(), &div});
	auto price_nondiv = MonteCarlo_EUdiv_helper(S0, r, vol, T, payoff, zero_divs, npath);

	// Compute Delta
	auto St_list = MonteCarlo_EUdiv_getSt(S0, r, vol, T, divs, npath);
	auto St_nondiv_list = MonteCarlo_EUdiv_getSt(S0, r, vol, T, zero_divs, npath);
	vector<double> delta(St_list.size());
	vector<double> delta_nondiv(St_list.size());
	double disc = exp(-r*T);
	for(size_t i=0; i<St_list.size(); i++){
		delta[i] = -(payoff(St_list[i])>0) * disc * St_nondiv_list[i] / S0 * 0.99;
		delta_nondiv[i] = -(payoff(St_nondiv_list[i])>0) * disc * St_nondiv_list[i] / S0;
	}
	auto delta_cv = Control_Variance(delta, delta_nondiv, BS_delta);
	*output_delta = accumulate(delta.begin(), delta.end(), 0.0) / static_cast<double>(delta.size());
	*output_delta_cv = accumulate(delta_cv.begin(), delta_cv.end(), 0.0) / static_cast<double>(delta_cv.size());

	auto price_cv = Control_Variance(price_div, price_nondiv, BS_Price);
	return accumulate(price_cv.begin(), price_cv.end(), 0.0) / static_cast<double>(price_cv.size());
}

