#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include "RandomGenerator.h"

#define ALL_ELE(V) V.begin(), V.end()

const auto PUT_PRICE = 6.61665465864117;

typedef std::function<RDN_GEN()> Factory;

void priceEU_PUT() {
	const double S = 50.0;
	const double K = 55.0;
	const double vol = 0.3;
	const double r = 0.04;
	const double T = 0.5;
	auto payoff = [=](double spot) { return exp(-r*T) * std::max(0.0, K - spot); };
	auto mc_price = [=](double z) {return S * exp((r - 0.5*pow(vol, 2)) * T + vol*sqrt(T)*z); };

	auto runMC = [payoff, mc_price](auto normal) {
		auto size = normal.size();
		std::vector<double> put(size);
		std::transform(normal.begin(), normal.begin()+size, put.begin(), [=](double z) {return payoff(mc_price(z)); });
		std::cout << size;
		return std::accumulate(ALL_ELE(put), 0.0) / size;
	};

	//std::cout << "The Inverse Transform method: " << runMC(getInvTransNormalGen(), 10000) << std::endl;
	//std::cout << "The AR method: " << runMC(getARNormalGen(), 10000) << std::endl;
	//std::cout << "The Box-Muller method: " << runMC(getBoxMullerNormalGen(), 10000) << std::endl;

	auto tryMehtod = [runMC](std::vector<double> v, auto listgen) {
		int k = 1;
		while (k <= 512) {
			std::vector<double> vv(v.begin(), v.begin() + k * 10000);
			double tmp_price = runMC(listgen(vv));
			std::cout << "\t" << tmp_price << "\t" << abs(PUT_PRICE - tmp_price) << std::endl;
			k *= 2;
		}
	};

	std::vector<double> v(5120000);
	std::generate(v.begin(), v.end(), getLinCongUniGen());
	std::cout << "----------------- The Inverse Transform method ----------------" << std::endl;
	tryMehtod(v, getInverselist);
	std::cout << "----------------- The AR method ----------------" << std::endl;
	tryMehtod(v, getARlist);
	std::cout << "----------------- The Box-Muller method ----------------" << std::endl;
	tryMehtod(v, getBMlist);
}

void HW6() {
	// Price EU option
	std::cout.precision(12);
	priceEU_PUT();
}

void HW7() {
	const double S = 50.0;
	const double K = 55.0;
	const double vol = 0.3;
	const double r = 0.04;
	const double T = 0.5;
	auto payoff = [=](double spot) { return exp(-r*T) * std::max(0.0, K - spot); };
	auto mc_price = [=](double z) {return S * exp((r - 0.5*pow(vol, 2)) * T + vol*sqrt(T)*z); };
	auto runMC = [payoff, mc_price](auto normal) {
		auto size = normal.size();
		std::vector<double> put(size);
		std::transform(normal.begin(), normal.begin()+size, put.begin(), [=](double z) {return payoff(mc_price(z)); });
		std::cout << size;
		return std::accumulate(ALL_ELE(put), 0.0) / size;
	};

	// Control variate
	//auto runMCCV = [payoff, mc_price, runMC, r, T, S](auto normal) {
	//	int size = normal.size();
	//	std::vector<double> stock_price(size);
	//	std::vector<double> option_price(size);
	//	std::transform(normal.begin(), normal.end(), stock_price.begin(), mc_price);
	//	std::transform(stock_price.begin(), stock_price.end(), option_price.begin(), payoff);
	//	auto s_mean = std::accumulate(ALL_ELE(stock_price), 0.0) / size;
	//	auto v_mean = std::accumulate(ALL_ELE(option_price), 0.0) / size;
	//	double up = 0.0, down = 0.0;
	//	for (int i = 0; i < size; ++i) {
	//		up += (stock_price[i] - s_mean) * (option_price[i] - v_mean);
	//		down += pow(stock_price[i] - s_mean, 2.0);
	//	}
	//	auto b = up / down;
	//	std::vector<double> c_cv(size);
	//	auto FVs0 = exp(r*T) * S;
	//	for (int i = 0; i < size; ++i) {
	//		c_cv[i] = option_price[i] - b * (stock_price[i] - FVs0);
	//	}
	//	return std::accumulate(ALL_ELE(c_cv), 0.0) / size;
	//};

	//std::vector<double> Zs(5120000);
	//std::generate(ALL_ELE(Zs), getBoxMullerNormalGen());
	//int k = 1;
	//while(k<=512){
	//	std::vector<double> vv(Zs.begin(), Zs.begin() + k * 10000);
	//	double tmp_price = runMCCV(vv);
	//	std::cout << k*10000 << "\t" << tmp_price << "\t" << abs(PUT_PRICE - tmp_price) << std::endl;
	//	k *= 2;
	//}

	//Antithetic Variables
	auto runMCAV = [payoff, mc_price](auto normal) {
		auto size = normal.size();
		std::vector<double> put(size);
		std::vector<double> put2(size);
		std::transform(normal.begin(), normal.end(), put.begin(), [=](double z) {return payoff(mc_price(z)); });
		std::transform(normal.begin(), normal.end(), put2.begin(), [=](double z) {return payoff(mc_price(-z)); });
		return (std::accumulate(ALL_ELE(put), 0.0)
			+ std::accumulate(ALL_ELE(put2), 0.0)) / (2 * size);
	};
	std::vector<double> Zs(5120000);
	std::generate(ALL_ELE(Zs), getInvTransNormalGen());
	int k = 1;
	while(k<=512){
		std::vector<double> vv(Zs.begin(), Zs.begin() + k * 10000);
		double tmp_price = runMCAV(vv);
		std::cout << k*10000 << "\t" << tmp_price << "\t" << abs(PUT_PRICE - tmp_price) << std::endl;
		k *= 2;
	}


	//// Moment Matching
	//auto runMCMM = [payoff, mc_price, r, T, S](auto normal) {
	//	int size = normal.size();
	//	std::vector<double> stock_price(size);
	//	std::transform(normal.begin(), normal.end(), stock_price.begin(), mc_price);
	//	auto s_mean = std::accumulate(ALL_ELE(stock_price), 0.0) / size;
	//	auto FVs0 = exp(r*T) * S;
	//	auto factor = FVs0 / s_mean;
	//	std::transform(ALL_ELE(stock_price), stock_price.begin(), [=](auto k) {return k*factor; });
	//	std::vector<double> v_mm(size);
	//	std::transform(ALL_ELE(stock_price), v_mm.begin(), payoff);
	//	return std::accumulate(ALL_ELE(v_mm), 0.0) / size;
	//};

	//std::vector<double> Zs(5120000);
	//std::generate(ALL_ELE(Zs), getBoxMullerNormalGen());
	//int k = 1;
	//while(k<=512){
	//	std::vector<double> vv(Zs.begin(), Zs.begin() + k * 10000);
	//	double tmp_price = runMCMM(vv);
	//	std::cout << k*10000 << "\t" << tmp_price << "\t" << abs(PUT_PRICE - tmp_price) << std::endl;
	//	k *= 2;
	//}


/*	auto runMCCVMM = [payoff, mc_price, r, T, S](auto normal) {
		int size = normal.size();
		std::vector<double> stock_price(size);
		std::transform(normal.begin(), normal.end(), stock_price.begin(), mc_price);
		auto s_mean = std::accumulate(ALL_ELE(stock_price), 0.0) / size;
		auto FVs0 = exp(r*T) * S;
		auto factor = FVs0 / s_mean;
		std::transform(ALL_ELE(stock_price), stock_price.begin(), [=](auto k) {return k*factor; });
		std::vector<double> v_mm(size);
		std::transform(ALL_ELE(stock_price), v_mm.begin(), payoff);
		auto v_mean = std::accumulate(ALL_ELE(v_mm), 0.0) / size;
		double up = 0.0, down = 0.0;
		for (int i = 0; i < size; ++i) {
			up += (stock_price[i] - FVs0) * (v_mm[i] - v_mean);
			down += pow(stock_price[i] - FVs0, 2.0);
		}
		auto b = up / down;
		std::vector<double> c_cv(size);
		for (int i = 0; i < size; ++i) {
			c_cv[i] = v_mm[i] - b * (stock_price[i] - FVs0);
		}
		std::cout << b << " " << std::accumulate(ALL_ELE(c_cv), 0.0)
			- std::accumulate(ALL_ELE(v_mm), 0.0) << std::endl;
		return std::accumulate(ALL_ELE(c_cv), 0.0) / size;
	};

	std::vector<double> Zs(5120000);
	std::generate(ALL_ELE(Zs), getBoxMullerNormalGen());
	int k = 1;
	while(k<=512){
		std::vector<double> vv(Zs.begin(), Zs.begin() + k * 10000);
		double tmp_price = runMCCVMM(vv);
		std::cout << k*10000 << "\t" << tmp_price << "\t" << abs(PUT_PRICE - tmp_price) << std::endl;
		k *= 2;
	}*/
}

int main() {

	//std::vector<double> v(100);
	//std::generate(v.begin(), v.end(), getLinCongUniGen());
	//std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cout, "\n"));
	//std::cout << std::accumulate(v.begin(), v.end(), 0.0) / v.size() << std::endl;
	//std::cout << std::accumulate(v.begin(), v.end(), 0.0, [](auto a, auto b) {return a+(b)*(b);}) / v.size() << std::endl;

	std::cout.precision(12);

	HW7();



	return 0;
}
