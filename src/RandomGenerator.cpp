#include "RandomGenerator.h"
#include <algorithm>
#include <math.h>
#include <vector>
//#include <iostream>
//using std::cout;
//using std::endl;

class LinCongUniGen{
public:
	LinCongUniGen(int _x0, int _a, int _c, unsigned int _k) 
		: x_last(_x0), a(_a), c(_c), k(_k){}
	int a, c;
	uBigInt x_last;
	unsigned int k;
	double operator()() {
		uBigInt x_new = (a * x_last + c) % k;
		x_last = x_new;
		return double(x_new) / double(k);
	}
};

//RDN_GEN getLinCongUniGen(int x0, int a, int c, unsigned int k)
//{
//	return [=]() {
//		static uBigInt x_last = x0;
//		uBigInt x_new = (a * x_last + c) % k;
//		x_last = x_new;
//		return double(x_new) / double(k);
//	};
//}

RDN_GEN getLinCongUniGen(int x0, int a, int c, unsigned int k) {
	return LinCongUniGen(x0, a, c, k);
}

double inv_cdf(const double& quantile) {
	// This is the Beasley-Springer-Moro algorithm which can 
	// be found in Glasserman [2004]. We won't go into the
	// details here, so have a look at the reference for more info
	static double a[4] = { 2.50662823884,
		-18.61500062529,
		41.39119773534,
		-25.44106049637 };

	static double b[4] = { -8.47351093090,
		23.08336743743,
		-21.06224101826,
		3.13082909833 };

	static double c[9] = { 0.3374754822726147,
		0.9761690190917186,
		0.1607979714918209,
		0.0276438810333863,
		0.0038405729373609,
		0.0003951896511919,
		0.0000321767881768,
		0.0000002888167364,
		0.0000003960315187 };

	if (quantile >= 0.5 && quantile <= 0.92) {
		double num = 0.0;
		double denom = 1.0;

		for (int i = 0; i < 4; i++) {
			num += a[i] * pow((quantile - 0.5), 2 * i + 1);
			denom += b[i] * pow((quantile - 0.5), 2 * i + 2);
		}
		return num / denom;

	}
	else if (quantile > 0.92 && quantile < 1) {
		double num = 0.0;

		for (int i = 0; i < 9; i++) {
			num += c[i] * pow((log(-log(1 - quantile))), i);
		}
		return num;

	}
	else {
		return -1.0*inv_cdf(1 - quantile);
	}
}

RDN_GEN getInvTransNormalGen(RDN_GEN uniform)
{
	return [=]() {
		return inv_cdf(uniform());
	};
}

std::vector<double> getInverselist(std::vector<double> v)
{
	std::vector<double> normal(v.size());
	std::transform(v.begin(), v.end(), normal.begin(), inv_cdf);
	return normal;
}

RDN_GEN getARNormalGen(RDN_GEN uniform)
{
	return [=]() {
		double u1 = uniform();
		double u2 = uniform();
		double u3 = uniform();
		double X = -log(u1);
		while (u2 > exp(-0.5*pow(X - 1.0, 2.0))) {
			u1 = uniform();
			u2 = uniform();
			u3 = uniform();
			X = -log(u1);
		}
		if (u3 <= 0.5)
			X = -X;
		return X;
	};
}

std::vector<double> getARlist(std::vector<double> v) {
	auto iter = v.begin();
	auto stop = v.end();
	auto getuni = [&](auto& x) {
		if (iter == stop)	return false;
		else
		{
			x = *iter;
			iter++;
			return true;
		}
	};
	auto get3uni = [=](auto &x, auto &y, auto &z) {
		return getuni(x) && getuni(y) && getuni(z);
	};

	std::vector<double> normal;
	double u1, u2, u3;
	bool flag = true;
	while (get3uni(u1, u2, u3)) {
		double X = -log(u1);
		while (u2 > exp(-0.5*pow(X - 1.0, 2.0))) {
			if (!get3uni(u1, u2, u3)) {
				flag = false;
				break;
			}
			X = -log(u1);
		}
		if (u3 <= 0.5)
			X = -X;
		if(flag)
			normal.push_back(X);

	}
	return normal;
}


RDN_GEN getBoxMullerNormalGen(RDN_GEN uniform)
{
	return [=] {
		static bool skip = false;
		static double Z2;
		if (skip) {
			skip = false;
			return Z2;
		}
		double X = 2.0;
		double u1, u2;
		while (X > 1) {
			u1 = 2 * uniform() - 1;
			u2 = 2 * uniform() - 1;
			X = pow(u1, 2) + pow(u2, 2);
		}
		double Y = sqrt(-2.0 * log(X) / X);
		double Z1 = u1*Y;
		Z2 = u2*Y;
		skip = true;
		return Z1;
	};
}

std::vector<double> getBMlist(std::vector<double> v)
{
	auto iter = v.begin();
	auto stop = v.end();
	auto getuni = [&](auto& x) {
		if (iter == stop)	return false;
		else
		{
			x = (*iter) * 2 - 1;
			iter++;
			return true;
		}
	};
	auto get2uni = [=](auto& x, auto& y) {
		return getuni(x) && getuni(y);
	};
	std::vector<double> normal;
	double flag = true;
	while (flag) {
		double X = 2.0;
		double u1, u2;
		while (X > 1) {
			if (!get2uni(u1, u2)) {
				flag = false;
				break;
			}
			X = pow(u1, 2) + pow(u2, 2);
		}
		double Y = sqrt(-2.0 * log(X) / X);
		if (flag) {
			normal.push_back(u1*Y);
			normal.push_back(u2*Y);
		}
	}
	return normal;
}
