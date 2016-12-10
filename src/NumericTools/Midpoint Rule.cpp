#include "Midpoint Rule.h"

#include <iostream>
using std::cout;
using std::endl;

double I_midpoint(double a, double b, int n, F f_int)
{
	double h = (b - a) / n;
	double result = 0.0;
	for (int i = 1; i <= n; i++) {
		result += f_int(a + (i - 0.5)*h);
	}
	result *= h;
	return result;
}

double Newton(double x0, F f, F ff, double tol_approx, double tol_consec)
{
	double x_new = x0;
	double x_old = x0 - 1;
	cout << "x0 = " << x0 << "\t";
	cout << "f(x0) = " << f(x0) << endl;
	while (abs(f(x_new)) > tol_approx || abs(x_new - x_old) > tol_consec) {
		x_old = x_new;
		x_new = x_old - f(x_old) / ff(x_old);
		cout << "x = " << std::fixed << x_new  << '\t'
			<< "f(x) = " << std::fixed << f(x_new) << '\t' 
			<< endl;
	}
	return x_new;
}

double secant(double x_1, double x0, F f, double tol_approx, double tol_consec)
{
	double x_new = x0;
	double x_old = x_1;
	double x_oldest;
	cout << "x0 = " << x0 << "\t";
	cout << "f(x0) = " << f(x0) << endl;
	while (abs(f(x_new)) > tol_approx || abs(x_new - x_old) > tol_consec) {
		x_oldest = x_old;
		x_old = x_new;
		x_new = x_old - f(x_old) * (x_old - x_oldest) / (f(x_old) - f(x_oldest));
		cout << "x = " << std::fixed << x_new << '\t'
			<< "f(x) = " << std::fixed << f(x_new) << '\t'
			<< endl;
	}
	return x_new;
}

double cum_dist_normal(double t)
{
	//double z = abs(t);
	//double y = 1.0 / (1.0 + 0.2316419 * z);
	//double a1 = 0.319381530;
	//double a2 = -0.356563782;
	//double a3 = 1.781477937;
	//double a4 = -1.821255978;
	//double a5 = 1.330274429;
	//double pi = 3.1415926535897;
	//double m = 1 - exp(-t*t / 2)
	//		* (a1*y + a2*pow(y, 2) + a3*pow(y, 3) + a4*pow(y, 4) + a5*pow(y, 5))
	//		/ sqrt(2 * pi);
	//return t > 0 ? m : 1 - m;
	F pdf = [](double x) {return 1.0 / sqrt(atan(1.0)*8.0)*exp(-0.5*pow(x, 2.0)); };
	if (t >= 0) return 0.5 + I_tol(0, t, 1e-12, pdf);
	else return 0.5 - I_tol(t, 0, 1e-12, pdf);
}

double bondPrice(int n, std::vector<double> t_cash_flow, std::vector<double> v_cash_flow, F r_zero)
{
	double B = 0;
	double disc;
	for (int i = 0; i < n; i++) {
		disc = exp(-t_cash_flow[i] * r_zero(t_cash_flow[i]));
		B += v_cash_flow[i] * disc;
	}
	return B;
}

double bisection(double a, double b, F f, double tol_int, double tol_approx)
{
	double xl = a;
	double xr = b;
	double xm;
	auto max = [](double a, double b) {return a > b ? a : b; };
	while (max(abs(f(xl)), abs(f(xr))) > tol_approx
		|| xr - xl > tol_int) {
		xm = (xl + xr) / 2.0;
		if (f(xl)*f(xm) < 0)
			xr = xm;
		else
			xl = xm;
	}
	return xr;
}

std::tuple<double, double> duration_n_convexity(int n, std::vector<double> t_cash_flow, std::vector<double> v_cash_flow, double y)
{
	double b = 0;
	double d = 0;
	double c = 0;
	double disc;
	for (int i = 0; i < n; i++) {
		disc = exp(-t_cash_flow[i] * y);
		b += v_cash_flow[i] * disc;
		d += t_cash_flow[i] * v_cash_flow[i] * disc;
		c += pow(t_cash_flow[i], 2.0) * v_cash_flow[i] * disc;
	}
	d /= b;
	c /= b;
	return std::make_tuple(d, c);
}

F bootstrap(std::vector<double> x, std::vector<double> y)
{
	F fit = [=](double t) {
		if (t == x[0]) 
			return y[0];
		int len = y.size();
		for (int i = 1; i < len; i++) {
			if (t <= x[i]) {
				double xl = x[i - 1];
				double xr = x[i];
				double yl = y[i - 1];
				double yr = y[i];
				return ((xr - t)*yl + (t - xl)*yr) / (xr - xl);
			}
		}
		return 0.0;
	};
	return fit;
}

EU_Option::EU_Option(OptionType type, double S, double K, double r, double q, double T, double sigma)
	: type{ type }, S{ S }, K{ K },
	r{ r }, q{ q }, T{ T }, sigma{ sigma }, delta(0), gamma(0)
{
	double d1 = (log(S / K) + (r - q + 0.5*pow(sigma, 2.0))*T) / (sigma * sqrt(T));
	double d2 = d1 - sigma*sqrt(T);
	if (type == EU_CALL) {
		price = S * exp(-q*T)*cum_dist_normal(d1) - K*exp(-r*T)*cum_dist_normal(d2);
		delta = exp(-q*T)*cum_dist_normal(d1);
	}
	else {
		price = -S * exp(-q*T)*cum_dist_normal(-d1) + K*exp(-r*T)*cum_dist_normal(-d2);
		delta = -exp(-q*T)*cum_dist_normal(-d1);
	}
	gamma = 1.0 / sqrt(atan(1)*8.0) * exp(-0.5*pow(d1, 2.0)) * exp(-q*T) / (S*sigma*sqrt(T));
}

double EU_Option::getDelta()
{
	return delta;
}

double EU_Option::getGamma()
{
	return gamma;
}

double EU_Option::getPrice()
{
	return price;
}

