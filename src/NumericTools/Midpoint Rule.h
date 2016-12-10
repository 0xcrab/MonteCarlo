#pragma once

#include <functional>
#include <vector>
#include <tuple>

using std::vector;

typedef std::function<double(double)> F;

enum OptionType { EU_CALL, EU_PUT };

class EU_Option {
public:
	EU_Option(OptionType type, double S, double K, double r, double q, double T, double sigma);
	double getDelta();
	double getGamma();
	double getPrice();

private:
	/*
	C : price of the call
	S : spot price of the underlying asset
	K : strike price of the option
	T : maturity of the option
	*/
	OptionType type;
	double S, K, r, q, T, sigma;
	double delta, gamma;
	double price;
};

/* This function returns the numerical definite integral
	within [a, b] of f_int().
	a,b : defines the interval
	n : number of partition intervals
	f_int : routine evaluating f(x)
	*/
double I_midpoint(double a, double b, int n, F f_int);


/* This function returns the approximate solution for
	f(x) = 0 using Newton's Method.
	x0 : initial guess
	f(x) : given function
	ff(x) : derivative of f(x)
	tol_approx : largest admissible value of |f(x)| when
		solution is found
	tol_consec : largest admissible distance between two
		consecutive approximations when solution is found
	*/
double Newton(double x0, F f, F ff, double tol_approx, double tol_consec);


/* This function returns the approximate solution for
f(x) = 0 using Secant Method.
x0 : initial guess
f(x) : given function
tol_approx : largest admissible value of |f(x)| when
solution is found
tol_consec : largest admissible distance between two
consecutive approximations when solution is found
*/
double secant(double x_1, double x0, F f, double tol_approx, double tol_consec);

/* This function returns the cumulative distribution of Z
	t : real number
	*/
double cum_dist_normal(double t);

/* This function computes the bond price given the zero rate curve
	n : number of cash flows
	**** Be careful for the inputed cash flow time, 
	**** if input in fraction form. 2/12 will produce 0.
	**** Always use 2.0/12.0
	t_cash_flow : vector of cash flow dates ( of size n ) 
	v_cahs_flow : vector of cash flows ( of size n )
	r_zero(t) : zero rate corresponding to time t
	*/
double bondPrice(int n, std::vector<double> t_cash_flow,
	std::vector<double> v_cash_flow, F r_zero);

/* This function use bisection method to find the root of give
	function
	a,b : left and right endpoints of the initial interval
	f : given function
	tol_int : largest admissible size of active interval when solution is found
	tol_approx : largest admissible value of |f(x)| when solution is found
	*/
double bisection(double a, double b, F f, double tol_int, double tol_approx);

/* This function compute the duration and convexity of given bond
	T : bond maturity
	n : number of cash flows
	t_cash_flow : vector of cash flow dates (of size n)
	v_cash_flow : vector of cash flows (of size n)
	y : yield of the bond
	*/
std::tuple<double, double> duration_n_convexity(int n, std::vector<double> t_cash_flow,
	std::vector<double> v_cash_flow, double y);


/* This function returns a interpolation function by bootstrap using given points
	x : x coordinates
	y : y coordinates
	*/
F bootstrap(std::vector<double> x, std::vector<double> y);

/* This function calculate EU options price and greeks
	S : underlying asset price
	K : strike price
	r : annual risk-free rate
	q : continuously compounded divident rate
	sigma : implied volatility
	T : maturity
	*/



/* This funtion computes implied volatility corresponding 
	to a given option price
	C : price of the call
	S : spot price of the underlying asset
	K : strike price of the option
	T : maturity of the option*/


double I_Simpson(double left, double right, int length, F f_int);  //function using Simpson rule to calculate numeric integration
double I_tol(double a, double b, double tol, F f_int);  // function to calculate numeric integration with tolerance
//double Bond_ins(int n, vector<double> t_cash, vector<double> v_cash, F r_inst, vector<double> tols); //function to calculate bond price for instantaneous rates