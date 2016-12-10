#include "Midpoint Rule.h"

#include <iostream>
#include <iomanip>

using std::endl;
using std::cout;

double I_tol(double a, double b,double tol, F f_int) {
	
	int n = 4;
	double I_old = I_Simpson(a, b, n, f_int);
	n = 8;
	double I_new = I_Simpson(a, b, n, f_int);

	while (abs(I_old - I_new) > tol) {
		I_old = I_new;
		n = 2 * n;
		I_new = I_Simpson(a, b, n, f_int);
	}

	return I_new;
}