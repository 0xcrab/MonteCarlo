#include "Midpoint Rule.h"

double I_Simpson(double left,double right, int length, F f_int) {
	double h = (right - left) / length;
	double result = f_int(left) / 6 + f_int(right) / 6;
	for (int i = 1; i < length; i++) {
		result = result + f_int(left + i*h) / 3;
	}
	for (int i = 1; i < (length+1); i++) {
		result = result + 2 * f_int(left + (i - 0.5)*h) / 3;
	}
	result = h*result;
	return result;
}