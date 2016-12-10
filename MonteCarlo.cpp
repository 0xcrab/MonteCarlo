#include "MonteCarlo.hpp"
#include <vector>
#include <algorithm>
using namespace std;

double MonteCarlo_EUdiv(double S0, double r, double vol,
		std::map<double, Divident> divs,
		size_t size){
	vector<double> v(size);
	generate(v.begin(), v.end(), getLinCongUniGen());
	
}
