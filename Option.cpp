#include "Option.h"
#include <algorithm>
#include <iostream>
using namespace std;

Option::Option(OptionType type, double T, double K, double S0)
	: _type(type), _T(T), _K(K), _S0(S0)
{
	if (type == OptionType::CALL)
		_payoff = [=](double St) {return std::max(0.0, St - K); };
	else
		_payoff = [=](double St) {return std::max(0.0, K - St); };
}

double EuOption::payoff(double t, double St) const
{
	return t == _T ? _payoff(St) : 0;
}

double AmOption::payoff(double t, double St) const
{
	return _payoff(St);
}

void DividentOption::addDivident(const Divident& d)
{
	div_list.push_back(&d);
	// If a fixed div is added, we need to alter S0
	if (d.getType() == DividentType::FIXED)
		_S0 -= d.getDivident(_S0, 0);
}

Divident_Proportional::Divident_Proportional(double _time, double _div)
{
	time = _time;
	div = _div;
	type = DividentType::PROPORTIONAL;
}

double Divident_Proportional::getDivident(double St, double t) const
{
	if (t >= time)
		return St * div;
	return 0.0;
}

Divident_Fixed::Divident_Fixed(double _time, double _div, double _rf)
{
	time = _time;
	div = _div;
	rf = _rf;
	type = DividentType::FIXED;
}

double Divident_Fixed::getDivident(double St, double t) const
{
	// S_div = S_nondiv + PV(D), if t <= time
	if (t <= time)
		return div * exp(-rf * (time - t));
	return 0.0;
}

double Divident::getTime() const
{
	return time;
}

DividentType Divident::getType() const
{
	return type;
}


void decide_change_S(double& S, double t, const Divident& d, const OptionType& op_type) {
	// Am. Put option may excercise after ex-div
	if (op_type == OptionType::PUT) {
		if (d.getType() == DividentType::PROPORTIONAL && d.getTime() <= t)
			S -= d.getDivident(S, t);
		else if (d.getType() == DividentType::FIXED && d.getTime() > t)
			S += d.getDivident(S, t);
	}
	// Am. Call option may exercise before ex-div
	else if (op_type == OptionType::CALL) {
		if (d.getType() == DividentType::PROPORTIONAL && d.getTime() < t)
			S -= d.getDivident(S, t);
		else if (d.getType() == DividentType::FIXED && d.getTime() >= t)
			S += d.getDivident(S, t);
	}
}

double EuDivOption::payoff(double t, double St) const
{
	if (t < _T) return 0.0;
	for (auto& i : div_list)
		decide_change_S(St, t, *i, _type);
	return _payoff(St);
}

double AmDivOption::payoff(double t, double St) const
{
	for (auto& i : div_list)
		decide_change_S(St, t, *i, _type);
	return _payoff(St);
}
