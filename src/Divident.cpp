#include "Divident.hpp"
#include <cmath>

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
