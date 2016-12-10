#ifndef DIVIDENT_HPP
#define DIVIDENT_HPP

enum class DividentType { PROPORTIONAL, FIXED };
class Divident {
public:
	virtual double getDivident(double St, double t) const {};
	double getTime() const;
	DividentType getType() const;
protected:
	double time, div;
	DividentType type;
};

// Proportional divident
class Divident_Proportional : public Divident {
public:
	Divident_Proportional(double _time, double _div);
	double getDivident(double St, double t) const override;
};

// Fixed divident
class Divident_Fixed : public Divident {
public:
	Divident_Fixed(double _time, double _div, double _rf);
	double getDivident(double St, double t) const override;
private:
	double rf;
};


#endif
