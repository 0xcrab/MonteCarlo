#pragma once

#include "Security.h"
#include "BinaryTreePricing.h"
#include <functional>
#include <vector>

enum class OptionType { CALL, PUT };

class Option : public Security {
public:
	Option(OptionType type, double T, double K, double S0);
	virtual double payoff(double T, double St) const = 0;
	template<typename BTree>
	OptionPrice price(double r, double q, double sigma, int depth) const;

protected:
	Option() = default;
	OptionType _type;
	double _S0;
	double _T;
	double _K;
	std::function<double(double)> _payoff;
};

// implement price/greeks method
#include "OptionPrice.tpp"

class EuOption : public virtual Option {
public:
	using Option::Option;
	double payoff(double t, double St) const override;
protected:
	EuOption() = default;
};

class AmOption : public virtual Option {
public:
	using Option::Option;
	double payoff(double t, double St) const override;
};

//struct Divident {
//	enum DIV_TYPE { PROPORTIONAL, FIXED };
//	DIV_TYPE type;
//	double time, div;
//	double getDivident(double St, double t) const{
//		double amount;
//		switch (type) {
//		case PROPORTIONAL:	amount = St * div;
//		case FIXED:			amount = div * exp(-r;
//		}
//		return amount;
//	}
//};
enum class DividentType { PROPORTIONAL, FIXED };
class Divident {
public:
	virtual double getDivident(double St, double t) const = 0;
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

class DividentOption : public virtual Option{
public:
	// If a fixed div is added, we need to alter S0
	void addDivident(const Divident& d);
	//void addDivident(const Divident& d, const Divident& args...);
	//TODO: override a payoff funtion adapting to dividents

protected:
	DividentOption() = default;
	std::vector<const Divident*> div_list;
};

class EuDivOption : public EuOption, public DividentOption {
public:
	EuDivOption(OptionType type, double T, double K, double S0)
		:Option(type, T, K, S0){}
	double payoff(double t, double St) const override;
};

class AmDivOption : public AmOption, public DividentOption {
public:
	AmDivOption(OptionType type, double T, double K, double S0)
		:Option(type, T, K, S0) {}
	double payoff(double t, double St) const override;
};