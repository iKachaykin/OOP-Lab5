#pragma once
#include "Function.h"
#include "MathVector.h"

class Polinom :
	public Function {
private:
	int polinom_pow;
	MathVector polinom_koeffs;
public:
	Polinom(int);
	Polinom(double *, int);
	Polinom(MathVector);
	template<class init_list>
	Polinom(init_list);
	Polinom(const Polinom &);
	Polinom & operator=(const Polinom &);
	double operator()(double);
	int power() const;
	std::string to_string() const;
	bool operator==(const Polinom &);
	bool operator!=(const Polinom &);
	Polinom operator+(const Polinom &);
	Polinom operator-(const Polinom &);
	Polinom operator*(const Polinom &);
	Polinom operator/(const Polinom &);
	Polinom operator%(const Polinom &);
	Polinom operator-();
	Polinom & operator+=(const Polinom &);
	Polinom & operator-=(const Polinom &);
	Polinom & operator*=(const Polinom &);
	Polinom & operator/=(const Polinom &);
	int Sturm_method(double, double);
	Polinom create_n_order_derivative(int);
	double n_order_deriv(double, unsigned int) override;
	~Polinom();
};

template<class init_list>
inline Polinom::Polinom(init_list) : Function(MathVector(init_list)){
	polinom_koeffs = MathVector(init_list);
	polinom_pow = polinom_koeffs.dimension() - 1;
}

std::ostream& operator<<(std::ostream &, Polinom &);