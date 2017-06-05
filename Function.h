#pragma once
#include "InterfaceFunction.h"
#include "MathVector.h"

class Function :
	public InterfaceFunction {
public:
	typedef double(*math_func)(double);			//макрос математической функции одной переменной
	typedef double(*math_n_func)(double, unsigned int);	//макрос последовательности математических функций одной переменной
	const double pi = 3.141592653589793, e = 2.718281828459045, CONST_EPS = 0.0001, GOLDEN_SECTION = (1 + sqrt(5)) / 2, H_INTEGR_CONST = 0.01;
private:
	MathVector hyp_polinom_vector;
	math_func this_appropriate_func;
protected:
	Function(MathVector);
public:
	Function(math_func);
	Function(const Function &);
	Function & operator=(const Function &);
	double operator()(double);
	virtual double n_order_deriv(double, unsigned int);
	double argmin_golden_section(double, double);
	double n_deriv_argmin_golden_section(unsigned int, double, double);
	double abs_argmin_golden_section(double, double);
	double abs_n_deriv_argmin_golden_section(unsigned int, double, double);
	double argmax_golden_section(double, double);
	double n_deriv_argmax_golden_section(unsigned int, double, double);
	double abs_argmax_golden_section(double, double);
	double abs_n_deriv_argmax_golden_section(unsigned int, double, double);
	short int signum(double, double, double);
	short int n_deriv_signum(unsigned int, double, double, double);
	void equation_sol_bisection(double, double, double, double &, int &, double &, int = 10000, double = 0.0001);
	void equation_sol_iteration(double, double, double, double &, int &, double &, int = 10000, double = 0.0001);
    void equation_sol_chord(double, double, double, double &, int &, double &, int = 10000, double = 0.0001);
	void equation_sol_tangent(double, double, double, double &, int &, double &, int = 10000, double = 0.0001);
	void integral_Simpson(double, double, double &, double &);
	~Function();
};