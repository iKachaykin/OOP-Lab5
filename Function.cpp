#include "Function.h"
#include <vector>

Function::Function(math_func this_appropriate_func) {
	this->this_appropriate_func = this_appropriate_func;
}

Function::Function(MathVector hyp_polinom_vector) {
	this->hyp_polinom_vector = hyp_polinom_vector;
	this_appropriate_func = NULL;
}

Function::Function(const Function& other) {
	this_appropriate_func = other.this_appropriate_func;
}

Function & Function::operator=(const Function& other) {
	this_appropriate_func = other.this_appropriate_func;
	return *this;
}

double Function::operator()(double x) {
	if(this_appropriate_func)
		return this_appropriate_func(x);
	double result = 0, tmp_comp = 1;
	for (int i = hyp_polinom_vector.dimension() - 1; i >= 0; i--) {
		result += hyp_polinom_vector[i] * tmp_comp;
		tmp_comp *= x;
	}
	return result;
}

double Function::n_order_deriv(double x, unsigned int order) {
	if (order == 0)
		return (*this)(x);
	double eps = order <= 4 ? CONST_EPS * pow(10, order - 1) : 0.1;
	std::vector<double> func_and_derivs_vals;
	for (int i = 0; i <= order; i++) {
		func_and_derivs_vals.push_back((*this)(x + ((int)order - 2 * i) * eps));
	}
	for (int order_index = 0; order_index < order; order_index++) {
		for (int i = 0; i < order - order_index; i++)
			func_and_derivs_vals[i] = (func_and_derivs_vals[i] - func_and_derivs_vals[i + 1]) / 2 / eps;
		func_and_derivs_vals.pop_back();
	}
	return func_and_derivs_vals[0];
}

double Function::argmin_golden_section(double left_interval_border, double right_interval_border) {
	if (left_interval_border > right_interval_border) {
		double tmp = left_interval_border;
		left_interval_border = right_interval_border;
		right_interval_border = tmp;
	}
	double x_1 = 0, x_2 = 0;
	while (fabs(right_interval_border - left_interval_border) > CONST_EPS) {
		x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		if ((*this)(x_1) >= (*this)(x_2))
			left_interval_border = x_1, x_1 = x_2, x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		else
			right_interval_border = x_2, x_2 = x_1, x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
	}
	return (left_interval_border + right_interval_border) / 2;
}

double Function::n_deriv_argmin_golden_section(unsigned int order, double left_interval_border, double right_interval_border) {
	if (left_interval_border > right_interval_border) {
		double tmp = left_interval_border;
		left_interval_border = right_interval_border;
		right_interval_border = tmp;
	}
	double x_1 = 0, x_2 = 0;
	while (fabs(right_interval_border - left_interval_border) > CONST_EPS) {
		x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		if (n_order_deriv(x_1, order) >= n_order_deriv(x_2, order))
			left_interval_border = x_1, x_1 = x_2, x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		else
			right_interval_border = x_2, x_2 = x_1, x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
	}
	return (left_interval_border + right_interval_border) / 2;
}

double Function::abs_argmin_golden_section(double left_interval_border, double right_interval_border) {
	if (left_interval_border > right_interval_border) {
		double tmp = left_interval_border;
		left_interval_border = right_interval_border;
		right_interval_border = tmp;
	}
	double x_1 = 0, x_2 = 0;
	while (fabs(right_interval_border - left_interval_border) > CONST_EPS) {
		x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		if (fabs((*this)(x_1)) >= fabs((*this)(x_2)))
			left_interval_border = x_1, x_1 = x_2, x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		else
			right_interval_border = x_2, x_2 = x_1, x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
	}
	return (left_interval_border + right_interval_border) / 2;
}

double Function::abs_n_deriv_argmin_golden_section(unsigned int order, double left_interval_border, double right_interval_border) {
	if (left_interval_border > right_interval_border) {
		double tmp = left_interval_border;
		left_interval_border = right_interval_border;
		right_interval_border = tmp;
	}
	double x_1 = 0, x_2 = 0;
	while (fabs(right_interval_border - left_interval_border) > CONST_EPS) {
		x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		if (fabs(n_order_deriv(x_1, order)) >= fabs(n_order_deriv(x_2, order)))
			left_interval_border = x_1, x_1 = x_2, x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		else
			right_interval_border = x_2, x_2 = x_1, x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
	}
	return (left_interval_border + right_interval_border) / 2;
}

double Function::argmax_golden_section(double left_interval_border, double right_interval_border) {
	if (left_interval_border > right_interval_border) {
		double tmp = left_interval_border;
		left_interval_border = right_interval_border;
		right_interval_border = tmp;
	}
	double x_1 = 0, x_2 = 0;
	while (fabs(right_interval_border - left_interval_border) > CONST_EPS) {
		x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		if ((*this)(x_1) <= (*this)(x_2))
			left_interval_border = x_1, x_1 = x_2, x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		else
			right_interval_border = x_2, x_2 = x_1, x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
	}
	return (left_interval_border + right_interval_border) / 2;
}

double Function::n_deriv_argmax_golden_section(unsigned int order, double left_interval_border, double right_interval_border) {
	if (left_interval_border > right_interval_border) {
		double tmp = left_interval_border;
		left_interval_border = right_interval_border;
		right_interval_border = tmp;
	}
	double x_1 = 0, x_2 = 0;
	while (fabs(right_interval_border - left_interval_border) > CONST_EPS) {
		x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		if (n_order_deriv(x_1, order) <= n_order_deriv(x_2, order))
			left_interval_border = x_1, x_1 = x_2, x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		else
			right_interval_border = x_2, x_2 = x_1, x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
	}
	return (left_interval_border + right_interval_border) / 2;
}

double Function::abs_argmax_golden_section(double left_interval_border, double right_interval_border) {
	if (left_interval_border > right_interval_border) {
		double tmp = left_interval_border;
		left_interval_border = right_interval_border;
		right_interval_border = tmp;
	}
	double x_1 = 0, x_2 = 0;
	while (fabs(right_interval_border - left_interval_border) > CONST_EPS) {
		x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		if (fabs((*this)(x_1)) <= fabs((*this)(x_2)))
			left_interval_border = x_1, x_1 = x_2, x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		else
			right_interval_border = x_2, x_2 = x_1, x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
	}
	return (left_interval_border + right_interval_border) / 2;
}

double Function::abs_n_deriv_argmax_golden_section(unsigned int order, double left_interval_border, double right_interval_border) {
	if (left_interval_border > right_interval_border) {
		double tmp = left_interval_border;
		left_interval_border = right_interval_border;
		right_interval_border = tmp;
	}
	double x_1 = 0, x_2 = 0;
	while (fabs(right_interval_border - left_interval_border) > CONST_EPS) {
		x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		if (fabs(n_order_deriv(x_1, order)) <= fabs(n_order_deriv(x_2, order)))
			left_interval_border = x_1, x_1 = x_2, x_2 = left_interval_border + (right_interval_border - left_interval_border) / GOLDEN_SECTION;
		else
			right_interval_border = x_2, x_2 = x_1, x_1 = right_interval_border - (right_interval_border - left_interval_border) / GOLDEN_SECTION;
	}
	return (left_interval_border + right_interval_border) / 2;
}

short int Function::signum(double left_interval_border, double right_interval_border, double step = 0.01) {
	bool pos_sgn = false, neg_sgn = false;
	for (double x = left_interval_border; x <= right_interval_border; x += step) {
		if ((*this)(x) > CONST_EPS)
			pos_sgn = true;
		else if ((*this)(x) < -CONST_EPS)
			neg_sgn = true;
	}
	return (pos_sgn && neg_sgn) ? 2 : (pos_sgn ? 1 : (neg_sgn ? -1 : 0));
}

short int Function::n_deriv_signum(unsigned int order, double left_interval_border, double right_interval_border, double step = 0.01) {
	bool pos_sgn = false, neg_sgn = false;
	for (double x = left_interval_border; x <= right_interval_border; x += step) {
		if (n_order_deriv(x, order) > CONST_EPS)
			pos_sgn = true;
		else if (n_order_deriv(x, order) < -CONST_EPS)
			neg_sgn = true;
	}
	return (pos_sgn && neg_sgn) ? 2 : (pos_sgn ? 1 : (neg_sgn ? -1 : 0));
}



void Function::equation_sol_bisection(double func_val, double left_interval_border, double right_interval_border, double &solve, int &iters, double &abs_err, int iters_lim_num, double eps) {
	if (((*this)(left_interval_border) - func_val) * ((*this)(right_interval_border) - func_val) >= 0)
		throw 1;
	solve = (left_interval_border + right_interval_border) / 2;
	iters = 0;
	while (iters < iters_lim_num && (right_interval_border - left_interval_border) / 2 > eps) {
		if (((*this)(solve) - func_val) * ((*this)(left_interval_border) - func_val) < 0)
			right_interval_border = (left_interval_border + right_interval_border) / 2;
		else
			left_interval_border = (left_interval_border + right_interval_border) / 2;
		solve = (left_interval_border + right_interval_border) / 2;
		iters++;
	}
	abs_err = (right_interval_border - left_interval_border) / 2;
}

void Function::equation_sol_iteration(double func_val, double left_interval_border, double right_interval_border, double &solve, int &iters, double &abs_err, int iters_lim_num, double eps) {
	if (((*this)(left_interval_border) - func_val) * ((*this)(right_interval_border) - func_val) > 0)
		throw 1;
	else if (n_deriv_signum(1, left_interval_border, right_interval_border) == 2)
		throw 2;
	double iter_koeff = 2 * n_deriv_signum(1, left_interval_border, right_interval_border) /
		(fabs(n_order_deriv(abs_n_deriv_argmax_golden_section(1, left_interval_border, right_interval_border), 1)) +
			fabs(n_order_deriv(abs_n_deriv_argmin_golden_section(1, left_interval_border, right_interval_border), 1))),
		abs_err_koeff = fabs(n_order_deriv(abs_n_deriv_argmax_golden_section(1, left_interval_border, right_interval_border), 1)) * fabs(iter_koeff) - 1;
	solve = (left_interval_border + right_interval_border) / 2;
	do {
		abs_err = fabs(((*this)(solve) - func_val) * iter_koeff * abs_err_koeff / (1 - abs_err_koeff));
		solve = solve - ((*this)(solve) - func_val) * iter_koeff;
		iters++;
	} while (fabs((*this)(solve) - func_val) > eps && iters < iters_lim_num && abs_err > eps);
}

void Function::equation_sol_chord(double func_val, double left_interval_border, double right_interval_border, double &solve, int &iters, double &abs_err, int iters_lim_num, double eps) {
	if (n_deriv_signum(1, left_interval_border, right_interval_border) == 2 || n_deriv_signum(2, left_interval_border, right_interval_border) == 2)
		throw 1;
	iters = 0;
	solve = ((*this)(right_interval_border) - func_val) * n_deriv_signum(2, left_interval_border, right_interval_border) > 0 ? left_interval_border : right_interval_border;
	if (fabs(solve - left_interval_border) < CONST_EPS)
		do
		{
			solve = right_interval_border - ((*this)(right_interval_border) - func_val) * (right_interval_border - solve) / ((*this)(right_interval_border) - (*this)(solve));
			abs_err = fabs(right_interval_border - ((*this)(right_interval_border) - func_val) * (right_interval_border - solve) / ((*this)(right_interval_border) - (*this)(solve)) - solve);
			iters++;
		} while (fabs((*this)(solve) - func_val) > eps && iters < iters_lim_num && abs_err > eps);
	else
		do
		{
			solve = left_interval_border - ((*this)(left_interval_border) - func_val) * (solve - left_interval_border) / ((*this)(solve) - (*this)(left_interval_border));
			abs_err = fabs(left_interval_border - ((*this)(left_interval_border) - func_val) * (solve - left_interval_border) / ((*this)(solve) - (*this)(left_interval_border)) - solve);
			iters++;
		} while (fabs((*this)(solve) - func_val) > eps && iters < iters_lim_num && abs_err > eps);
}

void Function::equation_sol_tangent(double func_val, double left_interval_border, double right_interval_border, double &solve, int &iters, double &abs_err, int iters_lim_num, double eps) {
	if (n_deriv_signum(1, left_interval_border, right_interval_border) == 2 || n_deriv_signum(2, left_interval_border, right_interval_border) == 2)
		throw 1;
	solve = ((*this)(right_interval_border) - func_val) * n_order_deriv(right_interval_border, 2) > 0 ? right_interval_border : left_interval_border;
	iters = 0;
	do {
		solve = solve - ((*this)(solve) - func_val) / n_order_deriv(solve, 1);
		abs_err = fabs(((*this)(solve) - func_val) / n_order_deriv(solve, 1));
		iters++;
	} while (fabs((*this)(solve) - func_val) > eps && iters < iters_lim_num && abs_err > eps);
}

void Function::integral_Simpson(double a, double b, double &integral, double &abs_delt) {
	integral = 0;
	double h = 0;
	int number_of_partition = (b - a) / H_INTEGR_CONST;
	if (number_of_partition % 2)
		number_of_partition++;
	h = (b - a) / number_of_partition;
	for (int i = 0; i < number_of_partition / 2; i++)
		integral += this_appropriate_func(a + 2 * i * h) + 4 * this_appropriate_func(a + (2 * i + 1) * h) + this_appropriate_func(a + (2 * i + 2) * h);
	integral *= h / 3;
	abs_delt = (b - a) / 180 * n_order_deriv(abs_n_deriv_argmax_golden_section(4, a, b), 4);
	for (int i = 0; i < 4; i++)
		abs_delt *= h;
}


Function::~Function() {}
