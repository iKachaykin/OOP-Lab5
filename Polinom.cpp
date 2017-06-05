#include "Polinom.h"
#include <string>
#include <vector>

Polinom::Polinom(int polinom_pow) : Function(MathVector(1, polinom_pow + 1)) {
	this->polinom_pow = polinom_pow;
	polinom_koeffs = MathVector(1, polinom_pow + 1);
}

Polinom::Polinom(double *polinom_koeffs, int polinom_pow) : Function(MathVector(polinom_koeffs, polinom_pow + 1)){
	//if (!polinom_koeffs[0] && polinom_pow)
	//	throw 1;
	this->polinom_pow = polinom_pow;
	this->polinom_koeffs = MathVector(polinom_koeffs, polinom_pow + 1);
}

Polinom::Polinom(MathVector polinom_koeffs) : Function(polinom_koeffs){
	//if (!polinom_koeffs[0] && polinom_koeffs.dimension() != 1)
	//	throw 1;
	this->polinom_koeffs = polinom_koeffs;
	polinom_pow = polinom_koeffs.dimension() - 1;
}

Polinom::Polinom(const Polinom &other) : Function(other.polinom_koeffs) {
	polinom_koeffs = other.polinom_koeffs;
	polinom_pow = other.polinom_pow;
}

Polinom & Polinom::operator=(const Polinom &other) {
	polinom_pow = other.polinom_pow;
	polinom_koeffs = other.polinom_koeffs;
	return *this;
}

double Polinom::operator()(double x) {
	double result = 0, tmp_comp = 1;
	for (int i = polinom_pow; i >= 0; i--) {
		result += polinom_koeffs[i] * tmp_comp;
		tmp_comp *= x;
	}
	return result;
}

int Polinom::power() const {
	return polinom_pow;
}

std::string Polinom::to_string() const {
	std::string result = "";
	int to_ins_i_val = round(polinom_koeffs[0]);
	double to_ins_d_val = polinom_koeffs[0];
	switch (polinom_pow) {
	case 0:
		result = fabs(to_ins_d_val - to_ins_i_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val);
		break;
	case 1:
		result = (fabs(polinom_koeffs[0] - 1) < EPS ? "" : fabs(polinom_koeffs[0] + 1) < EPS ? "-" : (fabs(to_ins_d_val - to_ins_i_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val) + " * ")) + "x";
		break;
	default:
		result = (fabs(polinom_koeffs[0] - 1) < EPS ? "" : fabs(polinom_koeffs[0] + 1) < EPS ? "-" : (fabs(to_ins_d_val - to_ins_i_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val) + " * ")) +  "x ^ " + std::to_string(polinom_pow);
		break;
	}
	for (int i = 1; i <= polinom_pow; i++) {
		to_ins_i_val = round(polinom_koeffs[i]);
		to_ins_d_val = polinom_koeffs[i];
		if (fabs(polinom_koeffs[i]) < EPS)
			continue;
		else if (fabs(polinom_koeffs[i] - 1) < EPS) {
			result += " + ";
			if (polinom_pow == i)
				result += "1";
		}
		else if (fabs(polinom_koeffs[i] + 1) < EPS) {
			result += " - ";
			if (polinom_pow == i)
				result += "1";
		}
		else if (polinom_koeffs[i] > EPS)
			result += " + " + (fabs(to_ins_d_val - to_ins_i_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val));
		else
			result += " - " + (fabs(to_ins_d_val - to_ins_i_val) < EPS ? std::to_string(-to_ins_i_val) : std::to_string(-to_ins_d_val));
		if (fabs(to_ins_d_val - to_ins_i_val) > EPS && i != polinom_pow)
			result += " * ";
		if (polinom_pow - i == 1)
			result += "x";
		else if (i != polinom_pow)
			result += "x ^ " + std::to_string(polinom_pow - i);
	}
	return result;
}

bool Polinom::operator==(const Polinom &other) {
	return polinom_koeffs == other.polinom_koeffs;
}

bool Polinom::operator!=(const Polinom &other) {
	return polinom_koeffs != other.polinom_koeffs;
}

Polinom Polinom::operator+(const Polinom &other) {
	int max_of_pow = polinom_pow >= other.polinom_pow ? polinom_pow : other.polinom_pow,
		this_index = polinom_pow, other_index = other.polinom_pow;
	Polinom result(max_of_pow);
	for (; other_index >= 0 && this_index >= 0; other_index--, this_index--)
		result.polinom_koeffs[(this_index >= other_index ? this_index : other_index)] = polinom_koeffs[this_index] + other.polinom_koeffs[other_index];
	for (int i = (this_index >= 0 ? this_index : other_index); i >= 0; i--)
		result.polinom_koeffs[i] = (this_index >= 0 ? polinom_koeffs[i] : other.polinom_koeffs[i]);
	return result;
}

Polinom Polinom::operator-(const Polinom &other) {
	int max_of_pow = polinom_pow >= other.polinom_pow ? polinom_pow : other.polinom_pow,
		this_index = polinom_pow, other_index = other.polinom_pow;
	Polinom result(max_of_pow);
	for (; other_index >= 0 && this_index >= 0; other_index--, this_index--)
		result.polinom_koeffs[(this_index >= other_index ? this_index : other_index)] = polinom_koeffs[this_index] - other.polinom_koeffs[other_index];
	for (int i = (this_index >= 0 ? this_index : other_index); i >= 0; i--)
		result.polinom_koeffs[i] = (this_index >= 0 ? polinom_koeffs[i] : -other.polinom_koeffs[i]);
	return result;
}

Polinom Polinom::operator*(const Polinom &other) {
	Polinom result(polinom_pow + other.polinom_pow);
	result.polinom_koeffs = MathVector((double)0, result.polinom_pow + 1);
	for (int i = result.polinom_pow; i >= 0; i--)
		for (int j = polinom_pow; j >= 0; j--)
			for (int k = other.polinom_pow; k >= 0; k--)
				if (j + k == i)
					result.polinom_koeffs[i] += polinom_koeffs[j] * other.polinom_koeffs[k];
	return result;
}

Polinom Polinom::operator/(const Polinom &other) {
	if (other.polinom_pow > polinom_pow)
		return Polinom(MathVector(1));
	std::vector<double> vect;
	MathVector result_v(polinom_pow - other.polinom_pow + 1), tmp_vect(polinom_pow - other.polinom_pow + 1);
	Polinom tmp_polinom(*this);
	for (int i = 0; i <= polinom_pow - other.polinom_pow; i++) {
		result_v[i] = tmp_polinom.polinom_koeffs[0] / other.polinom_koeffs[0];
		tmp_vect = MathVector(polinom_pow - other.polinom_pow + 1 - i);
		tmp_vect[0] = result_v[i];
		tmp_polinom = tmp_polinom - Polinom(tmp_vect) * other;
		vect.clear();
		for (int i = 0; i <= tmp_polinom.polinom_pow; i++)
			vect.push_back(tmp_polinom.polinom_koeffs[i]);
		vect.erase(vect.begin());
		tmp_polinom = Polinom(MathVector(vect));
	}
	return Polinom(result_v);
}

Polinom Polinom::operator%(const Polinom &other) {
	if (other.polinom_pow > polinom_pow)
		return Polinom(MathVector(1));
	std::vector<double> vect;
	MathVector result_v(polinom_pow - other.polinom_pow + 1), tmp_vect(polinom_pow - other.polinom_pow + 1);
	Polinom tmp_polinom(*this);
	for (int i = 0; i <= polinom_pow - other.polinom_pow; i++) {
		result_v[i] = tmp_polinom.polinom_koeffs[0] / other.polinom_koeffs[0];
		tmp_vect = MathVector(polinom_pow - other.polinom_pow + 1 - i);
		tmp_vect[0] = result_v[i];
		tmp_polinom = tmp_polinom - Polinom(tmp_vect) * other;
		if (tmp_polinom.polinom_pow) {
			vect.clear();
			for (int j = 0; j <= tmp_polinom.polinom_pow; j++)
				vect.push_back(tmp_polinom.polinom_koeffs[j]);
			vect.erase(vect.begin());
			tmp_polinom = Polinom(MathVector(vect));
		}
	}
	if (tmp_polinom.polinom_pow) {
		vect.clear();
		for (int j = 0; j <= tmp_polinom.polinom_pow; j++)
			vect.push_back(tmp_polinom.polinom_koeffs[j]);
		while (!*vect.begin() && vect.begin() < vect.end() - 1)
			vect.erase(vect.begin());
		tmp_polinom = Polinom(MathVector(vect));
	}
	return tmp_polinom;
}

Polinom Polinom::operator-() {
	Polinom result(*this);
	result.polinom_koeffs *= -1;
	return result;
}

Polinom & Polinom::operator+=(const Polinom &other) {
	*this = *this + other;
	return *this;
}

Polinom & Polinom::operator-=(const Polinom &other) {
	*this = *this - other;
	return *this;
}

Polinom & Polinom::operator*=(const Polinom &other) {
	*this = *this * other;
	return *this;
}

Polinom & Polinom::operator/=(const Polinom &other) {
	*this = *this / other;
	return *this;
}

int Polinom::Sturm_method(double left_interval_border, double right_interval_border) {
	int left_res = 0, right_res = 0;
	std::vector<Polinom> Sturm_system;
	Sturm_system.push_back(*this);
	Sturm_system.push_back(create_n_order_derivative(1));
	for (int i = 0; Sturm_system[i] % Sturm_system[i + 1] != Polinom(MathVector(1)); i++)
		Sturm_system.push_back(-(Sturm_system[i] % Sturm_system[i + 1]));
	for (int i = 0; i < Sturm_system.size() - 1; i++)
		if (Sturm_system[i](left_interval_border) * Sturm_system[i + 1](left_interval_border) < 0)
			left_res++;
	for (int i = 0; i < Sturm_system.size() - 1; i++)
		if (Sturm_system[i](right_interval_border) * Sturm_system[i + 1](right_interval_border) < 0)
			right_res++;
	return left_res - right_res;
}

Polinom Polinom::create_n_order_derivative(int order) {
	if (order < 0)
		throw 1;
	else if (order > polinom_pow)
		return Polinom(MathVector(1));
	Polinom result(polinom_pow - order);
	for (int i = 0; i <= result.polinom_pow; i++) {
		result.polinom_koeffs[i] = polinom_koeffs[i];
		for (int j = i; j < order + i; j++)
			result.polinom_koeffs[i] *= polinom_pow - j;
	}
	return result;
}

double Polinom::n_order_deriv(double x, unsigned int order) {
	return create_n_order_derivative(order)(x);
}

Polinom::~Polinom() {}

std::ostream & operator<<(std::ostream &os, Polinom &pol) {
	return os << pol.to_string();
}

//Polinom Polinom::operator-() {
//	return (*this) * (double)-1;
//}
