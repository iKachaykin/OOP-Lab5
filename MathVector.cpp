#include "MathVector.h"
#include "Matrix.h"
#include <cmath>
#include <string>

MathVector::MathVector(){
	dim = CONST_VECTOR_DIM;
	vect = new double[dim];
	for (int i = 0; i < dim; i++)
		vect[i] = 0;
}

MathVector::MathVector(int dim) {
	if (dim <= 0)
		throw 1;
	this->dim = dim;
	vect = new double[dim];
	for (int i = 0; i < dim; i++)
		vect[i] = 0;
}

MathVector::MathVector(double *other_vect, int other_dim) {
	if (other_dim <= 0)
		throw 1;
	dim = other_dim;
	vect = new double[dim];
	for (int i = 0; i < dim; i++)
		vect[i] = other_vect[i];
}

MathVector::MathVector(double default_val, int dim) {
	if (dim <= 0)
		throw 1;
	this->dim = dim;
	vect = new double[dim];
	for (int i = 0; i < dim; i++)
		vect[i] = default_val;
}

MathVector::MathVector(const MathVector &other) {
	dim = other.dim;
	vect = new double[dim];
	for (int i = 0; i < dim; i++)
		vect[i] = other.vect[i];
}

MathVector::MathVector(const Matrix &matrix) {
	if (matrix.rows_num != 1 && matrix.cols_num != 1)
		throw 1;
	else if (matrix.rows_num == 1) {
		dim = matrix.cols_num;
		vect = new double[matrix.cols_num];
		for (int i = 0; i < matrix.cols_num; i++)
			vect[i] = matrix.mat[0][i];
	}
	else {
		dim = matrix.rows_num;
		vect = new double[matrix.rows_num];
		for (int i = 0; i < matrix.rows_num; i++)
			vect[i] = matrix.mat[i][0];
	}
}

MathVector & MathVector::operator=(const MathVector &other) {
	delete[] vect;
	dim = other.dim;
	vect = new double[dim];
	for (int i = 0; i < dim; i++)
		vect[i] = other.vect[i];
	return *this;
}

MathVector & MathVector::operator=(const Matrix &matrix) {
	delete[] vect;
	if (matrix.rows_num != 1 && matrix.cols_num != 1)
		throw 1;
	else if (matrix.rows_num == 1) {
		dim = matrix.cols_num;
		vect = new double[matrix.cols_num];
		for (int i = 0; i < matrix.cols_num; i++)
			vect[i] = matrix.mat[0][i];
	}
	else {
		dim = matrix.rows_num;
		vect = new double[matrix.rows_num];
		for (int i = 0; i < matrix.rows_num; i++)
			vect[i] = matrix.mat[i][0];
	}
	return *this;
}

bool MathVector::operator==(const MathVector &other) {
	if (dim != other.dim)
		return false;
	for (int i = 0; i < dim; i++)
		if (fabs(vect[i] - other.vect[i]) >= EPS)
			return false;
	return true;
}

bool MathVector::operator!=(const MathVector &other) {
	if (dim != other.dim)
		return true;
	for (int i = 0; i < dim; i++)
		if (fabs(vect[i] - other.vect[i]) >= EPS)
			return true;
	return false;
}

double MathVector::cube_norm() const {
	double norm = fabs(vect[0]);
	for (int i = 1; i < dim; i++)
		if (fabs(vect[i]) > norm)
			norm = fabs(vect[i]);
	return norm;
}

double MathVector::oct_norm() const {
	double norm = 0;
	for (int i = 0; i < dim; i++)
		norm += fabs(vect[i]);
	return norm;
}

double MathVector::Euclid_norm() const {
	double norm = 0;
	for (int i = 0; i < dim; i++)
		norm += vect[i] * vect[i];
	return sqrt(norm);
}

std::string MathVector::to_string() const {
	std::string result = "(";
	int to_ins_i_val;
	double to_ins_d_val;
	for (int i = 0; i < dim - 1; i++) {
		to_ins_i_val = round(vect[i]);
		to_ins_d_val = vect[i];
		result += (fabs(to_ins_i_val - to_ins_d_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val)) + "; ";
	}
	to_ins_i_val = round(vect[dim - 1]);
	to_ins_d_val = vect[dim - 1];
	result += (fabs(to_ins_i_val - to_ins_d_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val)) + ")";
	return result;
}

int MathVector::dimension() const {
	return dim;
}

MathVector::operator Matrix() {
	return Matrix(&vect, 1, dim);
}

double MathVector::operator[](int index) const {
	if (index >= 0 && index < dim)
		return vect[index];
	else
		throw 1;
}


MathVector operator*(double koeff, const MathVector &right) {
	MathVector result(right);
	for (int i = 0; i < result.dim; i++)
		result.vect[i] *= koeff;
	return result;
}

double & MathVector::operator[](int index) {
	if (index >= 0 && index < dim)
		return vect[index];
	else
		throw 1;
}

MathVector MathVector::operator+(const MathVector &other) {
	if (dim != other.dim)
		throw 1;
	MathVector result(dim);
	for (int i = 0; i < dim; i++)
		result.vect[i] = vect[i] + other.vect[i];
	return result;
}

MathVector MathVector::operator-(const MathVector &other) {
	if (dim != other.dim)
		throw 1;
	MathVector result(dim);
	for (int i = 0; i < dim; i++)
		result.vect[i] = vect[i] - other.vect[i];
	return result;
}

double MathVector::operator*(const MathVector &other) {
	double comp = 0;
	for (int i = 0; i < dim; i++)
		comp += vect[i] * other.vect[i];
	return comp;
}

MathVector MathVector::operator*(double koeff) {
	MathVector result(*this);
	for (int i = 0; i < dim; i++)
		result.vect[i] *= koeff;
	return result;
}

MathVector MathVector::operator*(const Matrix &right_matrix) {
	if (dim != right_matrix.rows_num)
		throw 1;
	MathVector result_vector(right_matrix.cols_num);
	for (int vect_i = 0; vect_i < result_vector.dim; vect_i++)
		for (int sum_i = 0; sum_i < right_matrix.rows_num; sum_i++)
			result_vector.vect[vect_i] += vect[sum_i] * right_matrix.mat[sum_i][vect_i];
	return result_vector;
}

MathVector MathVector::operator/(double koeff) {
	if (fabs(koeff) < EPS)
		throw 1;
	MathVector result(*this);
	for (int i = 0; i < dim; i++)
		result.vect[i] /= koeff;
	return result;
}

//MathVector MathVector::operator-() {
//	return *this * -1;
//}

MathVector & MathVector::operator+=(const MathVector &other) {
	if (dim != other.dim)
		throw 1;
	for (int i = 0; i < dim; i++)
		vect[i] += other.vect[i];
	return *this;
}

MathVector & MathVector::operator-=(const MathVector &other) {
	if (dim != other.dim)
		throw 1;
	for (int i = 0; i < dim; i++)
		vect[i] -= other.vect[i];
	return *this;
}

MathVector & MathVector::operator*=(double koeff) {
	for (int i = 0; i < dim; i++)
		vect[i] *= koeff;
	return *this;
}

MathVector & MathVector::operator/=(double koeff) {
	if (fabs(koeff) < EPS)
		throw 1;
	for (int i = 0; i < dim; i++)
		vect[i] /= koeff;
	return *this;
}

MathVector::~MathVector() {
	delete[] vect;
}

