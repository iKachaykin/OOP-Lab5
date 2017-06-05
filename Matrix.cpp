#include "Matrix.h"
#include "SquareMatrix.h"
#include <string>

Matrix::Matrix() {
	rows_num = CONST_ROWS;
	cols_num = CONST_COLS;
	mat = new double*[rows_num];
	for (int i = 0; i < rows_num; i++)
		mat[i] = new double[cols_num];
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			mat[i][j] = 0;
}

Matrix::Matrix(int rows_num, int cols_num) {
	if (rows_num <= 0 || cols_num <= 0)
		throw 1;
	this->rows_num = rows_num;
	this->cols_num = cols_num;
	mat = new double*[rows_num];
	for (int i = 0; i < rows_num; i++)
		mat[i] = new double[cols_num];
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			mat[i][j] = 0;
}

Matrix::Matrix(double default_val, int rows_num, int cols_num) {
	if (rows_num <= 0 || cols_num <= 0)
		throw 1;
	this->rows_num = rows_num;
	this->cols_num = cols_num;
	mat = new double*[rows_num];
	for (int i = 0; i < rows_num; i++)
		mat[i] = new double[cols_num];
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			mat[i][j] = default_val;
}

Matrix::Matrix(double **mat, int rows_num, int cols_num) {
	if (rows_num <= 0 || cols_num <= 0)
		throw 1;
	this->rows_num = rows_num;
	this->cols_num = cols_num;
	this->mat = new double*[rows_num];
	for (int i = 0; i < rows_num; i++)
		this->mat[i] = new double[cols_num];
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			this->mat[i][j] = mat[i][j];
}

Matrix::Matrix(const MathVector &vector) {
	rows_num = 1;
	cols_num = vector.dim;
	mat = new double*[1];
	mat[0] = new double[cols_num];
	for (int i = 0; i < cols_num; i++)
		mat[0][i] = vector.vect[i];
}

Matrix::Matrix(const Matrix &other) {
	rows_num = other.rows_num;
	cols_num = other.cols_num;
	mat = new double*[rows_num];
	for (int i = 0; i < rows_num; i++)
		mat[i] = new double[cols_num];
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			mat[i][j] = other.mat[i][j];
}

Matrix & Matrix::operator=(const MathVector &vector)
{
	for (int i = 0; i < rows_num; i++)
		delete[] mat[i];
	delete[] mat;
	rows_num = 1;
	cols_num = vector.dim;
	mat = new double*[1];
	mat[0] = new double[cols_num];
	for (int i = 0; i < cols_num; i++)
		mat[0][i] = vector.vect[i];
	return *this;
}

Matrix & Matrix::operator=(const Matrix &other) {
	for (int i = 0; i < rows_num; i++)
		delete[] mat[i];
	delete[] mat;
	rows_num = other.rows_num;
	cols_num = other.cols_num;
	mat = new double*[rows_num];
	for (int i = 0; i < rows_num; i++)
		mat[i] = new double[cols_num];
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			mat[i][j] = other.mat[i][j];
	return *this;
}

bool Matrix::operator==(const Matrix &other) {
	if (rows_num != other.rows_num || cols_num != other.cols_num)
		return false;
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			if (fabs(mat[i][j] - other.mat[i][j]) >= 0)
				return false;
	return true;
}

bool Matrix::operator!=(const Matrix &other) {
	if (rows_num != other.rows_num || cols_num != other.cols_num)
		return true;
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			if (fabs(mat[i][j] - other.mat[i][j]) >= 0)
				return true;
	return false;
}

double * Matrix::operator[](int index) {
	if (index < 0 || index >= rows_num)
		throw 1;
	return mat[index];
}

double Matrix::min_matrix_norm() const {
	if (cube_norm() <= oct_norm() && cube_norm() <= Euclid_norm())
		return cube_norm();
	else if (oct_norm() <= cube_norm() && oct_norm() <= Euclid_norm())
		return oct_norm();
	else 
		return Euclid_norm();
}

double Matrix::cube_norm() const {
	double norm = 0, hyp_max = 0;
	for (int j = 0; j < cols_num; j++)
		norm += fabs(mat[0][j]);
	for (int i = 1; i < rows_num; i++) {
		hyp_max = 0;
		for (int j = 0; j < cols_num; j++)
			hyp_max += fabs(mat[i][j]);
		if (hyp_max > norm)
			norm = hyp_max;
	}
	return norm;
}

double Matrix::oct_norm() const {
	double norm = 0, hyp_max = 0;
	for (int i = 0; i < rows_num; i++)
		norm += fabs(mat[i][0]);
	for (int j = 1; j < cols_num; j++) {
		hyp_max = 0;
		for (int i = 0; i < rows_num; i++)
			hyp_max += fabs(mat[i][j]);
		if (hyp_max > norm)
			norm = hyp_max;
	}
	return norm;
}

double Matrix::Euclid_norm() const {
	double norm = 0;
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			norm += mat[i][j] * mat[i][j];
	return sqrt(norm);
}

std::string Matrix::to_string() const {
	std::string result = rows_num == 1 ? "(" : "/";
	int to_ins_i_val;
	double to_ins_d_val;
	for (int i = 0; i < rows_num; i++) {
		if (i != 0 && i != rows_num - 1)
			result += "|";
		else if (i == rows_num - 1)
			result += "\\";
		for (int j = 0; j < cols_num; j++) {
			to_ins_i_val = round(mat[i][j]);
			to_ins_d_val = mat[i][j];
			result += (fabs(to_ins_i_val - to_ins_d_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val)) + (j != cols_num - 1 ? "\t\t" : "\t");
		}
		if (i != 0 && i != rows_num - 1)
			result += "|\n";
		else if (i == 0)
			result += "\\\n";
	}
	result += rows_num == 1 ? ")\n\n" : "/\n\n";
	return result;
}

int Matrix::rows() const {
	return rows_num;
}

int Matrix::cols() const {
	return cols_num;
}

void Matrix::transpose() {
	double **t_matrix = new double*[cols_num];
	for (int i = 0; i < cols_num; i++)
		t_matrix[i] = new double[rows_num];
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			t_matrix[j][i] = mat[i][j];
	for (int i = 0; i < rows_num; i++)
		delete[] mat[i];
	delete[] mat;
	mat = t_matrix;
}

Matrix Matrix::create_transposed() const {
	Matrix result(cols_num, rows_num);
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			result.mat[j][i] = mat[i][j];
	return result;
}

Matrix Matrix::operator+(const Matrix &other) {
	if (rows_num != other.rows_num || cols_num != other.cols_num)
		throw 1;
	Matrix result(rows_num, cols_num);
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			result.mat[i][j] = mat[i][j] + other.mat[i][j];
	return result;
}

Matrix Matrix::operator-(const Matrix &other) {
	if (rows_num != other.rows_num || cols_num != other.cols_num)
		throw 1;
	Matrix result(rows_num, cols_num);
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			result.mat[i][j] = mat[i][j] - other.mat[i][j];
	return result;
}

Matrix Matrix::operator*(const Matrix &other) {
	if (cols_num != other.rows_num)
		throw 1;
	Matrix result(rows_num, other.cols_num);
	for (int i = 0; i < result.rows_num; i++)
		for (int j = 0; j < result.cols_num; j++)
			for (int sum_index = 0; sum_index < cols_num; sum_index++)
				result.mat[i][j] += mat[i][sum_index] * other.mat[sum_index][j];
	return result;
}

Matrix Matrix::operator*(double koeff) {
	Matrix result(*this);
	for (int i = 0; i < result.rows_num; i++)
		for (int j = 0; j < result.cols_num; j++)
			result.mat[i][j] *= koeff;
	return result;
}

Matrix Matrix::operator/(double koeff) {
	Matrix result(*this);
	for (int i = 0; i < result.rows_num; i++)
		for (int j = 0; j < result.cols_num; j++)
			result.mat[i][j] /= koeff;
	return result;
}

//Matrix Matrix::operator-() {
//	return *this * -1;
//}

MathVector Matrix::operator*(const MathVector &right_vector) {
	if (cols_num != right_vector.dimension())
		throw 1;
	MathVector result_vector(rows_num);
	for (int vect_i = 0; vect_i < rows_num; vect_i++)
		for (int sum_i = 0; sum_i < cols_num; sum_i++)
			result_vector.vect[vect_i] += mat[vect_i][sum_i] * right_vector.vect[sum_i];
	return result_vector;
}

Matrix & Matrix::operator+=(const Matrix &other) {
	if (rows_num != other.rows_num || cols_num != other.cols_num)
		throw 1;
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			mat[i][j] += other.mat[i][j];
	return *this;
}

Matrix & Matrix::operator-=(const Matrix &other) {
	if (rows_num != other.rows_num || cols_num != other.cols_num)
		throw 1;
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			mat[i][j] -= other.mat[i][j];
	return *this;
}

Matrix & Matrix::operator*=(const Matrix &other) {
	*this = *this * other;
	return *this;
}

Matrix & Matrix::operator*=(double koeff) {
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			mat[i][j] *= koeff;
	return *this;
}

Matrix & Matrix::operator/=(double koeff) {
	for (int i = 0; i < rows_num; i++)
		for (int j = 0; j < cols_num; j++)
			mat[i][j] /= koeff;
	return *this;
}

Matrix::operator SquareMatrix() {
	if (rows_num != cols_num)
		throw 1;
	return SquareMatrix(mat, rows_num);
}

Matrix::operator MathVector() {
	return MathVector(*this);
}


Matrix::~Matrix() {
	for (int i = 0; i < rows_num; i++)
		delete[] mat[i];
	delete[] mat;
}

Matrix operator*(double koeff, const Matrix &right) {
	Matrix result(right);
	for (int i = 0; i < result.rows_num; i++)
		for (int j = 0; j < result.cols_num; j++)
			result.mat[i][j] *= koeff;
	return result;
}
