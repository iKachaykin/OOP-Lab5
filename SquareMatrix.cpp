#include "SquareMatrix.h"


SquareMatrix::SquareMatrix() : Matrix() {}

SquareMatrix::SquareMatrix(int matrix_dim) : Matrix(matrix_dim, matrix_dim) {
	this->matrix_dim = matrix_dim;
}

SquareMatrix::SquareMatrix(double default_val, int matrix_dim) : Matrix(default_val, matrix_dim, matrix_dim) {
	this->matrix_dim = matrix_dim;
}

SquareMatrix::SquareMatrix(double **mat, int matrix_dim) : Matrix(mat, matrix_dim, matrix_dim) {
	this->matrix_dim = matrix_dim;
}

SquareMatrix::SquareMatrix(const Matrix &other) : Matrix(other) {
	if (other.cols() != other.rows())
		throw 1;
	matrix_dim = other.rows();
}

SquareMatrix::SquareMatrix(const SquareMatrix &other) : Matrix(other) {
	matrix_dim = other.matrix_dim;
}

SquareMatrix & SquareMatrix::operator=(const Matrix &other) {
	if (other.cols() != other.rows())
		throw 1;
	Matrix tmp(other);
	for (int i = 0; i < matrix_dim; i++)
		delete[] mat[i];
	delete[] mat;
	matrix_dim = tmp.rows();
	mat = new double*[matrix_dim];
	for (int i = 0; i < matrix_dim; i++)
		mat[i] = new double[matrix_dim];
	for (int i = 0; i < matrix_dim; i++)
		for (int j = 0; j < matrix_dim; j++)
			mat[i][j] = tmp[i][j];
	return *this;
}

bool SquareMatrix::symmetry() const {
	for (int i = 0; i < matrix_dim; i++)
		for (int j = i + 1; j < matrix_dim; j++)
			if (fabs(mat[i][j] - mat[j][i]) >= EPS)
				return false;
	return true;
}

bool SquareMatrix::non_zero_corner_minors() const {
	SquareMatrix tmp_matrix(mat, 1);
	for (int i = 1; i < matrix_dim; i++, tmp_matrix = SquareMatrix(mat, i))
		if (fabs(tmp_matrix.determinant()) < EPS)
			return false;
	return true;
}

bool SquareMatrix::diagonal_dominating() const {
	double abs_hyp_max_in_row = 0;
	for (int i = 0; i < matrix_dim; i++) {
		abs_hyp_max_in_row = fabs(mat[i][0]);
		for (int j = 1; j < matrix_dim; j++)
			if (abs_hyp_max_in_row < fabs(mat[i][j]))
				abs_hyp_max_in_row = fabs(mat[i][j]);
		if (fabs(abs_hyp_max_in_row - mat[i][i]) > EPS)
			return false;
	}
	return true;
}

int SquareMatrix::dimension() const {
	return matrix_dim;
}

SquareMatrix SquareMatrix::create_reversed() const {
	double tmp = 0;
	SquareMatrix result(matrix_dim), this_copy(*this);
	for (int i = 0; i < matrix_dim; i++)
		result.mat[i][i] = 1;
	for (int main_diag_i = 0, abs_max_in_col_i; main_diag_i < matrix_dim; main_diag_i++) {
		if (fabs(this_copy.mat[main_diag_i][main_diag_i]) < EPS) {
			abs_max_in_col_i = main_diag_i;
			for (int i = main_diag_i + 1; i < matrix_dim; i++) {
				if (fabs(this_copy.mat[i][main_diag_i]) > fabs(this_copy.mat[abs_max_in_col_i][main_diag_i]))
					abs_max_in_col_i = i;
			}
			if (fabs(this_copy.mat[abs_max_in_col_i][main_diag_i]) < EPS)
				throw 1;
			for (int j = 0; j < matrix_dim; j++) {
				tmp = this_copy.mat[abs_max_in_col_i][j];
				this_copy.mat[abs_max_in_col_i][j] = this_copy.mat[main_diag_i][j];
				this_copy.mat[main_diag_i][j] = tmp;
				tmp = result.mat[abs_max_in_col_i][j];
				result.mat[abs_max_in_col_i][j] = result.mat[main_diag_i][j];
				result.mat[main_diag_i][j] = tmp;
			}
		}
		for (int i = main_diag_i + 1; i < matrix_dim; i++) {
			tmp = this_copy.mat[i][main_diag_i];
			for (int j = main_diag_i; j < matrix_dim; j++)
				this_copy.mat[i][j] -= this_copy.mat[main_diag_i][j] / this_copy.mat[main_diag_i][main_diag_i] * tmp;
			for(int j = 0; j < matrix_dim; j++)
				result.mat[i][j] -= result.mat[main_diag_i][j] / this_copy.mat[main_diag_i][main_diag_i] * tmp;
		}
	}
	for (int main_diag_i = matrix_dim - 1; main_diag_i >= 0; main_diag_i--) {
		for (int i = main_diag_i - 1; i >= 0; i--) {
			tmp = this_copy.mat[i][main_diag_i];
			for (int j = 0; j < matrix_dim; j++) {
				this_copy.mat[i][j] -= this_copy.mat[main_diag_i][j] / this_copy.mat[main_diag_i][main_diag_i] * tmp;
				result.mat[i][j] -= result.mat[main_diag_i][j] / this_copy.mat[main_diag_i][main_diag_i] * tmp;
			}
		}
	}
	for (int i = 0; i < matrix_dim; i++) {
		tmp = this_copy.mat[i][i];
		for (int j = 0; j < matrix_dim; j++) {
			this_copy.mat[i][j] /= tmp;
			result.mat[i][j] /= tmp;
		}
	}
	return result;
}


double SquareMatrix::determinant() const {
	double det = 1, tmp = 0;
	int transp = 0;
	SquareMatrix this_copy(*this);
	for (int main_diag_i = 0, abs_max_in_row_i; main_diag_i < matrix_dim; main_diag_i++) {
		if (fabs(this_copy.mat[main_diag_i][main_diag_i]) < EPS) {
			abs_max_in_row_i = main_diag_i;
			for (int j = main_diag_i + 1; j < matrix_dim; j++) {
				if (fabs(this_copy.mat[main_diag_i][j]) > fabs(this_copy.mat[main_diag_i][abs_max_in_row_i]))
					abs_max_in_row_i = j;
			}
			if (fabs(this_copy.mat[main_diag_i][abs_max_in_row_i]) < EPS)
				return 0;
			transp++;
			for (int i = 0; i < matrix_dim; i++) {
				tmp = this_copy.mat[i][abs_max_in_row_i];
				this_copy.mat[i][abs_max_in_row_i] = this_copy.mat[i][main_diag_i];
				this_copy.mat[i][main_diag_i] = tmp;
			}
		}
		for (int i = main_diag_i + 1; i < matrix_dim; i++) {
			tmp = this_copy.mat[i][main_diag_i];
			for (int j = main_diag_i; j < matrix_dim; j++)
				this_copy.mat[i][j] -= this_copy.mat[main_diag_i][j] / this_copy.mat[main_diag_i][main_diag_i] * tmp;
		}
	}
	for (int main_diag_i = 0; main_diag_i < matrix_dim; main_diag_i++)
		det *= this_copy.mat[main_diag_i][main_diag_i];
	det *= pow(-1, transp);
	return det;
}

SquareMatrix::~SquareMatrix() {}
