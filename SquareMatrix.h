#pragma once
#include "Matrix.h"
class SquareMatrix :
	public Matrix {
private:
	int matrix_dim;
public:
	SquareMatrix();
	SquareMatrix(int);
	SquareMatrix(double, int);
	SquareMatrix(double **, int);
	SquareMatrix(const Matrix&);
	SquareMatrix(const SquareMatrix&);
	SquareMatrix& operator=(const Matrix&);
	bool symmetry() const;
	bool non_zero_corner_minors() const;
	bool diagonal_dominating() const;
	int dimension() const;
	SquareMatrix create_reversed() const;
	double determinant() const;
	~SquareMatrix();
};

