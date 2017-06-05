#pragma once
#include "InterfaceMatrix.h"
class Matrix;
class SquareMatrix;

#include "MathVector.h"

Matrix operator*(double, const Matrix&);

class Matrix :
	public InterfaceMatrix {
private:
	const int CONST_ROWS = 3, CONST_COLS = 3;
	int rows_num, cols_num;
protected:
	double **mat;
public:
	friend class MathVector;
	friend Matrix operator*(double, const Matrix&);
	Matrix();
	Matrix(int, int);
	Matrix(double, int, int);
	Matrix(double **, int, int);
	Matrix(const MathVector &);
	Matrix(const Matrix &);
	Matrix& operator=(const MathVector &);
	Matrix& operator=(const Matrix &);
	bool operator==(const Matrix &);
	bool operator!=(const Matrix &);
	double * operator[](int);
	double min_matrix_norm() const;
	double cube_norm() const;
	double oct_norm() const;
	double Euclid_norm() const;
	std::string to_string() const;
	int rows() const;
	int cols() const;
	void transpose();
	Matrix create_transposed() const;
	operator SquareMatrix();
	operator MathVector();
	Matrix operator+(const Matrix &);
	Matrix operator-(const Matrix &);
	Matrix operator*(const Matrix &);
	Matrix operator*(double);
	Matrix operator/(double);
	//Matrix operator-();
	MathVector operator*(const MathVector &);
	Matrix & operator+=(const Matrix &);
	Matrix & operator-=(const Matrix &);
	Matrix & operator*=(const Matrix &);
	Matrix & operator*=(double);
	Matrix & operator/=(double);
	virtual ~Matrix();
};

