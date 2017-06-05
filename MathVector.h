#pragma once
#include "InterfaceVectorSpaceElement.h"
#include <iostream>
#include <vector>
class MathVector;
class Matrix;

MathVector operator*(double, const MathVector &);

class MathVector :
	public InterfaceVectorSpaceElement {
private:
	const int CONST_VECTOR_DIM = 8;
	double *vect;
	int dim;
public:
	friend class Matrix;
	friend MathVector operator*(double, const MathVector &);
	MathVector();
	MathVector(int);
	template<class list>
	MathVector(list);
	MathVector(double*, int);
	MathVector(double, int);
	MathVector(const MathVector &);
	MathVector(const Matrix &);
	MathVector & operator=(const MathVector &);
	MathVector & operator=(const Matrix &);
	bool operator==(const MathVector &);
	bool operator!=(const MathVector &);
	double cube_norm() const;
	double oct_norm() const;
	double Euclid_norm() const;
	std::string to_string() const;
	int dimension() const;
	operator Matrix();
	double operator[](int) const;
	double &operator[](int);
	MathVector operator+(const MathVector &);
	MathVector operator-(const MathVector &);
	double operator*(const MathVector &);
	MathVector operator*(double);
	MathVector operator*(const Matrix &);
	MathVector operator/(double);
	//MathVector operator-();
	MathVector& operator+=(const MathVector &);
	MathVector& operator-=(const MathVector &);
	MathVector& operator*=(double);
	MathVector& operator/=(double);
	~MathVector();
};

template<class list>
inline MathVector::MathVector(list lst) {
	dim = lst.size();
	vect = new double[dim];
	for (int i = 0; i < dim; i++)
		vect[i] = lst.at(i);
}
