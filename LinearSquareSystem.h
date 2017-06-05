#pragma once
#include "SquareMatrix.h"
#include "MathVector.h"
#include "InterfaceLinearSquareSystem.h"
class LinearSquareSystem :
	public InterfaceLinearSquareSystem {
private:
	SquareMatrix system_matrix;
	MathVector free_koeffs_vector, solution_vector;
	double err;
	int sys_dim;
	bool is_determined;
public:
	LinearSquareSystem(SquareMatrix, MathVector);
	LinearSquareSystem(double **, int);
	LinearSquareSystem(double **, double *, int);
	std::string to_string();
	MathVector solution() const;
	bool determined() const;
	int dimension() const;
	double error() const;
	MathVector free_vector();
	SquareMatrix matrix();
	Matrix extended_matrix();
	Matrix swaped_system_equations_extended_matrix();
	void optimal_solution_method();
	void matrix_method();
	void Cramer_method();
	void Gauss_solution_method();
	void LU_decomposition_solution_method();
	void square_root_solution_method();
	bool simple_iteration_solution_method();
	bool Seidel_solution_method();
	~LinearSquareSystem();
};

