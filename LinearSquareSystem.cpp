#include "LinearSquareSystem.h"
#include <string>

LinearSquareSystem::LinearSquareSystem(SquareMatrix system_matrix, MathVector free_koeffs_vector) {
	if (system_matrix.dimension() != free_koeffs_vector.dimension())
		throw 1;
	this->system_matrix = system_matrix;
	this->free_koeffs_vector = free_koeffs_vector;
	sys_dim = system_matrix.dimension();
	is_determined = fabs(system_matrix.determinant()) >= EPS;
	if (is_determined)
		optimal_solution_method();
}

LinearSquareSystem::LinearSquareSystem(double **extended_system_matrix, int sys_dim) {
	this->sys_dim = sys_dim;
	system_matrix = SquareMatrix(extended_system_matrix, sys_dim);
	free_koeffs_vector = MathVector(sys_dim);
	for (int i = 0; i < sys_dim; i++)
		free_koeffs_vector[i] = extended_system_matrix[i][sys_dim];
	is_determined = fabs(system_matrix.determinant()) >= EPS;
	if (is_determined)
		optimal_solution_method();
}

LinearSquareSystem::LinearSquareSystem(double **system_matrix, double *free_koeffs_vector, int sys_dim) {
	this->sys_dim = sys_dim;
	this->system_matrix = SquareMatrix(system_matrix, sys_dim);
	this->free_koeffs_vector = MathVector(free_koeffs_vector, sys_dim);
	is_determined = fabs(this->system_matrix.determinant()) >= EPS;
	if (is_determined)
		optimal_solution_method();
}

std::string LinearSquareSystem::to_string(){
	std::string result = "Заданная СЛАУ:\n";
	double to_ins_d_val;
	int to_ins_i_val;
	for (int i = 0, first_non_zero_i; i < sys_dim; i++) {
		first_non_zero_i = -1;
		for (int j = 0; j < sys_dim; j++)
			if (fabs(system_matrix[i][j]) > EPS) {
				first_non_zero_i = j;
				break;
			}
		for (int j = 0; j < sys_dim && first_non_zero_i >= 0; j++) {
			to_ins_d_val = system_matrix[i][j];
			to_ins_i_val = round(system_matrix[i][j]);
			if (j == first_non_zero_i)
				result += (fabs(system_matrix[i][j] - 1) < EPS ? "" : fabs(system_matrix[i][j] + 1) < EPS ? "-" : (fabs(to_ins_d_val - to_ins_i_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val) + " * ")) 
					+ "x" + std::to_string(j + 1);
			else if(fabs(system_matrix[i][j] - 1) < EPS)
				result += " + x" + std::to_string(j + 1);
			else if (fabs(system_matrix[i][j] + 1) < EPS)
				result += " - x" + std::to_string(j + 1);
			else if (system_matrix[i][j] > EPS)
				result += " + " + (fabs(to_ins_d_val - to_ins_i_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val) + " * ") + "x" + std::to_string(j + 1);
			else if (system_matrix[i][j] < -EPS)
				result += " - " + (fabs(to_ins_d_val - to_ins_i_val) < EPS ? std::to_string(-to_ins_i_val) : std::to_string(-to_ins_d_val) + " * ") + "x" + std::to_string(j + 1);
		}
		to_ins_d_val = free_koeffs_vector[i];
		to_ins_i_val = round(free_koeffs_vector[i]);
		result += " = " + (fabs(to_ins_d_val - to_ins_i_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val)) + "\n\n";
	}
	result += "Решение:";
	if (!is_determined)
		result += " отсутствует или не единственное.";
	else
		for (int i = 0; i < sys_dim; i++) {
			to_ins_d_val = solution_vector[i];
			to_ins_i_val = round(solution_vector[i]);
			result += "\n\nx" + std::to_string(i + 1) + " = " + (fabs(to_ins_d_val - to_ins_i_val) < EPS ? std::to_string(to_ins_i_val) : std::to_string(to_ins_d_val)) + ";";
		}
	return result += "\n\nПогрешность вычисление не превышает: " + std::to_string(err) + "\n\n";
}

MathVector LinearSquareSystem::solution() const {
	if (!is_determined)
		throw 1;
	return solution_vector;
}

bool LinearSquareSystem::determined() const {
	return is_determined;
}

int LinearSquareSystem::dimension() const {
	return sys_dim;
}

double LinearSquareSystem::error() const {
	return err;
}

MathVector LinearSquareSystem::free_vector() {
	return free_koeffs_vector;
}

SquareMatrix LinearSquareSystem::matrix() {
	return system_matrix;
}

Matrix LinearSquareSystem::extended_matrix() {
	Matrix result(sys_dim, sys_dim + 1);
	for (int i = 0; i < sys_dim; i++)
		for (int j = 0; j < sys_dim; j++)
			result[i][j] = system_matrix[i][j];
	for (int i = 0; i < sys_dim; i++)
		result[i][sys_dim] = free_koeffs_vector[i];
	return result;
}

Matrix LinearSquareSystem::swaped_system_equations_extended_matrix() {
	double tmp = 0;
	Matrix result(extended_matrix());
	for (int main_diag_i = 0, abs_lower_max_i; main_diag_i < result.rows(); main_diag_i++) {
		abs_lower_max_i = main_diag_i;
		for (int i = main_diag_i + 1; i < result.rows(); i++)
			if (fabs(result[i][main_diag_i]) > fabs(result[abs_lower_max_i][main_diag_i]))
				abs_lower_max_i = i;
		if (main_diag_i != abs_lower_max_i)
			for (int j = 0; j <  result.cols(); j++) {
				tmp = result[main_diag_i][j];
				result[main_diag_i][j] = result[abs_lower_max_i][j];
				result[abs_lower_max_i][j] = tmp;
			}

	}
	return result;
}

void LinearSquareSystem::optimal_solution_method() {
	if (!is_determined)
		throw 1;
	if (!Seidel_solution_method()) {
		Gauss_solution_method();
		err = 0;
	}
}

void LinearSquareSystem::matrix_method() {
	if (!is_determined)
		throw 1;
	solution_vector = system_matrix.create_reversed() * free_koeffs_vector;
}

void LinearSquareSystem::Cramer_method() {
	if (!is_determined)
		throw 1;
	solution_vector = MathVector(sys_dim);
	double tmp_det = system_matrix.determinant();
	SquareMatrix tmp_system_matrix(system_matrix);
	for (int j = 0; j < sys_dim; j++) {
		tmp_system_matrix = system_matrix;
		for (int i = 0; i < sys_dim; i++)
			tmp_system_matrix[i][j] = free_koeffs_vector[i];
		solution_vector[j] = tmp_system_matrix.determinant() / tmp_det;
	}
}

void LinearSquareSystem::Gauss_solution_method() {
	if (!is_determined)
		throw 1;
	double tmp = 0;
	Matrix tmp_extended_matrix(extended_matrix());
	for (int main_diag_i = 0, abs_lower_max_i; main_diag_i < sys_dim; main_diag_i++, abs_lower_max_i = main_diag_i) {
		if (fabs(tmp_extended_matrix[main_diag_i][main_diag_i]) < EPS) {
			abs_lower_max_i = main_diag_i;
			for (int i = main_diag_i + 1; i < sys_dim; i++)
				if (fabs(tmp_extended_matrix[i][main_diag_i]) > fabs(tmp_extended_matrix[abs_lower_max_i][main_diag_i]))
					abs_lower_max_i = i;
			if (main_diag_i != abs_lower_max_i)
				for (int j = 0; j < sys_dim + 1; j++) {
					tmp = tmp_extended_matrix[main_diag_i][j];
					tmp_extended_matrix[main_diag_i][j] = tmp_extended_matrix[abs_lower_max_i][j];
					tmp_extended_matrix[abs_lower_max_i][j] = tmp;
				}
		}
		if (fabs(tmp_extended_matrix[main_diag_i][main_diag_i]) < EPS)
			throw 1;
		for (int i = main_diag_i + 1; i < sys_dim; i++) {
			tmp = tmp_extended_matrix[i][main_diag_i];
			for (int j = main_diag_i; j < sys_dim + 1; j++)
				tmp_extended_matrix[i][j] -= tmp_extended_matrix[main_diag_i][j] / tmp_extended_matrix[main_diag_i][main_diag_i] * tmp;
		}
	}
	for (int sol_i = sys_dim - 1; sol_i >= 0; sol_i--) {
		solution_vector[sol_i] = tmp_extended_matrix[sol_i][sys_dim];
		for (int i = sys_dim - 1; i > sol_i; i--)
			solution_vector[sol_i] -= solution_vector[i] * tmp_extended_matrix[sol_i][i];
		solution_vector[sol_i] /= tmp_extended_matrix[sol_i][sol_i];
	}
}

void LinearSquareSystem::LU_decomposition_solution_method() {
	if (!is_determined || !system_matrix.non_zero_corner_minors())
		throw 1;
	SquareMatrix lower_matrix(sys_dim), upper_matrix(sys_dim); 
	MathVector tmp_vector(sys_dim); 
	for (int i = 0; i < sys_dim; i++) {
		tmp_vector[i] = 0;
		for (int j = 0; j < sys_dim; j++) {
			lower_matrix[i][j] = 0;
			if (i != j)
				upper_matrix[i][j] = 0;
			else
				upper_matrix[i][j] = 1;
		}
	}
	for (int i = 0; i < sys_dim; i++) {
		for (int j = i; j < sys_dim; j++) {
			lower_matrix[j][i] = system_matrix[j][i];
			for (int sum_index = 0; sum_index < i; sum_index++)
				lower_matrix[j][i] -= lower_matrix[j][sum_index] * upper_matrix[sum_index][i];
		}
		for (int j = i + 1; j < sys_dim; j++) {
			upper_matrix[i][j] = system_matrix[i][j];
			for (int sum_index = 0; sum_index < i; sum_index++)
				upper_matrix[i][j] -= lower_matrix[i][sum_index] * upper_matrix[sum_index][j];
			upper_matrix[i][j] /= lower_matrix[i][i];
		}
	}
	for (int i = 0; i < sys_dim; i++) {
		tmp_vector[i] = free_koeffs_vector[i];
		for (int j = 0; j < i; j++)
			tmp_vector[i] -= tmp_vector[j] * lower_matrix[i][j];
		tmp_vector[i] /= lower_matrix[i][i];
	}
	for (int i = sys_dim - 1; i >= 0; i--) {
		solution_vector[i] = tmp_vector[i];
		for (int j = sys_dim - 1; j > i; j--)
			solution_vector[i] -= solution_vector[j] * upper_matrix[i][j];
	}
}

void LinearSquareSystem::square_root_solution_method() {
	if (!is_determined || !system_matrix.non_zero_corner_minors() || !system_matrix.symmetry())
		throw 1;
	SquareMatrix sqrt_matrix(sys_dim);
	MathVector tmp_vector(sys_dim);
	int *sgn_vector = new int[sys_dim];
	for (int i = 0; i < sys_dim; i++)
		for (int j = 0; j < sys_dim; j++)
			sqrt_matrix[i][j] = 0;
	for (int j = 0; j < sys_dim; j++) {
		for (int i = 0; i < j; i++) {
			sqrt_matrix[i][j] = system_matrix[i][j];
			for (int k = 0; k < i; k++)
				sqrt_matrix[i][j] -= sqrt_matrix[k][i] * sqrt_matrix[k][j] * sgn_vector[k];
			sqrt_matrix[i][j] /= sqrt_matrix[i][i] * sgn_vector[i];
		}
		sqrt_matrix[j][j] = system_matrix[j][j];
		for (int k = 0; k < j; k++)
			sqrt_matrix[j][j] -= sgn_vector[k] * sqrt_matrix[k][j] * sqrt_matrix[k][j];
		sgn_vector[j] = sqrt_matrix[j][j];
		sqrt_matrix[j][j] = sqrt(fabs(sqrt_matrix[j][j]));
		sgn_vector[j] = sgn_vector[j] > 0 ? 1 : sgn_vector[j] < 0 ? -1 : 0;
	}
	for (int i = 0; i < sys_dim; i++) {
		tmp_vector[i] = free_koeffs_vector[i];
		for (int j = 0; j < i; j++)
			tmp_vector[i] -= tmp_vector[j] * sqrt_matrix[j][i];
		tmp_vector[i] /= sqrt_matrix[i][i];
	}
	for (int i = sys_dim - 1; i >= 0; i--) {
		solution_vector[i] = tmp_vector[i];
		for (int j = sys_dim - 1; j > i; j--)
			solution_vector[i] -= solution_vector[j] * sqrt_matrix[i][j] * sgn_vector[i];
		solution_vector[i] /= sqrt_matrix[i][i] * sgn_vector[i];
	}
	delete[] sgn_vector;
}

bool LinearSquareSystem::simple_iteration_solution_method() {
	if (!is_determined)
		throw 1;
	bool stop = false;
	unsigned short int norm_num = 0;
	double tmp = 0; 
	MathVector betta_vector(sys_dim), tmp_vector(sys_dim);
	SquareMatrix alpha_matrix(sys_dim);	
	Matrix c_matrix = system_matrix.diagonal_dominating() ? extended_matrix() : swaped_system_equations_extended_matrix();
	for (int i = 0; i < sys_dim; i++) {
		for (int j = 0; j < sys_dim; j++) {
			if (i == j)
				alpha_matrix[i][j] = 0;
			else
				alpha_matrix[i][j] = -c_matrix[i][j] / c_matrix[i][i];
		}
		betta_vector[i] = free_koeffs_vector[i] / c_matrix[i][i];
	}
	if (alpha_matrix.min_matrix_norm() >= 1)
		return false;
	for (int i = 0; i < sys_dim; i++)
		solution_vector[i] = betta_vector[i];
	if (alpha_matrix.cube_norm() <= alpha_matrix.oct_norm() && alpha_matrix.cube_norm() <= alpha_matrix.Euclid_norm())
		norm_num = 0;
	else if (alpha_matrix.oct_norm() <= alpha_matrix.cube_norm() && alpha_matrix.oct_norm() <= alpha_matrix.Euclid_norm())
		norm_num = 1;
	else
		norm_num = 2;
	do {
		for (int i = 0; i < sys_dim; i++) {
			tmp_vector[i] = betta_vector[i];
			for (int j = 0; j < sys_dim; j++)
				tmp_vector[i] += alpha_matrix[i][j] * solution_vector[j];
		}
		for (int i = 0; i < sys_dim; i++)
			solution_vector[i] -= tmp_vector[i];
		err = alpha_matrix.min_matrix_norm() /
			(1 - alpha_matrix.min_matrix_norm());
		err *= !norm_num ? solution_vector.cube_norm() : norm_num == 1 ? solution_vector.oct_norm() : solution_vector.Euclid_norm();
		if (err < EPS)
			stop = true;
		for (int i = 0; i < sys_dim; i++)
			solution_vector[i] = tmp_vector[i];

	} while (!stop);
	return true;
}

bool LinearSquareSystem::Seidel_solution_method() {
	if (!is_determined)
		throw 1;
	bool stop = false;
	unsigned short int norm_num = 0;
	double tmp = 0;
	MathVector betta_vector(sys_dim), tmp_vector(sys_dim);
	SquareMatrix alpha_matrix(sys_dim);
	Matrix c_matrix = system_matrix.diagonal_dominating() ? extended_matrix() : swaped_system_equations_extended_matrix();
	for (int i = 0; i < sys_dim; i++) {
		for (int j = 0; j < sys_dim; j++) {
			if (i == j)
				alpha_matrix[i][j] = 0;
			else
				alpha_matrix[i][j] = -c_matrix[i][j] / c_matrix[i][i];
		}
		betta_vector[i] = free_koeffs_vector[i] / c_matrix[i][i];
	}
	if (alpha_matrix.min_matrix_norm() >= 1)
		return false;
	for (int i = 0; i < sys_dim; i++)
		solution_vector[i] = betta_vector[i];
	if (alpha_matrix.cube_norm() <= alpha_matrix.oct_norm() && alpha_matrix.cube_norm() <= alpha_matrix.Euclid_norm())
		norm_num = 0;
	else if (alpha_matrix.oct_norm() <= alpha_matrix.cube_norm() && alpha_matrix.oct_norm() <= alpha_matrix.Euclid_norm())
		norm_num = 1;
	else
		norm_num = 2;
	do {
		for (int i = 0; i < sys_dim; i++) {
			tmp_vector[i] = betta_vector[i];
			for (int j = i; j < sys_dim; j++)
				tmp_vector[i] += alpha_matrix[i][j] * solution_vector[j];
			for (int j = 0; j < i; j++)
				tmp_vector[i] += alpha_matrix[i][j] * tmp_vector[j];
		}
		for (int i = 0; i < sys_dim; i++)
			solution_vector[i] -= tmp_vector[i];
		err = alpha_matrix.min_matrix_norm() /
			(1 - alpha_matrix.min_matrix_norm());
		err *= !norm_num ? solution_vector.cube_norm() : norm_num == 1 ? solution_vector.oct_norm() : solution_vector.Euclid_norm();
		if (err < EPS)
			stop = true;
		for (int i = 0; i < sys_dim; i++)
			solution_vector[i] = tmp_vector[i];

	} while (!stop);
	return true;
}

LinearSquareSystem::~LinearSquareSystem() {}
