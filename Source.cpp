#define _CRT_SECURE_NO_WARNINGS
#include "Polinom.h"
#include "Function.h"
#include <cmath>
#include <ctime>
#include "SquareMatrix.h"
#include "LinearSquareSystem.h"
#include <string>

double func(double x) { return (x - 2) * (x - 1); }

//double normal_density(double x, double shift, double scale) {
//	return 1 / sqrt(2 * pi) / scale * exp(-(x - shift) * (x - shift) / 2 / scale / scale);
//}
//
//double concrete_normal_density(double x) {
//	return normal_density(x, 2, 2);
//}

//double concrete_normal_func(double x) {
//	double res_val = 0, delt;
//	integral_Simpson(concrete_normal_density, -20, x, res_val, delt);
//	return res_val;
//}

using namespace std;

void main() {
	srand(time(0));
	setlocale(LC_ALL, "ru");
	//int power = 3;
	//double **matrix = new double*[power];
	//for (int i = 0; i < power; i++)
	//	matrix[i] = new double[power + 1];
	//for (int i = 0; i < power; i++)
	//	for (int j = i; j < power; j++)
	//		if (i != j)
	//			matrix[i][j] = matrix[j][i] = rand() % 21 - 10;
	//		else
	//			matrix[i][j] = 10;
	//for(int i = 0; i < power; i++)
	//	matrix[i][power] = rand() % 21 - 10;
	//LinearSquareSystem sys(matrix, power);
	//cout << sys;
	//const int sample_size = 20;
	//double sol, delt, iters, uniform_distr_sample[sample_size] 
	//{ 0.0842, 0.0937, 0.1009, 0.1180, 0.1280, 0.1964, 0.3106, 0.3754, 0.6357, 0.6489,
	//  0.6606, 0.7379, 0.7652, 0.8015, 0.8345, 0.8526, 0.8868, 0.9852, 0.9901, 0.9959 }, 
	//normal_distr_sample[sample_size];
	//for (int i = 0; i < sample_size; i++)
	//	equation_sol_bisection(concrete_normal_func, uniform_distr_sample[i], -6, 10, normal_distr_sample[i], iters, delt);
	//cout << "Выборка нормального распределения:\n";
	//for (int i = 0; i < sample_size; i++)
	//	cout << normal_distr_sample[i] << "\n";
	//cout << "Значение функций распределения в данных точках выборки:\n";
	//for (int i = 0; i < sample_size; i++)
	//	cout << concrete_normal_func(normal_distr_sample[i]) << "\n";
	//for (int i = 0; i < power; i++)
	//	delete[] matrix[i];
	//delete[] matrix;
	int iters = 0, power = rand() % 11;
	double sol = 0, delt = 0;
	cout << "Input power: ";
	cin >> power;
	Polinom pol(power);
	MathVector tmp_vect(power + 1);
	for (int i = 0; i <= power; i++) {
		cout << "Input a" << i << ": ";
		cin >> tmp_vect[i];
	}
	pol = Polinom(tmp_vect);
	cout << "Inputed polinom:\n" << pol << "\n";
	cout << "derives:\n";
	for (int i = 0; i <= power; i++)
		cout << "Order: " << i << "\nDeriv:" << pol.create_n_order_derivative(i) << "\n";
	cout << "Number of real solves in interval (-500;500): " << pol.Sturm_method(-500, 0) << "\n";
	system("pause");
}