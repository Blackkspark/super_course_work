#include <iostream>
#include <cmath>
#include <fstream>
#include <thread>
#include <string>


using namespace std;

double C = 5;
double a = 0.05;

const int iterations_x = 10;
const int iterations_t = 1000;


double approximated_solution(double w_minus_1, double w, double w_plus_1, double x_move, double t_move)
{
	return w + t_move * (pow(w, -2. / 3.) * ((w_minus_1 - 2 * w + w_plus_1) / pow(x_move, 2))
		- (w_plus_1 - w_minus_1) / (3 * x_move * pow(w, 5. / 3.)));
}


double correct_solution(double x, double t)
{
	double mult1 = pow((C - 4 * a * t), (double)3 / 2);
	double mult2 = mult1 - pow(x, 2);

	return mult1 * pow(mult2, (double)-3 / 2);
}


void write_to_file(string name, double matrix[iterations_t][iterations_x], double dx, double dt)
{
	ofstream file(name);
	double x = 0.0;
	double t = 0.0;

	file << '{';
	for (int i = 0; i < iterations_t; ++i)
	{
		x = 0.0;
		for (int j = 0; j < iterations_x; ++j)
		{
			file << '{' << x << ", " << t << ", " << matrix[i][j] << '}' << ',';
			x += dx;
		}
		t += dt;
	}
	file << '}';

	file.close();
}


void calc_correct_matrix(double matrix[iterations_t][iterations_x], double iter_grid_of_x, double iter_grid_of_t, double move_grid_x, double move_grid_t)
{
	for (int i = 0; i < iterations_t; ++i)
	{
		iter_grid_of_x = 0.0;
		for (int j = 0; j < iterations_x; ++j)
		{
			matrix[i][j] = correct_solution(iter_grid_of_x, iter_grid_of_t);
			iter_grid_of_x += move_grid_x;
		}
		iter_grid_of_t += move_grid_t;
	}
}


void calc_partial_derivatives(double matrix[iterations_t][iterations_x], double iter_grid_of_x, double iter_grid_of_t, double move_grid_x, double move_grid_t)
{
	for (int j = 0; j < iterations_x; ++j)
	{
		matrix[0][j] = correct_solution(iter_grid_of_x, 0);
		iter_grid_of_x += move_grid_x;
	}

	for (int i = 0; i < iterations_t; ++i)
	{
		matrix[i][0] = correct_solution(0, iter_grid_of_t);
		matrix[i][iterations_x - 1] = correct_solution(1, iter_grid_of_t);
		iter_grid_of_t += move_grid_t;
	}
}


void successively_calc(double calc_matrix[iterations_t][iterations_x], double correct_matrix[iterations_t][iterations_x], double move_grid_x, double move_grid_t, double &error, int &error_place_i, int &error_place_j)
{
	for (int i = 1; i < iterations_t; ++i)
	{
		for (int j = 1; j < iterations_x - 1; ++j)
		{
			calc_matrix[i][j] = approximated_solution(calc_matrix[i - 1][j - 1], calc_matrix[i - 1][j],
				calc_matrix[i - 1][j + 1], move_grid_x, move_grid_t);
			if (fabs(calc_matrix[i][j] - correct_matrix[i][j]) > error)
			{
				error = fabs(calc_matrix[i][j] - correct_matrix[i][j]);
				error_place_i = i;
				error_place_j = j;
			}
		}
	}
}


void parallel_calc(double calc_matrix[iterations_t][iterations_x], double correct_matrix[iterations_t][iterations_x], double move_grid_x, double move_grid_t, double &error, int &error_place_i, int &error_place_j)
{
	for (int i = 1; i < iterations_t; ++i)
	{
#pragma omp parallel for
		for (int j = 1; j < iterations_x - 1; ++j)
		{
			calc_matrix[i][j] = approximated_solution(calc_matrix[i - 1][j - 1], calc_matrix[i - 1][j],
				calc_matrix[i - 1][j + 1], move_grid_x, move_grid_t);
			if (fabs(calc_matrix[i][j] - correct_matrix[i][j]) > error)
			{
				error = fabs(calc_matrix[i][j] - correct_matrix[i][j]);
				error_place_i = i;
				error_place_j = j;
			}
		}
	}
}



int main()
{
	setlocale(0, "rus");

	double x_iter_grid = 0.0, t_iter_grid = 0.0;
	double x_move = 1. / (iterations_x - 1);
	double t_move = 1. / (iterations_t - 1);


	double correct_matrix[iterations_t][iterations_x];
	calc_correct_matrix(correct_matrix, x_iter_grid, t_iter_grid, x_move, t_move);
	write_to_file("super_correct_result.txt", correct_matrix, x_move, t_move);


	x_iter_grid = 0;
	t_iter_grid = 0;
	double approximated_matrix_sequential[iterations_t][iterations_x];
	calc_partial_derivatives(approximated_matrix_sequential, x_iter_grid, t_iter_grid, x_move, t_move);


	double error_value = 0.0;
	int error_i_pos = 0;
	int error_j_pos = 0;
	successively_calc(approximated_matrix_sequential, correct_matrix, x_move, t_move, error_value, error_i_pos, error_j_pos);
	cout << "Последовательное исчисление:" << endl;
	cout << "Значение абсолютной погрешности: " << error_value << endl;
	cout << "Относительная погрешность: " <<
		fabs((approximated_matrix_sequential[error_i_pos][error_j_pos] - correct_matrix[error_i_pos][error_j_pos])
			/ correct_matrix[error_i_pos][error_j_pos]) * 100 << "%" << endl << endl;
	write_to_file("super_approximated_result_sic.txt", approximated_matrix_sequential, x_move, t_move);


	double approximated_matrix_parallel[iterations_t][iterations_x];
	calc_partial_derivatives(approximated_matrix_parallel, x_iter_grid, t_iter_grid, x_move, t_move);
	error_value = 0.0;
	error_i_pos = 0;
	error_j_pos = 0;
	parallel_calc(approximated_matrix_parallel, correct_matrix, x_move, t_move, error_value, error_i_pos, error_j_pos);
	cout << "Паралельное исчисление:" << endl;
	cout << "Значение абсолютной погрешности: " << error_value << endl;
	cout << "Относительная погрешность: " <<
		fabs((approximated_matrix_parallel[error_i_pos][error_j_pos] - correct_matrix[error_i_pos][error_j_pos])
			/ correct_matrix[error_i_pos][error_j_pos]) * 100 << "%" << endl << endl;
	write_to_file("super_approximated_result_paral.txt", approximated_matrix_parallel, x_move, t_move);


	system("pause");
	return 0;
}