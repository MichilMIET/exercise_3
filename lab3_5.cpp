#include <stdio.h>
#include <ctime>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <chrono>
#include <iostream>

using namespace std;
using namespace chrono;

// ���������� ����� � �������� ���������� �������
const int MATRIX_SIZE = 1500;
/// ������� InitMatrix() ��������� ���������� � �������� 
/// ��������� ���������� ������� ���������� ����������
/// matrix - �������� ������� ����
void InitMatrix(double** matrix)
{
	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		matrix[i] = new double[MATRIX_SIZE + 1];
	}

	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		for (int j = 0; j <= MATRIX_SIZE; ++j)
		{
			matrix[i][j] = rand() % 2500 + 1;
		}
	}
}

/// ������� SerialGaussMethod() ������ ���� ������� ������ 
/// matrix - �������� ������� �������������� ���������, �������� � ����,
/// ��������� ������� ������� - �������� ������ ������ ���������
/// rows - ���������� ����� � �������� �������
/// result - ������ ������� ����
auto SerialGaussMethod(double** matrix, const int rows, double* result)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	int k;
	double koef;

	// ������ ��� ������ ������
	for (k = 0; k < rows; ++k)
	{
		//
		for (int i = k + 1; i < rows; ++i)
		{
			koef = -matrix[i][k] / matrix[k][k];

			for (int j = k; j <= rows; ++j)
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	// �������� ��� ������ ������
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for (k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];

		//
		for (int j = k + 1; j < rows; ++j)
		{
			result[k] -= matrix[k][j] * result[j];
		}

		result[k] /= matrix[k][k];
	}
	return (t2 - t1);
}

auto ParallellGaussMethod(double** matrix, const int rows, double* result)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	int k;
	double koef;

	// ������ ��� ������ ������
	for (k = 0; k < rows; ++k)
	{
		//
		cilk_for(int i = k + 1; i < rows; ++i)
		{
			koef = -matrix[i][k] / matrix[k][k];

			for (int j = k; j <= rows; ++j)
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	// �������� ��� ������ ������
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for (k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];

		//
		cilk_for(int j = k + 1; j < rows; ++j)
		{
			result[k] -= matrix[k][j] * result[j];
		}

		result[k] /= matrix[k][k];
	}
	return (t2 - t1);
}


int main()
{
	srand((unsigned)time(0));

	int i;
	// ���-�� ����� � �������, ���������� � �������� �������
	const int test_matrix_lines = 4;

	double** test_matrix = new double* [MATRIX_SIZE];

	// ���� �� �������
	/*for (i = 0; i < test_matrix_lines; ++i)
	{
		// (test_matrix_lines + 1)- ���������� �������� � �������� �������,
		// ��������� ������� ������� ������� ��� ������ ����� ���������, �������� � ����
		test_matrix[i] = new double[test_matrix_lines + 1];
	}*/

	// ������ ������� ����
	double* result = new double[MATRIX_SIZE];

	// ������������� �������� �������
	/*test_matrix[0][0] = 2; test_matrix[0][1] = 5;  test_matrix[0][2] = 4;  test_matrix[0][3] = 1;  test_matrix[0][4] = 20;
	test_matrix[1][0] = 1; test_matrix[1][1] = 3;  test_matrix[1][2] = 2;  test_matrix[1][3] = 1;  test_matrix[1][4] = 11;
	test_matrix[2][0] = 2; test_matrix[2][1] = 10; test_matrix[2][2] = 9;  test_matrix[2][3] = 7;  test_matrix[2][4] = 40;
	test_matrix[3][0] = 3; test_matrix[3][1] = 8;  test_matrix[3][2] = 9;  test_matrix[3][3] = 2;  test_matrix[3][4] = 37;
	*/
	InitMatrix(test_matrix);
	duration<double> duration = SerialGaussMethod(test_matrix, MATRIX_SIZE, result);
	cout << endl;
	cout << "Duration is: " << duration.count() << " seconds" << endl; // ������� ����
	duration = ParallellGaussMethod(test_matrix, MATRIX_SIZE, result);
	cout << "Parrell Duration is: " << duration.count() << " seconds" << endl; // ������� ����
	for (i = 0; i < MATRIX_SIZE; ++i)
	{
		delete[]test_matrix[i];
	}
	delete[] result;
	system("pause");
	return 0;
}