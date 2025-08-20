// ConsoleApplication25.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

double*** Inverter_pivotChecker(double** a, double** e, int n, int shift) {
	double*** p = new double** [2];
	p[0] = a; p[1] = e;
	double* temp = new double[n];
	if (a[shift][shift] != 0) {
		return p;
	}
	else {
		for (int k = shift + 1; k < n; k++)
			if (a[k][shift] != 0) {
				for (int j = 0; j < n; j++) {
					temp[j] = a[shift][j];
					a[shift][j] = a[k][j];
					a[k][j] = temp[j];
				}
				for (int j = 0; j < n; j++) {
					temp[j] = e[shift][j];
					e[shift][j] = e[k][j];
					e[k][j] = temp[j];
				}
				p[0] = a; p[1] = e;
				delete temp;
				return p;
			}
	}
	return p;
}

double** Inverter(double** a, int n) {
	double** inv = new double* [n];
	double*** p;
	double t, t2;
	for (int i = 0; i < n; i++)
		inv[i] = new double[n];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i == j)
				inv[i][j] = 1.0;
			else
				inv[i][j] = 0.0;

	for (int i = 0; i < n; i++) {
		p = Inverter_pivotChecker(a, inv, n, i);
		a = p[0]; inv = p[1];
		t = a[i][i];
		for (int j = 0; j < n; j++) {
			a[i][j] /= t;
			inv[i][j] /= t;
		}
		for (int k = i + 1; k < n; k++) {
			t = -a[k][i] / a[i][i];
			for (int j = 0; j < n; j++) {
				a[k][j] += t * a[i][j];
				inv[k][j] += t * inv[i][j];
			}
		}
	}
	for (int k = n - 1; k >= 0; k--) {
		for (int i = 0; i < k; i++) {
			t = -a[i][k];
			for (int j = 0; j < n; j++) {
				a[i][j] += a[k][j] * t;
				inv[i][j] += inv[k][j] * t;
			}
		}

	}
	return inv;
}

double** matrix_multiplicator(double** u, double** v, int n) {
	double** product = new double* [n];
	for (int i = 0; i < n; i++)
		product[i] = new double[n];
	double sum = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			sum = 0;
			for (int k = 0; k < n; k++) {
				sum += u[i][k] * v[k][j];
			}
			product[i][j] = sum;
		}
	}
	return product;
}

double** Matrix_Adder(double** u, double** v, int n) {
	double** sum = new double* [n];
	for (int i = 0; i < n; i++) {
		sum[i] = new double[n];
	}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			sum[i][j] = u[i][j] + v[i][j];
	return sum;
}

int main() {
	int n = 51;
	double dx = 1.0 / (n - 1), dy = 1.0 / (n - 1);
	double** A = new double* [n], ** B = new double* [n], ** C = new double* [n], * D = new double[n],
		*** alpha = new double** [n], ** beta = new double* [n], ** P = new double* [n];
	double** inverse_matrix = new double* [n];

	for (int i = 0; i < n; i++) {
		A[i] = new double[n];
		B[i] = new double[n];
		C[i] = new double[n];
		P[i] = new double[n];
		beta[i] = new double[n];
		inverse_matrix[i] = new double[n];
		alpha[i] = new double* [n];
		for (int j = 0; j < n; j++) {
			alpha[i][j] = new double[n];
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			P[i][j] = 0.0;
			A[i][j] = 0.0;
			B[i][j] = 0.0;
			C[i][j] = 0.0;
			beta[i][j] = 0.0;
		}
		D[i] = 0.0;
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				A[i][j] = 1.0 / (dx * dx);
				C[i][j] = A[i][j];
				B[i][j] = -2.0 / (dx * dx) - 2.0 / (dx * dx);
				if (i != n - 1)
					B[i + 1][j] = 1.0 / (dy * dy);
				if (j != n - 1)
					B[i][j + 1] = 1.0 / (dy * dy);
			}
		}
	}


	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				alpha[k][i][j] = 0.0;
			}
		}
	}

	for (int i = n / 3; i < 2 * n / 3; i++) {
		beta[1][i] = 1.0;
		P[n - 1][i] = 1.0;
	}

	double sum;
	double* beta_vector = new double[n];

	for (int h = 1; h < n - 1; h++) {

		inverse_matrix = Inverter(Matrix_Adder(B, matrix_multiplicator(C, alpha[h], n), n), n);

		//alpha
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				sum = 0;
				for (int k = 0; k < n; k++)
				{
					sum += -inverse_matrix[i][k] * A[k][j];
				}
				alpha[h + 1][i][j] = sum;
			}
		}

		//beta
		for (int i = 0; i < n; i++) {
			sum = 0;
			for (int j = 0; j < n; j++) {
				sum += C[i][j] * beta[h][j];
			}
			beta_vector[i] = D[i] - sum;
		}

		for (int i = 0; i < n; i++) {
			sum = 0;
			for (int j = 0; j < n; j++) {
				sum += inverse_matrix[i][j] * beta_vector[j];
			}
			beta[h + 1][i] = sum;
		}
	}

	//прогонка
	for (int h = n - 2; h >= 0; h--) {
		for (int i = 0; i < n; i++) {
			sum = 0;
			for (int j = 0; j < n; j++) {
				sum += alpha[h + 1][i][j] * P[h + 1][j];
			}
			P[h][i] = sum + beta[h + 1][i];
		}
	}

	ofstream fout("Answer.dat", ios::out);
	fout << "VARIABLES=\"X\",\"Y\",\"P\"" << endl;
	fout << "ZONE I=" << n << "J=" << n << "F=point" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << i * dx << "\t" << j * dy << "\t" << P[i][j] << "\t" << endl;
		}
	}

	system("pause");
	return 0;
}


// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
