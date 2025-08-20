// ConsoleApplication22.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

const int n = 51;
int main() {
	double dx = 0.02, dy = 0.02, dt = 1.0 / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))),
		eps = pow(10, -6), a_sqr = 1.0, max, t=0;
	double Unew[n][n], U[n][n], Uhalf[n][n];
	double A, B, C, D[n][n];
	double alpha[n][n], beta[n][n];
	int n_it = 0;

	A = -a_sqr / (2.0 * dx * dx); B = a_sqr / (dx * dx) + 1.0 / dt; C = -a_sqr / (2.0 * dx * dx);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			U[i][j] = 0.0;
			if (j >= n / 3.0 && j < 2 * n / 3.0) {
				U[0][j] = 1.0;
				U[n - 1][j] = 1.0;
			}
		}
	}

	do {

			
		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				D[i][j] = U[i][j] / dt + a_sqr * (U[i + 1][j] - 2 * U[i][j] + U[i - 1][j]) / (2.0 * dx * dx) + a_sqr * (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]) / (dy * dy);
			}
		}

		for (int i = 0; i < n; i++) {
			if (i >= n / 3.0 && i < 2 * n / 3.0) {
				Uhalf[n - 1][i] = 1.0;
				Uhalf[0][i] = 1.0;
			}
			else {
				Uhalf[n - 1][i] = 0.0;
				Uhalf[0][i] = 0.0;
			}
		}

		for (int j = 1; j < n - 1; j++) {
			alpha[1][j] = 0.0;
			beta[1][j] = Uhalf[0][j];
			for (int i = 1; i < n - 1; i++) {
				alpha[i + 1][j] = -(A / (B + alpha[i][j] * C));
				beta[i + 1][j] = (D[i][j] - C * beta[i][j]) / (B + C * alpha[i][j]);
			}
			for (int i = n - 2; i >= 0; i--) {
				Uhalf[i][j] = alpha[i + 1][j] * Uhalf[i + 1][j] + beta[i + 1][j];
			}
		}

		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				D[i][j] = Uhalf[i][j] / dt - a_sqr * (U[i][j + 1] - 2.0 * U[i][j] + U[i][j - 1]) / (2.0 * dy * dy);
			}
		}

		for (int i = 0; i < n; i++) {
			Unew[i][0] = 0.0;
			Unew[i][n - 1] = 0.0;
			Unew[0][i] = Uhalf[0][i];
			Unew[n - 1][i] = Uhalf[n - 1][i];
		}

		for (int i = 1; i < n - 1; i++) {
			alpha[i][1] = 0.0;
			beta[i][1] = Unew[i][0];
			for (int j = 1; j < n - 1; j++) {
				alpha[i][j + 1] = -(A / (B + alpha[i][j] * C));
				beta[i][j + 1] = (D[i][j] - C * beta[i][j]) / (B + C * alpha[i][j]);
			}
			for (int j = n - 2; j >= 0; j--) {
				Unew[i][j] = alpha[i][j + 1] * Unew[i][j + 1] + beta[i][j + 1];
			}
		}

		max = -1;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (max < abs(Unew[i][j] - U[i][j])) {
					max = abs(Unew[i][j] - U[i][j]);
				}
				U[i][j] = Unew[i][j];
			}
		}
		t++;
		n_it++;
	}  
		while (max > eps);
		cout << "Maxsimalnaya razniza mezhdu Unew i U : " << max << endl;
		cout << "Kol-vo iterazii:" << n_it << endl;
	fstream fout("lab7.dat", ios::out);
	fout << "Variables=\"X\",\"Y\",\"T\",t\"" << endl;
	fout << "Zone I=" << n << ", J=" << n << ", F=POINT" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << i * dx << "\t" << j * dy << "\t" << t<< "\t" << Unew[i][j] << endl;
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
