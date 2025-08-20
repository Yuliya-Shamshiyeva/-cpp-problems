// ConsoleApplication21.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;
int main() {
	const int n = 51;
	double dx = 0.02, dy = 0.02, dt = 1.0 / (2 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
	double Uold[n][n], Unew[n][n];
	double max, e = pow(10, -6);
	int it = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Uold[i][j] = 0.0;
			0
		}
	}
	for (int i = 0; i < n / 3; i++) {
		Uold[n - 1][i] = 1.0;
		Uold[n - 1][i + 2 * n / 3] = 1.0;
	}
	do {
		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				Unew[i][j] = dt * ((Uold[i + 1][j] - 2 * Uold[i][j] + Uold[i -
					1][j]) / (dx * dx) + (Uold[i][j + 1] - 2 * Uold[i][j] + Uold[i][j - 1]) / (dy * dy)) +
					Uold[i][j];
			}
		}
		for (int i = 0; i < n; i++) {
			Unew[0][i] = 0;
			Unew[n - 1][i] = 0;
			Unew[i][0] = 0;
			Unew[i][n - 1] = 0;
		}
		for (int i = 0; i < n / 3; i++) {
			Unew[n - 1][i] = 1.0;
			Unew[n - 1][i + 2 * n / 3] = 1.0;
		}
		max = -1.0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (max < abs(Unew[i][j] - Uold[i][j])) {
					max = abs(Unew[i][j] - Uold[i][j]);
				}
				Uold[i][j] = Unew[i][j];
			}
		}
		it++;
	} while (max > e);
	cout << "Number of iterations: " << it << endl;
	cout << "Max difference: " << max << endl;
	fstream fout("Heat_eq_exp.dat", ios::out);
	fout << "Variables=\"X\",\"Y\",\"U\"" << endl;
	fout << "Zone I=" << n << ", J=" << n << ", F=POINT" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << i * dx << "\t" << j * dy << "\t" << Uold[i][j] << endl;
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
