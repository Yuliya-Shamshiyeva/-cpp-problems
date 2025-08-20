// ConsoleApplication23.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <fstream>
#include <cmath>
#define n 31
#define m 31
using namespace std;
int main() {
	int i, j, it = 0;
	double Lx = 3, Ly = 3;
	double dx, dy, newP[n][m], oldP[n][m], f[n][m], eps = 0.00001, max;
	dx = Lx / (n - 1); dy = Ly / (m - 1);
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			oldP[i][j] = 0.0;
			newP[i][j] = 0.0;
			f[i][j] = 0.0;
		}
	}
	do {
		for (i = 0; i < n / 3; i++) {
			oldP[0][i] = 1.0;
			newP[0][i] = 1.0;
			oldP[i + 2 * n / 3][n - 1] = 1.0;
			newP[i + 2 * n / 3][n - 1] = 1.0;
		}
		for (i = 0; i < n - 1; i++) {
			oldP[i][0] = 0.0;
			newP[i][0] = 0.0;
			oldP[n - 1][i] = 0.0;
			newP[n - 1][i] = 0.0;
		}
		for (i = 0; i < n - 1; i++) {
			for (j = 0; j < m - 1; j++) {
				newP[i][j] = 0.25 * (oldP[i + 1][j] + newP[i - 1][j] + oldP[i][j + 1] + newP[i][j - 1] + f[i][j] * dx * dx);
			}
		}
		for (i = 0; i < n / 3; i++) {
			oldP[0][i] = 1.0;
			newP[0][i] = 1.0;
			oldP[i + 2 * n / 3][n - 1] = 1.0;
			newP[i + 2 * n / 3][n - 1] = 1.0;
		}
		for (i = 0; i < n - 1; i++) {
			oldP[i][0] = 0.0;
			newP[i][0] = 0.0;
			oldP[n - 1][i] = 0.0;
			newP[n - 1][i] = 0.0;
		}
		max = 0.0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				if (max < fabs(newP[i][j] - oldP[i][j])) { max = fabs(newP[i][j] - oldP[i][j]); }
			}
		}
		for (i = 0; i < n; i++) {
			for (j = 0; j < m; j++) {
				oldP[i][j] = newP[i][j];
			}
		}
		it++;
	} while (max > eps);
	setlocale(LC_ALL, "Russian");
	fstream fout("ФизпроцессыЛаб9.dat", ios::out);
	fout << "Variables=\"X\",\"Y\",\"P\"" << endl;
	fout << "Zone I=" << n << ", J=" << m << ", F=POINT" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			fout << i * dx << "\t" << j * dy << "\t" << newP[i][j] << endl;
		}
	}
	cout << "Максимальная разница: " << max << endl;
	cout << "Количество итераций:" << it << endl;
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
