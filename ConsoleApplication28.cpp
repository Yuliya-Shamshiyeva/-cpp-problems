// ConsoleApplication28.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
int main() {
	const int n = 51;
	int itP = 0, it = 0;;
	double h = 1.0 / (n - 1), dt = h * h / 2.0, Re = 50.0, ro = 1.0, eps = pow(10, -6);
	double Uin[n], Us[n][n], U[n][n], Un[n][n], Vs[n][n], V[n][n], Vn[n][n], P[n][n], Pn[n][n];
	double max_p, max_u, max_v;

	//Начальные условия
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			P[i][j] = 0.0;
			U[i][j] = 0.0;
			V[i][j] = 0.0;
		}
		Uin[i] = 0;
	}

	//Граничные условия
	for (int i = 2 * n / 3; i < n; i++) {
		Uin[i] = 1.0;
		U[0][i] = 1.0;
		P[0][i] = 1.0;
	}

	do {
		//1. Us, Vs 
		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				Us[i][j] = U[i][j] + dt * ((-U[i][j] * (U[i + 1][j] - U[i - 1][j]) - V[i][j] * (U[i][j + 1] - U[i][j - 1])) / (2 * h) +
					(U[i + 1][j] - 2.0 * U[i][j] + U[i - 1][j] + U[i][j + 1] - 2.0 * U[i][j] + U[i][j - 1]) / (Re * h * h));
				Vs[i][j] = V[i][j] + dt * ((-U[i][j] * (V[i + 1][j] - V[i - 1][j]) - V[i][j] * (V[i][j + 1] - V[i][j - 1])) / (2 * h) +
					(V[i + 1][j] - 2.0 * V[i][j] + V[i - 1][j] + V[i][j + 1] - 2.0 * V[i][j] + V[i][j - 1]) / (Re * h * h));
			}
		}

		//Граничные условия 
		for (int i = 0; i < n; i++) {
			//Стены
			Us[0][i] = Uin[i];
			Us[n - 1][i] = 0.0;
			Us[i][0] = 0.0;
			Us[i][n - 1] = 0.0;

			Vs[0][i] = 0.0;
			Vs[n - 1][i] = 0.0;
			Vs[i][0] = 0.0;
			Vs[i][n - 1] = 0.0;

			//Выход
			if (i < n / 3 || i >= 2 * n / 3) {
				Us[n - 1][i] = Us[n - 2][i];
				Vs[n - 1][i] = Vs[n - 2][i];
			}
		}


		//2. Давление
		do {
			for (int i = 1; i < n - 1; i++) {
				for (int j = 1; j < n - 1; j++) {
					Pn[i][j] = (P[i + 1][j] + P[i - 1][j] + P[i][j + 1] + P[i][j - 1] -
						ro * h * (Us[i + 1][j] - Us[i - 1][j] + Vs[i][j + 1] - Vs[i][j - 1]) / (2 * dt)) / 4.0;
				}
			}

			//Граничные условия 
			for (int i = 0; i < n; i++) {
				//Стены
				Pn[0][i] = Pn[1][i];
				Pn[i][0] = Pn[i][1];
				Pn[i][n - 1] = Pn[i][n - 2];
				Pn[n - 1][i] = Pn[n - 2][i];

				//Вход
				if (i >= 2 * n / 3) {
					Pn[0][i] = 1.0;
				}

				//Выход
				if (i < n / 3 || i >= 2 * n / 3) {
					Pn[n - 1][i] = 0.0;
				}
			}

			max_p = -1.0;
			//Максимальная разница по P
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (max_p < abs(Pn[i][j] - P[i][j])) {
						max_p = abs(Pn[i][j] - P[i][j]);
					}
					P[i][j] = Pn[i][j];
				}
			}

			itP++;
		} while (max_p > eps);

		//3. Un, Vn
		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				Un[i][j] = Us[i][j] - dt * (Pn[i + 1][j] - Pn[i - 1][j]) / (ro * 2 * h);
				Vn[i][j] = Vs[i][j] - dt * (Pn[i][j + 1] - Pn[i][j - 1]) / (ro * 2 * h);
			}
		}

		//Граничные условия 
		for (int i = 0; i < n; i++) {
			//Стены
			Un[0][i] = Uin[i];
			Un[n - 1][i] = 0.0;
			Un[i][0] = 0.0;
			Un[i][n - 1] = 0.0;

			Vn[0][i] = 0.0;
			Vn[n - 1][i] = 0.0;
			Vn[i][0] = 0.0;
			Vn[i][n - 1] = 0.0;

			//Выход
			if (i < n / 3 || i >= 2 * n / 3) {
				Un[n - 1][i] = Un[n - 2][i];
				Vn[n - 1][i] = Vn[n - 2][i];
			}
		}

		max_u = -1.0; max_v = -1.0;
		//Максимальная разница по U,V
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (max_u < abs(Un[i][j] - U[i][j])) {
					max_u = abs(Un[i][j] - U[i][j]);
				}
				U[i][j] = Un[i][j];
				if (max_v < abs(Vn[i][j] - V[i][j])) {
					max_v = abs(Vn[i][j] - V[i][j]);
				}
				V[i][j] = Vn[i][j];
			}
		}

		//Присваивание значений 
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				U[i][j] = Un[i][j];
				V[i][j] = Vn[i][j];
			}
		}
		cout << max_u << " " << max_v << endl;
		it++;
	} while (max_u > eps || max_v > eps);

	fstream fout("навье.dat", ios::out);
	fout << "Variables=\"X\",\"Y\",\"U\",\"V\",\"P\"" << endl;
	fout << "Zone I=" << n << ", J=" << n << ", F=POINT" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << i * h << "\t" << j * h << "\t" << U[i][j] << "\t" << V[i][j] << "\t" << P[i][j] << endl;
		}
	}

	//тут наши вывод как обычно пишешь свой, а то я всегда по разному пишу 
	cout << "Number of iter: " << it << endl;
	cout << "Number of iter Pressure: " << itP << endl;
	cout << "Max difference:\n";
	cout << "by U: " << max_u << ",\t by V: " << max_v << endl;

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
