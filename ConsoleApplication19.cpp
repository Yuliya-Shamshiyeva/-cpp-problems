// ConsoleApplication19.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main() {
	int const N = 31, M = 31;
	int k = 0, kk = 0;
	double dx = 1.0 / (N - 1), dy = 1.0 / (M - 1), eps = pow(10, -5), Re = 10.0, dt = 0.000025, udiff, pdiff, vdiff;
	double** unew = new double* [N],
		** uold = new double* [N],
		** uhalf = new double* [N],
		** vnew = new double* [N],
		** vold = new double* [N],
		** vhalf = new double* [N],
		** Pold = new double* [N],
		** Pnew = new double* [N];
	for (int i = 0; i < N; i++) {
		unew[i] = new double[M];
		uold[i] = new double[M];
		uhalf[i] = new double[M];
		vnew[i] = new double[M];
		vold[i] = new double[M];
		vhalf[i] = new double[M];
		Pold[i] = new double[M];
		Pnew[i] = new double[M];

	}
	for (int i = 0; i <= N - 1; i++) {
		for (int j = 0; j <= M - 1; j++) {
			uold[i][j] = 0;
			vold[i][j] = 0;
			unew[i][j] = 0;
			vnew[i][j] = 0;
			Pold[i][j] = 0;
			Pnew[i][j] = 0;
			uhalf[i][j] = 0;
			vhalf[i][j] = 0;
		}
	}
	for (int i = 0; i <= (N - 1) / 3; i++) {
		uold[0][i] = 1.0;
		unew[0][i] = 1.0;
		uhalf[0][i] = 1.0;

	}
	for (int i = (N - 1) / 3; i < N; i++) {
		uold[0][i] = 1.0;
		unew[0][i] = 1.0;
		uhalf[0][i] = 1.0;


	}
	for (int i = 2 * (N - 1) / 3; i < 30; i++) {
		unew[i][0] = unew[i][1];
		uold[i][0] = uold[i][1];
		uhalf[i][0] = uhalf[i][1];
		unew[i][30] = unew[i][29];
		uold[i][30] = uold[i][29];
		uhalf[i][30] = uhalf[i][29];
	}
	do {
		for (int i = 0; i <= (N - 1) / 3; i++) {
			uold[0][i] = 1.0;
			unew[0][i] = 1.0;
			uhalf[0][i] = 1.0;

		}
		for (int i = (N - 1) / 3; i < N; i++) {
			uold[0][i] = 1.0;
			unew[0][i] = 1.0;
			uhalf[0][i] = 1.0;


		}
		for (int i = 2 * (N - 1) / 3; i < 30; i++) {
			unew[i][0] = unew[i][1];
			uold[i][0] = uold[i][1];
			uhalf[i][0] = uhalf[i][1];
			unew[i][30] = unew[i][29];
			uold[i][30] = uold[i][29];
			uhalf[i][30] = uhalf[i][29];
		}

		for (int i = 1; i < N - 1; i++) {
			for (int j = 1; j < M - 1; j++) {
				if (i >= 10 && i <= 20 && j >= 10 && j <= 20) {
					uhalf[i][j] = 0.;
					vhalf[i][j] = 0.;
				}
				else {
					uhalf[i][j] = uold[i][j] - dt * (uold[i][j] * (uold[i + 1][j] - uold[i - 1][j]) / (2.0 * dx) + vold[i][j] * (uold[i][j + 1] - uold[i][j - 1]) / (2.0 * dy) - ((uold[i + 1][j] - 2.0 * uold[i][j] + uold[i - 1][j]) / (dx * dx) + (uold[i][j + 1] - 2.0 * uold[i][j] + uold[i][j - 1]) / (dy * dy)) / Re);
					vhalf[i][j] = vold[i][j] - dt * (uold[i][j] * (vold[i + 1][j] - vold[i - 1][j]) / (2.0 * dx) + vold[i][j] * (vold[i][j + 1] - vold[i][j - 1]) / (2.0 * dy) - ((vold[i + 1][j] - 2.0 * vold[i][j] + vold[i - 1][j]) / (dx * dx) + (vold[i][j + 1] - 2.0 * vold[i][j] + vold[i][j - 1]) / (dy * dy)) / Re);
				}
			}
		}
		k = 0;
		do {
			for (int i = 0; i < N; i++) {
				Pold[i][0] = Pold[i][1];
				Pnew[i][0] = Pnew[i][1];
				Pold[i][N - 1] = Pold[i][N - 2];
				Pnew[i][N - 1] = Pnew[i][N - 2];
				Pold[N - 1][i] = Pold[N - 2][i];
				Pnew[N - 1][i] = Pnew[N - 2][i];
				Pold[0][i] = Pold[1][i];
				Pnew[0][i] = Pnew[1][i];
			}
			for (int i = 20; i < N; i++) {
				Pold[i][0] = 0.0;
				Pnew[i][0] = 0.0;
				Pold[i][30] = 0.0;
				Pnew[i][30] = 0.0;
			}


			for (int i = 1; i < N - 1; i++) {
				for (int j = 1; j < M - 1; j++) {
					if (i >= 10 && i <= 20 && j >= 10 && j <= 20) {
						Pold[i][j] = 0.;
						Pnew[i][j] = 0.;
					}
					else {
						for (int l = 10; l <= 20; l++) {
							Pold[l][10] = Pold[l][9];
							Pold[l][20] = Pold[l][21];
							Pold[10][l] = Pold[9][l];
							Pold[20][l] = Pold[21][l];
							Pnew[l][10] = Pnew[l][9];
							Pnew[l][20] = Pnew[l][21];
							Pnew[10][l] = Pnew[9][l];
							Pnew[20][l] = Pnew[21][l];
						}


						Pnew[i][j] = ((Pold[i + 1][j] + Pold[i - 1][j]) / (dx * dx) + (Pold[i][j + 1] + Pold[i][j - 1]) / (dy * dy) - ((uhalf[i + 1][j] - uhalf[i][j]) / dx + (vhalf[i][j + 1] - vhalf[i][j]) / dy) / dt) / (2 * (1.0 / (dx * dx) + (1.0 / (dy * dy))));
					}
				}
			}
			for (int i = 20; i < 31; i++) {
				Pold[i][0] = 0.0;
				Pnew[i][0] = 0.0;
				Pold[i][30] = 0.0;
				Pnew[i][30] = 0.0;
			}


			pdiff = 0.0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++) {
					if (fabs(Pnew[i][j] - Pold[i][j]) > pdiff)
						pdiff = fabs(Pnew[i][j] - Pold[i][j]);
				}
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++) {
					Pold[i][j] = Pnew[i][j];
				}


			k++;
		} while (pdiff > eps);

		for (int i = 1; i < N - 1; i++) {
			for (int j = 1; j < N - 1; j++) {
				if (i >= 10 && i <= 20 && j >= 10 && j <= 20) {
					unew[i][j] = 0.;
					vnew[i][j] = 0.;
				}
				else {
					unew[i][j] = uhalf[i][j] - dt * (Pnew[i][j] - Pnew[i - 1][j]) / dx;
					vnew[i][j] = vhalf[i][j] - dt * (Pnew[i][j] - Pnew[i][j - 1]) / dy;
				}
			}
		}
		udiff = 0.0;
		vdiff = 0.0;
		for (int i = 0; i <= N - 1; i++)
			for (int j = 0; j <= N - 1; j++) {
				if (fabs(unew[i][j] - uold[i][j]) > udiff)
					udiff = fabs(unew[i][j] - uold[i][j]);
				if (fabs(vnew[i][j] - vold[i][j]) > vdiff)
					vdiff = fabs(vnew[i][j] - vold[i][j]);
			}

		for (int i = 0; i <= N - 1; i++)
			for (int j = 0; j <= N - 1; j++) {
				uold[i][j] = unew[i][j];
				vold[i][j] = vnew[i][j];
			}
		kk++;
		cout << udiff << "\t" << vdiff << endl;
	} while (udiff > eps || vdiff > eps);

	cout << "iterations: " << kk << endl;
	ofstream fout("NavieS.dat", ios::out);
	fout << "Variables=\"X\",\"Y\",\"U\",\"V\",\"P\"" << endl;
	fout << "Zone I=" << N << ", J=" << M << ", F=POINT" << endl;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			fout << i * dx << "\t" << j * dy << "\t" << unew[i][j] << "\t" << vnew[i][j] << "\t" << Pnew[i][j] << endl;
			cout << i * dx << "\t" << j * dy << "\t" << unew[i][j] << "\t" << vnew[i][j] << "\t" << Pnew[i][j] << endl;
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
