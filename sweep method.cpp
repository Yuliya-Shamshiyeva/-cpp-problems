// ConsoleApplication10.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

float f(float x, float y, float t) {
	return 0;
}

int main(int argc, char** argv) {
	int N = 10, M = 10, P = 10;//размерность
	float u0 = 0, psi1 = 1, psi2 = 1, psi3 = 0, psi4 = 0, Lx = 1, Ly = 1, T = 1, A, B, C;
	float u[25][25][25], F[25][25], alpha[25], beta[25], x[20], y[20], t[21], hx = Lx / N, hy = Ly / M, tau = T / (2 * P), a = 1;
	x[0] = 0;
	y[0] = 0;
	t[0] = 0;
	for (int i = 1; i <= N; i++) //х
		x[i] = x[i - 1] + hx;
	for (int i = 1; i <= M; i++) //y
		y[i] = y[i - 1] + hy;
	for (int i = 1; i <= 2 * P; i++) //t (время)
		t[i] = t[i - 1] + tau;

	//Начальные условия
	for (int i = 0; i <= N; i++)
		for (int j = 0; j <= M; j++)
			u[i][j][0] = u0;
	//Краевые условия
	for (int k = 0; k <= 2 * P; k++) {
		for (int i = 0; i <= M; i++)
			if (i >= M / 2) {
				u[0][i][k] = 1;
			}
			else {
				u[0][i][k] = 0;
			}
		for (int i = 0; i <= M; i++)
			if (i >= M / 2) {
				u[N][i][k] = 0;
			}
			else {
				u[N][i][k] = 1;
			}
		for (int i = 0; i <= N; i++) {
			if (i <= N / 4) {
				u[i][0][k] = 1;
			}
			if ((i >= N / 4) && (i <= 2 * N / 4)) {
				u[i][0][k] = 0;
			}
			if ((i >= 2* N / 4) && (i <= 3 * N / 4)) {
				u[i][0][k] = 1;
			}
			if  (i >= 3 * N / 4) {
				u[i][0][k] = 0;
			}
		}
		for (int i = 0; i <= N; i++) {
			if (i <= N / 4) {
				u[i][M][k] = 0;
			}
			if ((i >= N / 4) && (i <= 2 * N / 4)) {
				u[i][M][k] = 1;
			}
			if ((i >= 2 * N / 4) && (i <= 3 * N / 4)) {
				u[i][M][k] = 0;
			}
			if (i >= 3 * N / 4) {
				u[i][M][k] = 1;
			}
		}
	}
	for (int k = 0; k <= 2 * P; k += 2) {
		//Иницмализируем оэфициентф матрицы по х
		A = -1 / (hx*hx);
		B = 1 / tau + 2 / (hx*hx);
		C = -1 / (hx*hx);
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				F[i][k * 2] = u[i][j][k * 2] / tau; 
				//ПРОГОНКА 1
		for (int j = 1; j < M - 1; j++) {
			alpha[1] = 0;
			beta[1] = 0;
			if (j >= M / 2) {
				alpha[1] = 1;
				beta[1] = 1;
			}
			for (int i = 2; i < N; i++) {
				alpha[i] = -C / (A*alpha[i - 1] + B);
				beta[i] = (F[i - 1][k * 2] - A * beta[i - 1]) / (A*alpha[i - 1] + B);
			}
			u[N - 1][j][k + 1] = (F[N - 1][k] - A * beta[N - 1]) / (C + A * alpha[N - 1]);
			for (int i = N - 1; i > 0; i--) {
				u[i][j][k + 1] = alpha[i + 1] * u[i + 1][j][k + 1] + beta[i + 1]; //puoluchaem znacheniya na sloe k+1/2
			}
		}
		//Иницмализируем оэфициентф матрицы по у
		A = -1 / (hy*hy);
		B = 1 / tau + 2 / (hy*hy);
		C = -1 / (hy*hy);
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				F[j][k] = u[i][j][k + 1] / tau; 

				//ПРОГОНКА 2

		for (int i = 1; i < N - 1; i++) {
			alpha[1] = 0;
			beta[1] = 0;
			for (int i = 0; i <= N; i++) {
				if (i <= N / 4) {
					u[i][0][k] = 1;
				}
				if ((i >= N / 4) && (i <= 2 * N / 4)) {
					u[i][0][k] = 0;
				}
				if ((i >= 2 * N / 4) && (i <= 3 * N / 4)) {
					u[i][0][k] = 1;
				}
				if (i >= 3 * N / 4) {
					u[i][0][k] = 0;
				}
			}
			for (int i = 0; i <= N; i++) {
				
				if ((i >= N / 4) && (i <= 2 * N / 4)) {
					alpha[1] = 0;
					beta[1] = 1;
				}
				
				if (i >= 3 * N / 4) {
					alpha[1] = 0;
					beta[1] = 1;
				}
			}
			if ((i >= N / 4) && (i <= 2 * N / 4)) {
				alpha[1] = 0;
				beta[1] = 1;
			}
			for (int j = 2; j < M; j++) {
				alpha[j] = -C / (A*alpha[j - 1] + B);
				beta[j] = (F[j - 1][k] - A * beta[j - 1]) / (A*alpha[j - 1] + B);
			}
			u[i][M - 1][k + 2] = (F[M - 1][k] - A * beta[M - 1]) / (C + A * alpha[M - 1]);
			for (int j = M - 1; j > 0; j--) {
				u[i][j][k + 2] = alpha[j + 1] * u[i + 1][j][k + 2] + beta[j + 1]; 
			}
		}
	}
	//Ввыводим на экран результат
	for (int k = 0; k <= 2 * P; k += 2) {
		cout <<"t="<<t[k]<<endl;
		for (int j = M; j >= 0; j--) {
			for (int i = 0; i <= N; i++)
				cout <<left<< setw(10)<< u[i][j][k]<<"\t";
			cout <<endl;
		}
	}
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
