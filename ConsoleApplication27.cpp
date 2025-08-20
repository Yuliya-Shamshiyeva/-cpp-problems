
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
int main() {
	fstream fout("lab13.dat", ios::out);  ЗАПИСЬ ДЛЯ ПО ВРЕМЕНИ
	fout << "Variables=\"X\",\"Y\",\"Z\",\"U\"" << endl;

	int n = 51, it = 0;   
	double dx = 0.02, dy = 0.02, dz = 0.02, dt = 0.001,   ШАГИ 
		A, B, C, max, a_sqr = 1.0, eps = 0.000001,   ЕСЛИ 3 МЕРНЫЙ ВЕКТОР Д ТРЕХМЕРНЫЙ МАССИВ
		*** D = new double** [n], *** alpha = new double** [n], *** beta = new double** [n],

		*** Unew = new double** [n], *** U = new double** [n], *** Uone = new double** [n], *** Utwo = new double** [n],
		НОВЫЙ СЛОЙ     ТЕКУЩИЙ СЛОЙ       1\3 ПО СЛОЮ       2.3 ПО СЛОЮ
		** Ubup = new double* [n],БАУНДАРИ АП ЭТО ГДЕ НЕ НУЛИВЫЕ ЗНАЧЕНИЯ НА ГРАНИЦАХ СВЕРХУ - СЛЕВА ** Ubleft = new double* [n]; ГРАНИЦЫ 
	

		ИНИЦИАЛИЗИРОВАНИЕ ТРЕХМЕРНЫХ МАССИВОВ
	for (int i = 0; i < n; i++) {
		D[i] = new double* [n];
		alpha[i] = new double* [n];
		beta[i] = new double* [n];
		Unew[i] = new double* [n];
		U[i] = new double* [n];
		Uone[i] = new double* [n];
		Utwo[i] = new double* [n]; 
		Ubup[i] = new double[n];
		Ubleft[i] = new double[n];
		for (int j = 0; j < n; j++) {
			D[i][j] = new double[n];
			alpha[i][j] = new double[n];
			beta[i][j] = new double[n];
			Unew[i][j] = new double[n];
			U[i][j] = new double[n];
			Uone[i][j] = new double[n];
			Utwo[i][j] = new double[n];
		}
	}


	НАЧАЛЬНЫЕ УСЛОВИЯ
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				U[i][j][k] = 0.0;
			}
		}
	}


	ГРАНИЧНЫЕ УСЛОВИЯ
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i >= 0 && i < n / 3 && j >= 2 * n / 3 && j < n) {
				Ubup[i][j] = 1.0;    ВЕРХ ГДЕ КВАДРАТ СВЕРХУ
			}
			else {
				Ubup[i][j] = 0.0;   ГДЕ НЕТ КВАДРАТА СВЕРХУ
			}
			if (i >= 0 && i < n / 3 && j >= 0 && j < n / 3) {  СЛЕВА
				Ubleft[i][j] = 1.0;
			}
			else {
				Ubleft[i][j] = 0.0;   ПУСТО СЛЕВА
			}
			U[i][j][n - 1] = Ubup[i][j];  ВЕРХНЯЯ ПОСЛЕДНЯЯ СТЕНКА
			U[i][n - 1][j] = Ubleft[i][j];  ЛЕВАЯ
		}
	}

	do {  ВСЕ ДИРЕХЛЕ ПОЭТОМУ СРАЗУ ЗАДАЕМ ГРАНИЧНЫЕ УСЛОВИЯ
		max = -1.0;
		for (int i = 0; i < n; i++) {
			 СТЕНКИ ЗАНОВО
			for (int j = 0; j < n; j++) {
				Unew[0][i][j] = 0.0;
				Unew[n - 1][i][j] = 0.0; 
				Unew[i][0][j] = 0.0;
				Unew[i][n - 1][j] = Ubleft[i][j];
				Unew[i][j][0] = 0.0;
				Unew[i][j][n - 1] = Ubup[i][j];
			}
		}

		A = -a_sqr / (2.0 * dx * dx); C = A; B = 1.0 / dt + a_sqr / (dx * dx);

		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				for (int k = 1; k < n - 1; k++) {
					D[i][j][k] = U[i][j][k] / dt + a_sqr * (0.5 * (U[i + 1][j][k] - 2.0 * U[i][j][k] + U[i - 1][j][k]) / (dx * dx)
						+ (U[i][j + 1][k] - 2.0 * U[i][j][k] + U[i][j - 1][k]) / (dy * dy)
						+ (U[i][j][k + 1] - 2.0 * U[i][j][k] + U[i][j][k - 1]) / (dz * dz));
				}
			}
		}
		ОБРАТНАЯ ПРОГОНКА АЛЬФА=0 ТК НЕТ НЕЙМАНА РИСУНОЧКИ
		for (int j = 1; j < n - 1; j++) {
			for (int k = 1; k < n - 1; k++) {
				alpha[1][j][k] = 0.0;
				beta[1][j][k] = Unew[0][j][k];
				Uone[n - 1][j][k] = Unew[n - 1][j][k];

				for (int i = 1; i < n - 1; i++) {
					alpha[i + 1][j][k] = -(A / (B + alpha[i][j][k] * C));
					beta[i + 1][j][k] = (D[i][j][k] - C * beta[i][j][k]) / (B + C * alpha[i][j][k]);
				}
				ТУТ
				for (int i = n - 2; i >= 0; i--) {
					Uone[i][j][k] = alpha[i + 1][j][k] * Uone[i + 1][j][k] + beta[i + 1][j][k];
				}
			}
		}

		A = -a_sqr / (2.0 * dy * dy); C = A; B = 1.0 / dt + a_sqr / (dy * dy);

		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				for (int k = 1; k < n - 1; k++) {
					D[i][j][k] = Uone[i][j][k] / dt - a_sqr * 0.5 * (U[i][j + 1][k] - 2.0 * U[i][j][k] + U[i][j - 1][k]) / (dy * dy);
				}
			}
		}

		for (int i = 1; i < n - 1; i++) {
			for (int k = 1; k < n - 1; k++) {
				alpha[i][1][k] = 0.0;
				beta[i][1][k] = Unew[i][0][k];
				Utwo[i][n - 1][k] = Unew[i][n - 1][k];
				for (int j = 1; j < n - 1; j++) {
					alpha[i][j + 1][k] = -(A / (B + alpha[i][j][k] * C));
					beta[i][j + 1][k] = (D[i][j][k] - C * beta[i][j][k]) / (B + C * alpha[i][j][k]);
				}

				for (int j = n - 2; j >= 0; j--) {
					Utwo[i][j][k] = alpha[i][j + 1][k] * Utwo[i][j + 1][k] + beta[i][j + 1][k];
				}
			}
		}

		A = -a_sqr / (2.0 * dz * dz); C = A; B = 1.0 / dt + a_sqr / (dz * dz);

		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				for (int k = 1; k < n - 1; k++) {
					D[i][j][k] = Utwo[i][j][k] / dt - a_sqr * 0.5 * (U[i][j][k + 1] - 2.0 * U[i][j][k] + U[i][j][k - 1]) / (dz * dz);
				}
			}
		}

		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				alpha[i][j][1] = 0.0;
				beta[i][j][1] = Unew[i][j][0];
				for (int k = 1; k < n - 1; k++) {
					alpha[i][j][k + 1] = -(A / (B + alpha[i][j][k] * C));
					beta[i][j][k + 1] = (D[i][j][k] - C * beta[i][j][k]) / (B + C * alpha[i][j][k]);
				}
				for (int k = n - 2; k >= 0; k--) {
					Unew[i][j][k] = alpha[i][j][k + 1] * Unew[i][j][k + 1] + beta[i][j][k + 1];
				}
			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					if (max < abs(Unew[i][j][k] - U[i][j][k])) {
						max = abs(Unew[i][j][k] - U[i][j][k]);
					}
				}
			}
		}
		if (it % 50 == 0) {
			fout << "Zone T=""\"" << it << "\""", I=" << n << ", J=" << n << ", K=" << n << ", F=POINT" << endl;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					for (int k = 0; k < n; k++) {
						fout << i * dx << "\t" << j * dx << "\t" << k * dx << "\t" << U[i][j][k] << endl;
					}
				}
			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					U[i][j][k] = Unew[i][j][k];
				}
			}
		}
		it++;

	} while (max > eps);

	fout << "Zone T=""\"" << it << "\""", I=" << n << ", J=" << n << ", K=" << n << ", F=POINT" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				fout << i * dx << "\t" << j * dx << "\t" << k * dx << "\t" << U[i][j][k] << endl;
			}
		}
	}

	cout << "Max razniz " << max << endl;
	cout << "iteraz " << it << endl;

	system("pause");
	return 0;
}