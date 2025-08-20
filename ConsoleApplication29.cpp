
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

int main() {
	//Бюргерс
	int n= 51, it = 0, itP = 0 ;
	double dx = 0.02, dy = 0.02, dt = dx * dy / 2.0, ro = 1.0, eps = pow(10, -6);
	double maxu, maxv, maxp, Re = 50.0;//для стабильности
	

	double** u_new = new double* [n], ** v_new = new double* [n], ** u = new double* [n], ** v = new double* [n],** u_h = new double* [n], ** v_h = new double* [n],
		** A = new double* [n], ** B = new double* [n], ** C = new double* [n],** D_u = new double* [n], ** D_v = new double* [n],
		** alpha_u = new double* [n], ** alpha_v = new double* [n],** beta_u = new double* [n], ** beta_v = new double* [n],** P = new double* [n], ** Pn = new double* [n],
		** u_1 = new double* [n], ** v_1 = new double* [n], **Un = new double* [n],  ** Vn = new double* [n];
	for (int i = 0; i < n; i++) {
		u_new[i] = new double[n]; v_new[i] = new double[n];//временные
		u[i] = new double[n]; v[i] = new double[n];//конечные
		u_h[i] = new double[n]; v_h[i] = new double[n];//половинные
		A[i] = new double[n]; B[i] = new double[n]; C[i] = new double[n];
		D_u[i] = new double[n]; D_v[i] = new double[n];
		alpha_u[i] = new double[n]; alpha_v[i] = new double[n];
		beta_u[i] = new double[n]; beta_v[i] = new double[n];
		P[i] = new double[n];
		Pn[i] = new double[n];
		
		//промежуточные
		v_1[i] = new double[n];
		u_1[i] = new double[n];
		Vn[i] = new double[n];
		Un[i] = new double[n];
	}



	//Начальные условия ВЫХОДЫ ПО НЕЙМАНУ
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			//скорости
			u[i][j] = 0.0;//по х
			v[i][j] = 0.0;//по у
			P[i][j] = 0.0;//давление 
		}
		if (i >= 2 * n / 3) {
			u[i][0] = 1.0;//на входе
		}
	}

	do {
		//нач условия
		maxu = -10; maxv = -10;
		for (int i = 0; i < n; i++) {
			u_new[i][n - 1] = 0.0;
			v_new[i][n - 1] = 0.0;

			u_new[i][0] = 0.0;
			v_new[i][0] = 0.0;

			if (i < 2 * n / 3) {
				u_new[0][i] = 0.0;
			}
			else {
				u_new[0][i] = 1.0;
			}
			v_new[0][i] = 0.0;
		}

		//X
		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				A[i][j] = -1.0 / (2.0 * Re * dx * dx) + u[i][j] / (2.0 * dx);
				B[i][j] = 1.0 / dt + 1.0 / (Re * dx * dx) - u[i][j] / (2.0 * dx);
				C[i][j] = -1.0 / (2.0 * Re * dx * dx);

				D_u[i][j] = u[i][j] / dt + 0.5 * ((u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) / (Re * dx * dx) -
					u[i][j] * (u[i + 1][j] - u[i][j]) / dx) +
					(u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / (Re * dy * dy) -
					v[i][j] * (u[i][j + 1] - u[i][j]) / dy;

				D_v[i][j] = v[i][j] / dt + 0.5 * ((v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]) / (Re * dx * dx) -
					u[i][j] * (v[i + 1][j] - v[i][j]) / dx) +
					(v[i][j + 1] - 2 * v[i][j] + v[i][j - 1]) / (Re * dy * dy) -
					v[i][j] * (v[i][j + 1] - v[i][j]) / dy;
			}
		}

		for (int j = 1; j < n - 1; j++) {
			alpha_u[1][j] = 0.0;
			alpha_v[1][j] = 0.0;
			beta_u[1][j] = u_new[0][j];
			beta_v[1][j] = v_new[0][j];

			for (int i = 1; i < n - 1; i++) {
				alpha_u[i + 1][j] = -A[i][j] / (B[i][j] + C[i][j] * alpha_u[i][j]);
				beta_u[i + 1][j] = (D_u[i][j] - C[i][j] * beta_u[i][j]) / (B[i][j] + C[i][j] * alpha_u[i][j]);
				alpha_v[i + 1][j] = -A[i][j] / (B[i][j] + C[i][j] * alpha_u[i][j]);
				beta_v[i + 1][j] = (D_v[i][j] - C[i][j] * beta_v[i][j]) / (B[i][j] + C[i][j] * alpha_u[i][j]);
			}

			if (n / 3 <= j && j < 2 * n / 3) {
				u_h[n - 1][j] = 0.0;
				v_h[n - 1][j] = 0.0;
			}
			else {
				u_h[n - 1][j] = beta_u[n - 1][j] / (1.0 - alpha_u[n - 1][j]);
				v_h[n - 1][j] = beta_v[n - 1][j] / (1.0 - alpha_v[n - 1][j]);
			}

			for (int i = n - 2; i >= 0; i--) {
				u_h[i][j] = alpha_u[i + 1][j] * u_h[i + 1][j] + beta_u[i + 1][j];
				v_h[i][j] = alpha_v[i + 1][j] * v_h[i + 1][j] + beta_v[i + 1][j];
			}
		}

		for (int i = 1; i < n - 1; i++) {
			u_new[n - 1][i] = u_h[n - 1][i];
			v_new[n - 1][i] = v_h[n - 1][i];
		}

		//Y
		for (int i = 1; i < n - 1; i++) {
			//НАХОДИМ КОЭФ А Б С Д
			for (int j = 1; j < n - 1; j++) {
				A[i][j] = -1.0 / (2.0 * Re * dy * dy) + v[i][j] / (2.0 * dy);
				B[i][j] = 1.0 / dt + 1.0 / (Re * dy * dy) - v[i][j] / (2.0 * dy);
				C[i][j] = -1.0 / (2.0 * Re * dy * dy);

				D_u[i][j] = u_h[i][j] / dt - (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / (Re * 2.0 * dy * dy) +
					v[i][j] * (u[i][j + 1] - u[i][j]) / (2.0 * dy);

				D_v[i][j] = v_h[i][j] / dt - (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]) / (Re * 2.0 * dy * dy) +
					v[i][j] * (v[i][j + 1] - v[i][j]) / (2.0 * dy);
			}
		}

		for (int i = 1; i < n - 1; i++) {
			alpha_u[i][1] = 0.0;
			alpha_v[i][1] = 0.0;
			beta_u[i][1] = u_new[i][0];
			beta_v[i][1] = v_new[i][0];

			for (int j = 1; j < n - 1; j++) {
				alpha_u[i][j + 1] = -(A[i][j] / (B[i][j] + alpha_u[i][j] * C[i][j]));
				beta_u[i][j + 1] = (D_u[i][j] - C[i][j] * beta_u[i][j]) / (B[i][j] + C[i][j] * alpha_u[i][j]);
				alpha_v[i][j + 1] = -(A[i][j] / (B[i][j] + alpha_v[i][j] * C[i][j]));
				beta_v[i][j + 1] = (D_v[i][j] - C[i][j] * beta_v[i][j]) / (B[i][j] + C[i][j] * alpha_v[i][j]);
			};
			//обратная прогонка
			for (int j = n - 2; j >= 0; j--) {
				u_new[i][j] = alpha_u[i][j + 1] * u_new[i][j + 1] + beta_u[i][j + 1];
				v_new[i][j] = alpha_v[i][j + 1] * v_new[i][j + 1] + beta_v[i][j + 1];
			}
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (maxu < abs(u_new[i][j] - u[i][j])) {
					maxu = abs(u_new[i][j] - u[i][j]);
				}
				if (maxv < abs(v_new[i][j] - v[i][j])) {
					maxv = abs(v_new[i][j] - v[i][j]);
				}
				u_1[i][j] = u_new[i][j];
				v_1[i][j] = v_new[i][j];
			}
		}


		//2. Давление
		do {
			for (int i = 1; i < n - 1; i++) {
				for (int j = 1; j < n - 1; j++) {
					Pn[i][j] = (P[i + 1][j] + P[i - 1][j] + P[i][j + 1] + P[i][j - 1] -
						ro * dx * (u_1[i + 1][j] - u_1[i - 1][j] + v_1[i][j + 1] - v_1[i][j - 1]) / (2 * dt)) / 4.0; //из аппроксимации
				}
			}

			//Граничные условия 
			for (int i = 0; i < n; i++) {
				//Стены
				Pn[0][i] = Pn[1][i];
				Pn[i][0] = Pn[i][1];
				Pn[i][n - 1] = Pn[i][n - 2];
				Pn[n - 1][i] = Pn[n - 2][i];

				

				//Выход
				if (i < n / 3 || i >= 2 * n / 3) {
					Pn[n - 1][i] = 0.0;//Дирихле
				}
			}

			maxp = -1.0;
			//Максимальная разница по P
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (maxp < abs(Pn[i][j] - P[i][j])) {
						maxp = abs(Pn[i][j] - P[i][j]);
					}
					P[i][j] = Pn[i][j];
				}
			}

			itP++;
		} while (maxp > eps);

		//3. Un, Vn
		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				Un[i][j] = u_1[i][j] - dt * (Pn[i + 1][j] - Pn[i - 1][j]) / (ro * 2 * dx);
				Vn[i][j] = v_1[i][j] - dt * (Pn[i][j + 1] - Pn[i][j - 1]) / (ro * 2 * dy);
			}
		}

		//Граничные условия 
		for (int i = 0; i < n; i++) {
			//Стены
			Un[0][i] = 0.0;
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

		maxu = -1.0; maxv = -1.0;
		//Максимальная разница по U,V
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (maxu < abs(Un[i][j] - u[i][j])) {
					maxu = abs(Un[i][j] - u[i][j]);
				}
				u[i][j] = Un[i][j];
				if (maxv < abs(Vn[i][j] - v[i][j])) {
					maxv = abs(Vn[i][j] - v[i][j]);
				}
				v[i][j] = Vn[i][j];
			}
		}

		//Присваивание значений 
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				u[i][j] = Un[i][j];
				v[i][j] = Vn[i][j];
			}
		}
		cout << maxu << " " << maxv << endl;
		it++;
	} while (maxu > eps || maxv > eps);



	fstream fout("навье.dat", ios::out);
	fout << "Variables=\"X\",\"Y\",\"U\",\"V\",\"P\"" << endl;
	fout << "Zone I=" << n << ", J=" << n << ", F=POINT" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << i * dx << "\t" << j * dy<< "\t" << u[i][j] << "\t" << v[i][j] << "\t" << P[i][j] << endl;
		}
	}
	cout << "Количество итераций:" << it << endl;
	cout << "Количество итераций Pressure: " << itP << endl;
	cout << "Максимальная разница:" << endl;
	cout << "по U: " << maxu << ", по V: " << maxv << endl;


	system("pause");
	return 0;
}





