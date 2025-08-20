#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

int main() {

	fstream fout("lab11.dat", ios::out);
	fout << "Variables=\"X\",\"Y\",\"U\",\"V\"" << endl;
	int n = 51, it = 0;
	double dx = 0.02, dy = 0.02, dt = dx * dy / 4.0, eps = pow(10, -6);
	double maxu, maxv, Re = 50.0;

	double** u_new = new double* [n], ** v_new = new double* [n],
		** u = new double* [n], ** v = new double* [n],
		** u_h = new double* [n], ** v_h = new double* [n],
		** A = new double* [n], ** B = new double* [n], ** C = new double* [n],
		** D_u = new double* [n], ** D_v = new double* [n],
		** alpha_u = new double* [n], ** alpha_v = new double* [n],
		** beta_u = new double* [n], ** beta_v = new double* [n];

	for (int i = 0; i < n; i++) {
		u_new[i] = new double[n]; v_new[i] = new double[n];
		u[i] = new double[n]; v[i] = new double[n];
		u_h[i] = new double[n]; v_h[i] = new double[n];
		A[i] = new double[n]; B[i] = new double[n]; C[i] = new double[n];
		D_u[i] = new double[n]; D_v[i] = new double[n];
		alpha_u[i] = new double[n]; alpha_v[i] = new double[n];
		beta_u[i] = new double[n]; beta_v[i] = new double[n];
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			u[i][j] = 0.0;
			v[i][j] = 0.0;
		}
		if (i >= 2 * n / 3) {
			u[i][0] = 1.0;
		}
	}

	do {

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
				u[i][j] = u_new[i][j];
				v[i][j] = v_new[i][j];
			}
		}

		it++;

		if (it % 10000 == 0) {
			fout << "Zone T=""\"" << it << "\""", I=" << n << ", J=" << n << ", F=POINT" << endl;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					fout << i * dx << "\t" << j * dy << "\t" << u[i][j] << "\t" << v[i][j] << endl;
				}
			}
		}
	} while (maxu > eps || maxv > eps);

	fout << "Zone T=""\"" << it << "\""", I=" << n << ", J=" << n << ", F=POINT" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << i * dx << "\t" << j * dy << "\t" << u[i][j] << "\t" << v[i][j] << endl;
		}
	}

	cout << "Kol-vo iter:" << it << endl;
	cout << "Max razniza:" << endl;
	cout << "U: " << maxu << ",  V: " << maxv << endl;

	system("pause");
	return 0;
}


