
#include <iostream>
#include <fstream>

using namespace std;
const int n = 51;
int main() {
	double dx = 0.02, dy = 0.02;
	double max = -1, eps = pow(10, -6);
	double Pnew[n][n], P[n][n], f[n][n];
	double W = 1;
	int n_it = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			P[i][j] = 0.0;
			f[i][j] = 0.0;
			if (j >= n / 3.0 && j < 2 * n / 3.0) {
				P[0][j] = 1.0;
				P[n - 1][j] = 1.0;
			}
		}
	}
	double a;

	do {

		for (int i = 0; i < n; i++) {
			Pnew[i][0] = 0.0;
			Pnew[i][n - 1] = 0.0;
			if (i >= n / 3.0 && i < 2 * n / 3.0) {
				Pnew[n - 1][i] = 1.0;
				Pnew[0][i] = 1.0;
			}
			else {
				Pnew[n - 1][i] = 0.0;
				Pnew[0][i] = 0.0;
			}
		}

		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < n - 1; j++) {
				a = (-f[i][j] - (P[i + 1][j] + Pnew[i - 1][j]) / (dx * dx) - (P[i][j + 1] + Pnew[i][j - 1]) / (dy * dy)) / (-2.0 / (dx * dx) - 2.0 / (dy * dy));
				Pnew[i][j] = (a - (1 - 1 / W) * P[i][j]) * W;
			}
		}

		max = -1.0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (max < abs(Pnew[i][j] - P[i][j])) {
					max = abs(Pnew[i][j] - P[i][j]);
				}
				P[i][j] = Pnew[i][j];
			}
		}

		n_it++;

	} while (max > eps);

	fstream fout("lab10.dat", ios::out);
	fout << "Variables=\"X\",\"Y\",\"P\"" << endl;
	fout << "Zone I=" << n << ", J=" << n << ", F=POINT" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << i * dx << "\t" << j * dy << "\t" << Pnew[i][j] << endl;
		}
	}

	cout << "Max raznica: " << max << endl;
	cout << "Kol-vo iteracii:" << n_it << endl;

	system("pause");
	return 0;
}



