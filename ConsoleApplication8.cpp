#include "pch.h"
#include "conio.h"
#include "math.h"
#include "iostream"
using namespace std;

int i, n, k, n1;
double z;
double A[50][50];
double B[50];
double eps[50];
double C[50];
double et[50];

int main()
{
	cout << "Vvedite razmernost matrici: ";
	cin >> n1;

	cout << "Vvedite elementy matrici:" << endl;
	for (i = 0; i < n1; i++)
		for (k = 0; k < n1; k++)
			cin >> A[i][k];

	cout << "Vasha matrix A:" << endl;
	for (i = 0; i < n1; i++)
	{
		for (k = 0; k < n1; k++)
			cout << A[i][k] << "\t ";
		cout << endl;
	}

	cout << "Vvedite elementy vectora B:" << endl;
	for (i = 0; i < n1; i++)
		cin >> B[i];

	cout << "Vash vector B:" << endl;
	for (i = 0; i < n1; i++)
		cout << B[i] << endl;

	n = n1 - 1;
	eps[0] = -A[0][1] / A[0][0];
	et[0] = B[0] / A[0][0];

	for (i = 1; i < n; i++)
	{
		z = A[i][i] + A[i][i - 1] * eps[i - 1];
		eps[i] = -A[i][i + 1] / z;
		et[i] = (B[i] - A[i][i - 1] * et[i - 1]) / z;
	}

	C[n] = (B[n] - A[n][n - 1] * et[n - 1]) / (A[n][n] + A[n][n - 1] * eps[n - 1]);

	for (i = n - 1; i >= 0; i--)
		C[i] = eps[i] * C[i + 1] + et[i];

	cout << "Rezultat C:" << endl;
	for (i = 0; i < n1; i++)
		cout << C[i] << endl;

	return 0;
}