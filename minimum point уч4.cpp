// ConsoleApplication15.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cmath>
#include <iomanip>

#define eps 0.001

using namespace std;

double Fun(double x)
{
	return exp(x);
}
int F(int n)
{
	int f, f1 = 1, f2 = 1, m = 0;
	while (m < n - 1)
	{
		f = f1 + f2;
		f1 = f2;
		f2 = f;
		++m;
	}
	return f1;
}
void Fib(double a, double b)
{


	int delta = 0.001;

	double x1, x2, _x, xf1, xf2;
	int k = 0;
	int N = 0;
	double fn1 = 1, fn2 = 1, fn, f = (b - a) / eps;
	cout << "\tj\t\tFj" << endl;
	cout << "\t" << N << "\t\t" << fn1 << endl;
	while (fn1 < f)
	{
		fn = fn1 + fn2;
		fn1 = fn2;
		fn2 = fn;
		++N;
		cout << "\t" << N << "\t\t" << fn1 << endl;
	}

	x1 = a + (double)F(N - 2) / F(N) * (b - a);
	x2 = a + (double)F(N - 1) / F(N) * (b - a);
	xf1 = Fun(x1);
	xf2 = Fun(x2);
	cout << "j\ta \t\t\tb" << endl;
	cout << " " << "\t" << a << "\t\t\t" << b << endl;
P:
	if (xf1 >= xf2)
	{
		a = x1;
		x1 = x2;
		xf1 = xf2;
		x2 = a + (double)F(N - k - 1) / F(N - k) * (b - a);
		xf2 = Fun(x2);
	}
	else
	{
		b = x2;
		x2 = x1;
		xf2 = xf1;
		x1 = a + (double)F(N - k - 2) / F(N - k) * (b - a);
		xf1 = Fun(x1);
	}
	cout << N - k << "\t" << round(a * 100000) / 100000 << " \t\t" << round(b * 100000) / 100000 << endl;
	k++;

	if (fabs(b - a) <= eps)
	{
		_x = (a + b) / 2;
		cout << endl << "Interval sodergashiy tochku minimuma: [" << round(a * 100000) / 100000 << "; " << round(b * 100000) / 100000 << "]" << endl;
		cout << "Tochka minimuma : " << round(_x * 100000) / 100000 << endl;
		cout << "minJ: " << Fun(_x) << endl;
	}
	else
		goto P;
}

int main()
{
	double a;
	double b;
	cout << "Vvidite a i b: ";
	cin >> a >> b;
	
	Fib(a, b);
	system("pause");
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
