#include <iostream>
#include <iomanip>

using namespace std;

const double N = 100;
const double eps = 0.0000001;

double f(double x)
{
	return x * x * x * x + 2 * x * x + 1;
}

int main()
{
	double a, b, h;
	cout << "enter a: ";
	cin >> a;
	cout << "enter b: ";
	cin >> b;

	h = (b - a) / N;

	double yMin = f(a);
	double xMin = a;

	//double yMax = f(b);
	//double xMax = b;

	cout << "i" << "\t\t" << " u_i" << "\t\t" << " J_i" << endl;

	double x; int i;

	for (i = 1, x = a + (i - 1) * h; i < N + 1, x <= b; x += h, i++)
	{
		double y = f(x);
		if (y < yMin)
		{
			xMin = x;
			yMin = y;
		}
		cout << i << "\t\t" << x << "\t\t" << y << endl;
		cout << endl;
	}

	cout << setprecision(2) << xMin << " - minimum x" << endl;
	cout << yMin << " - minimum y" << endl;

	cout << endl;

	//cout « xMax « " is maximum x " « endl;
	//cout « yMax « " is maximum y " « endl;
}
