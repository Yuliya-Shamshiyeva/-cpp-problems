#include "pch.h"
#include <iostream>
#include <omp.h>
#include <conio.h>
const int nNumThreadss = 2;
using namespace std;
int main()
{
	long n = 1000000;
	long long sum = 0;
	long long sum1 = 0;
	long long sum2 = 0;
	omp_set_num_threads(nNumThreadss);
#pragma omp parallel reduction(+:sum)
	{
		if (omp_get_thread_num() == 0)
		{
			for (int i = 1; i < n / 2; i++)
			{
				sum1 += i;
			}
			cout << "[" << omp_get_thread_num() << "]: " << "sum: " << sum1 << endl << endl;
		}
		else
		{
			for (int i = n / 2; i <= n; i++)
			{
				sum2 += i;
			}
			cout << "[" << omp_get_thread_num() << "]: " << "sum: " << sum2 << endl << endl;
		}
	}
	sum = sum1 + sum2;
	cout << "summa: " << sum << endl;
}