#include <iostream> 
#include <cmath>
#include <fstream>
#include<math.h> 

using namespace std;

int main()

{

    const int N = 20;

    int T = 0;

    double c = 1, dt = 0.01, dx = 0.1052631579;

    double U[42], Un[42], Uno[42];

    for (int i = 0; i < N; i++)
    {
        U[i] = 0.0;
        Uno[i] = 0.0;
    }
    do
    {
        U[0] = 1.0;
        Uno[0] = 1.0;
        Un[0] = 1.0;
        U[N] = 0.0;
        Uno[N] = 0.0;
        Un[N] = 0.0;
        for (int i = 1; i < N - 1; i++)

        {
            Un[i] = Uno[i] - c * dt / dx * (U[i + 1] - U[i - 1]);
        }

        for (int i = 1; i < N - 1; i++)

        {

            Uno[i] = U[i];

            U[i] = Un[i];
        }

        T++;

    } while (T != 50);

    fstream fout("Данные2.dat", ios::out);
    
    fout << "\"X\",\"An\"" << endl;
    for (int i = 0; i <= N; i++) {

        fout << Un[i] <<  endl;
        
    }
    

}