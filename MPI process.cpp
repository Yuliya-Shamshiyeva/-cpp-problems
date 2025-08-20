#include <iostream>
#include <mpi.h>

using namespace std;

int main()
{
	MPI_Init(NULL, NULL); //инициализация MPI
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Status status;
	int N; //Объявляем числа N и K
	int K;

	if (world_rank == 0) //главный процесс
	{
		cout << "Total processors: " << world_size << endl;//Всего процессов
		cout << "N: ";
		cin >> N; //Ввод числа N
		for (int i = 1; i < world_size; i++) //Высылаем всем процессам значение N
		{
			MPI_Send(&N, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
		}
		K = world_rank * N;// Определяем K для соответствующего процесса
		MPI_Send(&K, 1, MPI_INT, world_size - 1, 3, MPI_COMM_WORLD);// Высылаем последнему процессу значение K
		MPI_Recv(&K, 1, MPI_INT, world_rank + 1, 3, MPI_COMM_WORLD, &status);//Принимаем значение K от следующего процесса
		cout << "[ rank: " << world_rank << " ] --> " << K;//Вывод полученного K

	}
	if ((world_rank != 0) && (world_rank != world_size - 1))//Рассматриваем процессы кроме главного и последнего
	{
		MPI_Recv(&N, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status); //Принимаем значение N от главного процесса
		K = world_rank * N;// Определяем K для соответствующего процесса
		MPI_Send(&K, 1, MPI_INT, world_rank - 1, 3, MPI_COMM_WORLD);//Передаем значение K предыдущему процессу
		MPI_Recv(&K, 1, MPI_INT, world_rank + 1, 3, MPI_COMM_WORLD, &status);//Принимаем значение K от следующего процесса
		cout << "[ rank: " << world_rank << " ] --> " << K;//Вывод полученного K
	}
	if (world_rank == world_size - 1)//Рассматриваем последний процесс
	{
		MPI_Recv(&N, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);//Принимаем значение N от главного процесса
		K = world_rank * N;// Определяем K для соответствующего процесса
		MPI_Send(&K, 1, MPI_INT, world_rank - 1, 3, MPI_COMM_WORLD);//Передаем значение K предыдущему процессу
		MPI_Recv(&K, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);//Принимаем значение K от главного процесса
		cout << "[ rank: " << world_rank << " ] --> " << K;//Вывод полученного K
	}
	cout << endl;
	MPI_Finalize();
}
