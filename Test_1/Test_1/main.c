#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "mpi.h"

#define N 300

int main(int argc, char* argv[]) {
    srand(time(NULL));
    int ProcNum, ProcRank, RecvRank;
    int vector[N];  // Объявляем и определяем одномерный пассив
    for (int i = 0; i < N; i++)
        vector[i] = rand();
    int min = INT_MAX;
    int ProcMin = INT_MAX;
    MPI_Status Status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    // Выполняется процесс 0
    if (ProcRank == 0) {
        for (int j = 0; j < N; j++)
            printf("%6d  ", vector[j]);
        for (int i = 0; i < N / ProcNum; i++) {
            if (vector[i] < min)
                min = vector[i];
        }
        printf("\n\n Minimum from process %d = %d", ProcRank, min);
        for (int i = 1; i < ProcNum; i++) {
            MPI_Recv(&ProcMin, 1, MPI_INT, i,
                MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
            if (min > ProcMin)
                min = ProcMin;
            printf("\n Minimum from process %d = %d", Status.MPI_SOURCE, ProcMin);
        }
        printf("\n\n Total minimum = %d\n", min);
    }
    // Выполняются остальные процессы
    else {
        if (ProcRank == ProcNum - 1) {
            for (int i = ProcRank * (N / ProcNum); i < N; i++)
                if (vector[i] < ProcMin)
                    ProcMin = vector[i];
        }
        else {
            for (int i = ProcRank * (N / ProcNum); i < ProcRank * (N / ProcNum) + N / ProcNum; i++)
            if (vector[i] < ProcMin)
                ProcMin = vector[i];
        }
        MPI_Send(&ProcMin, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}