#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "mpi.h"

#define N 100

/*int main(int argc, char* argv[]) {
    srand(time(NULL));
    int ProcNum, ProcRank;
    int vector[N];  // ќбъ€вл€ем и определ€ем одномерный пассив
    int min = INT_MAX, ProcMin = INT_MAX;
    MPI_Status Status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    // Process 0
    if (ProcRank == 0) {
        for (int i = 0; i < N; i++) {
            vector[i] = rand();
            printf("%6d  ", vector[i]);
        }
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
    // Others process
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
}*/

void DataInitialization(int *x);

int main(int argc, char* argv[]) {
    int vector[N], ProcMin = INT_MAX, TotalMin;
    int ProcRank, ProcNum;
    MPI_Status Status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    if (ProcRank == 0)
        DataInitialization(vector);
    MPI_Bcast(&vector, N, MPI_INT, 0, MPI_COMM_WORLD);
    int part = N / ProcNum;
    int i1 = part * ProcRank;
    int i2 = (ProcRank == ProcNum - 1)? N : part * (ProcRank + 1);
    for (; i1 < i2; i1++)
        if (ProcMin > vector[i1])
            ProcMin = vector[i1];

    if (ProcRank == 0) {  // Process 0
        TotalMin = ProcMin;
        /*for (int i = 1; i < ProcNum; i++) {
            MPI_Recv(&ProcMin, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &Status);
            if (TotalMin > ProcMin)
                TotalMin = ProcMin;
        }*/
    }
    else  // Others process
        MPI_Send(&ProcMin, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    MPI_Reduce(&ProcMin, &TotalMin, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

    // Print result
    if (ProcRank == 0)
        printf("\nTotal minimum = %d", TotalMin);
    MPI_Finalize();
    return 0;
}

void DataInitialization(int *x) {
    srand(time(NULL));
    for (int i = 0; i < N; i++) {
        x[i] = rand();
        printf("%6d  ", x[i]);
    }
}