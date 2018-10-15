#include <stdio.h>
#include "mpi.h"

#define N 10

int main(int argc, char* argv[]) {
    int ProcNum, ProcRank, RecvRank;
    int vector[N] = { 60, 30, 70, 50, 10, 15, 14, 13, 12, 11 };
    int min = 200;
    int ProcMin = 200;
    MPI_Status Status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    if (ProcRank == 0) {
        for (int i = 0; i < N / ProcNum; i++) {
            if (vector[i] < min)
                min = vector[i];
        }
        printf("\n Minimum from process %d = %d", ProcRank, min);
        for (int i = 1; i < ProcNum; i++) {
            MPI_Recv(&ProcMin, 1, MPI_INT, i,
                MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
            if (min > ProcMin)
                min = ProcMin;
            printf("\n Minimum from process %d = %d", Status.MPI_SOURCE, ProcMin);
        }
        printf("\n\n Total minimum = %d\n", min);
    }
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