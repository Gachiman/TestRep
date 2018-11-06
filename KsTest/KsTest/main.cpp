#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <time.h>

void Initialization(char* str, int size);  // String initialization
int FindCount(char* str, int size);  // Search count of words

int main(int argc, char *argv[]) {
    int count;  // Count of symbols
    char *str = nullptr;
    int ProcNum, ProcRank;
    int num, sum = 0;
    double t1, t2;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0)
        count = (argc != 1) ? atoi(argv[1]) : 100;
    MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int size = count / ProcNum;

    if (ProcRank == 0) {
        str = new char[count];
        Initialization(str, count);
        t1 = MPI_Wtime();
        for (int i = 1; i < ProcNum - 1; i++)
            MPI_Send(str + size * i, size, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        if (ProcNum > 1)
            MPI_Send(str + size * (ProcNum - 1), count - (ProcNum - 1) * size,
                MPI_CHAR, ProcNum - 1, 0, MPI_COMM_WORLD);
        num = FindCount(str, size);
        delete[] str;
    }
    else {
        if (ProcRank == ProcNum - 1)
            size = count - (ProcNum - 1) * size;
        str = new char[size];
        MPI_Recv(str, size, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        num = FindCount(str, size);
        delete[] str;
    }

    MPI_Reduce(&num, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Print result
    if (ProcRank == 0) {
        t2 = MPI_Wtime();
        std::cout << "Sum of word = " << sum << std::endl
            << "Time spent = " << t2 - t1 << std::endl;
    }

    MPI_Finalize();
    return 0;
}

void Initialization(char* str, int size) {
    char alphabet[27];
    srand(time(0));
    for (int i = 0, j = 'A'; i < 26; i++, j++)
        alphabet[i] = j;
    alphabet[26] = ' ';

    for (int i = 0; i < size; i++)
        str[i] = alphabet[rand() % 27];

    for (int i = 0; i < size; i++)
        std::cout << str[i];
    std::cout << std::endl;
}

int FindCount(char* str, int size) {
    int space = 0;
    for (int i = 0; i < size; i++) {
        if (str[i] == ' ')
            space++;
    }
    return space;
}
