// 31. Linear image filtering (vertical split). The Gaussian kernel 3x3.

#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <mpi.h>

void printImage(int* image, int width, int height);
void createKernel(double* kernel, int radius, int size, double sigma);
void processImage(int* originIm, int* checkImage, double* kernel, int height, int width, int size, int radius);
int Clamp(int value, int min, int max);
void equality_check(int * res, int * res2, int height, int width);

int main(int argc, char* argv[]) {
    srand(time(NULL));
    int ProcNum, ProcRank;
    int width, height, kernelRadius, size;
    double* kernel = nullptr;
    int* image = nullptr;
    int* checkImage = nullptr;  // For checking
    double sequentialTime, parallelTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    // Non-parallel algorithm

    if (ProcRank == 0) {
        // Input parametrs
        width = ((argc >= 2) && (atoi(argv[1]) > 0)) ? atoi(argv[1]) : 3;
        height = ((argc >= 3) && (atoi(argv[2]) > 0)) ? atoi(argv[2]) : 3;
        double sigma = (argc >= 4) ? atof(argv[3]) : 1.0;
        kernelRadius = ((argc >= 5) && (argv[4] > 0)) ? atoi(argv[4]) : 1;

        size = 2 * kernelRadius + 1;
        image = new int[width * height];
        kernel = new double[size * size];
        checkImage = new int[width * height];

        // Image initialization
        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++)
                image[i * width + j] = rand() % 100 + 1;
        std::cout << "Original image:\n";
        printImage(image, width, height);  // Print image

        createKernel(kernel, kernelRadius, size, sigma);  // Kernel initialization

        // Image for checking
        sequentialTime = MPI_Wtime();
        processImage(image, checkImage, kernel, height, width, size, kernelRadius);
        sequentialTime = MPI_Wtime() - sequentialTime;
        std::cout << "Image for checking:\n";
        printImage(checkImage, width, height); // Print image for checking
    }
    
    // Parallel algorithm
    
    MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kernelRadius, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (ProcRank != 0) {
        image = new int[width * height];
        kernel = new double[size * size];
    }
    MPI_Bcast(kernel, size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int* part = new int[ProcNum];
    if (ProcRank == 0) {
        for (int i = 0; i < ProcNum; i++)
            part[i] = width / ProcNum;
        int x = width % ProcNum;
        while (x > 0)
            part[x--]++;
    }

    MPI_Datatype mpi_LANE1;
    MPI_Type_vector(height, width / ProcNum, width, MPI_INT, &mpi_LANE1);
    MPI_Type_commit(&mpi_LANE1);

    MPI_Datatype mpi_LANE2;
    MPI_Type_vector(height, width / ProcNum + 1, width, MPI_INT, &mpi_LANE2);
    MPI_Type_commit(&mpi_LANE2);

    if (ProcRank == 0)
        parallelTime = MPI_Wtime();

    int remainder = width % ProcNum;
    int size1 = width / ProcNum;

    if (ProcRank == 0) {
        for (int i = 1; i < ProcNum - remainder; i++)
            MPI_Send(image + ((width / ProcNum) * i), 1, mpi_LANE1, i, 0, MPI_COMM_WORLD);
        for (int i = ProcNum - remainder; i < ProcNum; i++)
            MPI_Send(image + ((width / ProcNum) * i), 1, mpi_LANE2, i, 0, MPI_COMM_WORLD);
    }
    else {
        if (ProcRank  < ProcNum - remainder)
            MPI_Recv(image, 1, mpi_LANE1, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        else {
            MPI_Recv(image, 1, mpi_LANE2, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            size1 = width / ProcNum + 1;
        }
    }

    int* result = new int[height * size1];

    processImage(image, result, kernel, height, size1, size, kernelRadius);

    if (ProcRank != 0) printImage(result, size1, height);

    if (ProcRank == 0) {
        printImage(result, size1, height);
        int* res1 = new int[height * size1];
        int* res2 = new int[height * size1 + 1];
        for (int procs = 1; procs < ProcNum - remainder; procs++) {
            MPI_Recv(res1, height * size1, MPI_INT, procs, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < height; i++)
                for (int j = 0; j < size1; j++)
                    image[i * size1 + j + size1 * procs] = res1[i * size1 + j];
        }
        for (int procs = ProcNum - remainder; procs < ProcNum; procs++) {
            int size2 = ProcNum - remainder;
            MPI_Recv(res2, height * (size1 + 1), MPI_INT, procs, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < height; i++)
                for (int j = 0; j < (size1 + 1); j++)
                    image[i * width + j] = res2[i * (size1 + 1) + j];
        }
        for (int i = 0; i < height; i++)
            for (int j = 0; j < size1; j++)
                image[i * height + j] = result[i * size1 + j];
    }
    else
        MPI_Send(result, height * size1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    /*double tmp;
    for (int i = 0; i < height; i++)
        for (int j = 0; j < size1; j++) {
            tmp = 0;
            for (int y = -kernelRadius; y <= kernelRadius; y++)
                for (int x = -kernelRadius; x <= kernelRadius; x++) {
                    int idY = Clamp(i + y, 0, height - 1);
                    int idX = Clamp(j + x, 0, tmp - 1);
                    tmp += image[idY * height + idX] * kernel[(y + kernelRadius) * size + x + kernelRadius];
                }
            result[i * size1 + j] = Clamp(round(tmp), 0, 255);
        }*/

    /*int row = (width / ProcNum + 2);
    //int tmp = width - ((width / ProcNum) * (ProcNum - 2) + (width / ProcNum + 1));
    int tmp = width / ProcNum - 1;
    MPI_Datatype mpi_LANE;
    MPI_Type_vector(height, row, width, MPI_INT, &mpi_LANE);
    MPI_Type_commit(&mpi_LANE);

    if (ProcRank == 0)
        parallelTime = MPI_Wtime();

    if (ProcRank == 0)
        for (int i = 1; i < ProcNum; i++)
            MPI_Send((image + ((tmp - i) + (i - 1) * (width / ProcNum + 1))), 1, mpi_LANE, i, 0, MPI_COMM_WORLD);
    else
        MPI_Recv(image, 1, mpi_LANE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Type_free(&mpi_LANE);

    if (ProcRank == 0)
        tmp++;
    else
        tmp = row;

    int* result = new int[height * tmp];
    double qwerty;

    for (int i = 0; i < height; i++)
        for (int j = 0; j < tmp; j++) {
            qwerty = 0;
            for (int y = -kernelRadius; y <= kernelRadius; y++)
                for (int x = -kernelRadius; x <= kernelRadius; x++) {
                    int idY = Clamp(i + y, 0, height - 1);
                    int idX = Clamp(j + x, 0, tmp - 1);
                    qwerty += image[idY * height + idX] * kernel[(y + kernelRadius) * size + x + kernelRadius];
                }
            result[i * tmp + j] = Clamp(round(qwerty), 0, 255);
        }*/

    /*if (ProcRank == 0) {
        tmp--;
        int* res = new int[height * row];
        for (int procs = 1; procs < ProcNum; procs++) {
            MPI_Recv(res, height * row, MPI_INT, procs, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < height; i++)
                for (int j = 0; j < row; j++)
                    image[i * height + j + tmp + (row - 2) * ((procs - 1))] = res[i * row + j + 1];
        }
        for (int i = 0; i < height; i++)
            for (int j = 0; j < tmp; j++)
                image[i * height + j] = result[i * tmp + j];
        int tempbus = tmp;
        tmp++;
        if (ProcNum == 1)
            tempbus++;

        for (int i = 0; i < height; i++)
            for (int j = 0; j < tempbus; j++)
                image[i * height + j] = result[i * tmp + j];
    }
    else
        MPI_Send(result, height * row, MPI_INT, 0, 0, MPI_COMM_WORLD);*/

   /* MPI_Group firstGr, secondGr, WorldGroup;
    MPI_Comm firstComm, secondComm;

    int* ranks = new int[width % ProcNum];
    for (int i = 0; i < width % ProcNum; i++)
        ranks[i] = i + 1;

    MPI_Comm_group(MPI_COMM_WORLD, &WorldGroup);
    MPI_Group_incl(WorldGroup, width % ProcNum, ranks, &secondGr);
    MPI_Group_difference(WorldGroup, secondGr, &firstGr);

    MPI_Comm_create(MPI_COMM_WORLD, firstGr, &firstComm);
    MPI_Comm_create(MPI_COMM_WORLD, secondGr, &secondComm);

    int row = (width / ProcNum + 2);
    int tmp = width / ProcNum - 1;
    
    MPI_Datatype mpi_LANE1;
    MPI_Type_vector(height, row, width, MPI_INT, &mpi_LANE1);
    MPI_Type_commit(&mpi_LANE1);

    MPI_Datatype mpi_LANE2;
    MPI_Type_vector(height, row, width, MPI_INT, &mpi_LANE2);
    MPI_Type_commit(&mpi_LANE2);*/

    if (ProcRank == 0) {  // Print conclusion
        parallelTime = MPI_Wtime() - parallelTime;
        std::cout << "Image made with parallel algorithm:\n";
        printImage(image, width, height);
        equality_check(image, checkImage, height, width);  // Test for equality
        std::cout << std::fixed << "Non-parallel algorithm time: "
            << sequentialTime << "\nParallel algorithm time: "
            << parallelTime << std::endl;
    }

    MPI_Finalize();
    return 0;
}

void printImage(int* image, int width, int height) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++)
            std::cout << image[i * width + j] << ' ';
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void createKernel(double* kernel, int radius, int size, double sigma) {
    double norm = 0;
    for (int i = -radius; i <= radius; i++)
        for (int j = -radius; j <= radius; j++) {
            kernel[(i + radius) * size + (j + radius)] = exp(-(i * i + j * j) / (sigma * sigma));
            norm += kernel[(i + radius) * size + (j + radius)];
        }
    for (int i = 0; i < size; i++)  // Normalization
        for (int j = 0; j < size; j++)
            kernel[i * size + j] /= norm;
    // Print kernel
    std::cout << "Kernel:\n";
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++)
            std::cout << std::fixed << kernel[i * size + j] << ' ';
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void processImage(int* originIm, int* checkImage, double* kernel, int height, int width, int size, int radius) {
    double tmp;
    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++) {
            tmp = 0;
            for (int y = -radius; y <= radius; y++)
                for (int x = -radius; x <= radius; x++) {
                    int idY = Clamp(i + y, 0, height - 1);
                    int idX = Clamp(j + x, 0, width - 1);
                    tmp += originIm[idY * height + idX] * kernel[(y + radius) * size + x + radius];
                }
            checkImage[i * width + j] = Clamp(round(tmp), 0, 255);
        }
}

int Clamp(int value, int min, int max) {
    if (value < min)
        return min;
    if (value > max)
        return max;
    return value;
}

void equality_check(int * res, int * res2, int height, int width) {
    bool flag = false;
    for (int i = 0; i < height*width; i++)
        if (res2[i] - res[i] != 0) {
            flag = !flag;
            break;
        }
    if (!flag)
        std::cout << "The program works correctly.\n";
}
