#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define SM (64 / sizeof (double))

const size_t N = 1024;

void ZeroMatrix(double* A, size_t N) {
    for (size_t i = 0; i < N; ++i)
        for(size_t j = 0; j < N; ++j)
            A[i * N + j] = 0.0;
}

void RandomMatrix(double* A, size_t N) {
    srand(time(NULL));

    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            A[i * N + j] = rand() / RAND_MAX;
}

double CalcMatMulTime_ijk(double* A, double* B, double* C, size_t N) {
    struct timeval start, end;
    double r_time = 0.0;

    ZeroMatrix(&C[0], N);

    gettimeofday(&start, NULL);
    
    for (size_t i = 0; i < N; ++i)
        for(size_t j = 0; j < N; ++j)
            for(size_t k = 0; k < N; ++k)
                C[i * N + j] += A[i * N + k] * B[k * N + j];

    gettimeofday(&end, NULL);
    
    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;
    
    return r_time;
}

double CalcMatMulTime_jik(double* A, double* B, double* C, size_t N) {
    struct timeval start, end;
    double r_time = 0.0;

    ZeroMatrix(&C[0], N);

    gettimeofday(&start, NULL);

    for (size_t j = 0; j < N; ++j)
        for (size_t i = 0; i < N; ++i)
            for (size_t k = 0; k < N; ++k)
                C[i * N + j] += A[i * N + k] * B[k * N + j];

    gettimeofday(&end, NULL);
    
    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;
    
    return r_time;
}

double CalcMatMulTime_kij(double* A, double* B, double* C, size_t N) {
    struct timeval start, end;
    double r_time = 0.0;

    ZeroMatrix(&C[0], N);

    gettimeofday(&start, NULL);

    for (size_t k = 0; k < N; ++k) 
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j)
                C[i * N + j] += A[i * N + k] * B[k * N + j];

    gettimeofday(&end, NULL);

    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;

    return r_time;
}

double CalcMatMulTime_kij_opt(double* A, double* B, double* C, size_t N) {
    struct timeval start, end;
    double r_time = 0.0;

    size_t dummy = 0;

    ZeroMatrix(&C[0], N);

    gettimeofday(&start, NULL);

    for (size_t k = 0; k < N; ++k)
        for(size_t i = 0; i < N; ++i) {
            dummy = i * N;
            for(size_t j = 0; j < N; ++j)
                C[dummy + j] += A[dummy + k] * B[k * N + j];
        }

    gettimeofday(&end, NULL);

    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;

    return r_time;
}

double CalcMatMulTime_ijk_transposed(double* A, double* B, double* C, size_t N) {
    struct timeval start, end;
    double r_time = 0.0;

    size_t dummy = 0;
    size_t dummy_ext = 0;

    ZeroMatrix(&C[0], N);

    gettimeofday(&start, NULL);

    double* tmp = (double*) malloc(N * N * sizeof(double));
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            tmp[i * N + j] = B[j * N + i];
    for (size_t i = 0; i < N; ++i) {
        dummy_ext = i * N;
        for (size_t j = 0; j < N; ++j) {
            dummy = j * N;
            for (size_t k = 0; k < N; ++k)
                C[dummy_ext + j] += A[dummy_ext + k] * tmp[dummy + k];
        }
    }

    gettimeofday(&end, NULL);

    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;

    return r_time;
}

double CalcMatMulTime_ijk_transposed_sub(double* A, double* B, double* C, size_t N) {
    struct timeval start, end;
    double r_time = 0.0;

    size_t i, j, k, i2, j2, k2;
    double* rmul1, * rmul2, * rres;
    size_t dummy = 0;
    size_t dummy_ext = 0;

    ZeroMatrix(&C[0], N);

    gettimeofday(&start, NULL);

    for (i = 0; i < N; i += SM)
        for (j = 0; j < N; j += SM)
            for (k = 0; k < N; k += SM)
                for (i2 = 0, rres = &C[i * N + j],
                     rmul1 = &A[i * N + k]; i2 < SM;
                     ++i2, rres += N, rmul1 += N)
                    for (k2 = 0, rmul2 = &B[k * N + j];
                         k2 < SM; ++k2, rmul2 += N)
                        for (j2 = 0; j2 < SM; ++j2)
                            rres[j2] += rmul1[k2] * rmul2[j2];

    gettimeofday(&end, NULL);

    r_time = end.tv_sec - start.tv_sec + ((double) (end.tv_usec - start.tv_usec)) / 1000000;

    return r_time;
}

int main() {
    int NRuns = 5;
    size_t i, j, k;

    double* runtimes;
    double* A, * B, * C;

    A = (double*) malloc(N * N * sizeof(double));
    B = (double*) malloc(N * N * sizeof(double));
    C = (double*) malloc(N * N * sizeof(double));
    runtimes = (double*) malloc(NRuns * sizeof(double));

    RandomMatrix(&A[0], N);
    RandomMatrix(&B[0], N);

    // ijk ordering
    double average_runtime = 0.0;
    for (int n = 0; n < NRuns; ++n) {
        runtimes[n] = CalcMatMulTime_ijk(&A[0], &B[0], &C[0], N);
        printf("runtime %lf seconds\n", runtimes[n]);
        average_runtime += runtimes[n] / NRuns;
    }

    printf("average runtime ijk %lf seconds\n", average_runtime);
    printf("---------------------------------\n");

    // jik ordering
    average_runtime = 0.0;
    for (int n = 0; n < NRuns; ++n) {
        runtimes[n] = CalcMatMulTime_jik(&A[0], &B[0], &C[0], N);
        printf("runtime %lf seconds\n", runtimes[n]);
        average_runtime += runtimes[n] / NRuns;
    }

    printf("average runtime jik %lf seconds\n", average_runtime);
    printf("---------------------------------\n");

    // kij ordering
    average_runtime = 0.0;
    for (int n = 0; n < NRuns; ++n) {
        runtimes[n] = CalcMatMulTime_kij(&A[0], &B[0], &C[0], N);
        printf("runtime %lf seconds\n", runtimes[n]);
        average_runtime += runtimes[n] / NRuns;
    }
    printf("average runtime kij %lf seconds\n", average_runtime);
    printf("---------------------------------\n");

    // kij ordering naive optimization (useless for -O3)
    average_runtime = 0.0;
    for (int n = 0; n < NRuns; ++n) {
        runtimes[n] = CalcMatMulTime_kij_opt(&A[0], &B[0], &C[0], N);
        printf("runtime %lf seconds\n", runtimes[n]);
        average_runtime += runtimes[n] / NRuns;
    }
    printf("average runtime kij opt %lf seconds\n", average_runtime);
    printf("---------------------------------\n");

    // custom optimizations
    average_runtime = 0.0;
    for (int n = 0; n < NRuns; ++n) {
        runtimes[n] = CalcMatMulTime_ijk_transposed(&A[0], &B[0], &C[0], N);
        printf("runtime %lf seconds\n", runtimes[n]);
        average_runtime += runtimes[n] / NRuns;
    }
    printf("average runtime ijk transposed %lf seconds\n", average_runtime);
    printf("---------------------------------\n");

    average_runtime = 0.0;
    for (int n = 0; n < NRuns; ++n) {
        runtimes[n] = CalcMatMulTime_ijk_transposed_sub(&A[0], &B[0], &C[0], N);
        printf("runtime %lf seconds\n", runtimes[n]);
        average_runtime += runtimes[n] / NRuns;
    }
    printf("average runtime ijk transposed submatrix %lf seconds\n", average_runtime);
    printf("---------------------------------\n");

    free(A); 
    free(B);
    free(C);
    return 0;
}
