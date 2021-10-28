#include "matrix.h"

void random_matrix(double *matrix, int N) {
    unsigned int seed = (unsigned) time(NULL);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            matrix[i * N + j] = (i == j) ? 0 : rand_r(&seed) & 1;
        }
    }
}

void print_vector(double *vector, int N) {
    for (int i = 0; i < N; ++i) {
        printf("%lf ", vector[i]);
    }
    printf("\n");
}

void print_matrix(double *matrix, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            printf("%.0lf ", matrix[i * N + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void matrix_vector_mult(double *matrix, double *vector,
                        double *result, int N) {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N,
                1.0, matrix, N, vector, 1, 0.0, result, 1);
}

void matrix_matrix_mult(double *mat_1, double *mat_2,
                        double *result, int N) {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N,
                1.0, mat_1, N, mat_2, N, 0.0, result, N);
}

void matrix_power(double *matrix, double *result, int N, int power) {
    /*
     * https://www.hackerearth.com/practice/notes/matrix-exponentiation-1/
     */

    // result = identity matrix
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            result[i * N + j] = (i == j) ? 1 : 0;
        }
    }

    double *buffer = (double *) malloc(N * N * sizeof(double));
    double *matrix_pow = (double *) malloc(N * N * sizeof(double));

    // initialize: matrix_pow = matrix
    memcpy(matrix_pow, matrix, N * N * sizeof(double));

    while (power > 0) {
        if (power & 1) {
            // result *= matrix_pow
            matrix_matrix_mult(matrix_pow, result, buffer, N);
            memcpy(result, buffer, N * N * sizeof(double));
        }

        // matrix_pow *= matrix_pow
        matrix_matrix_mult(matrix_pow, matrix_pow, buffer, N);
        memcpy(matrix_pow, buffer, N * N * sizeof(double));

        power >>= 1;
    }

    free(buffer);
    free(matrix_pow);
}