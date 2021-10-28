#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>
#include <time.h>

void random_matrix(double *matrix, int N);
void print_vector(double *vector, int N);
void print_matrix(double *matrix, int N);
void matrix_vector_mult(double *matrix, double *vector,
                        double *result, int N);
void matrix_matrix_mult(double *mat_1, double *mat_2,
                        double *result, int N);
void matrix_power(double *matrix, double *result,
                  int N, int power);