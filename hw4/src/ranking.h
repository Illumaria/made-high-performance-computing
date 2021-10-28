#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <omp.h>
#include "matrix.h"

double l1_distance(double *vec_1, double *vec_2, int N);
void naive_ranking(double *matrix, double *result, int N);
void pagerank(double *matrix, double *result, int N,
              double tolerance, double damping_factor, int max_iter);