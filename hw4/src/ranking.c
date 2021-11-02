#include "ranking.h"

double l1_distance(double *vec_1, double *vec_2, int N) {
    double result = 0;

    int i;

#pragma omp parallel private(i) shared(result, vec_1, vec_2, N)
    {
#pragma omp for reduction(+:result)
        for (i = 0; i < N; ++i) {
            result += fabs(vec_1[i] - vec_2[i]);
        }
    }

    return result;
}

void naive_ranking(double *matrix, double *result, int N) {
    int i, j;

    double total_links = 0;

#pragma omp parallel private(i, j) shared(result, total_links, N)
    {
#pragma omp for reduction(+:total_links)
        for (j = 0; j < N; ++j) {
            result[j] = 0;

            for (i = 0; i < N; ++i) {
                result[j] += matrix[i * N + j];
            }

            total_links += result[j];
        }

#pragma omp for
        for (i = 0; i < N; ++i) {
            result[i] /= total_links;
        }
    }
}

void pagerank(double *matrix, double *pagerank_vector, int N,
              double tolerance, double damping_factor, int max_iter) {
    /*
     * Iterative method from
     * https://en.wikipedia.org/wiki/PageRank#Computation
     */

    int i, j;
    double *fixed_matrix = (double *) malloc(N * N * sizeof(double));

#pragma omp parallel private(i, j) shared(matrix, fixed_matrix, N)
    {
#pragma omp for
        for (i = 0; i < N; ++i) {
            double num_outbound_links = 0;

            for (j = 0; j < N; ++j) {
                num_outbound_links += matrix[i * N + j];
            }

            for (j = 0; j < N; ++j) {
                fixed_matrix[j * N + i] = (num_outbound_links == 0) ? 1.0 / N : matrix[i * N + j] / num_outbound_links;
            }
        }
    }

    double *prev_pagerank_vector = (double *) malloc(N * sizeof(double));

    // initial probability distribution
    for (i = 0; i < N; ++i) {
        prev_pagerank_vector[i] = 1.0 / N;
    }

    for (i = 0; i < max_iter; ++i) {
        // X = M * R
        matrix_vector_mult(fixed_matrix, prev_pagerank_vector, pagerank_vector, N);

        // R' = d * X + (1 - d) / N
        for (j = 0; j < N; ++j) {
            pagerank_vector[j] = damping_factor * pagerank_vector[j] + (1 - damping_factor) / N;
        }

        // check for convergence
        double l1_dist = l1_distance(prev_pagerank_vector, pagerank_vector, N);
        if (l1_dist < tolerance) break;

        memcpy(prev_pagerank_vector, pagerank_vector, N * sizeof(double));
    }

    free(fixed_matrix);
    free(prev_pagerank_vector);
}
