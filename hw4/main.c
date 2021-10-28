#include "src/matrix.h"
#include "src/ranking.h"

int main(int argc, char *argv[]) {
    if (argc < 3 || atoi(argv[1]) < 1 || atoi(argv[2]) < 1) {
        return 0;
    }

    int N = atoi(argv[1]);
    int power = atoi(argv[2]);
    printf("N = %d, power = %d\n", N, power);

    double *matrix = (double *) malloc(N * N * sizeof(double));
    random_matrix(matrix, N);
    printf("Generated adjacency matrix:\n");
    print_matrix(matrix, N);

    double *powered = (double *) malloc(N * N * sizeof(double));
    matrix_power(matrix, powered, N, power);
    printf("Generated adjacency matrix powered by %d:\n", power);
    print_matrix(powered, N);

    double *naive_ranking_result = (double *) malloc(N * sizeof(double));
    naive_ranking(matrix, naive_ranking_result, N);

    double *pagerank_result = (double *) malloc(N * sizeof(double));
    double tolerance = 1e-6;
    double damping_factor = 0.85;
    int num_iterations = 1e6;
    pagerank(matrix, pagerank_result, N, tolerance, damping_factor, num_iterations);

    printf("Ranking results:\n");
    for (int i = 0; i < N; ++i) {
        printf("%2d: naive_ranking=%.4lf, pagerank=%.4lf\n",
               i, naive_ranking_result[i], pagerank_result[i]);
    }

    free(matrix);
    free(powered);
    free(naive_ranking_result);
    free(pagerank_result);

    return 0;
}