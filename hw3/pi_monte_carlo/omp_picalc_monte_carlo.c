#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

const unsigned int NUM_POINTS = 100000000;
const unsigned int RADIUS = 1.0;

bool is_inside_circle(double x, double y, double radius) {
    return (x * x + y * y) < radius;
}

void seedThreads(unsigned int* seeds) {
    unsigned int seed, thread_id;

#pragma omp parallel private(seed, thread_id)
    {
        thread_id = omp_get_thread_num();
        seed = (unsigned) time(NULL);
        seeds[thread_id] = (seed & 0xFFFFFFF0) | (thread_id + 1);
        printf("thread %d has seed %u\n", thread_id, seeds[thread_id]);
    }
}

int main (int argc, char *argv[]) {
    unsigned int seeds[omp_get_num_threads()];
    seedThreads(seeds);

    unsigned int correct_points = 0;

    double start = omp_get_wtime();

#pragma omp parallel for shared(seeds) reduction(+:correct_points)
    for (size_t i = 0; i < NUM_POINTS; ++i) {
        double x = ((double)rand_r(&seeds[omp_get_thread_num()]) / (double)RAND_MAX);
        double y = ((double)rand_r(&seeds[omp_get_thread_num()]) / (double)RAND_MAX);

        correct_points += is_inside_circle(x, y, RADIUS);
    }

    double end = omp_get_wtime();

    printf("pi = %f (time: %f s, points: %u)\n",
           (double)(4.0 * correct_points / NUM_POINTS), end - start, NUM_POINTS);

    return 0;
}