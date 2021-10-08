#include <chrono>
#include <random>

#include "src/matrix.h"

using task::Matrix;

double RandomDouble() {
  static std::mt19937 rand(std::random_device{}());

  std::uniform_real_distribution<double> dist{-10., 10.};
  return dist(rand);
}

Matrix RandomMatrix(size_t rows, size_t cols) {
  Matrix temp(rows, cols);
  for (size_t row = 0; row < rows; ++row) {
    for (size_t col = 0; col < cols; ++col) {
      temp[row][col] = RandomDouble();
    }
  }
  return temp;
}

int main(int argc, char** argv) {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    size_t N = std::stoi(argv[1]);
    std::cout << "N: " << N << std::flush;

    auto matrix_1 = RandomMatrix(N, N);
    auto matrix_2 = RandomMatrix(N, N);

    auto start = high_resolution_clock::now();
    auto result = matrix_1 * matrix_2;
    auto end = high_resolution_clock::now();
    duration<double, std::milli> classic_time = end - start;
    std::cout << "; CLASSIC: " << classic_time.count() << " ms" << std::flush;

    start = high_resolution_clock::now();
    auto result_strassen = strassen_product(matrix_1, matrix_2);
    end = high_resolution_clock::now();
    duration<double, std::milli> strassen_time = end - start;
    std::cout << "; STRASSEN: " << strassen_time.count() << " ms\n";

    return 0;
}