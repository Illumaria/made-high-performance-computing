#pragma once

#include <iostream>
#include <vector>

namespace task {

const double EPS = 1e-6;
const size_t MATRIX_SIZE_EFFICIENT_BOUND = 64;

class OutOfBoundsException : public std::exception {};
class SizeMismatchException : public std::exception {};

class Matrix {
 public:
  Matrix();
  Matrix(size_t rows, size_t cols);
  Matrix(const Matrix& copy);
  Matrix& operator=(const Matrix& a);

  const size_t getRows() const;
  const size_t getCols() const;

  double& get(size_t row, size_t col);
  const double& get(size_t row, size_t col) const;
  void set(size_t row, size_t col, const double& value);
  void resize(size_t new_rows, size_t new_cols);

  double* operator[](size_t row);
  double* operator[](size_t row) const;

  Matrix& operator+=(const Matrix& a);
  Matrix& operator-=(const Matrix& a);
  Matrix& operator*=(const Matrix& a);
  Matrix& operator*=(const double& number);

  Matrix operator+(const Matrix& a) const;
  Matrix operator-(const Matrix& a) const;
  Matrix operator*(const Matrix& a) const;
  Matrix operator*(const double& a) const;

  Matrix operator-() const;
  Matrix operator+() const;

  double det() const;
  void transpose();
  Matrix transposed() const;
  double trace() const;

  std::vector<double> getRow(size_t row);
  std::vector<double> getColumn(size_t column);

  bool operator==(const Matrix& a) const;
  bool operator!=(const Matrix& a) const;

  friend Matrix operator*(const double&, const Matrix&);
  friend std::ostream& operator<<(std::ostream&, const Matrix&);
  friend std::istream& operator>>(std::istream&, Matrix&);

 private:
  double** ptr;
  size_t rows, cols;
};

Matrix operator*(const double& a, const Matrix& b);

std::ostream& operator<<(std::ostream& output, const Matrix& matrix);
std::istream& operator>>(std::istream& input, Matrix& matrix);

Matrix strassen_product(const Matrix& left, const Matrix& right);

}  // namespace task