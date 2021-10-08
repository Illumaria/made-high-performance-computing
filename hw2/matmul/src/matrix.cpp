#include "matrix.h"

using namespace task;

Matrix::Matrix() : rows(1), cols(1), ptr(new double*[1]) {
  ptr[0] = new double[1];
  ptr[0][0] = 0;
}

Matrix::Matrix(size_t rows, size_t cols)
    : rows(rows), cols(cols), ptr(new double*[rows]) {
  for (size_t i = 0; i < rows; ++i) ptr[i] = new double[cols];

  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) ptr[i][j] = 0;
}

Matrix::Matrix(const Matrix& copy)
    : rows(copy.rows), cols(copy.cols), ptr(new double*[copy.rows]) {
  for (size_t i = 0; i < rows; ++i) ptr[i] = new double[cols];

  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) ptr[i][j] = copy[i][j];
}

Matrix& Matrix::operator=(const Matrix& a) {
  if (&a != this) {
    if (rows != a.rows || cols != a.cols) this->resize(a.rows, a.cols);

    for (size_t i = 0; i < rows; ++i)
      for (size_t j = 0; j < cols; ++j) ptr[i][j] = a[i][j];
  }

  return *this;
}

double& Matrix::get(size_t row, size_t col) {
  if (row < 0 || row >= rows || col < 0 || col >= cols)
    throw OutOfBoundsException();
  return ptr[row][col];
}

const double& Matrix::get(size_t row, size_t col) const {
  if (row < 0 || row >= rows || col < 0 || col >= cols)
    throw OutOfBoundsException();
  return ptr[row][col];
}

void Matrix::set(size_t row, size_t col, const double& value) {
  if (row < 0 || row >= rows || col < 0 || col >= cols)
    throw OutOfBoundsException();
  ptr[row][col] = value;
}

void Matrix::resize(size_t new_rows, size_t new_cols) {
  Matrix result(new_rows, new_cols);

  for (size_t i = 0; i < new_rows; ++i)
    for (size_t j = 0; j < new_cols; ++j)
      if (i >= rows || j >= cols)
        result[i][j] = 0;
      else
        result[i][j] = ptr[i][j];

  for (size_t i = 0; i < rows; ++i) delete[] ptr[i];
  delete[] ptr;

  ptr = result.ptr;
  rows = result.rows;
  cols = result.cols;
}

double* Matrix::operator[](size_t index) { return ptr[index]; }

double* Matrix::operator[](size_t index) const { return ptr[index]; }

Matrix& Matrix::operator+=(const Matrix& a) {
  if (rows != a.rows || cols != a.cols) throw SizeMismatchException();

  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) ptr[i][j] = ptr[i][j] + a[i][j];

  return *this;
}

Matrix& Matrix::operator-=(const Matrix& a) {
  if (rows != a.rows || cols != a.cols) throw SizeMismatchException();

  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) ptr[i][j] = ptr[i][j] - a[i][j];

  return *this;
}

Matrix& Matrix::operator*=(const Matrix& a) {
  if (cols != a.rows) throw SizeMismatchException();

  Matrix result(rows, a.cols);

  for (size_t r = 0; r < cols; ++r)
    for (size_t i = 0; i < rows; ++i)
      for (size_t j = 0; j < a.cols; ++j) result[i][j] += ptr[i][r] * a[r][j];

  for (size_t i = 0; i < rows; ++i) delete[] ptr[i];
  delete[] ptr;
  ptr = result.ptr;
  rows = result.rows;
  cols = result.cols;

  return *this;
}

Matrix& Matrix::operator*=(const double& number) {
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) ptr[i][j] *= number;

  return *this;
}

Matrix Matrix::operator+(const Matrix& a) const {
  if (rows != a.rows || cols != a.cols) throw SizeMismatchException();

  Matrix result(rows, cols);
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) result[i][j] = ptr[i][j] + a[i][j];

  return result;
}

Matrix Matrix::operator-(const Matrix& a) const {
  if (rows != a.rows || cols != a.cols) throw SizeMismatchException();

  Matrix result(rows, cols);
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) result[i][j] = ptr[i][j] - a[i][j];

  return result;
}

Matrix Matrix::operator*(const Matrix& a) const {
  if (cols != a.rows) throw SizeMismatchException();

  Matrix result(rows, a.cols);

  for (size_t r = 0; r < cols; ++r)
    for (size_t i = 0; i < rows; ++i)
      for (size_t j = 0; j < a.cols; ++j) result[i][j] += ptr[i][r] * a[r][j];

  return result;
}

Matrix Matrix::operator*(const double& number) const {
  Matrix result(*this);
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) result[i][j] *= number;

  return result;
}

Matrix Matrix::operator-() const {
  Matrix result(rows, cols);
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) result[i][j] = -ptr[i][j];

  return result;
}

Matrix Matrix::operator+() const {
  Matrix result(*this);

  return result;
}

double Matrix::det() const {
  if (rows != cols) throw SizeMismatchException();

  Matrix temp(*this);

  double coeff, result = 1.0;
  for (size_t row_n = 0; row_n < rows - 1; ++row_n) {
    coeff = temp[row_n][row_n];
    result *= coeff;
    for (size_t col_n = row_n; col_n < cols; ++col_n)
      temp[row_n][col_n] /= coeff;

    for (size_t next_row_n = row_n + 1; next_row_n < rows; ++next_row_n) {
      coeff = temp[next_row_n][row_n] / temp[row_n][row_n];

      for (size_t col_n = 0; col_n < cols; ++col_n)
        temp[next_row_n][col_n] -= coeff * temp[row_n][col_n];
    }
  }

  result *= temp[rows - 1][cols - 1];

  return result;
}

void Matrix::transpose() {
  Matrix result(cols, rows);

  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) result[j][i] = ptr[i][j];

  for (size_t i = 0; i < rows; ++i) delete[] ptr[i];
  delete[] ptr;
  ptr = result.ptr;
  size_t temp = rows;
  rows = cols;
  cols = temp;
}

Matrix Matrix::transposed() const {
  Matrix result(cols, rows);

  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) result[j][i] = ptr[i][j];

  return result;
}

double Matrix::trace() const {
  if (rows != cols) throw SizeMismatchException();

  double result;
  for (size_t i = 0; i < rows; ++i) result += ptr[i][i];

  return result;
}

std::vector<double> Matrix::getRow(size_t row) {
  if (row < 0 || row >= rows) throw OutOfBoundsException();

  std::vector<double> result(cols);
  for (size_t i = 0; i < cols; ++i) result[i] = ptr[row][i];

  return result;
}

std::vector<double> Matrix::getColumn(size_t column) {
  if (column < 0 || column >= cols) throw OutOfBoundsException();

  std::vector<double> result(rows);
  for (size_t i = 0; i < rows; ++i) result[i] = ptr[i][column];

  return result;
}

bool Matrix::operator==(const Matrix& a) const {
  if (rows != a.rows || cols != a.cols) return false;

  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j)
      if ((ptr[i][j] - a[i][j]) * (ptr[i][j] - a[i][j]) > EPS * EPS)
        return false;

  return true;
}

bool Matrix::operator!=(const Matrix& a) const { return !(*this == a); }

const size_t Matrix::getRows() const { return rows; }

const size_t Matrix::getCols() const { return cols; }

Matrix task::operator*(const double& a, const Matrix& b) {
  Matrix result(b);
  return result *= a;
}

std::ostream& task::operator<<(std::ostream& output, const Matrix& matrix) {
  for (size_t i = 0; i < matrix.getRows(); ++i) {
    for (size_t j = 0; j < matrix.getCols(); ++j) output << matrix[i][j] << " ";
    output << "\n";
  }

  return output;
}

std::istream& task::operator>>(std::istream& input, Matrix& matrix) {
  size_t rows, cols;
  input >> rows >> cols;
  if (rows < 0 || cols < 0) return input;
  matrix.resize(rows, cols);

  double value;
  for (size_t i = 0; i < rows; ++i)
    for (size_t j = 0; j < cols; ++j) {
      input >> value;
      matrix[i][j] = value;
    }

  return input;
}

Matrix task::strassen_product(const Matrix& left, const Matrix& right) {
  if (left.getCols() != left.getRows()) throw SizeMismatchException();
  if (right.getCols() != right.getRows()) throw SizeMismatchException();
  if (left.getCols() != right.getRows()) throw SizeMismatchException();

  size_t n = left.getCols();

  if (n <= MATRIX_SIZE_EFFICIENT_BOUND)
    return left * right;

  // pad sub-matrices with odd dimentions
  size_t m = n % 2 ? n / 2 + 1 : n / 2;

  Matrix a11(m, m), a12(m, m), a21(m, m), a22(m, m);
  Matrix b11(m, m), b12(m, m), b21(m, m), b22(m, m);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i < m && j < m) {  // top left part
          a11[i][j] = left[i][j];
          b11[i][j] = right[i][j];
      } else if (i < m && j >= m) {  // top right part
          a12[i][j - m] = left[i][j];
          b12[i][j - m] = right[i][j];
      } else if (i >= m && j < m) {  // bottom left part
          a21[i - m][j] = left[i][j];
          b21[i - m][j] = right[i][j];
      } else if (i >= m && j >= m) {  // bottom right part
          a22[i - m][j - m] = left[i][j];
          b22[i - m][j - m] = right[i][j];
      }
    }
  }

  auto p11 = a11 + a22;
  auto p12 = b11 + b22;
  auto p1 = strassen_product(p11, p12);

  auto p21 = a21 + a22;
  auto p22 = b11;
  auto p2 = strassen_product(p21, p22);

  auto p31 = a11;
  auto p32 = b12 - b22;
  auto p3 = strassen_product(p31, p32);

  auto p41 = a22;
  auto p42 = b21 - b11;
  auto p4 = strassen_product(p41, p42);

  auto p51 = a11 + a12;
  auto p52 = b22;
  auto p5 = strassen_product(p51, p52);

  auto p61 = a21 - a11;
  auto p62 = b11 + b12;
  auto p6 = strassen_product(p61, p62);

  auto p71 = a12 - a22;
  auto p72 = b21 + b22;
  auto p7 = strassen_product(p71, p72);

  auto c11 = p1 + p4 - p5 + p7;
  auto c12 = p3 + p5;
  auto c21 = p2 + p4;
  auto c22 = p1 - p2 + p3 + p6;

  Matrix c(n, n);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i < m && j < m) {  // top left part
          c[i][j] = c11[i][j];
      } else if (i < m && j >= m) {  // top right part
          c[i][j] = c12[i][j - m];
      } else if (i >= m && j < m) {  // bottom left part
          c[i][j] = c21[i - m][j];
      } else if (i >= m && j >= m) {  // bottom right part
          c[i][j] = c22[i - m][j - m];
      }
    }
  }

  return c;
}