#ifndef EIGEN_UTIL_H
#define EIGEN_UTIL_H

#include <memory>
#include <array>
#include <algorithm>
#include <vector>

using namespace std;

array<array<double, DIM>, DIM> get_transpose(const array<array<double, DIM>, DIM>&);
array<array<double, DIM>, DIM> multiply_matrices(const array<array<double, DIM>, DIM>&, const array<array<double, DIM>, DIM>&);
double get_determinant(const array<array<double, DIM>, DIM>&);
array<double, DIM> multiply_matrix_by_vector(const array<array<double, DIM>, DIM>&, const array<double, DIM>&);
array<array<double, DIM>, DIM> get_inverse(const array<array<double, DIM>, DIM>&);
void print_matrix(const array<array<double, DIM>, DIM>&);
array<double, DIM> get_roots(const array<double, 3>&);
array<double, 3> get_characteristic_polynomial(const array<array<double, DIM>, DIM>&);
array<double, DIM> multiply_vector_by_scalar(const array<double, DIM>&, const double&);
bool is_zero_matrix(const array<array<double, DIM>, DIM>&);
void verify_eigenvalues(const array<array<double, DIM>, DIM>&, const array<double, DIM>&);
void verify_eigenvectors_verbose(const array<array<double, DIM>, DIM>&, const array<array<double, DIM>, DIM>&, const array<double, DIM>&);
void collect_eigenvector_errors(const array<array<double, DIM>, DIM>&, const array<array<double, DIM>, DIM>&, const array<double, DIM>&, shared_ptr<vector<double>>&);
void verify_eigenvectors(const array<array<double, DIM>, DIM>&, const array<array<double, DIM>, DIM>&, const array<double, DIM>&);
array<double, DIM> jacobi(const array<array<double, DIM>, DIM>&, const array<double, DIM>&, array<double, DIM>&, const double&);
array<array<double, DIM>, DIM> find_eigenvectors(const array<array<double, DIM>, DIM>&, const array<double, DIM>&);
array<double, DIM> get_sym_eigenvalues(const array<array<double, DIM>, DIM>&);
#endif