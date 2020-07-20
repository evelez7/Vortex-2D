#ifndef UTIL_H
#define UTIL_H

#include <memory>
#include <vector>
#include <tuple>

typedef std::vector<std::vector<double>> two_d_matrix;
typedef std::shared_ptr<std::vector<double>> vector_ptr_double;
typedef std::vector<std::tuple<double, std::vector<double>>> tuple_vector;

std::shared_ptr<two_d_matrix> get_rotation_matrix(double);

std::shared_ptr<std::vector<double>> get_unit_vector(int const&);

std::tuple<double, double> get_r_tuple(double (*)(double));

std::tuple<std::vector<double>, std::vector<double>> get_interpolation_limits(std::shared_ptr<std::vector<double>>, double (*)(double));

bool vectors_are_equal(std::vector<double> const&, std::vector<double> const&);

bool limit_reached(std::vector<double> const&, std::vector<double> const&);

std::shared_ptr<std::vector<double>> divide_vector_by_scalar(std::shared_ptr<std::vector<double>> const&, double);

std::shared_ptr<std::vector<double>> multiply_vector_by_scalar(std::shared_ptr<std::vector<double>> const&, double);

std::shared_ptr<std::vector<double>> subtract_vectors(std::shared_ptr<std::vector<double>> const&, std::shared_ptr<std::vector<double>> const&);

std::shared_ptr<std::vector<double>> add_vectors(std::shared_ptr<std::vector<double>> const&, std::shared_ptr<std::vector<double>> const&);

std::shared_ptr<std::vector<double>> multiply_matrix_by_vector(std::shared_ptr<two_d_matrix> const&, std::shared_ptr<std::vector<double>> const&);

void floor_vector(std::shared_ptr<std::vector<double>>&);

void increment_vector(std::shared_ptr<std::vector<double>> &, double);

void print_vector(std::vector<double> const&);

void print_tuple_vector(std::shared_ptr<tuple_vector> const&);

bool vector_is_present_in_matrix(std::shared_ptr<two_d_matrix> const&, std::vector<double> const&);

std::shared_ptr<tuple_vector> compute_input_set(int, double, double, double);

std::shared_ptr<two_d_matrix> compute_points_to_consider(std::shared_ptr<std::vector<double>> const&, double (*)(double));

void write_to_file(double (*)(double), int, double, std::shared_ptr<two_d_matrix> const&);
#endif