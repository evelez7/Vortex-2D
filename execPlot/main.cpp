#include "util.h"
#include "w.h"
//#include "VisitWriter.H"
#include "CutoffKernel.H"
#include "ParticleSet.H"
#include "ParticleVelocities.H"
#include "ParticleWriter.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_Point.H"
#include "Proto_RK4.H"
#include "Proto_WriteBoxData.H"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#define DIMENSIONS = 2

int N_array[2] = {32, 64};
double (*W_array[2])(double) = {W_2, W_4};
std::vector<double> time_array{0.0, 1.0 / 6.0, 1.0 / 4.0};

// forward declarations
void I_p(std::shared_ptr<std::vector<std::vector<double>>> &, std::tuple<double, std::vector<double>> const &, double const &, double (*)(double), int);
double f(double);
std::shared_ptr<tuple_vector> compute_input_set(int, double, double, double);
std::shared_ptr<std::vector<std::vector<double>>> work(double (*)(double), int, double, double, double);
void write_to_file(double (*)(double), int, double, std::shared_ptr<std::vector<std::vector<double>>> const &);

int main(int argc, char **argv) {
    for (auto W_ptr : W_array) {
        for (auto N : N_array) {
            for (auto time : time_array) {
                double h = 1.0 / (2.0 * static_cast<double>(N));
                double h_g = 2* h;
                auto interpolation_results = work(W_ptr, N, time, h, h_g);
                auto new_cube = Box::Cube(N);
                auto interpolation_box_data = BoxData<double>(new_cube);

                for (int i = 0; i <= N - 1; ++i) {
                    for (int j = 0; j <= N - 1; ++j) {
                        interpolation_box_data(Point(i, j)) = interpolation_results->at(i).at(j);
                    }
                }
                std::string file_name = "inteprolation_data";
                std::string var_name = "./output/" + get_W_name(W_ptr) + "_" +std::to_string(N) + "_" + std::to_string(time);
                WriteData(interpolation_box_data, 0, h_g, file_name, var_name);
            }
        }
    }

    return 0;
}

std::shared_ptr<std::vector<std::vector<double>>> work(double (*W_pointer)(double), int N, double time, double h, double h_g) {
    auto answers = std::make_shared<std::vector<std::vector<double>>>();
    answers->resize(N);
    if (answers->size() != N) {
        std::cout << "!!!" << std::endl;
        throw std::exception();
    }

    for (int i = 0; i < answers->size(); ++i) {
        answers->at(i).resize(N);
    }

    auto input_set = compute_input_set(N, h, h_g, time);
    for (auto tuple : *input_set) {
        I_p(answers, tuple, h_g, W_pointer, N);
    }
    return answers;
}

void I_p(std::shared_ptr<std::vector<std::vector<double>>> &answers, std::tuple<double, std::vector<double>> const &input_tuple, double const &h_g, double (*w_pointer)(double), int N) {
    auto x_k = std::make_shared<std::vector<double>>(std::get<1>(input_tuple));
    auto left_hand_corner = divide_vector_by_scalar(x_k, h_g);
    floor_vector(left_hand_corner);
    auto grid_points = compute_points_to_consider(left_hand_corner, w_pointer);

    double f_k = std::get<0>(input_tuple);
    for (auto j : *grid_points) {
        auto x_bar = multiply_vector_by_scalar(std::make_shared<std::vector<double>>(j), h_g);
        auto z = subtract_vectors(x_bar, x_k);
        double W_value = W(z, h_g, w_pointer);
        double interpolation_value = W_value * f_k;
        auto x = j.at(0);
        auto y = j.at(1);
        if (x < 0.0 || y < 0.0 || x >= N || y >= N) {
            continue;
        }

        answers->at(x).at(y) += interpolation_value;
    }
}

double f(double h) {
    return pow(h, 2.0);
}