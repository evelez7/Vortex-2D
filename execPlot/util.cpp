#include "util.h"
#include "w.h"
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>

std::shared_ptr<two_d_matrix> create_grid(std::shared_ptr<std::vector<double>> const&);
void add_to_grid(std::shared_ptr<two_d_matrix>, std::shared_ptr<std::vector<double>>);

std::map<double, std::string> time_map = {
    {0.0, "0"},
    {1.0 / 6.0, "0.1667"},
    {1.0 / 4.0, "0.25"}};

std::map<double (*)(double), std::tuple<int, int>> r_map = {
    {W_2, std::make_tuple(0, 1)},
    {W_3, std::make_tuple(-1, 2)},
    {W_4, std::make_tuple(-1, 2)},
    {W_6, std::make_tuple(-2, 3)}
};

std::shared_ptr<two_d_matrix> get_rotation_matrix(double time) {
    auto rotation_matrix = std::make_shared<two_d_matrix>();

    std::vector<double> first_row{cos(M_PI * time), sin(M_PI * time)};
    std::vector<double> second_row{-1 * sin(M_PI * time), cos(M_PI * time)};

    rotation_matrix->push_back(first_row);
    rotation_matrix->push_back(second_row);

    return rotation_matrix;
}

std::tuple<double, double> get_r_tuple(double (*W_pointer)(double)) {
    return r_map[W_pointer];
}

std::shared_ptr<std::vector<double>> get_unit_vector(int const& dimensions) {
    auto new_unit_vector = std::make_shared<std::vector<double>>();

    for (int i = 0; i < dimensions; ++i) {
        new_unit_vector->push_back(1);
    }

    return new_unit_vector;
}

std::tuple<std::vector<double>, std::vector<double>> get_interpolation_limits(std::shared_ptr<std::vector<double>> i, double (*w_pointer)(double)) {
    auto r_double = r_map[w_pointer];
    std::vector<double> vec_low, vec_high;

    auto unit_vector = get_unit_vector(i->size());
    for (int j = 0; j < i->size(); ++j) {
        vec_low.push_back(i->at(j) + (std::get<0>(r_double) * unit_vector->at(j)));
        vec_high.push_back(i->at(j) + (std::get<1>(r_double) * unit_vector->at(j)));
    }
    return std::make_tuple(vec_low, vec_high); 
}

bool vectors_are_equal(std::vector<double> const& first, std::vector<double> const& second) {
    bool equal = true;  

    if (first.size() != second.size()) {
        throw std::exception();
    }

    for (int i = 0; i < first.size(); ++i) {
        if (first.at(i) != second.at(i)) {
            equal = false;
        }
    }

    return equal;
}

bool vector_is_present_in_matrix(std::shared_ptr<two_d_matrix> const& matrix, std::vector<double> const& vec) {
    if (std::find(matrix->begin(), matrix->end(), vec) != matrix->end()) {
        return true;
    }
    return false;
}

bool limit_reached(std::vector<double> const& first, std::vector<double> const& second) {
    bool limit_passed = false;

    for (int i = 0; i < first.size(); ++i) {
        if (first.at(i) > second.at(i)) {
            limit_passed = true;
        }
    }

    return limit_passed || vectors_are_equal(first, second);

}

std::shared_ptr<std::vector<double>> divide_vector_by_scalar(std::shared_ptr<std::vector<double>> const& vec, double scalar) {
    auto new_vector = std::make_shared<std::vector<double>>();
    for (int i = 0; i < vec->size(); ++i) {
        new_vector->push_back(vec->at(i) / scalar);
    }

    return new_vector;
}

std::shared_ptr<std::vector<double>> multiply_vector_by_scalar(std::shared_ptr<std::vector<double>> const& vec, double scalar) {
    auto new_vector = std::make_shared<std::vector<double>>();
    for (int i = 0; i < vec->size(); ++i) {
        new_vector->push_back(vec->at(i)*scalar);
    }

    return new_vector;
}

std::shared_ptr<std::vector<double>> subtract_vectors(std::shared_ptr<std::vector<double>> const& first, std::shared_ptr<std::vector<double>> const& second) {
    auto new_vector = std::make_shared<std::vector<double>>();
    if (first->size() != second->size()) {
        std::cout << "SUBTRACT VECTORS EXCEPTION" << std::endl;
        throw std::exception();
    }
    for (int i = 0; i < first->size(); ++i) {
       new_vector->push_back(first->at(i) - second->at(i)); 
    }
    return new_vector;
}

std::shared_ptr<std::vector<double>> add_vectors(std::shared_ptr<std::vector<double>> const& first, std::shared_ptr<std::vector<double>> const& second) {
    auto new_vector = std::make_shared<std::vector<double>>();
    if (first->size() != second->size()) {
        std::cout << "SUBTRACT VECTORS EXCEPTION" << std::endl;
        throw std::exception();
    }
    for (int i = 0; i < first->size(); ++i) {
        new_vector->push_back(first->at(i) + second->at(i));
    }
    return new_vector;
}

std::shared_ptr<std::vector<double>> multiply_matrix_by_vector(std::shared_ptr<two_d_matrix> const& matrix, std::shared_ptr<std::vector<double>> const& vec) {
    auto new_vector = std::make_shared<std::vector<double>>();
    if (vec->size() != matrix->size()) {
        throw std::exception();
    }

    for (auto row : *matrix) {
        double sum = 0;
        for (int i = 0; i < row.size(); ++i) {
           sum +=  row.at(i) * vec->at(i);
        }
        new_vector->push_back(sum);
    }
    return new_vector;
}

void increment_vector(std::shared_ptr<std::vector<double>> & to_increment, double value) {
    for (int i = 0; i < to_increment->size(); ++i) {
        to_increment->at(i)+=value;
    }
}

void floor_vector(std::shared_ptr<std::vector<double>> & vector_to_floor) {
    for (int i = 0; i < vector_to_floor->size(); ++i) {
        vector_to_floor->at(i) = std::floor(vector_to_floor->at(i));
    }
}

void print_vector(std::vector<double> const& vec) {
    std::cout << "(";
    for (int i = 0; i < vec.size(); ++i) {
        std::cout << vec.at(i);
        if (i != vec.size() - 1) {
            std::cout << ",";
        }
    }
    std::cout << ")" << std::endl;
}

void print_tuple_vector(std::shared_ptr<tuple_vector> const& to_print) {
    for (int i = 0; i < to_print->size(); ++i) {
        std::cout << "f_k: " << std::get<0>(to_print->at(i)) << "   vector: ";
        auto vec = std::get<1>(to_print->at(i));
        for (auto entry : vec) {
            std::cout << entry << " ";
        }
        std::cout<< " size: " << to_print->size() << std::endl;
    }
}

std::shared_ptr<two_d_matrix> create_grid(std::shared_ptr<std::vector<double>> const &left_hand_corner) {
    auto grid = std::make_shared<two_d_matrix>();

    grid->push_back(*left_hand_corner);

    grid->push_back(std::vector<double>{left_hand_corner->at(0) + 1, left_hand_corner->at(1)});
    grid->push_back(std::vector<double>{left_hand_corner->at(0) + 1, left_hand_corner->at(1) + 1});
    grid->push_back(std::vector<double>{left_hand_corner->at(0), left_hand_corner->at(1) + 1});
    return grid;
}

void add_to_grid(std::shared_ptr<two_d_matrix> grid, std::shared_ptr<std::vector<double>> to_add) {
    grid->push_back(*to_add);

    std::vector<double> new_vector_x_1{to_add->at(0) + 1, to_add->at(1)};
    if (!vector_is_present_in_matrix(grid, new_vector_x_1)) {
        grid->push_back(new_vector_x_1);
    }

    std::vector<double> new_vector_x_1_y_1{to_add->at(0) + 1, to_add->at(1) + 1};
    if (!vector_is_present_in_matrix(grid, new_vector_x_1_y_1)) {
        grid->push_back(new_vector_x_1_y_1);
    }

    std::vector<double> new_vector_y_1{to_add->at(0), to_add->at(1) + 1};
    if (!vector_is_present_in_matrix(grid, new_vector_y_1)) {
        grid->push_back(new_vector_y_1);
    }
}

std::shared_ptr<tuple_vector> compute_input_set(int N, double h, double h_g, double time) {
    auto input_set = std::make_shared<std::vector<std::tuple<double, std::vector<double>>>>();
    for (int i_count = 1; i_count <= N-1; ++i_count) {
        for (int j_count = 1; j_count <= N-1; ++j_count) {
            auto unit_vector = get_unit_vector(2);
            std::vector<double> i_raw{static_cast<double>(i_count), static_cast<double>(j_count)};
            auto i = std::make_shared<std::vector<double>>(i_raw);
            auto parens_left = multiply_vector_by_scalar(i, h);
            auto parens_right = multiply_vector_by_scalar(unit_vector, 0.25);
            auto parens = subtract_vectors(parens_left, parens_right);
            auto product = multiply_matrix_by_vector(get_rotation_matrix(time), parens);
            auto right_side = multiply_vector_by_scalar(unit_vector, 0.5);

            auto new_vector = add_vectors(product, right_side);
            auto new_f_k = (pow(h, 2.0));
            input_set->emplace_back(new_f_k, *new_vector);
        }
    }
    return input_set;
}

std::shared_ptr<two_d_matrix> compute_points_to_consider(std::shared_ptr<std::vector<double>> const &initial_point, double (*W_pointer)(double)) {
    if (initial_point->size() != 2) {
        throw std::exception();
    }
    auto points_to_consider = create_grid(initial_point);
    auto r_tuple = get_r_tuple(W_pointer);
    double r_low = std::get<0>(r_tuple);
    double r_high = std::get<1>(r_tuple);

    int initial = static_cast<int>(r_low);
    int end = static_cast<int>(r_high);
    for (int i = initial; i <= end; ++i) {
        std::vector<double> new_vector{initial_point->at(0) + i, initial_point->at(1)};
        if (vector_is_present_in_matrix(points_to_consider, new_vector)) {
            continue;
        }
        add_to_grid(points_to_consider, std::make_shared<std::vector<double>>(new_vector));
    }

    for (int i = initial; i <= end; ++i) {
        std::vector<double> new_vector{initial_point->at(0), initial_point->at(1) + i};
        if (vector_is_present_in_matrix(points_to_consider, new_vector)) {
            continue;
        }
        add_to_grid(points_to_consider, std::make_shared<std::vector<double>>(new_vector));
    }

    for (int i = initial; i <= end; ++i) {
        std::vector<double> new_vector{initial_point->at(0) + i, initial_point->at(1) + i};
        if (vector_is_present_in_matrix(points_to_consider, new_vector)) {
            continue;
        }
        add_to_grid(points_to_consider, std::make_shared<std::vector<double>>(new_vector));
    }
    return points_to_consider;
}

void write_to_file(double (*W_script)(double), int N, double time, std::shared_ptr<two_d_matrix> const &answers) {
    std::ofstream output_file("/Users/evelez/code/lbl/2D-Particles/output/" + get_W_name(W_script) + "-" + time_map[time] + "-" + std::to_string(N) + ".txt");
    if (output_file.is_open()) {
        for (auto row : *answers) {
            for (int i = 0; i < N; ++i) {
                output_file << row.at(i);
                if (i != N - 1) {
                    output_file << ',';
                }
            }
            output_file << '\n';
        }
        output_file.flush();
        output_file.close();
    }
}