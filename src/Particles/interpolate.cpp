#include "w.h"
#include "interpolate.h"
#include <array>
#include <cmath>
#include <iterator>
#include <string>
#include <tuple>
#include <vector>

using namespace std;
using namespace Proto;

auto w_ptr = W_2;

std::map<double (*)(double), std::tuple<int, int>> r_map = {
    {W_2, std::make_tuple(0, 1)},
    {W_3, std::make_tuple(-1, 2)},
    {W_4, std::make_tuple(-1, 2)},
    {W_6, std::make_tuple(-2, 3)}};


bool point_is_present_in_vector(const Point &point_to_find, const shared_ptr<vector<Point>> &vec) {
    if (find(vec->begin(), vec->end(), point_to_find) != vec->end()) {
        return true;
    }
    return false;
}

void complete_grid(shared_ptr<vector<Point>> &points_to_consider, const Point &to_add) {
  if (point_is_present_in_vector(to_add, points_to_consider)) {
    // if the point is present, then its points to consider must already be within
    return;
  } else {
    points_to_consider->push_back(to_add);
  }

    Point new_1(to_add[0] + 1, to_add[1]);
    if (!point_is_present_in_vector(new_1, points_to_consider)) {
        points_to_consider->push_back(new_1);
    }

    Point new_2(to_add[0], to_add[1] + 1);
    if (!point_is_present_in_vector(to_add, points_to_consider)) {
      points_to_consider->push_back(new_2);
    }

    Point new_3(to_add[0] + 1, to_add[1] + 1);
    if (!point_is_present_in_vector(to_add, points_to_consider)) {
      points_to_consider->push_back(new_3);
    }
}

shared_ptr<vector<Point>> compute_points_to_consider(const Point &original_point) {
  auto r_tuple = r_map[w_ptr];
  auto r_high = get<1>(r_tuple);
  auto r_low = get<0>(r_tuple);

  auto points_to_consider = make_shared<vector<Point>>();
  
  for (int i = r_low; i < r_high; ++i) {
    Point new_1(original_point[0] + i, original_point[1]);
    complete_grid(points_to_consider, new_1);
    Point new_2(original_point[0], original_point[1] + i);
    complete_grid(points_to_consider, new_2);
    Point new_3(original_point[0] + i, original_point[1] + i);
    complete_grid(points_to_consider, new_3);
  }

  return points_to_consider;
}

// m_dx is h_g
array<array<double, DIM>, DIM> interpolate(const BoxData<double> G_i[DIM][DIM], const Point &original_i, const Particle &x_k, const double &h_g) {
    array<array<double, DIM>, DIM> G_k;
    Point new_i(floor(original_i[0]/h_g), floor(original_i[0]/h_g));
    auto grid_points = compute_points_to_consider(new_i);
    for (auto i : *grid_points) {

        Point x_bar(i[0] * h_g, i[1] * h_g);
        Point z(x_k.m_x[0] - x_bar[0], x_k.m_x[1] - x_bar[1]);
        // w_ptr = W_2 function
        double W_value = W(z, h_g, w_ptr);

        // Multiply G_i by W_value and incrememnt G_k matrix by G_i
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                G_k[i][j] += G_i[i][j](Point(0, 0)) * W_value;
            }
        }
    }

    return G_k;
}
