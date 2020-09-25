#include "interpolate.h"
#include "w.h"
#include <array>
#include <cmath>
#include <iterator>
#include <string>
#include <tuple>
#include <array>
#include <vector>
#include <algorithm>

using namespace std;
using namespace Proto;

auto w_ptr = W_2;

std::map<double (*)(double), std::tuple<int, int>> r_map = {
    {W_2, std::make_tuple(0, 1)},
    {W_3, std::make_tuple(-1, 2)},
    {W_4, std::make_tuple(-1, 2)},
    {W_6, std::make_tuple(-2, 3)}};

bool point_is_present_in_vector(const Point& point_to_find,
                                const vector<Point>& vec) {
  if (find(vec.begin(), vec.end(), point_to_find) != vec.end()) {
    return true;
  }
  return false;
}

void complete_grid(vector<Point>& points_to_consider, const Point& to_add) {
  if (!point_is_present_in_vector(to_add, points_to_consider)) {
    // if the point is present, then its points to consider must already be within
    points_to_consider.push_back(to_add);
    return;
  }

  Point new_1(to_add[0] + 1, to_add[1]);
  if (!point_is_present_in_vector(new_1, points_to_consider)) {
    points_to_consider.push_back(new_1);
  }

  Point new_2(to_add[0], to_add[1] + 1);
  if (!point_is_present_in_vector(to_add, points_to_consider)) {
    points_to_consider.push_back(new_2);
  }

  Point new_3(to_add[0] + 1, to_add[1] + 1);
  if (!point_is_present_in_vector(to_add, points_to_consider)) {
    points_to_consider.push_back(new_3);
  }
}

vector<Point> compute_points_to_consider(const Point &original_point) {
  auto r_tuple = r_map[w_ptr];
  auto r_low = get<0>(r_tuple);
  auto r_high = get<1>(r_tuple);

  vector<Point> points_to_consider;
  points_to_consider.push_back(original_point);
  // points_to_consider.push_back(Point(original_point[0] + 1, original_point[1]));
  // points_to_consider.push_back(Point(original_point[0], original_point[1] + 1));
  // points_to_consider.push_back(Point(original_point[0] + 1, original_point[1] + 1));

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
array<array<double, DIM>, DIM> interpolate(const BoxData<double> G_i[DIM][DIM],
                                           const Particle &x_k,
                                           const double &h_g,
                                           const double& h_p)
{
  // cout << "BEGIN INTERPOLATION FUNC" << endl;
  array<array<double, DIM>, DIM> G_k;
  for (int i=0; i<DIM; ++i)
    for (int j=0; j<DIM;++j)
      G_k[i][j] = 0;
  // multidimensional interpolation algorithm
  Point left_hand_corner(floor(x_k.m_alpha[0] / h_g), floor(x_k.m_alpha[1] / h_g));
  auto grid_points = compute_points_to_consider(left_hand_corner);
  for (auto i_pos : grid_points) {
    Point x_bar(i_pos[0] * h_g, i_pos[1] * h_g);
    Point z(x_k.m_alpha[0] - x_bar[0], x_k.m_alpha[1] - x_bar[1]);
    // w_ptr = W_2 function
    double W_value = W(z, h_g, w_ptr);

    // Multiply G_i by W_value and incrememnt G_k matrix by G_i
    for (int i = 0; i < DIM; ++i) {
      for (int j = 0; j < DIM; ++j) {
        double G_i_val = G_i[i][j](i_pos);
        G_k[i][j] += G_i_val * W_value;
      }
    }
  }

  return G_k;
}

/**
 * BEGIN ARRAY SECTION
 */


bool point_is_present_in_vector(const array<double, DIM>& point_to_find,
                                const vector<array<double, DIM>>& vec) {
  if (find(vec.begin(), vec.end(), point_to_find) != vec.end()) {
    return true;
  }
  return false;
}

void complete_grid(vector<array<double, DIM>>& points_to_consider, const array<double, DIM>& to_add) {
  if (!point_is_present_in_vector(to_add, points_to_consider)) {
    // if the point is present, then its points to consider must already be within
    points_to_consider.push_back(to_add);
    return;
  }

  array<double, DIM> new_1 = {to_add[0] + 1, to_add[1]};
  if (!point_is_present_in_vector(new_1, points_to_consider)) {
    points_to_consider.push_back(new_1);
  }

  array<double, DIM> new_2 = {to_add[0], to_add[1] + 1};
  if (!point_is_present_in_vector(to_add, points_to_consider)) {
    points_to_consider.push_back(new_2);
  }

  array<double, DIM> new_3 = {to_add[0] + 1, to_add[1]};
  if (!point_is_present_in_vector(to_add, points_to_consider)) {
    points_to_consider.push_back(new_3);
  }
}

vector<array<double, DIM>> compute_points_to_consider(const array<double, DIM> &original_point) {
  auto r_tuple = r_map[w_ptr];
  auto r_low = get<0>(r_tuple);
  auto r_high = get<1>(r_tuple);

  vector<array<double, DIM>> points_to_consider;
  points_to_consider.push_back(original_point);
  array<double, DIM> original_2 = {original_point[0] + 1, original_point[1]};
  points_to_consider.push_back(original_2);
  array<double, DIM> original_3 = {original_point[0], original_point[1] + 1};
  points_to_consider.push_back(original_3);
  array<double, DIM> original_4 = {original_point[0] + 1, original_point[1] + 1};
  points_to_consider.push_back(original_4);

  for (int i = r_low; i <= r_high; ++i) {
    array<double, DIM> new_1 = {original_point[0] + i, original_point[1]};
    complete_grid(points_to_consider, new_1);
    array<double, DIM> new_2 = {original_point[0], original_point[1] + i};
    complete_grid(points_to_consider, new_2);
    array<double, DIM> new_3 = {original_point[0] + i, original_point[1] + i};
    complete_grid(points_to_consider, new_3);
  }

  return points_to_consider;
}

// m_dx is h_g
array<array<double, DIM>, DIM> interpolate_array(const BoxData<double> G_i[DIM][DIM],
                                           const Particle &x_k,
                                           const double &h_g,
                                           const double& h_p)
{
  // cout << "BEGIN INTERPOLATION FUNC" << endl;
  array<array<double, DIM>, DIM> G_k;
  for (int i=0; i<DIM; ++i)
  {
    for (int j=0; j<DIM;++j)
    {
      G_k[i][j] = 0;
    }
  }
  // multidimensional interpolation algorithm
  array<double, DIM> left_hand_corner = {floor(x_k.m_alpha[0] / h_g), floor(x_k.m_alpha[1] / h_g)};
  // cout << "lhc: " << left_hand_corner[0] << "," << left_hand_corner[1] << endl;
  auto grid_points = compute_points_to_consider(left_hand_corner);
  for (auto i_pos : grid_points) {
    // cout << "i_pos: " << i_pos[0] << "," << i_pos[1] << endl;
    array<double, DIM> x_bar = {i_pos[0] * h_g, i_pos[1] * h_g};
    array<double, DIM> z = {x_k.m_alpha[0] - x_bar[0], x_k.m_alpha[1] - x_bar[1]};
    // cout << "z: " << z[0] << "," << z[1] << endl;

    // w_ptr = W_2 function
    double W_value = W(z, h_g, w_ptr);

    // Multiply G_i by W_value and incrememnt G_k matrix by G_i
    for (int i = 0; i < DIM; ++i) {
      for (int j = 0; j < DIM; ++j) {
        Point point_accessor(i_pos[0], i_pos[1]);
        double G_i_val = G_i[i][j](point_accessor);
        G_k[i][j] += G_i_val * W_value;
      }
    }

  }

  return G_k;
}