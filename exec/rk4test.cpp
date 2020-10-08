#include "Proto_RK4.H"
#include "Particle.H"
#include "SingleParticleVelocity.H"
#include "util.h"
#include <iostream>
#include <cmath>
#include <array>
using namespace std;


array<array<double, DIM>, DIM> real_solution_part(const double x)
{
  array<array<double, DIM>, DIM> real;
  real[0][0] = cos(x);
  real[0][1] = sin(x);
  real[1][0] = -sin(x);
  real[1][1] = cos(x);
  return real;
}

array<array<double, DIM>, DIM> real_solution_dx(const double x)
{
  array<array<double, DIM>, DIM> real;
  real[0][0] = -sin(x);
  real[0][1] = cos(x);
  real[1][0] = -cos(x);
  real[1][1] = -sin(x);
  return real;
}
array<array<double, DIM>, DIM> matrix_subtraction(array<array<double, DIM>, DIM> exact, array<array<double, DIM>, DIM> estimate)
{
  array<array<double, DIM>, DIM> error;
  for (int i=0; i<DIM; ++i) {

    for (int j=0; j<DIM;++j)
      error[i][j] = exact[i][j]-estimate[i][j];
  }
  return error;
}

double get_norm(array<array<double, DIM>, DIM> to_get)
{
  double max = -INFINITY;
  for (int i=0; i <DIM; ++i)
  {
    double sum = 0;
    for (int j=0; j<DIM; ++j)
      sum+= to_get[i][j];

    if (sum > max)
      max = sum;
  }
  return max;
}

double error_check(array<array<double, DIM>, DIM> exact, array<array<double, DIM>, DIM> estimate)
{
  auto difference = matrix_subtraction(estimate, exact);
  double difference_norm = get_norm(difference);
  auto exact_norm = get_norm(exact);
  auto estimate_norm=get_norm(estimate);
  return difference_norm/exact_norm;
}

int main(int argc, char** argv)
{
  double time = 0.25;
  array<double, 3> dt = {0.025, 0.0125, 0.00625};
  // array<double, 1> dt = {0};

  Proto::RK4<Particle, SingleParticleVelocity, DX> solver;
  DX a_shift;
  Particle a_particle;

  for (int i=0; i < dt.size(); ++i)
  {
    cout << "=======NEW ADVANCE at dt= " << dt[i] << "============" << endl;
    solver.advance(time, dt[i], a_particle);
    cout << "PARTICLE AFTER ADVANCE" << endl;
    print_matrix_here(a_particle.m_gradx);
    cout << "REAL SOLUTION PARTICLE" << endl;
    print_matrix_here(real_solution_part(dt[i]));
    cout << "REAL SOLUTION DX" << endl;
    print_matrix_here(real_solution_dx(dt[i]));
    cout << "PARTICLE ERROR" << endl;
    cout << error_check(real_solution_part(dt[i]), a_particle.m_gradx) << endl;
    // a_particle.m_gradx[0][0]=1;
    // a_particle.m_gradx[0][1]=0;
    // a_particle.m_gradx[1][0]=0;
    // a_particle.m_gradx[1][1]=1;
  }

  return 0;
}
