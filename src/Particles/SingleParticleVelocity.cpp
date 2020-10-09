
#include <iostream>
#include <cassert>
#include <cmath>
#include <array>
#include "util.h"
#include "interpolate.h"
#include "SingleParticleVelocity.H"

using namespace std;

array<array<double, DIM>, DIM> subtraction_here(array<array<double, DIM>, DIM> exact, array<array<double, DIM>, DIM> estimate)
{
  array<array<double, DIM>, DIM> error;
  for (int i=0; i<DIM; ++i)
    for (int j=0; j<DIM;++j)
      error[i][j] = exact[i][j]-estimate[i][j];

  return error;
}
array<array<double, DIM>, DIM> example = {
  {
    {{0., 1.}},
    {{-1., 0.}}
  }
};
array<array<double, DIM>, DIM> real_solution(const double x)
{
  array<array<double, DIM>, DIM> real;
  real[0][0] = -sin(x);
  real[0][1] = cos(x);
  real[1][0] = -cos(x);
  real[1][1] = -sin(x);
  return real;
}
SingleParticleVelocity::SingleParticleVelocity(){};

void SingleParticleVelocity::operator()(DX &a_shift, double a_time, double a_dt, Particle &a_particle)
{
  for (int i=0; i<DIM; ++i)
  {
    for (int j=0; j<DIM; ++j)
    {
      a_shift.m_gradx[i][j]=0;
      for (int k=0; k<DIM; ++k)
        a_shift.m_gradx[i][j] += real_solution(a_time)[i][k] * a_particle.m_gradx[k][j];
        // a_shift.m_gradx[i][j] += real_solution(a_dt)[i][k] * a_particle.m_gradx[k][j];
      a_shift.m_gradx[i][j] *= a_dt;
    }
  }
}
