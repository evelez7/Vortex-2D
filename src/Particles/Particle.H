#ifndef _PARTICLE_H_
#define _PARTICLE_H_
#include <cstdio>
#include <cmath>
#include <cassert>
#include <array>
#include <iostream>
using namespace std;
class Particle;
class DX
{
public:
  DX()
  {
    for(unsigned int i=0; i<DIM; i++)
      {
	m_x[i]=0;
      }
  }
  array<double, DIM> m_x;
  array<array<double, DIM >, DIM > m_gradx {
    {
      {{0, 0}},
      {{0, 0}}
    }
  };
  inline void increment(double a_scale, const DX& a_rhs)
  {
    for(unsigned int i=0; i<DIM; i++)
      {
	      m_x[i]+=a_scale*a_rhs.m_x[i];
        for (unsigned int j=0; j<DIM; ++j)
        {
          m_gradx[i][j] += a_scale*a_rhs.m_gradx[i][j];
        }
      }
  }
  inline void operator*=(double a_scale)
  {
    for(unsigned int i=0; i<DIM; i++)
      {
	      m_x[i]*=a_scale;
        for (unsigned int j=0; j<DIM; ++j)
        {
          m_gradx[i][j] *= a_scale;
        }
      }
  }

  inline void init(const Particle &a_state)
  {
    for (int i=0; i<DIM; ++i)
      for (int j=0; j<DIM; ++j)
        m_gradx[i][j]=0;
  }
};
class Particle
{
public:
  array<double, DIM> m_x;
  double strength;
  array<double, DIM> m_alpha;
  array<array<double, DIM>, DIM> m_gradx {
    {
      {{1, 0}},
      {{0, 1}}
    }
  };
  inline void increment(const DX& a_shift)
  {
    for(unsigned int i=0; i<DIM; i++)
      {
	      m_x[i]+=a_shift.m_x[i];
        for (unsigned int j= 0; j<DIM; ++j)
        {
          m_gradx[i][j] += a_shift.m_gradx[i][j];
        }
      }
  }
};

#endif
