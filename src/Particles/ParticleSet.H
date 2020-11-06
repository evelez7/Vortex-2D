#ifndef _PARTICLESET_H_
#define _PARTICLESET_H_

#include <vector>
#include "Particle.H"
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_WriteBoxData.H"
#include "Hockney.H"

using namespace std;
using namespace Proto;

class ParticleSet;
class ParticleShift
{
public:
  vector<DX> m_particles;
  ParticleShift(){};
  ~ParticleShift(){};
  void init(const ParticleSet& a_particles);
  void increment(double a_scale,
                 const ParticleShift& a_rhs);
  void operator*=(double a_scale);
  void setToZero();
};
class ParticleSet
{
public:
  ParticleSet(
              shared_ptr<ConvKernel>& a_kerptr,
              Box& a_box,
              double& a_dx,
              double& a_hp,
              array<double, DIM>& a_lowCorner,
              int a_M);

  ParticleSet(){};
  ~ParticleSet(){};
  vector<Particle> m_particles;
  double m_dx;
  double m_hp;
  Box m_box;
  int RK_count = 0;
  int file_count = 0;
  array<double,DIM> m_lowCorner;
  Hockney m_hockney;
  void increment(const ParticleShift& a_shift);

};


#endif
