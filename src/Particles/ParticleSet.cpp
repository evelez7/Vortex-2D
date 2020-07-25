#include "ParticleSet.H"
#include <cassert>
ParticleSet::ParticleSet(shared_ptr<ConvKernel>& a_kerptr,Box& a_box,double& a_dx, double& a_hp, array<double, DIM>& a_lowCorner,int a_M)
{
  m_box = a_box;
  m_dx = a_dx;
  m_hp = a_hp;
  for (int l = 0; l < DIM; l++)
    {
      m_lowCorner[l] = a_lowCorner[l];
    }
  m_hockney.define(a_kerptr,m_dx,a_M);
};
void ParticleShift::increment(double a_scale, const ParticleShift& a_rhs)
{
  assert(m_particles.size() == a_rhs.m_particles.size());
  for(unsigned int i=0; i<m_particles.size(); i++)
    {
      m_particles[i].increment(a_scale, a_rhs.m_particles[i]);
    }
}

void ParticleShift::operator*=(double a_scale)
{
  for(unsigned int i=0; i<m_particles.size(); i++)
    {
      m_particles[i]*=a_scale;
    }
}

void ParticleShift::setToZero()
{
  for(unsigned int i=0; i<m_particles.size(); i++)
    {
      m_particles[i]*=0;
    }
}

void ParticleShift::init(const ParticleSet& a_rhs)
{
  m_particles.resize(0);
  m_particles.resize(a_rhs.m_particles.size());
}

void ParticleSet::increment(const ParticleShift& a_rhs)
{
  assert(m_particles.size() == a_rhs.m_particles.size());
  for(unsigned int i=0; i<m_particles.size(); i++)
    {
      m_particles[i].increment(a_rhs.m_particles[i]);
    }
}



