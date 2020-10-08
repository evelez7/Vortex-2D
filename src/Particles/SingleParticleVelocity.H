#include "Particle.H"

class SingleParticleVelocity
{
public:
  SingleParticleVelocity();
  void operator()(DX &, double a_time, double a_dt, Particle &);
};