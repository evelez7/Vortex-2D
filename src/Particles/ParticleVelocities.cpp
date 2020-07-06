#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Hockney.H"
#include "ParticleSet.H"
#include "ParticleVelocities.H"
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_WriteBoxData.H"
//#include "VisitWriter.H"
#include "Proto_Timer.H"
using namespace std;
using namespace Proto;

ParticleVelocities:: ParticleVelocities(){};
void ParticleVelocities::operator()
  (ParticleShift& a_k,const double& a_time,const double& a_dt,
   ParticleSet& a_state)
{
  PR_TIMERS("ParticleVelocities");
  PR_TIMER("getRHS",t1);
  PR_TIMER("Hockney",t2);
  PR_TIMER("Fields",t3);
  BoxData<double> rhs(a_state.m_box);
  Point ipos;
  double xpos[DIM];
  double weight;
  Point e0 = Point::Basis(0);
  Point e1 = Point::Basis(1);
  rhs.setVal(0.);
  PR_START(t1);
  
  int N = a_state.m_box.high()[1];
  
  const vector<Particle >& oldPart = a_state.m_particles;
  vector<DX >& dPart = a_k.m_particles;
  for (int k = 0; k < a_state.m_particles.size(); k++)
    {
      for (int l = 0; l < DIM; l++)
        {
          double newpos = oldPart[k].m_x[l] + dPart[k].m_x[l];
          ipos[l] = newpos/a_state.m_dx;
          xpos[l] = (newpos - ipos[l]*a_state.m_dx)/a_state.m_dx;
        }
      //Point pt(ipos);
      int kind = ipos[0] + ipos[1]*(N + 1);
      for (int l0=0; l0 < DIM;l0++)
        {
          for (int l1=0;l1 < DIM ; l1++)
            {
              rhs(ipos + e0*l0 + e1*l1) += (1.-xpos[0] + (2*xpos[0] - 1.)*l0)*
                (1.- xpos[1] + (2*xpos[1] - 1.)*l1)*oldPart[k].strength;
            }
        }
      
    }
  PR_STOP(t1);
  PR_START(t2);
  a_state.m_hockney.convolve(rhs);
  PR_STOP(t2);
  PR_START(t3);
  Box gradBox = a_state.m_box.grow(-1);
  BoxData<double > field[DIM];
  for (int l=0 ; l < DIM; l++)
    {
      field[l].define(gradBox);
      Point evec = Point::Basis(l);
      for (auto it = gradBox.begin();!it.done();++it)
        {
          Point pt = *it;
          field[l](pt) = (rhs(pt + evec) - rhs(pt - evec))/(2*a_state.m_dx);
        }
    }
  
  for (int k = 0; k < a_state.m_particles.size(); k++)
    {
      for (int l = 0; l < DIM; l++)
        {
          double newpos = a_state.m_particles[k].m_x[l] + a_k.m_particles[k].m_x[l];
          ipos[l] = newpos/a_state.m_dx;
          xpos[l] = (newpos - ipos[l]*a_state.m_dx)/a_state.m_dx;
        }
      Point pt(ipos);
      assert(a_state.m_box.contains(pt));
      double vel[DIM];
      for (int dir = 0 ; dir < DIM ; dir++)
        {
          int dirperp = (dir + 1)%DIM;
          int sdir = 1-2*dir;
          double field00 = field[dirperp](pt);
          double field10 = field[dirperp](pt + e0);
          double field01 = field[dirperp](pt + e1);
          double field11 = field[dirperp](pt + e0 + e1);
         
          vel[dir] =
            (field00*(1.-xpos[0])*(1. - xpos[1]) +
             field10*(xpos[0])*(1.-xpos[1]) +
             field01*(xpos[1])*(1.-xpos[0]) +
             field11*xpos[0]*xpos[1])*sdir;
          dPart[k].m_x[dir] = vel[dir]*a_dt;
        }
      double velmag = sqrt(vel[0]*vel[0] + vel[1]*vel[1]);
      // cout << "velmag at " << k << " = " << velmag <<endl; 
      //cout << "vel at particle " << k << " = " << vel[0] << " , " << vel[1] << endl;
    }
  PR_STOP(t3);
};
