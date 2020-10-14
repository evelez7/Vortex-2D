#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <string>
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_WriteBoxData.H"
#include "ParticleSet.H"
#include "ParticleVelocities.H"
#include "CutoffKernel.H"
//#include "VisitWriter.H"
#include "ParticleWriter.H"
#include "Proto_RK4.H"
using namespace std;
void outField(ParticleSet& a_state)
{
  BoxData<double> field(a_state.m_box);
  Point ipos;
  double xpos[DIM];
  double weight;
  Point e0 = Point::Basis(0);
  Point e1 = Point::Basis(1);
  field.setVal(0.);

  int N = a_state.m_box.high()[1];

  const vector<Particle >& oldPart = a_state.m_particles;
  for (int k = 0; k < a_state.m_particles.size(); k++)
    {
      for (int l = 0; l < DIM; l++)
        {
          double newpos = oldPart[k].m_x[l];
          ipos[l] = newpos/a_state.m_dx;
          xpos[l] = (newpos - ipos[l]*a_state.m_dx)/a_state.m_dx;
        }
      int kind = ipos[0] + ipos[1]*(N + 1);
      for (int l0=0; l0 < DIM;l0++)
        {
          for (int l1=0;l1 < DIM ; l1++)
            {
              field(ipos + e0*l0 + e1*l1) += (1.-xpos[0] + (2*xpos[0] - 1.)*l0)*
                (1.- xpos[1] + (2*xpos[1] - 1.)*l1)*oldPart[k].strength;
            }
        }

    }
  a_state.m_hockney.convolve(field);
  WriteBoxData("field",field);

}
void outVort(ParticleSet& p, int a_coarsenFactor,unsigned int a_nstep)
{
  int coarsenFactor = a_coarsenFactor;
  Box bx = p.m_box.coarsen(coarsenFactor);
  double h = p.m_dx*coarsenFactor;
  BoxData<double> outVort(bx);
  int ipos[DIM];
  array<double,DIM> xpos;
  double weight;
  Point e0 = Point::Basis(0);
  Point e1 = Point::Basis(1);
  //
  outVort.setVal(0.);
  for (int k = 0; k < p.m_particles.size(); k++)
    {
      for (int l = 0; l < DIM; l++)
        {
          double newpos = p.m_particles[k].m_x[l];
          newpos -= .5*h;
          ipos[l] = newpos/h;
          xpos[l] = (newpos - ipos[l]*h)/h;
        }
      Point pt(ipos);
      assert(p.m_box.contains(pt));
      for (int l0=0; l0 < DIM;l0++)
        {
          for (int l1=0;l1 < DIM ; l1++)
            {
              outVort(pt+e0*l0 + e1*l1) +=
                (1.-xpos[0] + (2*xpos[0] - 1.)*l0)*
                (1.- xpos[1] + (2*xpos[1] - 1.)*l1)*p.m_particles[k].strength/coarsenFactor/coarsenFactor;
            }
        }
    }
  string filename = string("vorticity") ;
  WriteData(outVort,   a_nstep, h, filename,filename);
};
int main(int argc, char* argv[])
{
  unsigned int M;
  unsigned int N;
  cout << "input log_2(number of grid points)" << endl;
  cin >> M;
  int test = 4;
  cout << "input particle refinement factor" << endl;
  unsigned int cfactor;
  cin >> cfactor;
  cout << "enter stopping time" << endl;
  double timeStop;
  cin >> timeStop;
  ParticleSet p;
  N = Power(2,M);
  PR_TIMER_SETFILE(to_string(N) + "proto.time.table");
  PR_TIMERS("main");
  double h = 1./N;
  double hp = h/cfactor; //pow(h,4./3.);
  int Np = 1./hp;
  hp = 1./Np;
  double delta = h;
  int pcfactor = 4/cfactor;
  if (pcfactor < 1 ) pcfactor = 1;
  cout << "number of particles per cell = " << h*h/hp/hp << endl;
  shared_ptr<CutoffKernel> cutkptr =
    shared_ptr<CutoffKernel>(new CutoffKernel(h,delta));
  shared_ptr<ConvKernel> convkerptr = dynamic_pointer_cast<ConvKernel>(cutkptr);
  p.m_hockney.define(convkerptr,h,M);
  p.m_dx = h;
  p.m_hp = hp;
  array<double,DIM> lowCorner;
  if (test == 1)
  {
    p.m_particles.resize(1);
      Particle& particle = p.m_particles[0];
      particle.m_x[0] = .49;
      particle.m_x[1] = .24;
      particle.strength = 1./h/h;
  }
  else if (test == 2)
    {
      p.m_particles.resize(2);
      Particle& particle = p.m_particles[0];
      particle.m_x[0] = .5;
      particle.m_x[1] = .25;
      particle.strength = 0.;
      Particle& particle2 = p.m_particles[1];
      particle2.m_x[0] = .5;
      particle2.m_x[1] = .5;
      particle2.strength = 1./h/h;
    }
else if (test == 3)
    {
      p.m_particles.resize(2);
      Particle& particle = p.m_particles[0];
      particle.m_x[0] = .5;
      particle.m_x[1] = .25;
      particle.strength = 1./h/h;
      Particle& particle2 = p.m_particles[1];
      particle2.m_x[0] = .5;
      particle2.m_x[1] = .75;
      particle2.strength = 1./h/h;
    }
  else
    {
      array<double,DIM> xp;
      for (int i = 0;i < Np;i++)
        {
          xp[0] = i*hp;
          for (int j = 0; j< Np; j++)
            {
              xp[1] = j*hp;
              double dist1 = sqrt(pow(xp[0] - .375,2) + pow(xp[1] - .5,2));
              double dist2 = sqrt(pow(xp[0] - .625,2) + pow(xp[1] - .5,2));
              if ((dist1 < .12 ) | (dist2 < .12))
            // if (dist1 < .1125 )
                {
                  Particle part;
                  part.m_x[0] = xp[0];
                  part.m_x[1] = xp[1];
                  part.strength = hp*hp/h/h;
                  part.m_alpha[0] = xp[0];
                  part.m_alpha[1] = xp[1];
                  p.m_particles.push_back(part);
                }
            }
        }
    }
  double dx = 1./N;
  cout << "number of particles = " << p.m_particles.size() << endl;
  Point low = Point::Zeros();
  Point high = Point::Unit()*N;
  Box bx(low,high);
  p.m_box=bx;
  lowCorner[0] = 0.;
  lowCorner[1] = 0.;
  ParticleShift kIn,kOut;
  kIn.init(p);
  kOut.init(p);
  kIn.setToZero();
  ParticleVelocities pv;
  double time = 0.;
  double dt = 140*.025/N;
  // double dt = 0;
  int m = 5000;
  // int m = 1;

  RK4<ParticleSet,ParticleVelocities,ParticleShift> integrator;
#if ANIMATION
  outVort(p,pcfactor,0);
  PWrite(&p);
#endif
  for(int i=0; i<m; i++)
    {
      integrator.advance(time, dt, p);
      time = time + dt;
#if ANIMATION
      outVort(p,pcfactor,i);
      PWrite(&p);
#endif
      if (time >= timeStop)
        {
          break;
        }
    }
  if (((test == 1) || (test == 2) || (test == 3)))
    {
      Particle part0 = p.m_particles[0];
      double radius = sqrt(pow(part0.m_x[0] - .5,2) + pow(part0.m_x[1] - .5,2));
      cout << "test = "<<test << ": radial position of first particle = " <<
        radius << endl;
      cout << "Cartesian position of first particle = " << part0.m_x[0] << " , " << part0.m_x[1] << endl;
    }
  outField(p);
  PR_TIMER_REPORT();
}
