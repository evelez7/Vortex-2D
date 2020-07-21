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
typedef array<array<double,DIM> , DIM> matrix;
typedef array<double,DIM> vec;
matrix transpose(const matrix& a_m)
{
  matrix retval;
  for (int i = 0; i < DIM; i++)
    {
      for (int j = 0; j < DIM; j++)
        {
          retval[i][j] = a_m[j][i];
        }
    }
  return retval; 
}
vec prod(const matrix& a_m,const vec& a_v)
{
  vec retval;
  for (int i = 0; i < DIM; i++)
    {
      retval[i] = 0;
      for (int j = 0; j < DIM; j++)
        {
          retval[i] += a_v[j]*a_m[i][j];
        }
    }
  return retval;
}
matrix prod(const matrix& a_m1,const matrix& a_m2)
{
  matrix retval;
  for (int i = 0; i < DIM; i++)
    {
      for (int j = 0; j < DIM; j++)
        {
          retval[i][j] = 0;
          for (int k = 0; k < DIM; k++)
            {
              retval[i][j] += a_m1[i][k]*a_m2[k][j];
            }
        }
    }
  return retval;
}
matrix rotMat(const double& t)
{
  matrix retval;
  retval[0][0] = cos(M_PI*t);
  retval[1][1] = cos(M_PI*t);
  retval[0][1] = -sin(M_PI*t);
  retval[1][0] = sin(M_PI*t);
  return retval;
}
matrix strainMat(const double& t)
{
  matrix retval;
  retval[0][0] = 1/t;
  retval[1][1] = t;
  retval[0][1] = 0;
  retval[1][0] = 0;
  return retval;
}
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
void outVort(ParticleSet& p, int a_coarsenFactor,double a_angle)
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
  string filename = string("rotate") + to_string(a_angle);
  WriteData(outVort, 0, h, filename,filename);
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
  cout << "enter angle, alpha" << endl;
  double angle,alpha;
  cin >> angle;
  cin >> alpha;
  ParticleSet p;
  N = Power(2,M);
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
      vec xp;
      matrix rot,sym,mat;
      rot = rotMat(angle);
      sym = strainMat(alpha);
      mat = prod(rot,prod(prod(rot,sym),transpose(rot)));
      for (int i = 1;i < Np/4;i++)
        {
          xp[0] = i*hp;// + sqrt(3)*hp/2;
          for (int j = 1; j< Np/4; j++)
            {
              xp[1] = j*hp;// + sqrt(2)*hp/2;
              vec dist;
              dist[0] = xp[0] - .125;
              dist[1] = xp[1] - .125;
              Particle part;
              part.m_x = prod(mat,dist);
              part.m_x[0] += .5;
              part.m_x[1] += .5;
              //part.m_x[0] = mat[0][0]*dist0 + mat[0][1]*dist1 + .5;
              //part.m_x[1] = mat[1][0]*dist0 + mat[1][1]*dist1 + .5;
              part.strength = hp*hp/h/h;
              p.m_particles.push_back(part);
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
  outVort(p,pcfactor,angle);
}
