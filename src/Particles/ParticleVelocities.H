#ifndef _PARTICLEVELOCITIES_H_
#define _PARTICLEVELOCITIES_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_WriteBoxData.H"
//#include "VisitWriter.H"
#include "Proto_Timer.H"
#include "Hockney.H"
#include "ParticleSet.H"
using namespace std;
using namespace Proto;
class ParticleVelocities
{
public:
  ParticleVelocities();
  void operator()(ParticleShift& a_k, 
                     const double& a_time, const double& dt, 
                     ParticleSet& a_state);
};
#endif
