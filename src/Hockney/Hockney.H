#ifndef _HOCKNEY_
#define _HOCKNEY_
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
#include "FFTMD.H"
#include "ConvKernel.H"
#include <array>
#include <memory>
using namespace std;
using namespace Proto;
class Hockney
{
public:
  Hockney();
  Hockney(shared_ptr<ConvKernel>& a_kerPtr,const double& a_h,int a_M);
  void define(shared_ptr<ConvKernel>& a_kerPtr,const double& a_h,int a_M);
  void convolve(BoxData<double>& a_rhs);
  ~Hockney(){};
protected:
  double m_h;
  int m_M,m_N;
  FFTMD m_fftmd;
  shared_ptr<ConvKernel> m_kerPtr;
  bool m_isDefined;
};
#endif
  
