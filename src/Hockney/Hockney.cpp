#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include "PowerItoI.H"
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_WriteBoxData.H"
//#include "VisitWriter.H"
#include "Proto_Timer.H"
#include "FFT1DW.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "Hockney.H"
#include "ConvKernel.H"

using namespace std;
using namespace Proto;
void copyReal(BoxData<complex<double> >& a_cxarray,BoxData<double >& a_real)
{
  Box d = a_real.box()&a_cxarray.box();
  for (auto it = d.begin();!it.done();++it)
    {
      Point pt = *it;
      a_real(pt) = real(a_cxarray(pt));
    }
};
void copyImag(BoxData<complex<double> >& a_cxarray,BoxData<double >& a_imag)
{
  Box d = a_imag.box();
  for (auto it = d.begin();!it.done();++it)
    {
      Point pt = *it;
      a_imag(pt) = imag(a_cxarray(pt));
    }
};
Hockney::Hockney()
{
  m_isDefined = false;
};
Hockney::Hockney(shared_ptr<ConvKernel>& a_kerPtr,const double& a_h,int a_M)
{ 
  shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
  shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
  m_fftmd.define(p_fft);
  m_kerPtr = a_kerPtr;
  m_h = a_h;
  m_M = a_M;
  m_N = Power(2,a_M);
};
void Hockney::define(shared_ptr<ConvKernel>& a_kerPtr,const double& a_h,int a_M)
{ 
  shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
  shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
  m_fftmd.define(p_fft);
  m_kerPtr = a_kerPtr;
  m_h = a_h;
  m_M = a_M;
  m_N = Power(2,a_M);
};
void Hockney::convolve(BoxData<double>& a_rhs)
{
  Box rhsDomain = a_rhs.box();
  Point low = rhsDomain.low();
  Point high = rhsDomain.high();

  assert(low == Point::Zeros());
  assert(high == Point::Ones()*m_N);

  low = high*(-1);
  Box ddomain(low,high);
  BoxData<complex<double> > rhsDouble(ddomain);
  complex<double> zero(0.,0.);
  rhsDouble.setVal(zero);
  double scale = 1./pow(m_N*1.,DIM*2)/4;
  for (auto it = rhsDomain.begin();!it.done();++it)
    {
      Point pt = *it;
      rhsDouble(pt).real(a_rhs(pt));
    }
  BoxData<complex<double> > kernel(ddomain);
  BoxData<double > realOut(ddomain);
  m_kerPtr->getKernel(kernel,m_h);
  m_fftmd.forwardCCcen(rhsDouble);
  m_fftmd.forwardCCcen(kernel);
  for (auto it = ddomain.begin();!it.done();++it)
    {
      Point pt = *it;
      rhsDouble(pt) *= kernel(pt);
    }
  m_fftmd.inverseCCcen(rhsDouble);
  a_rhs.setVal(0.);
  Box bx(rhsDomain.low(),rhsDomain.high() - Point::Ones());
   for (auto it = bx.begin();!it.done();++it)
     {
       Point pt = *it;
       a_rhs(pt) = real(rhsDouble(pt))*scale;
    }
}


