#ifndef _CONVKERNEL_H_
#define  _CONVKERNEL_H_

#include <assert.h>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_WriteBoxData.H"

#include "PowerItoI.H"
using namespace std;
using namespace Proto;
class ConvKernel
{
public:
  /// Default constructor.
  ConvKernel(){};
  /// Destructor.
  virtual ~ConvKernel(){};
  /// Defines a default-constructed RectMDArray.  
  virtual void getKernel(BoxData<complex<double> >& a_srcArray,double& a_h) = 0; 
private:
};  
#endif
