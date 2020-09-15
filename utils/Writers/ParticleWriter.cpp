#include "ParticleWriter.H"
#include "Proto_VisitWriter.H"
#include "ParticleSet.H"
#include "eigen_util.h"
#include <iostream>
#include <array>
#include <algorithm>
#include <tuple>
#include <memory>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
static FILE *fp = NULL;
static int useBinary = 0;
static int numInColumn = 0;

using namespace std;

const char* PWrite(const ParticleSet* a_array)
{
  static int fileCount = 0;
  static char nameBuffer[10];
  if(a_array == NULL)
    {
      return nameBuffer;
    }
  sprintf(nameBuffer, "PART.%d",fileCount);
  PWrite(nameBuffer, a_array);
  fileCount++;
  return nameBuffer;
}

void PWrite(const char* a_filename, const ParticleSet* a_p)
{
  if(a_filename == NULL || a_p == NULL)
    {
      return;
    }
  vector<vector<double> > vars(6);
  unsigned int size = a_p->m_particles.size();
  std::vector<double> x(3*size);
  vars[0] = std::vector<double>(size);
  vars[1] = std::vector<double>(size);
  vars[2] = std::vector<double>(size);
  vars[3] = std::vector<double>(size);
  vars[4] = std::vector<double>(size);
  vars[5] = std::vector<double>(size);
  for(unsigned int i=0; i<size; i++)
    {
      const Particle& p = a_p->m_particles[i];
      // equation 30
      auto A_t_A = multiply_matrices(get_transpose(p.m_gradx), p.m_gradx); // the symmetric and positive definite matrix
      //equation 30-32
      auto eigenvalues = get_sym_eigenvalues(A_t_A); // get the eigenvalues from the symmetric part of the polar decomp
      double eigen_product = eigenvalues[0] * eigenvalues[1];
      auto grad_det = get_determinant(p.m_gradx);

      double max_eigenvalue = -INFINITY;
      // find greatest eigenvalue from the decomp
      for (auto elem : eigenvalues)
        if (elem > max_eigenvalue)
          max_eigenvalue = elem;

      vars[0][i] = p.strength;
      vars[1][i] = p.m_alpha[0];
      vars[2][i] = p.m_alpha[1];
      vars[3][i] = max_eigenvalue;
      vars[4][i] = eigen_product;
      vars[5][i] = grad_det;
      x[i*3] = p.m_x[0];
      x[i*3+1] = p.m_x[1];
#if DIM==3
      x[i*3+2] = p.m_x[2];
#else
      x[i*3+2] = 0.0;
#endif
    }
  double* varPtr[6];
  varPtr[0] = &vars[0][0];
  varPtr[1] = &vars[1][0];
  varPtr[2] = &vars[2][0];
  varPtr[3] = &vars[3][0];
  varPtr[4] = &vars[4][0];
  varPtr[5] = &vars[5][0];
  int vardim[6] = {1,1,1,1,1, 1};
  const char* const varnames[] = {"strength","alpha1","alpha2","max_eigenvalue", "eigen_product", "grad_det"};

  write_point_mesh(a_filename, size,
		   &(x[0]), 6, vardim,
                   varnames, varPtr);
}
inline void write_point_mesh(const char* filename, int npts, double *pts,
                      int nvars, int *vardim, const char * const *varnames,
                      double **vars)
{
    FILE* fp = vtk_open_file(filename);
    int numInColumn = 0;

    int i;
    char  str[128];
    int  *centering = NULL;

    write_header(fp);

    write_string("DATASET UNSTRUCTURED_GRID\n",fp);
    sprintf(str, "POINTS %d double\n", npts);
    write_string(str,fp);
    for (i = 0 ; i < 3*npts ; i++)
    {
        write_double(pts[i],fp,numInColumn);
    }

    new_section(fp,numInColumn);
    sprintf(str, "CELLS %d %d\n", npts, 2*npts);
    write_string(str,fp);
    for (i = 0 ; i < npts ; i++)
    {
        write_int(1,fp,numInColumn);
        write_int(i,fp,numInColumn);
        end_line(fp,numInColumn);
    }
    centering = (int *) malloc(nvars*sizeof(int));
    for (i = 0 ; i < nvars ; i++)
        centering[i] = 1;
    write_variables(nvars, vardim, centering, varnames, vars, npts, npts, fp, numInColumn);
    free(centering);

    vtk_close_file(fp,numInColumn);
}
