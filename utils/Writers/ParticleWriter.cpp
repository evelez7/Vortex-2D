#include "ParticleWriter.H"
#include "Proto_VisitWriter.H"
#include "ParticleSet.H"

#include <cstdio>
#include <cstdlib>
#include <cstring>
static FILE *fp = NULL;
static int useBinary = 0;
static int numInColumn = 0;

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
  vector<vector<double> > vars(3);
  unsigned int size = a_p->m_particles.size();
  std::vector<double> x(3*size);
  vars[0] = std::vector<double>(size);
  vars[1] = std::vector<double>(size);
  vars[2] = std::vector<double>(size);
  for(unsigned int i=0; i<size; i++)
    {
      const Particle& p = a_p->m_particles[i];
      vars[0][i] = p.strength;
      vars[1][i] = p.m_alpha[0];
      vars[2][i] = p.m_alpha[1];
      x[i*3] = p.m_x[0];
      x[i*3+1] = p.m_x[1];
#if DIM==3
      x[i*3+2] = p.m_x[2];
#else
      x[i*3+2] = 0.0;
#endif
    }
  double* varPtr[3];
  varPtr[0] = &vars[0][0];
  varPtr[1] = &vars[1][0];
  varPtr[2] = &vars[2][0];
  int vardim[3] = {1,1,1};
  const char* const varnames[] = {"strength","alpha1","alpha2"}; 

  write_point_mesh(a_filename, size,
		   &(x[0]), 3, vardim,
                   varnames, varPtr);
}
inline void write_point_mesh(const char* filename, int npts, double *pts,
                      int nvars, int *vardim, const char * const *varnames,
                      double **vars)
{
    FILE* fp = vtk_open_file(filename);
    int numInColumn = 0;
     
    int   i;
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
