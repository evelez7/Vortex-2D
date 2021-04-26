#ifndef _PARTICLEWRITE_H_
#define _PARTICLEWRITE_H_
#include "ParticleSet.H"
#include <cstring>

void PWrite(const char* filename, const ParticleSet* a_particles, int fileCount);

const char* PWrite(const ParticleSet* a_particles);

inline void write_point_mesh(const char *filename, int npts,
                      double *pts, int nvars, int *vardim,
                      const char * const *varnames, double **vars);

#endif
