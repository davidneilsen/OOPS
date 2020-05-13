#ifndef OUTPUT_H
#define OUTPUT_H
#include <oopsconfig.h>
#include <cstdio>

#ifdef USE_SDF
#include "sdf.h"
#endif

namespace output{
  void output_data(char *name, double *v, double *r, int size, double time);
  void output_vtk(char *name, double *v, const double *r, int size, double time, int iter);
};

#endif
