//#include "nrutil.h"
#ifndef __SIMANN_H
#define __SIMANN_H
#include <ran2.h>
//#include "vegas.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  double Tin, Tfin, dt;
  int nIter;
} annealing_par;

void help();

double annealing(double *par, int nPar, double *parErr, void *addPar,
                 double (*fEn)(double *, void *), annealing_par *anPar,
                 int *ranMode);

// double chi2(double *par, void *extraPar);
#endif
