/* this is the header file for montecarlo integration routine vegas */

#pragma once
#include <default.h>
#include <nrutil.h>
#include <ran2.h>

typedef struct {
  double *regn;
  int ndim;
  unsigned long ncall;
  int init, itmx, nprn;
  void *addPar;
} vegas_par;

void vegas(double (*fxn)(double[], double wgt, void *), double *tgral,
           double *sd, double *chi2a, vegas_par *vP);

void rebin(double rc, int nd, double r[], double xin[], double xi[]);
