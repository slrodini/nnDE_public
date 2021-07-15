#include <adam.h>
#include <gsl/gsl_multimin.h>

double minimize(double *par, int nPar, void *addPar,
                void (*fnc)(void *, double *, double *));

double minimize2(double *par, int nPar, void *addPar,
                 void (*fnc)(void *, double *, double *));
