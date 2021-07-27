#include <default.h>

typedef struct
{
  double alpha, beta1, beta2, dfToll;
  int maxIt, changeLearn;
} adamPar;

double adamax(double *par, int nPar, void *addPar,
              void (*fgrad)(void *, double *, double *), adamPar *adPar);
double adam2(double *par, int nPar, void *addPar,
             void (*fgrad)(void *, double *, double *), adamPar *adPar);
