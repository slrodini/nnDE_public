#include <minimize.h>

double run_adam(double *par, int nPar, void *addPar,
                void (*fnc)(void *, double *, double *))
{
  // Set up multiple pass of ADAM (using ADAM2 in this case)
  // To speed up the first iterations and then refine the result
  adamPar aP;
  aP.alpha = 0.1;
  aP.beta1 = 0.9;
  aP.beta2 = 0.92;
  aP.maxIt = 1e+3;
  aP.dfToll = 1e-6;
  aP.changeLearn = 400;
  double chiMin = adam2(par, nPar, addPar, fnc, &aP);

  aP.alpha = 0.01;
  aP.beta1 = 0.9;
  aP.beta2 = 0.99;
  aP.maxIt = 1e+3;
  aP.dfToll = 1e-6;
  aP.changeLearn = 400;
  chiMin = adam2(par, nPar, addPar, fnc, &aP);

  aP.alpha = 0.001;
  aP.beta1 = 0.9;
  aP.beta2 = 0.999;
  aP.maxIt = 1e+4;
  aP.changeLearn = 5000;

  aP.dfToll = 1e-6;
  return adam2(par, nPar, addPar, fnc, &aP);
}

double minimize(double *par, int nPar, void *addPar,
                void (*fnc)(void *, double *, double *))
{

  double adamRes = run_adam(par, nPar, addPar, fnc);

  return adamRes;
}

double minimize2(double *par, int nPar, void *addPar,
                 void (*fnc)(void *, double *, double *))
{
  adamPar aP;
  aP.alpha = 0.001;
  aP.beta1 = 0.9;
  aP.beta2 = 0.999;
  aP.maxIt = 1e+4;
  aP.changeLearn = 3000;

  aP.dfToll = 0.5 * 1e-4;
  return adam2(par, nPar, addPar, fnc, &aP);
}
