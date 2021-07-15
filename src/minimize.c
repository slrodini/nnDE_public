#include <minimize.h>

static double stepsize = 0.1;
static double lineeps = 0.05;

double run_adam(double *par, int nPar, void *addPar,
                void (*fnc)(void *, double *, double *)) {
  adamPar aP;
  aP.alpha = 0.1;
  aP.beta1 = 0.9;
  aP.beta2 = 0.92;
  aP.maxIt = 1e+3;
  aP.dfToll = 1e-1;
  aP.changeLearn = 400;
  double chiMin = adam2(par, nPar, addPar, fnc, &aP);

  aP.alpha = 0.01;
  aP.beta1 = 0.9;
  aP.beta2 = 0.99;
  aP.maxIt = 1e+3;
  aP.dfToll = 1e-2;
  aP.changeLearn = 400;
  chiMin = adam2(par, nPar, addPar, fnc, &aP);

  aP.alpha = 0.001;
  aP.beta1 = 0.9;
  aP.beta2 = 0.999;
  aP.maxIt = 1e+4;
  aP.changeLearn = 5000;

  aP.dfToll = 1e-5;
  return adam2(par, nPar, addPar, fnc, &aP);
}

typedef struct {
  void *addPar;
  void (*fnc)(void *, double *, double *);
  int nPar;
  double *par;
  double *c2, *grad;

} gslMinimPar;

double my_f(const gsl_vector *v, void *params) {
  gslMinimPar *gP = (gslMinimPar *)params;
  for (int i = 0; i < gP->nPar; i++) {
    gP->par[i] = gsl_vector_get(v, i);
  }
  gP->fnc(gP->addPar, gP->c2, gP->grad);
  return *(gP->c2);
}

/* The gradient of f, df = (df/dx, df/dy). */
void my_df(const gsl_vector *v, void *params, gsl_vector *df) {
  gslMinimPar *gP = (gslMinimPar *)params;
  for (int i = 0; i < gP->nPar; i++) {
    gP->par[i] = gsl_vector_get(v, i);
  }
  gP->fnc(gP->addPar, gP->c2, gP->grad);
  for (int i = 0; i < gP->nPar; i++)
    gsl_vector_set(df, i, gP->grad[i]);
}

/* Compute both f and df together. */
void my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df) {
  *f = my_f(x, params);
  my_df(x, params, df);
}

double run_bfgs2(double *par, int nPar, void *addPar,
                 void (*fnc)(void *, double *, double *)) {
  int iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  /* Position of the minimum (1,2), scale factors
     10,20, height 30. */
  gsl_vector *x, *x_grad;
  gsl_multimin_function_fdf my_func;

  double grad[nPar];
  double c2 = 0.0;

  gslMinimPar gP;
  gP.addPar = addPar;
  gP.nPar = nPar;
  gP.par = par;
  gP.grad = grad;
  gP.c2 = &c2;
  gP.fnc = fnc;

  my_func.n = nPar;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = (void *)&gP;

  x = gsl_vector_alloc(nPar);
  x_grad = gsl_vector_alloc(nPar);
  for (int i = 0; i < nPar; i++)
    gsl_vector_set(x, i, par[i]);
  // x->data = par;

  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc(T, nPar);

  gsl_multimin_fdfminimizer_set(s, &my_func, x, stepsize, lineeps);

  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);
    printf("bfgs2 Iter: %4d  cur_chi2: %.4e\n", iter, s->f);

    if (status) {
      break;
    }

    // status = gsl_multimin_test_gradient(s->gradient, 1e-4);
    /*
    x_grad = gsl_multimin_fdfminimizer_gradient(s);
    double sum = 0.0;
    for (int i = 0; i < nPar; i++) {
      sum += pow(gsl_vector_get(x_grad, i), 2);
    }
    sum = sqrt(sum);

    if (sum > 1e-8) {
      status = GSL_CONTINUE;
    } else {
      status = GSL_SUCCESS;
    }
    */

  } while (status == GSL_CONTINUE && iter < 1e+3);
  // printf("Exit status: %d\n", status);
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return c2;
}

double minimize(double *par, int nPar, void *addPar,
                void (*fnc)(void *, double *, double *)) {

  double adamRes = run_adam(par, nPar, addPar, fnc);
  // double gslRes = run_bfgs2(par, nPar, addPar, fnc);
  // stepsize = 0.001;
  // lineeps = 0;
  // gslRes = run_bfgs2(par, nPar, addPar, fnc);

  // printf("Chi2 after adam2: %e\n", adamRes);
  // printf("Chi2 after bfgs2: %e\n", gslRes);
  return adamRes;
}

double minimize2(double *par, int nPar, void *addPar,
                 void (*fnc)(void *, double *, double *)) {
  adamPar aP;
  aP.alpha = 0.001;
  aP.beta1 = 0.9;
  aP.beta2 = 0.999;
  aP.maxIt = 1e+4;
  aP.changeLearn = 3000;

  aP.dfToll = 0.5 * 1e-4;
  return adam2(par, nPar, addPar, fnc, &aP);
  double gslRes = run_bfgs2(par, nPar, addPar, fnc);
  return gslRes;
}
