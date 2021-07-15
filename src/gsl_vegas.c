#include <gsl_vegas.h>

double vegas(double (*futil)(double *, size_t, void *), double *xl, double *xu,
             int ndim, void *addPar) {
  double res, err;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = {futil, ndim, addPar};

  size_t calls = 5000;

  gsl_rng_env_setup();

  T = gsl_rng_ranlxd2;
  r = gsl_rng_alloc(T);

  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(ndim);

  gsl_monte_vegas_integrate(&G, xl, xu, ndim, 10000, r, s, &res, &err);
  // display_results("vegas warm-up", res, err);

  // printf("converging...\n");

  do {
    gsl_monte_vegas_integrate(&G, xl, xu, ndim, calls / 5, r, s, &res, &err);

  } while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);

  // display_results("vegas final", res, err);

  gsl_monte_vegas_free(s);

  gsl_rng_free(r);

  return res;
}
