#pragma once
#include <default.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <stdlib.h>

double vegas(double (*futil)(double *, size_t, void *), double *xl, double *xu,
             int ndim, void *addPar);
