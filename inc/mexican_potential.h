#include <gsl/gsl_sf_hyperg.h>

double potential(double x) { return pow(x, 4) - pow(x, 2) + 1; }
double ground_state_energy() { return 0.0; }
double exact(double x) { return 0.0; }

// Works with 70 70 70 70 & 1000 points of int
