#include <gsl/gsl_sf_hyperg.h>

#define kl 2
#define ll 3
#define alpha 1

double potential(double x) {
  return 0.5 * pow(alpha, 2) *
         (kl * (kl - 1) / pow(sin(alpha * x), 2) +
          ll * (ll - 1) / pow(cos(alpha * x), 2));
}
double ground_state_energy() { return 0.5 * pow(alpha, 2) * pow(kl + ll, 2); }
double exact(double x) {
  return pow(cos(alpha * x), ll) * pow(sin(alpha * x), kl) /
         sqrt(3.0 * M_PI / 512.0);
}

#undef V0
#undef kl
#undef ll
#undef alpha
