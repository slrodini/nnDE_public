#include <activation.h>

#define aC 1.1745
#define mC (2.0 / 3.0)
#define bC 0.1

double act_sigma(double x) { return aC * tanh(x * mC) + bC * x; }
double act_sigma_d(double x) {
  if (fabs(x) > 100) {
    return 0.0;
  }
  return (2.0 * aC * mC) / (1.0 + cosh(2.0 * mC * x)) + bC;
}
double act_sigma_d2(double x) {
  if (fabs(x) > 100) {
    return 0.0;
  }
  return (-8.0 * aC * pow(mC, 2) * sinh(mC * x)) /
         (3.0 * cosh(mC * x) + cosh(3.0 * mC * x));
}
double act_sigma_d3(double x) {
  if (fabs(x) > 100) {
    return 0;
  }

  return 2.0 *
         (-2.0 * aC * pow(mC, 3) * pow(1.0 / cosh(mC * x), 4) +
          aC * pow(mC, 3) * cosh(2.0 * mC * x) * pow(1.0 / cosh(mC * x), 4));
}

double act_sin(double x) { return sin(x); }
double act_sin_d(double x) { return cos(x); }
double act_sin_d2(double x) { return -sin(x); }
double act_sin_d3(double x) { return -cos(x); }
#undef aC
#undef mC
#undef bC
