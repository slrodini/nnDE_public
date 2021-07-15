#define ee 0.5
double potential(double x) { return 0.5 * x * x - ee * x; }
double exact(double x) {
  return exp(-(x - ee) * (x - ee) * 0.5) / sqrt(sqrt(M_PI));
}
double ground_state_energy() { return 0.5 - 0.5 * ee * ee; }
#undef ee

// Works
