#define lambda 1.0

double potential(double x) { return 0.5 * x * x + lambda * pow(x, 4); }
double exact(double x) { return exp(-x * x * 0.5) / sqrt(sqrt(M_PI)); }
double ground_state_energy() { return 0.5 + lambda * 3.0 / 4.0; }
#undef lambda
// Works
