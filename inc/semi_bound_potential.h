// double potential(double x) { return 2.0 *(exp(-2.0 * x) -  exp(-x)) + 1.0; }
// //works with 50 50 nX10000 from -3 to 20
double potential(double x) { return exp(-2.0 * x) - 2.0 * exp(-x) + 1.0; }
double ground_state_energy() {
  double lambda = sqrt(2.0);
  return 1 - pow(lambda - 0.5, 2) / 2.0;
}
double exact(double x) {
  double lambda = sqrt(2.0);
  double z = 2 * lambda * exp(-x);
  return pow(z, lambda - 0.5) * exp(-z / 2.);
}
// works nx 3000, 30, 30 from -3 to 15
