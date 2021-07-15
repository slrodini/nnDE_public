double potential(double x) { return 0.5 * x * x; }
double exact(double x) { return exp(-x * x * 0.5) / sqrt(sqrt(M_PI)); }
double ground_state_energy() { return 0.5; }

// Works
