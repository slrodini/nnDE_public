#include <simAnnealing.h>
extern long idum;

void help() {
  printf("annealing(double *par, int nPar, double *parErr, void *addPar,double "
         "(*fEn)(double *, void *),annealing_par *anPar, int *ranMode, const "
         "char fileName[]);\n annealing_many(double *par, int nPar, double "
         "*parErr, void *addPar,double (*fEn)(double *, void *), annealing_par "
         "*anPar);\n ");
}

void copy_state(double *src, double *dst, int n) {
  for (int i = 0; i < n; i++) {
    dst[i] = src[i];
  }
}

double ran2U(long *id) { return 2.0 * ran2(id) - 1.0; }

double annealing(double *par, int nPar, double *parErr, void *addPar,
                 double (*fEn)(double *, void *), annealing_par *anPar,
                 int *ranMode) {
  double T = anPar->Tin;
  double Tfin = anPar->Tfin;
  int nIter = anPar->nIter; // Iteration at fixed T
  double dt = anPar->dt;

  double best_x[nPar];
  double pErr[nPar];
  if (parErr == NULL) {
    for (int i = 0; i < nPar; i++)
      pErr[i] = 1.0;
  } else {
    for (int i = 0; i < nPar; i++)
      pErr[i] = parErr[i];
  }
  double E, new_E, best_E;

  E = fEn(par, addPar);
  copy_state(par, best_x, nPar);

  best_E = E;

  int count = 0;

  double ranVal[nIter], ranValVar[nIter];
  int indexes[nIter];
  double (*rFunc[nPar])(long *);

  int save = 0;

  if (ranMode != NULL) {
    for (int i = 0; i < nPar; i++) {
      rFunc[i] = (ranMode[i] == 0) ? ran2U : ran2N;
    }
  } else {
    for (int i = 0; i < nPar; i++) {
      rFunc[i] = ran2U;
    }
  }

  while (1) {
    for (int i = 0; i < nIter; i++) {
      indexes[i] = (int)(nPar * ran2(&idum));
      ranVal[i] = ran2(&idum);
      ranValVar[i] = rFunc[indexes[i]](&idum) * pErr[indexes[i]];
    }
    printf("Iteration: %d Temperature: %.6e Energy: %.6e\n", count, T, E);

    for (int i = 0; i < nIter; i++) {
      // copy_state(x, new_x, nPar);

      // step (uniform)
      int index = indexes[i];
      double u = ranValVar[i];
      // new_x[index] = u * 2 * parErr[index] - parErr[index] + new_x[index];
      par[index] = u + par[index];

      // new_E = fEn(new_x, addPar);
      new_E = fEn(par, addPar);
      if (new_E <= best_E) {
        // copy_state(new_x, best_x, nPar);
        copy_state(par, best_x, nPar);
        best_E = new_E;
      }
      // printf("%d\t%d\t%.6e\t%.6e\n", i, index, E, new_E);
      if (new_E == E) {
        par[index] = par[index] - u;
      } else if (new_E < E || ranVal[i] < exp((-new_E + E) / (T))) {
        // mi sembra superfluo
        // if (new_E < best_E) {
        // copy_state(new_x, best_x, n);
        //}
        // copy_state(new_x, x, nPar);
        E = new_E;
      } else {
        par[index] = par[index] - u;
      }
      // printf("%d\t%d\t%.6e\t%.6e\n", 0, i, T, E);

    } // end equal-T cycles
    count++;

    T = T * dt;
    if (T < Tfin)
      break;
  }
  copy_state(best_x, par, nPar);
  return best_E;
}
