#include <activation.h>
#include <default.h>
#include <minimize.h>
#include <network.h>
#include <ran2.h>
#include <simAnnealing.h>
#include <time.h>
long idum = -2;

typedef struct {
  int nX;
  double *x;
  multilayerD *net;
} minim_par;

#include <semi_bound_potential.h>

double xmin = -2.0;
double xmax = 11.5;

double get(multilayerD *net) { return multilD_get(0, net); }
double get_grad(int j, multilayerD *net) { return multilD_get_grad(0, j, net); }
double get_d(multilayerD *net) { return multilD_get_d(0, 0, net); }
double get_grad_d(int j, multilayerD *net) {
  return multilD_get_grad_d(0, 0, j, net);
}
double get_d2(multilayerD *net) { return multilD_get_d2(0, 0, 0, net); }
double get_grad_d2(int j, multilayerD *net) {
  return multilD_get_grad_d2(0, 0, 0, j, net);
}

void Full_chi2(void *addPar, double *c2, double *grad) {
  minim_par *mp = (minim_par *)addPar;
  multilayerD *net = mp->net;
  double res = 0.0;

  double psi2_grad[net->nPar];
  double hpsi_grad[net->nPar];

  for (int j = 0; j < net->nPar; j++) {
    grad[j] = 0.0;
    psi2_grad[j] = 0.0;
    hpsi_grad[j] = 0.0;
  }
  double norm = 0.0;
  double en = 0.0;
  double dx = mp->x[1] - mp->x[0];

  double div = (double)mp->nX + 1.0;
  for (int i = 0; i < mp->nX; i++) {
    double x = mp->x[i];
    multilD_FullEvaluate(net, &x);

    double psi = get(net);

    double h_psi = -0.5 * get_d2(net) + potential(x) * psi;
    if (isnan(h_psi) != 0) {
      printf("err %lf \t %lf \t %lf\n", get_d2(net), potential(x), x);
      exit(-1);
    }

    en += psi * h_psi * dx;
    norm += psi * psi * dx;
    for (int j = 0; j < net->nPar; j++) {

      psi2_grad[j] += 2.0 * psi * get_grad(j, net) * dx;
      hpsi_grad[j] +=
          psi * dx *
          (-0.5 * get_grad_d2(j, net) + potential(x) * get_grad(j, net));
      hpsi_grad[j] += h_psi * get_grad(j, net) * dx;
    }
  }
  en /= norm;
  double diff = norm - 1.0;

  // printf("%e\n", norm);
  *c2 = pow(en, 2) + pow(diff, 2); //

  for (int j = 0; j < net->nPar; j++) {
    grad[j] = 2.0 * en * (hpsi_grad[j] / norm - en * psi2_grad[j] / norm);
    grad[j] += 2.0 * diff * psi2_grad[j];
  }
}

double test_solution(void *addPar, double en) {
  minim_par *mp = (minim_par *)addPar;
  multilayerD *net = mp->net;
  double res = 0.0;
  double n = 1e+5;
  double dx = (xmax - xmin) / (n - 1);
  for (int i = 0; i < n; i++) {
    double x = xmin + i * dx;
    multilD_FullEvaluate(net, &x);

    double psi = multilD_get(0, net);
    double h_psi = -0.5 * multilD_get_d2(0, 0, 0, net) + potential(x) * psi;

    res += pow(h_psi - en * psi, 2);
  }
  return res / n;
}

double get_energy(multilayerD *net) {
  double res = 0.0;
  double norm = 0.0;
  double en = 0.0;
  double n = 1e+5;
  double dx = (xmax - xmin) / (n - 1);

  for (int i = 0; i < n; i++) {
    double x = xmin + i * dx;
    multilD_FullEvaluate(net, &x);

    double psi = multilD_get(0, net);
    double h_psi = -0.5 * multilD_get_d2(0, 0, 0, net) + potential(x) * psi;

    en += psi * h_psi * dx;
    norm += psi * psi * dx;
  }
  en /= norm;
  printf("Norm: %lf\n", norm);
  return en;
}

double get_norm(multilayerD *net) {
  double res = 0.0;
  double norm = 0.0;
  double n = 1e+5;
  double dx = (xmax - xmin) / (n - 1);

  for (int i = 0; i < n; i++) {
    double x = xmin + i * dx;
    multilD_FullEvaluate(net, &x);

    double psi = multilD_get(0, net);

    norm += psi * psi * dx;
  }
  return norm;
}

int main() {
  int nI = 1;
  int nO = 1;
  int arch[8] = {nI, 10, 10, 10, 10, 10, 10, nO};
  int nL = 8;
  int nPar = multilD_getNpar(nL, arch);
  multilayerD net = multilD_init_net(nL, arch);
  // multilD_set_act_one(&net, 0, act_map, act_map_d, act_map_d2, act_map_d3);
  // multilD_load_net(&net, "networkPar.dat");
  minim_par mp;
  mp.net = &net;

  mp.nX = 1000;
  mp.x = (double *)malloc(sizeof(double) * mp.nX);
  for (int i = 0; i < mp.nX; i++) {
    mp.x[i] = xmin + (xmax - xmin) * (double)(i) / ((double)mp.nX - 1.0);
  }

  double chiMin = minimize(net.par, nPar, (void *)(&mp), Full_chi2);
  double energy = get_energy(&net);

  printf("Final enrgy: %lf\t Expected: %lf\n", energy, ground_state_energy());
  printf("Check sol: %e\n", test_solution((void *)(&mp), energy));
  FILE *fp = fopen("harmOsc_nl7_10xLayer.dat", "w");
  double norm = get_norm(&net);
  double x;
  for (int i = 0; i < 1000; i++) {
    x = xmin + (xmax - xmin) * (double)(i) / (999.0);
    multilD_FullEvaluate(&net, &x);
    double res = multilD_get(0, &net);
    fprintf(fp, "%e\t%e\t%e\n", x, res / sqrt(norm), exact(x));
  }

  multilD_save_net(&net, "networkPar.dat");

  multilD_free_net(&net);

  return 0;
}
