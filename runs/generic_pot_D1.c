#include <default.h>
#include <minimize.h>
#include <network.h>
#include <ran2.h>
#include <simAnnealing.h>
#include <time.h>
long idum = -2;

typedef struct
{
  int nX;
  double *x;
  multilayerD1 *net;
} minim_par;

#include <harmonic_potential.h>

double xmin = -4.0;
double xmax = +4.0;

double get(multilayerD1 *net) { return multilD1_get(0, net); }
double get_grad(int j, multilayerD1 *net) { return multilD1_get_grad(0, j, net); }
double get_d(multilayerD1 *net) { return multilD1_get_d(0, 0, net); }
double get_grad_d(int j, multilayerD1 *net)
{
  return multilD1_get_grad_d(0, 0, j, net);
}

void Full_chi2(void *addPar, double *c2, double *grad)
{
  minim_par *mp = (minim_par *)addPar;
  multilayerD1 *net = mp->net;
  double res = 0.0;

  double psi2_grad[net->nPar];
  double psi_h_psi_grad[net->nPar];

  for (int j = 0; j < net->nPar; j++)
  {
    grad[j] = 0.0;
    psi2_grad[j] = 0.0;
    psi_h_psi_grad[j] = 0.0;
  }
  double norm = 0.0;
  double en = 0.0;
  double dx = mp->x[1] - mp->x[0];

  double div = (double)mp->nX + 1.0;
  for (int i = 0; i < mp->nX; i++)
  {
    double x = mp->x[i];
    multilD1_FullEvaluate(net, &x);

    double psi = get(net);

    double psi_h_psi = 0.5 * get_d(net) * get_d(net) + potential(x) * psi * psi;

    en += psi_h_psi * dx;
    norm += psi * psi * dx;
    for (int j = 0; j < net->nPar; j++)
    {
      double temp = 2.0 * psi * get_grad(j, net) * dx;
      psi2_grad[j] += temp;
      psi_h_psi_grad[j] += temp * potential(x);
      psi_h_psi_grad[j] += get_grad_d(j, net) * get_d(net) * dx;
    }
  }
  en /= norm;
  double diff = norm - 1.0;

  // printf("%e\n", norm);
  *c2 = en + pow(diff, 2);

  for (int j = 0; j < net->nPar; j++)
  {
    grad[j] = (psi_h_psi_grad[j] / norm - en * psi2_grad[j] / norm);
    grad[j] += 2.0 * diff * psi2_grad[j];
  }
}

double get_energy(multilayerD1 *net)
{
  double res = 0.0;
  double norm = 0.0;
  double en = 0.0;
  double n = 1e+5;
  double dx = (xmax - xmin) / (n - 1);

  for (int i = 0; i < n; i++)
  {
    double x = xmin + i * dx;
    multilD1_FullEvaluate(net, &x);

    double psi = multilD1_get(0, net);
    double psi_h_psi = 0.5 * get_d(net) * get_d(net) + potential(x) * psi * psi;
    en += psi_h_psi * dx;
    norm += psi * psi * dx;
  }
  en /= norm;
  printf("Norm: %lf\n", norm);
  return en;
}

double get_norm(multilayerD1 *net)
{
  double res = 0.0;
  double norm = 0.0;
  double n = 1e+5;
  double dx = (xmax - xmin) / (n - 1);

  for (int i = 0; i < n; i++)
  {
    double x = xmin + i * dx;
    multilD1_FullEvaluate(net, &x);

    double psi = multilD1_get(0, net);

    norm += psi * psi * dx;
  }
  return norm;
}

int main()
{
  clock_t begin = clock();
  int nI = 1;
  int nO = 1;
  int arch[8] = {nI, 10, 10, 10, 10, 10, 10, nO};
  int nL = 8;
  int nPar = multilD1_getNpar(nL, arch);
  multilayerD1 net = multilD1_init_net(nL, arch);

  minim_par mp;
  mp.net = &net;

  mp.nX = 10000;
  mp.x = (double *)malloc(sizeof(double) * mp.nX);
  for (int i = 0; i < mp.nX; i++)
  {
    mp.x[i] = xmin + (xmax - xmin) * (double)(i) / ((double)mp.nX - 1.0);
  }

  double chiMin = minimize(net.par, nPar, (void *)(&mp), Full_chi2);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Time elapsed in s: %lf\n", time_spent);

  double energy = get_energy(&net);

  printf("Final enrgy: %lf\t Expected: %lf\n", energy, ground_state_energy());

  multilD1_free_net(&net);

  return 0;
}
