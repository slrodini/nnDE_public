#include <default.h>
#include <minimize.h>
#include <network.h>
#include <ran2.h>
#include <time.h>
long idum = -2;

typedef struct {
  int nX;
  double **x;
  multilayerD *net;
} minim_par;

double exact(double *x);

void Full_chi2(void *addPar, double *c2, double *grad) {
  minim_par *mp = (minim_par *)addPar;
  multilayerD *net = mp->net;
  double res = 0.0;

  for (int j = 0; j < net->nPar; j++) {
    grad[j] = 0.0;
  }

  double div = (double)mp->nX + 1.0;
  for (int i = 0; i < mp->nX; i++) {
    multilD_FullEvaluate(net, mp->x[i]);

    double fx = multilD_get_d(0, 0, net);
    double fy = multilD_get_d(0, 1, net);

    double temp = fx * fx + fy * fy - 1.0;
    res += temp * temp;

    for (int j = 0; j < net->nPar; j++) {
      grad[j] += 4.0 * temp *
                 (fx * multilD_get_grad_d(0, 0, j, net) +
                  fy * multilD_get_grad_d(0, 1, j, net));
    }
  }

  double x[2];
  x[0] = 0.0;
  x[1] = 0.0;
  multilD_FullEvaluate(net, x);
  double temp = multilD_get(0, net);
  res += temp * temp;
  for (int j = 0; j < net->nPar; j++) {
    grad[j] += 2.0 * temp * multilD_get_grad(0, j, net);
  }

  x[0] = 1.0;
  x[1] = 0.0;
  multilD_FullEvaluate(net, x);
  temp = multilD_get(0, net) - 1.0 / sqrt(2.0);
  res += temp * temp;
  for (int j = 0; j < net->nPar; j++) {
    grad[j] += 2.0 * temp * multilD_get_grad(0, j, net);
  }

  x[0] = 0.0;
  x[1] = 1.0;
  multilD_FullEvaluate(net, x);
  temp = multilD_get(0, net) - 1.0 / sqrt(2.0);
  res += temp * temp;
  for (int j = 0; j < net->nPar; j++) {
    grad[j] += 2.0 * temp * multilD_get_grad(0, j, net);
  }

  *c2 = res;
}

int main() {
  int nI = 2;
  int nO = 1;
  int arch[4] = {nI, 20, 10, nO};
  int nL = 4;
  int nPar = multilD_getNpar(nL, arch);
  multilayerD net = multilD_init_net(nL, arch);

  minim_par mp;
  mp.net = &net;

  mp.nX = 100;
  mp.x = (double **)malloc(sizeof(double *) * mp.nX);
  for (int i = 0; i < mp.nX; i++) {
    mp.x[i] = (double *)malloc(sizeof(double) * nI);
    int i1 = i % 10;
    int i2 = (i - i1) / 10;
    mp.x[i][0] = 1.0 * (double)i1 / 10.0;
    mp.x[i][1] = 1.0 * (double)i2 / 10.0;
  }

  double chiMin = minimize(net.par, nPar, (void *)(&mp), Full_chi2);
  double grad[nPar];
  double energy;
  Full_chi2((void *)(&mp), &energy, grad);
  FILE *fp = fopen("data.dat", "w");
  double x[2];
  for (int i = 0; i < 10000; i++) {
    int i1 = i % 100;
    int i2 = (i - i1) / 100;
    x[0] = 1.0 * (double)i1 / 100.0;
    x[1] = 1.0 * (double)i2 / 100.0;

    multilD_Evaluate(&net, x);
    double res = multilD_get(0, &net);
    fprintf(fp, "%e\t%e\t%e\t%e\n", x[0], x[1], res, exact(x));
  }

  multilD_save_net(&net, "networkPar.dat");

  multilD_free_net(&net);

  return 0;
}

double exact(double *x) { return (x[0] + x[1]) / sqrt(2.0); }
