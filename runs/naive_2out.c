#include <default.h>
#include <minimize.h>
#include <network.h>
#include <ran2.h>
#include <time.h>
long idum = -2;

typedef struct {
  int nX;
  double *x;
  multilayerD *net;
} minim_par;

void exact(double x, double *res);

void Full_chi2(void *addPar, double *c2, double *grad) {
  minim_par *mp = (minim_par *)addPar;
  multilayerD *net = mp->net;
  double res = 0.0;

  for (int j = 0; j < net->nPar; j++) {
    grad[j] = 0.0;
  }

  double div = (double)mp->nX + 1.0;
  for (int i = 0; i < mp->nX; i++) {
    double x = mp->x[i];
    multilD_FullEvaluate(net, &x);

    double temp1 = multilD_get_d2(0, 0, 0, net) + x * multilD_get_d(0, 0, net);
    temp1 += 4.0 * sin(2 * x) + cos(x) - 2.0 * x * cos(2 * x) + x * sin(x);

    double temp2 = x * multilD_get_d2(1, 0, 0, net) - x * exp(x * 0.5) / 8.0;

    res += temp1 * temp1 + temp2 * temp2;

    for (int j = 0; j < net->nPar; j++) {
      grad[j] += 2.0 * temp1 *
                 (multilD_get_grad_d2(0, 0, 0, j, net) +
                  x * multilD_get_grad_d(0, 0, j, net));
      grad[j] += 2.0 * temp2 * (x * multilD_get_grad_d2(1, 0, 0, j, net));
    }
  }

  double x = 0.0;
  multilD_FullEvaluate(net, &x);
  double temp1 = multilD_get(0, net) - 1.0;
  double temp2 = multilD_get_d(0, 0, net) - 2.0;
  double temp3 = multilD_get(1, net);
  double temp4 = multilD_get_d(1, 0, net) - 0.25;

  res += temp1 * temp1;
  res += temp2 * temp2;
  res += temp3 * temp3;
  res += temp4 * temp4;
  for (int j = 0; j < net->nPar; j++) {
    grad[j] += 2.0 * temp1 * (multilD_get_grad(0, j, net));
    grad[j] += 2.0 * temp2 * (multilD_get_grad_d(0, 0, j, net));
    grad[j] += 2.0 * temp3 * (multilD_get_grad(1, j, net));
    grad[j] += 2.0 * temp4 * (multilD_get_grad_d(1, 0, j, net));
  }

  *c2 = res;
}

double sin_act(double x) { return sin(x); }
double sin_act_d(double x) { return cos(x); }
double sin_act_d2(double x) { return -sin(x); }
double sin_act_d3(double x) { return -cos(x); }

int main() {
  int nI = 1;
  int nO = 2;
  int arch[6] = {nI, 10, 10, 10, 10, nO};
  int nL = 6;
  int nPar = multilD_getNpar(nL, arch);
  multilayerD net = multilD_init_net(nL, arch);
  // multilD_set_act_one(&net, 1, sin_act, sin_act_d, sin_act_d2, sin_act_d3);
  // multilD_set_act_one(&net, 2, sin_act, sin_act_d, sin_act_d2, sin_act_d3);

  minim_par mp;
  mp.net = &net;

  mp.nX = 100;
  mp.x = (double *)malloc(sizeof(double) * mp.nX);
  for (int i = 0; i < mp.nX; i++) {
    mp.x[i] = 2.0 * (double)i / ((double)mp.nX - 1.0);
  }

  double chiMin = minimize(net.par, nPar, (void *)(&mp), Full_chi2);
  double grad[nPar];
  double energy;
  Full_chi2((void *)(&mp), &energy, grad);
  FILE *fp = fopen("data.dat", "w");
  double x, ex[2];
  for (int i = 0; i < 10000; i++) {
    x = 2.0 * (double)i / 9999.0;
    exact(x, ex);

    multilD_Evaluate(&net, &x);
    double re = multilD_get(0, &net);
    double im = multilD_get(1, &net);
    fprintf(fp, "%e\t%e\t%e\t%e\t%e\n", x, re, im, ex[0], ex[1]);
  }

  multilD_save_net(&net, "networkPar.dat");

  multilD_free_net(&net);

  return 0;
}

void exact(double x, double *res) {
  res[0] = sin(2 * x) + cos(x);
  res[1] = exp(0.5 * x) * 0.5 - 0.5;
}
