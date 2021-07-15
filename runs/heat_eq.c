#include <default.h>
#include <minimize.h>
#include <network.h>
#include <ran2.h>
#include <time.h>
long idum = -2;

typedef struct {
  int nD, nBC, nT;
  double **xt;
  double **xbc;
  double **tbc;
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

  double div = (double)mp->nD + 1.0;
  for (int i = 0; i < mp->nD; i++) {
    multilD_FullEvaluate(net, mp->xt[i]);

    double f_t = multilD_get_d(0, 1, net);
    double f_xx = multilD_get_d2(0, 0, 0, net);

    double temp1 = f_t - f_xx / 4.0;

    res += pow(temp1, 2);

    for (int j = 0; j < net->nPar; j++) {
      double loc_val = 2.0 * temp1 *
                       (multilD_get_grad_d(0, 1, j, net) -
                        0.25 * multilD_get_grad_d2(0, 0, 0, j, net));
      grad[j] += loc_val;
    }
  }

  for (int i = 0; i < mp->nBC; i++) {
    multilD_FullEvaluate(net, mp->xbc[i]);

    double bc = exact(mp->xbc[i]);
    double temp1 = multilD_get(0, net) - bc;

    res += pow(temp1, 2);

    for (int j = 0; j < net->nPar; j++) {
      grad[j] += 2.0 * temp1 * multilD_get_grad(0, j, net);
    }
  }

  for (int i = 0; i < mp->nT; i++) {
    mp->tbc[i][0] = -4;
    multilD_FullEvaluate(net, mp->tbc[i]);

    double temp1 = multilD_get(0, net);

    res += pow(temp1, 2);

    for (int j = 0; j < net->nPar; j++) {
      grad[j] += 2.0 * temp1 * multilD_get_grad(0, j, net);
    }

    mp->tbc[i][0] = 4;
    multilD_FullEvaluate(net, mp->tbc[i]);

    temp1 = multilD_get(0, net);
    res += pow(temp1, 2);

    for (int j = 0; j < net->nPar; j++) {
      grad[j] += 2.0 * temp1 * multilD_get_grad(0, j, net);
    }
  }

  *c2 = res;
}

double sin_act(double x) { return sin(x); }
double sin_act_d(double x) { return cos(x); }
double sin_act_d2(double x) { return -sin(x); }
double sin_act_d3(double x) { return -cos(x); }

int main() {
  int nI = 2;
  int nO = 1;
  int arch[6] = {nI, 30, 30, 30, 30, nO};
  int nL = 6;
  int nPar = multilD_getNpar(nL, arch);
  multilayerD net = multilD_init_net(nL, arch);
  // multilD_set_act_one(&net, 3, sin_act, sin_act_d, sin_act_d2, sin_act_d3);
  // multilD_set_act_one(&net, 9, sin_act, sin_act_d, sin_act_d2, sin_act_d3);

  minim_par mp;
  mp.net = &net;

  int grid = 32;
  mp.nD = grid * grid;
  mp.xt = (double **)malloc(sizeof(double *) * mp.nD);
  for (int i = 0; i < mp.nD; i++) {
    mp.xt[i] = (double *)malloc(sizeof(double) * nI);
    int i1 = i % grid;
    int i2 = (i - i1) / grid;
    mp.xt[i][1] = 1.0 + 1.0 * (double)i1 / ((double)grid - 1);
    mp.xt[i][0] = -4.0 + 8.0 * (double)i2 / ((double)grid - 1);
    // printf("%lf\t%lf\n", mp.xt[i][0], mp.xt[i][1]);
  }

  mp.nBC = 200;
  mp.xbc = (double **)malloc(sizeof(double *) * mp.nBC);
  for (int i = 0; i < mp.nBC; i++) {
    mp.xbc[i] = (double *)malloc(sizeof(double) * nI);
    mp.xbc[i][1] = 1.0;
    mp.xbc[i][0] = -4.0 + 8.0 * (double)i / ((double)mp.nBC - 1.0);
  }

  mp.nT = 200;
  mp.tbc = (double **)malloc(sizeof(double *) * mp.nT);
  for (int i = 0; i < mp.nT; i++) {
    mp.tbc[i] = (double *)malloc(sizeof(double) * nI);
    mp.tbc[i][1] = 1.0 + 1.0 * (double)i / ((double)mp.nT - 1);
    mp.tbc[i][0] = -4.0;
  }

  double chiMin = minimize(net.par, nPar, (void *)(&mp), Full_chi2);

  FILE *fp = fopen("data.dat", "w");
  double x[2];
  double norm = 0.0;
  grid = 100;
  for (int i = 0; i < grid * grid; i++) {
    int i1 = i % grid;
    int i2 = (i - i1) / grid;
    x[1] = 1.0 + 1.0 * (double)i1 / ((double)grid - 1);
    x[0] = -4.0 + 8.0 * (double)i2 / ((double)grid - 1);
    double ex = exact(x);
    double f_t = multilD_get_d(0, 1, &net);
    double f_xx = multilD_get_d2(0, 0, 0, &net);
    norm += fabs(f_t - f_xx / 4.0);

    multilD_Evaluate(&net, x);
    double re = multilD_get(0, &net);
    fprintf(fp, "%e\t%e\t%e\t%e\n", x[0], x[1], re, ex);
  }
  fclose(fp);

  printf("%e\n", norm / (grid * grid));

  fp = fopen("data_bc.dat", "w");
  for (int i = 0; i < 10000; i++) {

    x[0] = -4.0 + 8.0 * (double)i / 9999.0;
    x[1] = 1.0;
    double ex = exact(x);

    multilD_Evaluate(&net, x);
    double re = multilD_get(0, &net);
    fprintf(fp, "%e\t%e\t%e\t%e\n", x[0], x[1], re, ex);
  }
  fclose(fp);

  multilD_save_net(&net, "networkPar.dat");

  multilD_free_net(&net);

  return 0;
}

double exact(double *x) {
  static double k = 0.25;
  return exp(-x[0] * x[0] / (4.0 * k * x[1])) / sqrt(4.0 * k * M_PI * x[1]);
}
