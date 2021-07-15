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

void exact(double *x, double *res);

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

    double re_psi_t = multilD_get_d(0, 1, net);
    double im_psi_t = multilD_get_d(1, 1, net);
    double re_hpsi = -0.5 * multilD_get_d2(0, 0, 0, net) +
                     0.5 * pow(mp->xt[i][0], 2) * multilD_get(0, net);
    double im_hpsi = -0.5 * multilD_get_d2(1, 0, 0, net) +
                     0.5 * pow(mp->xt[i][0], 2) * multilD_get(1, net);

    double temp1 = re_psi_t - im_hpsi;
    double temp2 = im_psi_t + re_hpsi;

    res += pow(temp1, 2) + pow(temp2, 2);

    for (int j = 0; j < net->nPar; j++) {
      grad[j] += 2.0 * temp1 *
                 (multilD_get_grad_d(0, 1, j, net) -
                  (-0.5) * multilD_get_grad_d2(1, 0, 0, j, net) -
                  0.5 * pow(mp->xt[i][0], 2) * multilD_get_grad(1, j, net));
      grad[j] += 2.0 * temp2 *
                 (multilD_get_grad_d(1, 1, j, net) +
                  (-0.5) * multilD_get_grad_d2(0, 0, 0, j, net) +
                  0.5 * pow(mp->xt[i][0], 2) * multilD_get_grad(0, j, net));
    }
  }

  for (int i = 0; i < mp->nBC; i++) {
    multilD_FullEvaluate(net, mp->xbc[i]);

    double bc[2];
    exact(mp->xbc[i], bc);

    double temp1 = multilD_get(0, net) - bc[0];
    double temp2 = multilD_get(1, net) - bc[1];

    res += pow(temp1, 2) + pow(temp2, 2);

    for (int j = 0; j < net->nPar; j++) {
      grad[j] += 2.0 * temp1 * multilD_get_grad(0, j, net);
      grad[j] += 2.0 * temp2 * multilD_get_grad(1, j, net);
    }
  }

  for (int i = 0; i < mp->nT; i++) {
    mp->tbc[i][0] = -6;
    multilD_FullEvaluate(net, mp->tbc[i]);

    double temp1 = multilD_get(0, net);
    double temp2 = multilD_get(1, net);

    res += pow(temp1, 2) + pow(temp2, 2);

    for (int j = 0; j < net->nPar; j++) {
      grad[j] += 2.0 * temp1 * multilD_get_grad(0, j, net);
      grad[j] += 2.0 * temp2 * multilD_get_grad(1, j, net);
    }

    mp->tbc[i][0] = 6;
    multilD_FullEvaluate(net, mp->tbc[i]);

    temp1 = multilD_get(0, net);
    temp2 = multilD_get(1, net);

    res += pow(temp1, 2) + pow(temp2, 2);

    for (int j = 0; j < net->nPar; j++) {
      grad[j] += 2.0 * temp1 * multilD_get_grad(0, j, net);
      grad[j] += 2.0 * temp2 * multilD_get_grad(1, j, net);
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
  int nO = 2;
  int arch[6] = {nI, 20, 20, 20, 20, nO};
  int nL = 6;
  int nPar = multilD_getNpar(nL, arch);
  multilayerD net = multilD_init_net(nL, arch);
  // multilD_set_act_one(&net, 1, sin_act, sin_act_d, sin_act_d2, sin_act_d3);
  // multilD_set_act_one(&net, 2, sin_act, sin_act_d, sin_act_d2, sin_act_d3);

  minim_par mp;
  mp.net = &net;

  int grid = 20;
  mp.nD = grid * grid * 2;
  mp.xt = (double **)malloc(sizeof(double *) * mp.nD);
  for (int i = 0; i < grid * grid; i++) {
    mp.xt[i] = (double *)malloc(sizeof(double) * nI);
    mp.xt[i + grid * grid] = (double *)malloc(sizeof(double) * nI);
    int i1 = i % grid;
    int i2 = (i - i1) / grid;
    mp.xt[i][1] = 1.0 * (double)i1 / ((double)grid - 1);
    mp.xt[i][0] = -6.0 + 12.0 * (double)i2 / ((double)grid - 1);
    // printf("%lf\t%lf\n", mp.xt[i][0], mp.xt[i][1]);
  }

  mp.nBC = 100;
  mp.xbc = (double **)malloc(sizeof(double *) * mp.nBC);
  for (int i = 0; i < mp.nBC; i++) {
    mp.xbc[i] = (double *)malloc(sizeof(double) * nI);
    mp.xbc[i][1] = 0.0;
    mp.xbc[i][0] = -6.0 + 12.0 * (double)i / ((double)mp.nBC - 1.0);
  }

  mp.nT = 100;
  mp.tbc = (double **)malloc(sizeof(double *) * mp.nT);
  for (int i = 0; i < mp.nT; i++) {
    mp.tbc[i] = (double *)malloc(sizeof(double) * nI);
    mp.tbc[i][1] = 1.0 * (double)i / ((double)mp.nT - 1);
    mp.tbc[i][0] = -6.0;
  }

  for (int j = 0; j < 20; j++) {
    for (int i = grid * grid; i < 2 * grid * grid; i++) {
      // int i1 = i % grid;
      // int i2 = (i - i1) / grid;
      mp.xt[i][1] = ran2(&idum);
      mp.xt[i][0] = 12.0 * ran2(&idum) - 6.0;
      // printf("%lf\t%lf\n", mp.xt[i][0], mp.xt[i][1]);
    }

    if (j == 0) {
      minimize(net.par, nPar, (void *)(&mp), Full_chi2);
    } else {
      minimize2(net.par, nPar, (void *)(&mp), Full_chi2);
    }
  }

  // Save res
  FILE *fp = fopen("data.dat", "w");

  double x[2], ex[2];
  double norm = 0.0;
  grid = 100;
  for (int i = 0; i < grid * grid; i++) {
    int i1 = i % grid;
    int i2 = (i - i1) / grid;
    x[1] = 1.0 * (double)i1 / ((double)grid - 1);
    x[0] = -6.0 + 12.0 * (double)i2 / ((double)grid - 1);

    exact(x, ex);

    multilD_Evaluate(&net, x);
    double re = multilD_get(0, &net);
    double im = multilD_get(1, &net);
    fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\n", x[0], x[1], re, im, ex[0], ex[1]);

    norm += fabs(re - ex[0]) + fabs(im - ex[1]);
  }
  fclose(fp);

  printf("%e\n", norm / (grid * grid * 2));

  fp = fopen("data_bc.dat", "w");
  for (int i = 0; i < 10000; i++) {

    x[0] = -4.0 + 8.0 * (double)i / 9999.0;
    x[1] = 0.0;
    exact(x, ex);

    multilD_Evaluate(&net, x);
    double re = multilD_get(0, &net);
    double im = multilD_get(1, &net);
    fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\n", x[0], x[1], re, im, ex[0], ex[1]);
  }
  fclose(fp);

  multilD_save_net(&net, "networkPar.dat");

  multilD_free_net(&net);

  return 0;
}

void exact(double *x, double *res) {
  static double temp1 = 0.0;
  static double phi = 0.5;
  temp1 = exp(-x[0] * x[0] * 0.5) / sqrt(sqrt(M_PI));
  res[0] = (cos(0.5 * x[1] + phi) + cos(1.5 * x[1] + phi) * sqrt(2.0) * x[0]) *
           temp1;
  res[1] = (-sin(0.5 * x[1] + phi) - sin(1.5 * x[1] + phi) * sqrt(2.0) * x[0]) *
           temp1;
}
