#include <default.h>
#include <minimize.h>
#include <network.h>
#include <ran2.h>
#include <time.h>
long idum = -2;

typedef struct {
  int nT;
  double *t;
  trilayer *net;
} minim_par;

double re_exact(double t) {
  double gammaj = 2.0;
  double omegaj2 = 3.0;
  double eEoM = 1.0;
  double omega = sqrt(omegaj2 + 0.1);
  return -eEoM *
         (cos(omega * t) * (omegaj2 - omega * omega) +
          sin(omega * t) * 2 * gammaj * omega) /
         (pow(omegaj2 - omega * omega, 2) +
          4 * gammaj * gammaj * omega * omega);
}

double im_exact(double t) {
  double gammaj = 2.0;
  double omegaj2 = 3.0;
  double eEoM = 1.0;
  double omega = sqrt(omegaj2 + 0.1);
  return -eEoM *
         (-sin(omega * t) * (omegaj2 - omega * omega) +
          cos(omega * t) * 2 * gammaj * omega) /
         (pow(omegaj2 - omega * omega, 2) +
          4 * gammaj * gammaj * omega * omega);
}

void Full_chi2(void *addPar, double *c2, double *grad) {
  minim_par *mp = (minim_par *)addPar;
  trilayer *net = mp->net;
  double res = 0.0;

  for (int j = 0; j < net->nPar; j++) {
    grad[j] = 0.0;
  }
  double gammaj = 2.0;
  double omegaj2 = 3.0;
  double eEoM = 1.0;
  double omega = sqrt(omegaj2 + 0.1);
  for (int j = 0; j < mp->nT; j++) {
    double norm = 0.0;

    // mp->t[j];

    tril_FullEvaluate(net, &(mp->t[j]));

    double r_re = net->res[0];
    double r_im = net->res[1];
    // norm += psi_re * psi_re + psi_im * psi_im;

    double d2r_re = tril_get_d2(0, 0, net);
    double d2r_im = tril_get_d2(1, 0, net);

    double d1r_re = tril_get_d(0, 0, net);
    double d1r_im = tril_get_d(1, 0, net);

    double T1 = d2r_re + 2.0 * gammaj * d1r_re + omegaj2 * r_re +
                eEoM * cos(mp->t[j] * omega);

    double T2 = d2r_im + 2.0 * gammaj * d1r_im + omegaj2 * r_im -
                eEoM * sin(mp->t[j] * omega);

    res += pow(T1, 2.0) + pow(T2, 2.0);
    for (int k = 0; k < net->nPar; k++) {
      grad[k] += T1 * (tril_get_grad(0, k, net) * omegaj2 +
                       2.0 * gammaj * tril_get_grad_d(0, 0, k, net) +
                       tril_get_grad_d2(0, 0, k, net));
      grad[k] += T2 * (tril_get_grad(1, k, net) * omegaj2 +
                       2.0 * gammaj * tril_get_grad_d(1, 0, k, net) +
                       tril_get_grad_d2(1, 0, k, net));
    }
  }

  res /= ((double)mp->nT + 1);
  for (int k = 0; k < net->nPar; k++) {
    grad[k] /= ((double)mp->nT + 1);
  }
  double ti = 0.0;
  tril_FullEvaluate(net, &ti);
  double r_re_0 = re_exact(0.0);
  double r_im_0 = im_exact(0.0);

  double T1 = net->res[0] - r_re_0;
  double T2 = net->res[1] - r_im_0;
  double T3 = tril_get_d(0, 0, net) - omega * r_im_0;
  double T4 = tril_get_d(1, 0, net) + omega * r_re_0;
  res += pow(T1, 2) + pow(T2, 2) + pow(T3, 2) + pow(T4, 2);
  for (int k = 0; k < net->nPar; k++) {
    grad[k] += T1 * tril_get_grad(0, k, net);
    grad[k] += T2 * tril_get_grad(1, k, net);
    grad[k] += T3 * tril_get_grad_d(0, 0, k, net);
    grad[k] += T4 * tril_get_grad_d(1, 0, k, net);
  }

  *c2 = res;
}

int main() {
  int nI = 1;
  int nO = 2;
  int nH1 = 30;
  int nH2 = 20;
  int nH3 = 10;
  int nPar = tril_getNpar(nI, nH1, nH2, nH3, nO);
  printf("%d\n", nPar);
  trilayer net = tril_init_net(nI, nH1, nH2, nH3, nO);

  minim_par mp;
  mp.net = &net;
  mp.nT = 100;
  mp.t = (double *)malloc(sizeof(double) * mp.nT);

  for (int i = 0; i < mp.nT; i++) {
    mp.t[i] = 1.0 * (double)i / ((double)mp.nT - 1.0);
  }

  adamPar aP;
  aP.alpha = 0.002;
  aP.beta1 = 0.9;
  aP.beta2 = 0.997;
  aP.dfToll = 1e-4;
  aP.maxIt = 1e+4;

  double chiMin = adam2(net.par, net.nPar, (void *)(&mp), Full_chi2, &aP, NULL);
  printf("Final chi: %e\n", chiMin);

  FILE *fp = fopen("data.dat", "w");
  for (int i = 0; i < 1000; i++) {
    double ti = (double)i / 999.0;

    tril_Evaluate(&net, &ti);
    double re = net.res[0];
    double im = net.res[1];
    fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n", ti, re, im, re_exact(ti),
            im_exact(ti));
  }

  FILE *fp_n = fopen("networkPar.dat", "w");
  for (int i = 0; i < nPar; i++) {
    fprintf(fp_n, "%lf\n", net.par[i]);
  }
  fclose(fp_n);

  tril_free_net(&net);

  return 0;
}
