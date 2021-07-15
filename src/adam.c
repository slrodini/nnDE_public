#include <adam.h>
#define EPS 1e-8

double adamax(double *par, int nPar, void *addPar,
              void (*fgrad)(void *, double *, double *), adamPar *adPar) {

  double alpha = adPar->alpha;
  double beta1 = adPar->beta1;
  double beta2 = adPar->beta2;
  int maxIt = adPar->maxIt;
  double gradToll = adPar->dfToll;

  double mt[nPar], ut[nPar], gt[nPar]; //, locPar[nPar];

  for (int i = 0; i < nPar; i++) {
    mt[i] = 0;
    ut[i] = 0;
    gt[i] = 0;
    // locPar[i] = par[i];
  }

  double t = 0;
  double oldChi = 0.0;
  fgrad(addPar, &oldChi, gt);

  double gradNorm;

  while (1) {
    t++;
    if ((int)t % adPar->changeLearn == 0) {
      alpha /= 1.5;
      beta2 = pow(beta2, 2.0 / 3.0);
    }
    // Compute the gradient

    fgrad(addPar, &oldChi, gt);

    gradNorm = 0.0;
    for (int i = 0; i < nPar; i++) {
      gradNorm += fabs(gt[i]);
    }

    printf("Iteration: %d\tcurrent_chi2: %.4e\tcurrent_grad: %.4e\n", (int)t,
           oldChi, gradNorm);

    double temp = alpha / (1.0 - pow(beta1, t));
    for (int i = 0; i < nPar; i++) {
      mt[i] = beta1 * mt[i] + (1.0 - beta1) * gt[i];
      ut[i] = (beta2 * ut[i] > fabs(gt[i])) ? beta2 * ut[i] : fabs(gt[i]);
      double change = -temp * mt[i] / (ut[i] + EPS);
      // locPar[i] = locPar[i] + change;
      par[i] = par[i] + change;
    }

    if (t >= maxIt || gradNorm <= gradToll)
      break;
  }

  // for (int i = 0; i < nPar; i++)
  //  par[i] = locPar[i];

  return oldChi;
}

double adam2(double *par, int nPar, void *addPar,
             void (*fgrad)(void *, double *, double *), adamPar *adPar) {

  double alpha = adPar->alpha;
  double beta1 = adPar->beta1;
  double beta2 = adPar->beta2;
  int maxIt = adPar->maxIt;
  double gradToll = adPar->dfToll;

  double mt[nPar], vt[nPar], gt[nPar]; //, locPar[nPar];

  for (int i = 0; i < nPar; i++) {
    mt[i] = 0;
    vt[i] = 0;
    gt[i] = 0;
    // locPar[i] = par[i];
  }
  double t = 0;
  double oldChi = 0.0;
  fgrad(addPar, &oldChi, gt);

  double gradNorm = 0;

  while (1) {
    t++;

    if ((int)t % adPar->changeLearn == 0) {
      alpha /= 1.5;
      beta2 = pow(beta2, 2.0 / 3.0);
    }

    fgrad(addPar, &oldChi, gt);
    gradNorm = 0.0;
    for (int i = 0; i < nPar; i++) {
      gradNorm += gt[i] * gt[i];
    }
    gradNorm = sqrt(gradNorm);
    printf("adam2 Iter: %4d  cur_chi2: %.4e  cur_grad: %.4e\n", (int)t, oldChi,
           gradNorm);

    double temp1 = 1.0 / (1.0 - pow(beta1, t));
    double temp2 = 1.0 / (1.0 - pow(beta2, t));

    for (int i = 0; i < nPar; i++) {

      mt[i] = beta1 * mt[i] + (1.0 - beta1) * gt[i];
      vt[i] = beta2 * vt[i] + (1.0 - beta2) * gt[i] * gt[i];
      double mhat = mt[i] * temp1;
      double vhat = vt[i] * temp2;

      double change = -alpha * mhat / (sqrt(fabs(vhat)) + EPS);

      // if (isnan(change) != 0) {
      //  printf("adam\n");
      //  exit(-1);
      //}

      // locPar[i] = locPar[i] + change;
      par[i] = par[i] + change;
    }

    if (t >= maxIt || gradNorm <= gradToll)
      break;
  }

  // for (int i = 0; i < nPar; i++)
  //  par[i] = locPar[i];
  return oldChi;
}
