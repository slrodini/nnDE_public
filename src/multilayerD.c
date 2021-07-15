#include <network.h>
#include <ran2.h>
extern long idum;

// Debug stuff
/*
void check(double val, const char name[]) {
  if (isnan(val) != 0) {
    printf("%s\n", name);
    exit(-1);
  }
}

void check2(double val, const char name[], int l) {
  if (isnan(val) != 0) {
    printf("%s %d\n", name, l);
    exit(-1);
  }
}
*/
// end debug

//#define aC 1.7159
//#define mC (2.0 / 3.0)
//#define bC 0.001

#define aC 1.1745
#define mC (2.0 / 3.0)
#define bC 0.1

static double loc_sigma(double x) { return aC * tanh(x * mC) + bC * x; }
static double loc_sigma_d(double x) {
  if (fabs(x) > 100) {
    return 0.0;
  }
  return (2.0 * aC * mC) / (1.0 + cosh(2.0 * mC * x)) + bC;
}
static double loc_sigma_d2(double x) {
  if (fabs(x) > 100) {
    return 0.0;
  }
  return (-8.0 * aC * pow(mC, 2) * sinh(mC * x)) /
         (3.0 * cosh(mC * x) + cosh(3.0 * mC * x));
}
static double loc_sigma_d3(double x) {
  if (fabs(x) > 100) {
    return 0;
  }

  return 2.0 *
         (-2.0 * aC * pow(mC, 3) * pow(1.0 / cosh(mC * x), 4) +
          aC * pow(mC, 3) * cosh(2.0 * mC * x) * pow(1.0 / cosh(mC * x), 4));
}

static double loc_id(double x) { return x; }
static double loc_id_d(double x) { return 1.0; }
static double loc_id_d2(double x) { return 0.0; }
static double loc_id_d3(double x) { return 0.0; }

static double loc_cubic(double x) { return x * x * x * x; }
static double loc_cubic_d(double x) { return 4 * x * x * x; }
static double loc_cubic_d2(double x) { return 12 * x * x; }
static double loc_cubic_d3(double x) { return 24.0 * x; }

int multilD_getNpar(int nL, int *arch) {
  int nPar = 0;
  for (int i = nL - 1; i >= 1; i--) {
    nPar += arch[i] * (arch[i - 1] + 1);
  }
  return nPar;
}

void multilD_set_act(multilayerD *net, multiL_act *activ, multiL_act *d_activ,
                     multiL_act *d2_activ, multiL_act *d3_activ) {
  for (int i = 0; i < net->nL; i++) {
    if (activ[i] == NULL || d_activ[i] == NULL || d2_activ[i] == NULL ||
        d3_activ[i] == NULL) {
      net->act_fun[i] = loc_sigma;
      net->d_act_fun[i] = loc_sigma_d;
      net->d2_act_fun[i] = loc_sigma_d2;
      net->d3_act_fun[i] = loc_sigma_d3;
    } else {
      net->act_fun[i] = activ[i];
      net->d_act_fun[i] = d_activ[i];
      net->d2_act_fun[i] = d2_activ[i];
      net->d3_act_fun[i] = d3_activ[i];
    }
  }
}

void multilD_set_act_one(multilayerD *net, int i, multiL_act activ,
                         multiL_act d_activ, multiL_act d2_activ,
                         multiL_act d3_activ) {

  net->act_fun[i] = activ;
  net->d_act_fun[i] = d_activ;
  net->d2_act_fun[i] = d2_activ;
  net->d3_act_fun[i] = d3_activ;
}

multilayerD multilD_init_net(int nL, int *arch) {
  multilayerD net;
  net.nI = arch[0];
  net.nO = arch[nL - 1];
  net.nL = nL;

  net.input = (double *)calloc(arch[0], sizeof(double));
  // here I create a local copy of the architecture instead of copying the
  // reference
  net.arch = (int *)malloc(sizeof(int) * nL);
  net.offsetW = (int *)calloc(nL, sizeof(int));
  net.offsetB = (int *)calloc(nL, sizeof(int));
  net.maxNH = 0;
  for (int i = 0; i < nL; i++) {
    net.arch[i] = arch[i];
    if (net.arch[i] > net.maxNH)
      net.maxNH = arch[i];
  }

  net.offsetW[0] = 0;
  net.offsetB[0] = arch[1] * arch[0];
  for (int i = 1; i < nL - 1; i++) {
    net.offsetW[i] = net.offsetB[i - 1] + arch[i];
    net.offsetB[i] = net.offsetW[i] + arch[i] * arch[i + 1];
  }
  net.res_d = (double *)calloc(net.nI * net.nO, sizeof(double));
  net.res_d2 = (double **)malloc(net.nO * sizeof(double *));
  for (int i = 0; i < net.nO; i++) {
    net.res_d2[i] = (double *)calloc(net.nI * net.nI, sizeof(double));
  }

  net.nPar = multilD_getNpar(nL, arch);
  // net.res = (double *)calloc(nO, sizeof(double));
  net.par = (double *)calloc(net.nPar, sizeof(double));
  net.parGrad = (double *)calloc(net.nPar * net.nO, sizeof(double));
  net.parGrad_d = (double *)calloc(net.nPar * net.nO * net.nI, sizeof(double));
  net.parGrad_d2 =
      (double *)calloc(net.nPar * net.nO * net.nI * net.nI, sizeof(double));
  for (int i = 0; i < net.nPar; i++) {
    net.par[i] = 0.5 * ran2(&idum);
  }

  // init the activation functions to default
  net.act_fun = (multiL_act *)malloc(sizeof(multiL_act) * nL);
  net.d_act_fun = (multiL_act *)malloc(sizeof(multiL_act) * nL);
  net.d2_act_fun = (multiL_act *)malloc(sizeof(multiL_act) * nL);
  net.d3_act_fun = (multiL_act *)malloc(sizeof(multiL_act) * nL);
  // input and output layers with the identity function
  net.act_fun[0] = loc_id;
  net.d_act_fun[0] = loc_id_d;
  net.d2_act_fun[0] = loc_id_d2;
  net.d3_act_fun[0] = loc_id_d3;
  for (int i = 1; i < nL - 1; i++) {
    net.act_fun[i] = loc_sigma;
    net.d_act_fun[i] = loc_sigma_d;
    net.d2_act_fun[i] = loc_sigma_d2;
    net.d3_act_fun[i] = loc_sigma_d3;
  }
  net.act_fun[nL - 1] = loc_id;
  net.d_act_fun[nL - 1] = loc_id_d;
  net.d2_act_fun[nL - 1] = loc_id_d2;
  net.d3_act_fun[nL - 1] = loc_id_d3;

  net.li = (double **)malloc(sizeof(double *) * nL);
  net.s_li = (double **)malloc(sizeof(double *) * nL);
  net.sd_li = (double **)malloc(sizeof(double *) * nL);
  net.sd2_li = (double **)malloc(sizeof(double *) * nL);
  net.sd3_li = (double **)malloc(sizeof(double *) * nL);
  net.ds_li_dx = (double **)malloc(sizeof(double *) * nL);
  net.dsd_li_dx = (double **)malloc(sizeof(double *) * nL);

  net.d2s_li_dx2 = (double **)malloc(sizeof(double *) * nL);
  net.d2sd_li_dx2 = (double **)malloc(sizeof(double *) * nL);
  // net.d2sd2_li_dx2 = (double **)malloc(sizeof(double *) * nL);

  for (int i = 0; i < nL; i++) {
    net.li[i] = (double *)malloc(arch[i] * sizeof(double));
    net.s_li[i] = (double *)malloc(arch[i] * sizeof(double));
    net.sd_li[i] = (double *)malloc(arch[i] * sizeof(double));
    net.sd2_li[i] = (double *)malloc(arch[i] * sizeof(double));
    net.sd3_li[i] = (double *)malloc(arch[i] * sizeof(double));

    // To access them input + (nI)*j, j = 0,...,arch[l]-1
    net.ds_li_dx[i] = (double *)malloc(sizeof(double) * arch[i] * net.nI);
    net.dsd_li_dx[i] = (double *)malloc(sizeof(double) * arch[i] * net.nI);

    // To access them inputA + (nI)*inputB + (nI*nI)*j, j = 0,...,arch[l]-1
    net.d2s_li_dx2[i] =
        (double *)malloc(sizeof(double) * arch[i] * net.nI * net.nI);
    net.d2sd_li_dx2[i] =
        (double *)malloc(sizeof(double) * arch[i] * net.nI * net.nI);
  }
  return net;
}

void multilD_save_net(multilayerD *net, const char fileName[]) {
  FILE *fp = fopen(fileName, "w");
  int nPar = net->nPar;
  for (int i = 0; i < nPar; i++) {
    fprintf(fp, "%lf\n", net->par[i]);
  }
  fclose(fp);
}

void multilD_load_net(multilayerD *net, const char fileName[]) {
  FILE *fp = fopen(fileName, "r");
  int nPar = net->nPar;
  for (int i = 0; i < nPar; i++) {
    fscanf(fp, "%lf", (net->par + i));
  }
  fclose(fp);
}

void multilD_free_net(multilayerD *net) {
  if (net->input != NULL) {
    free(net->input);
  }

  free(net->arch);

  free(net->offsetW);
  free(net->offsetB);
  // free(net->res);

  free(net->par);
  free(net->parGrad);
  free(net->parGrad_d);
  free(net->parGrad_d2);

  free(net->res_d);

  for (int i = 0; i < net->nO; i++) {
    free(net->res_d2[i]);
  }
  free(net->res_d2);

  free(net->act_fun);
  free(net->d_act_fun);
  free(net->d2_act_fun);
  free(net->d3_act_fun);

  for (int i = 0; i < net->nL; i++) {
    free(net->li[i]);
    free(net->s_li[i]);
    free(net->sd_li[i]);
    free(net->sd2_li[i]);
    free(net->sd3_li[i]);
    free(net->ds_li_dx[i]);
    free(net->dsd_li_dx[i]);
    free(net->d2s_li_dx2[i]);
    free(net->d2sd_li_dx2[i]);
  }
  free(net->li);
  free(net->s_li);
  free(net->sd_li);
  free(net->sd2_li);
  free(net->sd3_li);
  free(net->ds_li_dx);
  free(net->dsd_li_dx);
  free(net->d2s_li_dx2);
  free(net->d2sd_li_dx2);
}

static double Wjk(int j, int k, int L, multilayerD *net) {
  double temp = net->par[net->offsetW[L - 1] + k + (net->arch[L - 1]) * j];
  // check(temp, "wjk");
  return temp;
}

static double Bj(int j, int L, multilayerD *net) {
  double temp = net->par[net->offsetB[L - 1] + j];
  // check(temp, "bj");
  return temp;
}

void multilD_Evaluate(multilayerD *net) {

  for (int i = 0; i < net->nI; i++) {
    net->li[0][i] = net->input[i];
    net->s_li[0][i] = net->act_fun[0](net->input[i]);
    net->sd_li[0][i] = net->d_act_fun[0](net->input[i]);
    net->sd2_li[0][i] = net->d2_act_fun[0](net->input[i]);
    net->sd3_li[0][i] = net->d3_act_fun[0](net->input[i]);

    for (int a = 0; a < net->nI; a++) {
      net->ds_li_dx[0][a + i * (net->nI)] = (i == a ? net->sd_li[0][i] : 0.0);
      net->dsd_li_dx[0][a + i * (net->nI)] = (i == a ? net->sd2_li[0][i] : 0.0);

      for (int b = 0; b <= a; b++) {
        double temp1 = 0.0;
        double temp2 = 0.0;
        if (a == b && a == i) {
          temp1 = net->sd2_li[0][i];
          temp2 = net->sd3_li[0][i];
        }
        int index1 = a + b * (net->nI) + i * (net->nI) * (net->nI);
        int index2 = b + a * (net->nI) + i * (net->nI) * (net->nI);
        net->d2s_li_dx2[0][index1] = temp1;
        net->d2sd_li_dx2[0][index1] = temp2;

        net->d2s_li_dx2[0][index2] = temp1;
        net->d2sd_li_dx2[0][index2] = temp2;
      }
    }
  }
  double temp[net->nI];
  for (int L = 1; L < net->nL; L++) {
    for (int j = 0; j < net->arch[L]; j++) {
      net->li[L][j] = 0.0;
      for (int k = 0; k < net->arch[L - 1]; k++) {
        net->li[L][j] += net->s_li[L - 1][k] * Wjk(j, k, L, net);
      }
      net->li[L][j] += Bj(j, L, net);
      // check(net->li[L][j], "li");
      net->s_li[L][j] = net->act_fun[L](net->li[L][j]);
      net->sd_li[L][j] = net->d_act_fun[L](net->li[L][j]);
      net->sd2_li[L][j] = net->d2_act_fun[L](net->li[L][j]);
      net->sd3_li[L][j] = net->d3_act_fun[L](net->li[L][j]);
      // check2(net->sd3_li[L][j], "sd3_1", L);

      for (int a = 0; a < net->nI; a++) {
        temp[a] = 0.0;
        for (int i = 0; i < net->arch[L - 1]; i++) {
          temp[a] +=
              Wjk(j, i, L, net) * net->ds_li_dx[L - 1][a + i * (net->nI)];
        }
        net->ds_li_dx[L][a + j * (net->nI)] = temp[a] * net->sd_li[L][j];
        net->dsd_li_dx[L][a + j * (net->nI)] = temp[a] * net->sd2_li[L][j];
        for (int b = 0; b <= a; b++) {
          double temp1 = 0.0;
          for (int i = 0; i < net->arch[L - 1]; i++) {
            temp1 += Wjk(j, i, L, net) *
                     net->d2s_li_dx2[L - 1][b + a * (net->nI) +
                                            i * (net->nI) * (net->nI)];
          }
          int index1 = b + a * (net->nI) + j * (net->nI) * (net->nI);
          int index2 = a + b * (net->nI) + j * (net->nI) * (net->nI);
          net->d2s_li_dx2[L][index1] =
              temp[a] * temp[b] * net->sd2_li[L][j] + temp1 * net->sd_li[L][j];
          net->d2sd_li_dx2[L][index1] =
              temp[a] * temp[b] * net->sd3_li[L][j] + temp1 * net->sd2_li[L][j];

          net->d2s_li_dx2[L][index2] =
              temp[a] * temp[b] * net->sd2_li[L][j] + temp1 * net->sd_li[L][j];
          net->d2sd_li_dx2[L][index2] =
              temp[a] * temp[b] * net->sd3_li[L][j] + temp1 * net->sd2_li[L][j];
          // check(net->sd3_li[L][j], "sd3");
          // check(net->d2s_li_dx2[L][index1], "1.1");
          // check(net->d2sd_li_dx2[L][index1], "1.2");
          // check(net->d2s_li_dx2[L][index2], "2.1");
          // check(net->d2sd_li_dx2[L][index2], "2.2");
        }
      }
    }
  }
}
void multilD_EvaluateParGradient(multilayerD *net) {
  int count = 0;
  double sigma_vec[net->maxNH];
  double sigmaTemp_vec[net->maxNH];

  double sigmaDot_vec[net->maxNH][net->nI];
  double sigmaDotTemp_vec[net->maxNH][net->nI];

  double sigmaDot2_vec[net->maxNH][(net->nI) * (net->nI)];
  double sigmaDot2Temp_vec[net->maxNH][(net->nI) * (net->nI)];

  for (int i = 0; i < net->nO; i++) {
    for (int f = 0; f < net->maxNH; f++) {
      sigma_vec[f] = 0.0;
      for (int a = 0; a < net->nI; a++) {
        sigmaDot_vec[f][a] = 0.0;
        for (int b = 0; b < net->nI; b++) {
          sigmaDot2_vec[f][b + a * (net->nI)] = 0.0;
        }
      }
    }
    sigma_vec[i] = 1.0;
    // ok that sigmaDot = 0

    for (int L = net->nL - 1; L > 0; L--) {
      // gradient of the output w.r.t. the parameters
      for (int j = 0; j < net->arch[L]; j++) {
        int index = net->offsetB[L - 1] + i * (net->nPar) + j;
        net->parGrad[index] = sigma_vec[j] * net->sd_li[L][j];

        // net->parGrad[index] = net->parGrad[index] > 0 ? 1.0 : 1.0;
        // count++;
        for (int k = 0; k < net->arch[L - 1]; k++) {
          index = net->offsetW[L - 1] + i * (net->nPar) + k +
                  (net->arch[L - 1]) * j;
          net->parGrad[index] =
              sigma_vec[j] * net->sd_li[L][j] * net->s_li[L - 1][k];
          // net->parGrad[index] = net->parGrad[index] > 0 ? 1.0 : 1.0;
          // count++;
        }

        for (int a = 0; a < net->nI; a++) {
          // a * (net->nPar * net->nO) + i * net->nPar;
          int nCycle = a * (net->nPar * net->nO) + i * net->nPar;
          index = nCycle + net->offsetB[L - 1] + j;
          double foo = sigmaDot_vec[j][a] * net->sd_li[L][j] +
                       sigma_vec[j] * net->dsd_li_dx[L][a + j * (net->nI)];
          net->parGrad_d[index] = foo;
          for (int k = 0; k < net->arch[L - 1]; k++) {
            index = nCycle + net->offsetW[L - 1] + k + (net->arch[L - 1]) * j;
            net->parGrad_d[index] = foo * net->s_li[L - 1][k] +
                                    sigma_vec[j] * net->sd_li[L][j] *
                                        net->ds_li_dx[L - 1][a + k * (net->nI)];
            // net->parGrad[index] = net->parGrad[index] > 0 ? 1.0 : 1.0;
            // count++;
          }

          for (int b = 0; b <= a; b++) {
            int nCycle1 = b * (net->nPar * net->nO * net->nI) +
                          a * (net->nPar * net->nO) + i * net->nPar;
            int nCycle2 = a * (net->nPar * net->nO * net->nI) +
                          b * (net->nPar * net->nO) + i * net->nPar;
            int index1 = nCycle1 + net->offsetB[L - 1] + j;
            int index2 = nCycle2 + net->offsetB[L - 1] + j;
            foo = sigmaDot2_vec[j][b + a * (net->nI)] * net->sd_li[L][j];
            foo +=
                sigma_vec[j] * net->d2sd_li_dx2[L][b + a * (net->nI) +
                                                   j * (net->nI) * (net->nI)];
            foo += sigmaDot_vec[j][b] * net->dsd_li_dx[L][a + j * (net->nI)];
            foo += sigmaDot_vec[j][a] * net->dsd_li_dx[L][b + j * (net->nI)];
            net->parGrad_d2[index1] = foo;
            net->parGrad_d2[index2] = foo;
            for (int k = 0; k < net->arch[L - 1]; k++) {
              index1 =
                  nCycle1 + net->offsetW[L - 1] + k + (net->arch[L - 1]) * j;
              index2 =
                  nCycle2 + net->offsetW[L - 1] + k + (net->arch[L - 1]) * j;
              double foo2 = 0.0;
              foo2 += foo * net->s_li[L - 1][k];
              foo2 += sigmaDot_vec[j][a] * net->sd_li[L][j] *
                      net->ds_li_dx[L - 1][b + k * (net->nI)];
              foo2 += sigmaDot_vec[j][b] * net->sd_li[L][j] *
                      net->ds_li_dx[L - 1][a + k * (net->nI)];
              foo2 += sigma_vec[j] * net->dsd_li_dx[L][a + j * (net->nI)] *
                      net->ds_li_dx[L - 1][b + k * (net->nI)];
              foo2 += sigma_vec[j] * net->dsd_li_dx[L][b + j * (net->nI)] *
                      net->ds_li_dx[L - 1][a + k * (net->nI)];

              foo2 += sigma_vec[j] * net->sd_li[L][j] *
                      net->d2s_li_dx2[L - 1][b + a * (net->nI) +
                                             k * (net->nI) * (net->nI)];

              net->parGrad_d2[index1] = foo2;
              net->parGrad_d2[index2] = foo2;
              // net->parGrad[index] = net->parGrad[index] > 0 ? 1.0 : 1.0;
              // count++;
            }
          }
        }
      }

      for (int l = 0; l < net->arch[L - 1]; l++) {
        sigmaTemp_vec[l] = 0.0;
        for (int a = 0; a < net->nI; a++) {
          sigmaDotTemp_vec[l][a] = 0.0;
          for (int b = 0; b < net->nI; b++) {
            sigmaDot2Temp_vec[l][b + a * (net->nI)] = 0.0;
          }
        }

        for (int j = 0; j < net->arch[L]; j++) {
          sigmaTemp_vec[l] +=
              sigma_vec[j] * Wjk(j, l, L, net) * net->sd_li[L][j];
          for (int a = 0; a < net->nI; a++) {
            sigmaDotTemp_vec[l][a] +=
                sigmaDot_vec[j][a] * Wjk(j, l, L, net) * net->sd_li[L][j] +
                sigma_vec[j] * Wjk(j, l, L, net) *
                    net->dsd_li_dx[L][a + j * (net->nI)];
            for (int b = 0; b <= a; b++) {
              int index1 = b + a * (net->nI) + j * (net->nI) * (net->nI);
              double temp_loc = 0.0;
              temp_loc += sigmaDot2_vec[j][b + a * (net->nI)] *
                          Wjk(j, l, L, net) * net->sd_li[L][j];
              temp_loc += sigma_vec[j] * Wjk(j, l, L, net) *
                          net->d2sd_li_dx2[L][index1];
              temp_loc += sigmaDot_vec[j][b] * Wjk(j, l, L, net) *
                          net->dsd_li_dx[L][a + j * (net->nI)];
              temp_loc += sigmaDot_vec[j][a] * Wjk(j, l, L, net) *
                          net->dsd_li_dx[L][b + j * (net->nI)];
              sigmaDot2Temp_vec[l][b + a * (net->nI)] += temp_loc;
              sigmaDot2Temp_vec[l][a + b * (net->nI)] +=
                  (b == a ? 0.0 : temp_loc);
            }
          }
        }
      }
      for (int j = 0; j < net->arch[L - 1]; j++) {
        sigma_vec[j] = sigmaTemp_vec[j];
        for (int a = 0; a < net->nI; a++) {
          sigmaDot_vec[j][a] = sigmaDotTemp_vec[j][a];
          for (int b = 0; b < net->nI; b++) {
            sigmaDot2_vec[j][b + a * (net->nI)] =
                sigmaDot2Temp_vec[j][b + a * (net->nI)];
          }
        }
      }
    }
    // gradient w.r.t. the input as sigma L=0
    for (int a = 0; a < net->nI; a++) {
      net->res_d[a + (net->nI) * i] = sigma_vec[a];
      for (int b = 0; b < net->nI; b++) {
        net->res_d2[i][a + b * (net->nI)] =
            sigmaDot_vec[b][a] * net->sd_li[0][b] +
            sigma_vec[b] * net->dsd_li_dx[0][a + b * (net->nI)];
      }
    }
  }
}

void multilD_FullEvaluate(multilayerD *net, double *x) {
  if (x != NULL) {
    for (int i = 0; i < net->nI; i++) {
      net->input[i] = x[i];
    }
  }
  multilD_Evaluate(net);
  multilD_EvaluateParGradient(net);
}

double multilD_get_grad(int out, int par, multilayerD *net) {
  return net->parGrad[out * net->nPar + par];
}

double multilD_get_grad_d(int out, int i, int par, multilayerD *net) {
  return net->parGrad_d[i * (net->nPar * net->nO) + out * net->nPar + par];
}

double multilD_get_grad_d2(int out, int i1, int i2, int par, multilayerD *net) {
  return net->parGrad_d2[i2 * (net->nPar * net->nO * net->nI) +
                         i1 * (net->nPar * net->nO) + out * net->nPar + par];
}

double multilD_get(int out, multilayerD *net) {
  return net->s_li[net->nL - 1][out];
}

double multilD_get_d(int out, int i, multilayerD *net) {
  return net->res_d[i + out * (net->nI)];
}

double multilD_get_d2(int out, int i1, int i2, multilayerD *net) {
  return net->res_d2[out][i1 + i2 * (net->nI)];
}

#undef aC
#undef mC
#undef bC
